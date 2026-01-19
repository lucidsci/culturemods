"""
discrete column simulation

TODO - uncertainty eval based on height positional error, volume, concentration, D, and temperature
"""

import logging

try:
    from kinetics import media_vol_to_height
except:
    from culturemods.kinetics import media_vol_to_height

from rxd_fipy_1d import SimulationConfig, run_simulation, combine_results

# Module-level logger
logger = logging.getLogger(__name__)


def configure_logging(level=logging.INFO, format_string=None):
    """
    Configure logging for rxdiffusion module.

    Parameters
    ----------
    level : int
        Logging level (e.g., logging.DEBUG, logging.INFO)
    format_string : str, optional
        Custom format string for log messages
    """
    if format_string is None:
        format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(format_string))

    logger.setLevel(level)
    if not logger.handlers:
        logger.addHandler(handler)


class ReactionDiffusion1DParams:
    Nz = 100 # Number of spatial steps
    Nt = 50000  # Number of time steps
    dt = 2  # time step (seconds per step)
    L = 3.1 #media height mm
    D = 3e-3 #2e-3 #diffion coefficient mm2/s for 20-25C;  for 37C use ~3e-3
    C0=200
    Cs=200

    def validate(self):
        # Stability condition (for explicit method)
        dz = self.calc_dz()
        D_max = (dz**2) / (2 * self.D)
        logger.debug(f"Validating params: Nz={self.Nz}, Nt={self.Nt}, dt={self.dt:.4f}s, "
                     f"L={self.L:.2f}mm, D={self.D:.2e}mm2/s, dz={dz:.4f}mm")
        logger.debug(f"Stability check: dt={self.dt:.4f}s, D_max={D_max:.4f}s")
        if self.dt > D_max:
            logger.error(f"dt={self.dt} exceeds stability limit D_max={D_max}")
            raise Exception(f'dt too big - D max is {D_max}')
            self.dt = 0.5 * D_max  # Adjust time step for stability
        logger.debug("Parameter validation passed")

    def depth_to_z_index(self, depth_mm):
        dz = self.calc_dz()
        return int (depth_mm // dz)

    def calc_dz(self):
        #mm2 per discretized simulated height
        return self.L / (self.Nz - 1)

    def dz_volume_uL(self):
        #volume represented by each discreteized simulated column height step
        # since we are using mm2 for diffusion/flux and mm for height
        # the volume in uL is simply 1 mm2 times the discretized vertical step
        # dz
        return self.calc_dz() #dz in mm3 (uL)

    def characteristic_time(self):
        # characteristic time is driven by dimensionalized length (depth/height)
        # and diffusion constant D for our model
        return  self.L^2 / self.D

    def dimensionless_reaction_rate(self, rate_per_L_per_hour):
        #convert reaction rate in mols/L/hour to a dimensionless
        # reaction (consumption) rate (e.g. q*)
        # q* = q * L^2 / D = q * T
        # since L is in mm and D is mm2/s, we need the dimensionalized
        # rate q in units of mm2/s
        # a liter is 1e6 mm3 in 3D and 1e6 mm2 in 2D (cross sectional
        # area units for 1D cylindrical model).
        rate_per_mm2_per_sec  = rate_per_L_per_hour * 1e6 / 3600
        T = self.characteristic_time()
        q_star = rate_per_mm2_per_sec * T
        return q_star

    def normalized_reaction_rate(self, rate_mols_per_L_per_hour):
        #convert reaction rate in mols/L/hour to discretized units of
        # umols/uL/step
        rate_umolar_per_s = rate_mols_per_L_per_hour / 3600  * 1e6
        dC_dT = rate_umolar_per_s * self.dt * self.dz_volume_uL()
        return dC_dT #umolar per step

class ConstantRateUniformHeightProfile:
    first_order = False
    def __init__(self, reaction_consumption_rate=1):
        #reaction rate in mols/L/hour
        self.rate = reaction_consumption_rate

    def reaction_at_time(self, t, c):
        #reaction rate in mols/L/hour
        return self.rate

    def __repr__(self):
        return f'{self.rate}'


class FirstOrderRateProfile:
    first_order = True
    def __init__(self, k=1):
        #k is reaction rate cosntant in 1/s
        self.k = k

    def reaction_at_time(self, t, c):
        return c*self.k

    def __repr__(self):
        return f'k={self.k}'

class MichaelisMentenRateProfile:
    def __init__(self, v_max=10, K_m=100):
        self.v_max = v_max
        self.K_m = K_m

    #FIXME - get normalized rate k for first order sim
    def reaction_at_time(self, t, c):
        return (self.v_max * c) / (self.K_m + c)

    def __repr__(self):
        return f'MM V_max={self.v_max} K_m={self.K_m}'

#simulates concentration over time based on 1D diffusion and reaction model
# profile of consumption/production uniformly throughout all heights
class ReactionDiffusion1DModel:
    def __init__(self, params: ReactionDiffusion1DParams):
        logger.info(f"Initializing ReactionDiffusion1DModel: L={params.L:.2f}mm, "
                    f"Nz={params.Nz}, Nt={params.Nt}, dt={params.dt:.4f}s")
        logger.debug(f"Physical params: D={params.D:.2e}mm2/s, C0={params.C0}uM, Cs={params.Cs}uM")
        self.params = params

        #FIXME - let fipy model handle this as it controls whether
        # we use implicit, explicit, or crank-nicholson so
        # we don't know enough to determine if time step is too large
        #self.params.validate()

    def normalize_rate(self, rate_profile):
        """
        Convert rate profile to normalized rate k

        The FiPy model uses normalized concentration C* = C/C_sat, so the
        consumption rate k must be normalized: k = R / C_sat where R is in uM/s.

        For zero-order kinetics (constant rate):
            rate is in mols/L/hour -> R (uM/s) = rate * 1e6 / 3600
            k = R / C_sat

        For first-order kinetics:
            k is already in 1/s units (or similar)
        """
        p = self.params
        if rate_profile.first_order:
            # First-order rate constant, use directly
            k = rate_profile.k
            logger.debug(f"First-order rate profile: k={k} (1/s)")
            return k
        else:
            # Zero-order: convert mols/L/hour to normalized rate
            # rate (mols/L/hour) -> R (uM/s) = rate * 1e6 / 3600
            R_uM_per_s = rate_profile.rate * 1e6 / 3600
            # Normalize by C_sat to get dimensionless k
            k = R_uM_per_s / p.Cs if p.Cs > 0 else 0
            logger.debug(f"Zero-order rate profile: rate={rate_profile.rate} mols/L/hr -> "
                        f"R={R_uM_per_s:.4f} uM/s -> k={k:.6f} (normalized)")
            return k

    def run_fipy_result(self, rate_profile=ConstantRateUniformHeightProfile(), record_every=1):
        """
        Run simulation using FiPy solver and return full SimulationResult.

        Parameters
        ----------
        rate_profile : rate profile object
            Must have a `rate` or `k` attribute for the consumption rate.
        record_every : int
            Record concentration profile every N steps

        Returns
        -------
        SimulationResult
            Full result object with config, points, and final profile
        """

        logger.info(f"Starting FiPy simulation (result mode): profile={rate_profile}")

        p = self.params
        k = self.normalize_rate(rate_profile)


        config = SimulationConfig(
            D=p.D,
            C_air=p.Cs,
            L=p.L,
            nz=p.Nz,
            k=k,
            dt=p.dt,
            steps=p.Nt,
            C_initial_fraction=p.C0 / p.Cs if p.Cs > 0 else 1.0,
            first_order_reaction = rate_profile.first_order
        )

        total_time_hrs = (p.Nt * p.dt) / 3600
        logger.info(f"FiPy config: k={k:.6f}, steps={p.Nt}, total_time={total_time_hrs:.2f}hrs, "
                   f"record_every={record_every}")

        result = run_simulation(config, record_every=record_every, verbose=False)

        # Log final state
        if result.final_profile is not None:
            c_top = result.final_profile[-1] * p.Cs  # Last index is top (air interface)
            c_bottom = result.final_profile[0] * p.Cs  # First index is bottom
            logger.info(f"Simulation complete: final C_top={c_top:.2f}uM, C_bottom={c_bottom:.2f}uM")

        return result

    def run(self, rate_profile=ConstantRateUniformHeightProfile(), **kwargs):
        """
        Run simulation and return SimulationResult.

        Parameters
        ----------
        rate_profile : rate profile object
            Consumption rate profile
        **kwargs
            Additional arguments passed to run_fipy_result (e.g., record_every)

        Returns
        -------
        SimulationResult
            Full result object with config, points, and final profile
        """
        logger.debug(f"run() called with profile={rate_profile}, kwargs={kwargs}")
        return self.run_fipy_result(rate_profile, **kwargs)



def calc_parameterized_constant_rate(rates, volumes=[100], downsample_factor=100, T_minutes=60):
    """Run simulations with constant (zero-order) consumption rates."""
    logger.info(f"calc_parameterized_constant_rate: rates={rates}")
    profiles = [ConstantRateUniformHeightProfile(rate) for rate in rates]
    return calc_parameterized_profiles(profiles, volumes, downsample_factor, T_minutes)


def calc_parameterized_first_order(ks, volumes=[100], downsample_factor=100, T_minutes=60):
    """Run simulations with first-order consumption kinetics."""
    logger.info(f"calc_parameterized_first_order: k values={ks}")
    profiles = [FirstOrderRateProfile(k) for k in ks]
    return calc_parameterized_profiles(profiles, volumes, downsample_factor, T_minutes)


def calc_parameterized_profiles(rate_profiles, volumes=[100],  downsample_factor=100, duration_minutes=30):
    """
    Run simulations across rate profiles, volumes

    Parameters
    ----------
    rate_profiles : list
        List of rate profile objects
    volumes : list
        Media volumes in uL
    downsample_factor : int
        Record every N steps

    Returns
    -------
    list of dict
        Data points with rate, time, concentration, etc.
    """
    dt = 1 #seconds
    Nt = int(duration_minutes * 60 / dt)

    total_sims = len(volumes) * len(rate_profiles)
    logger.info(f"Starting parameterized simulation sweep: {len(rate_profiles)} profiles x "
               f"{len(volumes)} volumes = {total_sims} simulations")
    logger.info(f"Parameters: volumes={volumes}uL, Nt={Nt}, downsample_factor={downsample_factor}")
    logger.debug(f"Rate profiles: {rate_profiles}")

    dfs = []
    sim_count = 0

    for media_vol in volumes:
        media_height = media_vol_to_height(media_vol)
        for rate_profile in rate_profiles:
            sim_count += 1
            logger.info(f"Running simulation {sim_count}/{total_sims}: "
                       f"volume={media_vol}uL, profile={rate_profile}")

            params = ReactionDiffusion1DParams()
            params.L = media_vol_to_height(media_vol)
            params.Nt = Nt

            #target a strandard dz of around every 10 microns
            L_um = round(params.L * 1000)
            params.Nz  = L_um/10 - 1

            dz = params.calc_dz()

            logger.info(f"dz {dz*1000} microns")

            logger.debug(f"Media height for {media_vol}uL: {params.L:.2f}mm")

            model = ReactionDiffusion1DModel(params)

            sim_result = model.run(rate_profile, record_every=downsample_factor)
            df = sim_result.to_dataframe()



            df['profile'] = str(rate_profile)
            df['media_height'] = media_height
            df['media_vol'] = media_vol

            df['z'] = df.z_idx * dz

            dfs.append(df)

            logger.debug(f"Simulation {sim_count} complete: {len(sim_result.points)} points")

    return pd.concat(dfs)


if __name__ == '__main__':
    # Configure logging for demo
    configure_logging(level=logging.INFO)

    # Plot results
    import matplotlib.pyplot as plt
    from kinetics import media_vol_to_height

    logger.info("Starting rxdiffusion demo")

    #example reaction rates in mols/L/hr
    #ALGAL_REACTION_RATES = list(range(1,45, 8))
    ALGAL_REACTION_RATES = [5e-3, 1e-2, 2e-2, 3e-1]

    ENZYMATIC_REACTION_RATES = [1e-5*v for v in range(1, 200, 10)]#, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3]

    #mid range - peak rate reaches zero around 15 minutes, low rate stablizes at non-zero
    MID_RANGE_REACTION_RATES = [1e-5*v for v in range(10, 600, 10)]
    FULL_RANGE_REACTION_RATES = [1e-5*v for v in range(5, 600, 2)]


    import pandas as pd
    import seaborn as sns
    import numpy as np

    VOLS = [300]
    #df_all = calc_parameterized_constant_rate(FULL_RANGE_REACTION_RATES,
    #                             volumes=[100, 200, 300], downsample_factor=10)
    BACTERIAL_RATES = [1e-4* v for v in range(1,30,5 )]
    df_all = calc_parameterized_constant_rate(BACTERIAL_RATES,
                                 volumes=[100, 200, 300], downsample_factor=10)
    df_all['t_mins'] = df_all['t_s'] / 60
    df_all['z_um'] = df_all.z.apply(lambda z: round(z*1000))

    probe_height = 1
    #find closest z in simulation to probe height
    zs = df_all.z.unique()
    probe_z = zs[np.abs(zs - probe_height).argmin()]
    df_at_bottom = df_all[df_all.z == 0]
    df_at_pos = df_all[df_all.z == probe_z]
    ax = sns.relplot(x='t_mins', y='C', hue='profile',  data=df_at_bottom, kind='line',  row='media_vol')
    plt.ylim(0, 210)
    plt.show()

    if False:
        pts = calc_parameterized_constant_rate(ALGAL_REACTION_RATES,
                                 volumes=[300], downsample_factor=60*15, Nt=3600*10*12)
        df_all = pd.DataFrame(pts)
        df_all['t_mins'] = df_all['t_seconds'] / 60
        df_all['t_hrs'] = df_all['t_seconds'] / 3600

        ax = sns.relplot(x='t_hrs', y='c_at_z', hue='rate', col='height', data=df_all, kind='line',  row='media_vol')
        plt.ylim(0, 210)
        plt.show()

        def _plot_by_rates(rates, volume=300, height=1, downsample_factor=100, Nt=60000):
            pts = calc_parameterized_constant_rate(rates, volumes=[volume])
            df_all = pd.DataFrame(pts)
            df_all['t_mins'] = df_all['t_seconds'] / 60

            sns.lineplot(x='t_mins', y='c_at_z', hue='profile', data=df_all)
            plt.ylim(0, None)
            plt.show()

        ALGAL_PRODUCTION_RATES = [-5e-3, -1e-2, -2e-2, -3e-1]


    #_plot_by_rates(ENZYMATIC_REACTION_RATES, volume)
    #df_rs =
    if False:
        fig, ax = plt.subplots(figsize=(6, 4))
        dt_plot = 600
        for _i in range(0, Nz, 10):
            ax.plot([i*dt/60 for i in range(0, Nt, dt_plot)], [ct[z_i] for ct in C_results[::dt_plot]])

        ax.set_xlabel('Time')
        ax.set_ylabel('C')
#ax.invert_yaxis()
        ax.legend()
        ax.set_title('1D Diffusion with Consumption Term')
        plt.show()
