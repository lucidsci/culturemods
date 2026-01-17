"""
discrete column simulation

TODO - uncertainty eval based on height positional error, volume, concentration, D, and temperature
"""

try:
    from kinetics import media_vol_to_height
except:
    from culturemods.kinetics import media_vol_to_height

from rxd_fipy_1d import SimulationConfig, run_simulation, combine_results


class ReactionDiffusion1DParams:
    Nz = 100 # Number of spatial steps
    Nt = 50000  # Number of time steps
    dt = .15  # time step (seconds per step)
    L = 3.1 #media height mm
    D = 3e-3 #2e-3 #diffion coefficient mm2/s for 20-25C;  for 37C use ~3e-3
    C0=200
    Cs=200

    def validate(self):
        # Stability condition (for explicit method)
        dz = self.calc_dz()
        D_max = (dz**2) / (2 * self.D)
        if self.dt > D_max:
            raise Exception(f'dt too big - D max is {D_max}')
            self.dt = 0.5 * D_max  # Adjust time step for stability

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
        rate_per_mm2_per_sec  = rate_per_L_per_hour * 1e6 * 3600
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
    def __init__(self, reaction_consumption_rate=1):
        #reaction rate in mols/L/hour
        self.rate = reaction_consumption_rate

    def reaction_at_time(self, t, c):
        #reaction rate in mols/L/hour
        return self.rate

    def __repr__(self):
        return f'{self.rate}'


class FirstOrderRateProfile:
    def __init__(self, k=1):
        self.k = k

    def reaction_at_time(self, t, c):
        return c*self.k

    def __repr__(self):
        return f'k={self.k}'

class MichaelisMentenRateProfile:
    def __init__(self, v_max=10, K_m=100):
        self.v_max = v_max
        self.K_m = K_m

    def reaction_at_time(self, t, c):
        return (self.v_max * c) / (self.K_m + c)

    def __repr__(self):
        return f'MM V_max={self.v_max} K_m={self.K_m}'

#simulates concentration over time based on 1D diffusion and reaction model
# profile of consumption/production uniformly throughout all heights
class ReactionDiffusion1DModel:
    def __init__(self, params: ReactionDiffusion1DParams):
        self.params = params
        self.params.validate()

    def _rate_to_fipy_k(self, rate_profile):
        """
        Convert rate profile to FiPy k parameter.

        The FiPy model uses normalized concentration C* = C/C_air, so the
        consumption rate k must be normalized: k = R / C_air where R is in µM/s.

        For zero-order kinetics (constant rate):
            rate is in mols/L/hour -> R (µM/s) = rate * 1e6 / 3600
            k = R / C_air

        For first-order kinetics:
            k is already in 1/s units (or similar)
        """
        p = self.params
        if hasattr(rate_profile, 'k'):
            # First-order rate constant, use directly
            return rate_profile.k
        elif hasattr(rate_profile, 'rate'):
            # Zero-order: convert mols/L/hour to normalized rate
            # rate (mols/L/hour) -> R (uM/s) = rate * 1e6 / 3600
            R_uM_per_s = rate_profile.rate * 1e6 / 3600
            # Normalize by C_air to get dimensionless k
            k = R_uM_per_s / p.Cs if p.Cs > 0 else 0
            return k
        else:
            raise ValueError(f"Unsupported rate profile type: {type(rate_profile)}")

    def run_fipy(self, rate_profile=ConstantRateUniformHeightProfile(), record_every=10):
        """
        Run simulation using FiPy solver.

        Parameters
        ----------
        rate_profile : rate profile object
            Must have a `rate` or `k` attribute for the consumption rate.
            Currently supports ConstantRateUniformHeightProfile and FirstOrderRateProfile.
        record_every : int
            Yield concentration profile every N steps

        Yields
        ------
        np.ndarray
            Concentration profile at each recorded time step (bottom to top)
        """

        p = self.params
        k = self._rate_to_fipy_k(rate_profile)

        # Create FiPy simulation config from params
        config = SimulationConfig(
            D=p.D,
            C_air=p.Cs,
            L=p.L,
            nz=p.Nz,
            k=k,
            dt=p.dt,
            steps=p.Nt,
            C_initial_fraction=p.C0 / p.Cs if p.Cs > 0 else 1.0
        )

        # Run FiPy simulation
        result = run_simulation(config, record_every=record_every, verbose=False)

        # Convert results to match run_fdm output format
        # FiPy stores z_idx=0 as bottom (left), z_idx=nz-1 as top (right/air)
        # run_fdm stores index 0 as top (air), index -1 as bottom
        # So we need to reverse the profile order
        df = result.to_dataframe()
        for step in range(0, p.Nt + 1, record_every):
            step_data = df[df['step'] == step].sort_values('z_idx')
            if len(step_data) > 0:
                # Reverse to match FDM convention (top=0, bottom=-1)
                # and convert from dimensionless to µM
                profile = step_data['C'].values[::-1]
                yield profile

    def run_fipy_result(self, rate_profile=ConstantRateUniformHeightProfile()):
        """
        Run simulation using FiPy solver and return full SimulationResult.

        Parameters
        ----------
        rate_profile : rate profile object
            Must have a `rate` or `k` attribute for the consumption rate.

        Returns
        -------
        SimulationResult
            Full result object with config, points, and final profile
        """
        from rxd_fipy_1d import SimulationConfig, run_simulation

        p = self.params
        k = self._rate_to_fipy_k(rate_profile)

        config = SimulationConfig(
            D=p.D,
            C_air=p.Cs,
            L=p.L,
            nz=p.Nz,
            k=k,
            dt=p.dt,
            steps=p.Nt,
            C_initial_fraction=p.C0 / p.Cs if p.Cs > 0 else 1.0
        )

        return run_simulation(config, record_every=1, verbose=False)

    def run(self, rate_profile=ConstantRateUniformHeightProfile(), **kwargs):
        """

        Parameters
        ----------
        rate_profile : rate profile object
            Consumption rate profile
        solver : str
            'fdm' for finite difference method, 'fipy' for FiPy solver
        **kwargs
            Additional arguments passed to the solver method

        Yields
        ------
        np.ndarray
            Concentration profile at each time step
        """
        return self.run_fipy_result(rate_profile, **kwargs)



def calc_parameterized_constant_rate(rates, volumes=[100], heights=[1], downsample_factor=100, Nt=60000, solver='fdm'):
    profiles = [ConstantRateUniformHeightProfile(rate) for rate in rates]
    return calc_parameterized_profiles(profiles, volumes, heights, downsample_factor, Nt, solver=solver)


def calc_parameterized_first_order(ks, volumes=[100], heights=[1], downsample_factor=100, Nt=60000, solver='fdm'):
    profiles = [FirstOrderRateProfile(k) for k in ks]
    return calc_parameterized_profiles(profiles, volumes, heights, downsample_factor, Nt, solver=solver)


def calc_parameterized_profiles(rate_profiles, volumes=[100], heights=[1], downsample_factor=100, Nt=60000, solver='fdm'):
    """
    Run simulations across rate profiles, volumes, and heights.

    Parameters
    ----------
    rate_profiles : list
        List of rate profile objects
    volumes : list
        Media volumes in µL
    heights : list
        Heights from bottom to probe (mm)
    downsample_factor : int
        Record every N steps
    Nt : int
        Number of time steps
    solver : str
        'fdm' for finite difference method, 'fipy' for FiPy solver

    Returns
    -------
    list of dict
        Data points with rate, time, concentration, etc.
    """
    sim_results = []
    for media_vol in volumes:
        for rate_profile in rate_profiles:
            params = ReactionDiffusion1DParams()
            params.L = media_vol_to_height(media_vol)
            params.Nt = Nt
            model = ReactionDiffusion1DModel(params)
            media_height = media_vol_to_height(media_vol)
            dz = params.L / (params.Nz - 1)
            profile_str = str(rate_profile)

            sim_result  = model.run(rate_profile, record_every=downsample_factor)
            sim_results.append(sim_result)

    return combine_results(sim_results)


if __name__ == '__main__':
    # Plot results
    import matplotlib.pyplot as plt
    from kinetics import media_vol_to_height

    #example reaction rates in mols/L/hr
    #ALGAL_REACTION_RATES = list(range(1,45, 8))
    ALGAL_REACTION_RATES = [5e-3, 1e-2, 2e-2, 3e-1]

    ENZYMATIC_REACTION_RATES = [1e-5*v for v in range(5, 30, 5)]#, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3]


    if False:
        RAPID_REACTION_RATES = [1e-3*v for v in range(1, 100, 10)]
        pts = []
        for media_vol in [100, 200, 300]:
            for rate in [5, 25, 50, 100, 200, 400]:
                rate_profile = ConstantRateUniformHeightProfile(rate)
                params = ReactionDiffusion1DParams()
                params.L = media_vol_to_height(media_vol)
                model = ReactionDiffusion1DModel(params)
                media_height = media_vol_to_height(media_vol)
                dz = params.L / (params.Nz - 1)
                depths = [0.5, 1, 1.5, 2, 2.5]
                for t_i, concentrations in enumerate(model.run_fdm(rate_profile)):
                    for d in depths:
                        #z_i = int(h // params.Nz)
                        z_i = int (d // dz)
                        t = t_i * params.dt

                        #height from bottom
                        h = media_height - d

                        rate = rate_profile.reaction_at_time(t)
                        if t_i % 100 == 0:
                            pt = {'rate': rate, 't_seconds': t,
                                  'c_at_z': concentrations[z_i], 'depth': d, 'height': h, 'media_vol': media_vol}
                            pts.append(pt)

    import pandas as pd
    import seaborn as sns

    df_all = calc_parameterized_constant_rate(ENZYMATIC_REACTION_RATES,
                                 volumes=[100], heights=[1,2], downsample_factor=60*15, Nt=3600*12)
    df_all['t_mins'] = df_all['t_seconds'] / 60
    df_all['t_hrs'] = df_all['t_seconds'] / 3600

    ax = sns.relplot(x='t_mins', y='c_at_z', hue='profile', col='height', data=df_all, kind='line',  row='media_vol')
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
            pts = calc_parameterized_constant_rate(rates, volumes=[volume], heights=[height], Nt=Nt)
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
