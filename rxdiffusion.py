"""
discrete column simulation

references:
http://hplgit.github.io/num-methods-for-PDEs/doc/pub/diffu/

TODO - uncertainty eval based on height positional error, volume, concentration, D, and temperature
TODO - first order reaction rate model
"""
import numpy as np
# Numerical solution using finite difference method

class ReactionDiffusion1DParams:
    Nz = 100  # Number of spatial steps
    Nt = 50000  # Number of time steps
    dt = .1  # Time step size
    L = 3.1*3 #media height
    D = 3.2e-3
    C0=200
    Cs=200
    
    def depth_to_z_index(self, depth_mm):
        dz = self.L / (self.Nz - 1) 
        return int (depth // dz)

def params_from_well_model(csat=200, c_0=200, media_volume_uL=100, D=3.2e-3):
    pass

def R(z, t):
    return 0.01 #* np.exp(-t)# * np.sin(0.01 * t)  # Example: depth-dependent + time-varying

class ConstantRateProfile:
    def __init__(self, reaction_consumption_rate=0.1):
        self.rate = reaction_consumption_rate

    def reaction_at_height_time(self, z, t, c):
        return self.rate

class FirstOrderRateProfile:
    def __init__(self, reaction_consumption_rate=0.1):
        self.rate = reaction_consumption_rate

    def reaction_at_height_time(self, z, t, c):
        return self.rate

#simulates concentration over time based on 1D diffusion and reaction model
# profile of consumption/production at heights over time
# default to constant OCR across all heights
class ReactionDiffusion1DModel:
    def __init__(self, params: ReactionDiffusion1DParams):
        self.params = params

    def run_fdm(self, rate_profile=ConstantRateProfile()):
        p = self.params
        D = p.D
        dz = p.L / (p.Nz - 1)  # Spatial step size
        dt = p.dt
        L = p.L
        C0 = p.C0
        Cs = p.Cs
        Nz = p.Nz
        Nt = p.Nt
        # Stability condition (for explicit method)
        D_max = (dz**2) / (2 * D)
        if dt > D_max:
            dt = 0.5 * D_max  # Adjust time step for stability
    
        # Create grids
        z = np.linspace(0, L, Nz)
        t = np.linspace(0, Nt * dt, Nt)
        C = np.ones(Nz) * C0


        for n in range(Nt):
            C_new = C.copy()
            for i in range(1, Nz - 1):  # Exclude boundaries
                d2C_dz2 = (C[i + 1] - 2 * C[i] + C[i - 1]) / dz**2
                reaction_rate = rate_profile.reaction_at_height_time(z[i], t[n], C[i])
                C_new[i] = C[i] + dt * (D * d2C_dz2 - reaction_rate)

            # Apply boundary conditions
            C_new[0] = Cs  # Fixed concentration at the top
            C_new[-1] = C_new[-2]  # No-flux at the bottom

            C = C_new  # Update concentration
            yield C.copy()
            #C_results.append(C.copy())

def calc_parameterized(rates, volumes=[100], heights=[1], downsample_factor=100, Nt=60000):
    pts = []
    for media_vol in volumes:
        for ocr in rates:
            rate = flux_units_convert(ocr)
            rate_profile = ConstantRateProfile(rate)
            params = ReactionDiffusion1DParams()
            params.L = media_vol_to_height(media_vol)
            params.Nt = Nt
            model = ReactionDiffusion1DModel(params)
            media_height = media_vol_to_height(media_vol)
            dz = params.L / (params.Nz - 1)
            for t_i, concentrations in enumerate(model.run_fdm(rate_profile)):
                for h in heights:
                    depth = media_height - h
                    z_i = int (depth // dz)
                    t = t_i * params.dt

                    if t_i % downsample_factor == 0:
                        pt = {'ocr': ocr, 'rate': rate, 't_seconds': t, 
                              'c_at_z': concentrations[z_i], 'depth': depth, 'height': h, 'media_vol': media_vol}
                        pts.append(pt)
    return pts
    

if __name__ == '__main__':
    # Plot results
    import matplotlib.pyplot as plt
    from kinetics import media_vol_to_height, flux_units_convert
    
    if False:
        RAPID_REACTION_RATES = [200, 300, 400]
        pts = []
        for media_vol in [100, 200, 300]:
            for ocr in [5, 25, 50, 100, 200, 400]:
                rate = flux_units_convert(ocr)
                rate_profile = ConstantRateProfile(rate)
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

                        if t_i % 100 == 0:
                            pt = {'ocr': ocr, 'rate': rate, 't_seconds': t, 
                                  'c_at_z': concentrations[z_i], 'depth': d, 'height': h, 'media_vol': media_vol}
                            pts.append(pt)
        
    import pandas as pd
    import seaborn as sns
    
    pts = calc_parameterized([5, 25, 50, 100, 200, 400])
    df_all = pd.DataFrame(pts)
    df_all['t_mins'] = df_all['t_seconds'] / 60

    ax = sns.relplot(x='t_mins', y='c_at_z', hue='ocr', col='height', data=df_all, kind='line',  row='media_vol')
    plt.ylim(0, 210)
    plt.show()
    
    def _plot_by_rates(rates, volume=300, height=1, downsample_factor=100, Nt=60000):
        pts = calc_parameterized(rates, volumes=[volume], heights=[height], Nt=Nt)
        df_all = pd.DataFrame(pts)
        df_all['t_mins'] = df_all['t_seconds'] / 60

        ax = sns.lineplot(x='t_mins', y='c_at_z', hue='ocr', data=df_all)
        plt.ylim(0, None)
        plt.show()


    #df_rs = 
    if False:
        fig, ax = plt.subplots(figsize=(6, 4))
        dt_plot = 600
        for z_i in range(0, Nz, 10):
            ax.plot([i*dt/60 for i in range(0, Nt, dt_plot)], [ct[z_i] for ct in C_results[::dt_plot]])

        ax.set_xlabel('Time')
        ax.set_ylabel('C')
#ax.invert_yaxis()
        ax.legend()
        ax.set_title('1D Diffusion with Consumption Term')
        plt.show()
