"""
1D Reaction-Diffusion model for oxygen transport in cell culture media.

Models oxygen diffusion in a multiwell plate as a 1D system with:
- Fixed concentration at air-liquid interface (top)
- Zero-flux boundary at bottom (sealed well)
- Consumption term representing cellular oxygen uptake
    - models consumption as occurring uniformly throughout the volume
        (e.g. for suspension cultures, enzymatic reactions, etc. - not for monolayer)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Literal

from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, ImplicitSourceTerm


@dataclass
class SimulationConfig:
    """Configuration parameters for oxygen diffusion simulation."""

    # Physical parameters
    D: float = 3e-3  # Oxygen diffusion coefficient (mm²/s) (use ~2.2e-3 for 25C)
    C_air: float = 200.0  # Oxygen concentration at air interface (µM)

    # Geometry
    L: float = 3.1  # Depth of well (mm)
    nz: int = 200  # Number of mesh points

    # Reaction rate (dimensionless)
    k: float = 1.0  # zero order consumption rate

    k1: float = 0  # first order reaction rate

    # flux across bottom of well.  normally, zero to model
    # no flux across the plastic/glass well bottom barrier
    # however, to model a cell monolayer we can use this
    flux_bottom: float = 0

    # whether we are modeling an open air well (O2 air at top of media)
    # or a closed system (e.g. sealed at top)
    top_constraint: Literal['open', 'sealed'] = 'open'

    # Time parameters
    dt: Optional[float] = None  # Time step (s), computed from T if None
    steps: int = 100  # Number of time steps

    # Initial condition
    C_initial_fraction: float = 1.0  # Initial concentration as fraction of C_air

    # FIXME - no longer needed
    first_order_reaction: bool = False  # whether reaction is zero or first order

    halt_on_C_zero: bool = True  # halt simulation when C bottom drops below zero

    def __post_init__(self):
        """Compute derived parameters."""
        self.dz = self.L / self.nz
        self.T = self.L**2 / self.D  # Characteristic diffusion time (s)
        if self.dt is None:
            self.dt = self.T / 3600
        self.total_time = self.dt * self.steps

    @property
    def damkohler(self) -> float:
        """Damköhler number: ratio of reaction rate to diffusion rate."""
        return self.L**2 * self.k / self.D

    @property
    def L_char(self) -> float:
        """Characteristic diffusion-reaction length scale (mm)."""
        if self.k > 0:
            return np.sqrt(self.D / self.k)
        return float('inf')


@dataclass
class SimulationResult:
    """Results from a simulation run."""

    config: SimulationConfig
    points: List[Dict[str, Any]] = field(default_factory=list)
    final_profile: Optional[np.ndarray] = None

    def __getstate__(self) -> Dict[str, Any]:
        """Serialize for pickle/JSON. Converts numpy arrays to lists."""
        from dataclasses import asdict
        config_dict = asdict(self.config)
        # Remove derived fields computed in __post_init__
        for key in ['dz', 'T', 'total_time']:
            config_dict.pop(key, None)
        return {
            'config': config_dict,
            'points': self.points,
            'final_profile': self.final_profile.tolist() if self.final_profile is not None else None,
        }

    def __setstate__(self, state: Dict[str, Any]):
        """Deserialize from pickle/JSON."""
        self.config = SimulationConfig(**state['config'])
        self.points = state['points']
        self.final_profile = np.array(state['final_profile']) if state['final_profile'] is not None else None

    def to_dataframe(self):
        """Convert results to a pandas DataFrame with dimensionalized units."""
        import pandas as pd
        df = pd.DataFrame(self.points)
        if len(df) == 0:
            return df

        # Add dimensionalized units
        df['C'] = df['c_star'] * self.config.C_air
        df['height_mm'] = df['z_idx'] * self.config.dz
        df['t_hrs'] = df['step'] * self.config.dt / 3600
        df['t_mins'] = df['step'] * self.config.dt / 60
        df['t_s'] = df['step'] * self.config.dt
        return df

    def get_profile_at_step(self, step: int):
        """Get concentration profile at a specific time step."""
        import pandas as pd
        df = self.to_dataframe()
        return df[df['step'] == step].sort_values('z_idx')


def create_mesh(config: SimulationConfig) -> Grid1D:
    """Create the 1D mesh for the simulation."""
    return Grid1D(dx=config.dz, nx=config.nz)


def run_simulation(config: SimulationConfig,
                   record_every: int = 1,
                   verbose: bool = False,
                   step_callback=None) -> SimulationResult:
    """
    Run a 1D oxygen diffusion simulation.

    Parameters
    ----------
    config : SimulationConfig
        Simulation parameters
    record_every : int
        Record concentration profile every N steps (1 = every step)
    verbose : bool
        Print progress information
    step_callback : callable, optional
        Function called with step number during simulation

    Returns
    -------
    SimulationResult
        Object containing simulation results
    """
    mesh = create_mesh(config)

    # Create concentration variable (dimensionless, normalized by C_air)
    C = CellVariable(name="oxygen concentration", mesh=mesh,
                     value=config.C_initial_fraction)

    # Boundary conditions:
    # "right" = top of media (air-liquid interface for open system)
    if config.top_constraint == 'open':
        C.constrain(1.0, mesh.facesRight)  # Fixed concentration at air interface
    else:
        # closed/sealed system - enforce no flux boundary at top
        C.faceGrad.constrain((0,), where=mesh.facesRight)

    # "left" = bottom of well
    # flux at bottom is either zero (to model sealed plate bottom only)
    # or positive to represent a cell monolayer consumption
    C.faceGrad.constrain((config.flux_bottom,), where=mesh.facesLeft)

    eq = (TransientTerm() == DiffusionTerm(coeff=config.D) - config.k - ImplicitSourceTerm(coeff=config.k1))

    if verbose:
        print(f"Running simulation with k={config.k}, Da={config.damkohler:.3f}")
        print(f"Characteristic time: {config.T:.1f} s")
        print(f"Characteristic length: {config.L_char:.2f} mm")

    # Storage for results
    result = SimulationResult(config=config)

    def _record_step(C):
        for z_idx, c in enumerate(C.value):
            result.points.append({
                'k': config.k,
                'step': step,
                'c_star': c,
                'z_idx': z_idx
            })

    # Time stepping
    for step in range(config.steps):
        if step % record_every == 0:
            _record_step(C)
            if step_callback is not None:
                step_callback(step)

        eq.solve(var=C, dt=config.dt)

        if C.value[0] <= 0:
            if config.halt_on_C_zero:
                # halt simulation when C at bottom would go below zero
                C.value[0] = 0
                _record_step(C)
                break

    result.final_profile = C.value.copy()
    return result


def run_parameter_sweep(base_config: SimulationConfig,
                        param_name: str,
                        param_values: List[float],
                        record_every: int = 1,
                        verbose: bool = False) -> List[SimulationResult]:
    """
    Run simulations across a range of parameter values.

    Parameters
    ----------
    base_config : SimulationConfig
        Base configuration to modify
    param_name : str
        Name of parameter to sweep (must be a field of SimulationConfig)
    param_values : list
        Values to sweep over
    record_every : int
        Record every N steps
    verbose : bool
        Print progress

    Returns
    -------
    list of SimulationResult
        Results for each parameter value
    """
    from dataclasses import asdict

    results = []
    for val in param_values:
        # Create new config with modified parameter
        config_dict = asdict(base_config)
        # Remove derived fields that are computed in __post_init__
        for derived in ['dz', 'T', 'total_time']:
            config_dict.pop(derived, None)
        config_dict[param_name] = val
        config = SimulationConfig(**config_dict)

        if verbose:
            print(f"Running {param_name}={val}")

        result = run_simulation(config, record_every=record_every, verbose=verbose)
        results.append(result)

    return results


def combine_results(results: List[SimulationResult]):
    """Combine multiple simulation results into a single DataFrame."""
    import pandas as pd

    dfs = []
    for result in results:
        df = result.to_dataframe()
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)
