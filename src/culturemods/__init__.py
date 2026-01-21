"""
Culturemods - Oxygen diffusion-reaction models for cell culture systems.

This package provides tools for modeling oxygen transport and consumption
in cell culture systems, including:
- Analytical solutions for monolayer cultures (kinetics)
- Numerical solutions using FiPy for suspension cultures (rxd_fipy_1d)
- Unit conversion utilities (conversions)
- CLI tools for quick calculations
- GUI for interactive simulation
"""

from culturemods.kinetics import (
    concentration,
    flux_units_convert,
    media_vol_to_height,
    steady_state_o2_at_pos,
    time_to_constrained,
)
from culturemods.rxd_fipy_1d import (
    SimulationConfig,
    SimulationResult,
    run_simulation,
    run_parameter_sweep,
    combine_results,
)

__version__ = "0.1.0"

__all__ = [
    # kinetics
    "concentration",
    "flux_units_convert",
    "media_vol_to_height",
    "steady_state_o2_at_pos",
    "time_to_constrained",
    # rxd_fipy_1d
    "SimulationConfig",
    "SimulationResult",
    "run_simulation",
    "run_parameter_sweep",
    "combine_results",
]
