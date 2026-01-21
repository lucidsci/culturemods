"""
Unit conversion utilities for oxygen diffusion-reaction models.

Provides conversions between various units used in cell culture modeling:
- fmol, pmol, µmol, mol
- mm², mm³, L
- seconds, minutes, hours
"""


def fmols_per_mm2_per_s_to_mols_per_liter_per_hr(rate_fmols_per_mm2_per_s):
    """
    Convert flux from fmol/mm²/s to mol/L/hr.

    For a 1D cylindrical model:
    - 1 L = 1e6 mm³ (3D) = 1e6 mm² (cross-sectional area in 1D)
    - 1 fmol = 1e-15 mol
    - 1 hour = 3600 seconds
    """
    return 1e-15 * rate_fmols_per_mm2_per_s * 3600 * 1e6


def mols_per_liter_per_hr_to_fmols_per_mm2(rate_mols_per_L_per_hour):
    """
    Convert reaction rate from mol/L/hr to fmol/mm²/s.

    For a 1D cylindrical model:
    - 1 L = 1e6 mm³ (3D) = 1e6 mm² (cross-sectional area in 1D)
    - 1 mol = 1e15 fmol
    """
    rate_fmols_per_mm2_per_sec = 1e15 * rate_mols_per_L_per_hour / 1e6 / 3600
    return rate_fmols_per_mm2_per_sec


def fmols_per_mm3_to_micromolar(fmols_per_mm3):
    """Convert concentration from fmol/mm³ to µM (1 fmol/mm³ = 1 nM)."""
    return 1e3 * fmols_per_mm3


def rate_umolar_per_second_to_molar_per_hour(rate_umolar_per_s):
    """Convert rate from µM/s to M/hr."""
    return 3600 * rate_umolar_per_s / 1e6


def rate_mols_per_L_per_hour_umolar_per_s(rate):
    """Convert rate from mol/L/hr to µM/s."""
    return rate * 1e6 / 3600


def rate_pmols_per_L_per_minute_to_umolar_per_s(rate):
    """Convert rate from pmol/L/min to µM/s."""
    return rate * 60 / 1e3


def flux_fmols_per_mm2_per_s_to_umolar_per_s(flux, length_mm=1.0):
    """
    Convert flux in fmol/mm²/s to volumetric rate in µM/s.

    For a 1D model, flux (per area) becomes concentration rate (per volume)
    when divided by a characteristic length.

    1 fmol/mm³ = 0.001 µM (since 1 µM = 1000 fmol/mm³)
    So fmol/mm²/s ÷ mm = fmol/mm³/s → × 0.001 = µM/s

    Parameters
    ----------
    flux : float
        Flux in fmol/mm²/s
    length_mm : float
        Characteristic length (mm)

    Returns
    -------
    float
        Rate in µM/s
    """
    return flux / length_mm * 0.001
