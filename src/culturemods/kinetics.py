"""
Analytical solution for gas diffusion in mammalian cell culture systems.

Based on "Kinetics of Gas Diffusion in Mammalian Cell Culture Systems"
McLimans, Blumenson, Tunnah. 1968. Equations from Section 23.
"""
import math

# Constants
MONOLAYER_THICKNESS_MM = 0
DEFAULT_HEIGHT_MM = 3.1
O2_DIFFUSION_37C = 3.2e-3  # mm²/s


def H(n, a, s=MONOLAYER_THICKNESS_MM):
    """Calculate H_n from eq 19 where a is fluid height and s is monolayer thickness."""
    return math.pi * (2 * n + 1) / (2 * (a - s))


def summation(x, t, a=DEFAULT_HEIGHT_MM, s=MONOLAYER_THICKNESS_MM, iterations=100, D=O2_DIFFUSION_37C):
    """Calculate summand from eq 23 for the given n."""
    sigma = 0
    for n in range(iterations):
        H_n = H(n, a, s)
        sigma += math.pow(H_n, -2) * math.cos(H_n * x) * math.exp(-1 * math.pow(H_n, 2) * D * t)
    return sigma


def concentration(x, t, Q=5, c_initial=200, media_height=DEFAULT_HEIGHT_MM):
    """
    Calculate concentration over time at given height/depth based on Eq. 23.

    Parameters
    ----------
    x : float
        Position above bottom (mm)
    t : float
        Time (seconds)
    Q : float
        Flux in µmol/mm²/s
    c_initial : float
        Initial O2 concentration (µM)
    media_height : float
        Height of media column (mm)

    Returns
    -------
    float
        O2 concentration at position x and time t (µM)
    """
    a = media_height
    s = MONOLAYER_THICKNESS_MM
    D = O2_DIFFUSION_37C

    c_x_t = c_initial - (a - s - x) * (Q / D) + (2 * Q / ((a - s) * D)) * summation(x, t, s=s, a=a, D=D)
    return c_x_t


def flux_units_convert(flux_fmols_mm_sq_per_sec):
    """
    Convert flux from fmol/mm²/s to µmol/mm²/s.

    Parameters
    ----------
    flux_fmols_mm_sq_per_sec : float
        Flux in fmol/mm²/s

    Returns
    -------
    float
        Flux in µmol/mm²/s
    """
    return flux_fmols_mm_sq_per_sec * 1e-3


def media_vol_to_height(vol_uL, well_radius_mm=3.2):
    """
    Convert media volume to height for a cylindrical well.

    Parameters
    ----------
    vol_uL : float
        Volume in µL
    well_radius_mm : float
        Well radius in mm (default 3.2mm for 96-well)

    Returns
    -------
    float
        Media height in mm
    """
    return vol_uL / (math.pi * well_radius_mm**2)


def time_to_constrained(ocr, media_vol_uL=100, csat=185):
    """
    Calculate time until O2 at bottom reaches zero (diffusion-limited).

    Parameters
    ----------
    ocr : float
        Oxygen consumption rate in fmol/mm²/s
    media_vol_uL : float
        Media volume in µL
    csat : float
        Saturated O2 concentration (µM)

    Returns
    -------
    int or None
        Time in seconds, or None if not reached within 8 hours
    """
    h = media_vol_to_height(media_vol_uL)
    q = ocr * 1e-3
    for t_s in range(8 * 3600):
        c = concentration(0, t_s, q, c_initial=csat, media_height=h)
        if c <= 0:
            return t_s
    return None


def max_flux(vol=100, sat=200, D=3.2e-03):
    """Calculate maximum sustainable flux for given volume."""
    h = media_vol_to_height(vol)
    dC = sat
    return 1e3 * D * dC / h


def steady_state_o2_at_pos(flux, D=3.2e-03, o2_sat=200, media_height=3.1, pos_above_bottom=0):
    """
    Calculate steady-state O2 concentration at a given position.

    Parameters
    ----------
    flux : float
        O2 flux in fmol/mm²/s
    D : float
        Diffusion coefficient (mm²/s)
    o2_sat : float
        Saturated O2 concentration (µM)
    media_height : float
        Media height (mm)
    pos_above_bottom : float
        Position above bottom (mm)

    Returns
    -------
    float
        Steady-state O2 concentration (µM)
    """
    dC = (media_height - pos_above_bottom) * flux / D / 1000
    return min(o2_sat - dC, o2_sat)
