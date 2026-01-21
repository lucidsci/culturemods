"""
Command-line interface for oxygen diffusion kinetics calculations.

Provides quick access to common calculations and visualizations.
"""
import math

import click
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from culturemods import kinetics

matplotlib.style.use('fivethirtyeight')


def media_vol_to_height(vol_uL, well_radius_mm=3.2):
    """Convert media volume to height for a cylindrical well."""
    return vol_uL / (math.pi * well_radius_mm**2)


@click.group()
def cmds():
    """Culturemods CLI - Oxygen diffusion calculations for cell culture."""
    pass


@cmds.command()
@click.argument('ocrs', nargs=-1, type=int)
@click.option('--vol', default=100, type=int, help='Media volume in µL')
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def o2_at_bottom_by_ocr(ocrs, vol=100, csat=185):
    """Plot O2 at bottom over time for different OCR values."""
    height = media_vol_to_height(vol)
    if len(ocrs) == 0:
        ocrs = [100, 150, 200, 250]

    print(f"OCR values: {ocrs}")
    for ocr in ocrs:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600 * 4, 1))
        cs = [kinetics.concentration(0, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t / 3600 for t in ts], cs, label=f"OCR={ocr}")
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (µM)")
        plt.ylim(0, csat)

    plt.legend()
    plt.title(f"O2 at bottom by OCR in {vol} µL")
    plt.show()


@cmds.command()
@click.argument('vols', nargs=-1, type=int)
@click.option('--ocr', default=100, type=int, help='OCR in fmol/mm²/s')
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def o2_at_bottom_by_vol(vols, ocr=100, csat=185):
    """Plot O2 at bottom over time for different media volumes."""
    if len(vols) == 0:
        vols = [50, 100, 150, 200]

    q = kinetics.flux_units_convert(ocr)
    for vol in vols:
        height = media_vol_to_height(vol)
        ts = list(range(0, 3600 * 4, 1))
        cs = [kinetics.concentration(0, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t / 3600 for t in ts], cs, label=f"vol={vol}µL")
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (µM)")
        plt.ylim(0, csat)

    plt.legend()
    plt.title(f"O2 at bottom by volume (OCR={ocr})")
    plt.show()


@cmds.command()
@click.option('--ocr', default=100, type=int, help='OCR in fmol/mm²/s')
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def o2_at_heights(ocr=100, csat=185):
    """Plot O2 at different heights over time."""
    heights_um = list(range(1000, 1500, 100))
    vol = 100
    q = kinetics.flux_units_convert(ocr)
    for h in heights_um:
        height = media_vol_to_height(vol)
        ts = list(range(0, 3600 * 4, 1))
        cs = [kinetics.concentration(h / 1000, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t / 3600 for t in ts], cs, label=f"position={h}µm")
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (µM)")
        plt.ylim(0, csat)

    plt.legend()
    plt.title(f"O2 at position (OCR={ocr})")
    plt.show()


@cmds.command()
@click.argument('vols', nargs=-1, type=int)
@click.option('--position', default=1.25, help='Position above bottom in mm')
@click.option('--ocr', default=100, type=int, help='OCR in fmol/mm²/s')
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def o2_at_position_by_vol(vols, position=1.25, ocr=100, csat=185):
    """Plot O2 at a specific position for different media volumes."""
    if len(vols) == 0:
        vols = [50, 100, 150, 200]

    q = kinetics.flux_units_convert(ocr)
    for vol in vols:
        height = media_vol_to_height(vol)
        ts = list(range(0, 3600 * 4, 1))
        cs = [kinetics.concentration(position, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t / 3600 for t in ts], cs, label=f"vol={vol}µL")
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (µM)")
        plt.ylim(0, csat)

    plt.legend()
    plt.title(f"O2 at {position}mm by volume (OCR={ocr})")
    plt.show()


@cmds.command()
@click.argument('ocrs', nargs=-1, type=int)
@click.option('--position', default=1.25, help='Position above bottom in mm')
@click.option('--vol', default=100, type=int, help='Media volume in µL')
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def o2_at_position_by_ocr(ocrs, position=1.25, vol=100, csat=185):
    """Plot O2 at a specific position for different OCR values."""
    height = media_vol_to_height(vol)
    for ocr in ocrs:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600 * 4, 1))
        cs = [kinetics.concentration(position, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t / 3600 for t in ts], cs, label=f"OCR={ocr}")
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (µM)")
        plt.ylim(0, csat)

    plt.legend()
    plt.title(f"O2 at {position}mm above cells by OCR in {vol} µL")
    plt.show()


def plot_gradient_evolution(ocr=100, vol=100, csat=185, tmax_hrs=1, delta_t_mins=15):
    """Plot O2 gradient evolution over time."""
    nz = 20
    media_height = media_vol_to_height(vol)

    for t in range(0, 3600 * tmax_hrs, delta_t_mins * 60):
        cs = []
        hs = []
        t_min = t // 60
        for z_i in range(nz):
            h = z_i * media_height / nz
            q = ocr * 1e-3
            c = kinetics.concentration(h, t, q, c_initial=csat, media_height=media_height)
            cs.append(c)
            hs.append(h)
        plt.plot(hs, cs, label=f"t={t_min} mins")
        plt.legend()
    plt.title(f"Gradient over Time (OCR={ocr} fmol/mm²/s, volume={vol}µL)")


def plot_time_to_diffusion_limit_vs_ocr(ocrs, vols=[100, 150, 200], csat=185):
    """Plot time to diffusion limit vs OCR for different volumes."""
    for vol in vols:
        times = [kinetics.time_to_constrained(ocr, vol) for ocr in ocrs]
        t_hrs = [t / 3600 if t is not None else None for t in times]
        plt.plot(ocrs, t_hrs, label=f"{vol} µL")
    plt.xlabel("OCR (fmol/mm²/s)")
    plt.ylabel("Time to Diffusion-Limit (hrs)")
    plt.legend()
    plt.title("Time to OCR Diffusion Limit")


def plot_time_to_diffusion_limit_vs_vol(vols, ocrs=[100, 200, 300, 400], csat=185):
    """Plot time to diffusion limit vs volume for different OCRs."""
    for ocr in ocrs:
        times = [kinetics.time_to_constrained(ocr, vol) for vol in vols]
        t_hrs = [t / 3600 if t is not None else None for t in times]
        plt.plot(vols, t_hrs, label=f"OCR={ocr}")
    plt.xlabel("Volume (µL)")
    plt.ylabel("Time to Diffusion-Limit (hrs)")
    plt.legend()
    plt.title("Time to OCR Diffusion Limit")


@cmds.command()
@click.argument('vols', nargs=-1, type=int)
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def time_to_diffusion_limit_vs_ocr(vols, csat=185):
    """Plot time to diffusion limit vs OCR for given volumes."""
    if len(vols) == 0:
        vols = [100, 150, 200]
        print(f"Defaulting to volumes {vols}")
    ocrs = list(range(50, 400, 10))
    print("Calculating...")
    plot_time_to_diffusion_limit_vs_ocr(ocrs=ocrs, vols=vols, csat=csat)
    plt.show()


@cmds.command()
@click.argument('ocrs', nargs=-1, type=int)
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def time_to_diffusion_limit_vs_vol(ocrs, csat=185):
    """Plot time to diffusion limit vs volume for given OCRs."""
    if len(ocrs) == 0:
        ocrs = [100, 200, 300, 400]
        print(f"Defaulting to OCRs {ocrs}")

    vols = list(range(50, 300, 10))
    print("Calculating...")
    plot_time_to_diffusion_limit_vs_vol(ocrs=ocrs, vols=vols, csat=csat)
    plt.show()


@cmds.command()
@click.argument('ocr', default=100)
@click.option('--vol', default=100, type=int, help='Media volume in µL')
@click.option('--csat', default=185, type=int, help='Saturated O2 concentration in µM')
def o2_gradient_over_time(ocr=100, vol=100, csat=185, tmax_hrs=1, delta_t_mins=15):
    """Plot O2 gradient evolution over time."""
    plot_gradient_evolution(ocr=ocr, vol=vol, csat=csat, tmax_hrs=tmax_hrs, delta_t_mins=delta_t_mins)
    plt.show()


if __name__ == '__main__':
    cmds()
