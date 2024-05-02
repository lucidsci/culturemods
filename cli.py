import kinetics
import click
import math
import matplotlib 
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.style.use('fivethirtyeight')

@click.group()
def cmds():
    pass

def media_vol_to_height(vol_uL, well_radius_mm=3.2):
    return vol_uL / (math.pi * well_radius_mm**2)

@cmds.command()
@click.argument('ocrs', nargs=-1, type=int)
@click.option('--vol', default=100, type=int)
@click.option('--csat', default=185, type=int)
def o2_at_bottom_by_ocr(ocrs, vol=100, csat=185):
    height = media_vol_to_height(vol)
    if len(ocrs) == 0:
        ocrs = [100, 150, 200, 250]

    print(ocrs)
    for ocr in ocrs:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600*4, 1))
        cs = [kinetics.concentration(0, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="OCR={}".format(ocr))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, csat)
    
    plt.legend()
    plt.title("O2 at bottom by OCR in {} uL".format(vol))
    plt.show()

@cmds.command()
@click.argument('ocrs', nargs=-1, type=int)
@click.option('--position', default=1.25)
@click.option('--vol', default=100, type=int)
@click.option('--csat', default=185, type=int)
def o2_at_position_by_ocr(ocrs, position=1.25, vol=100, csat=185):
    height = media_vol_to_height(vol)
    for ocr in ocrs:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600*4, 1))
        cs = [kinetics.concentration(position, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="OCR={}".format(ocr))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, csat)
    
    plt.legend()
    plt.title("O2 at {}mm above cells by OCR in {} uL".format(position, vol))
    plt.show()

def plot_o2_at_position_by_ocr(ocrs, position=1.25, vol=100, csat=185):
    height = media_vol_to_height(vol)
    for ocr in ocrs:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600*4, 1))
        cs = [kinetics.concentration(position, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="OCR={}".format(ocr))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, csat)
    
    plt.legend()
    plt.title("O2 at {}mm above cells by OCR in {} uL".format(position, vol))
    plt.show()

def plot_gradient_evolution(ocr=100, vol=100, csat=185, tmax_hrs=1, delta_t_mins=15):
    nz = 20
    media_height = media_vol_to_height(vol)

    for t in range(0, 3600*tmax_hrs, delta_t_mins*60):#[3600*hr for hr in range(5)]:
        cs = []
        hs = []
        t_min = t//60
        for z_i in range(nz):
            h = z_i * media_height / nz
            #q = flux_units_convert(ocr)
            #print(q)
            q = ocr * 1e-3  #fmtomols/mm3 to micromolar
            #ts = list(range(0, 3600*24, 10))
            c = kinetics.concentration(h, t, q, c_initial=csat, media_height=media_height)
            cs.append(c)
            hs.append(h)
        plt.plot(hs, cs, label="t={} mins".format(t_min))
        plt.legend()
        #plt.ylim(0, 200)
    plt.title("Gradient over Time (OCR={} fmols/mm2/s  media volume={}uL".format(ocr, vol))

def plot_gradient_evolution(ocr=100, vol=100, csat=185, tmax_hrs=1, delta_t_mins=15):
    nz = 20
    media_height = media_vol_to_height(vol)

    for t in range(0, 3600*tmax_hrs, delta_t_mins*60):#[3600*hr for hr in range(5)]:
        cs = []
        hs = []
        t_min = t//60
        for z_i in range(nz):
            h = z_i * media_height / nz
            #q = flux_units_convert(ocr)
            #print(q)
            q = ocr * 1e-3  #fmtomols/mm3 to micromolar
            #ts = list(range(0, 3600*24, 10))
            c = kinetics.concentration(h, t, q, c_initial=csat, media_height=media_height)
            cs.append(c)
            hs.append(h)
        plt.plot(hs, cs, label="t={} mins".format(t_min))
        plt.legend()
        #plt.ylim(0, 200)
    plt.title("Gradient over Time (OCR={} fmols/mm2/s  media volume={}uL".format(ocr, vol))

def plot_time_to_diffusion_limit_vs_ocr(ocrs, vols=[100, 150, 200], csat=185):
    for vol in vols:
        times = [kinetics.time_to_constrained(ocr, vol) for ocr in ocrs]
        t_hrs = [t/3600 if t is not None else None for t in times]
        plt.plot(ocrs, t_hrs, label="{} uL".format(vol))
    plt.xlabel("OCR (fmols/mm2/s)")
    plt.ylabel("Time to Diffusion-Limit (hrs)")
    plt.legend()
    plt.title("Time to OCR Diffusion Limit")

def plot_time_to_diffusion_limit_vs_vol(vols, ocrs=[100, 200, 300, 400], csat=185):
    for ocr in ocrs :
        times = [kinetics.time_to_constrained(ocr, vol) for vol in vols]
        t_hrs = [t/3600 if t is not None else None for t in times]
        plt.plot(vols, t_hrs, label="OCR={}".format(ocr))
    plt.xlabel("Volume (uL)")
    plt.ylabel("Time to Diffusion-Limit (hrs)")
    plt.legend()
    plt.title("Time to OCR Diffusion Limit")


@cmds.command()
@click.argument('vols', nargs=-1, type=int)
@click.option('--csat', default=185, type=int)
def time_to_diffusion_limit_vs_ocr(vols, csat=185):
    if len(vols) == 0:
        vols = [100, 150, 200]
        print("Defulating to volumes {}".format(vols))
    ocrs = list(range(50, 400, 10))
    print("Calulating...")
    plot_time_to_diffusion_limit_vs_ocr(ocrs=ocrs, vols=vols, csat=csat)
    plt.show()

@cmds.command()
@click.argument('ocrs', nargs=-1, type=int)
@click.option('--csat', default=185, type=int)
def time_to_diffusion_limit_vs_vol(ocrs, csat=185):
    if len(ocrs) == 0:
        ocrs = [100, 200, 300, 400]
        print("Defulating to ocrs {}".format(ocrs))

    vols = list(range(50, 300, 10))
    print("Calulating...")
    plot_time_to_diffusion_limit_vs_vol(ocrs=ocrs, vols=vols, csat=csat)
    plt.show()

@cmds.command()
@click.argument('ocr', default=100)
@click.option('--vol', default=100, type=int)
@click.option('--csat', default=185, type=int)
def o2_gradient_over_time(ocr=100, vol=100, csat=185, tmax_hrs=1, delta_t_mins=15):
    plot_gradient_evolution(ocr=ocr, vol=vol, csat=csat, tmax_hrs=tmax_hrs, delta_t_mins=delta_t_mins)
    plt.show()

if __name__ == '__main__':
    cmds()
    ocr = 180
    for vol in [60, 75, 100, 150, 200]:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600*6, 30))
        cs = [concentration(0, t, q, c_initial=185, media_height=media_vol_to_height(vol)) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="volume={} uL".format(vol))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, 200)
    
    plt.legend()
    plt.title("O2 at bottom by media volume with OCR={} fmols/mm2/s".format(ocr))
    plt.show()
    
    for vol in [60, 75, 100, 150, 200]:
        q = flux_units_convert(ocr)
        ts = list(range(0, 3600*6, 30))
        cs = [concentration(1.25, t, q, c_initial=185, media_height=media_vol_to_height(vol)) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="volume={} uL".format(vol))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, 200)
    
    plt.legend()
    plt.title("O2 at 1.25mm by media volume with OCR={} fmols/mm2/s".format(ocr))
    plt.show()

    def plot_gradient_evolution(ocr=100, media_height=3.1, csat=185, tmax_hrs=1, delta_t_mins=15):
        nz = 20

        for t in range(0, 3600*tmax_hrs, delta_t_mins*60):#[3600*hr for hr in range(5)]:
            cs = []
            hs = []
            t_min = t//60
            for z_i in range(nz):
                h = z_i * media_height / nz
                #q = flux_units_convert(ocr)
                #print(q)
                q = ocr * 1e-3  #fmtomols/mm3 to micromolar
                #ts = list(range(0, 3600*24, 10))
                c = concentration(h, t, q, c_initial=csat, media_height=media_height)
                cs.append(c)
                hs.append(h)
            plt.plot(hs, cs, label="t={} mins".format(t_min))
            plt.legend()
            #plt.ylim(0, 200)
            plt.title("Gradient over Time (OCR={} fmols/mm2/s) media height {}mm".format(ocr, media_height))

    plot_gradient_evolution()
    plt.show()

    ocrs = list(range(50, 400, 10))
    for vol in [100, 150, 200]:
        times = [time_to_constrained(ocr, vol) for ocr in ocrs]
        t_hrs = [t/3600 if t is not None else None for t in times]
        plt.plot(ocrs, t_hrs, label="{} uL".format(vol))
    plt.xlabel("OCR (fmols/mm2/s)")
    plt.ylabel("Time to Diffusion-Limit (hrs)")
    plt.legend()
    plt.title("Time to OCR Diffusion Limit")
