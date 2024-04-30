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
@click.argument('ocrs', nargs=-1)
@click.option('--vol', default=100)
@click.option('--csat', default=185)
def o2_at_bottom_by_ocr(ocrs, vol=100, csat=185):
    height = media_vol_to_height(vol)
    for ocr in [50, 100, 150, 200, 300]:
        q = kinetics.flux_units_convert(ocr)
        ts = list(range(0, 3600*4, 1))
        cs = [kinetics.concentration(0, t, q, c_initial=csat, media_height=height) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="OCR={}".format(ocr))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, 200)
    
    plt.legend()
    plt.title("O2 at bottom by OCR in {} uL".format(vol))
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
