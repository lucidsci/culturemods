"""
#from "Kinetics of Gas Diffusion in Mammalian Cell Culture Systems" - McLimans, Blumension, Tunnah
#equations 23
"""
import math

#from "Kinetics of Gas Diffusion in Mammalian Cell Culture Systems" - McLimans, Blumension, Tunnah. 1968. 
MONOLAYER_THICKNESS_MM=0
DEFAULT_HEIGHT_MM=3.1
#o2 diffusion coefficient
O2_DIFFUSION_37C = 3e-3 #in mm2/s

def H(n, a , s=MONOLAYER_THICKNESS_MM):
    #calculate  H_n from eq 19 where a is fluid height and s is monolayer thickness
    return math.pi * (2*n + 1 ) / (2 * (a - s))

def summation(x, t, a=DEFAULT_HEIGHT_MM, s=MONOLAYER_THICKNESS_MM, iterations=100, D=O2_DIFFUSION_37C):
    #calculate summand from eq 23 for the given n
    sigma = 0 
    for n in range(iterations):
        H_n = H(n, a, s)
        sigma += math.pow(H_n, -2) * math.cos(H_n * x) * math.exp(-1 * math.pow(H_n, 2) * D * t)
    
    return sigma



def concentration(x, t, Q=5, c_initial = 200, media_height=DEFAULT_HEIGHT_MM):
    """
    Calculate concentration over time at given height/depth based on Eq. 23 from 
    "Kinetics of Gas Diffusion in Mammalian Cell Culture Sytems".  McLimans, et. al.  1968.

    Q is flux in umols/mm2/sec
    c_initial is o2 concentration in micromolar
    media_height in centimeters
    """
    a = media_height #reference uses 'a' for fluid height
    s = MONOLAYER_THICKNESS_MM
    D = O2_DIFFUSION_37C
    
    #reference uses 'b' for c_initial - the intitial o2 concentration in fluid
    c_x_t = c_initial - (a-s-x)*(Q/D) + (2*Q / ((a - s)*D)) * summation(x, t, s=s, a=a, D=D)
    return c_x_t

def flux_units_convert(flux_fmols_mm_sq_per_sec):
    """
    convert flux from 'normal' of fmols/mm2/sec to
    units the model expects of umols/mm2/sec
    """
    q = flux_fmols_mm_sq_per_sec * 1e-3
    return q
def media_vol_to_height(vol_uL, well_radius_mm=3.2):
    return vol_uL / (math.pi * well_radius_mm**2)

def time_to_constrained(ocr, media_vol_uL=100, c_sat=185):
    h = media_vol_to_height(media_vol_uL)
    tmax_hrs = 24
    q = ocr*1e-3
    for t_s in range(8*3600):
        c = concentration(0, t_s, q, c_initial=c_sat, media_height=h)
        if c <= 0:
            return t_s

    return None

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    for ocr in [50, 100, 150, 200, 300]:
        q = flux_units_convert(ocr)
        ts = list(range(0, 3600*4, 1))
        cs = [concentration(0, t, q, c_initial=185, media_height=3.1) for t in ts]
        plt.plot([t/3600 for t in ts], cs, label="OCR={}".format(ocr))
        plt.xlabel("Time (hours)")
        plt.ylabel("O2 (micromolar)")
        plt.ylim(0, 200)
    
    plt.legend()
    plt.title("O2 at bottom by OCR in 100 uL")
    plt.show()
    
    ocr = 180
    for vol in [60, 75, 100, 150, 200]:
        q = flux_units_convert(ocr)
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

    def plot_gradient_evolution(ocr=100, media_height=3.1, c_sat=185, tmax_hrs=1, delta_t_mins=15):
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
                c = concentration(h, t, q, c_initial=c_sat, media_height=media_height)
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
