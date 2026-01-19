
def fmols_per_mm2_per_s_to_mols_per_liter_per_hr(rate_fmols_per_mm2_per_s):
    # a liter is 1e6 cubic mm and an hour is 3600 seconds
    # a fmols is 1e-15 mols
    # since we are using a 1D cyclindrical model, 1 m^3 = 1 m^2 (cross
    # sectional area) and 1 m^2 = 1e6 mm^2

    return 1e-15*rate_fmols_per_mm2_per_s * 3600 * 1e6

def mols_per_liter_per_hr_to_fmols_per_mm2(rate_mols_per_L_per_hour):
    #convert reaction rate in mols/L/hour to fmols/mm2/s
    # NOTE: a liter is 1e6 mm3 in 3D and 1e6 mm2 in 2D (cross sectional
    # area units for 1D cylindrical model).
    # 1 mol is 1e15 fmols
    rate_fmols_per_mm2_per_sec  = 1e15 * rate_mols_per_L_per_hour / 1e6 / 3600
    return rate_fmols_per_mm2_per_sec

def fmols_per_mm3_to_micromolar(fmols_per_mm3):
    # 1 fmols per mm^3 equals 1 nanomolar
    return 1e3*fmols_per_mm3

def rate_umolar_per_second_to_molar_per_hour(rate_umolar_per_s):
    return 3600*rate_umolar_per_s / 1e6

def rate_mols_per_L_per_hour_umolar_per_s(rate):
    return rate * 1e6 / 3600

def rate_pmols_per_L_per_minute_to_umolar_per_s(rate):
    return rate * 60 / 1e3

"1 uM = 1000 fmol mm^−3"
# 1 L in 1e6 mm^3
