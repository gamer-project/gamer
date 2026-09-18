import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

# load the command-line parameters
parser = argparse.ArgumentParser(description="Dust dust_dens evolution plotter")
args = parser.parse_args()

# physical constants
CONST_CM  = 1.0
CONST_AMU = 1.660539040e-24
M_H       = 1.6735575e-24
K_B       = 1.380649e-16
AMU_IN_G  = CONST_AMU
MYR_IN_S  = 3.15569252e13
GYR_IN_S  = MYR_IN_S * 1.0e3

PREFIX     = '../'

# Constants
A_UM  = 0.1
OMEGA = 2.5
MU    = 0.6

# Load all simulation-wide constant parameters from the first output file,
# to ensure consistency with the actual simulation setup
first_file = os.path.join(PREFIX, 'Data_000000')
with h5py.File(first_file, "r") as f:
    DUST_TO_GAS    = float(f["Info"]["InputTest"]["GrackleTest_DustToGasRatio"])
    K_MYR          = float(f["Info"]["InputTest"]["GrackleTest_ExpCoolCoeff"])
    UNIT_L         = float(f['Info']['InputPara']['Unit_L'])
    UNIT_T         = float(f['Info']['InputPara']['Unit_T'])
    UNIT_D         = float(f['Info']['InputPara']['Unit_D'])
    GAMMA          = float(f['Info']['InputPara']['Gamma'])
    T_over_mmw_min = float(f["Info"]["InputTest"]["GrackleTest_TempOverMMW_Min"])

print("Dust-to-gas ratio =", DUST_TO_GAS)
print("k [Myr^-1] =", K_MYR)
if K_MYR != 0:
    t_cool_myr = 1.0 / K_MYR
    print("t_cool [Myr] =", t_cool_myr)
else:
    t_cool_myr = 1.0

k = K_MYR / MYR_IN_S

# initial conditions
# Choose ONE density setting:
# Case A: amu case
gas_rho0_amu = 1.0
gas_rho0_cgs = gas_rho0_amu * AMU_IN_G
dust_rho0    = DUST_TO_GAS * gas_rho0_amu
# Case B: g/cm^3 case
# gas_rho0_cgs = 1.0e-28
# gas_rho0_amu = gas_rho0_cgs / AMU_IN_G
# dust_rho0    = DUST_TO_GAS * gas_rho0_amu

T0 = T_over_mmw_min * MU

# Load data
table     = np.loadtxt(os.path.join(PREFIX, "Record__Conservation"))
time      = table[:, 0]
dust_dens = table[:,47] * UNIT_D / (CONST_AMU / CONST_CM**3)


# units
UNIT_V = UNIT_L / UNIT_T
UNIT_E = UNIT_V ** 2
e_code = T0 / (UNIT_E * ((GAMMA-1) * MU * M_H / K_B))
e_phys = e_code * UNIT_E

# Configuration: output filenames
fileout  = "fig__DustDensity_plot"
fig_name = "Dust Density v.s Time"

# Functions
def internal_energy(t):
    return e_phys * np.exp(-k*t)

def tsp_e(e_t):
    gas_rho_cgs = gas_rho0_cgs
    const_1 = 0.17 * (A_UM / 0.1) * (1.0e-27 / gas_rho_cgs) * GYR_IN_S
    const_2 = ((10**6.3 * K_B) / ((GAMMA-1)*MU*M_H)) ** OMEGA
    tsp = const_1 * (const_2 / e_t**OMEGA + 1.0)
    return tsp

def drho_dt(t, dust_rho):
    e_t = internal_energy(t)
    tsp = tsp_e(e_t)
    return -3.0 / tsp * dust_rho


# normalized
dust_dens_norm = dust_dens / dust_dens[0]

# time normalized by cooling time
time_cool = time / t_cool_myr

# Plot
f, ax = plt.subplots(1, 1)
f.subplots_adjust(wspace=0.4)
ax.set_xlabel(r"Number of cooling times $(t/t_{\rm cool})$", fontsize="large")
ax.set_title(fig_name)
ax.plot(time_cool, dust_dens_norm, 'ro', lw=1, mec='none', ms=5.0, label='Numerical')

# Reference solution
if K_MYR == 0:
    tsp     = tsp_e(e_phys)
    tsp_myr = tsp / MYR_IN_S
    rho_ref = dust_rho0 * np.exp((-3/tsp_myr) * time)
    rho_ref_norm = rho_ref / rho_ref[0]
    ax.plot(time_cool, rho_ref_norm, 'b-', lw=1.5, label="Reference")
else:
    t_span = (0, time[-1]*MYR_IN_S)
    t_eval = time * MYR_IN_S
    sol = solve_ivp(drho_dt, t_span, [dust_rho0], t_eval=t_eval, rtol=1e-10, atol=1e-14)
    rho_ref = sol.y[0]
    rho_ref_norm = rho_ref / rho_ref[0]
    ax.plot((sol.t / MYR_IN_S) / t_cool_myr, rho_ref_norm, 'b-', label="Reference")

# Final point comparison
t_num_final = time[-1]
rho_num_final = dust_dens_norm[-1]
rho_ref_final = rho_ref_norm[-1]
rho_abs_err = abs(rho_num_final - rho_ref_final)
rho_rel_err = rho_abs_err / abs(rho_ref_final)
print("====================================")
print("Final point comparison")
print("t_final =", t_num_final, "Myr")
print("t_final/t_cool =", time_cool[-1])
print("rho_num_final_norm =", rho_num_final)
print("rho_ref_final_norm =", rho_ref_final)
print("final absolute normalized error =", rho_abs_err)
print("final relative error =", rho_rel_err * 100.0, "%")
print("====================================")

# Text on the figure
textstr = (
    rf"$\mathrm{{Final\ relative\ error}} = {rho_rel_err*100.0:.3f}\%$" "\n"
    rf"$\mathrm{{Final\ abs.\ norm.\ error}} = {rho_abs_err:.3e}$"
)
ax.text(
    0.97, 0.78, textstr,
    transform=ax.transAxes,
    fontsize=14,
    fontfamily="serif",
    math_fontfamily="stix",
    verticalalignment="top",
    horizontalalignment="right"
)

# Finalize figure
ax.set_yscale('linear')
ax.set_xlim(0, time_cool[-1]*1.05)
ax.set_ylabel(r'$\rho_{\rm dust}/\rho_{\rm dust,0}$', fontsize='large')
ax.legend()
plt.savefig(fileout + ".png", bbox_inches='tight', pad_inches=0.05, dpi=150)
# plt.show()
