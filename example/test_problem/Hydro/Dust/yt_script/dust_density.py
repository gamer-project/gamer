import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

# load the command-line parameters
parser = argparse.ArgumentParser(description="Dust dust_dens evolution plotter")
parser.add_argument("-option", type=str, required=True, choices=["edot_0", "edot_const"], help="Analytic option")
args = parser.parse_args()

# Constants
Const_cm  = 1.0
Const_amu = 1.660539040e-24
UNIT_L    = 3.08567758149e21
UNIT_M    = 1.9885e42
UNIT_D    = UNIT_M / UNIT_L**3
UNIT_T    = 3.15569252e13
a_um      = 0.1
omega     = 2.5

# physical constants
gamma    = 5.0/3.0
mu       = 0.62
m_H      = 1.6735575e-24
k_B      = 1.380649e-16
amu_in_g = 1.66054e-24
Gyr_in_s = 3.154e16
Myr_in_s = Gyr_in_s / 1e3

# initial conditions
dust_rho0 = 0.01
gas_rho0  = 0.1
T0        = 1031251.393124
k         = 0.4 / Myr_in_s

# units
UNIT_V = UNIT_L / UNIT_T
UNIT_E = UNIT_V ** 2
e_code = 2.157133263873971e-02  # code unit
e_phys = e_code * UNIT_E


# Configuration: output filenames
fileout  = "fig__DustDensity_plot"
fig_name = "Dust Density v.s Time"
prefix   = '../'

# Functions
def internal_energy(e_0, k, t):
    return e_0 * np.exp(-k*t)

def tsp_e(e_t):
    gas_rho_cgs = gas_rho0 * amu_in_g
    const_1 = 0.17 * (a_um / 0.1) * (1.0e-27 / gas_rho_cgs) * Gyr_in_s
    const_2 = ( (10**6.3 * k_B) / ((gamma-1)*mu*m_H) ) ** omega
    tsp = const_1 * (const_2 / e_t**omega + 1.0)
    return tsp

def drho_dt(t, dust_rho):
    e_t = internal_energy(e_phys, k, t)
    tsp = tsp_e(e_t)
    return -3.0 / tsp * dust_rho

# Load data
f = h5py.File('../Data_%06d'%(0), 'r')
UNIT_D  = f['Info']['InputPara']['Unit_D']
BoxSize = f['Info']['InputPara']['BoxSize']
table = np.loadtxt("../Record__Conservation")
time      = table[:, 0]
dust_dens = table[:,47] * UNIT_D / (Const_amu / Const_cm**3) / BoxSize**3

# sorting 
sort_idx = np.argsort(time)
time = time[sort_idx]
dust_dens = dust_dens[sort_idx]

# remove duplicates
time, uniq_idx = np.unique(time, return_index=True)
dust_dens = dust_dens[uniq_idx]

# Plot
f, ax = plt.subplots(1, 1)
f.subplots_adjust(wspace=0.4)

ax.set_xlabel("$\\mathrm{t\\ [Myr]}$", fontsize="large")
ax.set_title(fig_name)
ax.plot(time, dust_dens, 'ro', lw=1, mec='none', ms=5.0, label='Numerical')

# Analytic solution
rho_analytic = None
if args.option == "edot_0":
    gas_rho_cgs = gas_rho0 * amu_in_g
    const_1 = 0.17 * (a_um / 0.1) * (1e-27 / gas_rho_cgs) * Gyr_in_s
    const_2 = (10**6.3 / T0)**omega
    tsp     = const_1 * (const_2 + 1.0)
    tsp_myr = tsp / Myr_in_s
    rho_analytic = dust_rho0 * np.exp((-3/tsp_myr) * time)
    ax.plot(time, rho_analytic, 'b-', lw=1.5, label="Analytic")

elif args.option == "edot_const":
    t_span = (0, time[-1]*Myr_in_s)
    t_eval = time * Myr_in_s
    sol = solve_ivp(drho_dt, t_span, [dust_rho0], t_eval=t_eval)
    rho_analytic = sol.y[0]
    ax.plot(sol.t / Myr_in_s, sol.y[0], 'b-', label="Analytic")

# Finalize figure
ax.set_yscale('log')
ax.set_xlim(1.0e-1, 20)
ax.yaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]))
ax.set_ylabel('$\\mathrm{Dust\ dust_dens\ [amu/cm^3]}$', fontsize='large')
ax.legend()

plt.savefig(fileout + ".png", bbox_inches='tight', pad_inches=0.05, dpi=150)
plt.show()
