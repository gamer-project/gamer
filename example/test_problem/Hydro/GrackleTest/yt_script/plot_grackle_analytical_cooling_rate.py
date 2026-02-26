import h5py
import numpy as np
import matplotlib
# Use 'Agg' backend to avoid PuTTY X11 errors
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import RegularGridInterpolator
from mpl_toolkits.axes_grid1 import make_axes_locatable


# -------------------------------------------------------------------------
# user-specified parameters
grackle_path          = "../CloudyData_UVB=HM2012.h5"
z_now                 = 3.0              # Redshift
output_mode           = int(2)           # [0] lists all groups and datasets of "CloudyData_UVB=HM2012.h5"
                                         # [1] computes cooling and heating rates for specified (z, nH, T) using interpolation 
                                         # [2] output cooling and heating rates for a grid of (z, nH, T) values and saves to a new HDF5 file
nH_now                = 1e-3             # Hydrogen density (cm^-3); for option [1]
T_now                 = 1e4              # Temperature (K); for option [1]
mp                    = 1.6726e-24       # Proton mass (g); for option [2]
X_h                   = 0.716            # Hydrogen mass fraction for HM2012; for option [2]
# Axis Ranges
rho_range             = np.logspace(-29, -21, 100)  # g/cm^3; for option [2]
t_mu_range            = np.logspace(0, 8, 100)      # K (T/mu); for option [2]
font_size             = 14               # Font size for plots; for option [2]


# -------------------------------------------------------------------------
# Define the rate function in the GLOBAL scope so all modes can use it
def get_grackle_rates(file_path, z_target, nh_target, t_target):
    with h5py.File(file_path, 'r') as f:
        ds_cool = f['CoolingRates/Primordial/Cooling']
        ds_heat = f['CoolingRates/Primordial/Heating']

        # Extract Grids
        nh_grid = ds_cool.attrs['Parameter1'][:]
        z_grid  = ds_cool.attrs['Parameter2'][:]
        t_grid  = ds_cool.attrs['Temperature'][:]

        cool_data = ds_cool[:]
        heat_data = ds_heat[:]

        # Re-order and Clean Grids
        grids = [nh_grid, z_grid, t_grid]
        for i in range(len(grids)):
            sort_idx = np.argsort(grids[i])
            grids[i] = grids[i][sort_idx]
            cool_data = np.take(cool_data, sort_idx, axis=i)
            heat_data = np.take(heat_data, sort_idx, axis=i)

            mask = np.concatenate(([True], np.diff(grids[i]) > 0))
            grids[i] = grids[i][mask]
            cool_data = np.compress(mask, cool_data, axis=i)
            heat_data = np.compress(mask, heat_data, axis=i)

        cool_interp = RegularGridInterpolator(tuple(grids), cool_data)
        heat_interp = RegularGridInterpolator(tuple(grids), heat_data)

        # File bounds for clipping
        # log_nh: Parameter1, Redshift: Parameter2, Temp: Temperature grid
        log_nh = np.clip(np.log10(nh_target), grids[0].min(), grids[0].max())
        z_safe = np.clip(z_target, grids[1].min(), grids[1].max())
        t_safe = np.clip(t_target, grids[2].min(), grids[2].max())

        points = np.dstack((log_nh, z_safe, t_safe))
        return cool_interp(points), heat_interp(points)
def calculate_cooling_time(nH, T, lambda_net_norm):
    # Constants (CGS)
    kB = 1.3806e-16
    yr_to_sec = 3.154e7
    
    # Assumptions for fully ionized primordial (metal-free) gas at 10^4K
    # Mean molecular weight mu ~ 0.6
    mu = 0.6
    # Total number density of all species (H, He, e-)
    n_tot = nH / (mu * X_h) 
    
    # Internal Energy (erg/cm^3)
    energy_density = 1.5 * n_tot * kB * T
    
    # Volumetric Cooling Rate (erg/cm^3/s); lambda_net_norm is |Heating - Cooling|
    vol_cooling_rate = lambda_net_norm * (nH**2)
    
    if vol_cooling_rate == 0:
        return np.inf
        
    t_sec = energy_density / vol_cooling_rate
    return t_sec / (yr_to_sec * 1e6) # Convert to Myr


# -------------------------------------------------------------------------
if (output_mode == 0):
    with h5py.File(grackle_path, 'r') as f:
        f.visit(print)


# -------------------------------------------------------------------------
if (output_mode == 1):
    try:
        cooling_rate, heating_rate = get_grackle_rates(grackle_path, z_target=z_now, nh_target=nH_now, t_target=T_now)
        print("-" * 45)
        print(f"Metal-Free Gas Rates (z={z_now}, nH={nH_now:.2e}, T={T_now:.2e}K):")
        print(f"  Cooling Rate (Λ/n_H²): {cooling_rate[0,0]:.2e} erg*cm^3/s")
        print(f"  Heating Rate (Γ/n_H²): {heating_rate[0,0]:.2e} erg*cm^3/s")
        print(f"  Net: {'Heating' if heating_rate[0,0] > cooling_rate[0,0] else 'Cooling'}")
        print("-" * 45)
    except Exception as e:
        print(f"Error: {e}")
    net_rate = abs(heating_rate[0,0] - cooling_rate[0,0])
    t_cool_myr = calculate_cooling_time(nH_now, T_now, net_rate)
    print(f"Cooling Timescale: {t_cool_myr:.2e} Myr")


# -------------------------------------------------------------------------
if (output_mode == 2):
    RHO, TMU = np.meshgrid(rho_range, t_mu_range)
    net_rates = np.zeros_like(RHO)

    print("Calculating heating/cooling rates on grid")
    nh_val     = RHO / (mp / X_h)
    t_val      = TMU
    z_val      = z_now*np.ones_like(nh_val)
    cool, heat = get_grackle_rates(grackle_path, z_val, nh_val, t_val)
    net_rates  = (heat - cool) * (nh_val**2)

    print("Start plotting")
    # --- Refined Plotting ---
    fig, ax = plt.subplots(figsize=(10, 8)) # Slightly wider to accommodate colorbar

    # 1. Force the plot area to be a perfect square
    ax.set_box_aspect(1)

    # Normalization
    norm = colors.SymLogNorm(linthresh=1e-29, vmin=-1e-24, vmax=1e-24, base=10)

    # Plotting the mesh
    im = ax.pcolormesh(RHO, TMU, net_rates, norm=norm, cmap='bwr', shading='auto')

    # 2. Add black dashed line for net rate of zero
    ax.contour(RHO, TMU, net_rates, levels=[0], colors='black', 
            linestyles='dashed', linewidths=1.5)

    # Formatting Axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Density $\rho_{\rm gas}$ [g/cm$^3$]', fontsize=font_size)
    ax.set_ylabel(r'$T/\mu$ [K]', fontsize=font_size)

    # 3. Longer and thicker ticks pointing inward
    ax.tick_params(axis='both', which='major', direction='in', length=10, width=2, 
                top=True, right=True, labelsize=font_size)
    ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1.5, 
                top=True, right=True)

    # 4. Updated Colorbar Logic (Fixed Height)
    # 'shrink=0.82' matches the height of a square set_box_aspect(1) frame
    cbar = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.03) 
    cbar.set_label(r'Grackle Cooling Rate $\left( \frac{\text{erg}}{\text{cm}^3 \cdot \text{s}} \right)$', 
                fontsize=font_size, labelpad=15)

    # Formatting colorbar ticks
    tick_locs = [-1e-24, -1e-25, -1e-26, -1e-27, -1e-28, -1e-29, 0, 1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24]
    tick_labels = [r'$-10^{-24}$', r'$-10^{-25}$', r'$-10^{-26}$', r'$-10^{-27}$', r'$-10^{-28}$', r'$-10^{-29}$', 
                '0', 
                r'$10^{-29}$', r'$10^{-28}$', r'$10^{-27}$', r'$10^{-26}$', r'$10^{-25}$', r'$10^{-24}$']

    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(tick_labels, fontsize=font_size)
    cbar.ax.tick_params(direction='in', length=5, width=1.5) # Colorbar ticks inward too

    # 5. Thicker frame (spines)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
    for spine in cbar.ax.spines.values():
        spine.set_linewidth(2.0)

    # 6. Text annotations
    ax.text(0.03, 0.97, 
            f'max = {net_rates.max():.3e} erg/(cm$^3 \cdot$ s)\nmin = {net_rates.min():.3e} erg/(cm$^3 \cdot$ s)', 
            transform=ax.transAxes, ha='left', va='top', 
            fontsize=font_size*0.9, color='black', fontweight='bold')

    ax.text(0.97, 0.97, rf'$\mathbf{{z = {z_now:.2f}}}$', 
            transform=ax.transAxes, ha='right', va='top', 
            fontsize=font_size, color='black', fontweight='bold')

    plt.savefig(f"grackle_phase_plot_z={z_now}.png", dpi=300, bbox_inches='tight')
    print("Finished plotting")
