import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use 'Agg' backend to avoid PuTTY X11 errors
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import RegularGridInterpolator


# -------------------------------------------------------------------------
# Set user-specified parameters
GRACKLE_PATH = "../CloudyData_UVB=HM2012.h5"
Z_NOW        = 0.0                         # Redshift
NH_NOW       = 1.0e-3                      # Hydrogen density (cm^-3)
T_NOW        = 1.0e4                       # Temperature (K)
X_H          = 0.716                       # Hydrogen mass fraction for HM2012
RHO_RANGE    = np.logspace(-29, -21, 100)  # Denisty axis range (g/cm^3)
T_MU_RANGE   = np.logspace(  0,   8, 100)  # Temperature/mu axis range (K)


# -------------------------------------------------------------------------
# Interpolate the cooling and heating rates from the grackle table
def get_grackle_rates_and_mmw(file_path, z_target, nh_target, t_target):
    with h5py.File(file_path, 'r') as f:
        # Load Cooling, Heating, and Mean Molecular Weight directly from "HM2012.h5"
        ds_cool = f['CoolingRates/Primordial/Cooling']
        ds_heat = f['CoolingRates/Primordial/Heating']
        ds_mmw  = f['CoolingRates/Primordial/MMW']

        # Extract grids (shared across these datasets)
        nh_grid   = ds_cool.attrs['Parameter1'][:]
        z_grid    = ds_cool.attrs['Parameter2'][:]
        t_grid    = ds_cool.attrs['Temperature'][:]
        
        cool_data = ds_cool[:]
        heat_data = ds_heat[:]
        mmw_data  = ds_mmw[:]

        grids = [nh_grid, z_grid, t_grid]
        data_list = [cool_data, heat_data, mmw_data]

        # Sort and clean grids and data simultaneously
        for i in range(len(grids)):
            sort_idx = np.argsort(grids[i])
            grids[i] = grids[i][sort_idx]
            for j in range(len(data_list)):
                data_list[j] = np.take(data_list[j], sort_idx, axis=i)

            mask = np.concatenate(([True], np.diff(grids[i]) > 0))
            grids[i] = grids[i][mask]
            for j in range(len(data_list)):
                data_list[j] = np.compress(mask, data_list[j], axis=i)

        # Create interpolators for all three
        cool_interp = RegularGridInterpolator(tuple(grids), data_list[0])
        heat_interp = RegularGridInterpolator(tuple(grids), data_list[1])
        mmw_interp  = RegularGridInterpolator(tuple(grids), data_list[2])

        # Bounds clipping
        log_nh = np.clip(np.log10(nh_target), grids[0].min(), grids[0].max())
        z_safe = np.clip(z_target, grids[1].min(), grids[1].max())
        t_safe = np.clip(t_target, grids[2].min(), grids[2].max())
        points = np.dstack((log_nh, z_safe, t_safe))

        return cool_interp(points), heat_interp(points), mmw_interp(points)


# -------------------------------------------------------------------------
# Compute the cooling time
def calculate_cooling_time(nH, T, lambda_net_norm, mu_val):
    # Constants (CGS)
    kB         = 1.3806e-16
    yr_to_sec  = 3.154e7
    Myr_to_sec = 1.0e6 * yr_to_sec

    # Total number density of all species (H, He, e-)
    n_tot = nH / (mu_val * X_H)

    # Internal Energy (erg/cm^3)
    energy_density = 1.5 * n_tot * kB * T

    # Volumetric Cooling Rate (erg/cm^3/s); lambda_net_norm is |Heating - Cooling|
    vol_cooling_rate = lambda_net_norm * (nH**2)

    # Cooling time
    t_sec = energy_density / vol_cooling_rate # in sec
    t_Myr = t_sec / Myr_to_sec                # in Myr

    return t_Myr


# -------------------------------------------------------------------------
# List all groups and datasets of "CloudyData_UVB=HM2012.h5"
def inspect_grackle_data_file(grackle_path=GRACKLE_PATH):
    print("-" * 45)
    print("File structure of %s"%grackle_path)
    print("-" * 45)
    with h5py.File(grackle_path, 'r') as f:
        f.visit(print)
    print("-" * 45)
    print("")


# -------------------------------------------------------------------------
# Compute cooling and heating rates for specified (z, nH, T) using interpolation
def output_one_point():
    print("-" * 55)
    print(f"Metal-Free Gas Rates (z={Z_NOW}, nH={NH_NOW:.2e}, T={T_NOW:.2e}K)")
    print_one_point(*calculate_one_point())



def calculate_one_point(grackle_path=GRACKLE_PATH, z_now=Z_NOW, nH_now=NH_NOW, T_now=T_NOW):
    # Retrieve all 3 values from interpolations directly from grackle tables
    cool_rate, heat_rate, mu_val = get_grackle_rates(grackle_path, z_now, nH_now, T_now)
    # Extract scalar values from the interpolation arrays
    c_rate = cool_rate[0,0]
    h_rate = heat_rate[0,0]
    m_val  = mu_val[0,0]
    # Use the MU returned by Grackle instead of the hardcoded 0.6 global
    net_rate = abs(h_rate - c_rate)
    # Calculate cooling time using the table's MU
    t_cool_Myr = calculate_cooling_time(nH_now, T_now, net_rate, mu_val=m_val)
    return c_rate, h_rate, net_rate, t_cool_Myr, m_val


def calculate_one_point(grackle_path=GRACKLE_PATH, z_now=Z_NOW, nH_now=NH_NOW, T_now=T_NOW):
    cooling_rate, heating_rate, mmw = get_grackle_rates_and_mmw(grackle_path, z_target=z_now, nh_target=nH_now, t_target=T_now)
    # Use the MU returned by Grackle instead of the hardcoded 0.6 global
    net_rate   = abs(heating_rate[0,0] - cooling_rate[0,0])
    # Pass mmw[0,0] to the cooling time calculation for 100% consistency
    t_cool_Myr = calculate_cooling_time(nH_now, T_now, net_rate, mu_val=mmw[0,0])
    return cooling_rate, heating_rate, net_rate, t_cool_Myr


def print_one_point(cooling_rate, heating_rate, net_rate, t_cool_Myr):
    print("-" * 55)
    print(f"  Cooling Rate (Λ/n_H²): {cooling_rate[0,0]:.2e} erg*cm^3/s")
    print(f"  Heating Rate (Γ/n_H²): {heating_rate[0,0]:.2e} erg*cm^3/s")
    print(f"  Net                  : {'Heating' if heating_rate[0,0] > cooling_rate[0,0] else 'Cooling'}")
    print(f"  Timescale            : {t_cool_Myr:.2e} Myr")
    print("-" * 55)
    print("")


# -------------------------------------------------------------------------
# Output cooling and heating rates for a grid of (z, nH, T) values to a png file
def output_grid():
    plot_grid(Z_NOW, *calculate_grid())


def calculate_grid(grackle_path=GRACKLE_PATH, z_now=Z_NOW, rho_range=RHO_RANGE, t_mu_range=T_MU_RANGE):
    # Create the grid
    grid_rho, grid_t_mu = np.meshgrid(rho_range, t_mu_range)

    # Calculate heating/cooling rates on grid 
    mp     = 1.6726e-24           # Proton mass (g)
    nh_val = grid_rho / (mp / X_H)
    z_val  = z_now * np.ones_like(nh_val)

    # (1) Use a reasonable guess for MU to get an initial Temperature
    mu_guess = 0.6 
    t_guess  = grid_t_mu * mu_guess
    
    # (2) Get actual MMW from the table for these conditions
    _, _, mu_actual = get_grackle_rates_and_mmw(grackle_path, z_val, nh_val, t_guess)
    
    # (3) Use the self-consistent MU to find true T and Rates
    t_actual = grid_t_mu * mu_actual
    cool, heat, _ = get_grackle_rates_and_mmw(grackle_path, z_val, nh_val, t_actual)
    
    net_rates = (heat - cool) * (nh_val**2)
    return grid_rho, grid_t_mu, net_rates


def plot_grid(z_now, grid_rho, grid_t_mu, net_rates):
    # Create the figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Set the font size for plots
    font_size  = 14

    # Force the plot area to be a perfect square
    ax.set_box_aspect(1)

    # Normalize
    norm = colors.SymLogNorm(linthresh=1e-29, vmin=-1e-24, vmax=1e-24, base=10)

    # Plot the mesh
    im = ax.pcolormesh(grid_rho, grid_t_mu, net_rates, norm=norm, cmap='bwr', shading='auto')

    # Add black dashed line for net rate of zero
    ax.contour(grid_rho, grid_t_mu, net_rates, levels=[0], colors='black',
               linestyles='dashed', linewidths=1.5)

    # Format axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Density $\rho_{\rm gas}$ [g/cm$^3$]', fontsize=font_size)
    ax.set_ylabel(r'$T/\mu$ [K]', fontsize=font_size)

    # Add longer and thicker ticks pointing inward
    ax.tick_params(axis='both', which='major', direction='in', length=10, width=2,
                   top=True, right=True, labelsize=font_size)
    ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1.5,
                   top=True, right=True)

    # Update colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=1.0, pad=0.03)
    cbar.set_label(r'Grackle Cooling Rate $\left( \frac{\mathrm{erg}}{\mathrm{cm}^3 \cdot \mathrm{s}} \right)$',
                   fontsize=font_size, labelpad=15)

    # Format colorbar ticks
    tick_locs   = [-1e-24, -1e-25, -1e-26, -1e-27, -1e-28, -1e-29, 0, 1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24]
    tick_labels = [r'$-10^{-24}$', r'$-10^{-25}$', r'$-10^{-26}$', r'$-10^{-27}$', r'$-10^{-28}$', r'$-10^{-29}$',
                   '0',
                   r'$10^{-29}$', r'$10^{-28}$', r'$10^{-27}$', r'$10^{-26}$', r'$10^{-25}$', r'$10^{-24}$']
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(tick_labels, fontsize=font_size)
    cbar.ax.tick_params(direction='in', length=5, width=1.5)

    # Set a thicker frame (spines)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
    for spine in cbar.ax.spines.values():
        spine.set_linewidth(2.0)

    # Add text annotations
    ax.text(0.03, 0.97,
            f'max = {net_rates.max():.3e} erg/(cm$^3 \cdot$ s)\nmin = {net_rates.min():.3e} erg/(cm$^3 \cdot$ s)',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=font_size*0.9, color='black', fontweight='bold')
    ax.text(0.97, 0.97, rf'$\mathbf{{z = {z_now:.2f}}}$',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=font_size, color='black', fontweight='bold')

    # Save
    plt.savefig(f"grackle_phase_plot_z={z_now}.png", dpi=300, bbox_inches='tight')


# -------------------------------------------------------------------------
def main():
    inspect_grackle_data_file()
    output_one_point()
    output_grid()


if __name__ == '__main__':
    main()
