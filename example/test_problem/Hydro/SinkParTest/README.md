Compilation flags:
========================================
Enable : MODEL=HYDRO, EOS=EOS_USER, MHD, PARTICLE, GRAVITY, BAROTROPIC_EOS, STORE_POT_GHOST, STORE_PAR_ACC, STAR_FORMATION, FEEDBACK, PAR_NATT_INT_USER=1  

Simulation setup:
========================================
1. Units (cgs):  
    [L] = cm  
    [M] = g  
    [T] = s  
  
2. `STAR_FORMATION` parameters:  
    - `SF_CREATE_SINK_MIN_GAS_DENS`: The minimum gas density allowed to form sink particles (count/cm^3) [1.0e10].  
    - `SF_CREATE_SINK_ACC_RADIUS`: The accretion radius in cells at the highest refinement level [0.5*PATCH_SIZE].  
    - `SF_CREATE_SINK_MAX_NPAR_MPI`:  The maximum number of particles per MPI rank [100].  
  
3. `FEEDBACK` parameters:  
    - `FB_ACC`: Turn on particle accretion [0]. Note that `SF_CREATE_SINK_MIN_GAS_DENS` and `SF_CREATE_SINK_ACC_RADIUS` will be used.  
  
4. Create turbulence table:  
    To generate the turbulence table (`Tur_Table.dat`, even the turbulence Mach number is set to 0), run `python generate_TubulentVelocityField.py -n [n] -kmin [kmin] -seed [random seed]`, where  
    - `n`: It is the power law index for turbulent energy power spectrum $\propto k^{n}$ and $k$ is the wave number.  For Kolmogorov and Burgers turbulence, $n = -11/3$ and $-4$ respectively.  
    - `kmin`: It is the minimum wave number in the unit of base cell number.  
    - `random seed`: A random seed in order to get different turbulence realization.  
  
    Example usage: `python generate_TubulentVelocityField.py -n -4 -kmin 1 -seed 0` will generate a turbulence velocity field with energy power spectrum $\propto k^{-4}$.  
  
    The creation of trubulence is following O. Lomax et al., 2015, MNRAS, 449, 662.  

Description:
========================================
This simulation uses Boss & Bodenheimer test (BB test) to test the sink particle and accretion. BB test describes the collapse and fragmentation of a rotating core. The size and mass of the cloud is controlled by `R0` and `Core_Mass`, the angular velocity is set by `Omega0`. A m = 2 density perturbation can be add to the cloud and the strength of the perturbation is controlled by `Delta_Dens`. A density contrast of `Denss_Contrast` is used to scale the cloud density to the background density. 

A magnetic field with strength of `B0` initially along `z`-axis can be added, its inclination respect to the `z`-axis can be set by `theta_B`. A turbulence field described by  `Tur_Table.dat` (see Simulation Setup for how to generate it) is added to the initial velocity field and scale to a Mach number of `Mach_num`.

This simulation uses barotropic EOS, and the transition density is controlled by `rho_AD_SinkParTest` (see Note for more details).

When `Delta_Dens` > 0, the cloud will collapse into to two cores rotating with each other with their own disk, the sink particles will form at the center of the disk and accrete mass from the disk.

Note:
========================================
1. The sink particle and accretion criterias are mostly following C. Federrath et al., 2010, ApJ, 713, 269, with additional criteria from S. D. Clarke et al., 2017, MNRAS, 468, 2489.
2. A custom EOS is used (`CPU_EoS_Barotropic_SinkParTest.cpp`), which descripts a barotropic EOS: $$T(\rho) = T_0 (1 + (\frac{\rho}{\rho_{ad}})^{\gamma - 1})$$, where $\rho$ is gas density, $T_0$ is the initial temperature (`ISO_TEMP`), `gamma` is the adiabatic index (`GAMMA`) and $\rho_{ad}$ is the transition density (`rho_AD_SinkParTest` in `Input_TestProb`), which represents the approximate density at which the gas becomes optically thick (see Hirohiko Masunaga et al., 1998, ApJ, 495, 346, and Hirohiko Masunaga and Shu-ichiro Inutsuka, 2000, ApJ, 531, 350).
3. A script (`Analysis.py`) using `yt` is provided to conduct some analysis:
    - It plots the column density projected on `ProjPlane [xy/xz]` and overplot the sink particles on the maps (`./ColMap/Column_Density_*.png`).
    - The evolution of the two sink particles mass is recorded in `SinkMassEvo.png`.
    - The evolution of mass error is recorded in `MassError.png`.
    - Note that `tff` is the free-fall time calculated from the default parameters in `Input_TestProb`.