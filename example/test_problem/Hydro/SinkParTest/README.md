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
    To generate the turbulence table (`Tur_Table.dat`, even the turbulence Mach number is set to 0), run `python VG_turb.py [power law index/2] [kmin] [random seed]`, where  
    - `power law index`: It is the power law index ($\alpha$) for turbulent energy power spectrum $\propto k^{\alpha}$ and $k$ is the wave number.  
    - `kmin`: It is the minimum wave number in the unit  of base cell number.  
    - `random seed`: A random seed in order to get different turbulence realization.  
  
    Example usage: `python VG_turb.py -2 1 0` will generate a turbulence velocity field with energy power spectrum $\propto k^{-4}$.  
  
    The creation of trubulence is following O. Lomax et al., 2015, MNRAS, 449, 662.  
      
Note:
========================================
1. The sink particle and accretion criterias are mostly following Christoph Federrath et al., 2010, ApJ, 713, 269, with additional criteria from S. D. Clarke et al., 2017, MNRAS, 468, 2489.
2. A custom EOS is used (`CPU_EoS_Barotropic_SinkParTest.cpp`), which descripts a barotropic EOS: $$T(\rho) = T_0 (1 + (\frac{\rho}{\rho_{ad}})^{\gamma - 1})$$, where $\rho$ is gas density, $T_0$ is the initial temperature (`ISO_TEMP`), `gamma` is the adiabatic index (`GAMMA`) and $\rho_{ad}$ is the transition density, which represents the approximate density at which the gas becomes optically thick (see Hirohiko Masunaga et al., 1998, ApJ, 495, 346, and Hirohiko Masunaga and Shu-ichiro Inutsuka, 2000, ApJ, 531, 350).
3. A script (`PlotColMap.py`) using `yt` is provided to plot the column density projected on `ProjPlane [xy/xz]` and overplot the sink particles on the maps.
    - Note that `tff` is the free-fall time calculated from the default parameters in `Input_TestPro`.