# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[EOS=EOS_GAMMA | Installation: Simulation-Options#EOS]]
   - [[SUPPORT_FFTW | Installation: Simulation-Options#SUPPORT_FFTW]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[MHD | Installation: Simulation-Options#MHD]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. Resolution = 256^3
2. Run for initial conditions only


# Note
1. Calculate the power spectrum of a 3D Gaussian distribution of energy
2. A simple gnuplot script `plot.gpt` is attached
