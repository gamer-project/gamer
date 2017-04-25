
   *****************************
   ** HYDRO BH accretion test **
   *****************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. Units:
   (1) External (for Input__TestProb only):
       [L   ] = kpc
       [Dens] = g/cm^3
       [Mass] = Msun
       [Temp] = keV

   (2) Internal (for all other input files and internal usage):
       [L   ] = kpc
       [v   ] = 1e3 km/s
       [Dens] = 5.0e-25 g/cm^3

       -->
       [Time] = [L]/[v] ~ 0.978 Myr ~ 1.0 Myr
       [Mass] = [Dens]*[L]^3 ~ 7.4e6 Msun

2. Gravity dt is disabled --> set DT__GRAVITY to arbitrarily large

3. MIN_PRES is set to 1.0e-16

4. Dual-energy formalism is important for evolving pressure accurately in the
   high Mach number region (when Mach number is greater than ~10)
