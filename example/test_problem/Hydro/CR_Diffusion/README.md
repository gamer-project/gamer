# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[EOS=EOS_GAMMA_CR | Installation: Simulation-Options#EOS]]
   - [[FLOAT8 | Installation: Simulation-Options#FLOAT8]]
   - [[COSMIC_RAY | Installation: Simulation-Options#COSMIC_RAY]]
   - [[CR_DIFFUSION | Installation: Simulation-Options#CR_DIFFUSION]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. `CR_Diffusion_Type = 4`     (Blastwave)
2. `CR_Diffusion_Mag_Type = 4` (Suwa+ 2007)
3. 3D simulation


# Note
1. The analytical solution of `CR_Diffusion_Type={0,1,2,3}` must set all fluid to be fixed but cosmic rays.


## `CR_Diffusion_Type`
1. Gaussian distribution ball (`CR_Diffusion_Type = 0`)
   All of the fluid must be fixed, except for cosmic rays.

2. Step function ring         (`CR_Diffusion_Type = 1`)
   1. All fluids must be fixed, except for cosmic rays.
   2. This test problem is a 2D simulation.
   3. The analytical solution assumes the CRs diffuse along the azimuthal direction only.

3. Gaussian distribution ring (`CR_Diffusion_Type = 2`)
   1. All fluids must be fixed, except for cosmic rays.
   2. This test problem is a 2D simulation.
   3. The analytical solution assumes the CRs diffuse along the azimuthal direction only.

4. Gaussian distribution plane (`CR_Diffusion_Type = 3`)
   1. All fluids must be fixed, except for cosmic rays.
   2. The analytical solution assumes the CRs diffuse along the diagonal direction.

5. CR driven blast wave        (`CR_Diffusion_Type = 4`)
   1. This test problem is a 3D simulation.


## `CR_Diffusion_Mag_Type`
1. Uniform    (`CR_Diffusion_Mag_Type= 0`)
   Uniform magnetic field.
2. Circular   (`CR_Diffusion_Mag_Type= 1`)
   1. Only works for 2D simulation
   2. For x-y plane: `B = (-y/r, x/r, 0)` where `r = sqrt(x^2+y^2)`.
   3. Using `CR_Diffusion_MagX` as amplitude.
3. Random     (`CR_Diffusion_Mag_Type= 2`)
   1. NOT divergence free!
   2. `-1. < B_i < 1.`
4. Radial     (`CR_Diffusion_Mag_Type= 3`)
   1. NOT divergence free!
   2. Only works for 2D simulation
   3. For x-y plane: `B = (x/r, y/r, 0)` where `r = sqrt(x^2+y^2)`.
   4. Using `CR_Diffusion_MagX` as amplitude.
5. Suwa+ 2007 (`CR_Diffusion_Mag_Type= 4`)
   1. Follows [Suwa+ 2007, PASJ, 59, 771](https://doi.org/10.1093/pasj/59.4.771) setup
   2. Using `CR_Diffusion_MagX` as amplitude.


## Set simulation dimensions
The dimensions are controlled by `CR_Diffusion_G{X,Y,Z}`.

For example, the 1D simulation running on the y-axis should set the following:
```
CR_diffusion_GX 0
CR_diffusion_GY 1
CR_diffusion_GZ 0
```

For example, the 2D simulation running on the x-z plane should set the following:
```
CR_diffusion_GX 1
CR_diffusion_GY 0
CR_diffusion_GZ 1
```


## How to fix fluid
Not available right now.


## `yt_L1error`
1. The file structure should be like:
   <pre>
   ./ --+-- res_0032 --+-- Data_000000
        |              |
        |              |-- ...
        |
        |-- res_0064 --+-- Data_000000
        |              |
        |-- ...        |-- ...
   <\pre>
2. Uncomment or comment your simulation type and set the parameters.
