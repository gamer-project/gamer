# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--comoving | Installation:-Option-List#--comoving]]
  - [[--mhd | Installation:-Option-List#--mhd]]
  - [[--dual | Installation:-Option-List#--dual]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Adopt mass density as the refinement criteria in `Input__Flag_Rho`
- Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]])=`6`
- Isolated Poisson solver ([[OPT__BC_POT | Runtime-Parameters:-Gravity#OPT__BC_POT]]=`2`)


# Note
- This test problem will overwrite gas field to store the gravitational potential errors
  - gas momentum x --> absolute errors of potential
  - gas momentum y --> relative errors of potential

  Two output files will be generated:
  - `PotError.bin`: binary file similar to `Data_??????` but with the gas field overwritten
  - `PotError.txt`: text   file storing the overwritten data along the diagonal

- Two plot scripts are provided
  - `plot__pot_error_diagonal.gpt`: plot the text file `PotError.txt` using `gnuplot`
  - `plot__pot_error_slice.py`: plot a slice of potential error from `PotError.bin` using `yt`

- Currently supports 2 density profiles: NFW and Hernquist
  - Controlled by the parameter `Gra_DensProf`.
  - Note that errors in NFW can be very large since the total mass in NFW diverges and the adopted analytical
    solution assumes zero potential at infinity.
    - Numerical and analytical solutions will differ by a DC term.

- Set `Gra_NIterProf>0` to measure the average performance of the Poisson solver
  - The measured performance excludes the time for exchanging MPI buffer data.
  - When adopting the isolated Poisson solver (i.e., [[OPT__BC_POT | Runtime-Parameters:-Hydro#OPT__BC_POT]]=`2`), both "NCell" and "Cells/s"
    in `Record__PoissonPerformance` do not count the number of cells in the padded region.
