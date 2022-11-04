# Demonstration of phase mismatch after refinement and restriction

## Explanation 
This test problem demonstrates that refining and restricting the wave function in ELBDM leads to a mismatch of the phase on the coarse levels when the option `OPT__INT_PHASE` is on.
Further, it shows that the new option `OPT__RES_PHASE` remedies the mismatch immediately after restriction, but does not prevent the formation of artifacts at the coarse-fine boundary. 
It uses the `GaussianWavePacket` test with suitable initial parameters and refines a region at the center of the Gaussian wave packet. 

## Execution
To run it, compile a version of GAMER with MODEL = ELBDM and put the executable in the test problem directory.
Execute `./compare_restriction.sh` (possibly after making it executable via `chmod +x compare_restriction.sh`).

## Output
The script executes `3` runs for `2` time steps respectively and generates `5` plots.
Three plots titled `Data_%06d_Slice_z_density.png` are created from a run with the new option `OPT__RES_PHASE` turned on.
The grid is unrefined in timestep `0`, refinement and restriction are performed at the end of timestep `1` and the refined fields are evolved during timestep `2`.
The plot named `ComparisonOfRestrictionMethodsBeforeEvolution.png` compares the phase fields on level `0` immediately after refinement and restriction.
The refined fields are not evolved before creating this plot. 
The plots show the phase field on level `0`, their laplacian and the mismatch between the fields from a run with refinement and a run without refinement. 
The phase field after restriction with the standard restriction method in GAMER exhibits a significant mismatch with the run without refinement. 
In other words, refinement and restriction alter the values of the phase field on level `0`. 
The new restriction method does cause a similar mismatch of the phase field on level `0`. 
However, as the plot named `ComparisonOfRestrictionMethodsBeforeEvolution.png` shows, the phase fields on levels `0` and `1` still evolve differently.
The resulting mismatch is bigger than the one caused by the old restriction method.
