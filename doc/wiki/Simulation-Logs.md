This page lists all simulation log files recording information like
simulation performance, load-balancing, memory consumption, time-steps,
and runtime/compilation parameters. Click the filename for details.

## Log Files

| Filename | Description | Option(s) |
|:---|:---|:---|
| [[Record__Center \| Simulation-Logs:-Record__Center]] | Position of maximum density, minimum potential, and center of mass | [[OPT__RECORD_CENTER \| Runtime Parameters:-Miscellaneous#OPT__RECORD_CENTER]] |
| [[Record__Conservation \| Simulation-Logs:-Record__Conservation]] | Integrated values of conservative quantities | [[OPT__CK_CONSERVATION \| Runtime Parameters:-Miscellaneous#OPT__CK_CONSERVATION]] |
| [[Record__Dump \| Simulation-Logs:-Record__Dump]] | Physical time of each data dump | [[OPT__OUTPUT_TOTAL \| Runtime-Parameters:-Outputs#OPT__OUTPUT_TOTAL]], [[OPT__OUTPUT_PART \| Runtime-Parameters:-Outputs#OPT__OUTPUT_PART]], [[OPT__OUTPUT_USER\| Runtime-Parameters:-Outputs#OPT__OUTPUT_USER]] |
| [[Record__LoadBalance \| Simulation-Logs:-Record__LoadBalance]] | Load-balancing estimation | [[OPT__RECORD_LOAD_BALANCE \| Runtime-Parameters:-MPI-and-OpenMP#OPT__RECORD_LOAD_BALANCE]] |
| [[Record__MemInfo \| Simulation-Logs:-Record__MemInfo]] | Free CPU memory | [[OPT__RECORD_MEMORY \| Runtime Parameters:-Miscellaneous#OPT__RECORD_MEMORY]] |
| [[Record__NCorrUnphy \| Simulation-Logs:-Record__NCorrUnphy]] | Number of cells corrected by the 1st-order flux correction | [[OPT__1ST_FLUX_CORR \| Runtime-Parameters:-Hydro#OPT__1ST_FLUX_CORR]], [[OPT__RECORD_UNPHY \| Runtime Parameters:-Miscellaneous#OPT__RECORD_UNPHY]] |
| [[Record__Note \| Simulation-Logs:-Record__Note]] | Runtime parameters, compilation options, CPU/GPU specifications, OpenMP configuration, ... | [[OPT__RECORD_NOTE \| Runtime Parameters:-Miscellaneous#OPT__RECORD_NOTE]] |
| [[Record__PatchCount \| Simulation-Logs:-Record__PatchCount]] | Number of patches on each level in each MPI process | [[OPT__PATCH_COUNT \| Runtime-Parameters:-Refinement#OPT__PATCH_COUNT]] |
| [[Record__ParticleCount \| Simulation-Logs:-Record__ParticleCount]] | Number of particles on each level in each MPI process | [[OPT__PARRTICLE_COUNT \| Runtime-Parameters:-Refinement#OPT__PARTICLE_COUNT]] |
| [[Record__Performance \| Simulation-Logs:-Record__Performance]] | Code performance | [[OPT__RECORD_PERFORMANCE \| Runtime-Parameters:-Miscellaneous#OPT__RECORD_PERFORMANCE]] |
| [[Record__TimeStep \| Simulation-Logs:-Record__TimeStep]] | Time-step constraints | [[OPT__RECORD_DT \| Runtime-Parameters:-Timestep#OPT__RECORD_DT]] |
| [[Record__Timing \| Simulation-Logs:-Record__Timing]] | Detailed timing analysis of all major routines | [[--timing \| Installation:-Option-List#--timing]], [[--timing_solver \| Installation:-Option-List#--timing_solver]] |
| [[Record__TimingMPI_Rank* \| Simulation-Logs:-Record__TimingMPI_Rank*]] | MPI bandwidths achieved by various MPI calls | [[OPT__TIMING_MPI \| Runtime-Parameters:-Miscellaneous#OPT__TIMING_MPI]] |
| [[Record__DivB \| Simulation-Logs:-Record__DivB]] | Divergence-free error on the magnetic field | [[OPT__CK_DIVERGENCE_B \| Runtime-Parameters:-Miscellaneous#OPT__CK_DIVERGENCE_B]] |


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Simulation Snapshots | Outputs]]
