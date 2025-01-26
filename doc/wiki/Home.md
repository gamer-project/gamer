**GAMER** is a **G**PU-accelerated **A**daptive **ME**sh **R**efinement
code for astrophysics. It features extremely high performance and
parallel scalability and supports a rich set of physics modules.

***


### Physics Modules
* Hydrodynamics
* [Magnetohydrodynamics](https://iopscience.iop.org/article/10.3847/1538-4365/aac49e/meta)
* [Special relativistic hydrodynamics](https://academic.oup.com/mnras/article/504/3/3298/6224873)
* Self-gravity and external gravity
* Particles
* Chemistry and radiative processes with [GRACKLE](http://grackle.readthedocs.io/en/latest/index.html)
* General equation of state
* [Cosmic rays with anisotropic diffusion](https://iopscience.iop.org/article/10.3847/1538-4357/ad50c5#apjad50c5app2)
* Fuzzy (Wave) dark matter: [Nature Physics paper](http://www.nature.com/nphys/journal/v10/n7/covers/index.html), [code paper](https://arxiv.org/abs/2411.17288)

### Other Features
* Adaptive mesh refinement
* Adaptive timestep
* Hybrid MPI/OpenMP/GPU parallelization (also support a CPU-only mode)
* Load balancing with Hilbert space-filling curve
* Bitwise reproducibility
* [HDF5](https://support.hdfgroup.org/HDF5) output
* Data analysis with [yt](http://yt-project.org)
* In situ Python analysis with [libyt](https://github.com/yt-project/libyt)
* Source-term interface
* Feedback interface

### Upcoming Physics and Features
* Anisotropic diffusion
* Anisotropic/Braginskii viscosity

### Need Helps?
* Project leader: Hsi-Yu Schive (hyschive@phys.ntu.edu.tw)
* Mailing list: [GAMER Google Group](https://groups.google.com/forum/#!forum/gamer-amr)
* Live chat: [GAMER Slack](https://join.slack.com/t/gamer-project/shared_invite/enQtNTUwMDA5ODAwMTMzLTc3ZWY2MWE2YTlmMDI0MTQ4M2JjOTg2NmU4OWVkOGY1ZTI3MmY5NjUxOTk1ZjM5ZjNjOGViMGY3ZGExMDdiYzU)
* Code papers:
[GAMER-2](https://academic.oup.com/mnras/article/481/4/4815/5106358) <a name="CODE_PAPER"></a>,
[GAMER-MHD](http://iopscience.iop.org/article/10.3847/1538-4365/aac49e/meta) <a name="MHD_PAPER"></a> ,
[GAMER-SR](https://academic.oup.com/mnras/article/504/3/3298/6224873) <a name="SR_PAPER"></a>,
[GAMER-FDM](https://arxiv.org/abs/2411.17288) <a name="FDM_PAPER"></a>
