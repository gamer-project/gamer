# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--dual | Installation:-Option-List#--dual]]=`ENPY`
- Must disable
  - [[--comoving | Installation:-Option-List#--comoving]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Simulation setup
- Units

  - External (for `Input__TestProb` only)

    See the comment of each runtime variable

  - Internal (for all other input files and internal usage)
    | Unit  | [L] | [M]        | [T] | [V] = [L]/[T] | [D] = [M]/[D]^3  |
    |-------|-----|------------|-----|---------------|------------------|
    | Value | kpc | 1.0e9 Msun | Myr | ~ 9.8e2 km/s  | ~ 6.8e-23 g/cm^3 |

- Low-resolution default setup

  1. Download the low-resolution initial conditions and `Input_*` by executing
     ```bash
     sh download_ic_low_res.sh
     ```

  2. Default resolution ~ 80 pc (root grid 128^3; [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]=`7`)

- High-resolution setup

  1. Download the high-resolution initial conditions and `Input_*` by executing
     ```bash
     sh download_ic_high_res.sh
     ```

  2. Highest resolution ~ 20 pc (root grid 128^3; [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]=`9`)

  3. This setup reproduces the AGORA high-resolution run presented in the GAMER-2 paper [Schive et al. 2018](https://academic.oup.com/mnras/article/481/4/4815/5106358)


# Note
- Mimic the "AgoraGalaxy" test problem setup of Enzo.

  See:
  - [Nathan Goldbaum, et al., 2015, ApJ, 814, 131](https://dx.doi.org/10.1088/0004-637X/814/2/131) [(arXiv: 1510.08458)](https://arxiv.org/abs/1510.08458) and
  - [Ji-hoon Kim, et al., 2016, ApJ, 833, 202](https://dx.doi.org/10.3847/1538-4357/833/2/202) [(arXiv: 1610.03066)](https://arxiv.org/abs/1610.03066)

- Other references

  - [AGORA website](https://sites.google.com/site/santacruzcomparisonproject/)
  - [AGORA initial conditions](https://goo.gl/8JzbIJ)
  - [Enzo setup](https://bitbucket.org/enzo/enzo-dev/src/19f4a44e06f1c386573dc77b3608ba66b64d93bc/run/Hydro/Hydro-3D/AgoraGalaxy/?at=week-of-code)
  - [Goldbaum et al. 2016](https://arxiv.org/abs/1605.00646)
  - [yt hub](https://girder.hub.yt/#collection/5736481ddd9119000164acf1)
  - [CloudyData_UVB table](https://github.com/grackle-project/grackle_data_files/tree/main/input)

- The low-resolution initial condition files from `sh download_ic_low_res.sh` are the same initial conditions used
  in the `AgoraGalaxy` test problem of Enzo

- Some handy yt analysis scripts are put at `yt_script`


# First-time GRACKLE installation guide (on NTU clusters as an example)

Please refer to the README file at `./example/grackle/`
