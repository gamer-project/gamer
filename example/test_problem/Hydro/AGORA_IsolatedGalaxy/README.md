Compilation flags:
========================================
Enable : MODEL=HYDRO, PARTICLE, GRAVITY, DUAL_ENERGY
Disable: COMOVING


Simulation setup:
========================================
1. Units:
   (1) External (for Input__TestProb only):
       See the comment of each runtime variable

   (2) Internal (for all other input files and internal usage):
       [L] = kpc
       [M] = 1.0e9 Msun
       [T] = Myr
    -->[V] = [L]/[T] ~ 9.8e2 km/s
    -->[D] = [M]/[D]^3 ~ 6.8e-23 g/cm^3

2.  Low-resolution default setup
   (1) Download the low-resolution initial conditions and Input_* by executing "sh download_ic_low_res.sh"
   (2) Default resolution ~ 80 pc (root grid 128^3; MAX_LEVEL = 7)

3.  High-resolution setup
   (1) Download the high-resolution initial conditions and Input_* by executing "sh download_ic_high_res.sh"
   (2) Highest resolution ~ 20 pc (root grid 128^3; MAX_LEVEL = 9)
   (3) This setup reproduces the AGORA high-resolution run presented in the GAMER-2 paper (Schive et al. 2018)


Note:
========================================
1. Mimic the "AgoraGalaxy" test problem setup of Enzo
   --> Nathan Goldbaum, et al., 2015, ApJ, 814, 131 (arXiv: 1510.08458)
       Ji-hoon Kim, et al., 2016, ApJ, 833, 202 (arXiv: 1610.03066)
2. Other references:
       AGORA website           : https://sites.google.com/site/santacruzcomparisonproject/
       AGORA initial conditions: https://goo.gl/8JzbIJ
       Enzo setup              : https://bitbucket.org/enzo/enzo-dev/src/19f4a44e06f1c386573dc77b3608ba66b64d93bc/run/Hydro/Hydro-3D/AgoraGalaxy/?at=week-of-code
       Goldbaum et al. 2016    : https://arxiv.org/abs/1605.00646
       yt hub                  : https://girder.hub.yt/#collection/5736481ddd9119000164acf1
       CloudyData_UVB table    : https://github.com/grackle-project/grackle_data_files/tree/main/input
3. The low-resolution initial condition files from "sh download_ic_low_res.sh" are the same initial conditions used
   in the "AgoraGalaxy" test problem of Enzo
4. Some handy yt analysis scripts are put at "yt_script"


First-time GRACKLE installation guide (on NTU clusters as an example):
========================================
Please refer to the README file at "./example/grackle/"
