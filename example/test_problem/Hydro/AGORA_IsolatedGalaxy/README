Compilation flags:
========================================
Enable : MODEL=HYDRO, PARTICLE, GRAVITY, DUAL_ENERGY
Disable: COMOVING


Default setup:
========================================
1. Units:
   (1) External (for Input__TestProb only):
       See the comment of each runtime variable

   (2) Internal (for all other input files and internal usage):
       [L] = kpc
       [M] = 1.0e8 Msun
       [T] = Myr
    -->[V] = [L]/[T] ~ 9.8e2 km/s
    -->[D] = [M]/[D]^3 ~ 6.8e-23 g/cm^3

2. Default resolution ~ 80 pc (MAX_LEVEL = 7)


Note:
========================================
1. Mimic the "AgoraGalaxy" test problem setup of Enzo
   --> Nathan Goldbaum, et al., 2015, ApJ, 814, 131 (arXiv: 1510.08458)
       Ji-hoon Kim, et al., 2016, ApJ, 833, 202 (arXiv: 1610.03066)
2. Other references:
       AGORA website           : https://sites.google.com/site/santacruzcomparisonproject/
       AGORA initial condition : http://goo.gl/8JzbIJ
       Enzo setup              : https://bitbucket.org/enzo/enzo-dev/src/19f4a44e06f1c386573dc77b3608ba66b64d93bc/run/Hydro/Hydro-3D/AgoraGalaxy/?at=week-of-code
       Goldbaum et al. 2016    : https://arxiv.org/abs/1605.00646
       yt hub                  : https://girder.hub.yt/#collection/5736481ddd9119000164acf1
3. Run the script "download_ic.sh" to download the low-resolution initial condition files for this test
   --> It's the same script used in the "AgoraGalaxy" test problem of Enzo
4. Some handy yt analysis scripts are put at "yt_script"
