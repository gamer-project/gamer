> [!CAUTION]
> Please do not edit this file(page) manually since the workflow will overwrite your changes.
>
> This file(page) is automatically generated by the workflow `Update test problem wiki page` using the script `tool/wiki/sync_test_problem_pages.py`.
>
> The workflow is triggered by push changes to any of `example/test_problem/*/*/README.md` and `tool/wiki/sync_test_problem_pages.py`.


# `configure.py` options:
- Must enable
   - [[--model | Installation:-Option-List#--model]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--double | Installation:-Option-List#--double]]
   - [[--dual | Installation:-Option-List#--dual]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Units
   1. External (for `Input__TestProb` only)

      See the comment of each runtime variable

   2. Internal (for all other input files and internal usage)
      | [L] | [M]        | [T] | [V] = [L]/[T] | [D] = [M]/[D]^3  |
      |---  |---         |---  |---            |---               |
      | kpc | 1.0e8 Msun | Myr | ~ 9.8e2 km/s  | ~ 6.8e-23 g/cm^3 |

2. Default resolution ~ 80 pc ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]] = 7)


# Note
1. Mimic the "AgoraGalaxy" test problem setup of Enzo.

   See: [Nathan Goldbaum, et al., 2015, ApJ, 814, 131](https://dx.doi.org/10.1088/0004-637X/814/2/131) [(arXiv: 1510.08458)](https://arxiv.org/abs/1510.08458) and
   [Ji-hoon Kim, et al., 2016, ApJ, 833, 202](https://dx.doi.org/10.3847/1538-4357/833/2/202) [(arXiv: 1610.03066)](https://dx.doi.org/10.3847/1538-4357/833/2/202)
2. Other references
   - [AGORA website](https://sites.google.com/site/santacruzcomparisonproject/)
   - [AGORA initial condition](https://goo.gl/8JzbIJ)
   - [Enzo setup](https://bitbucket.org/enzo/enzo-dev/src/19f4a44e06f1c386573dc77b3608ba66b64d93bc/run/Hydro/Hydro-3D/AgoraGalaxy/?at=week-of-code)
   - [Goldbaum et al. 2016](https://arxiv.org/abs/1605.00646)
   - [yt hub](https://girder.hub.yt/#collection/5736481ddd9119000164acf1)
3. Run the script `download_ic.sh` to download the low-resolution initial condition files for this test
> [!NOTE]
> It's the same script used in the "AgoraGalaxy" test problem of Enzo.
4. Some handy yt analysis scripts are put at `yt_script`