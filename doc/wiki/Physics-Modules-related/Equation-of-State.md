This page describes various equation of states (EoS) supported by the
compilation option [[--eos | Installation:-Option-List#--eos]].

* [`EOS_GAMMA`](#EOS_GAMMA): constant-gamma EoS
* [`EOS_ISOTHERMAL`](#EOS_ISOTHERMAL): isothermal EoS
* [`EOS_COSMIC_RAY`](#EOS_COSMIC_RAY): cosmic-ray EoS
* [`EOS_TAUBMATHEWS`](#EOS_TAUBMATHEWS): special relativistic EoS
* [`EOS_USER`](#EOS_USER): user-specified EoS


## EOS_GAMMA
An ideal-gas EoS with a constant adiabatic index [[GAMMA | Runtime-Parameters:-Hydro#GAMMA]].


## EOS_ISOTHERMAL
An isothermal EoS with a constant sound speed set by
[[MOLECULAR_WEIGHT | Runtime-Parameters:-Hydro#MOLECULAR_WEIGHT]] and
[[ISO_TEMP | Runtime-Parameters:-Hydro#ISO_TEMP]].


## EOS_COSMIC_RAY
A cosmic-ray EoS with an adiabatic index for fluid [[GAMMA | Runtime-Parameters:-Hydro#GAMMA]]
and an effective adiabatic index for cosmic rays [[GAMMA_CR | Runtime-Parameters:-Hydro#GAMMA_CR]].
Must enable [[--cosmic_ray | Installation:-Option-List#--cosmic_ray]].


## EOS_TAUBMATHEWS
A special relativistic EoS with a variable gamma by [Taub 1948](https://ui.adsabs.harvard.edu/abs/1948PhRv...74..328T/abstract) and [Mathews 1971](https://ui.adsabs.harvard.edu/abs/1971ApJ...165..147M/abstract).
Must enable [[--srhd | Installation:-Option-List#--srhd]].


## EOS_USER
Follow the steps below to define your EoS when
[[adding a new simulation | Adding-New-Simulations]] named `NewProblem`.

1. Go to the new test problem folder and copy the EoS template.

    ```bash
    cd src/TestProblem/Hydro/NewProblem
    cp ../../../EoS/User_Template/CPU_EoS_User_Template.cpp CPU_EoS_NewProblem.cpp
    ln -s CPU_EoS_NewProblem.cpp GPU_EoS_NewProblem.cu # CPU/GPU share the same source file
    ```

2. Edit the EoS source file `CPU_EoS_NewProblem.cpp`.
    1. Rename `User_Template` as `NewProblem`. For example, with the `vim` editor
you can do `%s/User_Template/NewProblem/g` in the command line mode.

    2. [Optional] Set the auxiliary array `AuxArray[]` in `void EoS_SetAuxArray_NewProblem()`.
This array must be constant across the simulation period and domain and will be passed to
all EoS conversion functions (see below) automatically. The array size is set by `EOS_NAUX_MAX`
in `include/Macro.h` (default is 10).

    3. **Implement the following EoS conversion functions**.
        * `EoS_DensEint2Pres_NewProblem()`: convert gas mass density and internal energy density to gas pressure.
        * `EoS_DensPres2Eint_NewProblem()`: convert gas mass density and pressure to gas internal energy density.
        * `EoS_DensPres2CSqr_NewProblem()`: convert gas mass density and pressure to sound speed squared.
        * `EoS_DensEint2Temp_NewProblem()`: convert gas mass density and internal energy to gas temperature [OPTIONAL].
        * `EoS_DensTemp2Pres_NewProblem()`: convert gas mass density and temperature to gas pressure [OPTIONAL].
        * `EoS_General_NewProblem()`: general conversion between user-specified input and output variables [OPTIONAL].

> [!CAUTION]
> * All conversion functions must be thread-safe and not use any global variable.
> * When a conversion function fails, it is recommended to return `NAN`
in order to trigger auto-corrections such as [[OPT__1ST_FLUX_CORR | Runtime-Parameters:-Hydro#OPT__1ST_FLUX_CORR]]
and [[AUTO_REDUCE_DT | Runtime-Parameters:-Timestep#AUTO_REDUCE_DT]].

3. Edit the problem source file `Init_TestProb_Hydro_NewProblem.cpp` to enable this new EoS.

    1.  Put the following function prototype on the top of this file.

        ```C++
        void EoS_Init_NewProblem();
        ```

    2. Set the EoS function pointer in `Init_TestProb_Hydro_NewProblem()`.

    ```C++
    EoS_Init_Ptr = EoS_Init_NewProblem;
    ```

4. Make sure to set `EOS=EOS_USER` in the `Makefile`.



<br>

## Links
* [[Main page of Physics Modules | Physics-Modules]]
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Hydro |hydro]]
