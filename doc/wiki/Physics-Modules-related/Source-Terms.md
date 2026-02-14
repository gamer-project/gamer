This page describes the source terms.


## Compilation Options

Related options:


## Runtime Parameters
[[ [Runtime parameters] Source Terms | [Runtime-Parameters]-Source-Terms ]]

Other related parameters:
[[SRC_GPU_NPGROUP | [Runtime-Parameters]-GPU#SRC_GPU_NPGROUP]] &nbsp;


## Remarks

### Add User-defined Source Terms
Follow the steps below to define your source terms when
[[adding a new simulation | Adding-New-Simulations]] named `NewProblem`.

1. Go to the new test problem folder and copy the source terms template.

    ```bash
    cd src/TestProblem/Hydro/NewProblem
    cp ../../../SourceTerms/User_Template/CPU_Src_User_Template.cpp Src_NewProblem.cpp
    ```

2. Edit the source terms source file `Src_NewProblem.cpp`.
    1. Rename `User_Template` as `NewProblem`.

TBF.


<br>

## Links
* [[Main page of Physics Modules | Physics-Modules]]
* [[Main page of Runtime Parameters | Runtime Parameters]]
