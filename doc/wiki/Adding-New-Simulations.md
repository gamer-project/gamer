Follow the steps below to add a new simulation.
Mandatory steps are marked by &#x1F4CC;.

1. [Register a New Problem](#i-register-a-new-problem) &#x1F4CC;
2. [Copy the Problem Template](#ii-copy-the-problem-template) &#x1F4CC;
3. [Set Initial Conditions](#iii-set-initial-conditions) &#x1F4CC;
4. [Add Problem-specific Parameters](#iv-add-problem-specific-parameters)
5. [Add Problem-specific Grid Fields and Particle Attributes](#v-add-problem-specific-grid-fields-and-particle-attributes)
6. [Add Problem-specific Functionalities](#vi-add-problem-specific-functionalities)
   *  [Output](#output)
   *  [Initial Condition from Files - Grids](#initial-condition-from-files---grids)
   *  [Initial Condition from Files - Particles](#initial-condition-from-files---particles)
   *  [Work Before Output](#work-before-output)
   *  [Refinement Criteria](#refinement-criteria)
   *  [Work Before Refine](#work-before-refine)
   *  [Timestep Constraints](#timestep-constraints)
   *  [Boundary Conditions](#boundary-conditions)
   *  [Fields Resetting](#fields-resetting)
   *  [Diagnostics](#diagnostics)
   *  [Initialization Function](#initialization-function)
   *  [Finalize Function](#finalize-function)
   *  [External Acceleration](#external-acceleration)
   *  [External Potential](#external-potential)
   *  [Equation of State](#equation-of-state)
   *  [Feedback](#feedback)
7. [Add Problem-specific Validators](#vii-add-problem-specific-validators)
8. [Store Problem-specific Input Files](#viii-store-problem-specific-input-files)


## I. Register a New Problem

1. Add a new problem name and ID to `include/Typedef.h` with the
type `TestProbID_t`. The following example adds
`TESTPROB_HYDRO_NEW_PROBLEM = 123`:

    ```C++
    typedef int TestProbID_t;
    const TestProbID_t
       ...
       TESTPROB_HYDRO_NEW_PROBLEM = 123,
       ...
    ```

2. Edit `src/Init/Init_TestProb.cpp` to register a new problem.
The following example assumes that your new problem initialization
function is called `Init_TestProb_Hydro_NewProblem()`:

    1. Add the new function prototype on the top of this file.

        ```C++
        void Init_TestProb_Hydro_NewProblem();
        ```

    2. Invoke the new problem initializer.

        ```C++
        switch ( TESTPROB_ID )
        {
           ...
           case TESTPROB_HYDRO_NEW_PROBLEM : Init_TestProb_Hydro_NewProblem(); break;
           ...
        }
        ```


## II. Copy the Problem Template

1. Create a new problem directory and file(s).

    ```bash
    cd src/TestProblem/Hydro/
    cp -r ../Template NewProblem
    cd NewProblem
    mv Init_TestProb_Template.cpp Init_TestProb_Hydro_NewProblem.cpp
    ```

2. Edit `Init_TestProb_Hydro_NewProblem.cpp` (hereafter referred to as the
*problem source file*) to replace all strings
`Init_TestProb_Template` by `Init_TestProb_Hydro_NewProblem`.

> [!TIP]
> You can also add additional source files in `src/TestProblem/Hydro/NewProblem`.
They will be compiled automatically without the need to modify the `Makefile`.


## III. Set Initial Conditions
1. Grids IC &#8212; choose one of the following two methods:

    * Specify IC by editing the function `SetGridIC()` in the problem source file.
For details see
[[Setting IC from Analytical Functions &#8212; Grids | Initial-Conditions#IC-Func-Grids]].

    * Load IC from a uniform-mesh binary file. For details see
[[Setting IC from Files &#8212; Grids | Initial-Conditions#IC-File-Grids]].

2. Magnetic field IC &#8212; choose one of the following two methods
(only necessary when enabling
[[--mhd | Installation:-Option-List#--mhd]]):

    * Specify the magnetic field IC by editing the function `SetBFieldIC()`
or specify the vector potential IC by editing a function linking to the function pointer
`Init_BField_ByVecPot_User_Ptr` in the problem source file. For details see
[[Setting IC from Analytical Functions &#8212; Magnetic Field | Initial-Conditions#IC-Func-BField]].
    * Load either the magnetic field or vector potential IC from a
uniform-mesh binary file. For details see
[[Setting IC from Files &#8212; Magnetic Field | Initial-Conditions#IC-File-BField]].
        ${{\color{red}\textsf{Caution:\ magnetic\ field\ is\ not\ currently\ supported\ in\ this\ method.\}}}\$

3. Particles IC &#8212; choose one of the following two methods
(only necessary when enabling
[[--particle | Installation:-Option-List#--particle]]):

    * Specify a particle initialization function.
        1. Define the function `Par_Init_ByFunction_NewProblem()`.
For details see
[[Setting IC from Analytical Functions &#8212; Particles | Initial-Conditions#IC-Func-Particles]].

        2. Put the function prototype on the top of the problem source file.

            ```C++
            #ifdef PARTICLE
            void Par_Init_ByFunction_NewProblem( const long NPar_ThisRank, const long NPar_AllRank,
                                                 real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                                 real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                                 long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                                 long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
            #endif
            ```

        3. Set the function pointer `Par_Init_ByFunction_Ptr` in
`Init_TestProb_Hydro_NewProblem()`.

            ```C++
            #  ifdef PARTICLE
               Par_Init_ByFunction_Ptr = Par_Init_ByFunction_NewProblem;
            #  endif
            ```

    * Load IC from a particle binary file. For details see
[[Setting IC from Files &#8212; Particles | Initial-Conditions#IC-File-Particles]].


## IV. Add Problem-specific Parameters

The following example will load 4 problem-specific runtime parameters
`var_bool`, `var_double`, `var_int` and `var_str` from the input file
`Input__TestProb`.

1. Declare problem-specific global variables on the top of the
problem source file.

    ```C++
    static bool   var_bool;
    static double var_double;
    static int    var_int;
    static char   var_str[MAX_STRING];
    ```

2. Edit the function `SetParameter()` to load these parameters.

    ```C++
       const char FileName[] = "Input__TestProb";
       ReadPara_t *ReadPara  = new ReadPara_t;

    // add parameters in the following format:
    // *************************************************************************************************
    // ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,     DEFAULT,       MIN,            MAX          );
    // *************************************************************************************************
       ReadPara->Add( "var_bool",          &var_bool,     true,          Useless_bool,   Useless_bool );
       ReadPara->Add( "var_double",        &var_double,   1.0,           Eps_double,     NoMax_double );
       ReadPara->Add( "var_int",           &var_int,      2,             0,              5            );
       ReadPara->Add( "var_str",            var_str,      Useless_str,   Useless_str,    Useless_str  );

       ReadPara->Read( FileName );

       delete ReadPara;
    ```

> [!CAUTION]
> `VARIABLE`, `DEFAULT`, `MIN`, and `MAX` must have the same data type.
>
> Some handy constants (e.g., `Useless_bool`, `Eps_double`, `NoMin_int`, ...)
are defined in `include/ReadPara.h`. See [[Adding Parameters | Adding-Parameters]] for details.

3. [Optional] Edit `SetParameter()` to make a note of the values adopted
during the runtime.

    ```C++
    if ( MPI_Rank == 0 )
    {
       Aux_Message( stdout, "=============================================================================\n" );
       Aux_Message( stdout, "  test problem ID = %d\n",     TESTPROB_ID );
       Aux_Message( stdout, "  var_bool        = %d\n",     var_bool );
       Aux_Message( stdout, "  var_double      = %13.7e\n", var_double );
       Aux_Message( stdout, "  var_int         = %d\n",     var_int );
       Aux_Message( stdout, "  var_str         = %s\n",     var_str );
       Aux_Message( stdout, "=============================================================================\n" );
    }
    ```

4. Add these parameters to the input file `Input__TestProb`
(see [[Input__TestProb | Runtime-Parameters:-Input__TestProb]]
for the file format).
This file must be put in the same directory as the executable `gamer`
when running the code.

    ```
    # problem-specific runtime parameters
    var_bool       0             # boolean variable [1]
    var_double     123.0         # double precision variable (>0.0) [1.0]
    var_int        456           # integer variable (0-5) [2]
    var_str        my_string     # string variable
    ```


## V. Add Problem-specific Grid Fields and Particle Attributes

### Grid Fields

It takes 4 small steps to add a new grid field:

1. Set [[ --passive | Installation:-Option-List#--passive]]
to `N` (for `N` new fields) when generating the Makefile.

2. Declare a global integer variable on the top of the problem source
file to store the new field index. For example,

    ```C++
    static int NewFieldIdx = Idx_Undefined;
    ```

    Note that some field index variables have been pre-declared in
`include/Field.h` (e.g., `Idx_Metal` for the field `Metal` used by,
for example, [[ GRACKLE_METAL | Runtime-Parameters:-Chemistry-and-Radiation#GRACKLE_METAL ]]).
Whenever applicable, skip this step and use these pre-declared index variables
directly.

3. Define a function called, for example, `AddNewField_NewProblem()`
and invoke `AddField()` for each of the new field to set the field label
and get the field index. For example,

    ```C++
    void AddNewField_NewProblem()
    {
       if ( NewFieldIdx == Idx_Undefined )
          NewFieldIdx = AddField( "NewFieldLabel", NORMALIZE_YES, INTERP_FRAC_YES );
    }
    ```

    The field label `NewFieldLabel` will be used as the name of this field
in the output files, and the field index `NewFieldIdx` can be used to
access the field data (see the next step). The check `if ( NewFieldIdx == Idx_Undefined )`
is just to avoid redundant assignments to the same field index variable.

    The second parameter should be set to either `NORMALIZE_YES` or `NORMALIZE_NO`.
It controls whether the new field will be renormalized by the total gas density
after every update when enabling
[[ OPT__NORMALIZE_PASSIVE | Runtime-Parameters:-Hydro#OPT__NORMALIZE_PASSIVE ]].

    The third parameter should be set to either `INTERP_FRAC_YES` or `INTERP_FRAC_NO`.
It controls whether the new field will be converted to mass fraction during interpolation
when enabling
[[ OPT__INT_FRAC_PASSIVE_LR | Runtime-Parameters:-Hydro#OPT__INT_FRAC_PASSIVE_LR ]].

    One must also set the function pointer `Init_Field_User_Ptr` in the problem
initialization function `Init_TestProb_Hydro_NewProblem()`.

    ```C++
    Init_Field_User_Ptr = AddNewField_NewProblem;
    ```

> [!NOTE]
> The built-in field `Metal` with the field index `Idx_Metal`
will be added automatically when enabling [[ GRACKLE_METAL | Runtime-Parameters:-Chemistry-and-Radiation#GRACKLE_METAL ]].

4. Assign initial values to the new field in `SetGridIC()` using the corresponding
field index. For example,

    ```C++
    void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                    const int lv, double AuxArray[] )
    {
        ...
        fluid[NewFieldIdx] = 3.7; // it should be mass density instead of mass fraction
        ...
    }
    ```
> [!CAUTION]
> Assign mass density instead of mass fraction to `fluid[NewFieldIdx]`.

### Particle Attributes

Adding a new particle attribute is very similar to adding a new grid field.
So we only highlight the differences in each of the 4 steps above.

1. Set [[ --par_attribute_flt | Installation:-Option-List#--par_attribute_flt ]] and
[[ --par_attribute_int | Installation:-Option-List#--par_attribute_int ]]
instead when generating the Makefile.

2. Declare a global integer variable on the top of the problem source
file to store the new field index. For example,

    ```C++
    static int NewParAttFltIdx = Idx_Undefined;
    static int NewParAttIntIdx = Idx_Undefined;
    ```

Note that some particle attribute index variables have been pre-declared in
`include/Field.h` (e.g., `ParMetalFrac`, which represents the particle
metallicity fraction). Whenever applicable, skip this step and use these
pre-declared index variables directly.

3.  Define a function called, for example, `AddNewParticleAttribute_NewProblem()`
and invoke `AddParticleAttributeFlt()` and `AddParticleAttributeInt()` for each
of the new floating-point and integer attribute, respectively. For example,

    ```C++
    void AddNewParticleAttribute_NewProblem()
    {
       if ( NewParAttFltIdx == Idx_Undefined )
          NewParAttFltIdx = AddParticleAttributeFlt( "NewParFltAttLabel" );
       if ( NewParAttIntIdx == Idx_Undefined )
          NewParAttIntIdx = AddParticleAttributeInt( "NewParIntAttLabel" );
    }
    ```

    The attribute indices `NewParAttFltIdx` and `NewParAttIntIdx` can be used to access the particle
floating-point and integer attribute data, respectively (see the next step). One must also set the function pointer
`Par_Init_Attribute_User_Ptr` in the problem initialization function.

    ```C++
    Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_NewProblem;
    ```

4. Assign initial values to the new particle attribute by using the
corresponding attribute index to access the pointer arrays
`*AllAttributeFlt[PAR_NATT_FLT_TOTAL]` and `*AllAttributeInt[PAR_NATT_INT_TOTAL]` (see
[[Setting IC from Analytical Functions &#8212; Particles | Initial-Conditions#IC-Func-Particles]]).
For example,

    ```C++
    void Par_Init_ByFunction_NewProblem( const long NPar_ThisRank, const long NPar_AllRank,
                                         real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                         real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                         long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                         long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
    {
       ...
       for (long p=0; p<NPar_ThisRank; p++)   AllAttributeFlt[NewParAttFltIdx][p] = 1.0;
       for (long p=0; p<NPar_ThisRank; p++)   AllAttributeInt[NewParAttIntIdx][p] = 2;
       ...
    }
    ```


## VI. Add Problem-specific Functionalities

### Output
The following example illustrates the procedure to add a problem-specific
(i.e., user-specified) data output routine.

1. Define a new data output function called, for example, `Output_NewProblem()`.

    ```C++
    void Output_NewProblem()
    {
       ...
    }
    ```

2. Put its function prototype on the top of the problem source file.

    ```C++
    void Output_NewProblem();
    ```

3. Set the corresponding function pointer in the problem initialization function
`Init_TestProb_Hydro_NewProblem()`.

    ```C++
    Output_User_Ptr = Output_NewProblem;
    ```

4. Turn on the corresponding runtime option
[[OPT__OUTPUT_USER | Runtime-Parameters:-Outputs#OPT__OUTPUT_USER]].
when running the code.

    ```
    OPT__OUTPUT_USER     1
    ```

Other user-specified functionalities such as refinement criteria and
timestep constraints can be added in a similar way and are outlined below.


### Initial Condition from Files - Grids
* **Description:**
Provide a custom routine for [[Setting IC from Files - Grids | Initial-Conditions#IC-File-Grids]].
* **Prototype:**
`void Init_ByFile_NewProblem( real fluid_out[], const real fluid_in[], const int nvar_in,
                              const double x, const double y, const double z, const double Time,
                              const int lv, double AuxArray[] );`
* **Function Pointer:**
`Init_ByFile_User_Ptr`
* **Runtime Option:**
[[OPT__INIT | Runtime-Parameters:-Initial-Conditions#OPT__INIT]]=3
* **Example:**
`src/Init/Init_ByFile.cpp` &#8594; `Init_ByFile_Default()`

### Initial Condition from Files - Particles
* **Description:**
Provide a custom routine for [[Setting IC from Files - Particles | Initial-Conditions#IC-File-Particles]].
* **Prototype:**
`void Par_Init_ByFile_NewProblem();`
* **Function Pointer:**
`Par_Init_ByFile_User_Ptr`
* **Runtime Option:**
[[PAR_INIT | Runtime-Parameters:-Particles#PAR_INIT]]=3
* **Example:**
`src/Particle/Par_Init_ByFile.cpp` &#8594; `Par_Init_ByFile_Default()`
### Work Before Output
* **Description:**
Perform user-specified work before dumping simulation data (e.g., `Data_xxxxxx`).
One common usage is to set grid fields or particle attributes only useful
in post-processing.

* **Prototype:**
`void Output_UserWorkBeforeOutput_NewProblem();`
* **Function Pointer:**
`Output_UserWorkBeforeOutput_Ptr`
* **Runtime Option:**
None
* **Example:**
`src/Output/Output_UserWorkBeforeOutput.cpp`

### Refinement Criteria
* **Description:**
Add user-specified refinement criteria. See
[[OPT__FLAG_USER | Runtime-Parameters:-Refinement#OPT__FLAG_USER]]
for details.
* **Prototype:**
`bool Flag_NewProblem( const int, const int, const int, const int, const int, const double );`
* **Function Pointer:**
`Flag_User_Ptr`
* **Runtime Option:**
[[OPT__FLAG_USER | Runtime-Parameters:-Refinement#OPT__FLAG_USER]]
* **Example:**
`src/Refine/Flag_User.cpp`

### Work Before Refine
* **Description:**
Perform user-specified work before grid refinement.
* **Prototype:**
`void Flag_UserWorkBeforeFlag_NewProblem( const double Time, const int lv );`
* **Function Pointer:**
`Flag_UserWorkBeforeFlag_Ptr`
* **Runtime Option:**
None
* **Example:**
`src/Refine/Flag_UserWorkBeforeFlag.cpp`

### Timestep Constraints
* **Description:**
Add user-specified timestep constraints.
* **Prototype:**
`double Mis_GetTimeStep_User( const int, const double );`
* **Function Pointer:**
`Mis_GetTimeStep_User_Ptr`
* **Runtime Option:**
[[OPT__DT_USER | Runtime-Parameters:-Timestep#OPT__DT_USER]]
* **Example:**
`src/Miscellaneous/Mis_GetTimeStep_User.cpp`

### Boundary Conditions
* **Description:**
Add user-specified (i.e., inflow) boundary conditions.
* **Prototype:**
`void BC_NewProblem( real [], const double, const double, const double, const double, const int, double [] );`
* **Function Pointer:**
`BC_User_Ptr` for the cell-centered fluid variables and
`BC_BField_User_Ptr` for the face-centered magnetic field
* **Runtime Option:**
[[OPT__BC_FLU_* | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]] = 4
* **Example:**
`src/TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp`
&#8594; `BC()`

### Fields Resetting
* **Description:**
  Reset fluid and magnetic fields in specified regions after every sub-step.
This can be used to implement sinks and sources of fluid and magnetic fields.

  To reset the magnetic field, one can specify either the magnetic field directly
(via `MHD_ResetByUser_BField_Ptr`) or the vector potential
(via both `MHD_ResetByUser_BField_Ptr` and `MHD_ResetByNewProblem_VecPot`).
Using the vector potential is recommended since it reduces the divergence-free errors
(see the related restriction below).
Note that one still needs to define `MHD_ResetByUser_BField_Ptr` when using the vector potential
(see the example below).
* **Restriction:**
  * [[INIT_SUBSAMPLING_NCELL | Runtime-Parameters:-Initial-Conditions#INIT_SUBSAMPLING_NCELL]] has no effect on resetting the initial magnetic field.
  * To ensure that the reset magnetic field satisfies the divergence-free condition to the machine precision,
one must (i) use the vector potential and (ii) ensure that the reset fields do not touch any coarse-fine
AMR interfaces. Supporting resetting a divergence-free magnetic field across coarse-fine AMR interfaces
will be implemented in the future.
* **Prototype:**
`int Flu_ResetByNewProblem( real[], const double, const double, const double, const double, const double,
                            const double, const int, double[] );`,
`double MHD_ResetByNewProblem_BField( const double, const double, const double, const double,
                                      const double, const int, const char, double[] );`,
`double MHD_ResetByNewProblem_VecPot( const double, const double, const double, const double,
                                      const double, const int, const char, double[] );`
* **Function Pointer:**
`Flu_ResetByUser_Func_Ptr`, `MHD_ResetByUser_BField_Ptr`, `MHD_ResetByUser_VecPot_Ptr`
* **Runtime Option:**
[[OPT__RESET_FLUID | Runtime-Parameters:-Hydro#OPT__RESET_FLUID]], [[OPT__RESET_FLUID_INIT | Runtime-Parameters:-Hydro#OPT__RESET_FLUID_INIT]]
* **Example:**
`src/Fluid/Flu_ResetByUser.cpp`, `src/Model_Hydro/MHD_ResetByUser.cpp`, `src/TestProblem/Hydro/BlastWave/MHD_ResetByUser_BlastWave.cpp`

### Diagnostics
* **Description:**
Perform user-specified simulation diagnostics.
* **Prototype:**
`void Aux_Record_NewProblem();`
* **Function Pointer:**
`Aux_Record_User_Ptr`
* **Runtime Option:**
[[OPT__RECORD_USER | Runtime-Parameters:-Miscellaneous#OPT__RECORD_USER]]
* **Example:**
`src/Auxiliary/Aux_Record_User.cpp`

### Initialization Function
* **Description:**
These functions will be invoked at the end of the program initialization.
`Init_User_Ptr` and `Init_User_AfterPoisson_Ptr` are invoked before and after
the Poisson solver, respectively.

* **Prototype:**
   * `void Init_NewProblem();`
   * `void Init_AfterPoisson_NewProblem();`
* **Function Pointer:**
   * `Init_User_Ptr`
   * `Init_User_AfterPoisson_Ptr`
* **Runtime Option:**
None
* **Example:**
`src/Init/Init_User.cpp`

### Finalize Function
* **Description:**
This function will be invoked just before the program ends.
One common usage is to free problem-specific memory allocation.
* **Prototype:**
`void End_NewProblem();`
* **Function Pointer:**
`End_User_Ptr`
* **Runtime Option:**
None
* **Example:**
`src/TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp`
&#8594; `End_ClusterMerger()`

### External Acceleration
* **Description:**
Add external acceleration. See
[[External Acceleration/Potential | Gravity#external-accelerationpotential]]
for details.
* **Prototype:**
   * `void Init_ExtAcc_NewProblem();`
* **Function Pointer:**
   * `Init_ExtAcc_Ptr`
* **Runtime Option:**
[[OPT__EXT_ACC | Runtime-Parameters:-Gravity#OPT__EXT_ACC]]
* **Example:**
   * `src/TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp`
   * `src/TestProblem/Hydro/Plummer/ExtAcc_Plummer.cpp`
   * `src/SelfGravity/CPU_Gravity/CPU_ExtAcc_PointMass.cpp`

### External Potential
* **Description:**
Add external potential. See
[[External Acceleration/Potential | Gravity#external-accelerationpotential]]
for details.
* **Prototype:**
   * `void Init_ExtPot_NewProblem();`
* **Function Pointer:**
   * `Init_ExtPot_Ptr`
* **Runtime Option:**
[[OPT__EXT_POT | Runtime-Parameters:-Gravity#OPT__EXT_POT]]
* **Example:**
   * `src/TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp`
   * `src/TestProblem/Hydro/Plummer/ExtPot_Plummer.cpp`
   * `src/SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp`
   * `src/SelfGravity/CPU_Poisson/CPU_ExtPot_Tabular.cpp`

### Equation of State
* **Description:**
Add a user-specified equation of state. See [[here | equation-of-state]] for details.
* **Function Pointer:**
   * `EoS_Init_Ptr`
   * `Eos_End_Ptr`
* **Compilation Option:**
[[--eos | Installation:-Option-List#--eos]]
* **Example:**
   * `src/EoS/User_Template`
   * `src/EoS/Gamma`

### Feedback
* **Description:**
Add a user-specified feedback. See [[FB_USER | Runtime-Parameters:-Feedback#FB_USER]] for details.
* **Function Pointer:**
   * `FB_Init_User_Ptr`
* **Compilation Option:**
[[--feedback | Installation:-Option-List#--feedback]]
* **Runtime Option:**
[[FB_USER | Runtime-Parameters:-Feedback#FB_USER]]
* **Example:**
   * `src/TestProblem/Hydro/Plummer/FB_Plummer.cpp`


## VII. Add Problem-specific Validators

Validate the compilation flags and runtime parameters in the function
`Validate()` of the problem source file. This function will be invoked
during initialization to ensure that the adopted simulation configuration
is compatible with the test problem. For example,


```C++
void Validate()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );

// errors
#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test\n" );
#     endif

      if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
         Aux_Message( stderr, "WARNING : simulation domain is not cubic ??\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate

```


## VIII. Store Problem-specific Input Files

This step is necessary only if you want to add the new simulation
as one of the test problems in GAMER.

1. Create a directory to store problem-specific input files `Input__*`
and other relevant files such as README and analysis scripts.

    ```bash
    cd example/test_problem/Hydro
    cp -rL ../Template NewProblem
    ```

 2. Add and properly set all relevant input files `Input__*`.

3. Edit `README` to help conduct this test problem.

4. Add analysis scripts (if any).


<br>

## Remarks

### Quick Test Without Registering a New Problem
TBF.
