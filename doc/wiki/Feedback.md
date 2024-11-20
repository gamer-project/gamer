This page describes feedback from particles to grids and vice versa.
Please enable the compilation option [[--feedback | Installation:-Option-List#--feedback]].


## Compilation Options

Related options:
[[--particle | Installation:-Option-List#--particle]], &nbsp;
[[--feedback | Installation:-Option-List#--feedback]] &nbsp;


## Runtime Parameters

Parameters described on this page:
[FB_LEVEL](#FB_LEVEL), &nbsp;
[FB_RSEED](#FB_RSEED), &nbsp;
[FB_SNE](#FB_SNE), &nbsp;
[FB_USER](#FB_USER) &nbsp;

Other related parameters:

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="FB_LEVEL"></a>
* #### `FB_LEVEL` &ensp; (0 &#8804; input < [[--nlevel | Installation:-Option-List#--nlevel]]; <0 &#8594; set to [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL ]]) &ensp; [-1]
    * **Description:**
AMR level to apply feedback.
    * **Restriction:**
Must be [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL ]] for now.

<a name="FB_RSEED"></a>
* #### `FB_RSEED` &ensp; (&#8805;0) &ensp; [456]
    * **Description:**
Random seed used by feedback.
    * **Restriction:**

<a name="FB_SNE"></a>
* #### `FB_SNE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Supernova explosion feedback.
    * **Restriction:**
Not supported yet.

<a name="FB_USER"></a>
* #### `FB_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
User-defined feedback.
See [Add User-defined Feedback](#add-user-defined-feedback) for details.
    * **Restriction:**


## Remarks

### Add User-defined Feedback
Follow the steps below to define your feedback when
[[adding a new simulation | Adding-New-Simulations]] named `NewProblem`.

1. Go to the new test problem folder and copy the feedback template.

    ```bash
    cd src/TestProblem/Hydro/NewProblem
    cp ../../../Feedback/User_Template/FB_User_Template.cpp FB_NewProblem.cpp
    ```

2. Edit the feedback source file `FB_NewProblem.cpp`.
    1. Rename `User_Template` as `NewProblem`.

    2. Follow the example `src/TestProblem/Hydro/Plummer/FB_Plummer.cpp` to edit
       `FB_Init_NewProblem()`, `FB_End_NewProblem()`, and `FB_NewProblem()`.

3. Edit the problem source file `Init_TestProb_Hydro_NewProblem.cpp` to enable this new feedback.

    1.  Put the following function prototype on the top of this file.

        ```C++
        #ifdef FEEDBACK
        void FB_Init_NewProblem();
        #endif
        ```

    2. Set the feedback function pointer in `Init_TestProb_Hydro_NewProblem()`.

    ```C++
    #  ifdef FEEDBACK
    FB_Init_User_Ptr = FB_Init_NewProblem;
    #  endif
    ```

4. Make sure to enable [[--feedback | Installation:-Option-List#--feedback]] when generating `Makefile` and `FB_USER` in `Input__Parameter`.


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
