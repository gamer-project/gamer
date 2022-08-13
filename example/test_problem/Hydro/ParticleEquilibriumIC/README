Compilation flags:
========================================
Enable : MODEL=HYDRO, GRAVITY, PARTICLE, SUPPORT_GSL
Disable: COMOVING


Default setup:
========================================
OPT__FREEZE_FLUID    1


Note:
========================================
1. You can put many clouds of particles in the HYDRO mode,and observe their evolution through gravity.
   To assign the number of clouds, input the desired number in the variable "ParEqmIC_Cloud_Num" in "Init_TestProb".
   And input the file names of physical parameters in "Input__TestProb" for each cloud respectively.

2. Test the evolution of all sorts of models, like Plummer, NFW, Burkert, Jaffe, Hernquist, Einasto.
   To switch among different models, set "Cloud_Type" in your physical parameter files.
   Available options include:
   Plummer
   NFW
   Burkert
   Jaffe
   Hernquist
   Einasto
   Table
   If you want to initialize a model with a given tabular density profile, use the option "Table" and input
   your file name in "Density_Table_Name" in your physical parameter files; otherwise input "NONE".

3. Sometimes the spatial resolution is not high enough for some models like Jaffe. Please increase MAX_LEVEL in "Input__Parameter".

4. You can construct multiple clouds of particles at the same time.
   For example, replace Input__TestProb with Input_TestProb_Double will allow you to construct two clouds of particles.
   For details, you may check Input__TestProb for further instructions.

5. To generate example density or external potential tables, you may use the python files in "tool/table_maker/".
