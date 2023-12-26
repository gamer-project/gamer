#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double   HaloMerger_Background_Density;                       // background density in the box
static int      HaloMerger_Halo_InitMode;                            // initialization mode: 1=UM_IC real and imaginary parts, 2=UM_IC density-only, 3=density function
static int      HaloMerger_Halo_N;                                   // total number of halos
static char     HaloMerger_Halo_IC1_Filename[MAX_STRING];            // UM_IC file of halo 1
static double   HaloMerger_Halo_IC1_L_X;
static double   HaloMerger_Halo_IC1_L_Y;
static double   HaloMerger_Halo_IC1_L_Z;
static int      HaloMerger_Halo_IC1_N_X;
static int      HaloMerger_Halo_IC1_N_Y;
static int      HaloMerger_Halo_IC1_N_Z;
static double   HaloMerger_Halo_PosX1;                               // center coordinates of halo 1
static double   HaloMerger_Halo_PosY1;                               // center coordinates of halo 1
static double   HaloMerger_Halo_PosZ1;                               // center coordinates of halo 1
static double   HaloMerger_Halo_VelX1;                               // velocities of halo 1
static double   HaloMerger_Halo_VelY1;                               // velocities of halo 1
static double   HaloMerger_Halo_VelZ1;                               // velocities of halo 1
static double (*HaloMerger_Halo_Center)[3] = NULL;                  // center coordinates of each halo
static double (*HaloMerger_Halo_Velocity)[3] = NULL;                // center coordinates of each halo
static real    *IC1_Data = NULL;
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__INIT_RESTRICT )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__INIT_RESTRICT !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",                       &VARIABLE,                                   DEFAULT,          MIN,           MAX            );
// ********************************************************************************************************************************
   ReadPara->Add( "HaloMerger_Background_Density",         &HaloMerger_Background_Density,              0.0,              0.0,           NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_InitMode",              &HaloMerger_Halo_InitMode,                   1,                1,             3              );
   ReadPara->Add( "HaloMerger_Halo_N",                     &HaloMerger_Halo_N,                          1,                1,             3              );
   ReadPara->Add( "HaloMerger_Halo_IC1_Filename",           HaloMerger_Halo_IC1_Filename,               NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "HaloMerger_Halo_IC1_L_X",               &HaloMerger_Halo_IC1_L_X,                   -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_IC1_L_Y",               &HaloMerger_Halo_IC1_L_Y,                   -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_IC1_L_Z",               &HaloMerger_Halo_IC1_L_Z,                   -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_IC1_N_X",               &HaloMerger_Halo_IC1_N_X,                   -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_IC1_N_Y",               &HaloMerger_Halo_IC1_N_Y,                   -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_IC1_N_Z",               &HaloMerger_Halo_IC1_N_Z,                   -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_PosX1",                 &HaloMerger_Halo_PosX1,                     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_PosY1",                 &HaloMerger_Halo_PosY1,                     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_PosZ1",                 &HaloMerger_Halo_PosZ1,                     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_VelX1",                 &HaloMerger_Halo_VelX1,                      0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_VelY1",                 &HaloMerger_Halo_VelY1,                      0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_VelZ1",                 &HaloMerger_Halo_VelZ1,                      0.0,              NoMin_double,  NoMax_double   );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( HaloMerger_Halo_N > 0 )
   if ( HaloMerger_Halo_InitMode <= 2 )
   {
      if ( HaloMerger_Halo_IC1_L_X <= 0.0  ||  HaloMerger_Halo_IC1_L_Y <= 0.0  ||   HaloMerger_Halo_IC1_L_Z <= 0.0 )
         Aux_Error( ERROR_INFO, "Runtime parameter \"HaloMerger_Halo_IC1_L_X/Y/Z\" is not set !!\n" );
      if ( HaloMerger_Halo_IC1_N_X <= 0  ||  HaloMerger_Halo_IC1_N_Y <= 0  ||   HaloMerger_Halo_IC1_N_Z <= 0 )
         Aux_Error( ERROR_INFO, "Runtime parameter \"HaloMerger_Halo_IC1_N_X/Y/Z\" is not set !!\n" );
   }

// (2) set the problem-specific derived parameters
// (2-1) allocate memory
   if ( HaloMerger_Halo_N > 0 ){
      HaloMerger_Halo_Center           = new double [HaloMerger_Halo_N][3];
      HaloMerger_Halo_Velocity         = new double [HaloMerger_Halo_N][3];
   }

// (2-2) halo centers
   if ( HaloMerger_Halo_N > 0 )
   {
      if ( HaloMerger_Halo_PosX1 < 0 || HaloMerger_Halo_PosY1 < 0 || HaloMerger_Halo_PosZ1 < 0 ){
         for (int d=0; d<3; d++)    HaloMerger_Halo_Center[0][d] = amr->BoxCenter[d];
      }
      else{
         HaloMerger_Halo_Center[0][0] = HaloMerger_Halo_PosX1;
         HaloMerger_Halo_Center[0][1] = HaloMerger_Halo_PosY1;
         HaloMerger_Halo_Center[0][2] = HaloMerger_Halo_PosZ1;
      }

      HaloMerger_Halo_Velocity[0][0] = HaloMerger_Halo_VelX1;
      HaloMerger_Halo_Velocity[0][1] = HaloMerger_Halo_VelY1;
      HaloMerger_Halo_Velocity[0][2] = HaloMerger_Halo_VelZ1;

      if ( HaloMerger_Halo_InitMode == 1  &&  !Aux_CheckFileExist(HaloMerger_Halo_IC1_Filename) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", HaloMerger_Halo_IC1_Filename );

      if ( HaloMerger_Halo_InitMode == 1 )
      {
         const int IC_NVar = 2;
         FILE *File = fopen( HaloMerger_Halo_IC1_Filename, "rb" );
         fseek( File, 0, SEEK_END );
         const long NPoint3D   = (long)HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z;
         const long ExpectSize = NPoint3D*IC_NVar*sizeof(real);
         const long FileSize   = ftell( File );
         if ( FileSize != ExpectSize )
            Aux_Error( ERROR_INFO, "size of the UM_IC <%s> (%ld) != expect (%ld) !!\n",
                 HaloMerger_Halo_IC1_Filename, FileSize, ExpectSize );
         IC1_Data = new real [ IC_NVar*NPoint3D ];
         fseek( File, 0, SEEK_SET );
         fread( IC1_Data, sizeof(real), IC_NVar*NPoint3D, File );
         fclose( File );
      }
   }



// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID          = %d\n",           TESTPROB_ID );
      Aux_Message( stdout, "  background               = %13.e\n",        HaloMerger_Background_Density );
      Aux_Message( stdout, "  total number of halos    = %d\n",           HaloMerger_Halo_N );
      Aux_Message( stdout, "  halo initialization mode = %d\n",           HaloMerger_Halo_InitMode );
      Aux_Message( stdout, "  Halo info:\n" );
      Aux_Message( stdout, "  %7s  %20s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                   "ID", "IC_Filename", "Center_X", "Center_Y", "Center_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );
      for (int t=0; t<HaloMerger_Halo_N; t++)
      Aux_Message( stdout, "  %7d  %20s  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                   t, HaloMerger_Halo_IC1_Filename,
                   HaloMerger_Halo_Center[t][0], HaloMerger_Halo_Center[t][1], HaloMerger_Halo_Center[t][2], HaloMerger_Halo_Velocity[t][0], HaloMerger_Halo_Velocity[t][1], HaloMerger_Halo_Velocity[t][2] );
      Aux_Message( stdout, "=============================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
   double Dens = HaloMerger_Background_Density;
   double Real = sqrt(Dens);
   double Imag = 0.0;

   switch ( HaloMerger_Halo_InitMode )
   {
      case 1:
      {
         const double IC1_dh   = HaloMerger_Halo_IC1_L_X/HaloMerger_Halo_IC1_N_X;

         const double IC1_Range_x0   = HaloMerger_Halo_PosX1 - 0.5*HaloMerger_Halo_IC1_L_X;
         const double IC1_Range_y0   = HaloMerger_Halo_PosY1 - 0.5*HaloMerger_Halo_IC1_L_Y;
         const double IC1_Range_z0   = HaloMerger_Halo_PosZ1 - 0.5*HaloMerger_Halo_IC1_L_Z;

         const double IC1_Range_x1   = HaloMerger_Halo_PosX1 + 0.5*HaloMerger_Halo_IC1_L_X;
         const double IC1_Range_y1   = HaloMerger_Halo_PosY1 + 0.5*HaloMerger_Halo_IC1_L_Y;
         const double IC1_Range_z1   = HaloMerger_Halo_PosZ1 + 0.5*HaloMerger_Halo_IC1_L_Z;

         if ( IC1_Range_x1 > amr->BoxSize[0]  ||  IC1_Range_x0 < 0.0  ||  IC1_Range_y1 > amr->BoxSize[1]  ||  IC1_Range_y0 < 0.0  ||  IC1_Range_z1 > amr->BoxSize[2]  ||  IC1_Range_z0 < 0.0 )
            Aux_Error( ERROR_INFO, "IC1 is out of range !!\n" );

         if ( x >= IC1_Range_x0  &&  x <= IC1_Range_x1  &&  y >= IC1_Range_y0  &&  y <= IC1_Range_y1  &&  z >= IC1_Range_z0  &&  z <= IC1_Range_z1 )
         {
            const int Target_i = (int)floor( (x - IC1_Range_x0)/IC1_dh - 0.5 );
            const int Target_j = (int)floor( (y - IC1_Range_y0)/IC1_dh - 0.5 );
            const int Target_k = (int)floor( (z - IC1_Range_z0)/IC1_dh - 0.5 );
            double Real_000;
            double Real_001;
            double Real_010;
            double Real_011;
            double Real_100;
            double Real_101;
            double Real_110;
            double Real_111;
            double Imag_000;
            double Imag_001;
            double Imag_010;
            double Imag_011;
            double Imag_100;
            double Imag_101;
            double Imag_110;
            double Imag_111;

            for (int IC_k=0; IC_k<2; IC_k++)
            for (int IC_j=0; IC_j<2; IC_j++)
            for (int IC_i=0; IC_i<2; IC_i++)
            {
               const int i = Target_i + IC_i;
               const int j = Target_j + IC_j;
               const int k = Target_k + IC_k;

               if ( IC_k == 0 && IC_j== 0 && IC_i == 0 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_000 = 0.0;
                     Imag_000 = 0.0;
                  }
                  else
                  {
                     Real_000 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_000 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 0 && IC_j== 0 && IC_i == 1 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_001 = 0.0;
                     Imag_001 = 0.0;
                  }
                  else
                  {
                     Real_001 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_001 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 0 && IC_j== 1 && IC_i == 0 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_010 = 0.0;
                     Imag_010 = 0.0;
                  }
                  else
                  {
                     Real_010 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_010 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 0 && IC_j== 1 && IC_i == 1 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_011 = 0.0;
                     Imag_011 = 0.0;
                  }
                  else
                  {
                     Real_011 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_011 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 1 && IC_j== 0 && IC_i == 0 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_100 = 0.0;
                     Imag_100 = 0.0;
                  }
                  else
                  {
                     Real_100 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_100 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 1 && IC_j== 0 && IC_i == 1 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_101 = 0.0;
                     Imag_101 = 0.0;
                  }
                  else
                  {
                     Real_101 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_101 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 1 && IC_j== 1 && IC_i == 0 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_110 = 0.0;
                     Imag_110 = 0.0;
                  }
                  else
                  {
                     Real_110 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_110 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
               if ( IC_k == 1 && IC_j== 1 && IC_i == 1 )
               {
                  if ( i < 0 ||  i >= HaloMerger_Halo_IC1_N_X  ||  j < 0 ||  j >= HaloMerger_Halo_IC1_N_Y  ||  k < 0 ||  k >= HaloMerger_Halo_IC1_N_Z )
                  {
                     Real_111 = 0.0;
                     Imag_111 = 0.0;
                  }
                  else
                  {
                     Real_111 = IC1_Data[(long)0*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                     Imag_111 = IC1_Data[(long)1*HaloMerger_Halo_IC1_N_X*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_Z + (long)k*HaloMerger_Halo_IC1_N_Y*HaloMerger_Halo_IC1_N_X + (long)j*HaloMerger_Halo_IC1_N_X + (long)i];
                  }
               }
            }

            const double Interpolation_x0 =  IC1_Range_x0 + (Target_i + 0.5)*IC1_dh;
            const double Interpolation_y0 =  IC1_Range_y0 + (Target_j + 0.5)*IC1_dh;
            const double Interpolation_z0 =  IC1_Range_z0 + (Target_k + 0.5)*IC1_dh;
            const double Interpolation_x1 =  IC1_Range_x0 + (Target_i + 1.5)*IC1_dh;
            const double Interpolation_y1 =  IC1_Range_y0 + (Target_j + 1.5)*IC1_dh;
            const double Interpolation_z1 =  IC1_Range_z0 + (Target_k + 1.5)*IC1_dh;

            // linear interpolation

            double Real_z00 = Real_000 + ( Real_100- Real_000)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);
            double Real_z01 = Real_001 + ( Real_101- Real_001)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);
            double Real_z10 = Real_010 + ( Real_110- Real_010)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);
            double Real_z11 = Real_011 + ( Real_111- Real_011)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);

            double Real_zy0 = Real_z00 + ( Real_z10- Real_z00)/(Interpolation_y1 - Interpolation_y0)*(y - Interpolation_y0);
            double Real_zy1 = Real_z01 + ( Real_z11- Real_z01)/(Interpolation_y1 - Interpolation_y0)*(y - Interpolation_y0);

            double Real_zyx = Real_zy0 + ( Real_zy1- Real_zy0)/(Interpolation_x1 - Interpolation_x0)*(x - Interpolation_x0);

            double Imag_z00 = Imag_000 + ( Imag_100- Imag_000)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);
            double Imag_z01 = Imag_001 + ( Imag_101- Imag_001)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);
            double Imag_z10 = Imag_010 + ( Imag_110- Imag_010)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);
            double Imag_z11 = Imag_011 + ( Imag_111- Imag_011)/(Interpolation_z1 - Interpolation_z0)*(z - Interpolation_z0);

            double Imag_zy0 = Imag_z00 + ( Imag_z10- Imag_z00)/(Interpolation_y1 - Interpolation_y0)*(y - Interpolation_y0);
            double Imag_zy1 = Imag_z01 + ( Imag_z11- Imag_z01)/(Interpolation_y1 - Interpolation_y0)*(y - Interpolation_y0);

            double Imag_zyx = Imag_zy0 + ( Imag_zy1- Imag_zy0)/(Interpolation_x1 - Interpolation_x0)*(x - Interpolation_x0);

            // add velocity
            double Phase = ELBDM_ETA*( HaloMerger_Halo_Velocity[0][0]*x + HaloMerger_Halo_Velocity[0][1]*y + HaloMerger_Halo_Velocity[0][2]*z );

            Real += Real_zyx * cos(Phase) - Imag_zyx * sin(Phase);
            Imag += Real_zyx * sin(Phase) + Imag_zyx * cos(Phase);
         }

         break;
      }

      default:
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );
   }



#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Real;
   fluid[IMAG] = Imag;
   fluid[DENS] = SQR(Real) + SQR(Imag);
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else { // if ( amr->use_wave_flag[lv] )
   fluid[DENS] = SQR(Real) + SQR(Imag);
   fluid[PHAS] = 0.0;
   fluid[STUB] = 0.0;
   } // if ( amr->use_wave_flag[lv] ) ... else
#  endif

} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  End_HaloMerger
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_HaloMerger()
{

   if ( HaloMerger_Halo_N > 0 ){
      delete [] HaloMerger_Halo_Center;
      delete [] HaloMerger_Halo_Velocity;

      if ( HaloMerger_Halo_InitMode == 1 )
      {
          delete [] IC1_Data;
      }
   }

} // FUNCTION : End_HaloMerger
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_HaloMerger
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_HaloMerger()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr         = SetGridIC;
   End_User_Ptr                   = End_HaloMerger;
#  endif // if ( MODEL == ELBDM  &&  defined GRAVITY )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_HaloMerger
