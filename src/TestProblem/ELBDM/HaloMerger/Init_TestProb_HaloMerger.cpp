#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double   HaloMerger_Background_Density;                       // background density in the box

static int      HaloMerger_Halo_InitMode;                            // initialization mode: 1=UM_IC real and imaginary parts, 2=UM_IC density-only, 3=density function
static int      HaloMerger_Halo_N;                                   // total number of halos

static char     HaloMerger_Halo_1_IC_Filename[MAX_STRING];            // UM_IC file of halo 1
static double   HaloMerger_Halo_1_IC_L_X;
static double   HaloMerger_Halo_1_IC_L_Y;
static double   HaloMerger_Halo_1_IC_L_Z;
static int      HaloMerger_Halo_1_IC_N_X;
static int      HaloMerger_Halo_1_IC_N_Y;
static int      HaloMerger_Halo_1_IC_N_Z;
static int      HaloMerger_Halo_1_IC_FLOAT8;
static double   HaloMerger_Halo_1_PosX;                               // center coordinates of halo 1
static double   HaloMerger_Halo_1_PosY;                               // center coordinates of halo 1
static double   HaloMerger_Halo_1_PosZ;                               // center coordinates of halo 1
static double   HaloMerger_Halo_1_VelX;                               // velocities of halo 1
static double   HaloMerger_Halo_1_VelY;                               // velocities of halo 1
static double   HaloMerger_Halo_1_VelZ;                               // velocities of halo 1

static char     HaloMerger_Halo_2_IC_Filename[MAX_STRING];            // UM_IC file of halo 2
static double   HaloMerger_Halo_2_IC_L_X;
static double   HaloMerger_Halo_2_IC_L_Y;
static double   HaloMerger_Halo_2_IC_L_Z;
static int      HaloMerger_Halo_2_IC_N_X;
static int      HaloMerger_Halo_2_IC_N_Y;
static int      HaloMerger_Halo_2_IC_N_Z;
static int      HaloMerger_Halo_2_IC_FLOAT8;
static double   HaloMerger_Halo_2_PosX;                               // center coordinates of halo 2
static double   HaloMerger_Halo_2_PosY;                               // center coordinates of halo 2
static double   HaloMerger_Halo_2_PosZ;                               // center coordinates of halo 2
static double   HaloMerger_Halo_2_VelX;                               // velocities of halo 2
static double   HaloMerger_Halo_2_VelY;                               // velocities of halo 2
static double   HaloMerger_Halo_2_VelZ;                               // velocities of halo 2

static char   **IC_Data = NULL;
static char   (*HaloMerger_Halo_IC_Filename)[MAX_STRING] = NULL;
static double (*HaloMerger_Halo_IC_L)[3] = NULL;
static int    (*HaloMerger_Halo_IC_N)[3] = NULL;
static int     *HaloMerger_Halo_IC_FLOAT8 = NULL;
static double (*HaloMerger_Halo_Center)[3] = NULL;                  // center coordinates of each halo
static double (*HaloMerger_Halo_Velocity)[3] = NULL;                // center coordinates of each halo
static double (*HaloMerger_Halo_IC_dh)[3] = NULL;
static double (*HaloMerger_Halo_IC_Range_LCorner)[3] = NULL;
static double (*HaloMerger_Halo_IC_Range_RCorner)[3] = NULL;
// =======================================================================================


static void HaloMerger_Add_Velocity( double *RealPart, double *ImagPart,
                                     const double Velocity_X, const double Velocity_Y, const double Velocity_Z,
                                     const double Position_X, const double Position_Y, const double Position_Z );
static double Trilinear_Interpolation( const double Target_X, const double Target_Y, const double Target_Z,
                                       const double Ref_Value[2][2][2],
                                       const double Ref_X[2], const double Ref_Y[2], const double Ref_Z[2] );



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
   ReadPara->Add( "HaloMerger_Halo_1_IC_Filename",           HaloMerger_Halo_1_IC_Filename,               NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "HaloMerger_Halo_1_IC_L_X",              &HaloMerger_Halo_1_IC_L_X,                  -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_L_Y",              &HaloMerger_Halo_1_IC_L_Y,                  -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_L_Z",              &HaloMerger_Halo_1_IC_L_Z,                  -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_N_X",              &HaloMerger_Halo_1_IC_N_X,                  -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_1_IC_N_Y",              &HaloMerger_Halo_1_IC_N_Y,                  -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_1_IC_N_Z",              &HaloMerger_Halo_1_IC_N_Z,                  -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_1_IC_FLOAT8",           &HaloMerger_Halo_1_IC_FLOAT8,                0,                0,             1              );
   ReadPara->Add( "HaloMerger_Halo_1_PosX",                &HaloMerger_Halo_1_PosX,                    -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_PosY",                &HaloMerger_Halo_1_PosY,                    -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_PosZ",                &HaloMerger_Halo_1_PosZ,                    -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_VelX",                &HaloMerger_Halo_1_VelX,                     0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_VelY",                &HaloMerger_Halo_1_VelY,                     0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_VelZ",                &HaloMerger_Halo_1_VelZ,                     0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_Filename",          HaloMerger_Halo_2_IC_Filename,              NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "HaloMerger_Halo_2_IC_L_X",              &HaloMerger_Halo_2_IC_L_X,                  -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_L_Y",              &HaloMerger_Halo_2_IC_L_Y,                  -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_L_Z",              &HaloMerger_Halo_2_IC_L_Z,                  -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_N_X",              &HaloMerger_Halo_2_IC_N_X,                  -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_2_IC_N_Y",              &HaloMerger_Halo_2_IC_N_Y,                  -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_2_IC_N_Z",              &HaloMerger_Halo_2_IC_N_Z,                  -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_2_IC_FLOAT8",           &HaloMerger_Halo_2_IC_FLOAT8,                0,                0,             1              );
   ReadPara->Add( "HaloMerger_Halo_2_PosX",                &HaloMerger_Halo_2_PosX,                    -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_PosY",                &HaloMerger_Halo_2_PosY,                    -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_PosZ",                &HaloMerger_Halo_2_PosZ,                    -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_VelX",                &HaloMerger_Halo_2_VelX,                     0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_VelY",                &HaloMerger_Halo_2_VelY,                     0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_VelZ",                &HaloMerger_Halo_2_VelZ,                     0.0,              NoMin_double,  NoMax_double   );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
// (2-1) allocate memory
   if ( HaloMerger_Halo_N > 0 )
   {
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         IC_Data                          = new char*  [HaloMerger_Halo_N];
         HaloMerger_Halo_IC_Filename      = new char   [HaloMerger_Halo_N][MAX_STRING];
         HaloMerger_Halo_IC_L             = new double [HaloMerger_Halo_N][3];
         HaloMerger_Halo_IC_N             = new int    [HaloMerger_Halo_N][3];
         HaloMerger_Halo_IC_FLOAT8        = new int    [HaloMerger_Halo_N];
         HaloMerger_Halo_IC_dh            = new double [HaloMerger_Halo_N][3];
         HaloMerger_Halo_IC_Range_LCorner = new double [HaloMerger_Halo_N][3];
         HaloMerger_Halo_IC_Range_RCorner = new double [HaloMerger_Halo_N][3];
      }

      HaloMerger_Halo_Center              = new double [HaloMerger_Halo_N][3];
      HaloMerger_Halo_Velocity            = new double [HaloMerger_Halo_N][3];
   }

   if ( HaloMerger_Halo_N > 0 )
   {
      if ( HaloMerger_Halo_1_PosX < 0 || HaloMerger_Halo_1_PosY < 0 || HaloMerger_Halo_1_PosZ < 0 )
      {
         for (int d=0; d<3; d++)
            HaloMerger_Halo_Center[0][d] = amr->BoxCenter[d];
      }
      else
      {
         HaloMerger_Halo_Center[0][0] = HaloMerger_Halo_1_PosX;
         HaloMerger_Halo_Center[0][1] = HaloMerger_Halo_1_PosY;
         HaloMerger_Halo_Center[0][2] = HaloMerger_Halo_1_PosZ;
      }

      HaloMerger_Halo_Velocity[0][0] = HaloMerger_Halo_1_VelX;
      HaloMerger_Halo_Velocity[0][1] = HaloMerger_Halo_1_VelY;
      HaloMerger_Halo_Velocity[0][2] = HaloMerger_Halo_1_VelZ;

      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         if ( !Aux_CheckFileExist(HaloMerger_Halo_1_IC_Filename) )
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", HaloMerger_Halo_1_IC_Filename );

         if ( HaloMerger_Halo_1_IC_L_X <= 0.0  ||  HaloMerger_Halo_1_IC_L_Y <= 0.0  ||  HaloMerger_Halo_1_IC_L_Z <= 0.0 )
            Aux_Error( ERROR_INFO, "Runtime parameter \"HaloMerger_Halo_1_IC_L_X/Y/Z\" is not set !!\n" );

         if ( HaloMerger_Halo_1_IC_N_X <= 0  ||  HaloMerger_Halo_1_IC_N_Y <= 0  ||   HaloMerger_Halo_1_IC_N_Z <= 0 )
            Aux_Error( ERROR_INFO, "Runtime parameter \"HaloMerger_Halo_1_IC_N_X/Y/Z\" is not set !!\n" );

         strcpy( HaloMerger_Halo_IC_Filename[0], HaloMerger_Halo_1_IC_Filename);
         HaloMerger_Halo_IC_L[0][0]   = HaloMerger_Halo_1_IC_L_X;
         HaloMerger_Halo_IC_L[0][1]   = HaloMerger_Halo_1_IC_L_Y;
         HaloMerger_Halo_IC_L[0][2]   = HaloMerger_Halo_1_IC_L_Z;
         HaloMerger_Halo_IC_N[0][0]   = HaloMerger_Halo_1_IC_N_X;
         HaloMerger_Halo_IC_N[0][1]   = HaloMerger_Halo_1_IC_N_Y;
         HaloMerger_Halo_IC_N[0][2]   = HaloMerger_Halo_1_IC_N_Z;
         HaloMerger_Halo_IC_FLOAT8[0] = HaloMerger_Halo_1_IC_FLOAT8;
      }
   }

   // second halo
   if ( HaloMerger_Halo_N > 1 )
   {
      if ( HaloMerger_Halo_2_PosX < 0 || HaloMerger_Halo_2_PosY < 0 || HaloMerger_Halo_2_PosZ < 0 )
      {
         for (int d=0; d<3; d++)
            HaloMerger_Halo_Center[1][d] = 0.5*amr->BoxCenter[d];
      }
      else
      {
         HaloMerger_Halo_Center[1][0] = HaloMerger_Halo_2_PosX;
         HaloMerger_Halo_Center[1][1] = HaloMerger_Halo_2_PosY;
         HaloMerger_Halo_Center[1][2] = HaloMerger_Halo_2_PosZ;
      }

      HaloMerger_Halo_Velocity[1][0] = HaloMerger_Halo_2_VelX;
      HaloMerger_Halo_Velocity[1][1] = HaloMerger_Halo_2_VelY;
      HaloMerger_Halo_Velocity[1][2] = HaloMerger_Halo_2_VelZ;

      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         if ( !Aux_CheckFileExist(HaloMerger_Halo_2_IC_Filename) )
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", HaloMerger_Halo_2_IC_Filename );

         if ( HaloMerger_Halo_2_IC_L_X <= 0.0  ||  HaloMerger_Halo_2_IC_L_Y <= 0.0  ||  HaloMerger_Halo_2_IC_L_Z <= 0.0 )
            Aux_Error( ERROR_INFO, "Runtime parameter \"HaloMerger_Halo_2_IC_L_X/Y/Z\" is not set !!\n" );

         if ( HaloMerger_Halo_2_IC_N_X <= 0  ||  HaloMerger_Halo_2_IC_N_Y <= 0  ||   HaloMerger_Halo_2_IC_N_Z <= 0 )
            Aux_Error( ERROR_INFO, "Runtime parameter \"HaloMerger_Halo_2_IC_N_X/Y/Z\" is not set !!\n" );

         strcpy( HaloMerger_Halo_IC_Filename[1], HaloMerger_Halo_2_IC_Filename);
         HaloMerger_Halo_IC_L[1][0]   = HaloMerger_Halo_2_IC_L_X;
         HaloMerger_Halo_IC_L[1][1]   = HaloMerger_Halo_2_IC_L_Y;
         HaloMerger_Halo_IC_L[1][2]   = HaloMerger_Halo_2_IC_L_Z;
         HaloMerger_Halo_IC_N[1][0]   = HaloMerger_Halo_2_IC_N_X;
         HaloMerger_Halo_IC_N[1][1]   = HaloMerger_Halo_2_IC_N_Y;
         HaloMerger_Halo_IC_N[1][2]   = HaloMerger_Halo_2_IC_N_Z;
         HaloMerger_Halo_IC_FLOAT8[1] = HaloMerger_Halo_2_IC_FLOAT8;
      }
   }

// (1-3) check the runtime parameters
   if ( HaloMerger_Halo_N > 0 )
   {
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         for (int index_halo=0; index_halo<HaloMerger_Halo_N; index_halo++)
         {
            for (int d=0; d<3; d++)
            {
               HaloMerger_Halo_IC_dh[index_halo][d]            = HaloMerger_Halo_IC_L[index_halo][d]/HaloMerger_Halo_IC_N[index_halo][d];
               HaloMerger_Halo_IC_Range_LCorner[index_halo][d] = HaloMerger_Halo_Center[index_halo][d] - 0.5*HaloMerger_Halo_IC_L[index_halo][d];
               HaloMerger_Halo_IC_Range_RCorner[index_halo][d] = HaloMerger_Halo_Center[index_halo][d] + 0.5*HaloMerger_Halo_IC_L[index_halo][d];

               if ( HaloMerger_Halo_IC_Range_RCorner[index_halo][d] > amr->BoxSize[d]  ||
                    HaloMerger_Halo_IC_Range_LCorner[index_halo][d] < 0.0 )
                  Aux_Error( ERROR_INFO, "IC of halo (index = %d) is outside of simulation box in direction %d !!\n", index_halo, d );
            }

            // check whether the input halos overlap
            for (int index2_halo=0; index2_halo<index_halo; index2_halo++)
            {
               bool isOverlap = true;
               for (int d=0; d<3; d++)
               {
                  if (HaloMerger_Halo_IC_Range_RCorner[index_halo][d]  <= HaloMerger_Halo_IC_Range_LCorner[index2_halo][d] ||
                      HaloMerger_Halo_IC_Range_RCorner[index2_halo][d] <= HaloMerger_Halo_IC_Range_LCorner[index_halo][d] )
                     isOverlap = false;
               }
               if ( isOverlap )
                  Aux_Error( ERROR_INFO, "IC of halo (index = %d) overlaps with IC of halo (index = %d) !!\n", index_halo, index2_halo );
            }
         }
      }
   }


// (2) set the problem-specific derived parameters

// (2-2) read IC data
   if ( HaloMerger_Halo_N > 0 )
   {
      switch ( HaloMerger_Halo_InitMode )
      {
         case 1:
         case 2:
         {
            const long IC_NVar = ( HaloMerger_Halo_InitMode == 1 ) ? 2 : 1; // 1:(Real part & Imag part), 2:(density)

            for (int index_halo=0; index_halo<HaloMerger_Halo_N; index_halo++)
            {
                const long IC_NPoint3D   = (long)HaloMerger_Halo_IC_N[index_halo][0]*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][2];
                size_t load_data_size    = ( HaloMerger_Halo_IC_FLOAT8[index_halo] ) ? sizeof(double) : sizeof(float);
                IC_Data[index_halo]      = new char [ IC_NVar*IC_NPoint3D*load_data_size ];

                // open the file
                FILE *File = fopen( HaloMerger_Halo_IC_Filename[index_halo], "rb" );

                // check the file size
                fseek( File, 0, SEEK_END );
                const long ExpectSize = IC_NPoint3D*IC_NVar*load_data_size;
                const long FileSize   = ftell( File );
                if ( FileSize != ExpectSize )
                   Aux_Error( ERROR_INFO, "size of the IC <%s> (%ld) != expect (%ld) !!\n", HaloMerger_Halo_IC_Filename[index_halo], FileSize, ExpectSize );

                // read the data
                fseek( File, 0, SEEK_SET );
                fread( IC_Data[index_halo], load_data_size, IC_NVar*IC_NPoint3D, File );

                // close the file
                fclose( File );
            }

            break;
         }
         default:
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );
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
      Aux_Message( stdout, "  Background               =%13.6e\n",       HaloMerger_Background_Density );
      Aux_Message( stdout, "  Total number of halos    = %d\n",           HaloMerger_Halo_N );
      Aux_Message( stdout, "  Halo initialization mode = %d\n",           HaloMerger_Halo_InitMode );
      Aux_Message( stdout, "  Halo info:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                   "ID", "Center_X", "Center_Y", "Center_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );
      for (int index_halo=0; index_halo<HaloMerger_Halo_N; index_halo++)
      Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                   index_halo,
                   HaloMerger_Halo_Center[index_halo][0], HaloMerger_Halo_Center[index_halo][1], HaloMerger_Halo_Center[index_halo][2], HaloMerger_Halo_Velocity[index_halo][0], HaloMerger_Halo_Velocity[index_halo][1], HaloMerger_Halo_Velocity[index_halo][2] );
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
      Aux_Message( stdout, "  Halo IC info:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                   "ID", "IC_Filename", "IC_L_X", "IC_L_Y", "IC_L_Z", "IC_N_X", "IC_N_Y", "IC_N_Z", "FLOAT8" );
      for (int index_halo=0; index_halo<HaloMerger_Halo_N; index_halo++)
      Aux_Message( stdout, "  %7d  %13s  %13.6e  %13.6e  %13.6e  %13d  %13d  %13d  %13d\n",
                   index_halo, HaloMerger_Halo_IC_Filename[index_halo],
                   HaloMerger_Halo_IC_L[index_halo][0], HaloMerger_Halo_IC_L[index_halo][1], HaloMerger_Halo_IC_L[index_halo][2], HaloMerger_Halo_IC_N[index_halo][0], HaloMerger_Halo_IC_N[index_halo][1], HaloMerger_Halo_IC_N[index_halo][2], HaloMerger_Halo_IC_FLOAT8[index_halo] );
      }
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

   if ( HaloMerger_Halo_N > 0 )
   {
      switch ( HaloMerger_Halo_InitMode )
      {
         case 1:
         {
            for (int index_halo=0; index_halo<HaloMerger_Halo_N; index_halo++)
            {

               if ( x >= HaloMerger_Halo_IC_Range_LCorner[index_halo][0]  &&  x <= HaloMerger_Halo_IC_Range_RCorner[index_halo][0]  &&
                    y >= HaloMerger_Halo_IC_Range_LCorner[index_halo][1]  &&  y <= HaloMerger_Halo_IC_Range_RCorner[index_halo][1]  &&
                    z >= HaloMerger_Halo_IC_Range_LCorner[index_halo][2]  &&  z <= HaloMerger_Halo_IC_Range_RCorner[index_halo][2] )
               {
                  const int Target_i = (int)floor( (x - HaloMerger_Halo_IC_Range_LCorner[index_halo][0])/HaloMerger_Halo_IC_dh[index_halo][0] - 0.5 );
                  const int Target_j = (int)floor( (y - HaloMerger_Halo_IC_Range_LCorner[index_halo][1])/HaloMerger_Halo_IC_dh[index_halo][1] - 0.5 );
                  const int Target_k = (int)floor( (z - HaloMerger_Halo_IC_Range_LCorner[index_halo][2])/HaloMerger_Halo_IC_dh[index_halo][2] - 0.5 );

                  size_t load_data_size    = ( HaloMerger_Halo_IC_FLOAT8[index_halo] ) ? sizeof(double) : sizeof(float);
                  double Interpolation_Ref_Real[2][2][2];
                  double Interpolation_Ref_Imag[2][2][2];

                  for (int Int_k=0; Int_k<2; Int_k++){ const int k = Target_k + Int_k;
                  for (int Int_j=0; Int_j<2; Int_j++){ const int j = Target_j + Int_j;
                  for (int Int_i=0; Int_i<2; Int_i++){ const int i = Target_i + Int_i;

                     if ( i < 0 ||  i >= HaloMerger_Halo_IC_N[index_halo][0]  ||
                          j < 0 ||  j >= HaloMerger_Halo_IC_N[index_halo][1]  ||
                          k < 0 ||  k >= HaloMerger_Halo_IC_N[index_halo][2] )
                     {
                        Interpolation_Ref_Real[Int_k][Int_j][Int_i] = 0.0;
                        Interpolation_Ref_Imag[Int_k][Int_j][Int_i] = 0.0;
                     }
                     else
                     {
                        if ( HaloMerger_Halo_IC_FLOAT8[index_halo] ){
                          Interpolation_Ref_Real[Int_k][Int_j][Int_i] = (double)(*((double*)&IC_Data[index_halo][((long)0*HaloMerger_Halo_IC_N[index_halo][0]*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][2] + (long)k*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][0] + (long)j*HaloMerger_Halo_IC_N[index_halo][0] + (long)i)*load_data_size]));
                          Interpolation_Ref_Imag[Int_k][Int_j][Int_i] = (double)(*((double*)&IC_Data[index_halo][((long)1*HaloMerger_Halo_IC_N[index_halo][0]*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][2] + (long)k*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][0] + (long)j*HaloMerger_Halo_IC_N[index_halo][0] + (long)i)*load_data_size]));
                        }
                        else
                        {
                          Interpolation_Ref_Real[Int_k][Int_j][Int_i] = (double)(*((float*)&IC_Data[index_halo][((long)0*HaloMerger_Halo_IC_N[index_halo][0]*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][2] + (long)k*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][0] + (long)j*HaloMerger_Halo_IC_N[index_halo][0] + (long)i)*load_data_size]));
                          Interpolation_Ref_Imag[Int_k][Int_j][Int_i] = (double)(*((float*)&IC_Data[index_halo][((long)1*HaloMerger_Halo_IC_N[index_halo][0]*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][2] + (long)k*HaloMerger_Halo_IC_N[index_halo][1]*HaloMerger_Halo_IC_N[index_halo][0] + (long)j*HaloMerger_Halo_IC_N[index_halo][0] + (long)i)*load_data_size]));
                        }
                     }
                  }}}

                  const double Interpolation_Ref_X[2] = {HaloMerger_Halo_IC_Range_LCorner[index_halo][0] + (Target_i + 0.5)*HaloMerger_Halo_IC_dh[index_halo][0], HaloMerger_Halo_IC_Range_LCorner[index_halo][0] + (Target_i + 1.5)*HaloMerger_Halo_IC_dh[index_halo][0]};
                  const double Interpolation_Ref_Y[2] = {HaloMerger_Halo_IC_Range_LCorner[index_halo][1] + (Target_j + 0.5)*HaloMerger_Halo_IC_dh[index_halo][1], HaloMerger_Halo_IC_Range_LCorner[index_halo][1] + (Target_j + 1.5)*HaloMerger_Halo_IC_dh[index_halo][1]};
                  const double Interpolation_Ref_Z[2] = {HaloMerger_Halo_IC_Range_LCorner[index_halo][2] + (Target_k + 0.5)*HaloMerger_Halo_IC_dh[index_halo][2], HaloMerger_Halo_IC_Range_LCorner[index_halo][2] + (Target_k + 1.5)*HaloMerger_Halo_IC_dh[index_halo][2]};

                  // linear interpolation
                  double Real_halo = Trilinear_Interpolation( x, y, z, Interpolation_Ref_Real, Interpolation_Ref_X, Interpolation_Ref_Y, Interpolation_Ref_Z );
                  double Imag_halo = Trilinear_Interpolation( x, y, z, Interpolation_Ref_Imag, Interpolation_Ref_X, Interpolation_Ref_Y, Interpolation_Ref_Z );

                  // add velocity
                  HaloMerger_Add_Velocity( &Real_halo, &Imag_halo, HaloMerger_Halo_Velocity[index_halo][0], HaloMerger_Halo_Velocity[index_halo][1], HaloMerger_Halo_Velocity[index_halo][2], x, y, z );

                  // add the wavefunction to the box
                  Real += Real_halo;
                  Imag += Imag_halo;
               }
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );
      }
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
   fluid[PHAS] = SATAN2( Imag, Real );
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

      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         for (int index_halo=0; index_halo<HaloMerger_Halo_N; index_halo++)
         {
             delete [] IC_Data[index_halo];
         }

         delete [] IC_Data;
         delete [] HaloMerger_Halo_IC_Filename;
         delete [] HaloMerger_Halo_IC_L;
         delete [] HaloMerger_Halo_IC_N;
         delete [] HaloMerger_Halo_IC_FLOAT8;
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



#if ( MODEL == ELBDM )
//-------------------------------------------------------------------------------------------------------
// Function    :  HaloMerger_Add_Velocity
// Description :  Multiply the wave function by a plane wave wave function with a specific velocity
//
// Note        :  1.
//
// Parameter   :  RealPart
//                ImagPart
//                Velocity_X
//                Velocity_Y
//                Velocity_Z
//                Position_X
//                Position_Y
//                Position_Z
// Return      :  RealPart
//                ImagPart
//-------------------------------------------------------------------------------------------------------
void HaloMerger_Add_Velocity( double *RealPart, double *ImagPart,
                              const double Velocity_X, const double Velocity_Y, const double Velocity_Z,
                              const double Position_X, const double Position_Y, const double Position_Z )
{
   const double Real_Old = *RealPart;
   const double Imag_Old = *ImagPart;

   // Phase = kx = (m*v/hbar)*x = eta*v*x
   const double Phase = ELBDM_ETA*( Velocity_X*Position_X + Velocity_Y*Position_Y + Velocity_Z*Position_Z );

   // psi_new = psi_old * exp(i*Phase) = (R_old + i*I_old)*(cos(Phase) + i*sin(Phase))
   const double Real_New = Real_Old * cos(Phase) - Imag_Old * sin(Phase);
   const double Imag_New = Real_Old * sin(Phase) + Imag_Old * cos(Phase);

   // return the updated wave function
   *RealPart = Real_New;
   *ImagPart = Imag_New;

} // FUNCTION : HaloMerger_Add_Velocity
#  endif // if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Trilinear_Interpolation
// Description :  Linear interpolation the desired value from the 3D eigth corners
//
// Note        :  1. Ref_Value is in the order zyx (i.e. Ref_Value[z][y][z] )
//
// Parameter   :  Target_X
//                Target_Y
//                Target_Z
//                Ref_Value
//                Ref_X
//                Ref_Y
//                Ref_Z
// Return      :  Value_ZYX
//-------------------------------------------------------------------------------------------------------
double Trilinear_Interpolation( const double Target_X, const double Target_Y, const double Target_Z,
                                const double Ref_Value[2][2][2],
                                const double Ref_X[2], const double Ref_Y[2], const double Ref_Z[2] )
{
    double Value_Z00 = Ref_Value[0][0][0] + ( Ref_Value[1][0][0] - Ref_Value[0][0][0])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    double Value_Z01 = Ref_Value[0][0][1] + ( Ref_Value[1][0][1] - Ref_Value[0][0][1])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    double Value_Z10 = Ref_Value[0][1][0] + ( Ref_Value[1][1][0] - Ref_Value[0][1][0])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    double Value_Z11 = Ref_Value[0][1][1] + ( Ref_Value[1][1][1] - Ref_Value[0][1][1])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);

    double Value_ZY0 = Value_Z00          + ( Value_Z10          - Value_Z00         )/(Ref_Y[1] - Ref_Y[0])*(Target_Y - Ref_Y[0]);
    double Value_ZY1 = Value_Z01          + ( Value_Z11          - Value_Z01         )/(Ref_Y[1] - Ref_Y[0])*(Target_Y - Ref_Y[0]);

    double Value_ZYX = Value_ZY0          + ( Value_ZY1          - Value_ZY0         )/(Ref_X[1] - Ref_X[0])*(Target_X - Ref_X[0]);

    return Value_ZYX;

}
