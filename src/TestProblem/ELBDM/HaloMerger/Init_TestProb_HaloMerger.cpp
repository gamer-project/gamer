#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double   HaloMerger_Background_Density;                    // background density in the box

static int      HaloMerger_Halo_Num;                              // total number of halos
static int      HaloMerger_Halo_InitMode;                         // initialization mode:
                                                                  // 1=UM_IC real and imaginary parts, 2=UM_IC density-only, 3=function density profile

static double   HaloMerger_Halo_1_CenCoordX;                      // x coordinate of center of 1st halo
static double   HaloMerger_Halo_1_CenCoordY;                      // y coordinate of center of 1st halo
static double   HaloMerger_Halo_1_CenCoordZ;                      // z coordinate of center of 1st halo
static double   HaloMerger_Halo_1_VelocityX;                      // x component of bulk velocity of 1st halo
static double   HaloMerger_Halo_1_VelocityY;                      // y component of bulk velocity of 1st halo
static double   HaloMerger_Halo_1_VelocityZ;                      // z component of bulk velocity of 1st halo
static char     HaloMerger_Halo_1_IC_Filename[MAX_STRING];        // filename of IC for 1st halo (binary file in VZYX order)
static double   HaloMerger_Halo_1_IC_BoxLenX;                     // physical length of box in x-direction of IC for 1st halo
static double   HaloMerger_Halo_1_IC_BoxLenY;                     // physical length of box in y-direction of IC for 1st halo
static double   HaloMerger_Halo_1_IC_BoxLenZ;                     // physical length of box in z-direction of IC for 1st halo
static int      HaloMerger_Halo_1_IC_NCellsX;                     // number of cells of box in x-direction of IC for 1st halo
static int      HaloMerger_Halo_1_IC_NCellsY;                     // number of cells of box in y-direction of IC for 1st halo
static int      HaloMerger_Halo_1_IC_NCellsZ;                     // number of cells of box in z-direction of IC for 1st halo
static int      HaloMerger_Halo_1_IC_FLOAT8;                      // data precision of IC for 1st halo: 0=float, 1=double

static double   HaloMerger_Halo_2_CenCoordX;                      // x coordinate of center of 2nd halo
static double   HaloMerger_Halo_2_CenCoordY;                      // y coordinate of center of 2nd halo
static double   HaloMerger_Halo_2_CenCoordZ;                      // z coordinate of center of 2nd halo
static double   HaloMerger_Halo_2_VelocityX;                      // x component of bulk velocity of 2nd halo
static double   HaloMerger_Halo_2_VelocityY;                      // y component of bulk velocity of 2nd halo
static double   HaloMerger_Halo_2_VelocityZ;                      // z component of bulk velocity of 2nd halo
static char     HaloMerger_Halo_2_IC_Filename[MAX_STRING];        // filename of IC for 2nd halo (binary file in VZYX order)
static double   HaloMerger_Halo_2_IC_BoxLenX;                     // physical length of box in x-direction of IC for 2nd halo
static double   HaloMerger_Halo_2_IC_BoxLenY;                     // physical length of box in y-direction of IC for 2nd halo
static double   HaloMerger_Halo_2_IC_BoxLenZ;                     // physical length of box in z-direction of IC for 2nd halo
static int      HaloMerger_Halo_2_IC_NCellsX;                     // number of cells of box in x-direction of IC for 2nd halo
static int      HaloMerger_Halo_2_IC_NCellsY;                     // number of cells of box in y-direction of IC for 2nd halo
static int      HaloMerger_Halo_2_IC_NCellsZ;                     // number of cells of box in z-direction of IC for 2nd halo
static int      HaloMerger_Halo_2_IC_FLOAT8;                      // data precision of IC for 2nd halo: 0=float, 1=double

static double   HaloMerger_Halo_3_CenCoordX;                      // x coordinate of center of 3rd halo
static double   HaloMerger_Halo_3_CenCoordY;                      // y coordinate of center of 3rd halo
static double   HaloMerger_Halo_3_CenCoordZ;                      // z coordinate of center of 3rd halo
static double   HaloMerger_Halo_3_VelocityX;                      // x component of bulk velocity of 3rd halo
static double   HaloMerger_Halo_3_VelocityY;                      // y component of bulk velocity of 3rd halo
static double   HaloMerger_Halo_3_VelocityZ;                      // z component of bulk velocity of 3rd halo
static char     HaloMerger_Halo_3_IC_Filename[MAX_STRING];        // filename of IC for 3rd halo (binary file in VZYX order)
static double   HaloMerger_Halo_3_IC_BoxLenX;                     // physical length of box in x-direction of IC for 3rd halo
static double   HaloMerger_Halo_3_IC_BoxLenY;                     // physical length of box in y-direction of IC for 3rd halo
static double   HaloMerger_Halo_3_IC_BoxLenZ;                     // physical length of box in z-direction of IC for 3rd halo
static int      HaloMerger_Halo_3_IC_NCellsX;                     // number of cells of box in x-direction of IC for 3rd halo
static int      HaloMerger_Halo_3_IC_NCellsY;                     // number of cells of box in y-direction of IC for 3rd halo
static int      HaloMerger_Halo_3_IC_NCellsZ;                     // number of cells of box in z-direction of IC for 3rd halo
static int      HaloMerger_Halo_3_IC_FLOAT8;                      // data precision of IC for 3rd halo: 0=float, 1=double

static double (*HaloMerger_Halo_CenCoord)[3]             = NULL;  // center coordinates of each halo
static double (*HaloMerger_Halo_Velocity)[3]             = NULL;  // center coordinates of each halo
static char   **IC_Data = NULL;                                   // array to store the data read from IC
static char   (*HaloMerger_Halo_IC_Filename)[MAX_STRING] = NULL;  // IC filename of each halo
static double (*HaloMerger_Halo_IC_BoxLen)[3]            = NULL;  // length of box of IC of each halo
static int    (*HaloMerger_Halo_IC_NCells)[3]            = NULL;  // number of cells of IC of each halo
static int     *HaloMerger_Halo_IC_FLOAT8                = NULL;  // data precision of IC of each halo
static double (*HaloMerger_Halo_IC_dh)[3]                = NULL;  // grid sizes of each halo
static double (*HaloMerger_Halo_IC_Range_EdgeL)[3]       = NULL;  // left edge of the range of each halo
static double (*HaloMerger_Halo_IC_Range_EdgeR)[3]       = NULL;  // right edge of the range of each halo
// =======================================================================================


static void HaloMerger_Add_Velocity( double *RealPart, double *ImagPart,
                                     const double Velocity_X, const double Velocity_Y, const double Velocity_Z,
                                     const double Position_X, const double Position_Y, const double Position_Z );
static double Trilinear_Interpolation( const double Target_X, const double Target_Y, const double Target_Z,
                                       const double Ref_Value[2][2][2],
                                       const double Ref_X[2], const double Ref_Y[2], const double Ref_Z[2] );

static long IdxInIC( const int v, const int k, const int j, const int i, const int N_v, const int N_k, const int N_j, const int N_i );

static double Get_Value_from_Halo_IC_Data( const double x, const double y, const double z, const int v, const int index_halo );


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
   ReadPara->Add( "HaloMerger_Halo_Num",                   &HaloMerger_Halo_Num,                        1,                1,             3              );
   ReadPara->Add( "HaloMerger_Halo_InitMode",              &HaloMerger_Halo_InitMode,                   1,                1,             3              );
   ReadPara->Add( "HaloMerger_Halo_1_CenCoordX",           &HaloMerger_Halo_1_CenCoordX,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_CenCoordY",           &HaloMerger_Halo_1_CenCoordY,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_CenCoordZ",           &HaloMerger_Halo_1_CenCoordZ,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_VelocityX",           &HaloMerger_Halo_1_VelocityX,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_VelocityY",           &HaloMerger_Halo_1_VelocityY,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_VelocityZ",           &HaloMerger_Halo_1_VelocityZ,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_Filename",          HaloMerger_Halo_1_IC_Filename,              NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "HaloMerger_Halo_1_IC_BoxLenX",          &HaloMerger_Halo_1_IC_BoxLenX,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_BoxLenY",          &HaloMerger_Halo_1_IC_BoxLenY,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_BoxLenZ",          &HaloMerger_Halo_1_IC_BoxLenZ,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_1_IC_NCellsX",          &HaloMerger_Halo_1_IC_NCellsX,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_1_IC_NCellsY",          &HaloMerger_Halo_1_IC_NCellsY,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_1_IC_NCellsZ",          &HaloMerger_Halo_1_IC_NCellsZ,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_1_IC_FLOAT8",           &HaloMerger_Halo_1_IC_FLOAT8,                0,                0,             1              );
   ReadPara->Add( "HaloMerger_Halo_2_CenCoordX",           &HaloMerger_Halo_2_CenCoordX,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_CenCoordY",           &HaloMerger_Halo_2_CenCoordY,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_CenCoordZ",           &HaloMerger_Halo_2_CenCoordZ,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_VelocityX",           &HaloMerger_Halo_2_VelocityX,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_VelocityY",           &HaloMerger_Halo_2_VelocityY,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_VelocityZ",           &HaloMerger_Halo_2_VelocityZ,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_Filename",          HaloMerger_Halo_2_IC_Filename,              NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "HaloMerger_Halo_2_IC_BoxLenX",          &HaloMerger_Halo_2_IC_BoxLenX,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_BoxLenY",          &HaloMerger_Halo_2_IC_BoxLenY,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_BoxLenZ",          &HaloMerger_Halo_2_IC_BoxLenZ,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_2_IC_NCellsX",          &HaloMerger_Halo_2_IC_NCellsX,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_2_IC_NCellsY",          &HaloMerger_Halo_2_IC_NCellsY,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_2_IC_NCellsZ",          &HaloMerger_Halo_2_IC_NCellsZ,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_2_IC_FLOAT8",           &HaloMerger_Halo_2_IC_FLOAT8,                0,                0,             1              );
   ReadPara->Add( "HaloMerger_Halo_3_CenCoordX",           &HaloMerger_Halo_3_CenCoordX,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_CenCoordY",           &HaloMerger_Halo_3_CenCoordY,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_CenCoordZ",           &HaloMerger_Halo_3_CenCoordZ,               -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_VelocityX",           &HaloMerger_Halo_3_VelocityX,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_VelocityY",           &HaloMerger_Halo_3_VelocityY,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_VelocityZ",           &HaloMerger_Halo_3_VelocityZ,                0.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_IC_Filename",          HaloMerger_Halo_3_IC_Filename,              NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "HaloMerger_Halo_3_IC_BoxLenX",          &HaloMerger_Halo_3_IC_BoxLenX,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_IC_BoxLenY",          &HaloMerger_Halo_3_IC_BoxLenY,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_IC_BoxLenZ",          &HaloMerger_Halo_3_IC_BoxLenZ,              -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "HaloMerger_Halo_3_IC_NCellsX",          &HaloMerger_Halo_3_IC_NCellsX,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_3_IC_NCellsY",          &HaloMerger_Halo_3_IC_NCellsY,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_3_IC_NCellsZ",          &HaloMerger_Halo_3_IC_NCellsZ,              -1,                NoMin_int,     NoMax_int      );
   ReadPara->Add( "HaloMerger_Halo_3_IC_FLOAT8",           &HaloMerger_Halo_3_IC_FLOAT8,                0,                0,             1              );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) allocate memory
   if ( HaloMerger_Halo_Num > 0 )
   {
      HaloMerger_Halo_CenCoord            = new double [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_Velocity            = new double [HaloMerger_Halo_Num][3];

      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         IC_Data                          = new char*  [HaloMerger_Halo_Num];
         HaloMerger_Halo_IC_Filename      = new char   [HaloMerger_Halo_Num][MAX_STRING];
         HaloMerger_Halo_IC_BoxLen        = new double [HaloMerger_Halo_Num][3];
         HaloMerger_Halo_IC_NCells        = new int    [HaloMerger_Halo_Num][3];
         HaloMerger_Halo_IC_FLOAT8        = new int    [HaloMerger_Halo_Num];
         HaloMerger_Halo_IC_dh            = new double [HaloMerger_Halo_Num][3];
         HaloMerger_Halo_IC_Range_EdgeL   = new double [HaloMerger_Halo_Num][3];
         HaloMerger_Halo_IC_Range_EdgeR   = new double [HaloMerger_Halo_Num][3];
      }
   }

// (1-3) set the default values
// the 1st halo
   if ( HaloMerger_Halo_Num > 0 )
   {
      // set the center
      if ( HaloMerger_Halo_1_CenCoordX < 0.0  ||  HaloMerger_Halo_1_CenCoordY < 0.0  ||  HaloMerger_Halo_1_CenCoordZ < 0.0 )
      {
         for (int d=0; d<3; d++)
            HaloMerger_Halo_CenCoord[0][d] = amr->BoxCenter[d];
      }
      else if ( HaloMerger_Halo_1_CenCoordX > amr->BoxSize[0]  ||  HaloMerger_Halo_1_CenCoordY > amr->BoxSize[1]  ||  HaloMerger_Halo_1_CenCoordZ > amr->BoxSize[2] )
      {
         Aux_Error( ERROR_INFO, "HaloMerger_Halo_1_CenCoord [%13.6e, %13.6e, %13.6e] is outside of simulation box !!\n", HaloMerger_Halo_1_CenCoordX, HaloMerger_Halo_1_CenCoordY, HaloMerger_Halo_1_CenCoordZ );
      }
      else
      {
         HaloMerger_Halo_CenCoord[0][0] = HaloMerger_Halo_1_CenCoordX;
         HaloMerger_Halo_CenCoord[0][1] = HaloMerger_Halo_1_CenCoordY;
         HaloMerger_Halo_CenCoord[0][2] = HaloMerger_Halo_1_CenCoordZ;
      }

      // set the velocity
      HaloMerger_Halo_Velocity[0][0] = HaloMerger_Halo_1_VelocityX;
      HaloMerger_Halo_Velocity[0][1] = HaloMerger_Halo_1_VelocityY;
      HaloMerger_Halo_Velocity[0][2] = HaloMerger_Halo_1_VelocityZ;

      // set the IC related parameters
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         // check runtime parameters
         if ( !Aux_CheckFileExist(HaloMerger_Halo_1_IC_Filename) )
            Aux_Error( ERROR_INFO, "Halo_1 IC file \"%s\" does not exist !!\n", HaloMerger_Halo_1_IC_Filename );

         if ( HaloMerger_Halo_1_IC_BoxLenX <= 0.0  ||  HaloMerger_Halo_1_IC_BoxLenY <= 0.0  ||  HaloMerger_Halo_1_IC_BoxLenZ <= 0.0 )
            Aux_Error( ERROR_INFO, "HaloMerger_Halo_1_IC_BoxLen [%13.6e, %13.6e, %13.6e] is not set properly !!\n", HaloMerger_Halo_1_IC_BoxLenX, HaloMerger_Halo_1_IC_BoxLenY, HaloMerger_Halo_1_IC_BoxLenZ );

         if ( HaloMerger_Halo_1_IC_NCellsX <= 0  ||  HaloMerger_Halo_1_IC_NCellsY <= 0  ||  HaloMerger_Halo_1_IC_NCellsZ <= 0 )
            Aux_Error( ERROR_INFO, "HaloMerger_Halo_1_IC_NCells [%d, %d, %d] is not set properly !!\n", HaloMerger_Halo_1_IC_NCellsX, HaloMerger_Halo_1_IC_NCellsY, HaloMerger_Halo_1_IC_NCellsZ );

         // set the parameters
         strcpy( HaloMerger_Halo_IC_Filename[0], HaloMerger_Halo_1_IC_Filename );
         HaloMerger_Halo_IC_BoxLen[0][0]       = HaloMerger_Halo_1_IC_BoxLenX;
         HaloMerger_Halo_IC_BoxLen[0][1]       = HaloMerger_Halo_1_IC_BoxLenY;
         HaloMerger_Halo_IC_BoxLen[0][2]       = HaloMerger_Halo_1_IC_BoxLenZ;
         HaloMerger_Halo_IC_NCells[0][0]       = HaloMerger_Halo_1_IC_NCellsX;
         HaloMerger_Halo_IC_NCells[0][1]       = HaloMerger_Halo_1_IC_NCellsY;
         HaloMerger_Halo_IC_NCells[0][2]       = HaloMerger_Halo_1_IC_NCellsZ;
         HaloMerger_Halo_IC_FLOAT8[0]          = HaloMerger_Halo_1_IC_FLOAT8;
      }
   }

// the 2nd halo
   if ( HaloMerger_Halo_Num > 1 )
   {
      // set the center
      if ( HaloMerger_Halo_2_CenCoordX < 0.0  ||  HaloMerger_Halo_2_CenCoordY < 0.0  ||  HaloMerger_Halo_2_CenCoordZ < 0.0 )
      {
         for (int d=0; d<3; d++)
            HaloMerger_Halo_CenCoord[1][d] = 0.5*amr->BoxCenter[d];
      }
      else if ( HaloMerger_Halo_2_CenCoordX > amr->BoxSize[0]  ||  HaloMerger_Halo_2_CenCoordY > amr->BoxSize[1]  ||  HaloMerger_Halo_2_CenCoordZ > amr->BoxSize[2] )
      {
         Aux_Error( ERROR_INFO, "HaloMerger_Halo_2_CenCoord [%13.6e, %13.6e, %13.6e] is outside of simulation box !!\n", HaloMerger_Halo_2_CenCoordX, HaloMerger_Halo_2_CenCoordY, HaloMerger_Halo_2_CenCoordZ );
      }
      else
      {
         HaloMerger_Halo_CenCoord[1][0] = HaloMerger_Halo_2_CenCoordX;
         HaloMerger_Halo_CenCoord[1][1] = HaloMerger_Halo_2_CenCoordY;
         HaloMerger_Halo_CenCoord[1][2] = HaloMerger_Halo_2_CenCoordZ;
      }

      // set the velocity
      HaloMerger_Halo_Velocity[1][0] = HaloMerger_Halo_2_VelocityX;
      HaloMerger_Halo_Velocity[1][1] = HaloMerger_Halo_2_VelocityY;
      HaloMerger_Halo_Velocity[1][2] = HaloMerger_Halo_2_VelocityZ;

      // set the IC related parameters
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         // check runtime parameters
         if ( !Aux_CheckFileExist(HaloMerger_Halo_2_IC_Filename) )
            Aux_Error( ERROR_INFO, "Halo 2 IC file \"%s\" does not exist !!\n", HaloMerger_Halo_2_IC_Filename );

         if ( HaloMerger_Halo_2_IC_BoxLenX <= 0.0  ||  HaloMerger_Halo_2_IC_BoxLenY <= 0.0  ||  HaloMerger_Halo_2_IC_BoxLenZ <= 0.0 )
            Aux_Error( ERROR_INFO, "HaloMerger_Halo_2_IC_BoxLen [%13.6e, %13.6e, %13.6e] is not set properly !!\n", HaloMerger_Halo_2_IC_BoxLenX, HaloMerger_Halo_2_IC_BoxLenY, HaloMerger_Halo_2_IC_BoxLenZ );

         if ( HaloMerger_Halo_2_IC_NCellsX <= 0  ||  HaloMerger_Halo_2_IC_NCellsY <= 0  ||  HaloMerger_Halo_2_IC_NCellsZ <= 0 )
            Aux_Error( ERROR_INFO, "HaloMerger_Halo_2_IC_NCells [%d, %d, %d] is not set properly !!\n", HaloMerger_Halo_2_IC_NCellsX, HaloMerger_Halo_2_IC_NCellsY, HaloMerger_Halo_2_IC_NCellsZ );

         // set the parameters
         strcpy( HaloMerger_Halo_IC_Filename[1], HaloMerger_Halo_2_IC_Filename );
         HaloMerger_Halo_IC_BoxLen[1][0]       = HaloMerger_Halo_2_IC_BoxLenX;
         HaloMerger_Halo_IC_BoxLen[1][1]       = HaloMerger_Halo_2_IC_BoxLenY;
         HaloMerger_Halo_IC_BoxLen[1][2]       = HaloMerger_Halo_2_IC_BoxLenZ;
         HaloMerger_Halo_IC_NCells[1][0]       = HaloMerger_Halo_2_IC_NCellsX;
         HaloMerger_Halo_IC_NCells[1][1]       = HaloMerger_Halo_2_IC_NCellsY;
         HaloMerger_Halo_IC_NCells[1][2]       = HaloMerger_Halo_2_IC_NCellsZ;
         HaloMerger_Halo_IC_FLOAT8[1]          = HaloMerger_Halo_2_IC_FLOAT8;
      }
   }

// the 3rd halo
   if ( HaloMerger_Halo_Num > 2 )
   {
      // set the center
      if ( HaloMerger_Halo_3_CenCoordX < 0.0  ||  HaloMerger_Halo_3_CenCoordY < 0.0  ||  HaloMerger_Halo_3_CenCoordZ < 0.0 )
      {
         for (int d=0; d<3; d++)
            HaloMerger_Halo_CenCoord[2][d] = 1.5*amr->BoxCenter[d];
      }
      else if ( HaloMerger_Halo_3_CenCoordX > amr->BoxSize[0]  ||  HaloMerger_Halo_3_CenCoordY > amr->BoxSize[1]  ||  HaloMerger_Halo_3_CenCoordZ > amr->BoxSize[2] )
      {
         Aux_Error( ERROR_INFO, "HaloMerger_Halo_3_CenCoord [%13.6e, %13.6e, %13.6e] is outside of simulation box !!\n", HaloMerger_Halo_3_CenCoordX, HaloMerger_Halo_3_CenCoordY, HaloMerger_Halo_3_CenCoordZ );
      }
      else
      {
         HaloMerger_Halo_CenCoord[2][0] = HaloMerger_Halo_3_CenCoordX;
         HaloMerger_Halo_CenCoord[2][1] = HaloMerger_Halo_3_CenCoordY;
         HaloMerger_Halo_CenCoord[2][2] = HaloMerger_Halo_3_CenCoordZ;
      }

      // set the velocity
      HaloMerger_Halo_Velocity[2][0] = HaloMerger_Halo_3_VelocityX;
      HaloMerger_Halo_Velocity[2][1] = HaloMerger_Halo_3_VelocityY;
      HaloMerger_Halo_Velocity[2][2] = HaloMerger_Halo_3_VelocityZ;

      // set the IC related parameters
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         // check runtime parameters
         if ( !Aux_CheckFileExist(HaloMerger_Halo_3_IC_Filename) )
            Aux_Error( ERROR_INFO, "Halo 3 IC file \"%s\" does not exist !!\n", HaloMerger_Halo_3_IC_Filename );

         if ( HaloMerger_Halo_3_IC_BoxLenX <= 0.0  ||  HaloMerger_Halo_3_IC_BoxLenY <= 0.0  ||  HaloMerger_Halo_3_IC_BoxLenZ <= 0.0 )
            Aux_Error( ERROR_INFO, "HaloMerger_Halo_3_IC_BoxLen [%13.6e, %13.6e, %13.6e] is not set properly !!\n", HaloMerger_Halo_3_IC_BoxLenX, HaloMerger_Halo_3_IC_BoxLenY, HaloMerger_Halo_3_IC_BoxLenZ );

         if ( HaloMerger_Halo_3_IC_NCellsX <= 0  ||  HaloMerger_Halo_3_IC_NCellsY <= 0  ||  HaloMerger_Halo_3_IC_NCellsZ <= 0 )
            Aux_Error( ERROR_INFO, "HaloMerger_Halo_3_IC_NCells [%d, %d, %d] is not set properly !!\n", HaloMerger_Halo_3_IC_NCellsX, HaloMerger_Halo_3_IC_NCellsY, HaloMerger_Halo_3_IC_NCellsZ );

         // set the parameters
         strcpy( HaloMerger_Halo_IC_Filename[2], HaloMerger_Halo_3_IC_Filename );
         HaloMerger_Halo_IC_BoxLen[2][0]       = HaloMerger_Halo_3_IC_BoxLenX;
         HaloMerger_Halo_IC_BoxLen[2][1]       = HaloMerger_Halo_3_IC_BoxLenY;
         HaloMerger_Halo_IC_BoxLen[2][2]       = HaloMerger_Halo_3_IC_BoxLenZ;
         HaloMerger_Halo_IC_NCells[2][0]       = HaloMerger_Halo_3_IC_NCellsX;
         HaloMerger_Halo_IC_NCells[2][1]       = HaloMerger_Halo_3_IC_NCellsY;
         HaloMerger_Halo_IC_NCells[2][2]       = HaloMerger_Halo_3_IC_NCellsZ;
         HaloMerger_Halo_IC_FLOAT8[2]          = HaloMerger_Halo_3_IC_FLOAT8;
      }
   }

// (1-3) set the problem-specific derived parameters
   if ( HaloMerger_Halo_Num > 0 )
   {
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
         {
            for (int d=0; d<3; d++)
            {
               HaloMerger_Halo_IC_dh[index_halo][d]            = HaloMerger_Halo_IC_BoxLen[index_halo][d]/HaloMerger_Halo_IC_NCells[index_halo][d];
               HaloMerger_Halo_IC_Range_EdgeL[index_halo][d]   = HaloMerger_Halo_CenCoord[index_halo][d] - 0.5*HaloMerger_Halo_IC_BoxLen[index_halo][d];
               HaloMerger_Halo_IC_Range_EdgeR[index_halo][d]   = HaloMerger_Halo_CenCoord[index_halo][d] + 0.5*HaloMerger_Halo_IC_BoxLen[index_halo][d];

               // check whether the input halos cross the boundary
               if ( HaloMerger_Halo_IC_Range_EdgeR[index_halo][d] > amr->BoxSize[d]  ||
                    HaloMerger_Halo_IC_Range_EdgeL[index_halo][d] < 0.0 )
                  Aux_Error( ERROR_INFO, "The edge in direction-%d [%13.6e, %13.6e] of Halo_%d IC range is outside of simulation box !!\n",
                             d, HaloMerger_Halo_IC_Range_EdgeL[index_halo][d], HaloMerger_Halo_IC_Range_EdgeR[index_halo][d], index_halo+1 );
            }

            // check whether the input halos overlap
            for (int index2_halo=0; index2_halo<index_halo; index2_halo++)
            {
               bool isOverlap = true;

               for (int d=0; d<3; d++)
               {
                  if (HaloMerger_Halo_IC_Range_EdgeR[index_halo][d]  <= HaloMerger_Halo_IC_Range_EdgeL[index2_halo][d] ||
                      HaloMerger_Halo_IC_Range_EdgeR[index2_halo][d] <= HaloMerger_Halo_IC_Range_EdgeL[index_halo][d] )
                     isOverlap = false;
               }

               if ( isOverlap )
                  Aux_Error( ERROR_INFO, "Halo_%d IC range overlaps with the Halo_%d IC range !!\n", index_halo+1, index2_halo+1 );

            }
         }
      }
   }


// (2) Load the IC data

// (2-1) Halos
   if ( HaloMerger_Halo_Num > 0 )
   {
      switch ( HaloMerger_Halo_InitMode )
      {
         case 1:
         case 2:
         {
            const long IC_NVar = ( HaloMerger_Halo_InitMode == 1 ) ? 2 : 1; // 1:(Real part & Imag part), 2:(Density)

            for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
            {
                const long IC_NCells3D   = (long)HaloMerger_Halo_IC_NCells[index_halo][0]*HaloMerger_Halo_IC_NCells[index_halo][1]*HaloMerger_Halo_IC_NCells[index_halo][2];
                size_t load_data_size    = ( HaloMerger_Halo_IC_FLOAT8[index_halo] ) ? sizeof(double) : sizeof(float);
                IC_Data[index_halo]      = new char [ IC_NVar*IC_NCells3D*load_data_size ];

                // open the file
                FILE *File = fopen( HaloMerger_Halo_IC_Filename[index_halo], "rb" );

                // check the file size
                fseek( File, 0, SEEK_END );
                const long ExpectSize = IC_NVar*IC_NCells3D*load_data_size;
                const long FileSize   = ftell( File );
                if ( FileSize != ExpectSize )
                   Aux_Error( ERROR_INFO, "size of the IC <%s> (%ld) != expect (%ld) !!\n", HaloMerger_Halo_IC_Filename[index_halo], FileSize, ExpectSize );

                // load data from the file
                fseek( File, 0, SEEK_SET );
                fread( IC_Data[index_halo], 1, IC_NVar*IC_NCells3D*load_data_size, File );

                // close the file
                fclose( File );
            }

            break;
         }
         default:
            Aux_Error( ERROR_INFO, "unsupported halo initialization mode (%s = %d) !!\n",
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
      Aux_Message( stdout, "  background density       =%13.6e\n",        HaloMerger_Background_Density );
      Aux_Message( stdout, "  total number of halos    = %d\n",           HaloMerger_Halo_Num      );
      Aux_Message( stdout, "  halo initialization mode = %d\n",           HaloMerger_Halo_InitMode );
      Aux_Message( stdout, "  halo information:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                   "ID", "CenCoord_X", "CenCoord_Y", "CenCoord_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );
      for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
      Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                   index_halo+1,
                   HaloMerger_Halo_CenCoord[index_halo][0], HaloMerger_Halo_CenCoord[index_halo][1], HaloMerger_Halo_CenCoord[index_halo][2],
                   HaloMerger_Halo_Velocity[index_halo][0], HaloMerger_Halo_Velocity[index_halo][1], HaloMerger_Halo_Velocity[index_halo][2] );
      if ( HaloMerger_Halo_InitMode <= 2 )
      {
      Aux_Message( stdout, "  halo IC information:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %17s  %17s  %17s  %17s  %17s  %17s\n",
                   "ID", "IC_Filename", "IC_FLOAT8",
                         "IC_BoxLen_X", "IC_BoxLen_Y", "IC_BoxLen_Z",
                         "IC_NCells_X", "IC_NCells_Y", "IC_NCells_Z",
                         "IC_dh_X", "IC_dh_Y", "IC_dh_Z",
                         "IC_Range_EdgeL_X", "IC_Range_EdgeL_Y", "IC_Range_EdgeL_Z",
                         "IC_Range_EdgeR_X", "IC_Range_EdgeR_Y", "IC_Range_EdgeR_Z" );
      for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
      Aux_Message( stdout, "  %7d  %13s  %13d  %13.6e  %13.6e  %13.6e  %13d  %13d  %13d  %13.6e  %13.6e  %13.6e  %17.6e  %17.6e  %17.6e  %17.6e  %17.6e  %17.6e\n",
                   index_halo+1, HaloMerger_Halo_IC_Filename[index_halo], HaloMerger_Halo_IC_FLOAT8[index_halo],
                   HaloMerger_Halo_IC_BoxLen[index_halo][0], HaloMerger_Halo_IC_BoxLen[index_halo][1], HaloMerger_Halo_IC_BoxLen[index_halo][2],
                   HaloMerger_Halo_IC_NCells[index_halo][0], HaloMerger_Halo_IC_NCells[index_halo][1], HaloMerger_Halo_IC_NCells[index_halo][2],
                   HaloMerger_Halo_IC_dh[index_halo][0], HaloMerger_Halo_IC_dh[index_halo][1], HaloMerger_Halo_IC_dh[index_halo][2],
                   HaloMerger_Halo_IC_Range_EdgeL[index_halo][0], HaloMerger_Halo_IC_Range_EdgeL[index_halo][1], HaloMerger_Halo_IC_Range_EdgeL[index_halo][2],
                   HaloMerger_Halo_IC_Range_EdgeR[index_halo][0], HaloMerger_Halo_IC_Range_EdgeR[index_halo][1], HaloMerger_Halo_IC_Range_EdgeR[index_halo][2]);
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
   // (1) set the background
   double Dens = HaloMerger_Background_Density;
   double Real = sqrt(Dens);
   double Imag = 0.0;

   // (2) set the halos
   if ( HaloMerger_Halo_Num > 0 )
   {
      switch ( HaloMerger_Halo_InitMode )
      {
         case 1:
         case 2:
         {
            for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
            {
               // when the (x,y,z) is inside the range of this halo
               if ( x >= HaloMerger_Halo_IC_Range_EdgeL[index_halo][0]  &&  x <= HaloMerger_Halo_IC_Range_EdgeR[index_halo][0]  &&
                    y >= HaloMerger_Halo_IC_Range_EdgeL[index_halo][1]  &&  y <= HaloMerger_Halo_IC_Range_EdgeR[index_halo][1]  &&
                    z >= HaloMerger_Halo_IC_Range_EdgeL[index_halo][2]  &&  z <= HaloMerger_Halo_IC_Range_EdgeR[index_halo][2] )
               {
                  double Real_halo, Imag_halo;
                  if ( HaloMerger_Halo_InitMode == 1 )
                  {
                     Real_halo = Get_Value_from_Halo_IC_Data( x, y, z, 0, index_halo );
                     Imag_halo = Get_Value_from_Halo_IC_Data( x, y, z, 1, index_halo );
                  }
                  else if ( HaloMerger_Halo_InitMode == 2 )
                  {
                     double Dens_halo;
                     Dens_halo = Get_Value_from_Halo_IC_Data( x, y, z, 0, index_halo );
                     Real_halo = sqrt(Dens_halo);
                     Imag_halo = 0.0;
                  }
                  else
                      Aux_Error( ERROR_INFO, "unsupported halo initialization mode (%s = %d) !!\n",
                                 "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

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

   if ( HaloMerger_Halo_Num > 0 ){
      delete [] HaloMerger_Halo_CenCoord;
      delete [] HaloMerger_Halo_Velocity;

      if ( HaloMerger_Halo_InitMode <= 2 )
      {
         for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
         {
             delete [] IC_Data[index_halo];
         }

         delete [] IC_Data;
         delete [] HaloMerger_Halo_IC_Filename;
         delete [] HaloMerger_Halo_IC_BoxLen;
         delete [] HaloMerger_Halo_IC_NCells;
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
    // linear interpolation in z-direction
    double Value_Z00 = Ref_Value[0][0][0] + ( Ref_Value[1][0][0] - Ref_Value[0][0][0])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    double Value_Z01 = Ref_Value[0][0][1] + ( Ref_Value[1][0][1] - Ref_Value[0][0][1])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    double Value_Z10 = Ref_Value[0][1][0] + ( Ref_Value[1][1][0] - Ref_Value[0][1][0])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    double Value_Z11 = Ref_Value[0][1][1] + ( Ref_Value[1][1][1] - Ref_Value[0][1][1])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);

    // linear interpolation in y-direction
    double Value_ZY0 = Value_Z00          + ( Value_Z10          - Value_Z00         )/(Ref_Y[1] - Ref_Y[0])*(Target_Y - Ref_Y[0]);
    double Value_ZY1 = Value_Z01          + ( Value_Z11          - Value_Z01         )/(Ref_Y[1] - Ref_Y[0])*(Target_Y - Ref_Y[0]);

    // linear interpolation in x-direction
    double Value_ZYX = Value_ZY0          + ( Value_ZY1          - Value_ZY0         )/(Ref_X[1] - Ref_X[0])*(Target_X - Ref_X[0]);

    return Value_ZYX;

}



//-------------------------------------------------------------------------------------------------------
// Function    :  IdxInIC
// Description :  Get index in the IC_Data
//
// Note        :  1. in the orderi vzyx
//
// Parameter   :  v
//                k
//                j
//                i
//                N_v
//                N_k
//                N_j
//                N_i
// Return      :  Index_1D
//-------------------------------------------------------------------------------------------------------
long IdxInIC( const int v, const int k, const int j, const int i, const int N_v, const int N_k, const int N_j, const int N_i )
{
    const long Index_1D = (long)v*N_k*N_j*N_i + (long)k*N_j*N_i + (long)j*N_i + (long)i;
    return Index_1D;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  Get_Value_from_Halo_IC_Data
// Description :  Get the target value at (x,y,z) by reading the IC_Data and doing linear interpolation
//
// Note        :  1.
//
// Parameter   :  x
//                y
//                z
//                v
//                index_halo
// Return      :  Interpolated_Value
//-------------------------------------------------------------------------------------------------------
double Get_Value_from_Halo_IC_Data( const double x, const double y, const double z, const int v, const int index_halo )
{
    // 1. IC information
    const int IC_Nx          = HaloMerger_Halo_IC_NCells[index_halo][0];
    const int IC_Ny          = HaloMerger_Halo_IC_NCells[index_halo][1];
    const int IC_Nz          = HaloMerger_Halo_IC_NCells[index_halo][2];
    const int IC_NVar        = ( HaloMerger_Halo_InitMode == 1 ) ? 2 : 1; // 1:(Real part & Imag part), 2:(Density)
    const int IC_FLOAT8      = HaloMerger_Halo_IC_FLOAT8[index_halo];
    size_t load_data_size    = ( IC_FLOAT8 ) ? sizeof(double) : sizeof(float);

    double IC_Range_EdgeL[3];
    double IC_dh[3];
    for (int d=0; d<3; d++) IC_Range_EdgeL[d] = HaloMerger_Halo_IC_Range_EdgeL[index_halo][d];
    for (int d=0; d<3; d++) IC_dh[d]          = HaloMerger_Halo_IC_dh[index_halo][d];

    // Use eight corner from the IC to do linear interpolation

    // 2. Index of the bottom-left corner in the IC
    const int IntCorner000_ICIndex[3] = {(int)floor( (x - IC_Range_EdgeL[0])/IC_dh[0] - 0.5 ),
                                         (int)floor( (y - IC_Range_EdgeL[1])/IC_dh[1] - 0.5 ),
                                         (int)floor( (z - IC_Range_EdgeL[2])/IC_dh[2] - 0.5 )};

    // 3. Physical coordinates of the eight corners
    double IntCorner_Coord[3][2];
    for (int d=0; d<3; d++)
    for (int IntCorner_ID=0; IntCorner_ID<2; IntCorner_ID++)
       IntCorner_Coord[d][IntCorner_ID] = IC_Range_EdgeL[d] + (IntCorner000_ICIndex[d] + IntCorner_ID + 0.5)*IC_dh[d];

    // 4. Values at the eight corners
    double IntCorner_Value[2][2][2];
    for (int IntCorner_ID_k=0; IntCorner_ID_k<2; IntCorner_ID_k++){ const int IntCorner_ICIndex_k = IntCorner000_ICIndex[2] + IntCorner_ID_k;
    for (int IntCorner_ID_j=0; IntCorner_ID_j<2; IntCorner_ID_j++){ const int IntCorner_ICIndex_j = IntCorner000_ICIndex[1] + IntCorner_ID_j;
    for (int IntCorner_ID_i=0; IntCorner_ID_i<2; IntCorner_ID_i++){ const int IntCorner_ICIndex_i = IntCorner000_ICIndex[0] + IntCorner_ID_i;

       if ( IntCorner_ICIndex_i < 0  ||  IntCorner_ICIndex_i >= IC_Nx  ||
            IntCorner_ICIndex_j < 0  ||  IntCorner_ICIndex_j >= IC_Ny  ||
            IntCorner_ICIndex_k < 0  ||  IntCorner_ICIndex_k >= IC_Nz )
       {
          // set the value as zero when the corner is outside the IC file
          IntCorner_Value[IntCorner_ID_k][IntCorner_ID_j][IntCorner_ID_i] = (double)0.0;
       }
       else
       {
          if ( IC_FLOAT8 ) IntCorner_Value[IntCorner_ID_k][IntCorner_ID_j][IntCorner_ID_i] = (double)(*((double*)&IC_Data[index_halo][IdxInIC(v, IntCorner_ICIndex_k, IntCorner_ICIndex_j, IntCorner_ICIndex_i, IC_NVar, IC_Nz, IC_Ny, IC_Nx)*load_data_size]));
          else             IntCorner_Value[IntCorner_ID_k][IntCorner_ID_j][IntCorner_ID_i] = (double)(*( (float*)&IC_Data[index_halo][IdxInIC(v, IntCorner_ICIndex_k, IntCorner_ICIndex_j, IntCorner_ICIndex_i, IC_NVar, IC_Nz, IC_Ny, IC_Nx)*load_data_size]));
       }

    }}} // IntCorner_ICIndex_k, IntCorner_ICIndex_j, IntCorner_ICIndex_i

    // 5. 3D linear interpolation
    double Interpolated_Value = Trilinear_Interpolation( x, y, z, IntCorner_Value, IntCorner_Coord[0], IntCorner_Coord[1], IntCorner_Coord[2] );

    return Interpolated_Value;

}
