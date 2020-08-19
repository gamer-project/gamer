#include "GAMER.h"
#include "TestProb.h"

#if  ( MODEL == SR_HYDRO )

void Cartesian2Spherical( const double x[], double r[] );
void Spherical2Cartesian( const double r[], double x[] );
void CartesianRotate( double x[], double theta, double phi, bool inverse );


//  X²    Y²    Z²
// --- + --- + --- = 1
//  a²    b²    c²
// 
//
// problem-specific global variables
// =======================================================================================
static double  Blast_Dens_Bg;                      // background proper mass density
static double  Blast_Temp_Bg;                      // background temperature
static double  Blast_Dens_Src;                     // source density
static double  Blast_Temp_Src;                     // source temperature
static double  Blast_Radius_x;                     // length of source axis along x-axis (a)
static double  Blast_Radius_y;                     // length of source axis along y-axis (b)
static double  Blast_Radius_z;                     // length of source axis along z-axis (c)
static double  RotatedAngle[2];                    // angle of rotation in unit of degree ([0]: theta, [1]: phi)
// the angles `theta' and `phi' are defined as below:
// --> https://upload.wikimedia.org/wikipedia/commons/4/4f/3D_Spherical.svg

static double  Blast_Center[3];                    // coordinates of explosion center ([0]:x, [1]:y, [2]:z)

static int     NumSource;                          // number of blast wave sources
static char    Random_File[MAX_STRING];            // the name of file that contain information of source (e.g., position, angle of rotation)
static double *Random_Data               = NULL;   // pointer that point to file
static int     Random_Table_NBin;                  // number of columns in `Random_File`
static double *Table_x, *Table_y, *Table_z, *Table_theta,  *Table_phi, *Table_Temp;
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


#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != SR_HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( !OPT__INIT_RESTRICT )
      Aux_Error( ERROR_INFO, "OPT__INIT_RESTRICT must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == SR_HYDRO )
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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Blast_Dens_Bg",     &Blast_Dens_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Temp_Bg",     &Blast_Temp_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Dens_Src",    &Blast_Dens_Src,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Temp_Src",    &Blast_Temp_Src,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Radius_x",    &Blast_Radius_x,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Radius_y",    &Blast_Radius_y,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Radius_z",    &Blast_Radius_z,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Center_X",    &Blast_Center[0],       -1.0,          NoMin_double,     amr->BoxSize[0]   );
   ReadPara->Add( "Blast_Center_Y",    &Blast_Center[1],       -1.0,          NoMin_double,     amr->BoxSize[1]   );
   ReadPara->Add( "Blast_Center_Z",    &Blast_Center[2],       -1.0,          NoMin_double,     amr->BoxSize[2]   );
   ReadPara->Add( "NumSource",         &NumSource,              1,            1,                NoMax_int         );
   ReadPara->Add( "Rotated_Theta",     &RotatedAngle[0],        0.0,          0.0,              180.0             );
   ReadPara->Add( "Rotated_Phi",       &RotatedAngle[1],        0.0,          0.0,              360.0             );
   ReadPara->Add( "Random_File",       Random_File,             Useless_str,  Useless_str,      Useless_str       );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Blast_Center[d] < 0.0 )  Blast_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 0.5*BOX_SIZE;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                         = %d\n",     TESTPROB_ID      );
      Aux_Message( stdout, "  background mass density                 = %13.7e\n", Blast_Dens_Bg    );
      Aux_Message( stdout, "  background pressure                     = %13.7e\n", Blast_Temp_Bg    );
      Aux_Message( stdout, "  density ratio of center to background   = %13.7e\n", Blast_Dens_Src   );
      Aux_Message( stdout, "  pressure ratio of center to background  = %13.7e\n", Blast_Temp_Src   );
      Aux_Message( stdout, "  explosion radius a                      = %13.7e\n", Blast_Radius_x   );
      Aux_Message( stdout, "  explosion radius b                      = %13.7e\n", Blast_Radius_y   );
      Aux_Message( stdout, "  explosion radius c                      = %13.7e\n", Blast_Radius_z   );
      Aux_Message( stdout, "  explosion center x                      = %13.7e\n", Blast_Center[0]  );
      Aux_Message( stdout, "  explosion center y                      = %13.7e\n", Blast_Center[1]  );
      Aux_Message( stdout, "  explosion center z                      = %13.7e\n", Blast_Center[2]  );
      Aux_Message( stdout, "  rotated angle theta                     = %13.7e\n", RotatedAngle[0]  );
      Aux_Message( stdout, "  rotated angle phi                       = %13.7e\n", RotatedAngle[1]  );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );


   const bool RowMajor_No  = false;                // load data into the column-major order
   const bool AllocMem_Yes = true;                 // allocate memory for Merger_Prof1/2
   const int  NCol         = 6;                    // total number of columns to load
   const int  Col[NCol]    = {1, 2, 3, 4, 5, 6};   // target columns: (x, y, z, radius, density, temperature)
   
   if ( NumSource > 1 )   
     Random_Table_NBin = Aux_LoadTable( Random_Data, Random_File, NCol, Col, RowMajor_No, AllocMem_Yes );
              
   if ( Random_Data == NULL && NumSource > 1 )
    Aux_Error( ERROR_INFO, "Random_Data == NULL !!\n" );



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
   double Prim_Bg[5], Cons_Bg[5];
   double Prim_Src[5], Cons_Src[5];
   double RotatedCartesian[3];

   double Theta = RotatedAngle[0]*M_PI/180.0;
   double Phi   = RotatedAngle[1]*M_PI/180.0;

   Prim_Bg[0]  = Blast_Dens_Bg;
   Prim_Bg[1]  = 0.0;
   Prim_Bg[2]  = 0.0;
   Prim_Bg[3]  = 0.0;
   Prim_Bg[4]  = Blast_Temp_Bg*Blast_Dens_Bg;

   Prim_Src[0] = Blast_Dens_Src;
   Prim_Src[1] = 0.0;
   Prim_Src[2] = 0.0; 
   Prim_Src[3] = 0.0; 
   Prim_Src[4] = Blast_Dens_Src*Blast_Temp_Src;

   SRHydro_Pri2Con (Prim_Bg,  Cons_Bg,  GAMMA);


   bool InsideEllipsoid = false; 
// --> true/false: the current cell is inside/outside source

   for ( int i = 0; i < NumSource;i++ )
   {
	  if ( NumSource > 1 )
	  {
        Table_x          = Random_Data + 0*Random_Table_NBin;
        Table_y          = Random_Data + 1*Random_Table_NBin;
        Table_z          = Random_Data + 2*Random_Table_NBin;
        Table_theta      = Random_Data + 3*Random_Table_NBin;
        Table_phi        = Random_Data + 4*Random_Table_NBin;
        Table_Temp       = Random_Data + 5*Random_Table_NBin;

		Blast_Center[0]  = Table_x[i];
		Blast_Center[1]  = Table_y[i];
		Blast_Center[2]  = Table_z[i];
		Theta            = Table_theta[i];
	    Phi              = Table_phi[i];
		Prim_Src    [4]  = Blast_Dens_Src*Table_Temp[i];
	  }


//    1. rotate source
      RotatedCartesian[0] = x - Blast_Center[0];
      RotatedCartesian[1] = y - Blast_Center[1];
      RotatedCartesian[2] = z - Blast_Center[2];

//    2. rotate source
      CartesianRotate( RotatedCartesian, Theta, Phi, false );

//    3. check the current cell is inside source or not
      InsideEllipsoid  = SQR ( RotatedCartesian[0] / Blast_Radius_x ) 
                       + SQR ( RotatedCartesian[1] / Blast_Radius_y ) 
                       + SQR ( RotatedCartesian[2] / Blast_Radius_z ) < 1.0;

      if ( InsideEllipsoid == true )
      {
         SRHydro_Pri2Con (Prim_Src, Cons_Src, GAMMA);

         fluid[DENS] = (real) Cons_Src[0];
         fluid[MOMX] = (real) Cons_Src[1];
         fluid[MOMY] = (real) Cons_Src[2];
         fluid[MOMZ] = (real) Cons_Src[3];
         fluid[ENGY] = (real) Cons_Src[4];

         return;
      }
   }

   fluid[DENS] = (real) Cons_Bg[0];
   fluid[MOMX] = (real) Cons_Bg[1];
   fluid[MOMY] = (real) Cons_Bg[2];
   fluid[MOMZ] = (real) Cons_Bg[3];
   fluid[ENGY] = (real) Cons_Bg[4];

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == SR_HYDRO )




bool Flag_User_Blast( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{
   bool Flag = false;
   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
 
   const double dR[3]  = { Pos[0]-Blast_Center[0], Pos[1]-Blast_Center[1], Pos[2]-Blast_Center[2] };

   const double R      = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );
   

// use nested refinement to ensure tiny source can be properly refined at the beginning
// --> 1.8 is a approximated number of sqrt(3), sqrt(3)*dh is the diaganal length of cell at level `lv'
   Flag |= R < dh*1.8;
       
   if ( Flag ) return true;
   else        return false;
}


void End_Blast()
{
	  if ( NumSource > 1 )  delete [] Random_Data;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_BlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_BlastWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == SR_HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_User_Ptr            = Flag_User_Blast;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = End_Blast;
#  endif // #if ( MODEL == SR_HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BlastWave

#endif
