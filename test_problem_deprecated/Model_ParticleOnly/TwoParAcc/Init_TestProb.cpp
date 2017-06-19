#include "GAMER.h"

#ifdef PARTICLE



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_TwoParAcc( real fluid[], const double x, const double y, const double z, const double Time );
static void Output_TwoParAcc( FILE *FileOut, const int Count );


// global variables in the two particles force test
// =======================================================================================
int  TwoParAcc_Mode;          // 1/2/3 : 3D test/1D test/1D self-force test
int  TwoParAcc_RSeed;         // random seed for setting particle position
int  TwoParAcc_NTest;         // number of tests
real TwoParAcc_MinR;          // minimum distance between two particles
real TwoParAcc_MaxR;          // maximum distance between two particles
int  TwoParAcc_LogSample;     // true/false --> evenly sample in log/linear scale
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the two particles force test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "Two particles force";

// check
# if ( MODEL != ELBDM )
# error : ERROR : "MODEL != ELBDM" in the two particles force test !!
# endif

# ifndef PARTICLE
# error : ERROR : "PARTICLE must be ON" in the two particles force test !!
# endif

# ifndef STORE_PAR_ACC
# error : ERROR : "STORE_PAR_ACC must be ON" in the two particles force test !!
# endif

# ifndef GRAVITY
# error : ERROR : "GRAVITY must be ON" in the two particles force test !!
# endif

# ifdef COMOVING
# error : ERROR : "COMOVING must be OFF" in the two particles force test !!
# endif

# ifndef SERIAL
# error : ERROR : "Only support SERIAL" in the two particles force test !!
# endif

# ifdef TIMING
# error : ERROR : "TIMING must be OFF" in the two particles force test !!
# endif

   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "pleaset set parameter %s to %d in the %s test !!\n",
                 "BC_POT_ISOLATED", BC_POT_ISOLATED, TestProb );

   if ( amr->Par->NPar_Active_AllRank != 2 )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "amr->Par->NPar_Active_AllRank", amr->Par->NPar_Active_AllRank, TestProb );

   if ( amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "amr->Par->Init", amr->Par->Init, TestProb );

   if ( OPT__INIT != INIT_STARTOVER )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "OPT__INIT", OPT__INIT, TestProb );

   if ( OPT__EXTERNAL_POT )
      Aux_Error( ERROR_INFO, "please disable the parameter %s in the %s test !!\n",
                 "OPT__EXTERNAL_POT", TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_TwoParAcc;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set the test problem parameters


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, " Note: test mode                                    = %d\n",     TwoParAcc_Mode  );
      Aux_Message( stdout, "       random seed for setting particle position    = %d\n",     TwoParAcc_RSeed );
      Aux_Message( stdout, "       number of tests                              = %d\n",     TwoParAcc_NTest );
      Aux_Message( stdout, "       minimum radius                               = %13.7e\n", TwoParAcc_MinR  );
      Aux_Message( stdout, "       maximum radius                               = %13.7e\n", TwoParAcc_MaxR  );
      Aux_Message( stdout, "       evenly sample in log scale                   = %s\n",     (TwoParAcc_LogSample)?
                                                                                             "Yes":"No"      );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// ------------------------------------------------------
// follow the same initialization procedure in Init_GAMER
// ------------------------------------------------------

// set the GPU ID and several GPU parameters
#  ifdef GPU
#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
   CUAPI_SetDevice( OPT__GPUID_SELECT );

   CUAPI_Set_Default_GPU_Parameter( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP );
#  endif


// initialize the array recording the previous physical time as an arbitrary "negative" number
   for (int lv=0; lv<NLEVEL; lv++)     Time_Prev[lv] = -9999.9;


// verify the input parameters
   Aux_Check_Parameter();


// initialize the timer function (TIMING is useless in this test)
   /*
#  ifdef TIMING
   Aux_CreateTimer();
#  endif
   */


// load the tables of the flag criteria from the input files "Input__FLAG_XXX"
   Init_Load_FlagCriteria();


// allocate memory for several global arrays
   Init_MemAllocate();


// take note
   Aux_TakeNote();

#  ifdef GPU
   CUAPI_DiagnoseDevice();
#  endif


// initialize the k-space Green's function for the isolated BC.
   if (  OPT__BC_POT == BC_POT_ISOLATED  &&
         ( OPT__GRAVITY_TYPE == GRAVITY_SELF || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  )
   Init_GreenFuncK();


// initialize the base level
   Init_BaseLevel();

   ELBDM_Init_StartOver_AssignData( 0 );


// set the random seed
   srand( TwoParAcc_RSeed );


// open the output file
   FILE *FileOut = fopen( "Record__TwoParAcc", "a" );

   fprintf( FileOut, "#Count : ID of test\n" );
   fprintf( FileOut, "#x/y/z : Particle positions or acceleration directions\n" );
   fprintf( FileOut, "#0/1   : Particle ID\n" );
   fprintf( FileOut, "#r     : Distance between the two particles\n" );
   fprintf( FileOut, "#Acc   : Gravitational acceleration on particles\n" );
   fprintf( FileOut, "#Num   : Numerical result\n" );
   fprintf( FileOut, "#Ana   : Analytical result\n" );
   fprintf( FileOut, "#Ratio : Num/Ana\n" );
   fprintf( FileOut, "#------------------------------------------------------------------------------\n" );

   fprintf( FileOut, "#%9s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n",
            "[1]Count", "[2]x0", "[3]y0", "[4]z0", "[5]x1", "[6]y1", "[7]z1", "[8]r",
             "[9]Acc_Num_x0", "[10]Acc_Num_x1", "[11]Acc_Ana_x",
            "[12]Acc_Num_y0", "[13]Acc_Num_y1", "[14]Acc_Ana_y",
            "[15]Acc_Num_z0", "[16]Acc_Num_z1", "[17]Acc_Ana_z",
            "[18]Acc_Num0", "[19]Acc_Num1", "[20]Acc_Ana", "[21]Acc_Ratio0", "[22]Acc_Ratio1" );


// start the test
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Start the two particles force test\n" );

   for (int t=0; t<TwoParAcc_NTest; t++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Test %10d ...\n", t );

//    initialize particles
      Par_Init_ByFunction();

      Par_FindHomePatch_Base( BaseP );


//    set up higher levels
      for (int lv=0; lv<NLEVEL-1; lv++)
      {
         Flag_Real( lv, USELB_NO );

         Flag_Buffer( lv );

         Init_Refine( lv );

         ELBDM_Init_StartOver_AssignData( lv + 1 );
      }


//    get the total number of patches at all ranks
      for (int lv=0; lv<NLEVEL; lv++)     Mis_GetTotalPatchNumber( lv );


      if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
      {
//       average density is useless in this test
//       if ( AveDensity_Init <= 0.0 )    Poi_GetAverageDensity();


//       evaluate the gravitational potential
         for (int lv=0; lv<NLEVEL; lv++)
            Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false );
      } // if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )


//    record the number of patches
      if ( OPT__PATCH_COUNT > 0 )   Aux_PatchCount();


//    perform some checks if necessary
      Aux_Check();


      const bool StoreAcc_Yes    = true;
      const bool UseStoredAcc_No = false;

      for (int lv=0; lv<=MAX_LEVEL; lv++)
         Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY,
                             StoreAcc_Yes, UseStoredAcc_No );


//    output the results
      Output_TwoParAcc( FileOut, t );


//    delete refinement levels for the next run
      for (int PID=0; PID<amr->NPatchComma[0][1]; PID++)
      {
         if ( amr->patch[0][0][PID]->NPar != 0 )
            amr->patch[0][0][PID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[0], true );
      }

      for (int PID0=0; PID0<amr->NPatchComma[1][1]; PID0+=8)
         amr->patch[0][0][ amr->patch[0][1][PID0]->father ]->son = -1;

      for (int lv=1; lv<NLEVEL; lv++)  amr->Lvdelete( lv, OPT__REUSE_MEMORY==2 );

   } // for (int t=0; t<TwoParAcc_NTest; t++)


// terminate the program
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Two particles force test is done\n" );

   fclose( FileOut );

   End_GAMER();

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_TestProbSol_TwoParAcc
// Description :  Initialize the background density field as zero for the two particles force test
//
// Note        :  1. Currently particle test must work with the ELBDM model
//                2. Invoked by "ELBDM_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Par_TestProbSol_TwoParAcc( real *fluid, const double x, const double y, const double z, const double Time )
{

// set wave function as zero everywhere
   fluid[REAL] = 0.0;
   fluid[IMAG] = 0.0;
   fluid[DENS] = 0.0;

} // FUNCTION : Par_TestProbSol_TwoParAcc



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTestProbParameter
// Description :  Load parameters for the test problem
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadTestProbParameter()
{

   const char FileName[] = "Input__TestProb";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;
   int    temp_int;


   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &TwoParAcc_Mode,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &TwoParAcc_RSeed,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &TwoParAcc_NTest,          string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &TwoParAcc_MinR,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &TwoParAcc_MaxR,           string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &TwoParAcc_MinR,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &TwoParAcc_MaxR,           string );
#  endif // #ifdef FLOAT8 ... else ...

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &TwoParAcc_LogSample,      string );

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default test problem parameters
   if ( TwoParAcc_NTest <= 0 )
   {
      switch( TwoParAcc_Mode )
      {
         case 1 : TwoParAcc_NTest = 200;   break;
         case 2 : TwoParAcc_NTest = 100;   break;
         case 3 : TwoParAcc_NTest = 200;   break;
         Aux_Error( ERROR_INFO, "unsupported TwoParAcc_Mode (%d) !!\n", TwoParAcc_Mode );
      }

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "TwoParAcc_NTest", TwoParAcc_NTest );
   }

   if ( TwoParAcc_MinR <= 0.0 )
   {
      TwoParAcc_MinR = 0.1*amr->dh[0];

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "TwoParAcc_MinR", TwoParAcc_MinR );
   }

   if ( TwoParAcc_MaxR <= 0.0 )
   {
//    since the outmost base patches are not allowed to be refined for the isolated gravity solver,
//    we set MaxR = 0.5*L-3*PatchSize to avoid putting the second particle on that region
//    TwoParAcc_MaxR = 0.5*amr->BoxSize[0] - 3.0*PS1*amr->dh[0];

//    1.5 base level patch
      TwoParAcc_MaxR = 1.5*PS1*amr->dh[0];

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "TwoParAcc_MaxR", TwoParAcc_MaxR );
   }

   if ( TwoParAcc_LogSample < 0 )
   {
      switch( TwoParAcc_Mode )
      {
         case 1 : TwoParAcc_LogSample = 1;   break;
         case 2 : TwoParAcc_LogSample = 1;   break;
         case 3 : TwoParAcc_LogSample = 0;   break;
         Aux_Error( ERROR_INFO, "unsupported TwoParAcc_Mode (%d) !!\n", TwoParAcc_Mode );
      }

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "TwoParAcc_LogSample", TwoParAcc_LogSample );
   }

// check
   if ( TwoParAcc_Mode < 1  ||  TwoParAcc_Mode > 3 )
      Aux_Error( ERROR_INFO, "TwoParAcc_Mode (%d) must be 1/2/3 !!\n", TwoParAcc_Mode );

   if ( TwoParAcc_RSeed < 0 )    Aux_Error( ERROR_INFO, "TwoParAcc_RSeed (%d) < 0 !!\n", TwoParAcc_RSeed );
   if ( TwoParAcc_NTest < 0 )    Aux_Error( ERROR_INFO, "TwoParAcc_NTest (%d) < 0 !!\n", TwoParAcc_NTest );
   if ( TwoParAcc_MinR < 0.0 )   Aux_Error( ERROR_INFO, "TwoParAcc_MinR (%14.7e) < 0.0 !!\n", TwoParAcc_MinR );
   if ( TwoParAcc_MaxR < 0.0 )   Aux_Error( ERROR_INFO, "TwoParAcc_MaxR (%14.7e) < 0.0 !!\n", TwoParAcc_MaxR );

} // FUNCTION : LoadTestProbParameter



//------------------------------------------------------------------------------------------------------
// Function    :  Output_TwoParAcc
// Description :  Output the particle position
//
// Parameter   :  FileOut : Pointer to the file variable
//                Count   : ID of test
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_TwoParAcc( FILE *FileOut, const int Count )
{

   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Acc[3] = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };

   real dr[3], r, Acc_Num0[3], Acc_Num1[3], Acc_Ana[3], Acc_Num_Abs0, Acc_Num_Abs1;
   real Acc_Ana_Abs, Acc_Ratio0, Acc_Ratio1;


// calculate the analytical acceleration
   for (int d=0; d<3; d++)    dr[d] = Pos[d][1] - Pos[d][0];

   r = SQRT( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );

   for (int d=0; d<3; d++)    Acc_Ana[d] = NEWTON_G*Mass[1]*dr[d]/CUBE(r);

   Acc_Ana_Abs = SQRT( SQR(Acc_Ana[0]) + SQR(Acc_Ana[1]) + SQR(Acc_Ana[2]) );

   for (int d=0; d<3; d++)
   {
      Acc_Num0[d] = Acc[d][0];
      Acc_Num1[d] = Acc[d][1];
   }

   Acc_Num_Abs0 = SQRT( SQR(Acc_Num0[0]) + SQR(Acc_Num0[1]) + SQR(Acc_Num0[2]) );
   Acc_Num_Abs1 = SQRT( SQR(Acc_Num1[0]) + SQR(Acc_Num1[1]) + SQR(Acc_Num1[2]) );
   Acc_Ratio0   = Acc_Num_Abs0 / Acc_Ana_Abs;
   Acc_Ratio1   = Acc_Num_Abs1 / Acc_Ana_Abs;


// output
   fprintf( FileOut, "%10d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e",
            Count, Pos[0][0], Pos[1][0], Pos[2][0], Pos[0][1], Pos[1][1], Pos[2][1], r,
            Acc_Num0[0], Acc_Num1[0], Acc_Ana[0] );

   fprintf( FileOut, " %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
            Acc_Num0[1], Acc_Num1[1], Acc_Ana[1], Acc_Num0[2], Acc_Num1[2], Acc_Ana[2],
            Acc_Num_Abs0, Acc_Num_Abs1, Acc_Ana_Abs, Acc_Ratio0, Acc_Ratio1 );

} // FUNCTION : Output_TwoParAcc



#endif // #ifdef PARTICLE
