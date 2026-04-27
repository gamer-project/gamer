#include "GAMER.h"

static void AddNewField_ELBDM_UniformGranule( void );
static void Init_User_ELBDM_UniformGranule( void );
static void Do_CF( void );
static void End_UniformGranule( void );

void Aux_ComputeCorrelation( Profile_t *Correlation[], const Profile_t *prof_init[], const double Center[],
                             const double r_max_input, const double dr_min, const bool LogBin, const double LogBinRatio,
                             const bool RemoveEmpty, const long TVarBitIdx[], const int NProf, const int MinLv, const int MaxLv,
                             const PatchType_t PatchType, const double PrepTime, const double dr_min_prof );

// problem-specific global variables
// =======================================================================================
static FieldIdx_t Idx_Dens0 = Idx_Undefined;  // field index for storing the **initial** density
static bool     ComputeCorrelation;           // flag for computing correlation
static bool     ReComputeCorrelation;         // flag for recomputing correlation during restart; use the simulation time of RESTART as the initial time for computing time correlation; only available for RESTART
static char     FilePath_corr[MAX_STRING];    // folder for storing correlation function text files
static int      MinLv_corr;                   // compute correlation function from AMR levels MinLv to MaxLv
static int      MaxLv_corr;                   // compute correlation function from AMR levels MinLv to MaxLv
static int      StepInitial_corr;             // inital step for recording correlation function
static int      StepInterval_corr;            // interval for recording correlation function
static int      StepEnd_corr;                 // end step for recording correlation function
static double   Center_corr[3];               // center for computing correlation function (currently set to the center of mass)
static double   RadiusMax_prof;               // maximum radius for correlation function statistics (profile)
static double   dr_min_prof;                  // bin size of correlation function statistics (minimum size for logarithmic bin) (profile)
static bool     LogBin_prof;                  // logarithmic bin or not (profile)
static double   LogBinRatio_prof;             // ratio of adjacent log bins for logarithmic bin (profile)
static bool     RemoveEmpty_prof;             // remove bins without any samples; false: Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0 (profile)
static double   RadiusMax_corr;               // maximum radius for correlation function statistics (correlation)
static double   dr_min_corr;                  // bin size of correlation function statistics (minimum size for logarithmic bin) (correlation)
static bool     LogBin_corr;                  // logarithmic bin or not (correlation)
static double   LogBinRatio_corr;             // ratio of adjacent log bins for logarithmic bin (correlation)
static bool     RemoveEmpty_corr;             // remove bins without any samples; false: Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0 (correlation)

static Profile_t *Prof_Dens_initial = new Profile_t(); // pointer to save initial density profile
static Profile_t *Correlation_Dens  = new Profile_t(); // pointer to save density correlation function
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

#  if ( NCOMP_PASSIVE_USER != 1 )
   Aux_Error( ERROR_INFO, "must set NCOMP_PASSIVE_USER to 1 !!\n" );
#  endif

// only accept OPT__INIT == INIT_BY_RESTART or OPT__INIT == INIT_BY_FILE
   if ( OPT__INIT != INIT_BY_RESTART  &&  OPT__INIT != INIT_BY_FILE )
      Aux_Error( ERROR_INFO, "only support OPT__INIT == INIT_BY_RESTART or OPT__INIT == INIT_BY_FILE !!\n" );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

// check whether fluid boundary condition in Input__Parameter is set properly
   for ( int direction = 0; direction < 6; direction++ )
   {
      if ( OPT__BC_FLU[direction] != BC_FLU_PERIODIC )
         Aux_Error( ERROR_INFO, "must set periodic BC for fluid --> reset OPT__BC_FLU[%d] to 1 !!\n", direction );
   }
   if ( OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must set periodic BC for gravity --> reset OPT__BC_POT to 1 !!\n" );

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
//                3. Must NOT call any EoS routine here since it hasn't been initialized at this point
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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",          &VARIABLE,                DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "ComputeCorrelation",       &ComputeCorrelation,      false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "ReComputeCorrelation",     &ReComputeCorrelation,    false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "FilePath_corr",            FilePath_corr,            Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "MinLv_corr",               &MinLv_corr,              0,             0,                MAX_LEVEL         );
   ReadPara->Add( "MaxLv_corr",               &MaxLv_corr,              MAX_LEVEL,     0,                MAX_LEVEL         );
   ReadPara->Add( "StepInitial_corr",         &StepInitial_corr,        0,             0,                NoMax_int         );
   ReadPara->Add( "StepInterval_corr",        &StepInterval_corr,       1,             1,                NoMax_int         );
   ReadPara->Add( "StepEnd_corr",             &StepEnd_corr,            NoMax_int,     0,                NoMax_int         );
   ReadPara->Add( "RadiusMax_corr",           &RadiusMax_corr,          Eps_double,    NoMin_double,     NoMax_double      );
   ReadPara->Add( "dr_min_corr",              &dr_min_corr,             Eps_double,    NoMin_double,     NoMax_double      );
   ReadPara->Add( "LogBin_corr",              &LogBin_corr,             false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "LogBinRatio_corr",         &LogBinRatio_corr,        1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "RemoveEmpty_corr",         &RemoveEmpty_corr,        false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "dr_min_prof",              &dr_min_prof,             Eps_double,    NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if (ComputeCorrelation)
   {
       if ( !OPT__RECORD_USER )
          Aux_Error( ERROR_INFO, "OPT__RECORD_USER should be turned on to enable ComputeCorrelation !!\n" );

       if ( dr_min_corr <=Eps_double )          dr_min_corr = 1e-2*0.5*amr->BoxSize[0];
       if ( RadiusMax_corr<=Eps_double )        RadiusMax_corr = 0.5*amr->BoxSize[0];
       if ( LogBinRatio_corr<=1.0 )             LogBinRatio_corr = 2.0;

       if ( dr_min_prof <=Eps_double )          dr_min_prof = 0.25*dr_min_corr;
       RadiusMax_prof                           = RadiusMax_corr * 1.05;   // assigned by Test Problem
       LogBinRatio_prof                         = 1.0;                     // hard-coded by Test Problem (no effect)
       LogBin_prof                              = false;                   // hard-coded by Test Problem
       RemoveEmpty_prof                         = false;                   // hard-coded by Test Problem

       if ( MinLv_corr < 0 ) MinLv_corr = 0;
       if ( MaxLv_corr < MinLv_corr ) MaxLv_corr = MAX_LEVEL;
       if ( FilePath_corr == "\0" )  sprintf( FilePath_corr, "./" );
       else
       {
          FILE *file_checker = fopen(FilePath_corr, "r");
          if (!file_checker)
             Aux_Error( ERROR_INFO, "File path %s for saving correlation function text files does not exist!! Please create!!\n", FilePath_corr );
          else
             fclose(file_checker);
       }
   }

// (1-3) check the runtime parameters

// (2) set the problem-specific derived parameters

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = 50000;
   const double End_T_Default    = 2.0e-3;

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
      Aux_Message(    stdout, "==============================================================================\n" );
      Aux_Message(    stdout, "  test problem ID                             = %d\n"    , TESTPROB_ID            );
      Aux_Message(    stdout, "  compute correlation                         = %d\n"    , ComputeCorrelation     );
      if (ComputeCorrelation)
      {
         Aux_Message( stdout, "  radius maximum      (correlation)           = %13.7e\n", RadiusMax_corr         );
         Aux_Message( stdout, "  histogram bin size  (correlation)           = %13.7e\n", dr_min_corr            );
         Aux_Message( stdout, "  use logarithmic bin (correlation)           = %d\n"    , LogBin_corr            );
         Aux_Message( stdout, "  log bin ratio       (correlation)           = %13.7e\n", LogBinRatio_corr       );
         Aux_Message( stdout, "  remove empty bin    (correlation)           = %d\n"    , RemoveEmpty_corr       );
         Aux_Message( stdout, "  radius maximum      (profile, assigned)     = %13.7e\n", RadiusMax_prof         );
         Aux_Message( stdout, "  histogram bin size  (profile)               = %13.7e\n", dr_min_prof            );
         Aux_Message( stdout, "  use logarithmic bin (profile, hard-coded)   = %d\n"    , LogBin_prof            );
         Aux_Message( stdout, "  log bin ratio       (profile, no effect)    = %13.7e\n", LogBinRatio_prof       );
         Aux_Message( stdout, "  remove empty bin    (profile, hard-coded)   = %d\n"    , RemoveEmpty_prof       );
         Aux_Message( stdout, "  minimum level                               = %d\n"    , MinLv_corr             );
         Aux_Message( stdout, "  maximum level                               = %d\n"    , MaxLv_corr             );
         Aux_Message( stdout, "  folder for storing correlation text file    = %s\n"    , FilePath_corr          );
         if ( OPT__INIT == INIT_BY_RESTART )
         Aux_Message( stdout, "  re-compute correlation using restart time   = %d\n"    , ReComputeCorrelation   );
      }
      Aux_Message(    stdout, "==============================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_ELBDM_UniformGranule
// Description :  Store the initial density as Dens0
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void AddNewField_ELBDM_UniformGranule( void )
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   Idx_Dens0 = AddField( "Dens0", FIXUP_FLUX_NO, FIXUP_REST_NO, FLOOR_NO, NORMALIZE_NO, INTERP_FRAC_NO );
#  endif

} // FUNCTION : AddNewField_ELBDM_UniformGranule

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_ELBDM_UniformGranule
// Description :  Store the initial density
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Init_User_ELBDM_UniformGranule( void )
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
//    store the initial density in both Sg so that we don't have to worry about which Sg to be used
//    a. for restart and ReComputeCorrelation disabled, the initial density has already been loaded and we just need to copy the data to another Sg
      if ( ( OPT__INIT == INIT_BY_RESTART ) && ( !ReComputeCorrelation ) )
      {
         const real Dens0 = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i];

         amr->patch[ 1-amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
      }

//    b. for starting a new simulation, we must copy the initial density to both Sg
      else
      {
         const real Dens0 = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];

         amr->patch[   amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
         amr->patch[ 1-amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
      }
   }

   if (ComputeCorrelation)
   {
      const double InitialTime = Time[0];
      if ( MPI_Rank == 0 )  Aux_Message( stdout, "StepInitial_corr = %d ; StepInterval_corr = %d ; StepEnd_corr = %d\n", StepInitial_corr, StepInterval_corr, StepEnd_corr);

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "InitialTime = %13.6e \n", InitialTime );

//    determine the central position for computing correlation function
      if ( MPI_Rank == 0 )  Aux_Message( stdout, "Calculate the initial center-of-mass position:\n");

      double FinaldR;
      int    FinalNIter;
      double CoM_ref[3];
      double MaxR    = __FLT_MAX__;
      double MinWD   = 0.0;
      double TolErrR = amr->dh[MAX_LEVEL];
      int    MaxIter = 10;

      Extrema_t Max_Dens;
      Max_Dens.Field     = _DENS;
      Max_Dens.Radius    = __FLT_MAX__; // entire domain
      Max_Dens.Center[0] = amr->BoxCenter[0];
      Max_Dens.Center[1] = amr->BoxCenter[1];
      Max_Dens.Center[2] = amr->BoxCenter[2];

      Aux_FindExtrema( &Max_Dens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

      for (int d=0; d<3; d++) CoM_ref[d] = Max_Dens.Coord[d];

      Aux_FindWeightedAverageCenter( Center_corr, CoM_ref, MaxR, MinWD, _TOTAL_DENS, TolErrR, MaxIter, &FinaldR, &FinalNIter );

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "Initial center-of-mass position is ( %14.11e,%14.11e,%14.11e )\n", Center_corr[0], Center_corr[1], Center_corr[2] );

//    commpute density profile for passive field
      if ( MPI_Rank == 0 )  Aux_Message( stdout, "Calculate initial density profile:\n");

      const long TVar[] = {BIDX(Idx_Dens0)};
      Aux_ComputeProfile( &Prof_Dens_initial, Center_corr, RadiusMax_prof, dr_min_prof, LogBin_prof, LogBinRatio_prof, RemoveEmpty_prof, TVar, 1, MinLv_corr, MaxLv_corr, PATCH_LEAF, InitialTime, true );

      if ( MPI_Rank == 0 )
      {
         char Filename[MAX_STRING];
         sprintf( Filename, "%s/initial_profile.txt", FilePath_corr );
         FILE *output_initial_prof = fopen(Filename, "w");
         fprintf( output_initial_prof, "#%19s  %21s  %21s  %21s  %11s\n", "Radius", "Dens", "Dens_Sigma" , "Weighting", "Cell_Number");
         for (int b=0; b<Prof_Dens_initial->NBin; b++)
            fprintf( output_initial_prof, "%20.14e  %21.14e  %21.14e  %21.14e  %11ld\n",
                     Prof_Dens_initial->Radius[b], Prof_Dens_initial->Data[b], Prof_Dens_initial->Data_Sigma[b], Prof_Dens_initial->Weight[b], Prof_Dens_initial->NCell[b] );
         fclose(output_initial_prof);
      }
   }
#  endif // #if ( NCOMP_PASSIVE_USER > 0 )

} // FUNCTION : Init_User_ELBDM_UniformGranule

#endif // #if ( MODEL == ELBDM  && defined GRAVITY )

//-------------------------------------------------------------------------------------------------------
// Function    :  Do_CF
// Description :  Compute the correlation function
//
// Note        :  1. Linked to the function pointer "Aux_Record_User_Ptr"
//                2. Please turn on the runtime option "OPT__RECORD_USER"
//-------------------------------------------------------------------------------------------------------
static void Do_CF( void )
{
// Compute correlation if ComputeCorrelation flag is true
   if (ComputeCorrelation)
   {
      if ( (Step>=StepInitial_corr) && (((Step-StepInitial_corr)%StepInterval_corr)==0) && (Step<=StepEnd_corr) )
      {
         const long TVar[] = {_DENS};
         Aux_ComputeCorrelation( &Correlation_Dens, (const Profile_t**) &Prof_Dens_initial, Center_corr, RadiusMax_corr, dr_min_corr, LogBin_corr, LogBinRatio_corr,
                                 RemoveEmpty_corr, TVar, 1, MinLv_corr, MaxLv_corr, PATCH_LEAF, Time[0], dr_min_prof);

         if ( MPI_Rank == 0 )
         {
            char Filename[MAX_STRING];
            sprintf( Filename, "%s/correlation_function_t=%.4e.txt", FilePath_corr, Time[0] );
            FILE *output_correlation = fopen(Filename, "w");
            fprintf( output_correlation, "#%19s  %21s  %21s  %11s\n", "Radius", "Correlation_Function", "Weighting", "Cell_Number");
            for (int b=0; b<Correlation_Dens->NBin; b++)
                fprintf( output_correlation, "%20.14e  %21.14e  %21.14e  %11ld\n",
                         Correlation_Dens->Radius[b], Correlation_Dens->Data[b], Correlation_Dens->Weight[b], Correlation_Dens->NCell[b] );
            fclose(output_correlation);
         }
      }
   }  // if (ComputeCorrelation)
} // FUNCTION : Do_CF

//-------------------------------------------------------------------------------------------------------
// Function    :  End_UniformGranule
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//-------------------------------------------------------------------------------------------------------
void End_UniformGranule()
{

   delete Prof_Dens_initial;
   delete Correlation_Dens;
   Prof_Dens_initial = NULL;
   Correlation_Dens  = NULL;

} // FUNCTION : End_UniformGranule

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_UniformGranule
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_UniformGranule()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();

   Init_Field_User_Ptr = AddNewField_ELBDM_UniformGranule;
   Init_User_Ptr       = Init_User_ELBDM_UniformGranule;
   Aux_Record_User_Ptr = Do_CF;
   End_User_Ptr        = End_UniformGranule;
#  endif // #if ( MODEL == ELBDM )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_UniformGranule
