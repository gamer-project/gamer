#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_Parameter
// Description :  Load the initial values of simulation parameters from the file "Input__Parameter"
//-------------------------------------------------------------------------------------------------------
void Init_Load_Parameter() 
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_Parameter ... \n" );


   const char FileName[] = "Input__Parameter";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   int    temp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;


// simulation scale
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BOX_SIZE,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NX0_TOT[0],               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NX0_TOT[1],               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NX0_TOT[2],               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank_X[0],           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank_X[1],           string );
   
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank_X[2],           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OMP_NTHREAD,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &END_T,                    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%ld%s",  &END_STEP,                 string );

   getline( &input_line, &len, File );


// boundary condition
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_FLU[0] = (OptFluBC_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_FLU[1] = (OptFluBC_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_FLU[2] = (OptFluBC_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_FLU[3] = (OptFluBC_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_FLU[4] = (OptFluBC_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_FLU[5] = (OptFluBC_t)temp_int;

#  ifdef GRAVITY
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__BC_POT = (OptPotBC_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &GFUNC_COEFF0,             string );

#  else // #ifdef GRAVITY
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif // #ifdef GRAVITY ... else ...

   getline( &input_line, &len, File );


// particle
#  ifdef PARTICLE
   getline( &input_line, &len, File );
   sscanf( input_line, "%ld%s",  &amr->Par->NPar,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &amr->Par->Init,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &amr->Par->Interp,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &amr->Par->Integ,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   amr->Par->ImproveAcc = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   amr->Par->PredictPos = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &amr->Par->RemoveCell,     string );

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif // #ifdef PARTICLE ... else ...

   getline( &input_line, &len, File );


// cosmology simulations (COMOVING frame)
#  ifdef COMOVING
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &A_INIT,                   string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OMEGA_M0,                 string );

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif

   getline( &input_line, &len, File );


// time-step determination
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &DT__FLUID,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &DT__FLUID_INIT,           string );

   getline( &input_line, &len, File );
#  ifdef GRAVITY
   sscanf( input_line, "%lf%s",  &DT__GRAVITY,              string );
#  endif
   
   getline( &input_line, &len, File );
#  if ( MODEL == ELBDM )
   sscanf( input_line, "%lf%s",  &DT__PHASE,                string );
#  endif

   getline( &input_line, &len, File );
#  ifdef PARTICLE 
   sscanf( input_line, "%lf%s",  &DT__PARVEL,               string );
#  endif

   getline( &input_line, &len, File );
#  ifdef PARTICLE 
   sscanf( input_line, "%lf%s",  &DT__PARVEL_MAX,           string );
#  endif

   getline( &input_line, &len, File );
#  ifdef COMOVING
   sscanf( input_line, "%lf%s",  &DT__MAX_DELTA_A,          string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__ADAPTIVE_DT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RECORD_DT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__DT_USER = (bool)temp_int;

   getline( &input_line, &len, File );


// domain refinement
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &REGRID_COUNT,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &FLAG_BUFFER_SIZE,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MAX_LEVEL,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_RHO = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_RHO_GRADIENT = (bool)temp_int;

   getline( &input_line, &len, File );
#  if   ( MODEL == HYDRO ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_PRES_GRADIENT = (bool)temp_int;
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif

   getline( &input_line, &len, File );
#  if ( MODEL == ELBDM ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_ENGY_DENSITY = (bool)temp_int;
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_LOHNER_DENS = (bool)temp_int;

   getline( &input_line, &len, File );
#  if   ( MODEL == HYDRO ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_LOHNER_ENGY = (bool)temp_int;
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif

   getline( &input_line, &len, File );
#  if   ( MODEL == HYDRO ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_LOHNER_PRES = (bool)temp_int;
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_LOHNER_FORM = (OptLohnerForm_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_USER = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_REGION = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__PATCH_COUNT,         string );

   getline( &input_line, &len, File );
#  ifdef PARTICLE
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__PAR_LEVEL = (bool)temp_int;
#  endif

   getline( &input_line, &len, File );


// load balance
#  ifdef LOAD_BALANCE
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &LB_INPUT__WLI_MAX,        string );

#  else // #ifdef LOAD_BALANCE ... else ...
   getline( &input_line, &len, File );
#  endif // #ifdef LOAD_BALANCE ... else ...

   getline( &input_line, &len, File );


// fluid solvers in HYDRO
#  if   ( MODEL == HYDRO )
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &GAMMA,                    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &MINMOD_COEFF,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &EP_COEFF,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__LR_LIMITER = (LR_Limiter_t)temp_int;

   getline( &input_line, &len, File ); // skip two comment lines
   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__WAF_LIMITER = (WAF_Limiter_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CORR_UNPHY_SCHEME = (OptRSolver_t)temp_int;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif // MODEL

   getline( &input_line, &len, File );


// ELBDM solvers
#  if ( MODEL == ELBDM )
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ELBDM_MASS,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ELBDM_PLANCK_CONST,       string );

   getline( &input_line, &len, File );
#  ifdef QUARTIC_SELF_INTERACTION
   sscanf( input_line, "%lf%s",  &ELBDM_LAMBDA,             string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ELBDM_TAYLOR3_COEFF,      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   ELBDM_TAYLOR3_AUTO = (bool)temp_int;

#  else // #if ( MODEL == ELBDM ) ... else ...
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif // #if ( MODEL == ELBDM ) ... else ...

   getline( &input_line, &len, File );


// fluid solvers in both HYDRO/MHD/ELBDM
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &FLU_GPU_NPGROUP,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &GPU_NSTREAM,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FIXUP_FLUX = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FIXUP_RESTRICT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OVERLAP_MPI = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RESET_FLUID = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CORR_UNPHY = (bool)temp_int;

   getline( &input_line, &len, File );


// self-gravity
#  ifdef GRAVITY

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &NEWTON_G,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SOR_OMEGA,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &SOR_MAX_ITER,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &SOR_MIN_ITER,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MG_MAX_ITER,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MG_NPRE_SMOOTH,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MG_NPOST_SMOOTH,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &MG_TOLERATED_ERROR,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &POT_GPU_NPGROUP,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__GRA_P5_GRADIENT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__GRAVITY_TYPE = (OptGravityType_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__EXTERNAL_POT = (bool)temp_int;

#  else // #ifdef GRAVITY ... else ...

   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );

#  endif // #ifdef GRAVITY ... else ...

   getline( &input_line, &len, File );


// initialization
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INIT = (OptInit_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RESTART_HEADER = (OptRestartH_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__UM_START_LEVEL,      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__UM_START_NVAR,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__UM_START_DOWNGRADE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__UM_START_REFINE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__UM_FACTOR_5OVER3 = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INIT_RESTRICT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__GPUID_SELECT,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &INIT_SUBSAMPLING_NCELL,   string );

   getline( &input_line, &len, File );



// interpolation schemes
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INT_TIME = (bool)temp_int;

   getline( &input_line, &len, File );
#  if ( MODEL == ELBDM ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INT_PHASE = (bool)temp_int;
#  endif

   getline( &input_line, &len, File ); // skip one comment line

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLU_INT_SCHEME = (IntScheme_t)temp_int;

#  ifdef GRAVITY
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__POT_INT_SCHEME = (IntScheme_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RHO_INT_SCHEME = (IntScheme_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__GRA_INT_SCHEME = (IntScheme_t)temp_int;

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__REF_FLU_INT_SCHEME = (IntScheme_t)temp_int;

   getline( &input_line, &len, File );
#  ifdef GRAVITY
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__REF_POT_INT_SCHEME = (IntScheme_t)temp_int;
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &INT_MONO_COEFF,           string );

   getline( &input_line, &len, File );


// data dump
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_TOTAL = (OptOutputFormat_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_PART = (OptOutputPart_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_TEST_ERROR = (bool)temp_int;

   getline( &input_line, &len, File );
#  ifdef PARTICLE
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_PARTICLE = (bool)temp_int;
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_BASEPS = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_BASE = (bool)temp_int;

   getline( &input_line, &len, File );
#  ifdef GRAVITY
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_POT = (bool)temp_int;
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_MODE = (OptOutputMode_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OUTPUT_STEP,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_DT,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_PART_X,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_PART_Y,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_PART_Z,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &INIT_DUMPID,              string );

   getline( &input_line, &len, File );


// miscellaneous
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__VERBOSE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__TIMING_BALANCE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__TIMING_MPI = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RECORD_MEMORY = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RECORD_PERFORMANCE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__MANUAL_CONTROL = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RECORD_USER = (bool)temp_int;

   getline( &input_line, &len, File );


// simulation checks
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_REFINE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_PROPER_NESTING = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_CONSERVATION = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_RESTRICT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_FINITE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_PATCH_ALLOCATE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_FLUX_ALLOCATE = (bool)temp_int;

   getline( &input_line, &len, File );
#  if   ( MODEL == HYDRO )
   sscanf( input_line, "%d%s",   &OPT__CK_NEGATIVE,         string );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OPT__CK_MEMFREE,          string );

   getline( &input_line, &len, File );
#  ifdef PARTICLE
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_PARTICLE = (bool)temp_int;
#  endif


   fclose( File );

   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_Parameter ... done\n" );

} // FUNCTION : Init_Load_Parameter

