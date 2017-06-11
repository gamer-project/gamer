#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_SetDefaultParameter
// Description :  1. Set parameters to their default values if they are not specified in the Input__Parameter
//                2. Reset parameters which are either unsupported or useless
//-------------------------------------------------------------------------------------------------------
void Init_SetDefaultParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


// set parameters to their default values
// ------------------------------------------------------------------------------------------------------
// set the number of OpenMP threads
#  ifndef OPENMP
   if ( OMP_NTHREAD != 1  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 since \"OPENMP\" is disabled !!\n",
                   "OMP_NTHREAD" );

   OMP_NTHREAD = 1;
#  endif


// time-step factors
   if ( DT__FLUID < 0.0 )
   {
#     if   ( MODEL == HYDRO )
#     if   ( FLU_SCHEME == RTVD )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == WAF )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == MHM )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == MHM_RP )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == CTU )
      DT__FLUID = 0.50;
#     else
#     error : unsupported CPU hydro scheme
#     endif

#     elif  ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif  ( MODEL == ELBDM )
#     ifdef GRAVITY
      DT__FLUID = 0.20;                   // 1D k-max mode rotates 0.20*2*PI
#     else
#     ifdef LAPLACIAN_4TH
      DT__FLUID = SQRT(27.0)*M_PI/32.0;   // stability limit (~0.51)
#     else
      DT__FLUID = SQRT(3.0)*M_PI/8.0;     // stability limit (~0.68)
#     endif
#     endif // #ifdef GRAVITY ... else ...

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                        "DT__FLUID", DT__FLUID );
   } // if ( DT__FLUID < 0.0 )

   if ( DT__FLUID_INIT < 0.0 )
   {
      DT__FLUID_INIT = DT__FLUID;

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                        "DT__FLUID_INIT", DT__FLUID_INIT );
   }

#  ifdef GRAVITY
   if ( DT__GRAVITY < 0.0 )
   {
#     if   ( MODEL == HYDRO )
      DT__GRAVITY = 0.50;

#     elif  ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif  ( MODEL == ELBDM )
      DT__GRAVITY = 0.125;

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                        "DT__GRAVITY", DT__GRAVITY );
   } // if ( DT__GRAVITY < 0.0 )
#  endif

#  if ( defined PARTICLE  &&  !defined STORE_PAR_ACC )
   if ( DT__PARACC != 0.0 )
   {
      DT__PARACC = 0.0;    // disable it

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to %13.7e when \"%s\" is off !!\n",
                      "DT__PARACC", DT__PARACC, "STORE_PAR_ACC" );
   }
#  endif


// SOR parameters
#  ifdef GRAVITY
#  if   ( POT_SCHEME == SOR )
   Init_Set_Default_SOR_Parameter( SOR_OMEGA, SOR_MAX_ITER, SOR_MIN_ITER );
#  elif ( POT_SCHEME == MG  )
   Init_Set_Default_MG_Parameter( MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, MG_TOLERATED_ERROR );
#  endif
#  endif // GRAVITY


// set the GPU parameters to the default values when using CPU only (please set OMP_NTHREAD in advance)
#  ifndef GPU
   GPU_NSTREAM = 1;
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                      "GPU_NSTREAM", GPU_NSTREAM );

   if ( FLU_GPU_NPGROUP <= 0 )
   {
#     ifdef OPENMP
      FLU_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      FLU_GPU_NPGROUP = 1;
#     endif

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "FLU_GPU_NPGROUP", FLU_GPU_NPGROUP );
   }

#  ifdef GRAVITY
   if ( POT_GPU_NPGROUP <= 0 )
   {
#     ifdef OPENMP
      POT_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      POT_GPU_NPGROUP = 1;
#     endif

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "POT_GPU_NPGROUP", POT_GPU_NPGROUP );
   }
#  endif
#  endif // #ifndef GPU


// grid size in different refinement levels and box size and scale in different directions
   int NX0_Max;
   NX0_Max = ( NX0_TOT[0] > NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_Max = ( NX0_TOT[2] > NX0_Max    ) ? NX0_TOT[2] : NX0_Max;

   for (int lv=0; lv<NLEVEL; lv++)     amr->dh[lv] = BOX_SIZE / (double)( NX0_Max*(1<<lv) );

   for (int d=0; d<3; d++)
   {
      amr->BoxSize [d] = NX0_TOT[d]*amr->dh   [0];
      amr->BoxScale[d] = NX0_TOT[d]*amr->scale[0];
   }


// whether of not to allocate fluxes at the coarse-fine boundaries
#  if   ( MODEL == HYDRO )
   if ( OPT__FIXUP_FLUX )  amr->WithFlux = true;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( OPT__FIXUP_FLUX )  amr->WithFlux = true;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// (7) ELBDM parameters
#  if ( MODEL == ELBDM )
   ELBDM_ETA = ELBDM_MASS / ELBDM_PLANCK_CONST;

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "NOTE : ELBDM_ETA is set to %13.7e\n", ELBDM_ETA );

#     ifdef COMOVING
      const double JeansK = pow( 6.0*A_INIT*SQR(ELBDM_ETA), 0.25 );
      Aux_Message( stdout, "Note : initial Jean's wavenumber = %13.7e (corresponding to %13.7e Mpc/h)\n",
                   JeansK, 2.0*M_PI/JeansK );
#     endif
   }

   if ( ELBDM_TAYLOR3_AUTO )
   {
      ELBDM_TAYLOR3_COEFF = NULL_REAL;

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to %13.7e since ELBDM_TAYLOR3_AUTO is on\n",
                                        "ELBDM_TAYLOR3_COEFF", ELBDM_TAYLOR3_COEFF );
   }
#  endif // #if ( MODEL == ELBDM )


// interpolation scheme: fluid solver and fluid refinement
#  if   ( MODEL == HYDRO )
   if ( OPT__FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__FLU_INT_SCHEME = INT_CQUAD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_FLU_INT_SCHEME = INT_CQUAD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );
   }

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( OPT__FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__FLU_INT_SCHEME = INT_CQUAR;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_FLU_INT_SCHEME = INT_CQUAR;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );
   }

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

// interpolation scheme: monotonicity coefficient
   if ( INT_MONO_COEFF < 0.0 )
   {
      INT_MONO_COEFF = 2.0;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "INT_MONO_COEFF", INT_MONO_COEFF );
   }



// initial dump ID
   if ( INIT_DUMPID < 0 )  DumpID = 0;
   else                    DumpID = INIT_DUMPID;



// ResPower2 in the AMR_t structure
   int NBits0, NX0_TOT_Max;

   NX0_TOT_Max = ( NX0_TOT[0] > NX0_TOT[1]  ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_TOT_Max = ( NX0_TOT[2] > NX0_TOT_Max ) ? NX0_TOT[2] : NX0_TOT_Max;
   NBits0      = (int)log2( (double)NX0_TOT_Max );

   if (  ( NX0_TOT_Max & (NX0_TOT_Max-1) ) != 0  )    NBits0 ++;  // check if NX0_TOT_Max is a power of two

   for (int lv=0; lv<NLEVEL; lv++)  amr->ResPower2[lv] = NBits0 + lv;


// particle options
#  ifdef PARTICLE
   if ( OPT__BC_POT == BC_POT_ISOLATED  &&  amr->Par->RemoveCell < 0.0 )
   {
      switch ( amr->Par->Interp )
      {
//       set amr->Par->RemoveCell to the distance where potential extrapolation is required
         case ( PAR_INTERP_NGP ):   amr->Par->RemoveCell = 1.0;   break;
         case ( PAR_INTERP_CIC ):   amr->Par->RemoveCell = 1.5;   break;
         case ( PAR_INTERP_TSC ):   amr->Par->RemoveCell = 2.0;   break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "PAR_REMOVE_CELL", amr->Par->RemoveCell );
   }

// set the number of ghost zones for the interpolation scheme
   if ( amr->Par->GhostSize < 0 )
   {
      switch ( amr->Par->Interp )
      {
         case ( PAR_INTERP_NGP ): amr->Par->GhostSize = 0;  break;
         case ( PAR_INTERP_CIC ): amr->Par->GhostSize = 1;  break;
         case ( PAR_INTERP_TSC ): amr->Par->GhostSize = 1;  break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "amr->Par->GhostSize", amr->Par->GhostSize );
   }
#  endif // #ifdef PARTICLE


// coefficient of the Green's function at the origin (must be set after setting amr->Par->Interp )
#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_ISOLATED  &&  GFUNC_COEFF0 < 0.0 )
   {
#     ifdef PARTICLE
      switch ( amr->Par->Interp )
      {
         case ( PAR_INTERP_NGP ):   GFUNC_COEFF0 = 0.0;   break;
         case ( PAR_INTERP_CIC ):   GFUNC_COEFF0 = 4.0;   break;
         case ( PAR_INTERP_TSC ):   GFUNC_COEFF0 = 4.8;   break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }
#     else
      GFUNC_COEFF0 = 0.0;
#     endif

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "GFUNC_COEFF0", GFUNC_COEFF0 );
   }
#  endif


// 1st-order flux correction
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_NONE  &&  OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_NONE )
   {
      OPT__1ST_FLUX_CORR_SCHEME = RSOLVER_1ST_NONE;

      if ( MPI_Rank == 0 )  Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to %d since OPT__1ST_FLUX_CORR is off !!\n",
                                         "OPT__1ST_FLUX_CORR_SCHEME", OPT__1ST_FLUX_CORR_SCHEME );
   }
#  endif


// timing options
   if ( OPT__TIMING_BARRIER < 0 )
   {
#     ifdef TIMING
      if ( OPT__TIMING_BALANCE  ||  OPT__TIMING_MPI )
         OPT__TIMING_BARRIER = 1;
      else
#     endif
         OPT__TIMING_BARRIER = 0;

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                           "OPT__TIMING_BARRIER", OPT__TIMING_BARRIER );
   }


// physical time
   for (int lv=0; lv<NLEVEL; lv++)
   {
#     ifdef COMOVING
      Time     [lv] = A_INIT;          // will be overwritten during restart
#     else
      Time     [lv] = 0.0;             // will be overwritten during restart
#     endif
      Time_Prev[lv] = -__FLT_MAX__;    // initialized as arbitrary "negative" number to indicate that it has not been set yet

      amr->FluSgTime[lv][   amr->FluSg[lv] ] = Time[lv];
#     ifdef GRAVITY
      amr->PotSgTime[lv][   amr->PotSg[lv] ] = Time[lv];
#     endif

      amr->FluSgTime[lv][ 1-amr->FluSg[lv] ] = Time_Prev[lv];
#     ifdef GRAVITY
      amr->PotSgTime[lv][ 1-amr->PotSg[lv] ] = Time_Prev[lv];
#     endif
   }


// particle weight for load balance
#  if ( defined LOAD_BALANCE  &&  defined PARTICLE )
   if ( LB_INPUT__PAR_WEIGHT < 0.0 )
   {
      LB_INPUT__PAR_WEIGHT = 0.0;

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                           "LB_INPUT__PAR_WEIGHT", LB_INPUT__PAR_WEIGHT );
   }
#  endif


// OPT__CORR_AFTER_ALL_SYNC
   if ( OPT__CORR_AFTER_ALL_SYNC == CORR_DEFAULT )
   {
#     ifdef GAMER_DEBUG
      OPT__CORR_AFTER_ALL_SYNC = CORR_EVERY_STEP;
#     else
      OPT__CORR_AFTER_ALL_SYNC = CORR_BEFORE_DUMP;
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                           "OPT__CORR_AFTER_ALL_SYNC", OPT__CORR_AFTER_ALL_SYNC );
   }




// reset parameters and options which are either unsupported or useless
// ------------------------------------------------------------------------------------------------------
// general
// disable "OPT__ADAPTIVE_DT" (not supported yet)
   if ( OPT__ADAPTIVE_DT )
   {
      OPT__ADAPTIVE_DT = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since it is not supported yet !!\n",
                      "OPT__ADAPTIVE_DT" );
   }

// disable "OPT__OVERLAP_MPI" if "OVERLAP_MPI" is NOT turned on in the Makefile
#  ifndef OVERLAP_MPI
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since \"%s\" is off in the Makefile !!\n",
                      "OPT__OVERLAP_MPI", "OVERLAP_MPI" );
   }
#  endif

// disable "OPT__CK_FLUX_ALLOCATE" if no flux arrays are going to be allocated
   if ( OPT__CK_FLUX_ALLOCATE  &&  !amr->WithFlux )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since no flux is required !!\n",
                      "OPT__CK_FLUX_ALLOCATE" );
   }


// for the shared time-step integration
#  ifndef INDIVIDUAL_TIMESTEP
// disable "OPT__INT_TIME"
   if ( OPT__INT_TIME )
   {
      OPT__INT_TIME = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled in the shared time-step scheme !!\n",
                      "OPT__INT_TIME" );
   }
#  endif


// for serial mode
#  ifdef SERIAL
// reset the MPI rank
   if ( MPI_NRank != 1 )
   {
      MPI_NRank = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank" );
   }

   if ( MPI_NRank_X[0] != 1 )
   {
      MPI_NRank_X[0] = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank_X[0]");
   }

   if ( MPI_NRank_X[1] != 1 )
   {
      MPI_NRank_X[1] = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank_X[1]");
   }

   if ( MPI_NRank_X[2] != 1 )
   {
      MPI_NRank_X[2]= 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank_X[2]");
   }

// turn off "OPT__OVERLAP_MPI"
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled in the SERIAL mode !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif // ifdef SERIAL


// turn off "OPT__OVERLAP_MPI" if LOAD_BALANCE is off
#  ifndef LOAD_BALANCE
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since LOAD_BALANCE is NOT turned on !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif // #ifndef LOAD_BALANCE


// always turn on "OPT__VERBOSE" in the debug mode
#  ifdef GAMER_DEBUG
   if ( !OPT__VERBOSE )
   {
      OPT__VERBOSE = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is turned on automatically in the debug mode !!\n",
                      "OPT__VERBOSE" );
   }
#  endif


// for OpenMP
#  ifndef OPENMP
// turn off "OPT__OVERLAP_MPI" if OPENMP is not enabled
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since OPENMP is NOT turned on !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif

// for parallel mode
#  ifndef SERIAL
// check the level of MPI thread support
   int MPI_Thread_Status;
   MPI_Query_thread( &MPI_Thread_Status );
   if ( OPT__OVERLAP_MPI  &&  MPI_Thread_Status == MPI_THREAD_SINGLE )
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled if the level of MPI thread support == %s !!\n",
                      "OPT__OVERLAP_MPI", "MPI_THREAD_SINGLE" );
   }
#  endif


// for different modes
#  if ( MODEL != HYDRO  &&  MODEL != MHD  &&  MODEL != ELBDM )
// operations related to FLUX are useful in HYDRO/MHD/ELBDM only
   if ( OPT__FIXUP_FLUX )
   {
      OPT__FIXUP_FLUX = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is only supported in HYDRO/MHD/ELBDM and hence is disabled !!\n",
                      "OPT__FIXUP_FLUX" );
   }

   if ( OPT__CK_FLUX_ALLOCATE )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is only supported in HYDRO/MHD and hence is disabled !!\n",
                      "OPT__CK_FLUX_ALLOCATE" );
   }
#  endif // #if ( MODEL == HYDRO  &&  MODEL != MHD )

// turn off refinement criteria and checks related to density if "DENS" is not defined
#  ifndef DENS
   if ( OPT__FLAG_RHO )
   {
      OPT__FLAG_RHO = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since the variable DENS is not defined !!\n",
                      "OPT__FLAG_RHO" );
   }

   if ( OPT__FLAG_RHO_GRADIENT )
   {
      OPT__FLAG_RHO_GRADIENT = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since the variable DENS is not defined !!\n",
                      "OPT__FLAG_RHO_GRADIENT" );
   }

   if ( OPT__CK_REFINE )
   {
      OPT__CK_REFINE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is disabled since the variable DENS is not defined !!\n",
                      "OPT__CK_REFINE" );
   }
#  endif // #ifndef DENS

// conservation check is supported only in the models HYDRO, MHD, and ELBDM
#  if ( MODEL != HYDRO  &&  MODEL != MHD  &&  MODEL != ELBDM )
   if ( OPT__CK_CONSERVATION )
   {
      OPT__CK_CONSERVATION = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported in this MODEL and hence is disabled !!\n",
                      "OPT__CK_CONSERVATION" );
   }
#  endif

// set OPT__LR_LIMITER and OPT__WAF_LIMITER to NONE if they are useless (in HYDRO)
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
   if ( OPT__LR_LIMITER != LR_LIMITER_NONE )
   {
      OPT__LR_LIMITER = LR_LIMITER_NONE;

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stderr, "WARNING : parameter \"%s\" is useless in the adopted hydro scheme ", "OPT__LR_LIMITER" );
         Aux_Message( stderr, "and has been set to NONE !!\n" );
      }
   }
#  endif

#  if ( FLU_SCHEME != WAF )
   if ( OPT__WAF_LIMITER != WAF_LIMITER_NONE )
   {
      OPT__WAF_LIMITER = WAF_LIMITER_NONE;

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stderr, "WARNING : parameter \"%s\" is useless in the adopted hydro scheme ", "OPT__WAF_LIMITER" );
         Aux_Message( stderr, "and has been set to NONE !!\n" );
      }
   }
#  endif

#  ifndef GRAVITY
   if ( OPT__FLAG_JEANS )
   {
      OPT__FLAG_JEANS = false;

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stderr, "WARNING : parameter \"%s\" is useless when GRAVITY is off and hence is disabled !!\n",
                      "OPT__FLAG_JEANS" );
      }
   }
#  endif
#  endif // #if ( MODEL == HYDRO  ||  MODEL == MHD )

// operations related to FLUX are useful in ELBDM only if CONSERVE_MASS is on
#  if ( MODEL == ELBDM  &&  !defined CONSERVE_MASS )
   if ( OPT__FIXUP_FLUX )
   {
      OPT__FIXUP_FLUX = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is only supported if CONSERVE_MASS is on and hence is disabled !!\n",
                      "OPT__FIXUP_FLUX" );
   }

   if ( OPT__CK_FLUX_ALLOCATE )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is only supported if CONSERVE_MASS is on and hence is disabled !!\n",
                      "OPT__CK_FLUX_ALLOCATE" );
   }
#  endif // if ( MODEL == ELBDM  &&  !defined CONSERVE_MASS )


// currently OPT__OUTPUT_BASEPS is not supported if the self-gravity is disabled
#  ifndef GRAVITY
   if ( OPT__OUTPUT_BASEPS )
   {
      OPT__OUTPUT_BASEPS = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported when GRAVITY is off and hence is disabled !!\n",
                      "OPT__OUTPUT_BASEPS" );
   }
#  endif


// the option "OPT__RECORD_PERFORMANCE" must work with TIMING
#  ifndef TIMING
   if ( OPT__RECORD_PERFORMANCE )
   {
      OPT__RECORD_PERFORMANCE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported when TIMING is off and hence is disabled !!\n",
                      "OPT__RECORD_PERFORMANCE" );
   }
#  endif


// MPI_NRank_X is useless for restart if LOAD_BALANCE is on
#  ifdef LOAD_BALANCE
   if ( OPT__INIT == INIT_RESTART )
   {
      for (int d=0; d<3; d++)    MPI_NRank_X[d] = -1;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is useless during restart for LOAD_BALANCE\n", "MPI_NRANK_X" );
   }
#  endif


// timing options
#  ifndef TIMING
   if ( OPT__TIMING_BARRIER != 0 )
   {
      OPT__TIMING_BARRIER = 0;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported when TIMING is off and hence is disabled !!\n",
                      "OPT__TIMING_BARRIER" );
   }

   if ( OPT__TIMING_BALANCE )
   {
      OPT__TIMING_BALANCE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported when TIMING is off and hence is disabled !!\n",
                      "OPT__TIMING_BALANCE" );
   }

   if ( OPT__TIMING_MPI )
   {
      OPT__TIMING_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported when TIMING is off and hence is disabled !!\n",
                      "OPT__TIMING_MPI" );
   }
#  endif // #ifndef

#  ifndef LOAD_BALANCE
   if ( OPT__TIMING_MPI )
   {
      OPT__TIMING_MPI = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is not supported when LOAD_BALANCE is off and hence is disabled !!\n",
                      "OPT__TIMING_MPI" );
   }
#  endif


// OPT__UM_START_DOWNGRADE must be turned on for the isolated Poisson solver
#  ifdef GRAVITY
   if ( !OPT__UM_START_DOWNGRADE  &&  OPT__BC_POT == BC_POT_ISOLATED  &&  OPT__UM_START_LEVEL > 0 )
   {
      OPT__UM_START_DOWNGRADE = true;

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stderr, "WARNING : parameter \"%s\" must be turned on for the isolated BC. of gravity !!\n",
                      "OPT__UM_START_DOWNGRADE" );
         Aux_Message( stderr, "          --> It has been switched on automatically.\n" );
      }
   }
#  endif


// OPT__OUTPUT_TOTAL == OUTPUT_FORMAT_HDF5 is not supported if "SUPPORT_HDF5" is not defined
#  ifndef SUPPORT_HDF5
   if ( OPT__OUTPUT_TOTAL == OUTPUT_FORMAT_HDF5 )
   {
      OPT__OUTPUT_TOTAL = OUTPUT_FORMAT_CBINARY;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 2 (\"%s\") since \"%s\" is off !!\n",
                      "OPT__OUTPUT_TOTAL", "OUTPUT_FORMAT_CBINARY", "SUPPORT_HDF5" );
   }
#  endif


// always turn on "OPT__CK_PARTICLE" in the debug mode
#  ifdef DEBUG_PARTICLE
   if ( !OPT__CK_PARTICLE )
   {
      OPT__CK_PARTICLE = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is turned on automatically in the debug mode !!\n",
                      "OPT__CK_PARTICLE" );
   }
#  endif


// RemoveCell is useless for periodic B.C.
#  ifdef PARTICLE
   if ( OPT__BC_POT == BC_POT_PERIODIC  &&  amr->Par->RemoveCell >= 0.0 )
   {
      amr->Par->RemoveCell = -1.0;

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stderr, "WARNING : \"%s\" is useless in the periodic BC ", "PAR_REMOVE_CELL" );
         Aux_Message( stderr, "and has been reset to %14.7e !!\n", amr->Par->RemoveCell );
      }
   }
#  endif // #ifdef PARTICLE


// set particle initialization mode to PAR_INIT_BY_RESTART for restart
#  ifdef PARTICLE
   if ( OPT__INIT == INIT_RESTART  &&  amr->Par->Init != PAR_INIT_BY_RESTART )
   {
      amr->Par->Init = PAR_INIT_BY_RESTART;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : \"%s\" is reset set to \"%s\" for restart !!\n",
                      "PAR_INIT", "PAR_INIT_BY_RESTART" );
   }
#  endif


// always enable OPT__CORR_AFTER_ALL_SYNC in the debug mode
// --> but currently it is allowed to have "OPT__CORR_AFTER_ALL_SYNC == CORR_BEFORE_DUMP" in the debug mode if set by users
#  ifdef GAMER_DEBUG
   if ( OPT__CORR_AFTER_ALL_SYNC == CORR_NONE )
   {
      OPT__CORR_AFTER_ALL_SYNC = CORR_EVERY_STEP;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to %d in the debug mode !!\n",
                      "OPT__CORR_AFTER_ALL_SYNC", OPT__CORR_AFTER_ALL_SYNC );
   }
#  endif


// turn off OPT__NORMALIZE_PASSIVE if there is no passive scalars
#  if (  NCOMP_PASSIVE <= 0  ||  ( defined DUAL_ENERGY && NCOMP_PASSIVE == 1 )  )
   if ( OPT__NORMALIZE_PASSIVE )
   {
      OPT__NORMALIZE_PASSIVE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is turned off automatically since there is no passive scalar !!\n",
                      "OPT__NORMALIZE_PASSIVE" );
   }
#  endif


// turn off OPT__CK_NORMALIZE_PASSIVE when OPT__NORMALIZE_PASSIVE is off
#  if ( NCOMP_PASSIVE > 0 )
   if ( OPT__CK_NORMALIZE_PASSIVE  &&  !OPT__NORMALIZE_PASSIVE )
   {
      OPT__CK_NORMALIZE_PASSIVE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter \"%s\" is turned off automatically when \"%s\" is off !!\n",
                      "OPT__CK_NORMALIZE_PASSIVE", "OPT__NORMALIZE_PASSIVE" );
   }
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_SetDefaultParameter
