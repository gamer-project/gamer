#include "Copyright.h"
#include "GAMER.h"

extern double (*Mis_GetTimeStep_User_Ptr)( const double dTime_dt );




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep
// Description :  Estimate the evolution time-step (dt) and the physical time interval (dTime) at the target
//                refinement level
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the physical time interval (dTime) instead of the evolution time-step (dt)
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dTime, the physical time interval == "delta(scale_factor)"
//                       in the comoving coordinates, back to dt in EvolveLevel()
//                2. The function pointer "Mis_GetTimeStep_User_Ptr" points to "Mis_GetTimeStep_User()" by default
//                   but may be overwritten by various test problem initializers
//
// Parameter   :  lv             : Target refinement level
//                dTime_SyncFaLv : dt to synchronize lv and lv-1
//                                 --> Only used for OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE
//
// Return      :  dTime_min
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep( const int lv, const double dTime_SyncFaLv )
{

   const char  FileName[] = "Record__TimeStep";
   static bool FirstTime  = true;
   const int   NdTimeMax  = 20;

   char  (*dTime_Name)[MAX_STRING] = new char   [NdTimeMax][MAX_STRING];
   double *dTime                   = new double [NdTimeMax];

   int NdTime = 0;


// -1. return immediately if the target level has no patches
// =============================================================================================================
   if ( NPatchTotal[lv] == 0 )   return HUGE_NUMBER;



// 0. estimate the relation between the evolution time-step (dt) and the physical time interval (dTime)
// =============================================================================================================
   double dTime_dt;  // first derivative of dTime over dt
#  ifdef COMOVING
   dTime_dt = pow(  OMEGA_M0*pow( Time[lv], 3.0 ) + (1.0-OMEGA_M0)*pow( Time[lv], 6.0 ),  0.5  );
#  else
   dTime_dt = 1.0;
#  endif



// 1.1 CRITERION ONE : fluid solver
// =============================================================================================================
#  if   ( MODEL == HYDRO )
   dTime[NdTime] = dTime_dt * dt_InvokeSolver( DT_FLU_SOLVER, lv );
   sprintf( dTime_Name[NdTime++], "%s", "Hydro_CFL" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
// ELBDM_GetTimeStep_Fluid( dt1, dTime1, MinDtLv_Fluid, dt_dTime );
   sprintf( dTime_Name[NdTime++], "%s", "ELBDM_CFL" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 1.2 CRITERION TWO : gravitational acceleration
// =============================================================================================================
#  ifdef GRAVITY
#  if   ( MODEL == HYDRO )
   dTime[NdTime] = dTime_dt * dt_InvokeSolver( DT_GRA_SOLVER, lv );
   sprintf( dTime_Name[NdTime++], "%s", "Hydro_Acc" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
// ELBDM_GetTimeStep_Gravity( dt2, dTime2, MinDtLv_Gravity, MinDtVar_Gravity, dt_dTime );
   sprintf( dTime_Name[NdTime++], "%s", "ELBDM_Pot" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
#  endif // #ifdef GRAVITY


// 1.3 CRITERION THREE : maximum allowed variation of the expansion factor
// =============================================================================================================
#  ifdef COMOVING
   dTime[NdTime] = DT__MAX_DELTA_A * Time[lv];
   sprintf( dTime_Name[NdTime++], "%s", "Delta_A" );
#  endif


// 1.4 CRITERION FOUR : match the time of the next data dump
// =============================================================================================================
// DumpByTime : true --> dump data according to the physical time
#  ifdef PARTICLE
   const bool DumpData   = ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS ||
                             OPT__OUTPUT_PAR_TEXT );
#  else
   const bool DumpData   = ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS );
#  endif
   const bool DumpByTime = (  DumpData  &&  ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT || OPT__OUTPUT_MODE == OUTPUT_USE_TABLE )  )
                           ? true : false;

   if ( DumpByTime )
   {
      dTime[NdTime] = DumpTime - Time[lv];

      if ( dTime[NdTime] <= 0.0 )
      {
         Aux_Message( stderr, "ERROR : dTime (%20.14e) <= 0.0, something is wrong !!\n", dTime[NdTime] );
         Aux_Message( stderr, "        (DumpTime %20.14e, Time %20.14e, lv %d)\n", DumpTime, Time[lv], lv );
         MPI_Exit();
      }

      sprintf( dTime_Name[NdTime++], "%s", "Data_Dump" );
   }


// 1.5 CRITERION FIVE : match the program end time
// =============================================================================================================
   dTime[NdTime] = END_T - Time[lv];

   if ( dTime[NdTime] <= 0.0 )
   {
      Aux_Message( stderr, "ERROR : dTime (%20.14e) <= 0.0, something is wrong !!\n", dTime[NdTime] );
      Aux_Message( stderr, "        (END_T %20.14e, Time %20.14e, lv %d)\n", END_T, Time[lv], lv );
      MPI_Exit();
   }

   sprintf( dTime_Name[NdTime++], "%s", "End_Time" );


// 1.6 CRITERION SIX : user-defined criteria
// =============================================================================================================
   if ( OPT__DT_USER  &&  Mis_GetTimeStep_User_Ptr != NULL )
   {
      dTime[NdTime] = dTime_dt * Mis_GetTimeStep_User_Ptr( dTime_dt );
      sprintf( dTime_Name[NdTime++], "%s", "User" );
   }


// 1.7 CRITERION SEVEN : phase rotation ##ELBDM ONLY##
// =============================================================================================================
#  if ( MODEL == ELBDM )
   const bool ELBDM_PhaseDt = ( DT__PHASE != 0.0 ) ? true : false;
   int    MinDtLv_Phase;
   real   MinDtVar_Phase[NCOMP_FLUID];

   if ( ELBDM_PhaseDt )
   {
      ELBDM_GetTimeStep_Phase( dt7, dTime7, MinDtLv_Phase, MinDtVar_Phase, dt_dTime );
      sprintf( dTime_Name[NdTime++], "%s", "ELBDM_Phase" );
   }
#  endif


// 1.8 CRITERION EIGHT : particle evolution
// =============================================================================================================
   /*
#  ifdef PARTICLE
   real   MinDtVar_ParVelAcc[2];
   int    MinDtLv_ParVelAcc[2];

   Par_GetTimeStep_VelAcc( dt8, dTime8, MinDtLv_ParVelAcc, MinDtVar_ParVelAcc, dt_dTime );

   sprintf( dTime_Name[NdTime++], "%s", "Par_Vel" );
   if ( DT__PARACC > 0.0 )
   sprintf( dTime_Name[NdTime++], "%s", "Par_Acc" );
#  endif // #ifdef PARTICLE
   */



// 2. get the minimum time-step from all criteria
// =============================================================================================================
// 2.1 loop over all dt criteria
   double dTime_min = HUGE_NUMBER;

   for (int t=0; t<NdTime; t++)  dTime_min = fmin( dTime_min, dTime[t] );


// 2.2 synchronize with the parent level
   if ( OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE )
   {
      if ( lv > 0 )
      {
         if ( dTime_SyncFaLv <= 0.0 )
         {
            Aux_Message( stderr, "ERROR : dTime_SyncFaLv (%20.14e) <= 0.0, something is wrong !!\n", dTime_SyncFaLv );
            MPI_Exit();
         }

         if ( (1.0+DT__FLEXIBLE_RANGE)*dTime_min >= dTime_SyncFaLv )    dTime_min = dTime_SyncFaLv;
      }

      dTime[NdTime] = dTime_SyncFaLv;
      sprintf( dTime_Name[NdTime++], "%s", "Sync_FaLv" );
   }



// 3. record the dt info
// =============================================================================================================
   if ( OPT__RECORD_DT  &&  MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );

//    header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         fprintf( File, "#%3s  %8s  %8s  %13s  %13s  %13s", "Lv", "Step", "Counter", "TimeOld", "TimeNew", "dTime" );

#        ifdef COMOVING
         fprintf( File, "  %13s", "dTime_dt" );
#        endif

         for (int t=0; t<NdTime; t++)
         fprintf( File, "  %13s", dTime_Name[t] );

         fprintf( File, "\n" );
      }

//    dt info
      fprintf( File, "%4d  %8ld  %8ld  %13.7e  %13.7e  %13.7e",
               lv, Step, AdvanceCounter[lv], Time[lv], Time[lv]+dTime_min, dTime_min );

#     ifdef COMOVING
      fprintf( File, "  %13.7e", dTime_dt );
#     endif

      for (int t=0; t<NdTime; t++)
      fprintf( File, "  %13.7e", dTime[t] );

      fprintf( File, "\n" );

      fclose( File );
   } // if ( OPT__RECORD_DT  &&  MPI_Rank == 0 )



// 4. verify time-step
// =============================================================================================================
   if ( dTime_min <= 0.0  ||  !isfinite(dTime_min) )
      Aux_Error( ERROR_INFO, "incorrect time-step (dTime_min = %20.14e) !!\n", dTime_min );



   delete [] dTime_Name;
   delete [] dTime;
   FirstTime = false;

   return dTime_min;

} // FUNCTION : Mis_GetTimeStep
