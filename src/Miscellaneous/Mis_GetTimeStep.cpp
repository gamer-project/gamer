#include "Copyright.h"
#include "GAMER.h"

extern void (*Mis_GetTimeStep_User_Ptr)( double &dt, double &dTime, const double dt_dTime );




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep
// Description :  Estimate the evolution time-step (dt) and the physical time interval (dTime)
//
// Note        :  1. Physical coordinates : dTime == dt
//                   Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
//                2. The function pointer "Mis_GetTimeStep_User_Ptr" points to "Mis_GetTimeStep_User()" by default
//                   but may be overwritten by various test problem initializers
//
// Parameter   :  None
//
// Return      :  dTime_min
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep()
{

   const char FileName[] = "Record__TimeStep";
   static bool FirstTime = true;

   if ( MPI_Rank == 0  &&  FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      FirstTime = false;
   }


// 0. estimate the relation between the evolution time-step (dt) and the physical time interval (dTime)
// =============================================================================================================
   double dt_dTime;    // first derivative of dt over dTime
#  ifdef COMOVING
   dt_dTime = pow(  OMEGA_M0*pow( Time[0], 3.0 ) + (1.0-OMEGA_M0)*pow( Time[0], 6.0 ),  -0.5  );
#  else
   dt_dTime = 1.0;
#  endif



// 1.1 CRITERION ONE : fluid solver condition
// =============================================================================================================
   double dTime1, dt1;
   int    MinDtLv_Fluid;

#  if   ( MODEL == HYDRO )
   real MinDtVar_Fluid[5]; // 5: Rho, Vx/y/z, Cs
   Hydro_GetTimeStep_Fluid( dt1, dTime1, MinDtLv_Fluid, MinDtVar_Fluid, dt_dTime );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   ELBDM_GetTimeStep_Fluid( dt1, dTime1, MinDtLv_Fluid, dt_dTime );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL



// 1.2 CRITERION TWO : gravitation acceleration condition
// =============================================================================================================
#  ifdef GRAVITY
   double dTime2, dt2;
   int    MinDtLv_Gravity;

#  if   ( MODEL == HYDRO )
   real MinDtVar_Gravity;
   Hydro_GetTimeStep_Gravity( dt2, dTime2, MinDtLv_Gravity, MinDtVar_Gravity, dt_dTime );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real MinDtVar_Gravity[3];  // [0]: gravitational potential; [1]: lambda*rho (for self-interaction); [2] external potential
   ELBDM_GetTimeStep_Gravity( dt2, dTime2, MinDtLv_Gravity, MinDtVar_Gravity, dt_dTime );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
#  endif  // #ifdef GRAVITY



// 1.3 CRITERION THREE : maximum allowed variation of the expansion factor
// =============================================================================================================
#  ifdef COMOVING
   double dTime3, dt3;

   dTime3 = DT__MAX_DELTA_A * Time[0];
   dt3    = dTime3 * dt_dTime;
#  endif



// 1.4 CRITERION FOUR : fit the time of the next data dump
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

   double dTime4 = NULL_REAL;
   double dt4    = NULL_REAL;

   if ( DumpByTime )
   {
      dTime4 = DumpTime - Time[0];
      dt4    = dTime4 * dt_dTime;

      if ( dTime4 <= 0.0 )
      {
         Aux_Message( stderr, "ERROR : dTime4 (%20.14e) <= 0.0, something is wrong !!\n", dTime4 );
         Aux_Message( stderr, "        (DumpTime %20.14e, Time %20.14e)\n", DumpTime, Time[0] );
         MPI_Exit();
      }
   }



// 1.5 CRITERION FIVE : fit the program end time
// =============================================================================================================
   const double dTime5 = END_T - Time[0];
   const double dt5    = dTime5 * dt_dTime;

   if ( dTime5 <= 0.0 )
   {
      Aux_Message( stderr, "ERROR : dTime5 (%20.14e) <= 0.0, something is wrong !!\n", dTime5 );
      Aux_Message( stderr, "        (END_T %20.14e, Time %20.14e)\n", END_T, Time[0] );
      MPI_Exit();
   }



// 1.6 CRITERION SIX : user-defined criteria
// =============================================================================================================
   double dTime6 = NULL_REAL;
   double dt6    = NULL_REAL;

   if ( OPT__DT_USER  &&  Mis_GetTimeStep_User_Ptr != NULL )   Mis_GetTimeStep_User_Ptr( dt6, dTime6, dt_dTime );



// 1.7 CRITERION SEVEN : phase rotation condition ##ELBDM ONLY##
// =============================================================================================================
#  if ( MODEL == ELBDM )
   const bool ELBDM_PhaseDt = ( DT__PHASE != 0.0 ) ? true : false;
   double dTime7, dt7;
   int    MinDtLv_Phase;
   real   MinDtVar_Phase[NCOMP_FLUID];

   if ( ELBDM_PhaseDt )
   ELBDM_GetTimeStep_Phase( dt7, dTime7, MinDtLv_Phase, MinDtVar_Phase, dt_dTime );
#  endif


// 1.8 CRITERION EIGHT : particle evolution
// =============================================================================================================
#  ifdef PARTICLE
   double dTime8[2], dt8[2];
   real   MinDtVar_ParVelAcc[2];
   int    MinDtLv_ParVelAcc[2];

   Par_GetTimeStep_VelAcc( dt8, dTime8, MinDtLv_ParVelAcc, MinDtVar_ParVelAcc, dt_dTime );
#  endif // #ifdef PARTICLE



// 2. get the minimum time-step from all criteria
// =============================================================================================================
   double dTime_min = dTime1;

#  ifdef GRAVITY
   dTime_min= fmin( dTime_min, dTime2 );
#  endif

#  ifdef COMOVING
   dTime_min= fmin( dTime_min, dTime3 );
#  endif

   if ( DumpByTime )
   dTime_min= fmin( dTime_min, dTime4 );

   dTime_min= fmin( dTime_min, dTime5 );

   if ( OPT__DT_USER  &&  Mis_GetTimeStep_User_Ptr != NULL )
   dTime_min= fmin( dTime_min, dTime6 );

#  if ( MODEL == ELBDM )
   if ( ELBDM_PhaseDt )
   dTime_min= fmin( dTime_min, dTime7 );
#  endif

#  ifdef PARTICLE
   dTime_min= fmin( dTime_min, dTime8[0] );

   if ( DT__PARACC > 0.0 )
   dTime_min= fmin( dTime_min, dTime8[1] );
#  endif



// 3. estimate the evolution time-step (dt)
// =============================================================================================================
   const double dt_min = Mis_dTime2dt( Time[0], dTime_min );



// 4. record the information of time-step determination
// =============================================================================================================
   if ( OPT__RECORD_DT  &&  MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );

      fprintf( File, "Time = %12.6e, Step = %6ld -> %6ld, dt/dTime = %12.6e\n", Time[0], Step, Step+1, dt_dTime );
      fprintf( File, "------------------------------------------------------------------\n" );

#     if   ( MODEL == HYDRO )
      fprintf( File, "CFL Info  : Rho = %12.6e, Vx = %13.6e, Vy = %13.6e, Vz = %13.6e, Cs = %12.6e\n",
               MinDtVar_Fluid[0], MinDtVar_Fluid[1], MinDtVar_Fluid[2], MinDtVar_Fluid[3], MinDtVar_Fluid[4] );
#     elif ( MODEL == ELBDM )
#     ifdef GRAVITY
      if ( ELBDM_PhaseDt )
      fprintf( File, "Phase Info: Lap(Amp)/Amp = %13.6e, Vel^2 = %13.6e, Pot = %13.6e, dPhase_dt = %13.6e\n",
               MinDtVar_Phase[0], MinDtVar_Phase[1], MinDtVar_Phase[2],
               MinDtVar_Phase[0] + MinDtVar_Phase[1] + MinDtVar_Phase[2] );
#     else
      if ( ELBDM_PhaseDt )
      fprintf( File, "Phase Info: Lap(Amp)/Amp = %13.6e, Vel^2 = %13.6e, dPhase_dt = %13.6e\n",
               MinDtVar_Phase[0], MinDtVar_Phase[1], MinDtVar_Phase[0] + MinDtVar_Phase[1] );
#     endif // GRAVITY
#     else
#     warning : WARNING : DO YOU WANT TO PUT SOMETHING HERE FOR THE NEW MODEL ??
#     endif // MODEL

      fprintf( File, "Hydro     : dt = %12.6e, dTime = %12.6e, lv = %2d\n", dt1, dTime1, MinDtLv_Fluid );

#     ifdef GRAVITY
#     if   ( MODEL == HYDRO  ||  MODEL == MHD )
      fprintf( File, "Gravity   : dt = %12.6e, dTime = %12.6e, lv = %2d, MaxAcc = %13.6e\n",
               dt2, dTime2, MinDtLv_Gravity, MinDtVar_Gravity );

#     elif ( MODEL == ELBDM )
      fprintf( File, "Gravity   : dt = %12.6e, dTime = %12.6e, lv = %2d, Max(PotG) = %13.6e",
               dt2, dTime2, MinDtLv_Gravity, MinDtVar_Gravity[0] );
#     ifdef QUARTIC_SELF_INTERACTION
      fprintf( File, ", Max(PotS) = %13.6e", MinDtVar_Gravity[1] );
#     endif
      if ( OPT__EXTERNAL_POT )
      fprintf( File, ", Max(PotE) = %13.6e", MinDtVar_Gravity[2] );

      fprintf( File, "\n" );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
#     endif // #ifdef GRAVITY

#     if ( MODEL == ELBDM )
      if ( ELBDM_PhaseDt )
      fprintf( File, "Phase     : dt = %12.6e, dTime = %12.6e, lv = %2d\n", dt7, dTime7, MinDtLv_Phase );
#     endif

#     ifdef PARTICLE
      fprintf( File, "Particle  : dt = %12.6e, dTime = %12.6e, lv = %2d, MaxVel = %13.6e\n",
               dt8[0], dTime8[0], MinDtLv_ParVelAcc[0], MinDtVar_ParVelAcc[0] );

      if ( DT__PARACC > 0.0 )
      fprintf( File, "            dt = %12.6e, dTime = %12.6e, lv = %2d, MaxAcc = %13.6e\n",
               dt8[1], dTime8[1], MinDtLv_ParVelAcc[1], MinDtVar_ParVelAcc[1] );
#     endif

#     ifdef COMOVING
      fprintf( File, "Delta A   : dt = %12.6e, dTime = %12.6e\n", dt3, dTime3 );
#     endif

      if ( DumpByTime )
      fprintf( File, "Data Dump : dt = %12.6e, dTime = %12.6e\n", dt4, dTime4 );

      if ( dTime_min == dTime5 )
      fprintf( File, "End Time  : dt = %12.6e, dTime = %12.6e\n", dt5, dTime5 );

      if ( OPT__DT_USER  &&  Mis_GetTimeStep_User_Ptr != NULL )
      fprintf( File, "User      : dt = %12.6e, dTime = %12.6e\n", dt6, dTime6 );

      fprintf( File, "Minimum   : dt = %12.6e, dTime = %12.6e\n", dt_min, dTime_min );
      fprintf( File, "\n" );

      fclose( File );

   } // if ( OPT__RECORD_DT  &&  MPI_Rank == 0 )


// 5. verify time-step
// =============================================================================================================
   if ( dt_min <= 0.0  ||  dTime_min<= 0.0  ||  !isfinite(dt_min)  ||  !isfinite(dTime_min) )
      Aux_Error( ERROR_INFO, "incorrect time-step (dt = %20.14e, dTime = %20.14e) !!\n", dt_min, dTime_min );


   return dTime_min;

} // FUNCTION : Mis_GetTimeStep
