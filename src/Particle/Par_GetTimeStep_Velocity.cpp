#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_GetTimeStep_Velocity
// Description :  Estimate the evolution time-step and physical time interval by the maximum particle velocity
//                --> dt = DT__PARVEL*dh/v_max
//
// Note        :  Physical coordinates : dTime == dt
//                Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
// 
// Parameter   :  dt       : Time interval to advance solution
//                dTime    : Time interval to update physical time 
//                MinDtLv  : Refinement level determining the smallest time-step
//                MinDtVar : Maximum velocity determining the minimum time-step
//                dt_dTime : dt/dTime (== 1.0 if COMOVING is off)
//
// Return      :  dt, dTime, MinDtLv, MinDtVar
//-------------------------------------------------------------------------------------------------------
void Par_GetTimeStep_Velocity( double &dt, double &dTime, int &MinDtLv, real &MinDtVar, const double dt_dTime )
{

   const real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

   real  *MaxVel   = MinDtInfo_ParVel;    // "MinDtInfo_ParVel" is a global variable
   double dt_local = __FLT_MAX__;         // initialize it as an extremely large number
   double dt_min, dt_tmp;
   real   TempVel;
   long   ParID;


// get the maximum particle velocity at each level
   if ( !OPT__ADAPTIVE_DT )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
         MaxVel[lv] = __FLT_MIN__;

         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID      = amr->patch[0][lv][PID]->ParList[p];
            TempVel    = SQRT( SQR(Vel[0][ParID]) + SQR(Vel[1][ParID]) + SQR(Vel[2][ParID]) );
            MaxVel[lv] = MAX( TempVel, MaxVel[lv] ); 
         }
      }
   }


// get the time-step in one rank
   for (int lv=0; lv<NLEVEL; lv++)
   {
      dt_tmp  = amr->dh[lv] / MaxVel[lv];

//    return 2*dt for the individual time-step since at the base level each step actually includes two sub-steps
#     ifdef INDIVIDUAL_TIMESTEP
      dt_tmp *= double( 1<<(lv+1) );
#     endif

      if ( dt_tmp < dt_local )    
      {
         dt_local = dt_tmp;
         MinDtLv  = lv;
         MinDtVar = MaxVel[lv];
      }
   }


// get the minimum time-step from all ranks
   MPI_Allreduce( &dt_local, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


// verify the minimum time-step
   if ( dt_min == __FLT_MAX__ )
      Aux_Error( ERROR_INFO, "time-step estimation by particle velocity is incorrect (dt_min = %13.7e) !!\n", dt_min );


// gather the minimum time-step information from all ranks
   /*
#  ifndef SERIAL
   double *dt_AllRank     = new double [MPI_NRank];
   int    *MinDtLv_AllRank  = new int    [MPI_NRank];
   real   *MinDtVar_AllRank = new real   [MPI_NRank];
   
   MPI_Gather( &dt_local, 1, MPI_DOUBLE, dt_AllRank,       1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Gather( &MinDtLv,  1, MPI_INT,    MinDtLv_AllRank,  1, MPI_INT,    0, MPI_COMM_WORLD );
#  ifdef FLOAT8
   MPI_Gather( &MinDtVar, 1, MPI_DOUBLE, MinDtVar_AllRank, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Gather( &MinDtVar, 1, MPI_FLOAT,  MinDtVar_AllRank, 1, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   if ( MPI_Rank == 0 )
   {
      for (int Rank=0; Rank<MPI_NRank; Rank++)
      {
         if ( dt_AllRank[Rank] == dt_min )
         {
            MinDtLv  = MinDtLv_AllRank [Rank];
            MinDtVar = MinDtVar_AllRank[Rank];
            break;
         }

         if ( Rank == MPI_NRank-1 )    Aux_Message( stderr, "WARNING : no match of \"dt_min\" was found !!\n" );
      }
   }

   delete [] dt_AllRank;
   delete [] MinDtLv_AllRank;
   delete [] MinDtVar_AllRank;
#  endif // #ifndef SERIAL 
   */
#  ifndef SERIAL
#  error : ERROR : only SERIAL work here
#  endif

#  ifdef INDIVIDUAL_TIMESTEP
#  error : ERROR : INDIVIDUAL_TIMESTEP needs to be checked here
#  endif

#  ifdef COMOVING
#  error : ERROR : COMOVING needs to be checked here
#  endif


   dt    = DT__PARVEL * dt_min;
   dTime = dt / dt_dTime;

} // FUNCTION : Par_GetTimeStep_Velocity



#endif // #ifdef PARTICLE
