#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_GetTimeStep_VelAcc
// Description :  Estimate the evolution time-step and physical time interval by the maximum particle velocity
//                and acceleration
//                --> dt_vel = DT__PARVEL*dh/v_max, where v_max = max(vx,vy,vz,all_particles)
//                    dt_acc = DT__PARACC*(dh/a_max)^0.5, where a_max = max(accx,accy,accz,all_particles)
//
// Note        :  1. Physical coordinates : dTime == dt
//                   Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
//                2. Particle acceleration criterion is used only when DT__PARACC > 0.0
//
// Parameter   :  dt       : Time interval to advance solution
//                dTime    : Time interval to update physical time
//                MinDtLv  : Refinement level determining the smallest time-step
//                MinDtVar : Maximum velocity and acceleration determining the minimum time-step
//                           [0/1] ==> [velocity/acceleration]
//                dt_dTime : dt/dTime (== 1.0 if COMOVING is off)
//
// Return      :  dt, dTime, MinDtLv, MinDtVar
//-------------------------------------------------------------------------------------------------------
void Par_GetTimeStep_VelAcc( double dt[2], double dTime[2], int MinDtLv[2], real MinDtVar[2], const double dt_dTime )
{

#  ifdef COMOVING
#  error : ERROR : COMOVING needs to be checked here
#  endif


   const real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
#  ifdef STORE_PAR_ACC
   const real *Acc[3] = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };
#  else
   const real *Acc[3] = { NULL, NULL, NULL };
#  endif
   const bool  UseAcc = ( DT__PARACC > 0.0 );
#  ifdef OPENMP
   const int   NT     = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int   NT     = 1;
#  endif

#  ifndef STORE_PAR_ACC
   if ( UseAcc )
      Aux_Error( ERROR_INFO, "DT__PARACC (%14.7e) > 0.0 when STORE_PAR_ACC is off !!\n", DT__PARACC );
#  endif

   real  *MaxVel      = MinDtInfo_ParVelAcc[0];       // "MinDtInfo_ParVelAcc" is a global variable
   real  *MaxAcc      = MinDtInfo_ParVelAcc[1];       // useless if UseAcc == false
   double dt_local[2] = { __FLT_MAX__, __FLT_MAX__ }; // initialize as extremely large numbers

   double dt_min[2], dt_tmp[2];
   int    TID;
   long   ParID;

   real *MaxVel_OMP = new real [NT];
   real *MaxAcc_OMP = new real [NT];


// get the maximum particle velocity and acceleration at each level
   if ( !OPT__ADAPTIVE_DT )
   for (int lv=0; lv<NLEVEL; lv++)
   {
      MaxVel[lv] = (real)0.0;    // don't assign negative values since we assume it to be positive-definite
      MaxAcc[lv] = (real)0.0;

      for (int t=0; t<NT; t++)
      {
         MaxVel_OMP[t] = (real)0.0;
         MaxAcc_OMP[t] = (real)0.0;
      }

//###NOTE: OpenMP may not improve performance here
#     pragma omp parallel private( TID, ParID )
      {
#        ifdef OPENMP
         TID = omp_get_thread_num();
#        else
         TID = 0;
#        endif

#        pragma omp for schedule( runtime )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

            for (int d=0; d<3; d++)
            MaxVel_OMP[TID] = MAX( MaxVel_OMP[TID], FABS(Vel[d][ParID]) );

            if ( UseAcc )
            for (int d=0; d<3; d++)
            MaxAcc_OMP[TID] = MAX( MaxAcc_OMP[TID], FABS(Acc[d][ParID]) );
         }
      } // end of OpenMP parallel region

//    compare the maximum velocity and acceleration evaluated by different OMP threads
      for (int t=0; t<NT; t++)
      {
         MaxVel[lv] = MAX( MaxVel_OMP[t], MaxVel[lv] );
         if ( UseAcc )
         MaxAcc[lv] = MAX( MaxAcc_OMP[t], MaxAcc[lv] );
      }
   } // for (int lv=0; lv<NLEVEL; lv++)

   delete [] MaxVel_OMP;
   delete [] MaxAcc_OMP;


// get the time-step in one rank
   MinDtLv[0] = -1;  // indicating that MinDtVar cannot be obtained
   MinDtLv[1] = -1;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      dt_tmp[0]  =       amr->dh[lv] / MaxVel[lv];
      if ( UseAcc )
      dt_tmp[1]  = sqrt( amr->dh[lv] / MaxAcc[lv] );

//    return 2*dt for the individual time-step since at the base level each step actually includes two sub-steps
#     ifdef INDIVIDUAL_TIMESTEP
      dt_tmp[0] *= double( 1<<(lv+1) );
      if ( UseAcc )
      dt_tmp[1] *= double( 1<<(lv+1) );
#     endif

      if ( dt_tmp[0] < dt_local[0] )
      {
         dt_local[0] = dt_tmp[0];
         MinDtLv [0] = lv;
         MinDtVar[0] = MaxVel[lv];
      }

      if ( UseAcc )
      if ( dt_tmp[1] < dt_local[1] )
      {
         dt_local[1] = dt_tmp[1];
         MinDtLv [1] = lv;
         MinDtVar[1] = MaxAcc[lv];
      }
   } // for (int lv=0; lv<NLEVEL; lv++)


// get the minimum time-step from all ranks
   MPI_Allreduce( &dt_local[0], &dt_min[0], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

   if ( UseAcc )
   MPI_Allreduce( &dt_local[1], &dt_min[1], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


// verify the minimum time-step
   if ( dt_min[0] == __FLT_MAX__  &&  amr->Par->NPar_Active_AllRank > 0 )
   {
      Aux_Message( stderr, "WARNING : time-step estimation by particle velocity is incorrect (dt_min = %13.7e) !!\n", dt_min[0] );
      Aux_Message( stderr, "          --> Likely all particles have zero velocity\n" );

      if ( DT__PARVEL_MAX < 0.0 )
      Aux_Message( stderr, "          --> You might want to set DT__PARVEL_MAX properly\n" );
   }

   if ( UseAcc  &&  dt_min[1] == __FLT_MAX__  &&  amr->Par->NPar_Active_AllRank > 0 )
      Aux_Error( ERROR_INFO, "time-step estimation by particle acceleration is incorrect (dt_min = %13.7e) !!\n", dt_min[1] );



// gather the minimum time-step information from all ranks
#  ifndef SERIAL
   double (*dt_AllRank      )[2] = new double [MPI_NRank][2];
   int    (*MinDtLv_AllRank )[2] = new int    [MPI_NRank][2];
   real   (*MinDtVar_AllRank)[2] = new real   [MPI_NRank][2];

   MPI_Gather( dt_local, 2, MPI_DOUBLE, dt_AllRank,       2, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Gather( MinDtLv,  2, MPI_INT,    MinDtLv_AllRank,  2, MPI_INT,    0, MPI_COMM_WORLD );
#  ifdef FLOAT8
   MPI_Gather( MinDtVar, 2, MPI_DOUBLE, MinDtVar_AllRank, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Gather( MinDtVar, 2, MPI_FLOAT,  MinDtVar_AllRank, 2, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   if ( MPI_Rank == 0 )
   for (int t=0; t<2; t++)
   {
      for (int Rank=0; Rank<MPI_NRank; Rank++)
      {
         if ( dt_AllRank[Rank][t] == dt_min[t] )
         {
            MinDtLv [t] = MinDtLv_AllRank [Rank][t];
            MinDtVar[t] = MinDtVar_AllRank[Rank][t];
            break;
         }

         if ( Rank == MPI_NRank-1 )
            Aux_Message( stderr, "WARNING : no match of \"dt_min[%d]\" was found !!\n", t );
      }
   }

   delete [] dt_AllRank;
   delete [] MinDtLv_AllRank;
   delete [] MinDtVar_AllRank;
#  endif // #ifndef SERIAL


   dt[0] = DT__PARVEL * dt_min[0];
   if ( DT__PARVEL_MAX >= 0.0 )  dt[0] = MIN( dt[0], DT__PARVEL_MAX );
   dTime[0] = dt[0] / dt_dTime;

   if ( UseAcc )
   {
      dt[1] = DT__PARACC * dt_min[1];
      dTime[1] = dt[1] / dt_dTime;
   }

} // FUNCTION : Par_GetTimeStep_VelAcc



#endif // #ifdef PARTICLE
