#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"
#ifdef INTEL
#include <mathimf.h>
#endif

#if ( MODEL == HYDRO )

#if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
extern real CPU_PositivePres( const real Pres_In, const real Dens, const real _Dens );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetTimeStep_Fluid
// Description :  Estimate the evolution time-step and physical time interval from the fluid solver condition
//
// Note        :  1. Physical coordinates : dTime == dt
//                   Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
// 
// Parameter   :  dt       : Time interval to advance solution
//                dTime    : Time interval to update physical time 
//                MinDtLv  : Refinement level determining the smallest time-step
//                MinDtVar : Array to store the variables with the maximum speed (minimum time-step) at each level
//                dt_dTime : dt/dTime (== 1.0 if COMOVING is off)
//
// Return      :  dt, dTime, MinDtLv, MinDtVar
//-------------------------------------------------------------------------------------------------------
void Hydro_GetTimeStep_Fluid( double &dt, double &dTime, int &MinDtLv, real MinDtVar[], const double dt_dTime )
{

   real  *MaxCFL   = MinDtInfo_Fluid;  // "MinDtInfo_Fluid" is a global variable
   double dt_local = __FLT_MAX__;      // initialize it as an extremely large number
   double dt_min, dt_tmp; 
   real   MinDtVar_AllLv[NLEVEL][NCOMP];


// get the maximum CFL velocity ( sound speed + fluid velocity )
   if ( !OPT__ADAPTIVE_DT )   Hydro_GetMaxCFL( MaxCFL, MinDtVar_AllLv );


// get the time-step at the base level per sub-step in one rank
   for (int lv=0; lv<NLEVEL; lv++)
   {
      dt_tmp = amr->dh[lv] / (double)MaxCFL[lv];

//    return 2*dt for the individual time-step since at the base level each step actually includes two sub-steps
#     ifdef INDIVIDUAL_TIMESTEP
      dt_tmp *= double( 1<<(lv+1) );
#     endif

      if ( dt_tmp <= 0.0 )
         Aux_Error( ERROR_INFO, "dt_tmp = %14.7e <= 0 (Rank %d, Lv %d) !!\n", dt_tmp, MPI_Rank, lv );
      
      if ( dt_tmp < dt_local )    
      {
         dt_local = dt_tmp;
         MinDtLv  = lv;

         for (int v=0; v<NCOMP; v++)   MinDtVar[v] = MinDtVar_AllLv[lv][v];
      }
   } // for (int lv=0; lv<NLEVEL; lv++)


// get the minimum time-step from all ranks
   MPI_Allreduce( &dt_local, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


// verify the minimum time-step
   if ( dt_min == __FLT_MAX__ )
      Aux_Error( ERROR_INFO, "time-step estimation by %s is incorrect (dt_min = %13.7e) !!\n",
                 __FUNCTION__, dt_min );


// gather the minimum time-step information from all ranks
#  ifndef SERIAL
   double *dt_AllRank                = new double [MPI_NRank];
   int    *MinDtLv_AllRank           = new int    [MPI_NRank];
   real   (*MinDtVar_AllRank)[NCOMP] = new real   [MPI_NRank][NCOMP];
   
   MPI_Gather( &dt_local, 1,     MPI_DOUBLE, dt_AllRank,       1,     MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Gather( &MinDtLv,  1,     MPI_INT,    MinDtLv_AllRank,  1,     MPI_INT,    0, MPI_COMM_WORLD );
#  ifdef FLOAT8
   MPI_Gather( MinDtVar,  NCOMP, MPI_DOUBLE, MinDtVar_AllRank, NCOMP, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Gather( MinDtVar,  NCOMP, MPI_FLOAT,  MinDtVar_AllRank, NCOMP, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   if ( MPI_Rank == 0 )
   {
      for (int Rank=0; Rank<MPI_NRank; Rank++)
      {
         if ( dt_AllRank[Rank] == dt_min )
         {
            MinDtLv = MinDtLv_AllRank[Rank];
            for (int v=0; v<NCOMP; v++)   MinDtVar[v] = MinDtVar_AllRank[Rank][v];
            break;
         }

         if ( Rank == MPI_NRank-1 )    Aux_Message( stderr, "WARNING : no match of \"dt_min\" was found !!\n" );
      }
   }

   delete [] dt_AllRank;
   delete [] MinDtVar_AllRank;
   delete [] MinDtLv_AllRank;
#  endif // #ifndef SERIAL 


   dt    = ( Step == 0 ) ? DT__FLUID_INIT*dt_min : DT__FLUID*dt_min;
   dTime = dt / dt_dTime;

} // FUNCTION : Hydro_GetTimeStep_Fluid



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetMaxCFL
// Description :  Evaluate the maximum (fluid velocity + sound speed) for the time-step estimation
//
// Note        :  This function is also invoked in "Init_GAMER"
// 
// Parameter   :  MaxCFL         : Array to store the maximum speed at each level
//                MinDtVar_AllLv : Array to store the variables with the maximum speed at each level
//-------------------------------------------------------------------------------------------------------
void Hydro_GetMaxCFL( real MaxCFL[], real MinDtVar_AllLv[][NCOMP] )
{

   const real Gamma_m1 = GAMMA - (real)1.0;

#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   const bool PositivePres = true;
#  else
   const bool PositivePres = false;
#  endif

#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   real  *MaxCFL_OMP           = new real [NT];
   real (*MinDtVar_OMP)[NCOMP] = new real [NT][NCOMP];
   real Fluid[NCOMP], _Rho, Vx, Vy, Vz, Pres, Cs, MaxV, MaxCFL_candidate;
   int  TID = 0;  // thread ID


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the arrays MaxCFL and MaxCFL_OMP as extremely small numbers
      MaxCFL[lv] = __FLT_MIN__;
      for (int t=0; t<NT; t++)   MaxCFL_OMP[t] = __FLT_MIN__;

#     pragma omp parallel private( Fluid, _Rho, Vx, Vy, Vz, Pres, Cs, MaxV, MaxCFL_candidate, TID )
      {
#        ifdef OPENMP
         TID = omp_get_thread_num();
#        endif

//       loop over all patches
#        pragma omp for
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
//          if ( amr->patch[0][lv][PID]->son == -1 )  // only check the patches without son
            if ( true )
            {
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
                  for (int v=0; v<NCOMP; v++)   Fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

                 _Rho  = (real)1.0 / Fluid[DENS];
                  Vx   = FABS( Fluid[MOMX] )*_Rho;
                  Vy   = FABS( Fluid[MOMY] )*_Rho;
                  Vz   = FABS( Fluid[MOMZ] )*_Rho;
                  Pres = Hydro_GetPressure( Fluid, Gamma_m1, PositivePres );
                  Cs   = SQRT( GAMMA*Pres*_Rho );

#                 if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU  ||  FLU_SCHEME == WAF )
                  MaxV             = ( Vx > Vy   ) ? Vx : Vy;
                  MaxV             = ( Vz > MaxV ) ? Vz : MaxV;
                  MaxCFL_candidate = MaxV + Cs;

#                 elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
                  MaxV             = Vx + Vy + Vz;
                  MaxCFL_candidate = MaxV + (real)3.0*Cs;
#                 endif

                  if ( MaxCFL_candidate > MaxCFL_OMP[TID] )
                  {
                     MaxCFL_OMP[TID] = MaxCFL_candidate;

                     MinDtVar_OMP[TID][0] = Fluid[DENS];
                     MinDtVar_OMP[TID][1] = Vx;
                     MinDtVar_OMP[TID][2] = Vy;
                     MinDtVar_OMP[TID][3] = Vz;
                     MinDtVar_OMP[TID][4] = Cs;
                  }


//                verify the CFL number
                  if (  ! isfinite( MaxCFL_candidate )  )
                  {
                     Aux_Message( stderr, "ERROR : \"incorrect CFL evaluation (%14.7e)\" !!\n", MaxCFL_candidate);
                     Aux_Message( stderr, "        t %14.7e, Step %ld, Cs %14.7e, MaxV %14.7e, Rho %14.7e\n",
                                  Time[0], Step, Cs, MaxV, Fluid[DENS] );
                     Aux_Message( stderr, "        Vx %14.7e, Vy %14.7e, Vz %14.7e, Engy %14.7e Pres %14.7e\n",
                                  Vx, Vy, Vz, Fluid[ENGY], Pres );
                     Aux_Message( stderr, "        Rank %d, Lv %d, PID  %d, Coordinates (%20.14e,%20.14e,%20.14e)\n",
                                  MPI_Rank, lv, PID, amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv],  
                                                     amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv], 
                                                     amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv]  );
                     Aux_Message( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  
                                                                                           __FUNCTION__  );
                     Aux_Message( stderr, "        **********************************************************\n");
                     Aux_Message( stderr, "        ** Usually this error is caused by the negative density **\n");
                     Aux_Message( stderr, "        ** or pressure evaluated in the fluid solver            **\n");
                     Aux_Message( stderr, "        ** --> You might want to set MIN_PRES or MIN_PRES_DENS  **\n");
                     Aux_Message( stderr, "        **     in the header GAMER/include/CUFLU.h properly     **\n");
                     Aux_Message( stderr, "        **     or choose a smaller time-step                    **\n");
                     Aux_Message( stderr, "        **********************************************************\n");
                     MPI_Exit();
                  }
   
               } // i,j,k
//          } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // if ( true )
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // OpenMP parallel region


//    compare the maximum CFL evaluated by different OMP threads
      for (int t=0; t<NT; t++)
      {
         if ( MaxCFL_OMP[t] > MaxCFL[lv] )   
         {
            MaxCFL[lv] = MaxCFL_OMP[t];

            for (int v=0; v<NCOMP; v++)   MinDtVar_AllLv[lv][v] = MinDtVar_OMP[t][v];
         }
      }
   } // for (int lv=0; lv<NLEVEL; lv++)


   delete [] MaxCFL_OMP;
   delete [] MinDtVar_OMP;

} // FUNCTION : Hydro_GetMaxCFL



#endif // #if ( MODEL == HYDRO )
