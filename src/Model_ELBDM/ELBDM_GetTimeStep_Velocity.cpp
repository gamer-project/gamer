#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

static real GetMaxVelocity( const int lv );


//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Velocity
// Description :  Estimate the evolution time-step by via the CFL condition from the Hamilton-Jacobi equation
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. dt = 0.5 * dx * m/hbar / sum_i |v_i|
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Velocity( const int lv )
{

   real   MaxdS_dx;
   double dt;


// get the velocity dS/dx / dx as first derivative of phase
   MaxdS_dx  = GetMaxVelocity( lv );
   MaxdS_dx += (MaxdS_dx == 0) ? 1e-8 : 0;

// get the time-step
   dt = 1/MaxdS_dx*0.5/ELBDM_ETA * DT__VELOCITY;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Velocity



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxVelocity
// Description :  Evaluate the maximum first spatial derivative of phase for the time-step estimation
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Velocity()
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MaxdS_dx
//-------------------------------------------------------------------------------------------------------
real GetMaxVelocity( const int lv )
{
   // Maximum velocity for calculation of time step zero for wave patches
   // since we are only concerned with a velocity dependent time step criterion for the fluid solver
   //if ( amr -> use_wave_flag[lv] == true )
   //   return 0;


   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;
   const int  NGhost            = 1;               // number of ghost zones to calculate dS_dx
   const int  Size_Flu          = PS2 + 2*NGhost;  // size of the array Flu_Array
   const int  NPG               = 1;               // number of patch groups (must be ONE here)
   const int  NCOMP1            = 1;               // we only need to retrieve phase for velocity criterion
   real (*Flu_Array)[NCOMP1][Size_Flu][Size_Flu][Size_Flu] = NULL;

   real dS_dx, MaxdS_dx, _dh, _dh2, GradS[3];
   int im, ip, jm, jp, km, kp, I, J, K;


   MaxdS_dx = 0.0;
   _dh      = (real)1.0/amr->dh[lv];
   _dh2     = (real)0.5*_dh;



#  pragma omp parallel private( Flu_Array, dS_dx, GradS, \
                                im, ip, jm, jp, km, kp, I, J, K )
   {
      Flu_Array = new real [NPG][NCOMP1][Size_Flu][Size_Flu][Size_Flu];

//    loop over all patches
#     pragma omp for reduction( max:MaxdS_dx ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)
      {
//       prepare phase with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][0][0][0][0], NULL, NGhost, NPG, &PID0, _PHAS, _NONE,
                            INT_MINMOD1D, INT_NONE, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


//       evaluate dS_dt
         for (int k=NGhost; k<Size_Flu-NGhost; k++)    {  km = k - 1;    kp = k + 1;   K = k - NGhost;
         for (int j=NGhost; j<Size_Flu-NGhost; j++)    {  jm = j - 1;    jp = j + 1;   J = j - NGhost;
         for (int i=NGhost; i<Size_Flu-NGhost; i++)    {  im = i - 1;    ip = i + 1;   I = i - NGhost;
            GradS[0]   = _dh2 * ( Flu_Array[0][0][k ][j ][ip] - Flu_Array[0][0][k ][j ][im] );
            GradS[1]   = _dh2 * ( Flu_Array[0][0][k ][jp][i ] - Flu_Array[0][0][k ][jm][i ] );
            GradS[2]   = _dh2 * ( Flu_Array[0][0][kp][j ][i ] - Flu_Array[0][0][km][j ][i ] );

            dS_dx      = _dh * (FABS(GradS[0]) + FABS(GradS[1]) + FABS(GradS[2]));
            MaxdS_dx   = MAX( MaxdS_dx, dS_dx );

         }}} // i,j,k
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)

      delete [] Flu_Array;
   } // OpenMP parallel region


// get the maximum potential in all ranks
   real MaxdS_dx_AllRank;
#  ifdef FLOAT8
   MPI_Allreduce( &MaxdS_dx, &MaxdS_dx_AllRank, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#  else
   MPI_Allreduce( &MaxdS_dx, &MaxdS_dx_AllRank, 1, MPI_FLOAT,  MPI_MAX, MPI_COMM_WORLD );
#  endif


// check
   if ( MaxdS_dx_AllRank == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MaxdS_dx == 0.0 at lv %d !!\n", lv );


   return MaxdS_dx_AllRank;

} // FUNCTION : GetMaxPhaseDerivative



#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
