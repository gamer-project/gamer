#include "GAMER.h"

#if ( ELBDM_SCHEME == ELBDM_HYBRID )

static real GetMaxVelocity( const int lv, const bool ExcludeWaveCells);




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Hybrid_Velocity
// Description :  Estimate the evolution time-step via the CFL condition from the Hamilton-Jacobi equation
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical Coordinates : dt = physical time interval
//                       Comoving Coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving Coordinates, in Mis_GetTimeStep()
//                2. dt = 0.5 * dx * m/hbar / sum_i |v_i|
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Hybrid_Velocity( const int lv )
{

   const bool ExcludeWaveCells = true;

// get the velocity dS/dx * hbar/m as first derivative of phase
   const real MaxV = GetMaxVelocity( lv, ExcludeWaveCells );

// get the time-step
   double dt = 0.5 * amr->dh[lv] / MaxV;
   dt *= (Step==0) ? DT__HYBRID_VELOCITY_INIT : DT__HYBRID_VELOCITY;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Hybrid_Velocity




//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxVelocity
// Description :  Evaluate the maximum first spatial derivative of phase for the time step estimation
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Hybrid_Velocity()
//
// Parameter   :  lv               : Target refinement level
//                ExcludeWaveCells : Do not include cells with refined wave counterpart in time step estimation if set to True
//
// Return      :  Maximum velocity across all ranks
//-------------------------------------------------------------------------------------------------------
real GetMaxVelocity( const int lv, const bool ExcludeWaveCells )
{

// maximum velocity for calculation of time step is zero for wave patches
// since we are only concerned with a velocity-dependent time step criterion for the fluid solver
   if ( amr->use_wave_flag[lv] )
      return 0;


   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;
   const int  NGhost            = 1;               // number of ghost zones to calculate V
   const int  Size_Flu          = PS2 + 2*NGhost;  // size of the array Flu_Array
   const int  NPG               = 1;               // number of patch groups (must be ONE here)
   const int  NComp1            = 1;               // we only need to retrieve phase for velocity criterion
   real (*Flu_Array)[NComp1][Size_Flu][Size_Flu][Size_Flu] = NULL;

   real V, MaxV, _dh, _dh2, GradS[3];
   int im, ip, jm, jp, km, kp, I, J, K;

   MaxV = 0.0;
   _dh  = (real)1.0/amr->dh[lv];
   _dh2 = (real)0.5*_dh;

#  pragma omp parallel private( Flu_Array, V, GradS, im, ip, jm, jp, km, kp, I, J, K )
   {
      Flu_Array = new real [NPG][NComp1][Size_Flu][Size_Flu][Size_Flu];

//    loop over all patches
#     pragma omp for reduction( max:MaxV ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)
      {
//       prepare phase with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][0][0][0][0], NULL, NGhost, NPG, &PID0, _PHAS, _NONE,
                            INT_CQUAD, INT_NONE, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

//       evaluate dS_dt
         for (int k=NGhost; k<Size_Flu-NGhost; k++)    {  km = k - 1;    kp = k + 1;   K = k - NGhost;
         for (int j=NGhost; j<Size_Flu-NGhost; j++)    {  jm = j - 1;    jp = j + 1;   J = j - NGhost;
         for (int i=NGhost; i<Size_Flu-NGhost; i++)    {  im = i - 1;    ip = i + 1;   I = i - NGhost;

//          skip velocities of cells that have a wave counterpart on the refined levels
            long GID0                   = GlobalTree->PID2GID( PID0, lv );
            bool DoNotCalculateVelocity = false;
            if ( ExcludeWaveCells )
            {
//             in principle, we only need to check the single patch (I, J, K) belongs to
//             however, ELBDM_HasWaveCounterpart automatically returns False if (I, J, K) is outside of the patch GID0 + LocalID
               for (int LocalID=0; LocalID<8; LocalID++ )
               {
                  DoNotCalculateVelocity |= ELBDM_HasWaveCounterpart( I, J, K, GID0, GID0 + LocalID, *GlobalTree );
               }
            }

            if ( !DoNotCalculateVelocity ) {

               GradS[0] = _dh2*( Flu_Array[0][0][k ][j ][ip] - Flu_Array[0][0][k ][j ][im] );
               GradS[1] = _dh2*( Flu_Array[0][0][k ][jp][i ] - Flu_Array[0][0][k ][jm][i ] );
               GradS[2] = _dh2*( Flu_Array[0][0][kp][j ][i ] - Flu_Array[0][0][km][j ][i ] );

               V    = FABS( GradS[0] ) + FABS( GradS[1] ) + FABS( GradS[2] );
               MaxV = MAX( MaxV, V );
            }
         }}} // i,j,k
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)

      delete [] Flu_Array;
   } // OpenMP parallel region

// get the maximum potential in all ranks
   real MaxV_AllRank;
   MPI_Allreduce( &MaxV, &MaxV_AllRank, 1, MPI_GAMER_REAL, MPI_MAX, MPI_COMM_WORLD );

   return MaxV_AllRank / ELBDM_ETA;

} // FUNCTION : GetMaxVelocity



#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
