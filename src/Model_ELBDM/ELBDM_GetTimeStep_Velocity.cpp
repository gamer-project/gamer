#include "GAMER.h"

#if ( MODEL == ELBDM )

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
   dt = 1/MaxdS_dx*0.5/ELBDM_ETA;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Phase



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

   const real Eps               = 0;//1.0e-2;          // soften in the denominator for calculating dS_dx
   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;
   const int  NGhost            = 1;               // number of ghost zones to calculate dS_dt
   const int  Size_Flu          = PS2 + 2*NGhost;  // size of the array Flu_Array
   const int  NPG               = 1;               // number of patch groups (must be ONE here)
   const int  AMP               = DENS;            // array index to store amplitude

   real (*Flu_Array)[NCOMP_FLUID][Size_Flu][Size_Flu][Size_Flu] = NULL;

   real dS_dx, MaxdS_dx, _dh, _dh2, _dh_sqr, Re, Im, GradS[3];
   real _Dens, _Amp;                         // 0.5/Dens/dh, 1.0/amplitude
   int  im, ip, jm, jp, km, kp, I, J, K;


   MaxdS_dx = 0.0;
   _dh      = (real)1.0/amr->dh[lv];
   _dh2     = (real)0.5*_dh;
   _dh_sqr  = _dh*_dh;


#  pragma omp parallel private( Flu_Array, dS_dx, Re, Im, GradS, _Dens, _Amp, \
                                im, ip, jm, jp, km, kp, I, J, K )
   {
      Flu_Array = new real [NPG][NCOMP_FLUID][Size_Flu][Size_Flu][Size_Flu];

//    loop over all patches
#     pragma omp for reduction( max:MaxdS_dx ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)
      {
//       prepare real part with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][REAL][0][0][0], NULL, NGhost, NPG, &PID0, _REAL, _NONE,
                            INT_MINMOD1D, INT_NONE, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

//       prepare imag part with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][IMAG][0][0][0], NULL, NGhost, NPG, &PID0, _IMAG, _NONE,
                            INT_MINMOD1D, INT_NONE, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


//       evaluate amplitude with soften
         for (int k=0; k<Size_Flu; k++)
         for (int j=0; j<Size_Flu; j++)
         for (int i=0; i<Size_Flu; i++)
         {
            Re = Flu_Array[0][REAL][k][j][i];
            Im = Flu_Array[0][IMAG][k][j][i];

            Flu_Array[0][AMP][k][j][i] = SQRT( Re*Re + Im*Im + Eps );
         }

//       evaluate dS_dt
         for (int k=NGhost; k<Size_Flu-NGhost; k++)    {  km = k - 1;    kp = k + 1;   K = k - NGhost;
         for (int j=NGhost; j<Size_Flu-NGhost; j++)    {  jm = j - 1;    jp = j + 1;   J = j - NGhost;
         for (int i=NGhost; i<Size_Flu-NGhost; i++)    {  im = i - 1;    ip = i + 1;   I = i - NGhost;

            _Amp       = (real)1.0 / Flu_Array[0][AMP][k][j][i];
            _Dens      = _dh2*_Amp*_Amp;

            GradS[0]   = _Dens*(  Flu_Array[0][REAL][k][j][i]*( Flu_Array[0][IMAG][k ][j ][ip] -
                                                                Flu_Array[0][IMAG][k ][j ][im] )
                                 -Flu_Array[0][IMAG][k][j][i]*( Flu_Array[0][REAL][k ][j ][ip] -
                                                                Flu_Array[0][REAL][k ][j ][im] )  );
            GradS[1]   = _Dens*(  Flu_Array[0][REAL][k][j][i]*( Flu_Array[0][IMAG][k ][jp][i ] -
                                                                Flu_Array[0][IMAG][k ][jm][i ] )
                                 -Flu_Array[0][IMAG][k][j][i]*( Flu_Array[0][REAL][k ][jp][i ] -
                                                                Flu_Array[0][REAL][k ][jm][i ] )  );
            GradS[2]   = _Dens*(  Flu_Array[0][REAL][k][j][i]*( Flu_Array[0][IMAG][kp][j ][i ] -
                                                                Flu_Array[0][IMAG][km][j ][i ] )
                                 -Flu_Array[0][IMAG][k][j][i]*( Flu_Array[0][REAL][kp][j ][i ] -
                                                              Flu_Array[0][REAL][km][j ][i ] )  );

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



#endif // #if ( MODEL == ELBDM )
