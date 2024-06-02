#include "GAMER.h"

#if ( MODEL == ELBDM )

static real GetMaxPhaseDerivative( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Phase
// Description :  Estimate the evolution time-step by restricting the maximum phase rotation
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. dt = DT__PHASE*2*pi/Max(first temporal derivative of phase)
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Phase( const int lv )
{

   real   MaxdS_dt;
   double dt;

// get the maximum first temporal derivative of phase
   MaxdS_dt = GetMaxPhaseDerivative( lv );

// get the time-step
   dt = DT__PHASE*2.0*M_PI/MaxdS_dt;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Phase



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxPhaseDerivative
// Description :  Evaluate the maximum first temporal derivative of phase for the time-step estimation
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Phase()
//                2. Currently this function does not take into account the self-interaction potential
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MaxdS_dt
//-------------------------------------------------------------------------------------------------------
real GetMaxPhaseDerivative( const int lv )
{

   const real Eps               = 1.0e-2;          // soften in the denominator for calculating dS_dt
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

   real dS_dt, MaxdS_dt, _dh, _dh2, _dh_sqr, Re, Im, GradS[3];
   real LapAmp_Amp, Vel_Sqr;                 // laplacian(amplitude)/amplitude, -grad(phase)^2
   real _Dens, _Amp;                         // 0.5/Dens/dh, 1.0/amplitude
   int  im, ip, jm, jp, km, kp, I, J, K;

#  ifdef GRAVITY
   const real PotCoeff = -2.0*SQR(ELBDM_ETA);
   const int  Size_Pot = PS2;                // size of the array Pot_Array
   real Pot;                                 // -2.0*ELBDM_ETA^2*potential
   real (*Pot_Array)[Size_Pot][Size_Pot][Size_Pot] = NULL;
#  endif

   MaxdS_dt = 0.0;
   _dh      = (real)1.0/amr->dh[lv];
   _dh2     = (real)0.5*_dh;
   _dh_sqr  = _dh*_dh;


#  ifdef GRAVITY
#  pragma omp parallel private( Flu_Array, dS_dt, Re, Im, GradS, LapAmp_Amp, Vel_Sqr, _Dens, _Amp, Pot_Array, Pot, \
                                im, ip, jm, jp, km, kp, I, J, K )
#  else
#  pragma omp parallel private( Flu_Array, dS_dt, Re, Im, GradS, LapAmp_Amp, Vel_Sqr, _Dens, _Amp, \
                                im, ip, jm, jp, km, kp, I, J, K )
#  endif
   {
      Flu_Array = new real [NPG][NCOMP_FLUID][Size_Flu][Size_Flu][Size_Flu];
#     ifdef GRAVITY
      Pot_Array = new real [NPG]             [Size_Pot][Size_Pot][Size_Pot];
#     endif

//    loop over all patches
#     pragma omp for reduction( max:MaxdS_dt ) schedule( runtime )
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

//       prepare potential with no ghost zone
#        ifdef GRAVITY
         Prepare_PatchData( lv, Time[lv], &Pot_Array[0][0][0][0],       NULL,      0, NPG, &PID0, _POTE, _NONE,
                            INT_NONE,     INT_NONE, UNIT_PATCHGROUP, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
#        endif

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

            LapAmp_Amp = _dh_sqr*( Flu_Array[0][AMP][kp][j ][i ] + Flu_Array[0][AMP][km][j ][i ] +
                                   Flu_Array[0][AMP][k ][jp][i ] + Flu_Array[0][AMP][k ][jm][i ] +
                                   Flu_Array[0][AMP][k ][j ][ip] + Flu_Array[0][AMP][k ][j ][im] -
                                   (real)6.0*Flu_Array[0][AMP][k][j][i] ) * _Amp;
            Vel_Sqr    = -( GradS[0]*GradS[0] + GradS[1]*GradS[1] + GradS[2]*GradS[2] );
            dS_dt      = LapAmp_Amp + Vel_Sqr;

#           ifdef GRAVITY
            Pot        = PotCoeff*Pot_Array[0][K][J][I];
            dS_dt     += Pot;
#           endif

            dS_dt      = FABS( dS_dt );
            MaxdS_dt   = MAX( MaxdS_dt, dS_dt );

         }}} // i,j,k
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)

      delete [] Flu_Array;
#     ifdef GRAVITY
      delete [] Pot_Array;
#     endif
   } // OpenMP parallel region


// get the maximum potential in all ranks
   real MaxdS_dt_AllRank;
   MPI_Allreduce( &MaxdS_dt, &MaxdS_dt_AllRank, 1, MPI_GAMER_REAL, MPI_MAX, MPI_COMM_WORLD );


// check
   if ( MaxdS_dt_AllRank == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MaxdS_dt == 0.0 at lv %d !!\n", lv );


   return MaxdS_dt_AllRank*0.5/ELBDM_ETA;

} // FUNCTION : GetMaxPhaseDerivative



#endif // #if ( MODEL == ELBDM )
