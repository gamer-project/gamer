#include "GAMER.h"

#if ( MODEL == ELBDM )

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Flag_Interference
// Description :  Flag according to the interference criterion
//
// Note        :  1. Flag the input cell if all the interference criteria are met (minimum density, local extrema in phase and density fields, quantum pressure, phase curvature)
//                2. Size of the input array "Var" should be 2*(PS1+2)^3
//
// Parameter   :  i,j,k             : Indices of the target cell in the array "Var1D"
//                Var               : Input array holding the density and phase field
//                QPThreshold       : Refinement Threshold for quantum pressure
//                                    Refines when dimensionless quantum pressure / NDim exceeds QPThreshold
//                                    Ensure stability and correctness of hybrid scheme by refining regions where fluid approach produces high errors and fails to wave scheme
//                                    QPThreshold <= 0.03 avoids spurious halos and yields good agreement with wave-only simulations
//                DensThreshold     : Minimum density at which to check quantum pressure threshold
//                                    Should be set to zero by default, but can be used to avoid refinement of low-density regions with small density oscillations
//                                    Use values > 0 with care since they may lead to instability.                                    Use values > 0 with care since they may lead to instability.
//                LapPhaseThreshold : Refinement Threshold for second derivative of phase field in addition to QPThreshold
//                                    Should be set to zero by default, but can be used to avoid refinement of of regions with high quantum pressure without destructive interference
//                                    The motivation behind this flag is that destructive interference coincides with a high second derivative of the phase field
//                                    Use with care since it has been shown to lead to instability in some cases
//                OnlyAtExtrema     : Boolean flag indicating whether only extrema are refined
//                                    Should be set to False, but can be used to avoid refinement of regions with high quantum pressure without destructive interference
//                                    The motivation behind this flag is that destructive interference where the fluid scheme fails happens at extrema of the density and phase fields
//                                    Use with care since it has been shown to lead to instability in some cases
//
// Return      :  "true"  if the flag criterion is     fulfilled
//                "false" if the flag criterion is NOT fulfilled
//-------------------------------------------------------------------------------------------------------
bool ELBDM_Flag_Interference( const int i, const int j, const int k, const real Var1D[], const double QPThreshold, const double DensThreshold, const double LapPhaseThreshold, const bool OnlyAtExtrema )
{
// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   const int NGhost = 1;
   const int NCell  = PS1 + 2 * NGhost;   // size of the array Var

   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;

// convert the 1D array
   real (*Var)  [NCell][NCell][NCell] = ( real(*) [NCell][NCell][NCell] )  Var1D;

   kk = k + NGhost;   kkp = kk + 1;   kkm = kk - 1;
   jj = j + NGhost;   jjp = jj + 1;   jjm = jj - 1;
   ii = i + NGhost;   iip = ii + 1;   iim = ii - 1;

// check minimum density
   const bool DensCond = Var[0][kk][jj][ii] > DensThreshold;

   if ( !DensCond ) return false;

// check whether density and phase fields have local extrema
   const bool DChangeSignX = ( !OnlyAtExtrema ) || (( Var[0][kk ][jj ][iip] - Var[0][kk][jj][ii] ) * ( Var[0][kk][jj][ii] - Var[0][kk ][jj ][iim] ) < 0.0);
   const bool DChangeSignY = ( !OnlyAtExtrema ) || (( Var[0][kk ][jjp][ii ] - Var[0][kk][jj][ii] ) * ( Var[0][kk][jj][ii] - Var[0][kk ][jjm][ii ] ) < 0.0);
   const bool DChangeSignZ = ( !OnlyAtExtrema ) || (( Var[0][kkp][jj ][ii ] - Var[0][kk][jj][ii] ) * ( Var[0][kk][jj][ii] - Var[0][kkm][jj ][ii ] ) < 0.0);
   const bool SChangeSignX = ( !OnlyAtExtrema ) || (( Var[1][kk ][jj ][iip] - Var[1][kk][jj][ii] ) * ( Var[1][kk][jj][ii] - Var[1][kk ][jj ][iim] ) < 0.0);
   const bool SChangeSignY = ( !OnlyAtExtrema ) || (( Var[1][kk ][jjp][ii ] - Var[1][kk][jj][ii] ) * ( Var[1][kk][jj][ii] - Var[1][kk ][jjm][ii ] ) < 0.0);
   const bool SChangeSignZ = ( !OnlyAtExtrema ) || (( Var[1][kkp][jj ][ii ] - Var[1][kk][jj][ii] ) * ( Var[1][kk][jj][ii] - Var[1][kkm][jj ][ii ] ) < 0.0);

// compute second derivative of phase field
   const bool SCurvX       =  FABS( Var[1][kk ][jj ][iip] - 2 * Var[1][kk][jj][ii] + Var[1][kk ][jj ][iim] ) > LapPhaseThreshold;
   const bool SCurvY       =  FABS( Var[1][kk ][jjp][ii ] - 2 * Var[1][kk][jj][ii] + Var[1][kk ][jjm][ii ] ) > LapPhaseThreshold;
   const bool SCurvZ       =  FABS( Var[1][kkp][jj ][ii ] - 2 * Var[1][kk][jj][ii] + Var[1][kkm][jj ][ii ] ) > LapPhaseThreshold;

   real QP                 =  0;
   const real SqrtRhoC     =  SQRT(Var[0][kk][jj][ii]);

   if ( SChangeSignX && DChangeSignX )
   {
      QP +=  FABS( SQRT(Var[0][kk ][jj ][iip]) - 2 * SqrtRhoC + SQRT(Var[0][kk ][jj ][iim]) ) / SqrtRhoC;
   }
   if ( SChangeSignY && DChangeSignY )
   {
      QP +=  FABS( SQRT(Var[0][kk ][jjp][ii ]) - 2 * SqrtRhoC + SQRT(Var[0][kk ][jjm][ii ]) ) / SqrtRhoC;
   }
   if ( SChangeSignZ && DChangeSignZ )
   {
      QP +=  FABS( SQRT(Var[0][kkp][jj ][ii ]) - 2 * SqrtRhoC + SQRT(Var[0][kkm][jj ][ii ]) ) / SqrtRhoC;
   }

   return ( QP > QPThreshold ) || SCurvX || SCurvY || SCurvZ;
} // FUNCTION : ELBDM_Flag_Interference

#endif // #if ( MODEL == ELBDM )
