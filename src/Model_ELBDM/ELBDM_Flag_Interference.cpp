#include "GAMER.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Flag_Interference
// Description :  Flag according to the interference criterion
//
// Note        :  1. Flag the input cell if one of the following criteria is met:
//                   - minimum density + quantum pressure + local extremum in density field
//                   - second derivative in phase + local extremum in phase field
//                2. Size of the input array "Var1D" should be 2*(PS1+2)^3
//
// Parameter   :  i,j,k             : Indices of the target cell in the array "Var1D"
//                Var1D             : Input array holding the density and phase field
//                QPThreshold       : Refinement threshold for quantum pressure
//                                    Refines when dimensionless quantum pressure in any dimension exceeds QPThreshold + (DensThreshold and OnlyAtExtrema if set)
//                                    Ensure stability and correctness of hybrid scheme by refining regions where fluid approach produces high errors and fails to wave scheme
//                                    QPThreshold <= 0.03 avoids spurious halos and yields good agreement with wave-only simulations
//                DensThreshold     : Minimum density at which to check quantum pressure threshold
//                                    Should be set to zero by default, but can be used to avoid refinement of low-density regions with small density oscillations
//                                    Use values > 0 with care since they may lead to instability.
//                LapPhaseThreshold : Refinement threshold for second derivative of phase field
//                                    Refines when dimensionless second derivative of phase in any dimension exceeds LapPhaseThreshold
//                OnlyAtExtrema     : Boolean flag indicating whether only extrema are refined
//                                    Should be set to False, but can be used to avoid refinement of regions with high quantum pressure and phase curvature without destructive interference
//                                    The motivation behind this flag is that destructive interference where the fluid scheme fails happens at extrema of the density and phase fields
//                                    Use with care since it has been shown to lead to instability in some cases
//
// Return      :  "true"  if the flag criterion is     fulfilled
//                "false" if the flag criterion is NOT fulfilled
//-------------------------------------------------------------------------------------------------------
bool ELBDM_Flag_Interference( const int i, const int j, const int k, const real Var1D[], const double QPThreshold,
                              const double DensThreshold, const double LapPhaseThreshold, const bool OnlyAtExtrema )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif


   const int NGhost = 1;
   const int NCell  = PS1 + 2*NGhost;  // size of the array Var

   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;


// convert the 1D array
   real (*Var)[NCell][NCell][NCell] = ( real(*)[NCell][NCell][NCell] )Var1D;

   kk = k + NGhost;   kkp = kk + 1;   kkm = kk - 1;
   jj = j + NGhost;   jjp = jj + 1;   jjm = jj - 1;
   ii = i + NGhost;   iip = ii + 1;   iim = ii - 1;


// check minimum density
   const bool DensCond = ( Var[0][kk][jj][ii] > DensThreshold );


// check whether density and phase fields have local extrema
   const bool DensChangeSignX  = ( !OnlyAtExtrema )  ||  ( ( Var[0][kk ][jj ][iip] - Var[0][kk][jj][ii] ) * ( Var[0][kk][jj][ii] - Var[0][kk ][jj ][iim] ) < 0.0 );
   const bool DensChangeSignY  = ( !OnlyAtExtrema )  ||  ( ( Var[0][kk ][jjp][ii ] - Var[0][kk][jj][ii] ) * ( Var[0][kk][jj][ii] - Var[0][kk ][jjm][ii ] ) < 0.0 );
   const bool DensChangeSignZ  = ( !OnlyAtExtrema )  ||  ( ( Var[0][kkp][jj ][ii ] - Var[0][kk][jj][ii] ) * ( Var[0][kk][jj][ii] - Var[0][kkm][jj ][ii ] ) < 0.0 );
   const bool PhaseChangeSignX = ( !OnlyAtExtrema )  ||  ( ( Var[1][kk ][jj ][iip] - Var[1][kk][jj][ii] ) * ( Var[1][kk][jj][ii] - Var[1][kk ][jj ][iim] ) < 0.0 );
   const bool PhaseChangeSignY = ( !OnlyAtExtrema )  ||  ( ( Var[1][kk ][jjp][ii ] - Var[1][kk][jj][ii] ) * ( Var[1][kk][jj][ii] - Var[1][kk ][jjm][ii ] ) < 0.0 );
   const bool PhaseChangeSignZ = ( !OnlyAtExtrema )  ||  ( ( Var[1][kkp][jj ][ii ] - Var[1][kk][jj][ii] ) * ( Var[1][kk][jj][ii] - Var[1][kkm][jj ][ii ] ) < 0.0 );


// compute second derivative of phase field
   const bool LapPhaseX = PhaseChangeSignX  &&  ( FABS( Var[1][kk ][jj ][iip] - (real)2.0*Var[1][kk][jj][ii] + Var[1][kk ][jj ][iim] ) > LapPhaseThreshold );
   const bool LapPhaseY = PhaseChangeSignY  &&  ( FABS( Var[1][kk ][jjp][ii ] - (real)2.0*Var[1][kk][jj][ii] + Var[1][kk ][jjm][ii ] ) > LapPhaseThreshold );
   const bool LapPhaseZ = PhaseChangeSignZ  &&  ( FABS( Var[1][kkp][jj ][ii ] - (real)2.0*Var[1][kk][jj][ii] + Var[1][kkm][jj ][ii ] ) > LapPhaseThreshold );


// compute quantum pressure
   const real SqrtRhoC = SQRT(  MAX( Var[0][kk][jj][ii], TINY_NUMBER )  );

   const bool QPX = DensCond  &&  DensChangeSignX  &&  ( FABS( SQRT(Var[0][kk ][jj ][iip]) - (real)2.0*SqrtRhoC + SQRT(Var[0][kk ][jj ][iim]) ) / SqrtRhoC > QPThreshold );
   const bool QPY = DensCond  &&  DensChangeSignY  &&  ( FABS( SQRT(Var[0][kk ][jjp][ii ]) - (real)2.0*SqrtRhoC + SQRT(Var[0][kk ][jjm][ii ]) ) / SqrtRhoC > QPThreshold );
   const bool QPZ = DensCond  &&  DensChangeSignZ  &&  ( FABS( SQRT(Var[0][kkp][jj ][ii ]) - (real)2.0*SqrtRhoC + SQRT(Var[0][kkm][jj ][ii ]) ) / SqrtRhoC > QPThreshold );


   return ( QPX || QPY || QPZ || LapPhaseX || LapPhaseY || LapPhaseZ );

} // FUNCTION : ELBDM_Flag_Interference



#endif // #if ( MODEL == ELBDM )
