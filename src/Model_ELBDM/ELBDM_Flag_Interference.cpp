#include "GAMER.h"

#if ( MODEL == ELBDM )

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Flag_Interference
// Description :  Flag according to the interference criterion
//
// Note        :  1. Flag the input cell if the interference criteria are met (Interference, Phase jumps and additional check whether dB wavelength is resolved)
//                2. Size of the input array "Cond_Array" should be PATCH_SIZE^3
//
// Parameter   :  i,j,k       : Indices of the target cell in the arrays "Cond_Array"
//                Threshold   : Refinement Threshold for quantum pressure
//
// Return      :  "true"  if the flag criterion is     fulfilled
//                "false" if the flag criterion is NOT fulfilled
//-------------------------------------------------------------------------------------------------------
bool ELBDM_Flag_Interference( const int i, const int j, const int k, const real Cond_Array[], const double Threshold)
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   const int Idx = k*PS1*PS1 + j*PS1 + i;

   return ( Cond_Array[Idx] > Threshold );
} // FUNCTION : ELBDM_Flag_Interference




//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_for_Interference_Criterion
// Description :  Evaluate quantum pressure, phase jumps and whether dB wavelength will be resolved after refinement for the interference criterion
//
// Note        :  1. This function is called in "Flag_Real" before looping over all cells in the patch in order to
//                   achieve higher performance
//                2. Evaluate laplacian with second-order stencil and dB wavelength via first-order stencil
//                3. Do not take into account the physical size of each cell since criteria are dimensionless
//                4. The sizes of the arrays (Var1D, Temp1D, Cond1D) must be ( (PS1+2)^3, 2*(PS1+2)^3, 3*(PS1)^3 )
//
// Parameter   :  Var       : Array storing the input density and phase for UseWaveFlag == false and the real and imaginary parts for UseWaveFlag == true
//                Temp      : Array to store the intermediate variable sqrt(density) and the phase
//                Cond      : Array to store the output dimensionless quantum pressure, the curvature of the phase field as well as
//                            the maximum phase difference between neighbouring points
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Prepare_for_Interference_Criterion(const real *Var1D, real *Temp1D, real *Cond1D, const bool UseWaveFlag)
{

   const int NCell  = PS1 + 2;   // size of the arrays Var, Temp
   const int NCond  = PS1;       // size of the array  Cond

   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;

// convert the 1D arrays
   real (*Var)  [NCell][NCell][NCell] = ( real(*) [NCell][NCell][NCell] )  Var1D;
   real (*Temp) [NCell][NCell][NCell] = ( real(*) [NCell][NCell][NCell] )  Temp1D;
   real (*Cond) [NCond][NCond][NCond] = ( real(*) [NCond][NCond][NCond] )  Cond1D;

   for (int k=0; k<NCell; k++)    {
   for (int j=0; j<NCell; j++)    {
   for (int i=0; i<NCell; i++)    {
//    compute square root of density field
      Temp[0][k][j][i] = Var[0][k][j][i];

//    check whether to convert imaginary and real parts to the phase
      if ( UseWaveFlag ) {
         Temp[1][k][j][i] = SATAN2(Var[IMAG][k][j][i], Var[REAL][k][j][i]);
      } else {
         Temp[1][k][j][i] = Var[1][k][j][i];
      }
   }}} // k,j,i

   for (int k=0; k<NCond; k++)    {  kk = k + 1;   kkp = kk + 1;   kkm = kk - 1;
   for (int j=0; j<NCond; j++)    {  jj = j + 1;   jjp = jj + 1;   jjm = jj - 1;
   for (int i=0; i<NCond; i++)    {  ii = i + 1;   iip = ii + 1;   iim = ii - 1;

//    compute the dimensionless quantum pressure (divided by number of dimensions for normalisation) if phase and density fields have local extrema
      const bool DChangeSignX = ( Temp[0][kk ][jj ][iip] - Temp[0][kk][jj][ii] ) * ( Temp[0][kk][jj][ii] - Temp[0][kk ][jj ][iim] ) < 0;
      const bool DChangeSignY = ( Temp[0][kk ][jjp][ii ] - Temp[0][kk][jj][ii] ) * ( Temp[0][kk][jj][ii] - Temp[0][kk ][jjm][ii ] ) < 0;
      const bool DChangeSignZ = ( Temp[0][kkp][jj ][ii ] - Temp[0][kk][jj][ii] ) * ( Temp[0][kk][jj][ii] - Temp[0][kkm][jj ][ii ] ) < 0;
      const bool SChangeSignX = ( Temp[1][kk ][jj ][iip] - Temp[1][kk][jj][ii] ) * ( Temp[1][kk][jj][ii] - Temp[1][kk ][jj ][iim] ) < 0;
      const bool SChangeSignY = ( Temp[1][kk ][jjp][ii ] - Temp[1][kk][jj][ii] ) * ( Temp[1][kk][jj][ii] - Temp[1][kk ][jjm][ii ] ) < 0;
      const bool SChangeSignZ = ( Temp[1][kkp][jj ][ii ] - Temp[1][kk][jj][ii] ) * ( Temp[1][kk][jj][ii] - Temp[1][kkm][jj ][ii ] ) < 0;
      //const bool SCurvX     =  FABS( Temp[1][kk ][jj ][iip] - 2 * Temp[1][kk][jj][ii] + Temp[1][kk ][jj ][iim] );
      //const bool SCurvY     =  FABS( Temp[1][kk ][jjp][ii ] - 2 * Temp[1][kk][jj][ii] + Temp[1][kk ][jjm][ii ] );
      //const bool SCurvZ     =  FABS( Temp[1][kkp][jj ][ii ] - 2 * Temp[1][kk][jj][ii] + Temp[1][kkm][jj ][ii ] );
      //const real D          =  Temp[0][kk][jj][ii];
      Cond[0][k][j][i] =  0;

      if ( SChangeSignX && DChangeSignX )
         Cond[0][k][j][i] +=  FABS( Temp[0][kk ][jj ][iip] - 2 * Temp[0][kk ][jj ][ii] + Temp[0][kk ][jj ][iim] ) / Temp[0][kk][jj][ii];
      if ( SChangeSignY && DChangeSignY )
         Cond[0][k][j][i] +=  FABS( Temp[0][kk ][jjp][ii ] - 2 * Temp[0][kk ][jj ][ii] + Temp[0][kk ][jjm][ii ] ) / Temp[0][kk][jj][ii];
      if ( SChangeSignZ && DChangeSignZ )
         Cond[0][k][j][i] +=  FABS( Temp[0][kkp][jj ][ii ] - 2 * Temp[0][kk ][jj ][ii] + Temp[0][kkm][jj ][ii ] ) / Temp[0][kk][jj][ii];

//    flags for wave solver
      if (UseWaveFlag)   {
         Cond[1][k][j][i] = 0;

//    flags for fluid solver
      } else { // if (UseWaveFlag)

//       check second derivative of phase field (divided by number of dimensions for normalisation) to detect phase jumps
         Cond[1][k][j][i] =  FABS( Temp[1][kk ][jj ][iip] + Temp[1][kk ][jj ][iim]
                                 + Temp[1][kk ][jjp][ii ] + Temp[1][kk ][jjm][ii ]
                                 + Temp[1][kkp][jj ][ii ] + Temp[1][kkm][jj ][ii ]
                                 -  (real) 6.0 * Temp[1][kk][jj][ii])
                                 / ((real) 3.0);

      }  // if (UseWaveFlag) ... else

//    resolve de Broglie wavelength when vortices are absent
//    use quantum pressure combined with vortex threshold to determine whether vortex is present
      if ( Cond[0][k][j][i] < ELBDM_VORTEX_THRESHOLD ) {
         Cond[2][k][j][i] = FMAX(FMAX(FMAX(FMAX(FMAX(
            FABS(Temp[1][kk ][jj ][iip] - ELBDM_UnwrapPhase(Temp[1][kk ][jj ][iip], Temp[1][kk ][jj ][ii ])),
            FABS(Temp[1][kk ][jj ][ii ] - ELBDM_UnwrapPhase(Temp[1][kk ][jj ][ii ], Temp[1][kk ][jj ][iim]))),
            FABS(Temp[1][kk ][jjp][ii ] - ELBDM_UnwrapPhase(Temp[1][kk ][jjp][ii ], Temp[1][kk ][jj ][ii ]))),
            FABS(Temp[1][kk ][jj ][ii ] - ELBDM_UnwrapPhase(Temp[1][kk ][jj ][ii ], Temp[1][kk ][jjm][ii ]))),
            FABS(Temp[1][kkp][jj ][ii ] - ELBDM_UnwrapPhase(Temp[1][kkp][jj ][ii ], Temp[1][kk ][jj ][ii ]))),
            FABS(Temp[1][kk ][jj ][ii ] - ELBDM_UnwrapPhase(Temp[1][kk ][jj ][ii ], Temp[1][kkm][jj ][ii ]))) / (2 * M_PI);
      } else {
//       around the vortex, the velocity and the de Broglie wavelength diverge
//       disable de Broglie wavelength refinement there
         Cond[2][k][j][i] = 0;
      }
   }}} // k,j,i
} // FUNCTION : Prepare_for_Interference_Criterion



#endif // #if ( MODEL == ELBDM )
