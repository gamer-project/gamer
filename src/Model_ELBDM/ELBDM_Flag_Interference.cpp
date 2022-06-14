#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Flag_Interference
// Description :  Flag according to the interference criterion 
//
// Note        :  1. Flag the input cell if the level of interference as quantified by 
//                   M_int = dh^2 / dim * laplace sqrt(rho) / sqrt(rho)
//                2. Size of the input array "Cond_Array" should be PATCH_SIZE^3 
//
// Parameter   :  i,j,k       : Indices of the target cell in the arrays "Cond_Array"
//                Theshold    : Refinement Threshold for M_int
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

   const int    Idx       = k*PS1*PS1 + j*PS1 + i;
   const real   InterferenceCriterion   = Cond_Array[Idx];
   bool Flag = InterferenceCriterion > Threshold;

   return Flag;
} // FUNCTION : ELBDM_Flag_EngyDensity




//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_for_Interference_Criterion
// Description :  Evaluate quantum pressure for the interference criterium
//
// Note        :  1. This function is called in "Flag_Real" before looping over all cells in the patch in order to
//                   achieve higher performance
//                2. Evaluate laplacian with second-order stencil
//                3. Do not take into account the physical size of each cell since the Lohner error estimator
//                   is dimensionless
//                4. The sizes of the arrays (Dens1D, SqrtDens1D, QP1D) must be ( (PS1+2)^3, (PS1+2)^3, (PS1)^3 )
//
// Parameter   :  Dens1D       : Array storing the input density for the interference criterion
//                SqrtDens1D   : Array to store the intermediate variable sqrt(density)
//                QP1D         : Array to store the output dimensionless quantum pressures for the interference criterion
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Prepare_for_Interference_Criterion(const real *Dens1D, real *SqrtDens1D, real *QP1D)
{

   const int NCell  = PS1 + 2;   // size of the arrays Dens1D, SqrtDens1D
   const int NCurv  = PS1;       // size of the array  QP1D

   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;

// convert the 1D arrays
   real (*Dens)         [NCell ][NCell ] = ( real(*) [NCell ][NCell ] )  Dens1D;
   real (*SqrtDens)     [NCell ][NCell ] = ( real(*) [NCell ][NCell ] )  SqrtDens1D;
   real (*QP)           [NCurv ][NCurv ] = ( real(*) [NCurv ][NCurv ] )  QP1D;

   for (int k=0; k<NCell; k++)    {
   for (int j=0; j<NCell; j++)    {
   for (int i=0; i<NCell; i++)    {
      SqrtDens[k][j][i] = Dens[k][j][i];
   }}} // k,j,i

   for (int k=0; k<NCurv; k++)    {  kk = k + 1;   kkp = kk + 1;   kkm = kk - 1;
   for (int j=0; j<NCurv; j++)    {  jj = j + 1;   jjp = jj + 1;   jjm = jj - 1;
   for (int i=0; i<NCurv; i++)    {  ii = i + 1;   iip = ii + 1;   iim = ii - 1;

      QP[k][j][i] =  FABS(  SqrtDens[kk ][jj ][iip] + SqrtDens[kk ][jj ][iim] \
                         + SqrtDens[kk ][jjp][ii ] + SqrtDens[kk ][jjm][ii ] \
                         + SqrtDens[kkp][jj ][ii ] + SqrtDens[kkm][jj ][ii ] \
                         - (real) 6.0 * SqrtDens[kk ][jj ][ii])\
                    / ((real) 3.0 * SqrtDens[kk ][jj ][ii]);
   }}} // k,j,i


} // FUNCTION : Prepare_for_Interference_Criterion*/



#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
