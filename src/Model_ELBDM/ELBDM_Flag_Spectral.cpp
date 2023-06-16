#include "GAMER.h"

#if ( MODEL == ELBDM && WAVE_SCHEME == WAVE_GRAMFE )

#include <complex.h>

typedef std::complex<real> complex_type;


real Compute_Extension_Mass(const complex_type Row[FLU_NXT]) {
   real extensionMass = 0;
   for (size_t i = FLU_NXT; i < GRAMFE_FLU_NXT; ++i)
   {
      for (size_t j = 0; j < FLU_NXT; ++j)
      {
         psi = {0, 0};

         for (int order=0; order < GRAMFE_ORDER; order++) {
            psi += (complex_type) Gramfe_Extend[i][j] *  Row[j];
         }

         extensionMass += SQR(psi.real()) + SQR(psi.imag());
      }
   }
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_for_Spectral_Criterio
// Description :  Evaluate ratio of mass of wave function in patch group and mass of wave function in extension domain
//
// Note        :  1. This function is called once per patch group
//                4. The size of the arrays Var1D must be FLU_NXT^3
//
// Parameter   :  Var1D     : Array storing the input re & im
//                Cond      : Reference to floating point variable where mass ratio is stored
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Prepare_for_Spectral_Criterion(const real *Var1D, real& Cond1D)
{

   const int NCell  = PS1 + 2;   // size of the arrays Var, Temp
   const int NCond  = PS1;       // size of the array  Cond

   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;

// convert the 1D arrays
   real (*Var)  [NCell ][NCell ][NCell ] = ( real(*) [NCell ][NCell ][NCell ] )  Var1D;
   real  Row[FLU_NXT];

   real PhysicalMass = 0, ExtensionMass = 0;

   for (int k=0; k<NCell; k++)    {
   for (int j=0; j<NCell; j++)    {
   for (int i=0; i<NCell; i++)    {
      PhysicalMass += SQR(Var[REAL][k][j][i]) + SQR(Var[IMAG][k][j][i]);
   }}} // k,j,i

// x-direction
   for (int k=0; k<NCell; k++)    {
   for (int j=0; j<NCell; j++)    {
      for (int i=0; i<NCell; i++)    {
         Row[i] = complex_type(Var[REAL][k][j][i], Var[IMAG][k][j][i]);
      }
      ExtensionMass += Compute_Extension_Mass(Row);
   }} // k,j,i

// y-direction
   for (int k=0; k<NCell; k++)    {
   for (int j=0; j<NCell; j++)    {
      for (int i=0; i<NCell; i++)    {
         Row[i] = complex_type(Var[REAL][k][i][j], Var[IMAG][k][i][j]);
      }
      ExtensionMass += Compute_Extension_Mass(Row);
   }} // k,j,i


// z-direction
   for (int k=0; k<NCell; k++)    {
   for (int j=0; j<NCell; j++)    {
      for (int i=0; i<NCell; i++)    {
         Row[i] = complex_type(Var[REAL][i][j][k], Var[IMAG][i][j][k]);
      }
      ExtensionMass += Compute_Extension_Mass(Row);
   }} // k,j,i

   Cond1D = ExtensionMass / PhysicalMass;

} // FUNCTION : Prepare_for_Spectral_Criterion



#endif // #if ( MODEL == ELBDM && WAVE_SCHEME == WAVE_GRAMFE )
