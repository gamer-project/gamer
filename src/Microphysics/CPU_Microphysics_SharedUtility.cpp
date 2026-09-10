#ifndef __CUFLU_MICROSHARED__
#define __CUFLU_MICROSHARED__

#include "CUFLU.h"

#if defined( VISCOSITY ) || defined( CONDUCTION ) || defined( CR_DIFFUSION )

#ifdef __CUDACC__
#include "../Model_Hydro/GPU_Hydro/CUFLU_Shared_FluUtility.cu"
#endif

//-----------------------------------------------------------------------------------------
// Function    : minmod
// Description : Minmod slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real minmod( const real a, const real b )
{

   if      ( a > (real)0.0  &&  b > (real)0.0 )    return FMIN(a, b);
   else if ( a < (real)0.0  &&  b < (real)0.0 )    return FMAX(a, b);
   else                                            return (real)0.0;

} // FUNCTION : minmod

//-----------------------------------------------------------------------------------------
// Function    : MC_limiter
// Description : Monotonized central (MC) slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real MC_limiter( const real a, const real b )
{

   return minmod( (real)2.0*minmod(a,b), (real)0.5*(a+b) );

} // FUNCTION : MC_limiter

//-----------------------------------------------------------------------------------------
// Function    : van_leer 
// Description : vanLeer slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real van_leer( const real a, const real b )
{

   return ( a*b > (real)0.0 ) ? ( (real)2.0*a*b / (a+b) ) : (real)0.0;

} // FUNCTION : van_leer

//-----------------------------------------------------------------------------------------
// Function    : van_leer2 
// Description : vanLeer slope limiter called twice
//
// Parameter   : a, b, c, d : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real van_leer2( const real a, const real b, const real c, const real d )
{

   return van_leer( van_leer(a,b), van_leer(c,d) );

} // FUNCTION : van_leer2

GPU_DEVICE
real compute_temperature( const real ConVar[][ CUBE(FLU_NXT) ],
                          const real PriVar[][ CUBE(FLU_NXT) ],
                          const real   FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                          const int idx, const real MinTemp,
                          const long PassiveFloor, const EoS_t *EoS )
{
   const bool CheckMinTemp_Yes = true;
   real temp;
   real fluid[NCOMP_TOTAL];

   if ( PriVar == NULL )
   {
      real Emag;
#ifdef MHD
      const int size_ij = SQR( FLU_NXT );
      const int i       = idx % FLU_NXT;
      const int j       = idx % size_ij / FLU_NXT;
      const int k       = idx / size_ij;
      Emag = MHD_GetCellCenteredBEnergy( FC_B[0], FC_B[1], FC_B[2], FLU_NXT, FLU_NXT, FLU_NXT,
                                         i, j, k );
#else
      Emag = NULL_REAL;
#endif
      for (int v=0; v<NCOMP_TOTAL; v++) fluid[v] = ConVar[v][idx];
      temp = Hydro_Con2Temp( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ],
                             fluid[ENGY], fluid+NCOMP_FLUID, CheckMinTemp_Yes,
                             MinTemp, PassiveFloor, Emag, EoS->DensEint2Temp_FuncPtr,
                             EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr,
                             EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
   }
   else
   {
      for (int v=0; v<NCOMP_TOTAL; v++) fluid[v] = PriVar[v][idx];
      const real Eint = EoS->DensPres2Eint_FuncPtr( fluid[DENS], fluid[ENGY], fluid+NCOMP_FLUID,
                                                    EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                                                    EoS->Table );
      temp = EoS->DensEint2Temp_FuncPtr( fluid[DENS], Eint, fluid+NCOMP_FLUID, EoS->AuxArrayDevPtr_Flt,
                                         EoS->AuxArrayDevPtr_Int, EoS->Table );
      temp = Hydro_CheckMinTemp( temp, MinTemp );
   } // if ( PriVar == NULL ) ... else ...

   return temp;
} // FUNCTION : compute_temperature

#endif // #if defined( VISCOSITY ) || defined( CONDUCTION ) || defined( CR_DIFFUSION )

#endif // #ifndef __CUFLU_MICROSHARED__
