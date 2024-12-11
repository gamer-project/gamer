#include "GAMER.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_SetTaylor3Coeff
// Description :  Optimize the coefficient of the third-order term in the Taylor expansion
//
// Note        :  1. The optimization aims at minimizing the amplitude error for the smallest wavelength (k = kmax in 1D)
//                2. If the routine fails to find an optimized coefficient (which may be due to too large dt), it is
//                   set to 1.0/6.0 by default
//                3. Invoked by "CPU_FluidSolver" and "CUAPI_Asyn_FluidSolver"
//
// Parameter   :  dt    : Time interval to advance solution
//                dh    : Grid size
//                Eta   : Particle mass / Planck constant
//
// Return      :  Taylor3_Coeff
//-------------------------------------------------------------------------------------------------------
real ELBDM_SetTaylor3Coeff( const real dt, const real dh, const real Eta )
{

// check
   if ( dt <= 0.0 )  Aux_Error( ERROR_INFO, "dt = %14.7e <= 0.0 !!\n", dt );


   const real Alpha   = 0.5*dt/(Eta*dh*dh);
   const real Alpha2  = SQR( Alpha );
   const real Damping = 1.000157;   // dampling coefficient, which must be greater than 1

   real Temp, Taylor3_Coeff;

#  ifdef LAPLACIAN_4TH
   Temp = 81.0 - 576.0*Alpha2;

   if ( Temp >= 0.0 )
      Taylor3_Coeff = Damping * ( 9.0 - SQRT(Temp) ) / ( 256.0*Alpha2 );

   else
   {
      Taylor3_Coeff = real(1.0/6.0);

#     ifdef GAMER_DEBUG
      Aux_Message( stderr, "********************************************************************************\n" );
      Aux_Message( stderr, "WARNING : 81.0-576.0*Alpha2 = %12.5e < 0.0 (dt %12.5e, dh %12.5e, Eta %12.5e) !!\n",
                   Temp, dt, dh, Eta, __FUNCTION__ );
      Aux_Message( stderr, "        --> Program cannot find an optimized Taylor expansion coefficient.\n" );
      Aux_Message( stderr, "            Therefore, it is set to %12.5e by default.\n", Taylor3_Coeff );
      Aux_Message( stderr, "        --> Perhaps you should set DT__FLUID/DT__FLUID_INIT smaller ...\n" );
      Aux_Message( stderr, "        Rank <%d>, file <%s>, line <%d>, function <%s>\n",
                   MPI_Rank, __FILE__, __LINE__, __FUNCTION__ );
      Aux_Message( stderr, "********************************************************************************\n" );
#     endif
   }

#  else // #ifdef LAPLACIAN_4TH
   Temp = 1.0 - 4.0*Alpha2;

   if ( Temp >= 0.0 )
      Taylor3_Coeff = Damping * ( 1.0 - SQRT(Temp) ) / ( 16.0*Alpha2 );

   else
   {
      Taylor3_Coeff = real(1.0/6.0);

#     ifdef GAMER_DEBUG
      Aux_Message( stderr, "********************************************************************************\n" );
      Aux_Message( stderr, "ERROR : 1.0-4.0*Alpha2 = %12.5e < 0.0 (dt %12.5e, dh %12.5e, Eta %12.5e) !!\n",
                   Temp, dt, dh, Eta, __FUNCTION__ );
      Aux_Message( stderr, "        --> Program cannot find an optimized Taylor expansion coefficient.\n" );
      Aux_Message( stderr, "            Therefore, it is set to %12.5e by default.\n", Taylor3_Coeff );
      Aux_Message( stderr, "        --> Perhaps you should set DT__FLUID/DT__FLUID_INIT smaller ...\n" );
      Aux_Message( stderr, "        Rank <%d>, file <%s>, line <%d>, function <%s>\n",
                   MPI_Rank, __FILE__, __LINE__, __FUNCTION__ );
      Aux_Message( stderr, "********************************************************************************\n" );
#     endif
   }
#  endif // #ifdef LAPLACIAN_4TH ... else ...


   return Taylor3_Coeff;

} // FUNCTION : ELBDM_SetTaylor3Coeff



#endif // #if ( MODEL == ELBDM )
