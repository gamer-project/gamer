#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

#if ( MODEL == HYDRO )

#if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
extern real CPU_PositivePres( const real Pres_In, const real Dens, const real _Dens );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetPressure
// Description :  Evalute the fluid pressure
//
// Note        :  1. Currently only work with the adiabatic EOS 
//                2. Invoked by the functions "Hydro_GetTimeStep_Fluid", "Prepare_PatchData", "InterpolateGhostZone", 
//                   "Hydro_Aux_Check_Negative"
// 
// Parameter   :  Fluid          : Input fluid variables
//                Gamma_m1       : Gamma - 1.0
//                PositivePres   : Enforce the returned pressure to be positive by invoking the function "CPU_PositivePres"
//                                 --> useless if neither MIN_PRES_DENS nor MIN_PRES is defined
//
// Return      :  Pressure
//-------------------------------------------------------------------------------------------------------
real Hydro_GetPressure( const real Fluid[], const real Gamma_m1, const bool PositivePres )
{

   real _Dens, Pres;

  _Dens = (real)1.0 / Fluid[DENS];
   Pres = Gamma_m1*(  Fluid[ENGY] - (real)0.5*_Dens*( SQR(Fluid[MOMX]) + SQR(Fluid[MOMY]) + SQR(Fluid[MOMZ]) )  );
#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   if ( PositivePres )
   Pres = CPU_PositivePres( Pres, Fluid[DENS], _Dens );
#  endif

   return Pres;

} // FUNCTION : Hydro_GetPressure



#endif // #if ( MODEL == HYDRO )
