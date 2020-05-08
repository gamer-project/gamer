#include "CUFLU.h"

#if ( MODEL == HYDRO )



// internal function prototypes
// --> only necessary for GPU since they are included in Prototype.h for the CPU codes
//#ifdef __CUDACC__
//GPU_DEVICE
//static real Hydro_CheckMinPres( const real InPres, const real MinPres );
//#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres
// Description :  Convert gas mass density and internal energy density to pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//
// Parameter   :  Dens : Gas mass density
//                Eint : Gas internal energy density
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real EoS_DensEint2Pres( const real Dens, const real Eint )
{


} // FUNCTION : EoS_DensEint2Pres



#endif // #if ( MODEL == HYDRO )
