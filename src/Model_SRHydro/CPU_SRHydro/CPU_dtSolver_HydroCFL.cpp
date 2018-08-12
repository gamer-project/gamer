#include "GAMER.h"
#include "../../../include/CPU_prototypes.h"


#if ( !defined GPU  &&  MODEL == SR_HYDRO )

//-----------------------------------------------------------------------------------------
// Function    :  CPU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) required for the hydro solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
//
// Parameter   :  dt_Array  : Array to store the minimum dt in each target patch
//                Flu_Array : Array storing the prepared fluid data of each target patch
//                NPG       : Number of target patch groups
//                dh        : Grid size
//                Safety    : dt safety factor
//                Gamma     : Ratio of specific heats
//                MinPres   : Minimum allowed pressure
//
// Return      :  dt_Array
//-----------------------------------------------------------------------------------------
void CPU_dtSolver_HydroCFL( real dt_Array[], const real Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                            const int NPG, const real dh, const real Safety, const real Gamma, const real MinPres )
{
   const int  NPatch           = 8*NPG;


   for (int p=0; p<NPatch; p++)
   {
      dt_Array[p] = Safety*dh;

   } // for (int p=0; p<NPatch; p++)

} // FUNCTION : CPU_dtSolver_HydroCFL



#endif // #if ( !defined GPU  &&  MODEL == SR_HYDRO )
