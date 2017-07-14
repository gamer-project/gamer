#include "Copyright.h"
#include "GAMER.h"

#ifndef GPU

#if   ( MODEL == HYDRO )
void CPU_dtSolver_HydroCFL( real dt_Array[], const real Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                            const int NPG, const real dh, const real Safety, const real Gamma, const real MinPres );

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_dtSolver
// Description :  Invoke various CPU dt solvers
//
// Note        :  1. Invoked by InvokeSolver()
//                2. Corresponding preparation and closing steps are defined in dt_Prepare_XXX() and dt_Close()
//
// Parameter   :  TSolver   : Target dt solver
//                dt_Array  : Array to store the minimum dt in each target patch
//                Flu_Array : Array storing the prepared fluid data of each target patch
//                Pot_Array : Array storing the prepared potential data of each target patch
//                NPG       : Number of target patch groups
//                dh        : Grid size
//                Safety    : dt safety factor
//                Gamma     : Ratio of specific heats
//                MinPres   : Minimum allowed pressure
//                NewtonG   : Gravitational constant
//
// Return      :  dt_Array
//-------------------------------------------------------------------------------------------------------
void CPU_dtSolver( const Solver_t TSolver, real dt_Array[],
                   const real Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ], const real Pot_Array[][ CUBE(GRA_NXT) ],
                   const int NPG, const real dh, const real Safety, const real Gamma, const real MinPres,
                   const real NewtonG )
{

   switch ( TSolver )
   {
#     if   ( MODEL == HYDRO )
      case DT_SOLVER_HYDRO_CFL:
         CPU_dtSolver_HydroCFL( dt_Array, Flu_Array, NPG, dh, Safety, Gamma, MinPres );
      break;

#     ifdef GRAVITY
      case DT_SOLVER_HYDRO_GRAVITY:
         CPU_dtSolver_HydroGravity( dt_Array, Pot_Array, NPG, dh, Safety, NewtonG );
      break;
#     endif


#     elif ( MODEL == MHD )
#        warning : WAIT MHD !!!


#     elif ( MODEL == ELBDM )


#     else
#        error : ERROR : unsupported MODEL !!
#     endif // MODEL

      default :
         Aux_Error( ERROR_INFO, "unsupported \"TSolver\" (%d) !!\n", TSolver );
   } // switch ( TSolver )

} // FUNCTION : CPU_dtSolver



#endif // #ifndef GPU
