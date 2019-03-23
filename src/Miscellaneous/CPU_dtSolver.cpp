#include "GAMER.h"

#ifndef GPU


#ifdef GRAVITY
#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];
#endif


#if   ( MODEL == HYDRO )
void CPU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                            const int NPG, const real dh, const real Safety, const real Gamma, const real MinPres );
#ifdef GRAVITY
void CPU_dtSolver_HydroGravity( real g_dt_Array[],
                                const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                const double g_Corner_Array[][3],
                                const int NPatchGroup, const real dh, const real Safety, const bool P5_Gradient,
                                const OptGravityType_t GravityType, const double c_ExtAcc_AuxArray[],
                                const double ExtAcc_Time );
#endif

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#elif ( MODEL == SR_HYDRO )
void CPU_dtSolver_SRHydroCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                              const int NPG, const real dh, const real Safety, const real Gamma, const real MinPres );

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
// Parameter   :  TSolver      : Target dt solver
//                               --> DT_FLU_SOLVER : dt solver for fluid
//                                   DT_GRA_SOLVER : dt solver for gravity
//                dt_Array     : Array to store the minimum dt in each target patch
//                Flu_Array    : Array storing the prepared fluid data of each target patch
//                Pot_Array    : Array storing the prepared potential data of each target patch
//                Corner_Array : Array storing the physical corner coordinates of each patch
//                NPatchGroup  : Number of patch groups evaluated simultaneously by GPU
//                dh           : Cell size
//                Safety       : dt safety factor
//                Gamma        : Ratio of specific heats
//                MinPres      : Minimum allowed pressure
//                P5_Gradient  : Use 5-points stencil to evaluate the potential gradient
//                GravityType  : Types of gravity --> self-gravity, external gravity, both
//                ExtPot       : Add the external potential for ELBDM
//                TargetTime   : Target physical time
//
// Return      :  dt_Array
//-------------------------------------------------------------------------------------------------------
void CPU_dtSolver( const Solver_t TSolver, real dt_Array[], const real Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                   const real Pot_Array[][ CUBE(GRA_NXT) ], const double Corner_Array[][3],
                   const int NPatchGroup, const real dh, const real Safety, const real Gamma, const real MinPres,
                   const bool P5_Gradient, const OptGravityType_t GravityType, const bool ExtPot, const double TargetTime )
{

   switch ( TSolver )
   {
#     if   ( MODEL == HYDRO )
      case DT_FLU_SOLVER:
         CPU_dtSolver_HydroCFL( dt_Array, Flu_Array, NPatchGroup, dh, Safety, Gamma, MinPres );
      break;

#     ifdef GRAVITY
      case DT_GRA_SOLVER:
         CPU_dtSolver_HydroGravity( dt_Array, Pot_Array, Corner_Array, NPatchGroup, dh, Safety, P5_Gradient,
                                    GravityType, ExtAcc_AuxArray, TargetTime );
      break;
#     endif


#     elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#     elif   ( MODEL == SR_HYDRO )
      case DT_FLU_SOLVER:
         CPU_dtSolver_SRHydroCFL( dt_Array, Flu_Array, NPatchGroup, dh, Safety, Gamma, MinPres );
      break;

#     ifdef GRAVITY
      case DT_GRA_SOLVER:
         CPU_dtSolver_HydroGravity( dt_Array, Pot_Array, Corner_Array, NPatchGroup, dh, Safety, P5_Gradient,
                                    GravityType, ExtAcc_AuxArray, TargetTime );
      break;
#     endif



#     elif ( MODEL == ELBDM )


#     else
#        error : ERROR : unsupported MODEL !!
#     endif // MODEL

      default :
         Aux_Error( ERROR_INFO, "unsupported \"TSolver\" (%d) !!\n", TSolver );
   } // switch ( TSolver )

} // FUNCTION : CPU_dtSolver



#endif // #ifndef GPU
