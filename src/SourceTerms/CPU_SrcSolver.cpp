#include "GAMER.h"

#ifndef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_SrcSolver
// Description :  CPU solver for adding source terms
//
// Note        :  1. Input FLU_NIN_S fluid variables and NCOMP_MAG B field components
//                   --> FLU_NIN_S is fixed to NCOMP_TOTAL for now
//                   --> Should remove unused fields in the future
//                2. Use patches instead of patch groups as the basic unit
//                3. No ghost zones
//                   --> Should support ghost zones in the future
//
// Parameter   :  h_Flu_Array_In    : Host array storing the input fluid variables
//                h_Flu_Array_Out   : Host array to store the output fluid variables
//                h_Mag_Array_In    : Host array storing the input B field (for MHD only)
//                NPatchGroup       : Number of patch groups to be evaluated
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//-------------------------------------------------------------------------------------------------------
void CPU_SrcSolver( const real h_Flu_Array_In [][FLU_NIN_S ][ CUBE(PS1)      ],
                          real h_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)      ],
                    const real h_Mag_Array_In [][NCOMP_MAG ][ PS1P1*SQR(PS1) ],
                    const int NPatchGroup, const real dt, const real dh,
                    const double TimeNew, const double TimeOld,
                    const real MinDens, const real MinPres, const real MinEint )
{

// check
#  ifdef GAMER_DEBUG
   if ( h_Flu_Array_In  == NULL )   Aux_Error( ERROR_INFO, "h_Flu_Array_In = NULL !!\n" );
   if ( h_Flu_Array_Out == NULL )   Aux_Error( ERROR_INFO, "h_Flu_Array_Out = NULL !!\n" );
   if ( h_Mag_Array_In  == NULL )   Aux_Error( ERROR_INFO, "h_Mag_Array_In = NULL !!\n" );
#  endif

} // FUNCTION : CPU_SrcSolver



#endif // #ifndef GPU
