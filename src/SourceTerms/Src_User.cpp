#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Src_User( real fluid[], const double x, const double y, const double z, const double Time,
                      const int lv, double AuxArray[], const double dt );

// this function pointer may be overwritten by various test problem initializers
void (*Src_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                      const int lv, double AuxArray[], const double dt ) = Src_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_User
// Description :  User-defined source terms
//
// Note        :  1. Invoked by Src_AddSource() using the function pointer "Src_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "SRC_USER"
//
// Parameter   :  fluid    : Fluid array storing both the input and updated values
//                           --> Array size = NCOMP_TOTAL (so it includes both active and passive variables)
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//                dt       : Time interval to advance solution
//
// Return      :  fluid[]
//-------------------------------------------------------------------------------------------------------
void Src_User( real fluid[], const double x, const double y, const double z, const double Time,
               const int lv, double AuxArray[], const double dt )
{

// example
   /*
   const double CoolRate = 1.23; // set arbitrarily here
   double Ek, Eint;              // kinetic and internal energies

   Ek    = (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   Eint  = fluid[ENGY] - Ek;
   Eint -= CoolRate*dt;
   fluid[ENGY] = Ek + Eint;
   */

} // FUNCTION : Src_User
