#include "GAMER.h"

#if ( defined GRAVITY  &&  POT_SCHEME == SOR )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Set_Default_SOR_Parameter
// Description :  Set the SOR parameters to the default values
//
// Note        :  1. Work only when the corresponding input parameters are negative
//                2. The default values are determined empirically from the cosmological simulations
//
// Return      :  SOR_OMEGA, SOR_MAX_ITER, SOR_MIN_ITER
//-------------------------------------------------------------------------------------------------------
void Init_Set_Default_SOR_Parameter()
{

// reference to the optimum relaxation parameter: Yang & Gobbert, Appl. Math. Lett. 22, 325 (2009)
// --> Eq. [1.3] and the last paragraph in Sec. 3
// --> further confirmed with the test problem "Hydro/Gravity" using different PATCH_SIZE and POT_GHOST_SIZE
   const int    NCell           = PS1 + 2*POT_GHOST_SIZE;
   const double Default_Omega   = 2.0 / ( 1.0 + sin(M_PI/NCell) );
#  ifdef FLOAT8
   const int    Default_MaxIter = int( 13 + 10*NCell );  // add a factor of 2 buffer to the empirically determined maximum iteration
#  else
   const int    Default_MaxIter = int( 40 + 3*NCell );   // add a factor of 2 buffer to the empirically determined maximum iteration
#  endif
   const int    Default_MinIter = Default_MaxIter / 5;   // 20% of the maximum iteration (determined empirically)


   if ( SOR_OMEGA < 0.0 )
   {
      SOR_OMEGA = Default_Omega;

      PRINT_RESET_PARA( SOR_OMEGA, FORMAT_REAL, "" );
   }

   if ( SOR_MAX_ITER < 0 )
   {
      SOR_MAX_ITER = Default_MaxIter;

      PRINT_RESET_PARA( SOR_MAX_ITER, FORMAT_INT, "" );
   }

   if ( SOR_MIN_ITER < 0 )
   {
      SOR_MIN_ITER = Default_MinIter;

      PRINT_RESET_PARA( SOR_MIN_ITER, FORMAT_INT, "" );
   }

} // FUNCTION : Init_Set_Default_SOR_Parameter



#endif // #if ( defined GRAVITY  &&  POT_SCHEME == SOR )
