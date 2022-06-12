#include "GAMER.h"

#if ( defined GRAVITY  &&  POT_SCHEME == SOR )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Set_Default_SOR_Parameter
// Description :  Set the SOR parameters to the default values
//
// Note        :  1. Work only when the corresponding input parameters are negative
//                2. The default values are determined empirically from the cosmological simulations
//
// Parameter   :  SOR_Omega    : Over-relaxation parameter for SOR
//                SOR_Max_Iter : Maximum number of iterations for SOR
//                SOR_Min_Iter : Minimum number of iterations for SOR
//-------------------------------------------------------------------------------------------------------
void Init_Set_Default_SOR_Parameter( double &SOR_Omega, int &SOR_Max_Iter, int &SOR_Min_Iter )
{

// helper macro for printing warning messages
#  define FORMAT_INT    %- 21d
#  define FORMAT_FLT    %- 21.14e
#  define PRINT_WARNING( name, value, format, reason )                                                         \
   {                                                                                                           \
      if ( MPI_Rank == 0 )                                                                                     \
         Aux_Message( stderr, "WARNING : parameter [%-28s] is reset to [" EXPAND_AND_QUOTE(format) "] %s\n",   \
                      #name, value, reason );                                                                  \
   }


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

   if ( SOR_Omega < 0.0 )
   {
      SOR_Omega = Default_Omega;

      PRINT_WARNING( SOR_OMEGA, SOR_Omega, FORMAT_FLT, "" );
   }

   if ( SOR_Max_Iter < 0 )
   {
      SOR_Max_Iter = Default_MaxIter;

      PRINT_WARNING( SOR_MAX_ITER, SOR_Max_Iter, FORMAT_INT, "" );
   }

   if ( SOR_Min_Iter < 0 )
   {
      SOR_Min_Iter = Default_MinIter;

      PRINT_WARNING( SOR_MIN_ITER, SOR_Min_Iter, FORMAT_INT, "" );
   }


// remove symbolic constants and macros only used in this structure
#  undef FORMAT_INT
#  undef FORMAT_FLT
#  undef QUOTE

} // FUNCTION : Init_Set_Default_SOR_Parameter



#endif // #if ( defined GRAVITY  &&  POT_SCHEME == SOR )
