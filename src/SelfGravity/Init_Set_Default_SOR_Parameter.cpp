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


// check
#  if ( defined GRAVITY  &&  POT_GHOST_SIZE > 5 )
   if ( SOR_Omega < 0.0 )
      Aux_Error( ERROR_INFO, "function \"%s\" does not work for POT_GHOST_SIZE > 5 !!\n", __FUNCTION__ );
#  endif


   const double Default_Omega[5] = { 1.49, 1.57, 1.62, 1.65, 1.69 };    // for POT_GHOST_SIZE = [1,2,3,4,5]
#  ifdef FLOAT8
   const int    Default_MaxIter  = 110;
#  else
   const int    Default_MaxIter  = 60;
#  endif
   const int    Default_MinIter  = 10;

   if ( SOR_Omega < 0.0 )
   {
      SOR_Omega = Default_Omega[POT_GHOST_SIZE-1];

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
