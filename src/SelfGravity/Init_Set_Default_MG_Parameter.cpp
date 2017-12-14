#include "GAMER.h"

#if ( defined GRAVITY  &&  POT_SCHEME == MG )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Set_Default_MG_Parameter
// Description :  Set the multigrid parameters to the default values
//
// Note        :  1. Work only when the corresponding input parameters are negative
//                2. Default values are determined empirically
//
// Parameter   :  Max_Iter        : Maximum number of iterations
//                NPre_Smooth     : Number of pre-smoothing steps
//                NPos_tSmooth    : Number of post-smoothing steps
//                Tolerated_Error : Maximum tolerated error
//-------------------------------------------------------------------------------------------------------
void Init_Set_Default_MG_Parameter( int &Max_Iter, int &NPre_Smooth, int &NPost_Smooth, double &Tolerated_Error )
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


#  ifdef FLOAT8
   const int    Default_Max_Iter        = 20;
   const double Default_Tolerated_Error = 1.e-15;
#  else
   const int    Default_Max_Iter        = 10;
   const double Default_Tolerated_Error = 1.e-6;
#  endif
   const int    Default_NPre_Smooth     = 3;
   const int    Default_NPost_Smooth    = 3;

   if ( Max_Iter < 0 )
   {
      Max_Iter = Default_Max_Iter;

      PRINT_WARNING( MG_MAX_ITER, Max_Iter, FORMAT_INT, "" );
   }

   if ( NPre_Smooth < 0 )
   {
      NPre_Smooth = Default_NPre_Smooth;

      PRINT_WARNING( MG_NPRE_SMOOTH, NPre_Smooth, FORMAT_INT, "" );
   }

   if ( NPost_Smooth < 0 )
   {
      NPost_Smooth = Default_NPost_Smooth;

      PRINT_WARNING( MG_NPOST_SMOOTH, NPost_Smooth, FORMAT_INT, "" );
   }

   if ( Tolerated_Error < 0.0 )
   {
      Tolerated_Error = Default_Tolerated_Error;

      PRINT_WARNING( MG_TOLERATED_ERROR, Tolerated_Error, FORMAT_FLT, "" );
   }


// remove symbolic constants and macros only used in this structure
#  undef FORMAT_INT
#  undef FORMAT_FLT
#  undef QUOTE

} // FUNCTION : Init_Set_Default_MG_Parameter



#endif // #if ( defined GRAVITY  &&  POT_SCHEME == MG )
