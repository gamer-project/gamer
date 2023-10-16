#include "GAMER.h"

#if ( defined GRAVITY  &&  POT_SCHEME == MG )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Set_Default_MG_Parameter
// Description :  Set the multigrid parameters to the default values
//
// Note        :  1. Work only when the corresponding input parameters are negative
//                2. Default values are determined empirically
//
// Return      :  MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, MG_TOLERATED_ERROR
//-------------------------------------------------------------------------------------------------------
void Init_Set_Default_MG_Parameter()
{

#  ifdef FLOAT8
   const int    Default_Max_Iter        = 20;
   const double Default_Tolerated_Error = 1.e-15;
#  else
   const int    Default_Max_Iter        = 10;
   const double Default_Tolerated_Error = 1.e-6;
#  endif
   const int    Default_NPre_Smooth     = 3;
   const int    Default_NPost_Smooth    = 3;


   if ( MG_MAX_ITER < 0 )
   {
      MG_MAX_ITER = Default_Max_Iter;

      PRINT_RESET_PARA( MG_MAX_ITER, FORMAT_INT, "" );
   }

   if ( MG_NPRE_SMOOTH < 0 )
   {
      MG_NPRE_SMOOTH = Default_NPre_Smooth;

      PRINT_RESET_PARA( MG_NPRE_SMOOTH, FORMAT_INT, "" );
   }

   if ( MG_NPOST_SMOOTH < 0 )
   {
      MG_NPOST_SMOOTH = Default_NPost_Smooth;

      PRINT_RESET_PARA( MG_NPOST_SMOOTH, FORMAT_INT, "" );
   }

   if ( MG_TOLERATED_ERROR < 0.0 )
   {
      MG_TOLERATED_ERROR = Default_Tolerated_Error;

      PRINT_RESET_PARA( MG_TOLERATED_ERROR, FORMAT_REAL, "" );
   }

} // FUNCTION : Init_Set_Default_MG_Parameter



#endif // #if ( defined GRAVITY  &&  POT_SCHEME == MG )
