#include "Copyright.h"
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

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "MG_MAX_ITER", Default_Max_Iter );
   }

   if ( NPre_Smooth < 0 )
   {
      NPre_Smooth = Default_NPre_Smooth;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "MG_NPRE_SMOOTH", Default_NPre_Smooth );
   }

   if ( NPost_Smooth < 0 )
   {
      NPost_Smooth = Default_NPost_Smooth;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "MG_NPOST_SMOOTH", Default_NPost_Smooth );
   }

   if ( Tolerated_Error < 0.0 )
   {
      Tolerated_Error = Default_Tolerated_Error;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "MG_TOLERATED_ERROR", Default_Tolerated_Error );
   }

} // FUNCTION : Init_Set_Default_MG_Parameter



#endif // #if ( defined GRAVITY  &&  POT_SCHEME == MG )
