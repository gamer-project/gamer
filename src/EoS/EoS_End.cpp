#include "GAMER.h"

#if ( MODEL == HYDRO )



// prototypes of built-in EoS
#if   ( EOS == EOS_GAMMA )
// nothing to do
#elif ( EOS == EOS_ISOTHERMAL )
// nothing to do
#elif ( EOS == EOS_NUCLEAR )
# error : ERROR : EOS_NUCLEAR is NOT supported yet !!
#endif // # EOS

// this function pointer can be set by a test problem initializer for non-built-in EoS
void (*EoS_End_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_End
// Description :  Free the resources used by the EoS routines
//
// Note        :  1. Invoked by End_GAMER()
//                2. For a non-built-in EoS, "EoS_End_Ptr" can be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void EoS_End()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// set the function pointer for the built-in EoS
#  if   ( EOS == EOS_GAMMA )
// nothing to do
#  elif ( EOS == EOS_ISOTHERMAL )
// nothing to do
#  elif ( EOS == EOS_NUCLEAR )
#  error : ERROR : EOS_NUCLEAR is NOT supported yet !!
#  endif // # EOS


   if ( EoS_End_Ptr != NULL )    EoS_End_Ptr();


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : EoS_End



#endif // #if ( MODEL == HYDRO )
