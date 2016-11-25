#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize the particle position and velocity
//
// Note        :  1. Invoked by "Init_GAMER"
//                2. Assuming all particles are active initially
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass/Pos(X/Y/Z)/Vel(X/Y/Z)/Time
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );




   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
