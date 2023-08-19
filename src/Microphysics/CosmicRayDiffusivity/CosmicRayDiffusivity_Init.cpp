#include "GAMER.h"
#include <typeinfo>

#ifdef COSMICRAY

//-------------------------------------------------------------------------------------------------------
// Function    :  CosmicRay_Init
// Description :  Initialize cosmic ray
// 
// Note        : 
// 
// Parameter   :
// 
// Return      :
//-------------------------------------------------------------------------------------------------------
void CosmicRay_Init()
{
  if (MPI_Rank == 0)    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

  // iso

  // not iso

} // FUNCTION : CosmicRay_Init

#endif // #ifdef COSMICRAY
