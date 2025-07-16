#include "GAMER.h"


#ifdef SUPPORT_HYPRE
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_Free
// Description :  Destory the Hypre arrays
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Hypre_Free()
{

   HYPRE_CHECK_FUNC(   HYPRE_SStructGridDestroy( Hypre_grid )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructGraphDestroy( Hypre_graph )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructStencilDestroy( Hypre_stencil )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixDestroy( Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorDestroy( Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorDestroy( Hypre_b )   );

} // FUNCITON : Hypre_Free
#endif // #ifdef SUPPORT_HYPRE
