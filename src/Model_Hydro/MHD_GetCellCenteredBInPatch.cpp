#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBFieldInPatch
// Description :  Calculate the cell-centered magnetic field from the face-centered values of a given cell
//                in a given patch
//
// Note        :  1. Invoke MHD_GetCellCenteredBField() defined in CPU/CUFLU_Shared_FluUtility.cpp
//                2. Return all three components of the magnetic field
//
// Parameter   :  B_CC  : Cell-centered magnetic field to be returned
//                lv    : Target AMR level
//                PID   : Target patch index
//                i/j/k : Target array indices of the patch
//                MagSg : Sandglass of the magnetic field data
//
// Return      :  B_CC
//-------------------------------------------------------------------------------------------------------
void MHD_GetCellCenteredBFieldInPatch( real B_CC[], const int lv, const int PID, const int i, const int j, const int k,
                                       const int MagSg )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )   Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (total number of patches = %d) !!\n", PID, amr->num[lv] );

   if ( i < 0  ||  i >= PS1 )    Aux_Error( ERROR_INFO, "incorrect i = %d !!\n", i );
   if ( j < 0  ||  j >= PS1 )    Aux_Error( ERROR_INFO, "incorrect j = %d !!\n", j );
   if ( k < 0  ||  k >= PS1 )    Aux_Error( ERROR_INFO, "incorrect k = %d !!\n", k );

   if ( MagSg != 0  &&  MagSg != 1 )   Aux_Error( ERROR_INFO, "MagSg (%d) != 0/1 !!\n", MagSg );

   if ( amr->patch[MagSg][lv][PID]->magnetic == NULL )
      Aux_Error( ERROR_INFO, "magnetic == NULL (lv %d, PID %d, MagSg %d) !!\n", lv, PID, MagSg );
#  endif


// FC = face-centered
   const real *Bx_FC = amr->patch[MagSg][lv][PID]->magnetic[MAGX];
   const real *By_FC = amr->patch[MagSg][lv][PID]->magnetic[MAGY];
   const real *Bz_FC = amr->patch[MagSg][lv][PID]->magnetic[MAGZ];

   MHD_GetCellCenteredBField( B_CC, Bx_FC, By_FC, Bz_FC, PS1, PS1, PS1, i, j, k );

} // FUNCTION : MHD_GetCellCenteredBFieldInPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBEnergyInPatch
// Description :  Calculate the cell-centered magnetic energy (i.e., 0.5*B^2) from the face-centered values
//                of a given cell in a given patch
//
// Note        :  1. Invoke MHD_GetCellCenteredBEnergy() defined in CPU/CUFLU_Shared_FluUtility.cpp
//
// Parameter   :  lv    : Target AMR level
//                PID   : Target patch index
//                i/j/k : Target array indices of the patch
//                MagSg : Sandglass of the magnetic field data
//
// Return      :  0.5*B^2 at the center of the target cell
//-------------------------------------------------------------------------------------------------------
real MHD_GetCellCenteredBEnergyInPatch( const int lv, const int PID, const int i, const int j, const int k,
                                        const int MagSg )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )   Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (total number of patches = %d) !!\n", PID, amr->num[lv] );

   if ( i < 0  ||  i >= PS1 )    Aux_Error( ERROR_INFO, "incorrect i = %d !!\n", i );
   if ( j < 0  ||  j >= PS1 )    Aux_Error( ERROR_INFO, "incorrect j = %d !!\n", j );
   if ( k < 0  ||  k >= PS1 )    Aux_Error( ERROR_INFO, "incorrect k = %d !!\n", k );

   if ( MagSg != 0  &&  MagSg != 1 )   Aux_Error( ERROR_INFO, "MagSg (%d) != 0/1 !!\n", MagSg );

   if ( amr->patch[MagSg][lv][PID]->magnetic == NULL )
      Aux_Error( ERROR_INFO, "magnetic == NULL (lv %d, PID %d, MagSg %d) !!\n", lv, PID, MagSg );
#  endif


// FC = face-centered
   const real *Bx_FC = amr->patch[MagSg][lv][PID]->magnetic[MAGX];
   const real *By_FC = amr->patch[MagSg][lv][PID]->magnetic[MAGY];
   const real *Bz_FC = amr->patch[MagSg][lv][PID]->magnetic[MAGZ];

   return MHD_GetCellCenteredBEnergy( Bx_FC, By_FC, Bz_FC, PS1, PS1, PS1, i, j, k );

} // FUNCTION : MHD_GetCellCenteredBEnergyInPatch



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
