#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBField
// Description :  Calculate the cell-centered magnetic field from the face-centered values of a given cell
//                in a given patch
//
// Note        :  1. Use the central average operator
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
void MHD_GetCellCenteredBField( real B_CC[], const int lv, const int PID, const int i, const int j, const int k,
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
   const real (*B_FC)[PS1_P1*PS1*PS1] = amr->patch[MagSg][lv][PID]->magnetic;
   const int idx_Bx = IDX321_BX( i, j, k, PS1 );
   const int idx_By = IDX321_BY( i, j, k, PS1 );
   const int idx_Bz = IDX321_BZ( i, j, k, PS1 );

   B_CC[0] = (real)0.5*( B_FC[0][idx_Bx] + B_FC[0][ idx_Bx + 1        ] );
   B_CC[1] = (real)0.5*( B_FC[1][idx_By] + B_FC[1][ idx_By + PS1      ] );
   B_CC[2] = (real)0.5*( B_FC[2][idx_Bz] + B_FC[2][ idx_Bz + SQR(PS1) ] );

} // FUNCTION : MHD_GetCellCenteredBField



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBEnergy
// Description :  Calculate the cell-centered magnetic energy (i.e., 0.5*B^2) from the face-centered values
//                of a given cell in a given patch
//
// Note        :  1. Invoke MHD_GetCellCenteredBField()
//
// Parameter   :  lv    : Target AMR level
//                PID   : Target patch index
//                i/j/k : Target array indices of the patch
//                MagSg : Sandglass of the magnetic field data
//
// Return      :  0.5*B^2 at the center of the target cell
//-------------------------------------------------------------------------------------------------------
real MHD_GetCellCenteredBEnergy( const int lv, const int PID, const int i, const int j, const int k,
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
#  endif


// CC = cell-centered
   real B_CC[3], BEngy;

   MHD_GetCellCenteredBField( B_CC, lv, PID, i, j, k, MagSg );

   BEngy = (real)0.5*( SQR(B_CC[0]) + SQR(B_CC[1]) + SQR(B_CC[2]) );

   return BEngy;

} // FUNCTION : MHD_GetCellCenteredBEnergy



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
