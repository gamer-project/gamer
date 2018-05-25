#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBField
// Description :  Calculate the cell-centered magnetic field from the face-centered values of a given cell
//                in a given patch
//
// Note        :  1. Use the central average operator
//                2. Return all three components of the magnetic field
//                3. Always output the value stored in amr->MagSg[lv]
//
// Parameter   :  B_CC  : Cell-centered magnetic field to be returned
//                lv    : Target AMR level
//                PID   : Target patch index
//                i/j/k : Target array indices of the patch
//
// Return      :  B_CC
//-------------------------------------------------------------------------------------------------------
void MHD_GetCellCenteredBField( real B_CC[], const int lv, const int PID, const int i, const int j, const int k )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (total number of patches = %d) !!\n", PID, amr->num[lv] );

   if ( i < 0  ||  i >= PS1 )    Aux_Error( ERROR_INFO, "incorrect i = %d !!\n", i );
   if ( j < 0  ||  j >= PS1 )    Aux_Error( ERROR_INFO, "incorrect j = %d !!\n", j );
   if ( k < 0  ||  k >= PS1 )    Aux_Error( ERROR_INFO, "incorrect k = %d !!\n", k );

   if ( amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic == NULL )
         Aux_Error( ERROR_INFO, "magnetic == NULL (lv %d, PID %d, MagSg %d) !!\n", lv, PID, amr->MagSg[lv] );
#  endif


// FC = face-centered
   const real (*B_FC)[PS1_P1*PS1*PS1] = amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic;

   B_CC[0] = (real)0.5*(  B_FC[0][ IDX321_BX(i,j,k) ] + B_FC[0][ IDX321_BX(i+1,j,  k  ) ]  );
   B_CC[1] = (real)0.5*(  B_FC[1][ IDX321_BY(i,j,k) ] + B_FC[1][ IDX321_BY(i,  j+1,k  ) ]  );
   B_CC[2] = (real)0.5*(  B_FC[2][ IDX321_BZ(i,j,k) ] + B_FC[2][ IDX321_BZ(i,  j,  k+1) ]  );

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
//
// Return      :  0.5*B^2 at the center of the target cell
//-------------------------------------------------------------------------------------------------------
real MHD_GetCellCenteredBEnergy( const int lv, const int PID, const int i, const int j, const int k )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (total number of patches = %d) !!\n", PID, amr->num[lv] );

   if ( i < 0  ||  i >= PS1 )    Aux_Error( ERROR_INFO, "incorrect i = %d !!\n", i );
   if ( j < 0  ||  j >= PS1 )    Aux_Error( ERROR_INFO, "incorrect j = %d !!\n", j );
   if ( k < 0  ||  k >= PS1 )    Aux_Error( ERROR_INFO, "incorrect k = %d !!\n", k );
#  endif


// CC = cell-centered
   real B_CC[3], BEngy;

   MHD_GetCellCenteredBField( B_CC, lv, PID, i, j, k );

   BEngy = (real)0.5*( SQR(B_CC[0]) + SQR(B_CC[1]) + SQR(B_CC[2]) );

   return BEngy;

} // FUNCTION : MHD_GetCellCenteredBEnergy



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
