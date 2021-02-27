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



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredDivBInPatch
// Description :  Calculate the cell-centered divergence(magnetic field) from the face-centered values
//                of a given cell in a given patch
//
// Note        :  1. The MHD algorithms in GAMER should guarantee the divergence-free constraint to the
//                   machine precision
//                   --> This function is thus mainly for debugging purposes
//                2. To be more specific, this function computes the dimensionless quantity |div(B)*dh/(3*<|B|>)|
//                   --> Except when B=0 on all 6 faces, for which it simply returns 0
//
// Parameter   :  lv    : Target AMR level
//                PID   : Target patch index
//                i/j/k : Target array indices of the patch
//                MagSg : Sandglass of the magnetic field data
//
// Return      :  |div(B)*dh/(3*<|B|>)| at the center of the target cell
//-------------------------------------------------------------------------------------------------------
real MHD_GetCellCenteredDivBInPatch( const int lv, const int PID, const int i, const int j, const int k,
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


   const real (*B)[PS1P1*PS1*PS1] = amr->patch[MagSg][lv][PID]->magnetic;

   const int idx_BxL = IDX321_BX( i, j, k, PS1, PS1 );
   const int idx_ByL = IDX321_BY( i, j, k, PS1, PS1 );
   const int idx_BzL = IDX321_BZ( i, j, k, PS1, PS1 );

   real BxL, BxR, ByL, ByR, BzL, BzR, DivB, AmpB;

   BxL = B[MAGX][ idx_BxL            ];
   BxR = B[MAGX][ idx_BxL + 1        ];
   ByL = B[MAGY][ idx_ByL            ];
   ByR = B[MAGY][ idx_ByL + PS1      ];
   BzL = B[MAGZ][ idx_BzL            ];
   BzR = B[MAGZ][ idx_BzL + SQR(PS1) ];

   DivB = ( BxR - BxL ) + ( ByR - ByL ) + ( BzR - BzL );
   AmpB = FABS(BxR) + FABS(BxL) + FABS(ByR) + FABS(ByL) + FABS(BzR) + FABS(BzL);

// return 0 if B is zero on all 6 faces
   if ( AmpB != (real)0.0 )   DivB = FABS( DivB/AmpB );

   return DivB;

} // FUNCTION : MHD_GetCellCenteredDivBInPatch



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
