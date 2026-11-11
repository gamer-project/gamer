#include "GAMER.h"

#ifdef FEEDBACK



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Aux_CellPatchRelPos
// Description :  Return the relative position between a given cell index and its host patch group
//
// Note        :  1. Invoked by various feedback routines to check whether particles are too close to
//                   coarse-fine boundaries
//                2. Assume FB_GHOST_SIZE ghost zones surrounding a patch group
//                3. Use the sibling index 0-25 to specify the relative position
//
// Parameter   :  ijk : Cell indices
//
// Return      :  Sibling index 0-25 (-1 means the input cell is within the patch group)
//-------------------------------------------------------------------------------------------------------
int FB_Aux_CellPatchRelPos( const int ijk[] )
{

   const int SibID[3][3][3] = {  { {18, 10, 19}, {14,  4, 16}, {20, 11, 21} },
                                 { { 6,  2,  7}, { 0, -1,  1}, { 8,  3,  9} },
                                 { {22, 12, 23}, {15,  5, 17}, {24, 13, 25} }  };
   int RelPos[3];
   for (int d=0; d<3; d++)    RelPos[d] = ( ijk[d] < FB_GHOST_SIZE     ) ? 0 :
                                          ( ijk[d] < FB_GHOST_SIZE+PS2 ) ? 1 : 2;

   return SibID[ RelPos[2] ][ RelPos[1] ][ RelPos[0] ];

} // FUNCTION : FB_Aux_CellPatchRelPos



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Aux_PatchLocalID
// Description :  For a given cell index, return the patch local index in its host patch group + nearby patches
//
// Note        :  1. Invoked by various feedback routines to check whether particles are too close to
//                   non-leaf patches
//                2. Assume FB_GHOST_SIZE ghost zones surrounding a patch group
//                3. Use the patch local index 0-63 to specify the patch
//
// Parameter   :  ijk : Cell indices
//
// Return      :  Patch local index 0-63 (0-7 means the input cell is within the patch group)
//-------------------------------------------------------------------------------------------------------
int FB_Aux_PatchLocalID( const int ijk[] )
{

   const int PatchLocalID[4][4][4] = {  { {56, 40, 41, 57}, {48, 24, 26, 52}, {49, 25, 27, 53}, {58, 42, 43, 59} },
                                        { {32, 16, 17, 34}, { 8,  0,  1, 12}, { 9,  2,  4, 13}, {36, 20, 21, 38} },
                                        { {33, 18, 19, 35}, {10,  3,  6, 14}, {11,  5,  7, 15}, {37, 22, 23, 39} },
                                        { {60, 44, 45, 61}, {50, 28, 29, 54}, {51, 30, 31, 55}, {62, 46, 47, 63} }  };
   int RelPos[3];
   for (int d=0; d<3; d++)    RelPos[d] = ( ijk[d] < FB_GHOST_SIZE     ) ? 0 :
                                          ( ijk[d] < FB_GHOST_SIZE+PS1 ) ? 1 :
                                          ( ijk[d] < FB_GHOST_SIZE+PS2 ) ? 2 : 3;

   return PatchLocalID[ RelPos[2] ][ RelPos[1] ][ RelPos[0] ];

} // FUNCTION : FB_Aux_PatchLocalID



#endif // #ifdef FEEDBACK
