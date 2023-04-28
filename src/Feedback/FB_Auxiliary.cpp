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



#endif // #ifdef FEEDBACK
