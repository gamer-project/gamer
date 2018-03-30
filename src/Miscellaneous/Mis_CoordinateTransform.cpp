#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Scale2PhySize
// Description :  Scale --> corresponding physical size
//
// Parameter   :  Scale :  Input scale
//
// Return      :  Corresponding physical size
//-------------------------------------------------------------------------------------------------------
double Mis_Scale2PhySize( const int Scale )
{

   return Scale*amr->dh[TOP_LEVEL];

} // FUNCTION : Mis_Scale2PhySize



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Cell2PhySize
// Description :  Number of cells at the target level --> corresponding physical size
//
// Parameter   :  NCell :  Input number of cells
//                lv    :  Target level
//
// Return      :  Corresponding physical size
//-------------------------------------------------------------------------------------------------------
double Mis_Cell2PhySize( const int NCell, const int lv )
{

   return NCell*amr->dh[lv];

} // FUNCTION : Mis_Cell2PhySize



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Scale2Cell
// Description :  Scale --> corresponding number of cells at the target level
//
// Parameter   :  Scale :  Input scale
//                lv    :  Target level
//
// Return      :  Corresponding number of cells at lv
//-------------------------------------------------------------------------------------------------------
int Mis_Scale2Cell( const int Scale, const int lv )
{

#  ifdef GAMER_DEBUG
   if ( Scale % amr->scale[lv] != 0 )
      Aux_Message( stderr, "WARNING : input scale (%d) is not a multiple of the scale unit (lv %d --> %d) !!\n",
                   Scale, lv, amr->scale[lv] );
#  endif

   return ( Scale >> (TOP_LEVEL-lv) );

} // FUNCTION : Mis_Scale2Cell



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Cell2Scale
// Description :  Number of cells at the target level --> scale
//
// Parameter   :  NCell :  Input number of cells
//                lv    :  Target level
//
// Return      :  Corresponding scale
//-------------------------------------------------------------------------------------------------------
int Mis_Cell2Scale( const int NCell, const int lv )
{

   return NCell*amr->scale[lv];

} // FUNCTION : Mis_Cell2Scale
