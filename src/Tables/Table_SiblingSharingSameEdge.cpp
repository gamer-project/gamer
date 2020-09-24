#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_SiblingSharingSameEdge
// Description :  Return the three sibling direction indices sharing the same patch edge as the input one
//
// Note        :  1. Results are stored in SibID[]
//                2. Direction of edge indices (EdgeID) are defined in the same way as the sibling indices
//                3. If SibSibID != NULL, this function also returns the sibling directions **from the
//                   sibling patches associated with SibID[] to the target patch edge**
//
// Parameter   :  EdgeID   : Sibling index corresponding to the target patch edge
//                SibID    : Output array (see the notes above)
//                SibSibID : Output array (see the notes above)
//
// Return      :  SibID[], SibSibID[]
//-------------------------------------------------------------------------------------------------------
void TABLE_SiblingSharingSameEdge( const int EdgeID, int SibID[], int SibSibID[] )
{

#  ifdef GAMER_DEBUG
   if ( EdgeID < 6  ||  EdgeID >= 18 )
      Aux_Error( ERROR_INFO, "EdgeID (%d) not in the accepted range [6 ... 17] !!\n", EdgeID );
#  endif

   switch ( EdgeID )
   {
      case 6 :
            SibID[0] = 0;      SibID[1] = 2;      SibID[2] = 6;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 7;   SibSibID[1] = 8;   SibSibID[2] = 9; }
         break;

      case 7 :
            SibID[0] = 1;      SibID[1] = 2;      SibID[2] = 7;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 6;   SibSibID[1] = 9;   SibSibID[2] = 8; }
         break;

      case 8 :
            SibID[0] = 0;      SibID[1] = 3;      SibID[2] = 8;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 9;   SibSibID[1] = 6;   SibSibID[2] = 7; }
         break;

      case 9 :
            SibID[0] = 1;      SibID[1] = 3;      SibID[2] = 9;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 8;   SibSibID[1] = 7;   SibSibID[2] = 6; }
         break;

      case 10 :
            SibID[0] = 2;      SibID[1] = 4;      SibID[2] = 10;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 11;  SibSibID[1] = 12;  SibSibID[2] = 13; }
         break;

      case 11 :
            SibID[0] = 3;      SibID[1] = 4;      SibID[2] = 11;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 10;  SibSibID[1] = 13;  SibSibID[2] = 12; }
         break;

      case 12 :
            SibID[0] = 2;      SibID[1] = 5;      SibID[2] = 12;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 13;  SibSibID[1] = 10;  SibSibID[2] = 11; }
         break;

      case 13 :
            SibID[0] = 3;      SibID[1] = 5;      SibID[2] = 13;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 12;  SibSibID[1] = 11;  SibSibID[2] = 10; }
         break;

      case 14 :
            SibID[0] = 0;      SibID[1] = 4;      SibID[2] = 14;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 16;  SibSibID[1] = 15;  SibSibID[2] = 17; }
         break;

      case 15 :
            SibID[0] = 0;      SibID[1] = 5;      SibID[2] = 15;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 17;  SibSibID[1] = 14;  SibSibID[2] = 16; }
         break;

      case 16 :
            SibID[0] = 1;      SibID[1] = 4;      SibID[2] = 16;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 14;  SibSibID[1] = 17;  SibSibID[2] = 15; }
         break;

      case 17 :
            SibID[0] = 1;      SibID[1] = 5;      SibID[2] = 17;
         if ( SibSibID != NULL ) {
         SibSibID[0] = 15;  SibSibID[1] = 16;  SibSibID[2] = 14; }
         break;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "EdgeID", EdgeID );
         break;

   } // switch ( EdgeID )

} // FUNCTION : TABLE_SiblingSharingSameEdge
