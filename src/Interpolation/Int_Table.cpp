#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Int_Table
// Description :  Return the number of sides (6 or 26) and ghost zones for the queried interpolation scheme 
//
// Parameter   :  IntScheme : Interpolation scheme
//                            --> currently supported schemes include
//                                INT_MINMOD3D : MinMod-3D
//                                INT_MINMOD1D : MinMod-1D
//                                INT_VANLEER  : vanLeer
//                                INT_CQUAD    : conservative quadratic
//                                INT_QUAD     : quadratic
//                                INT_CQUAR    : conservative quartic
//                                INT_QUAR     : quartic
//
// Return      :  NSide, NGhost
//-------------------------------------------------------------------------------------------------------
void Int_Table( const IntScheme_t IntScheme, int &NSide, int &NGhost )
{

   switch ( IntScheme )
   {
      case INT_MINMOD3D :  NSide = 26;    NGhost = 1;    break;
      case INT_MINMOD1D :  NSide =  6;    NGhost = 1;    break;
      case INT_VANLEER  :  NSide = 26;    NGhost = 1;    break;
      case INT_CQUAD    :  NSide = 26;    NGhost = 1;    break;
      case INT_QUAD     :  NSide = 26;    NGhost = 1;    break;
      case INT_CQUAR    :  NSide = 26;    NGhost = 2;    break;
      case INT_QUAR     :  NSide = 26;    NGhost = 2;    break;

      default           :  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
                           exit(1);
   }

} // FUNCTION : Int_Table
