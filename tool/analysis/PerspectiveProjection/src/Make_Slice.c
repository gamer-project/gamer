#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Slice
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real Slice( int Idx1, int Idx2, int numCellXYZ[], real ***TargetQuantity, char CuttingPlane[] )
{
   if ( atoi(CuttingPlane) == atoi("x") )
   {
      int CenterIdx = 0.5 * numCellXYZ[0];
      return TargetQuantity[CenterIdx][Idx1][Idx2];
   }
   else if ( atoi(CuttingPlane) == atoi("y") )
   {
      int CenterIdy = 0.5 * numCellXYZ[1];
      return TargetQuantity[Idx1][CenterIdy][Idx2];
   }
   else if ( atoi(CuttingPlane) == atoi("z") )
   {
      int CenterIdz = 0.5 * numCellXYZ[2];
      return TargetQuantity[Idx1][Idx2][CenterIdz];
   }
   else
   {
      ERROR_EXIT( 0, "ERROR : Something wrong !!\n" );
      return 0.0;
   }
} // FUNCTION : Slice
