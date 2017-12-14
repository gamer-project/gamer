#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Idx1D2Idx3D
// Description :  Transform a 1D index to 3D indices 
//
// Note        :  Overloaded with different types 
//                --> Explicit template instantiation is put in the end of this file
//
// Parameter   :  Size  : Size in all three spatial directions
//                Idx1D : Input 1D index 
//                Idx3D : Output 3D index
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Mis_Idx1D2Idx3D( const int Size[], const T Idx1D, int Idx3D[] )
{

   const T SizeXY = (T)Size[0]*(T)Size[1];

   Idx3D[0] = Idx1D % (T)Size[0];
   Idx3D[1] = Idx1D % SizeXY / (T)Size[0];
   Idx3D[2] = Idx1D / SizeXY;

} // FUNCTION : Mis_Idx1D2Idx3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Idx3D2Idx1D
// Description :  Transform 3D indices to a 1D index
//
// Note        :  The macro "IDX321" does the same job
//                
// Parameter   :  Size  : Size in all three spatial directions
//                Idx3D : Output 3D index
//
// Return      :  1D index 
//-------------------------------------------------------------------------------------------------------
ulong Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] )
{

   return ( (ulong)Idx3D[2]*Size[1] + (ulong)Idx3D[1] )*Size[0] + (ulong)Idx3D[0];

} // FUNCTION : Mis_Idx3D2Idx1D



// explicit template instantiation
template void Mis_Idx1D2Idx3D <int>   ( const int Size[], const int   Idx1D, int Idx3D[] );
template void Mis_Idx1D2Idx3D <long>  ( const int Size[], const long  Idx1D, int Idx3D[] );
template void Mis_Idx1D2Idx3D <ulong> ( const int Size[], const ulong Idx1D, int Idx3D[] );
