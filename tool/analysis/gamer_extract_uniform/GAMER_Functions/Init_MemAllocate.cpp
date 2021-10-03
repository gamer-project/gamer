#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate
// Description :  Allocate memory for several global arrays
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate()
{

// a. allocate the BaseP
   const int NPatch1D[3] = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };
   BaseP = new int [ NPatch1D[0]*NPatch1D[1]*NPatch1D[2] ];

}



