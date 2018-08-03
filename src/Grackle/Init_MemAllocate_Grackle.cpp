#include "GAMER.h"

#ifdef SUPPORT_GRACKLE


// global variables for accessing h_Che_Array[]
// --> also used by Grackle_Prepare.cpp and Grackle_Close.cpp
// --> they are not declared in "Global.h" simply because they are only used by a few Grackle routines
int Che_NField   = NULL_INT;
int CheIdx_Dens  = Idx_Undefined;
int CheIdx_sEint = Idx_Undefined;
int CheIdx_Ek    = Idx_Undefined;
int CheIdx_e     = Idx_Undefined;
int CheIdx_HI    = Idx_Undefined;
int CheIdx_HII   = Idx_Undefined;
int CheIdx_HeI   = Idx_Undefined;
int CheIdx_HeII  = Idx_Undefined;
int CheIdx_HeIII = Idx_Undefined;
int CheIdx_HM    = Idx_Undefined;
int CheIdx_H2I   = Idx_Undefined;
int CheIdx_H2II  = Idx_Undefined;
int CheIdx_DI    = Idx_Undefined;
int CheIdx_DII   = Idx_Undefined;
int CheIdx_HDI   = Idx_Undefined;
int CheIdx_Metal = Idx_Undefined;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Grackle
// Description :  Allocate the CPU memory for the Grackle solver
//
// Note        :  1. Work even when GPU is enabled
//                2. Invoked by Init_MemAllocate()
//                3. Also set global variables for accessing h_Che_Array[]
//                   --> Declared on the top of this file
//
// Parameter   :  Che_NPG : Number of patch groups to be evaluated at a time
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Grackle( const int Che_NPG )
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


// set global variables related to h_Che_Array[]
   Che_NField   = 0;

   CheIdx_Dens  = Che_NField ++;
   CheIdx_sEint = Che_NField ++;
   CheIdx_Ek    = Che_NField ++;

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   CheIdx_e     = Che_NField ++;
   CheIdx_HI    = Che_NField ++;
   CheIdx_HII   = Che_NField ++;
   CheIdx_HeI   = Che_NField ++;
   CheIdx_HeII  = Che_NField ++;
   CheIdx_HeIII = Che_NField ++;
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   CheIdx_HM    = Che_NField ++;
   CheIdx_H2I   = Che_NField ++;
   CheIdx_H2II  = Che_NField ++;
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   CheIdx_DI    = Che_NField ++;
   CheIdx_DII   = Che_NField ++;
   CheIdx_HDI   = Che_NField ++;
   }

   if ( GRACKLE_METAL )
   CheIdx_Metal = Che_NField ++;


// allocate the input/output array for the Grackle solver
   for (int t=0; t<2; t++)
      h_Che_Array[t] = new real [ (long)Che_NField*(long)Che_NPG*(long)CUBE(PS2) ];

} // FUNCTION : Init_MemAllocate_Grackle



#endif // #ifdef SUPPORT_GRACKLE
