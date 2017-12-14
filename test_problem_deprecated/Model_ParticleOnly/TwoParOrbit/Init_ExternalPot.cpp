#include "GAMER.h"

#ifdef GRAVITY


#include "CUPOT.h"
real ExtPot_AuxArray[EXT_POT_NAUX_MAX];

extern real TwoParOrbit_M;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalPot 
// Description :  Initialize the external potential routines "CUPOT_ExternalPot.cu / CPU_ExternalPot.cpp"
//
// Note        :  Fill in the array "ExtPot_AuxArray" here 
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_ExternalPot()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// ExtPot_AuxArray has the size of EXT_POT_NAUX_MAX defined in CUPOT.h (default = 10)
   ExtPot_AuxArray[0] = (real)0.5*amr->BoxSize[0];
   ExtPot_AuxArray[1] = (real)0.5*amr->BoxSize[1];
   ExtPot_AuxArray[2] = (real)0.5*amr->BoxSize[2];
   ExtPot_AuxArray[3] = (real)0.25*NEWTON_G*TwoParOrbit_M;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ExternalPot



#endif // #ifdef GRAVITY
