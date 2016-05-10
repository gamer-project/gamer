#include "Copyright.h"
#include "GAMER.h"

#ifdef GRAVITY


#include "CUPOT.h"
double ExtPot_AuxArray[EXT_POT_NAUX_MAX];




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


   /*
   const double M  = 1.2e1;
   const double GM = NEWTON_G*M;

// ExtPot_AuxArray has the size of EXT_POT_NAUX_MAX defined in CUPOT.h (default = 10)
   ExtPot_AuxArray[0] = 0.5*amr->BoxSize[0];
   ExtPot_AuxArray[1] = 0.5*amr->BoxSize[1];
   ExtPot_AuxArray[2] = 0.5*amr->BoxSize[2];
   ExtPot_AuxArray[3] = GM;
   */


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ExternalPot



#endif // #ifdef GRAVITY
