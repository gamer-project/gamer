#include "GAMER.h"

#ifdef GRAVITY


#include "CUPOT.h"
real ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];

extern real TwoParOrbit_M;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc 
// Description :  Initialize the external potential routines "CUPOT_ExternalAcc.cu / CPU_ExternalAcc.cpp"
//
// Note        :  Fill in the array "ExtAcc_AuxArray" here 
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// ExtAcc_AuxArray has the size of EXT_ACC_NAUX_MAX defined in CUPOT.h (default = 10)
   ExtAcc_AuxArray[0] = (real)0.5*amr->BoxSize[0];
   ExtAcc_AuxArray[1] = (real)0.5*amr->BoxSize[1];
   ExtAcc_AuxArray[2] = (real)0.5*amr->BoxSize[2];
   ExtAcc_AuxArray[3] = (real)0.25*NEWTON_G*TwoParOrbit_M;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ExternalAcc



#endif // #ifdef GRAVITY
