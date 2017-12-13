#include "Copyright.h"
#include "GAMER.h"
#include "TestProb.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )



extern double Bondi_MassBH;
extern double Bondi_Soften_R;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc_Bondi
// Description :  Set the array "ExtAcc_AuxArray" used by the external acceration routines
//                "CUPOT_ExternalAcc.cu / CPU_ExternalAcc.cpp"
//
// Note        :  1. Enabled by the runtime option "OPT__GRAVITY_TYPE == 2/3"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc_Bondi()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// ExtAcc_AuxArray has the size of EXT_ACC_NAUX_MAX defined in CUPOT.h (default = 10)
   ExtAcc_AuxArray[0] = 0.5*amr->BoxSize[0];
   ExtAcc_AuxArray[1] = 0.5*amr->BoxSize[1];
   ExtAcc_AuxArray[2] = 0.5*amr->BoxSize[2];
   ExtAcc_AuxArray[3] = NEWTON_G*Bondi_MassBH;
   ExtAcc_AuxArray[4] = Bondi_Soften_R;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ExternalAcc_Bondi



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
