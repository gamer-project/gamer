#include "GAMER.h"



#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )

//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Sync
// Description :  Sync amr->use_wave_flag across all GAMER ranks
//
// Note        :  None
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Flag_Sync( const int lv )
{

   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

    int recv;
    int send = amr->use_wave_flag[lv];

    MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    amr->use_wave_flag[lv] = recv;
} // FUNCTION : Flag_Sync

#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )