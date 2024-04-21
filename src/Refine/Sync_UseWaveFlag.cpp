#include "GAMER.h"

#if ( ELBDM_SCHEME == ELBDM_HYBRID )




//-------------------------------------------------------------------------------------------------------
// Function    :  Sync_UseWaveFlag
// Description :  Synchronize amr->use_wave_flag[] across all MPI ranks
//
// Note        :  None
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Sync_UseWaveFlag( const int lv )
{

   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

    bool recv;
    bool send = amr->use_wave_flag[lv];

    MPI_Allreduce( &send, &recv, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD );

    amr->use_wave_flag[lv] = recv;

} // FUNCTION : Sync_UseWaveFlag



#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
