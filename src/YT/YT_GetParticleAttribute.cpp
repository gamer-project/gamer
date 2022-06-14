#include "GAMER.h"

#ifdef SUPPORT_LIBYT

void YT_GetPID(const long gid, int *level, int *PID);

#ifdef PARTICLE


//-------------------------------------------------------------------------------------------------------
// Function    :  Get_ParticleAttribute
// Description :  Get the particle attribute.
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. This function will be called if yt inline-analysis needs the data.
//
// Parameter   :  long  gid     : Grid GID
//                char *attr    : Request particle attribute
//                void *data    : Data pointer, write the results to here.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Get_ParticleAttribute(int list_len, long *list_gid, char *attr, yt_array *data_array){
    // Get attribute index in GAMER
    FieldIdx_t ParAttr_Idx = GetParticleAttributeIndex( attr, CHECK_ON );

    // loop through list_gid
    for(int lid=0; lid<list_len; lid++){
        // Parse level and pid from gid
        int level, PID0;
        YT_GetPID( list_gid[lid], &level, &PID0 );

        // write data
#ifdef  LIBYT_USE_PATCH_GROUP
        long  ParID;
        long  data_idx = 0;
        for (int i=0; i<8; i++){ // run through 8 patches
            for (int p=0; p<amr->patch[0][level][PID0 + i]->NPar; p++){ // run through particle data in one PID
                ParID = amr->patch[0][level][PID0 + i]->ParList[p];
                ((real *) data_array[lid].data_ptr)[data_idx] = amr->Par->Attribute[ParAttr_Idx][ParID];
                data_idx += 1;
            }
        }
#else
        long ParID;
        for(int p=0; p<amr->patch[0][level][PID0]->NPar; p++){
            ParID = amr->patch[0][level][PID0]->ParList[p];
            ((real *) data_array[lid].data_ptr)[p] = amr->Par->Attribute[ParAttr_Idx][ParID];
        }
#endif // #ifdef LIBYT_USE_PATCH_GROUP
    }
}

#endif // #ifdef PARTICLE

#endif // #ifdef SUPPORT_LIBYT
