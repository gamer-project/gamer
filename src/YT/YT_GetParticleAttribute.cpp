#include "GAMER.h"

#ifdef SUPPORT_LIBYT

#ifdef PARTICLE

#ifdef LIBYT_USE_PATCH_GROUP
//-------------------------------------------------------------------------------------------------------
// Function    :  Get_ParticleAttribute_PatchGroup
// Description :  Get the particle attribute in patch group.
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. This function will be called if yt inline-analysis needs the data.
//                4. Searching of GID can be optimized.
//                5. There should be 8 grids in a grid patch. Record all their pids, then merge data to
//                   raw_data.
//
// Parameter   :  long  gid     : Grid GID
//                char *attr    : Request particle attribute
//                void *data    : Data pointer, write the results to here.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Get_ParticleAttribute_PatchGroup(int list_len, long *list_gid, char *attr, yt_array *data_array){
    // Get attribute index in GAMER
    FieldIdx_t ParAttr_Idx = GetParticleAttributeIndex( attr, CHECK_ON );

    // loop through list_gid
    for(int lid=0; lid<list_len; lid++){
        // Parse level and pid from gid
        int level = 0;
        int PID0;
        for(int lv=1; lv<NLEVEL; lv++){
            if( list_gid[lid] < YT_GID_Offset[lv] ) break;
            level = lv;
        }
        for(int PID=0; PID<(amr->NPatchComma[level][1]); PID+=8){
            if ( amr->patch[0][level][PID]->libyt_GID == list_gid[lid] ){
                PID0 = PID;
                break;
            }
        }
        // write data
        long  ParID;
        long  data_idx = 0;
        for (int i=0; i<8; i++){ // run through 8 patches
            for (int p=0; p<amr->patch[0][level][PID0 + i]->NPar; p++){ // run through particle data in one PID
                ParID = amr->patch[0][level][PID0 + i]->ParList[p];
                ((real *) data_array[lid].data_ptr)[data_idx] = amr->Par->Attribute[ParAttr_Idx][ParID];
                data_idx += 1;
            }
        }
    }
}
#else // #ifdef LIBYT_USE_PATCH_GROUP
//-------------------------------------------------------------------------------------------------------
// Function    :  Get_ParticleAttribute
// Description :  Get the particle attribute
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. This function will be called if yt inline-analysis needs the data.
//                4. Searching of GID can be optimized.
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
        // parse level and pid.
        int level = 0;
        int pid;
        for(int lv=1; lv<NLEVEL; lv++){
            if( list_gid[lid] < YT_GID_Offset[lv] ) break;
            level = lv;
        }
        for(int PID=0; PID<(amr->NPatchComma[level][1]); PID++){
            if( amr->patch[0][level][PID]->libyt_GID == list_gid[lid] ){
                pid = PID;
                break;
            }
        }
        // write data.
        long ParID;
        for(int p=0; p<amr->patch[0][level][pid]->NPar; p++){
            ParID = amr->patch[0][level][pid]->ParList[p];
            ((real *) data_array[lid].data_ptr)[p] = amr->Par->Attribute[ParAttr_Idx][ParID];
        }
    }
}
#endif // #ifdef LIBYT_USE_PATCH_GROUP

#endif // #ifdef PARTICLE

#endif // #ifdef SUPPORT_LIBYT
