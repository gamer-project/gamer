#include "GAMER.h"

#ifdef SUPPORT_LIBYT

#ifdef PARTICLE

#ifdef LIBYT_USE_PATCH_GROUP
//-------------------------------------------------------------------------------------------------------
// Function    :  Get_ParticleAttribute
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
void Get_ParticleAttribute_PatchGroup(long gid, char *attr, void *raw_data){
    // Parse lv and pid from gid
    int level;
    int pid[8];
    int count = 0;
    for (int lv=0; lv<NLEVEL; lv++){
        for (int PID=0; PID<(amr->NPatchComma[lv][1]); PID+=8){
            if ( amr->patch[0][lv][PID]->libyt_GID == gid ){
                level = lv;
                pid[count] = PID;
                count += 1;
            }
            if (count >= 8) break;
        }
    }

    // Get attribute index in GAMER
    FieldIdx_t ParAttr_Idx = GetParticleAttributeIndex( attr, CHECK_ON );

    // Return the data
    real *data  = (real *) raw_data;
    long  ParID;
    long  data_idx = 0;
    for (int i=0; i<8; i++){ // run through pid
        for (int p=0; p<amr->patch[0][level][pid[i]]->NPar; p++){ // run through particle data in one pid
            ParID          = amr->patch[0][level][pid[i]]->ParList[p];
            data[data_idx] = amr->Par->Attribute[ParAttr_Idx][ParID];
            data_idx += 1;
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
void Get_ParticleAttribute(long gid, char *attr, void *raw_data){
    // Parse lv and pid from gid
    int level;
    int pid;
    for (int lv=0; lv<NLEVEL; lv++){
        for (int PID=0; PID<(amr->NPatchComma[lv][1]); PID++){
            if ( amr->patch[0][lv][PID]->libyt_GID == gid ){
                level = lv;
                pid = PID;
                break;
            }
        }
    }

    // Get attribute index in GAMER
    FieldIdx_t ParAttr_Idx = GetParticleAttributeIndex( attr, CHECK_ON );

    // Return the data
    real *data  = (real *) raw_data;
    long  ParID;
    for (int p=0; p<amr->patch[0][level][pid]->NPar; p++){
        ParID   = amr->patch[0][level][pid]->ParList[p];
        data[p] = amr->Par->Attribute[ParAttr_Idx][ParID];
    }
}
#endif // #ifdef LIBYT_USE_PATCH_GROUP

#endif // #ifdef PARTICLE

#endif // #ifdef SUPPORT_LIBYT
