#include "GAMER.h"


#ifdef SUPPORT_LIBYT

#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Get_ParticleAttribute
// Description :  Get the particle attribute
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
#endif

#endif