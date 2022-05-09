#include "GAMER.h"

#ifdef SUPPORT_LIBYT

//-------------------------------------------------------------------------------------------------------
// Function    :  YT_GetPID
// Description :  Get PID from YT passed in gid.
//
// Note        :  1. Get YT passed in gid, and then write its corresponding level and PID at *level and *PID.
//                2. Whether or not using LIBYT_USE_PATCH_GROUP would work.
//
// Parameter   :  gid      : gid YT passed in.
//                level    : fill in level of this gid.
//                PID      : fill in PID of this gid.
//
// Return      :  *level, *PID
//-------------------------------------------------------------------------------------------------------
void YT_GetPID(const long gid, int *level, int *PID){
#ifdef GAMER_DEBUG
    long NPatchAllLv = 0;
    for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];
    if ( gid < 0  ||  gid >= NPatchAllLv )  Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", gid, NPatchAllLv-1 );
#endif
    *level = 0;
    for(int lv=1; lv<NLEVEL; lv++){
#ifdef  LIBYT_USE_PATCH_GROUP
        if( gid < (YT_GID_Offset[lv] / 8) ) break;
#else
        if( gid <  YT_GID_Offset[lv] )      break;
#endif // #ifdef  LIBYT_USE_PATCH_GROUP
        *level = lv;
    }
#ifdef LIBYT_USE_PATCH_GROUP
    *PID = gid * 8 - YT_GID_Offset[*level];
#else
    *PID = gid - YT_GID_Offset[*level];
#endif // #ifdef  LIBYT_USE_PATCH_GROUP
}

#endif // #ifdef SUPPORT_LIBYT