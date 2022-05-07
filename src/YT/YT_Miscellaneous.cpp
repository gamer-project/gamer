#include "GAMER.h"
#include <iostream>

#ifdef SUPPORT_LIBYT

//-------------------------------------------------------------------------------------------------------
// Function    :  YT_GetPID
// Description :  Get PID from YT passed in gid.
//
// Note        :  1. Get YT passed in gid, and then write its corresponding level and PID at *level and *PID.
//                2. Whether or not using LIBYT_USE_PATCH_GROUP would work.
//
// Parameter   :  gid      : length of list_gid
//                level    : a list of grid id to prepare.
//                PID      : target field.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
inline void YT_GetPID(long &gid, int *level, int *PID){
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