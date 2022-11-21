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
void YT_GetPID(const long gid, int *level, int *PID) {
#   ifdef  LIBYT_USE_PATCH_GROUP
    LB_GetPID( 8 * gid, *level, *PID, YT_PatchCount ); 
#   else
    LB_GetPID(     gid, *level, *PID, YT_PatchCount );
#   endif // # ifdef  LIBYT_USE_PATCH_GROUP
}
#endif // #ifdef SUPPORT_LIBYT