#include "GAMER.h"

#ifdef SUPPORT_LIBYT

//-------------------------------------------------------------------------------------------------------
// Function    :  YT_GetPID
// Description :  Get PID from YT passed in gid.
//
// Note        :  1. Get YT passed in gid, and then write its corresponding level and PID at *level and *PID.
//                2. Whether or not using LIBYT_USE_PATCH_GROUP would work.
//                3. The global variable YT_GID_Offset must be filled with values before calling YT_GetPID
//
// Parameter   :  gid      : gid YT passed in.
//                level    : fill in level of this gid.
//                PID      : fill in PID of this gid.
//
// Return      :  *level, *PID
//-------------------------------------------------------------------------------------------------------
void YT_GetPID( const long gid, int *level, int *PID )
{
#   ifdef LIBYT_USE_PATCH_GROUP
    LB_GetPID( 8 * gid, *level, *PID, YT_GID_Offset );
#   else
    LB_GetPID(     gid, *level, *PID, YT_GID_Offset );
#   endif // #ifdef LIBYT_USE_PATCH_GROUP

} // FUNCTION : YT_GetPID
#endif // #ifdef SUPPORT_LIBYT
