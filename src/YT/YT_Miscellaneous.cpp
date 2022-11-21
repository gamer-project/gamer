#include "GAMER.h"

#ifdef SUPPORT_LIBYT


//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetPID
// Description :  Convert GID to local PID
//
// Note        :  1. Calculate PID and level from PID.
//
// Parameter   :  GID        : GID to convert
//                level      : fill in level of GID.
//                PID        : fill in PID of GID
//                GID_Offset : table with GID offsets on rank
//
// Return      :  *level, *PID
//-------------------------------------------------------------------------------------------------------
void LB_GetPID(long GID, int& level, int& PID, int* GID_Offset) {
#   ifdef GAMER_DEBUG
    long NPatchAllLv = 0;
    for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];
    if ( GID < 0  ||  GID >= NPatchAllLv )  Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", GID, NPatchAllLv-1 );
#   endif

   level = 0;

   for(int lv = 1; lv < NLEVEL; lv++) {
      if ( GID < GID_Offset[lv] )      
        break;
      level = lv;
   }

   PID = GID - GID_Offset[level];
}


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
void YT_GetPID(const long gid, int *level, int *PID) {
#   ifdef  LIBYT_USE_PATCH_GROUP
    LB_GetPID( 8 * gid, *level, *PID, YT_GID_Offset ); 
#   else
    LB_GetPID(     gid, *level, *PID, YT_GID_Offset );
#   endif // # ifdef  LIBYT_USE_PATCH_GROUP
}
#endif // #ifdef SUPPORT_LIBYT

