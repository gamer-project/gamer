#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_SameInterfaceB
// Description :  Ensure that adjacent patches on the same level (i.e., sibling patches) have *exactly*
//                the same B field on their shared interfaces
//                --> Even round-off errors are the same
//
// Note        :  1. Applied to both real and buffer patches
//                2. Invoked by Init_GAMER() and Main()
//                3. Controlled by the runtime option OPT__SAME_INTERFACE_B
//                4. Always use the B field on the +x/+y/+z sides to overwrite that on the -x/-y/-z sides
//                5. Mainly for debugging purposes since this consistency should already be guaranteed
//                   even when disabling OPT__SAME_INTERFACE_B
//
// Parameter   :  lv : AMR level
//
// Return      :  Longitudinal B field on the interfaces between sibling patches
//-------------------------------------------------------------------------------------------------------
void MHD_SameInterfaceB( const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );
#  endif


   const int MagSg = amr->MagSg[lv];

// iterate over all real and buffer patches
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->num[lv]; PID++)
   {
//    use the B field on the +x/+y/+z sides to overwrite that on the -x/-y/-z sides
//    --> skip s=1/3/5
      for (int s=0; s<6; s+=2)
      {
         const int SibPID = amr->patch[0][lv][PID]->sibling[s];

//       some buffer patches may have magnetic == NULL --> skip them
         if ( SibPID >= 0  &&
              amr->patch[MagSg][lv][   PID]->magnetic != NULL  &&
              amr->patch[MagSg][lv][SibPID]->magnetic != NULL     )
            MHD_CopyPatchInterfaceBField( lv, PID, s, MagSg );
      }
   } // for (int PID=0; PID<amr->num[lv]; PID++)

} // FUNCTION : MHD_SameInterfaceB



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
