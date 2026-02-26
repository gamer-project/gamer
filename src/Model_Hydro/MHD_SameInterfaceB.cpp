#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_SameInterfaceB
// Description :  Ensure that adjacent patches on the same level (i.e., sibling patches) have *exactly*
//                the same B field and B energy on their shared interfaces
//                --> Even round-off errors are the same
//
// Note        :  1. Applied to both real and buffer patches
//                2. Invoked by Init_GAMER() and Main()
//                3. Controlled by the runtime option OPT__SAME_INTERFACE_B
//                4. Always use the B field on the +x/+y/+z sides to overwrite that on the -x/-y/-z sides
//                5. Mainly for debugging purposes since this consistency should already be guaranteed
//                   even when disabling OPT__SAME_INTERFACE_B
//                6. The following two approaches are equivalent. We currently adopt approach (1).
//                   --> (1) Iterate over both real and buffer patches. In this case, there is no need to call
//                           Buf_GetBufferData() afterward to exchange buffer-patch data.
//                       (2) Iterate over real patches only. In this case, we must call Buf_GetBufferData()
//                           afterward to exchange Flu_ParaBuf buffer-patch data for both the magnetic field and energy
//                7. Possible optimizations:
//                   --> (1) When calling MHD_SameInterfaceB() in EvolveLevel(), iterate over real patches only,
//                           since Buf_GetBufferData() is always called in EvolveLevel() anyway.
//                       (2) Related to the above, when calling Buf_GetBufferData() before MHD_SameInterfaceB()
//                           in EvolveLevel(), we could exchange only 0 instead of Flu_ParaBuf data for the
//                           magnetic field and skip exchanging the energy, since buffer patches do not need to be
//                           updated inside MHD_SameInterfaceB() here.
//
// Parameter   :  lv    : AMR level
//                FluSg : Sandglass to access/store the fluid   data
//                MagSg : Sandglass to access/store the B field data
//
// Return      :  Longitudinal B field on the interfaces between sibling patches
//-------------------------------------------------------------------------------------------------------
void MHD_SameInterfaceB( const int lv, const int FluSg, const int MagSg )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );
#  endif

// start of OpenMP parallel region
#  pragma omp parallel
   {
      real ***Emag_old = NULL;

      Aux_AllocateArray3D( Emag_old, 3, PS1, PS1 );

// iterate over all real and buffer patches
#  pragma omp for schedule( runtime )
   for (int PID=0; PID<amr->num[lv]; PID++)
   {
//    use the B field on the +x/+y/+z sides to overwrite that on the -x/-y/-z sides
//    --> skip s=1/3/5
      for (int s=0; s<6; s+=2)
      {
         const int d      = s/2;
         const int SibPID = amr->patch[0][lv][PID]->sibling[s];

         if ( SibPID < 0 )   continue;

//       some buffer patches may have magnetic == NULL --> skip them
         if ( amr->patch[MagSg][lv][   PID]->magnetic == NULL  ||
              amr->patch[MagSg][lv][SibPID]->magnetic == NULL     )   continue;

         int ii, jj, kk;

         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            switch ( d )
            {
               case 0: ii = 0;  jj = i;  kk = j;  break;
               case 1: ii = j;  jj = 0;  kk = i;  break;
               case 2: ii = i;  jj = j;  kk = 0;  break;
            }

            Emag_old[d][j][i] = MHD_GetCellCenteredBEnergyInPatch( lv, PID, ii, jj, kk, MagSg );
         }

         MHD_CopyPatchInterfaceBField( lv, SibPID, s+1, MagSg );

         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            switch ( d )
            {
               case 0: ii = 0;  jj = i;  kk = j;  break;
               case 1: ii = j;  jj = 0;  kk = i;  break;
               case 2: ii = i;  jj = j;  kk = 0;  break;
            }

            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, ii, jj, kk, MagSg );

            amr->patch[FluSg][lv][PID]->fluid[ENGY][kk][jj][ii] += (Emag - Emag_old[d][j][i]);
         }
      } // for (int s=0; s<6; s+=2)
   } // for (int PID=0; PID<amr->num[lv]; PID++)

   Aux_DeallocateArray3D( Emag_old );

   } // end of OpenMP parallel region

} // FUNCTION : MHD_SameInterfaceB



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
