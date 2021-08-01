#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_CopyPatchInterfaceBField
// Description :  Copy the longitudinal B field from the common interface of one patch to its sibling patch
//
// Note        :  1. Copy data from [lv][PID] to [lv][PID]->sibling[SibID]
//                2. Invoked by Flu_FixUp_Restrict() and LB_GetBufferData()
//                3. Do nothing if the target sibling patch does not exist
//
// Parameter   :  lv    : AMR level
//                PID   : Source patch ID
//                SibID : Sibling direction (0~5)
//                MagSg : B field sandglass at level lv
//
// Return      :  Longitudinal B field of the patch [lv][PID]->sibling[SibID]
//-------------------------------------------------------------------------------------------------------
void MHD_CopyPatchInterfaceBField( const int lv, const int PID, const int SibID, const int MagSg )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (amr->num[%d] = %d) !!\n", PID, lv, amr->num[lv] );

   if ( SibID < 0  ||  SibID > 5 )
      Aux_Error( ERROR_INFO, "incorrect SibID = %d !!\n", SibID );

   if ( MagSg != 0  &&  MagSg != 1 )
      Aux_Error( ERROR_INFO, "incorrect MagSg = %d !!\n", MagSg );
#  endif


   const int MirrorSib[6]      = { 1, 0, 3, 2, 5, 4 };
   const int Bidx_offset[6]    = { 0, PS1, 0, SQR(PS1), 0, CUBE(PS1) };       // array offsets of the longitudinal B field on 6 faces
   const int Bidx_stride[3][2] = { PS1P1, PS1P1*PS1, 1, PS1P1*PS1, 1, PS1 };  // array strides along the transverse directions of Bx/y/z
   const int Bdir              = SibID/2;                                     // Bx/y/z = 0/1/2
   const int Bdidx_m           = Bidx_stride[Bdir][1];
   const int Bdidx_n           = Bidx_stride[Bdir][0];
   const int SibPID            = amr->patch[0][lv][PID]->sibling[SibID];


// do nothing if there is no sibling patch in the target direction
   if ( SibPID < 0 )
   {
#     ifdef GAMER_DEBUG
      Aux_Message( stderr, "WARNING : SibPID = %d < 0 (lv %d, PID %d, SibID %d) !!\n", SibPID, lv, PID, SibID );
#     endif

      return;
   }


// copy data
#  ifdef GAMER_DEBUG
   if ( amr->patch[MagSg][lv][   PID]->magnetic[Bdir] == NULL )
      Aux_Error( ERROR_INFO, "amr->patch[%d][%d][%d]->magnetic[%d] == NULL !!\n", MagSg, lv,    PID, Bdir );

   if ( amr->patch[MagSg][lv][SibPID]->magnetic[Bdir] == NULL )
      Aux_Error( ERROR_INFO, "amr->patch[%d][%d][%d]->magnetic[%d] == NULL !!\n", MagSg, lv, SibPID, Bdir );
#  endif

   const real    *MagPtr0 = amr->patch[MagSg][lv][   PID]->magnetic[Bdir] + Bidx_offset[           SibID  ];
         real *SibMagPtr0 = amr->patch[MagSg][lv][SibPID]->magnetic[Bdir] + Bidx_offset[ MirrorSib[SibID] ];

   for (int m=0; m<PS1; m++)
   {
      const real    *MagPtr =    MagPtr0 + m*Bdidx_m;
            real *SibMagPtr = SibMagPtr0 + m*Bdidx_m;

      for (int n=0; n<PS1; n++)  SibMagPtr[ n*Bdidx_n ] = MagPtr[ n*Bdidx_n ];
   }

} // FUNCTION : MHD_CopyPatchInterfaceBField



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
