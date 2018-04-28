#include "GAMER.h"

#ifdef LOAD_BALANCE



static void SetSiblingInSamePatchGroup( const int lv, const int PID0 );
static void SetSiblingInDiffPatchGroup( const int lv, const int PID0, const int SibPID0, const int SibID,
                                        const bool BothSide );
static void SetSiblingExternal( const int lv, const int NTarget0, const int *TargetPID0 );
static void SetSiblingExternal_CheckPeriodicity( int *Sibling, const int Boundary );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_SiblingSearch
// Description :  Construct the sibling patch relation
//
// Note        :  1. LB_PaddedCr1DList[] and LB_PaddedCr1DList_IdxTable[] at SonLv and FaLv must be properly prepared
//                2. SearchAllPID == true  --> Works on all patches at lv (including real, sibling-buffer
//                                             and father-buffer patches)
//                                == false --> Only works on PID0 recorded in TargetPID0
//
// Parameter   :  lv           : Target refinement level
//                SearchAllPID : Whether to search over all patches at lv or not
//                NInput       : Number of target patches (with LocalID==0) in "TargetPID0"
//                               (useful only if "SearchAllPID == false")
//                TargetPID0   : Lists recording all target patches (with LocalID==0)
//                               (useful only if "SearchAllPID == false")
//-------------------------------------------------------------------------------------------------------
void LB_SiblingSearch( const int lv, const bool SearchAllPID, const int NInput, int *TargetPID0 )
{

   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( !SearchAllPID  &&  NInput != 0  &&  TargetPID0 == NULL )
      Aux_Error( ERROR_INFO, "lv %d, NInput %d, TargetPID0 == NULL !!\n", lv, NInput );


   const int  NPatch              = amr->num[lv];
   const bool BothSide            = ( SearchAllPID ) ? false : true;             // construct relations in both side
   const int  NTarget0            = ( SearchAllPID ) ? NPatch/8 : NInput;
   const int  NSib                = 26;
   const int  NSearch_Max         = NSib*NTarget0;
   const int  Padded              = 1<<NLEVEL;
   const int  BoxNScale_Padded[3] = { amr->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                      amr->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                      amr->BoxScale[2]/PATCH_SIZE + 2*Padded };  //normalized and padded BoxScale
   const int  Scale2              = 2*amr->scale[lv];
   const long dr[3]               = { (long)Scale2,
                                      (long)Scale2*BoxNScale_Padded[0],
                                      (long)Scale2*BoxNScale_Padded[0]*BoxNScale_Padded[1] };

   int   NSearch, Count, *Match, PID0;
   long  Cr1D_Disp[26];

   ulong *SibCr1D_Search   = new ulong [ NSearch_Max ];
   ulong *SibCr1D          = new ulong [ NSearch_Max ];
   int   *SibCr1D_IdxTable = new  int  [ NSearch_Max ];


// nothing to do if there is no target patches
   if ( NTarget0 == 0 )
   {
      delete [] SibCr1D_Search;
      delete [] SibCr1D;
      delete [] SibCr1D_IdxTable;

      return;
   }


// 0. initialize all siblings as -1 and construct the target patch list with LocalID==0 (for SearchAllPID)
   if ( SearchAllPID )
   {
      for (int PID=0; PID<NPatch; PID++)
      for (int s=0; s<NSib; s++)
         amr->patch[0][lv][PID]->sibling[s] = -1;

      TargetPID0 = new int [NTarget0];

      for (int t=0; t<NTarget0; t++)   TargetPID0[t] = 8*t;
   }


// 1. construct the displacement matrix of padded 1D corner coordinates
   Count = 0;

   for (int k=-1; k<=1; k++)
   for (int j=-1; j<=1; j++)
   for (int i=-1; i<=1; i++)
      if ( i != 0  ||  j != 0  ||  k != 0 )  Cr1D_Disp[ Count++ ] = (long)i*dr[0] + (long)j*dr[1] + (long)k*dr[2];


// 2. construct the sorted sibling list (with duplicates)
   Count = 0;

   for (int t=0; t<NTarget0; t++)
   {
      PID0 = TargetPID0[t];

#     ifdef GAMER_DEBUG
      if ( PID0%8 != 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID0 %d is not a multiple of 8 !!\n", lv, PID0 );
#     endif

//###NOTE: Disp = i*dr[0] + j*dr[1] + k*dr[2] can be negative! But it's OK to conduct PaddedCr1D + (ulong)Disp
//         as long as we guarantee "PaddedCr1D + Disp >= 0"
//         --> ulong(Disp) = Disp + UINT_MAX + 1 (if Disp < 0; ==> reduced modulo)
//         --> PaddedCr1D + (ulong)Disp = PaddedCr1D + Disp + UINT_MAX + 1 = PaddedCr1D + Disp + UINT_MAX + 1 - (UINT_MAX + 1)
//                                      = PaddedCr1D + Disp
//             (because PaddedCr1D + Disp >= 0; ==> reduced modulo again)
      for (int s=0; s<NSib; s++)
         SibCr1D[ Count++ ] = amr->patch[0][lv][PID0]->PaddedCr1D + (ulong)Cr1D_Disp[s];
   }

   Mis_Heapsort( NSearch_Max, SibCr1D, SibCr1D_IdxTable );


// 3. prepare the searching list and remove the duplicates
   memcpy( SibCr1D_Search, SibCr1D, NSearch_Max*sizeof(ulong) );

   NSearch = ( NTarget0 > 0 ) ? 1 : 0;

   for (int t=1; t<NSearch_Max; t++)
      if ( SibCr1D_Search[t] != SibCr1D_Search[t-1] )   SibCr1D_Search[ NSearch ++ ] = SibCr1D_Search[t];


// 4. matching
   Match = new int [NSearch];
   Mis_Matching_int( NPatch, amr->LB->PaddedCr1DList[lv], NSearch, SibCr1D_Search, Match );


// 5. construct the sibling relation
   const int PGScale = PATCH_SIZE*Scale2;
   const int SibID[3][3][3] = {  { {18, 10, 19}, {14,  4, 16}, {20, 11, 21} },
                                 { { 6,  2,  7}, { 0, -1,  1}, { 8,  3,  9} },
                                 { {22, 12, 23}, {15,  5, 17}, {24, 13, 25} }  };
   ulong Cr1D;
   int   SibPID0, dID[3], Start = 0;
   int  *Cr1, *Cr2;

// 5.1 construct the sibling relation for patches within the same patch group
   for (int t=0; t<NTarget0; t++)   SetSiblingInSamePatchGroup( lv, TargetPID0[t] );


// 5.2 construct the sibling relation for patches in different patch groups
   for (int t=0; t<NSearch; t++)
   {
      if ( Match[t] != -1 )
      {
         Cr1D    = amr->LB->PaddedCr1DList         [lv][ Match[t] ];
         SibPID0 = amr->LB->PaddedCr1DList_IdxTable[lv][ Match[t] ];
         Cr2     = amr->patch[0][lv][SibPID0]->corner;

         for (int m=Start; m<NSearch_Max; m++)
         {
            if ( SibCr1D[m] == Cr1D )
            {
               do
               {
                  PID0 = TargetPID0[ SibCr1D_IdxTable[m] / NSib ];
                  Cr1  = amr->patch[0][lv][PID0]->corner;

                  for (int d=0; d<3; d++)    dID[d] = 1 + ( Cr2[d] - Cr1[d] ) / PGScale;

//                for NLEVEL == 1, buffer patch groups can have sibling PaddedCr1D map to wrong buffer
//                patch groups in the opposite direction (check the note for a more detailed explanation)
#                 if ( NLEVEL == 1 )
                  if (  dID[0]<0 || dID[0]>2 || dID[1]<0 || dID[1]>2  )
                  {
                     m++;
                     continue;
                  }
#                 endif

#                 ifdef GAMER_DEBUG
                  if (  ( NLEVEL != 1 && (dID[0]<0 || dID[0]>2 || dID[1]<0 || dID[1]>2) )
                        || dID[2]<0 || dID[2]>2 || ( dID[0]==1 && dID[1]==1 && dID[2]==1 )  )
                     Aux_Error( ERROR_INFO, "lv %d, PID0 %d, SibPID0 %d, incorrect dID[3]=(%d,%d,%d) !!\n",
                                lv, PID0, SibPID0, dID[0], dID[1], dID[2] );
#                 endif

                  SetSiblingInDiffPatchGroup( lv, PID0, SibPID0, SibID[ dID[2] ][ dID[1] ][ dID[0] ], BothSide );

                  m++;
               }
               while( m < NSearch_Max  &&  SibCr1D[m] == Cr1D );

               Start = m;
               break;
            } // if ( SibCr1D[m] == Cr1D )

#           ifdef GAMER_DEBUG
            if ( m == NSearch_Max - 1 )
               Aux_Error( ERROR_INFO, "lv %d, SibPID0 %d, Cr1D %lu has no matching patch !!\n",
                          lv, SibPID0, Cr1D );
#           endif

         } // for (int m=Start; m<NSearch_Max; m++)
      } // if ( Match[t] != -1 )
   } // for (int t=0; t<NSearch; t++)


// 5.3 set the sibling indices for the patches adjacent to the simulation domain (for non-periodic B.C. only)
   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC  ||
        OPT__BC_FLU[2] != BC_FLU_PERIODIC  ||
        OPT__BC_FLU[4] != BC_FLU_PERIODIC   )   SetSiblingExternal( lv, NTarget0, TargetPID0 );


// check results in debug mode
#  ifdef GAMER_DEBUG
   const int MirrorSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const int PScale        = PATCH_SIZE*amr->scale[lv];

   for (int PID=0; PID<NPatch; PID++)
   for (int s=0; s<NSib; s++)
   {
      int SibPID = amr->patch[0][lv][PID]->sibling[s];

      if ( SibPID >= 0 )
      {
//       check 1: PID's sibling's mirror-sibling = PID
         if ( amr->patch[0][lv][SibPID]->sibling[ MirrorSib[s] ] != PID )
            Aux_Error( ERROR_INFO, "lv %d, PID[%d]->Sib[%d] = %d != SibPID[%d]->MirrorSib[%d] = %d !!\n",
                       lv, PID, s, SibPID, SibPID, MirrorSib[s],
                       amr->patch[0][lv][SibPID]->sibling[ MirrorSib[s] ] );

//       check 2: PID's sibling has correct coordinates
         for (int d=0; d<3; d++)
         {
            if (    amr->patch[0][lv][   PID]->corner[d] + TABLE_01( s, 'x'+d, -PScale, 0, PScale )
                 != amr->patch[0][lv][SibPID]->corner[d]  )
               Aux_Error( ERROR_INFO, "lv %d, sibling %2d, PID %8d, SibPID %8d, dim %d: corner %8d + %8d != %8d\n",
                          lv, s, PID, SibPID, d, amr->patch[0][lv][PID]->corner[d],
                          TABLE_01( s, 'x'+d, -PScale, 0, PScale ), amr->patch[0][lv][SibPID]->corner[d] );
         }
      }
   } // for (PID, s)
#  endif // #ifdef GAMER_DEBUG


// free memory
   delete [] SibCr1D_Search;
   delete [] SibCr1D;
   delete [] SibCr1D_IdxTable;
   delete [] Match;
   if ( SearchAllPID )  delete [] TargetPID0;

} // FUNCTION : LB_SiblingSearch



//-------------------------------------------------------------------------------------------------------
// Function    :  SetSiblingInSamePatchGroup
// Description :  Construct the sibling patch relation for patches within the same patch group
//
// Note        :  Sibling relation of patches in different patch groups is constructed by
//                SetSiblingInDiffPatchGroup()
//
// Parameter   :  lv   : Target refinement level
//                PID0 : Index of the patch with LocalID==0
//-------------------------------------------------------------------------------------------------------
void SetSiblingInSamePatchGroup( const int lv, const int PID0 )
{

// check
#  ifdef GAMER_DEBUG
   if ( PID0%8 != 0 )   Aux_Error( ERROR_INFO, "lv %d, PID0 %d is not a multiple of 8 !!\n", lv, PID0 );
#  endif


   int PID[8];
   for (int LocalID=0; LocalID<8; LocalID++)    PID[LocalID] = PID0 + LocalID;

// LocalID == 0
   amr->patch[0][lv][ PID[0] ]->sibling[ 1] = PID[1];
   amr->patch[0][lv][ PID[0] ]->sibling[ 3] = PID[2];
   amr->patch[0][lv][ PID[0] ]->sibling[ 5] = PID[3];
   amr->patch[0][lv][ PID[0] ]->sibling[ 9] = PID[4];
   amr->patch[0][lv][ PID[0] ]->sibling[13] = PID[5];
   amr->patch[0][lv][ PID[0] ]->sibling[17] = PID[6];
   amr->patch[0][lv][ PID[0] ]->sibling[25] = PID[7];

// LocalID == 1
   amr->patch[0][lv][ PID[1] ]->sibling[ 0] = PID[0];
   amr->patch[0][lv][ PID[1] ]->sibling[ 8] = PID[2];
   amr->patch[0][lv][ PID[1] ]->sibling[15] = PID[3];
   amr->patch[0][lv][ PID[1] ]->sibling[ 3] = PID[4];
   amr->patch[0][lv][ PID[1] ]->sibling[24] = PID[5];
   amr->patch[0][lv][ PID[1] ]->sibling[ 5] = PID[6];
   amr->patch[0][lv][ PID[1] ]->sibling[13] = PID[7];

// LocalID == 2
   amr->patch[0][lv][ PID[2] ]->sibling[ 2] = PID[0];
   amr->patch[0][lv][ PID[2] ]->sibling[ 7] = PID[1];
   amr->patch[0][lv][ PID[2] ]->sibling[12] = PID[3];
   amr->patch[0][lv][ PID[2] ]->sibling[ 1] = PID[4];
   amr->patch[0][lv][ PID[2] ]->sibling[ 5] = PID[5];
   amr->patch[0][lv][ PID[2] ]->sibling[23] = PID[6];
   amr->patch[0][lv][ PID[2] ]->sibling[17] = PID[7];

// LocalID == 3
   amr->patch[0][lv][ PID[3] ]->sibling[ 4] = PID[0];
   amr->patch[0][lv][ PID[3] ]->sibling[16] = PID[1];
   amr->patch[0][lv][ PID[3] ]->sibling[11] = PID[2];
   amr->patch[0][lv][ PID[3] ]->sibling[21] = PID[4];
   amr->patch[0][lv][ PID[3] ]->sibling[ 3] = PID[5];
   amr->patch[0][lv][ PID[3] ]->sibling[ 1] = PID[6];
   amr->patch[0][lv][ PID[3] ]->sibling[ 9] = PID[7];

// LocalID == 4
   amr->patch[0][lv][ PID[4] ]->sibling[ 6] = PID[0];
   amr->patch[0][lv][ PID[4] ]->sibling[ 2] = PID[1];
   amr->patch[0][lv][ PID[4] ]->sibling[ 0] = PID[2];
   amr->patch[0][lv][ PID[4] ]->sibling[22] = PID[3];
   amr->patch[0][lv][ PID[4] ]->sibling[15] = PID[5];
   amr->patch[0][lv][ PID[4] ]->sibling[12] = PID[6];
   amr->patch[0][lv][ PID[4] ]->sibling[ 5] = PID[7];

// LocalID == 5
   amr->patch[0][lv][ PID[5] ]->sibling[10] = PID[0];
   amr->patch[0][lv][ PID[5] ]->sibling[19] = PID[1];
   amr->patch[0][lv][ PID[5] ]->sibling[ 4] = PID[2];
   amr->patch[0][lv][ PID[5] ]->sibling[ 2] = PID[3];
   amr->patch[0][lv][ PID[5] ]->sibling[16] = PID[4];
   amr->patch[0][lv][ PID[5] ]->sibling[ 7] = PID[6];
   amr->patch[0][lv][ PID[5] ]->sibling[ 1] = PID[7];

// LocalID == 6
   amr->patch[0][lv][ PID[6] ]->sibling[14] = PID[0];
   amr->patch[0][lv][ PID[6] ]->sibling[ 4] = PID[1];
   amr->patch[0][lv][ PID[6] ]->sibling[20] = PID[2];
   amr->patch[0][lv][ PID[6] ]->sibling[ 0] = PID[3];
   amr->patch[0][lv][ PID[6] ]->sibling[11] = PID[4];
   amr->patch[0][lv][ PID[6] ]->sibling[ 8] = PID[5];
   amr->patch[0][lv][ PID[6] ]->sibling[ 3] = PID[7];

// LocalID == 7
   amr->patch[0][lv][ PID[7] ]->sibling[18] = PID[0];
   amr->patch[0][lv][ PID[7] ]->sibling[10] = PID[1];
   amr->patch[0][lv][ PID[7] ]->sibling[14] = PID[2];
   amr->patch[0][lv][ PID[7] ]->sibling[ 6] = PID[3];
   amr->patch[0][lv][ PID[7] ]->sibling[ 4] = PID[4];
   amr->patch[0][lv][ PID[7] ]->sibling[ 0] = PID[5];
   amr->patch[0][lv][ PID[7] ]->sibling[ 2] = PID[6];

} // FUNCTION : SetSiblingInSamePatchGroup



//-------------------------------------------------------------------------------------------------------
// Function    :  SetSiblingInDiffPatchGroup
// Description :  Construct the sibling patch relation for patches in different patch groups
//
// Note        :  Sibling relation of patches in the same patch group is constructed by
//                SetSiblingInSamePatchGroup()
//
// Parameter   :  lv       : Target refinement level
//                PID0     : Index of the target patch with LocalID==0
//                SibPID0  : Index of the sibling patch with LocalID==0
//                SibID    : Sibling index
//                BothSide : Construct the relation between two nearby patches at the same time
//-------------------------------------------------------------------------------------------------------
void SetSiblingInDiffPatchGroup( const int lv, const int PID0, const int SibPID0, const int SibID,
                                 const bool BothSide )
{

#  ifdef GAMER_DEBUG
   if (    PID0%8 != 0 )   Aux_Error( ERROR_INFO, "lv %d,    PID0 %d is not a multiple of 8 !!\n", lv,    PID0 );
   if ( SibPID0%8 != 0 )   Aux_Error( ERROR_INFO, "lv %d, SibPID0 %d is not a multiple of 8 !!\n", lv, SibPID0 );
#  endif


   switch ( SibID )
   {
      case  0:
         amr->patch[0][lv][ PID0+0 ]->sibling[ 0] = SibPID0+1;
         amr->patch[0][lv][ PID0+0 ]->sibling[ 8] = SibPID0+4;
         amr->patch[0][lv][ PID0+0 ]->sibling[15] = SibPID0+6;
         amr->patch[0][lv][ PID0+0 ]->sibling[24] = SibPID0+7;

         amr->patch[0][lv][ PID0+2 ]->sibling[ 6] = SibPID0+1;
         amr->patch[0][lv][ PID0+2 ]->sibling[ 0] = SibPID0+4;
         amr->patch[0][lv][ PID0+2 ]->sibling[22] = SibPID0+6;
         amr->patch[0][lv][ PID0+2 ]->sibling[15] = SibPID0+7;

         amr->patch[0][lv][ PID0+3 ]->sibling[14] = SibPID0+1;
         amr->patch[0][lv][ PID0+3 ]->sibling[20] = SibPID0+4;
         amr->patch[0][lv][ PID0+3 ]->sibling[ 0] = SibPID0+6;
         amr->patch[0][lv][ PID0+3 ]->sibling[ 8] = SibPID0+7;

         amr->patch[0][lv][ PID0+5 ]->sibling[18] = SibPID0+1;
         amr->patch[0][lv][ PID0+5 ]->sibling[14] = SibPID0+4;
         amr->patch[0][lv][ PID0+5 ]->sibling[ 6] = SibPID0+6;
         amr->patch[0][lv][ PID0+5 ]->sibling[ 0] = SibPID0+7;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+1 ]->sibling[ 1] = PID0+0;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[ 7] = PID0+0;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[16] = PID0+0;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[19] = PID0+0;

            amr->patch[0][lv][ SibPID0+1 ]->sibling[ 9] = PID0+2;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[ 1] = PID0+2;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[21] = PID0+2;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[16] = PID0+2;

            amr->patch[0][lv][ SibPID0+1 ]->sibling[17] = PID0+3;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[23] = PID0+3;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[ 1] = PID0+3;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[ 7] = PID0+3;

            amr->patch[0][lv][ SibPID0+1 ]->sibling[25] = PID0+5;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[17] = PID0+5;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[ 9] = PID0+5;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[ 1] = PID0+5;
         }
      break;


      case  1:
         amr->patch[0][lv][ PID0+1 ]->sibling[ 1] = SibPID0+0;
         amr->patch[0][lv][ PID0+1 ]->sibling[ 9] = SibPID0+2;
         amr->patch[0][lv][ PID0+1 ]->sibling[17] = SibPID0+3;
         amr->patch[0][lv][ PID0+1 ]->sibling[25] = SibPID0+5;

         amr->patch[0][lv][ PID0+4 ]->sibling[ 7] = SibPID0+0;
         amr->patch[0][lv][ PID0+4 ]->sibling[ 1] = SibPID0+2;
         amr->patch[0][lv][ PID0+4 ]->sibling[23] = SibPID0+3;
         amr->patch[0][lv][ PID0+4 ]->sibling[17] = SibPID0+5;

         amr->patch[0][lv][ PID0+6 ]->sibling[16] = SibPID0+0;
         amr->patch[0][lv][ PID0+6 ]->sibling[21] = SibPID0+2;
         amr->patch[0][lv][ PID0+6 ]->sibling[ 1] = SibPID0+3;
         amr->patch[0][lv][ PID0+6 ]->sibling[ 9] = SibPID0+5;

         amr->patch[0][lv][ PID0+7 ]->sibling[19] = SibPID0+0;
         amr->patch[0][lv][ PID0+7 ]->sibling[16] = SibPID0+2;
         amr->patch[0][lv][ PID0+7 ]->sibling[ 7] = SibPID0+3;
         amr->patch[0][lv][ PID0+7 ]->sibling[ 1] = SibPID0+5;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+0 ]->sibling[ 0] = PID0+1;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[ 6] = PID0+1;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[14] = PID0+1;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[18] = PID0+1;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[ 8] = PID0+4;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[ 0] = PID0+4;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[20] = PID0+4;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[14] = PID0+4;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[15] = PID0+6;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[22] = PID0+6;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[ 0] = PID0+6;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[ 6] = PID0+6;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[24] = PID0+7;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[15] = PID0+7;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[ 8] = PID0+7;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[ 0] = PID0+7;
         }
      break;


      case  2:
         amr->patch[0][lv][ PID0+0 ]->sibling[ 2] = SibPID0+2;
         amr->patch[0][lv][ PID0+0 ]->sibling[ 7] = SibPID0+4;
         amr->patch[0][lv][ PID0+0 ]->sibling[12] = SibPID0+5;
         amr->patch[0][lv][ PID0+0 ]->sibling[23] = SibPID0+7;

         amr->patch[0][lv][ PID0+1 ]->sibling[ 6] = SibPID0+2;
         amr->patch[0][lv][ PID0+1 ]->sibling[ 2] = SibPID0+4;
         amr->patch[0][lv][ PID0+1 ]->sibling[22] = SibPID0+5;
         amr->patch[0][lv][ PID0+1 ]->sibling[12] = SibPID0+7;

         amr->patch[0][lv][ PID0+3 ]->sibling[10] = SibPID0+2;
         amr->patch[0][lv][ PID0+3 ]->sibling[19] = SibPID0+4;
         amr->patch[0][lv][ PID0+3 ]->sibling[ 2] = SibPID0+5;
         amr->patch[0][lv][ PID0+3 ]->sibling[ 7] = SibPID0+7;

         amr->patch[0][lv][ PID0+6 ]->sibling[18] = SibPID0+2;
         amr->patch[0][lv][ PID0+6 ]->sibling[10] = SibPID0+4;
         amr->patch[0][lv][ PID0+6 ]->sibling[ 6] = SibPID0+5;
         amr->patch[0][lv][ PID0+6 ]->sibling[ 2] = SibPID0+7;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+2 ]->sibling[ 3] = PID0+0;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[ 8] = PID0+0;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[11] = PID0+0;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[20] = PID0+0;

            amr->patch[0][lv][ SibPID0+2 ]->sibling[ 9] = PID0+1;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[ 3] = PID0+1;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[21] = PID0+1;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[11] = PID0+1;

            amr->patch[0][lv][ SibPID0+2 ]->sibling[13] = PID0+3;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[24] = PID0+3;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[ 3] = PID0+3;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[ 8] = PID0+3;

            amr->patch[0][lv][ SibPID0+2 ]->sibling[25] = PID0+6;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[13] = PID0+6;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[ 9] = PID0+6;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[ 3] = PID0+6;
         }
      break;


      case  3:
         amr->patch[0][lv][ PID0+2 ]->sibling[ 3] = SibPID0+0;
         amr->patch[0][lv][ PID0+2 ]->sibling[ 9] = SibPID0+1;
         amr->patch[0][lv][ PID0+2 ]->sibling[13] = SibPID0+3;
         amr->patch[0][lv][ PID0+2 ]->sibling[25] = SibPID0+6;

         amr->patch[0][lv][ PID0+4 ]->sibling[ 8] = SibPID0+0;
         amr->patch[0][lv][ PID0+4 ]->sibling[ 3] = SibPID0+1;
         amr->patch[0][lv][ PID0+4 ]->sibling[24] = SibPID0+3;
         amr->patch[0][lv][ PID0+4 ]->sibling[13] = SibPID0+6;

         amr->patch[0][lv][ PID0+5 ]->sibling[11] = SibPID0+0;
         amr->patch[0][lv][ PID0+5 ]->sibling[21] = SibPID0+1;
         amr->patch[0][lv][ PID0+5 ]->sibling[ 3] = SibPID0+3;
         amr->patch[0][lv][ PID0+5 ]->sibling[ 9] = SibPID0+6;

         amr->patch[0][lv][ PID0+7 ]->sibling[20] = SibPID0+0;
         amr->patch[0][lv][ PID0+7 ]->sibling[11] = SibPID0+1;
         amr->patch[0][lv][ PID0+7 ]->sibling[ 8] = SibPID0+3;
         amr->patch[0][lv][ PID0+7 ]->sibling[ 3] = SibPID0+6;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+0 ]->sibling[ 2] = PID0+2;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[ 6] = PID0+2;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[10] = PID0+2;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[18] = PID0+2;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[ 7] = PID0+4;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[ 2] = PID0+4;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[19] = PID0+4;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[10] = PID0+4;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[12] = PID0+5;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[22] = PID0+5;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[ 2] = PID0+5;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[ 6] = PID0+5;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[23] = PID0+7;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[12] = PID0+7;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[ 7] = PID0+7;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[ 2] = PID0+7;
         }
      break;


      case  4:
         amr->patch[0][lv][ PID0+0 ]->sibling[ 4] = SibPID0+3;
         amr->patch[0][lv][ PID0+0 ]->sibling[11] = SibPID0+5;
         amr->patch[0][lv][ PID0+0 ]->sibling[16] = SibPID0+6;
         amr->patch[0][lv][ PID0+0 ]->sibling[21] = SibPID0+7;

         amr->patch[0][lv][ PID0+1 ]->sibling[14] = SibPID0+3;
         amr->patch[0][lv][ PID0+1 ]->sibling[20] = SibPID0+5;
         amr->patch[0][lv][ PID0+1 ]->sibling[ 4] = SibPID0+6;
         amr->patch[0][lv][ PID0+1 ]->sibling[11] = SibPID0+7;

         amr->patch[0][lv][ PID0+2 ]->sibling[10] = SibPID0+3;
         amr->patch[0][lv][ PID0+2 ]->sibling[ 4] = SibPID0+5;
         amr->patch[0][lv][ PID0+2 ]->sibling[19] = SibPID0+6;
         amr->patch[0][lv][ PID0+2 ]->sibling[16] = SibPID0+7;

         amr->patch[0][lv][ PID0+4 ]->sibling[18] = SibPID0+3;
         amr->patch[0][lv][ PID0+4 ]->sibling[14] = SibPID0+5;
         amr->patch[0][lv][ PID0+4 ]->sibling[10] = SibPID0+6;
         amr->patch[0][lv][ PID0+4 ]->sibling[ 4] = SibPID0+7;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+3 ]->sibling[ 5] = PID0+0;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[12] = PID0+0;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[15] = PID0+0;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[22] = PID0+0;

            amr->patch[0][lv][ SibPID0+3 ]->sibling[17] = PID0+1;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[23] = PID0+1;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[ 5] = PID0+1;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[12] = PID0+1;

            amr->patch[0][lv][ SibPID0+3 ]->sibling[13] = PID0+2;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[ 5] = PID0+2;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[24] = PID0+2;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[15] = PID0+2;

            amr->patch[0][lv][ SibPID0+3 ]->sibling[25] = PID0+4;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[17] = PID0+4;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[13] = PID0+4;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[ 5] = PID0+4;
         }
      break;


      case  5:
         amr->patch[0][lv][ PID0+3 ]->sibling[ 5] = SibPID0+0;
         amr->patch[0][lv][ PID0+3 ]->sibling[17] = SibPID0+1;
         amr->patch[0][lv][ PID0+3 ]->sibling[13] = SibPID0+2;
         amr->patch[0][lv][ PID0+3 ]->sibling[25] = SibPID0+4;

         amr->patch[0][lv][ PID0+5 ]->sibling[12] = SibPID0+0;
         amr->patch[0][lv][ PID0+5 ]->sibling[23] = SibPID0+1;
         amr->patch[0][lv][ PID0+5 ]->sibling[ 5] = SibPID0+2;
         amr->patch[0][lv][ PID0+5 ]->sibling[17] = SibPID0+4;

         amr->patch[0][lv][ PID0+6 ]->sibling[15] = SibPID0+0;
         amr->patch[0][lv][ PID0+6 ]->sibling[ 5] = SibPID0+1;
         amr->patch[0][lv][ PID0+6 ]->sibling[24] = SibPID0+2;
         amr->patch[0][lv][ PID0+6 ]->sibling[13] = SibPID0+4;

         amr->patch[0][lv][ PID0+7 ]->sibling[22] = SibPID0+0;
         amr->patch[0][lv][ PID0+7 ]->sibling[12] = SibPID0+1;
         amr->patch[0][lv][ PID0+7 ]->sibling[15] = SibPID0+2;
         amr->patch[0][lv][ PID0+7 ]->sibling[ 5] = SibPID0+4;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+0 ]->sibling[ 4] = PID0+3;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[14] = PID0+3;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[10] = PID0+3;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[18] = PID0+3;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[11] = PID0+5;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[20] = PID0+5;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[ 4] = PID0+5;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[14] = PID0+5;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[16] = PID0+6;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[ 4] = PID0+6;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[19] = PID0+6;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[10] = PID0+6;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[21] = PID0+7;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[11] = PID0+7;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[16] = PID0+7;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[ 4] = PID0+7;
         }
      break;


      case  6:
         amr->patch[0][lv][ PID0+0 ]->sibling[ 6] = SibPID0+4;
         amr->patch[0][lv][ PID0+0 ]->sibling[22] = SibPID0+7;

         amr->patch[0][lv][ PID0+3 ]->sibling[18] = SibPID0+4;
         amr->patch[0][lv][ PID0+3 ]->sibling[ 6] = SibPID0+7;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+4 ]->sibling[ 9] = PID0+0;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[21] = PID0+0;

            amr->patch[0][lv][ SibPID0+4 ]->sibling[25] = PID0+3;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[ 9] = PID0+3;
         }
      break;


      case  7:
         amr->patch[0][lv][ PID0+1 ]->sibling[ 7] = SibPID0+2;
         amr->patch[0][lv][ PID0+1 ]->sibling[23] = SibPID0+5;

         amr->patch[0][lv][ PID0+6 ]->sibling[19] = SibPID0+2;
         amr->patch[0][lv][ PID0+6 ]->sibling[ 7] = SibPID0+5;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+2 ]->sibling[ 8] = PID0+1;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[20] = PID0+1;

            amr->patch[0][lv][ SibPID0+2 ]->sibling[24] = PID0+6;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[ 8] = PID0+6;
         }
      break;


      case  8:
         amr->patch[0][lv][ PID0+2 ]->sibling[ 8] = SibPID0+1;
         amr->patch[0][lv][ PID0+2 ]->sibling[24] = SibPID0+6;

         amr->patch[0][lv][ PID0+5 ]->sibling[20] = SibPID0+1;
         amr->patch[0][lv][ PID0+5 ]->sibling[ 8] = SibPID0+6;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+1 ]->sibling[ 7] = PID0+2;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[19] = PID0+2;

            amr->patch[0][lv][ SibPID0+1 ]->sibling[23] = PID0+5;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[ 7] = PID0+5;
         }
      break;


      case  9:
         amr->patch[0][lv][ PID0+4 ]->sibling[ 9] = SibPID0+0;
         amr->patch[0][lv][ PID0+4 ]->sibling[25] = SibPID0+3;

         amr->patch[0][lv][ PID0+7 ]->sibling[21] = SibPID0+0;
         amr->patch[0][lv][ PID0+7 ]->sibling[ 9] = SibPID0+3;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+0 ]->sibling[ 6] = PID0+4;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[18] = PID0+4;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[22] = PID0+7;
            amr->patch[0][lv][ SibPID0+3 ]->sibling[ 6] = PID0+7;
         }
      break;


      case 10:
         amr->patch[0][lv][ PID0+0 ]->sibling[10] = SibPID0+5;
         amr->patch[0][lv][ PID0+0 ]->sibling[19] = SibPID0+7;

         amr->patch[0][lv][ PID0+1 ]->sibling[18] = SibPID0+5;
         amr->patch[0][lv][ PID0+1 ]->sibling[10] = SibPID0+7;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+5 ]->sibling[13] = PID0+0;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[24] = PID0+0;

            amr->patch[0][lv][ SibPID0+5 ]->sibling[25] = PID0+1;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[13] = PID0+1;
         }
      break;


      case 11:
         amr->patch[0][lv][ PID0+2 ]->sibling[11] = SibPID0+3;
         amr->patch[0][lv][ PID0+2 ]->sibling[21] = SibPID0+6;

         amr->patch[0][lv][ PID0+4 ]->sibling[20] = SibPID0+3;
         amr->patch[0][lv][ PID0+4 ]->sibling[11] = SibPID0+6;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+3 ]->sibling[12] = PID0+2;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[22] = PID0+2;

            amr->patch[0][lv][ SibPID0+3 ]->sibling[23] = PID0+4;
            amr->patch[0][lv][ SibPID0+6 ]->sibling[12] = PID0+4;
         }
      break;


      case 12:
         amr->patch[0][lv][ PID0+3 ]->sibling[12] = SibPID0+2;
         amr->patch[0][lv][ PID0+3 ]->sibling[23] = SibPID0+4;

         amr->patch[0][lv][ PID0+6 ]->sibling[22] = SibPID0+2;
         amr->patch[0][lv][ PID0+6 ]->sibling[12] = SibPID0+4;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+2 ]->sibling[11] = PID0+3;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[20] = PID0+3;

            amr->patch[0][lv][ SibPID0+2 ]->sibling[21] = PID0+6;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[11] = PID0+6;
         }
      break;


      case 13:
         amr->patch[0][lv][ PID0+5 ]->sibling[13] = SibPID0+0;
         amr->patch[0][lv][ PID0+5 ]->sibling[25] = SibPID0+1;

         amr->patch[0][lv][ PID0+7 ]->sibling[24] = SibPID0+0;
         amr->patch[0][lv][ PID0+7 ]->sibling[13] = SibPID0+1;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+0 ]->sibling[10] = PID0+5;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[18] = PID0+5;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[19] = PID0+7;
            amr->patch[0][lv][ SibPID0+1 ]->sibling[10] = PID0+7;
         }
      break;


      case 14:
         amr->patch[0][lv][ PID0+0 ]->sibling[14] = SibPID0+6;
         amr->patch[0][lv][ PID0+0 ]->sibling[20] = SibPID0+7;

         amr->patch[0][lv][ PID0+2 ]->sibling[18] = SibPID0+6;
         amr->patch[0][lv][ PID0+2 ]->sibling[14] = SibPID0+7;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+6 ]->sibling[17] = PID0+0;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[23] = PID0+0;

            amr->patch[0][lv][ SibPID0+6 ]->sibling[25] = PID0+2;
            amr->patch[0][lv][ SibPID0+7 ]->sibling[17] = PID0+2;
         }
      break;


      case 15:
         amr->patch[0][lv][ PID0+3 ]->sibling[15] = SibPID0+1;
         amr->patch[0][lv][ PID0+3 ]->sibling[24] = SibPID0+4;

         amr->patch[0][lv][ PID0+5 ]->sibling[22] = SibPID0+1;
         amr->patch[0][lv][ PID0+5 ]->sibling[15] = SibPID0+4;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+1 ]->sibling[16] = PID0+3;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[19] = PID0+3;

            amr->patch[0][lv][ SibPID0+1 ]->sibling[21] = PID0+5;
            amr->patch[0][lv][ SibPID0+4 ]->sibling[16] = PID0+5;
         }
      break;


      case 16:
         amr->patch[0][lv][ PID0+1 ]->sibling[16] = SibPID0+3;
         amr->patch[0][lv][ PID0+1 ]->sibling[21] = SibPID0+5;

         amr->patch[0][lv][ PID0+4 ]->sibling[19] = SibPID0+3;
         amr->patch[0][lv][ PID0+4 ]->sibling[16] = SibPID0+5;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+3 ]->sibling[15] = PID0+1;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[22] = PID0+1;

            amr->patch[0][lv][ SibPID0+3 ]->sibling[24] = PID0+4;
            amr->patch[0][lv][ SibPID0+5 ]->sibling[15] = PID0+4;
         }
      break;


      case 17:
         amr->patch[0][lv][ PID0+6 ]->sibling[17] = SibPID0+0;
         amr->patch[0][lv][ PID0+6 ]->sibling[25] = SibPID0+2;

         amr->patch[0][lv][ PID0+7 ]->sibling[23] = SibPID0+0;
         amr->patch[0][lv][ PID0+7 ]->sibling[17] = SibPID0+2;

         if ( BothSide )
         {
            amr->patch[0][lv][ SibPID0+0 ]->sibling[14] = PID0+6;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[18] = PID0+6;

            amr->patch[0][lv][ SibPID0+0 ]->sibling[20] = PID0+7;
            amr->patch[0][lv][ SibPID0+2 ]->sibling[14] = PID0+7;
         }
      break;


      case 18:
         amr->patch[0][lv][ PID0+0 ]->sibling[18] = SibPID0+7;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+7 ]->sibling[25] = PID0+0;
      break;


      case 19:
         amr->patch[0][lv][ PID0+1 ]->sibling[19] = SibPID0+5;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+5 ]->sibling[24] = PID0+1;
      break;


      case 20:
         amr->patch[0][lv][ PID0+2 ]->sibling[20] = SibPID0+6;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+6 ]->sibling[23] = PID0+2;
      break;


      case 21:
         amr->patch[0][lv][ PID0+4 ]->sibling[21] = SibPID0+3;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+3 ]->sibling[22] = PID0+4;
      break;


      case 22:
         amr->patch[0][lv][ PID0+3 ]->sibling[22] = SibPID0+4;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+4 ]->sibling[21] = PID0+3;
      break;


      case 23:
         amr->patch[0][lv][ PID0+6 ]->sibling[23] = SibPID0+2;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+2 ]->sibling[20] = PID0+6;
      break;


      case 24:
         amr->patch[0][lv][ PID0+5 ]->sibling[24] = SibPID0+1;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+1 ]->sibling[19] = PID0+5;
      break;


      case 25:
         amr->patch[0][lv][ PID0+7 ]->sibling[25] = SibPID0+0;

         if ( BothSide )
            amr->patch[0][lv][ SibPID0+0 ]->sibling[18] = PID0+7;
      break;


      default:
         Aux_Error( ERROR_INFO, "lv %d, incorrect SibID (%d) !!\n", lv, SibID );
         exit(1);
   } // switch ( SibID )

} // FUNCTION : SetSiblingInDiffPatchGroup



//-------------------------------------------------------------------------------------------------------
// Function    :  SetSiblingExternal
// Description :  Set the sibling indices of patches adjacent to the simulation domain (for non-periodic B.C. only)
//
// Note        :  If a target sibling patch is an external patch (which lies outside the simulation domain),
//                the corresponding sibling index is set to "SIB_OFFSET_NONPERIODIC-Sibling", where Sibling
//                represents the sibling direction of the boundary region.
//
// Parameter   :  lv         : Target refinement level
//                NTarget0   : Number of target patches (with LocalID==0) in "TargetPID0"
//                TargetPID0 : Lists recording all target patches (with LocalID==0)
//-------------------------------------------------------------------------------------------------------
void SetSiblingExternal( const int lv, const int NTarget0, const int *TargetPID0 )
{

// check
#  ifdef GAMER_DEBUG
   if ( NTarget0 > 0  &&  TargetPID0 == NULL )
      Aux_Error( ERROR_INFO, "something is wrong ... (NTarget0 = %d) !!\n", NTarget0 );
#  endif


// maximum corner coordinates for the internal patches
   const int MaxIntCr[3] = { amr->BoxScale[0] - PATCH_SIZE*amr->scale[lv],
                             amr->BoxScale[1] - PATCH_SIZE*amr->scale[lv],
                             amr->BoxScale[2] - PATCH_SIZE*amr->scale[lv] };
   int *Sibling = NULL;
   int PID0, Cr[3];

#  pragma omp parallel for private( Cr, Sibling, PID0 ) schedule( runtime )
   for (int t=0; t<NTarget0; t++)
   {
      PID0 = TargetPID0[t];

//    no external patches should ever exist for the non-periodic B.C.
#     ifdef GAMER_DEBUG
      for (int d=0; d<3; d++)
      {
         if (  OPT__BC_FLU[2*d] != BC_FLU_PERIODIC  &&
               ( amr->patch[0][lv][PID0]->corner[d] < 0 || amr->patch[0][lv][PID0+7]->corner[d] > MaxIntCr[d] )  )
            Aux_Error( ERROR_INFO, "external patch (lv %d, PID0 %d, dir %d) is found for non-periodic B.C. !!\n",
                       lv, PID0, d );
      }
#     endif

//    skip all patch groups not adjacent to the simulation bondaries
      if ( amr->patch[0][lv][PID0]->corner[0] > 0  &&  amr->patch[0][lv][PID0+7]->corner[0] < MaxIntCr[0]  &&
           amr->patch[0][lv][PID0]->corner[1] > 0  &&  amr->patch[0][lv][PID0+7]->corner[1] < MaxIntCr[1]  &&
           amr->patch[0][lv][PID0]->corner[2] > 0  &&  amr->patch[0][lv][PID0+7]->corner[2] < MaxIntCr[2]     )
         continue;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         Sibling = amr->patch[0][lv][PID]->sibling;

//       apply periodicity to corner[] so that external and internal patches can be treated the same
         for (int d=0; d<3; d++)
         {
            Cr[d] = amr->patch[0][lv][PID]->corner[d];

            if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC )   Cr[d] = ( Cr[d] + amr->BoxScale[d] ) % amr->BoxScale[d];
         }

//       reset Sibling[] only if they are -1 to skip the periodic directions
//       sibling 0
         if ( Cr[0] == 0           )         SetSiblingExternal_CheckPeriodicity( Sibling+ 0,  0 );

//       sibling 1
         if ( Cr[0] == MaxIntCr[0] )         SetSiblingExternal_CheckPeriodicity( Sibling+ 1,  1 );

//       sibling 2
         if ( Cr[1] == 0           )         SetSiblingExternal_CheckPeriodicity( Sibling+ 2,  2 );

//       sibling 3
         if ( Cr[1] == MaxIntCr[1] )         SetSiblingExternal_CheckPeriodicity( Sibling+ 3,  3 );

//       sibling 4
         if ( Cr[2] == 0           )         SetSiblingExternal_CheckPeriodicity( Sibling+ 4,  4 );

//       sibling 5
         if ( Cr[2] == MaxIntCr[2] )         SetSiblingExternal_CheckPeriodicity( Sibling+ 5,  5 );

//       sibling 6
         if ( Cr[0] == 0 )
         {
            if ( Cr[1] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+ 6,  6 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+ 6,  0 );
         }
         else
            if ( Cr[1] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+ 6,  2 );

//       sibling 7
         if ( Cr[0] == MaxIntCr[0] )
         {
            if ( Cr[1] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+ 7,  7 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+ 7,  1 );
         }
         else
            if ( Cr[1] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+ 7,  2 );

//       sibling 8
         if ( Cr[0] == 0 )
         {
            if ( Cr[1] == MaxIntCr[1] )      SetSiblingExternal_CheckPeriodicity( Sibling+ 8,  8 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+ 8,  0 );
         }
         else
            if ( Cr[1] == MaxIntCr[1] )      SetSiblingExternal_CheckPeriodicity( Sibling+ 8,  3 );

//       sibling 9
         if ( Cr[0] == MaxIntCr[0] )
         {
            if ( Cr[1] == MaxIntCr[1] )      SetSiblingExternal_CheckPeriodicity( Sibling+ 9,  9 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+ 9,  1 );
         }
         else
            if ( Cr[1] == MaxIntCr[1] )      SetSiblingExternal_CheckPeriodicity( Sibling+ 9,  3 );

//       sibling 10
         if ( Cr[1] == 0 )
         {
            if ( Cr[2] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+10, 10 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+10,  2 );
         }
         else
            if ( Cr[2] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+10,  4 );

//       sibling 11
         if ( Cr[1] == MaxIntCr[1] )
         {
            if ( Cr[2] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+11, 11 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+11,  3 );
         }
         else
            if ( Cr[2] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+11,  4 );

//       sibling 12
         if ( Cr[1] == 0 )
         {
            if ( Cr[2] == MaxIntCr[2] )      SetSiblingExternal_CheckPeriodicity( Sibling+12, 12 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+12,  2 );
         }
         else
            if ( Cr[2] == MaxIntCr[2] )      SetSiblingExternal_CheckPeriodicity( Sibling+12,  5 );

//       sibling 13
         if ( Cr[1] == MaxIntCr[1] )
         {
            if ( Cr[2] == MaxIntCr[2] )      SetSiblingExternal_CheckPeriodicity( Sibling+13, 13 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+13,  3 );
         }
         else
            if ( Cr[2] == MaxIntCr[2] )      SetSiblingExternal_CheckPeriodicity( Sibling+13,  5 );

//       sibling 14
         if ( Cr[2] == 0 )
         {
            if ( Cr[0] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+14, 14 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+14,  4 );
         }
         else
            if ( Cr[0] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+14,  0 );

//       sibling 15
         if ( Cr[2] == MaxIntCr[2] )
         {
            if ( Cr[0] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+15, 15 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+15,  5 );
         }
         else
            if ( Cr[0] == 0 )                SetSiblingExternal_CheckPeriodicity( Sibling+15,  0 );

//       sibling 16
         if ( Cr[2] == 0 )
         {
            if ( Cr[0] == MaxIntCr[0] )      SetSiblingExternal_CheckPeriodicity( Sibling+16, 16 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+16,  4 );
         }
         else
            if ( Cr[0] == MaxIntCr[0] )      SetSiblingExternal_CheckPeriodicity( Sibling+16,  1 );

//       sibling 17
         if ( Cr[2] == MaxIntCr[2] )
         {
            if ( Cr[0] == MaxIntCr[0] )      SetSiblingExternal_CheckPeriodicity( Sibling+17, 17 );
            else                             SetSiblingExternal_CheckPeriodicity( Sibling+17,  5 );
         }
         else
            if ( Cr[0] == MaxIntCr[0] )      SetSiblingExternal_CheckPeriodicity( Sibling+17,  1 );

//       sibling 18
         if ( Cr[0] == 0 )
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+18, 18 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+18,  6 );
            }
            else
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+18, 14 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+18,  0 );
            }
         }
         else
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+18, 10 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+18,  2 );
            }
            else
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+18,  4 );
         }

//       sibling 19
         if ( Cr[0] == MaxIntCr[0] )
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+19, 19 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+19,  7 );
            }
            else
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+19, 16 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+19,  1 );
            }
         }
         else
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+19, 10 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+19,  2 );
            }
            else
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+19,  4 );
         }

//       sibling 20
         if ( Cr[0] == 0 )
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+20, 20 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+20,  8 );
            }
            else
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+20, 14 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+20,  0 );
            }
         }
         else
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+20, 11 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+20,  3 );
            }
            else
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+20,  4 );
         }

//       sibling 21
         if ( Cr[0] == MaxIntCr[0] )
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+21, 21 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+21,  9 );
            }
            else
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+21, 16 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+21,  1 );
            }
         }
         else
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+21, 11 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+21,  3 );
            }
            else
               if ( Cr[2] == 0 )             SetSiblingExternal_CheckPeriodicity( Sibling+21,  4 );
         }

//       sibling 22
         if ( Cr[0] == 0 )
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+22, 22 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+22,  6 );
            }
            else
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+22, 15 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+22,  0 );
            }
         }
         else
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+22, 12 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+22,  2 );
            }
            else
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+22,  5 );
         }

//       sibling 23
         if ( Cr[0] == MaxIntCr[0] )
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+23, 23 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+23,  7 );
            }
            else
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+23, 17 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+23,  1 );
            }
         }
         else
         {
            if ( Cr[1] == 0 )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+23, 12 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+23,  2 );
            }
            else
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+23,  5 );
         }

//       sibling 24
         if ( Cr[0] == 0 )
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+24, 24 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+24,  8 );
            }
            else
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+24, 15 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+24,  0 );
            }
         }
         else
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+24, 13 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+24,  3 );
            }
            else
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+24,  5 );
         }

//       sibling 25
         if ( Cr[0] == MaxIntCr[0] )
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+25, 25 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+25,  9 );
            }
            else
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+25, 17 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+25,  1 );
            }
         }
         else
         {
            if ( Cr[1] == MaxIntCr[1] )
            {
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+25, 13 );
               else                          SetSiblingExternal_CheckPeriodicity( Sibling+25,  3 );
            }
            else
               if ( Cr[2] == MaxIntCr[2] )   SetSiblingExternal_CheckPeriodicity( Sibling+25,  5 );
         }

      } // for (int PID=PID0; PID<PID0+8; PID++)
   } // for (int t=0; t<NTarget0; t++)

} // FUNCTION : SetSiblingExternal



//-------------------------------------------------------------------------------------------------------
// Function    :  SetSiblingExternal_CheckPeriodicity
// Description :  Set the sibling indices of patches adjacent to the simulation domain if the target
//                boundary direction is non-periodic
//
// Note        :  1. Invoked by SetSiblingExternal()
//                2. Do nothing if the target boundary direction is periodic
//                   --> A boundary direction is defined to be periodic only if its associated directions
//                       along x/y/z are all periodic
//                   --> For example, Boundary=26 is periodic only if x, y, and z directions are all periodic
//
// Parameter   :  Sibling  : Pointer of Sibling[] to be set
//                Boundary : Boundary direction (0-25)
//-------------------------------------------------------------------------------------------------------
void SetSiblingExternal_CheckPeriodicity( int *Sibling, const int Boundary )
{

   bool Periodic = true;

   for (int d=0; d<3; d++)
   {
      if (  TABLE_01( Boundary, 'x'+d, true, false, true )  &&  OPT__BC_FLU[2*d] != BC_FLU_PERIODIC  )
      {
         Periodic = false;
         break;
      }
   }

   if ( !Periodic )  *Sibling = SIB_OFFSET_NONPERIODIC - Boundary;

} // FUNCTION : SetSiblingExternal_CheckPeriodicity



#endif // #ifdef LOAD_BALANCE
