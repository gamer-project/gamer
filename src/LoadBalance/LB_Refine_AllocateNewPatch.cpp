#include "Copyright.h"
#include "GAMER.h"

#ifdef LOAD_BALANCE



void PrepareCData( const int FaLv, const int FaPID, real *const FaData,
                   const int FaSg_Flu, const int FaGhost_Flu, const int NSide_Flu,
                   const int FaSg_Pot, const int FaGhost_Pot, const int NSide_Pot,
                   const int BC_Face[], const int FluVarIdxList[] );
void LB_Refine_AllocateBufferPatch_Sibling( const int SonLv );
static int AllocateSonPatch( const int FaLv, const int *Cr, const int PScale, const int FaPID, real *CData,
                             const int CGhost_Flu, const int NSide_Flu, const int CGhost_Pot, const int NSide_Pot,
                             const int BC_Face[], const int FluVarIdxList[] );
static void DeallocateSonPatch( const int FaLv, const int FaPID, const int NNew_Real0, int NewSonPID0_Real[], 
                                const int SwitchIdx );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine_AllocateNewPatch
// Description :  Allocate/deallocate real son patches at FaLv+1
//
// Note        :  1. This function is invoked by the function "LB_Refine"
//                2. Home/Away : targeted patches at home/not at home
//                3. Input Cr1D and CData lists are unsorted
//                4. After invoking this function, some father-buffer patches at FaLv may become useless
//                   --> Currently we do not remove these patches for the consideration of better performance
//                5. This function will also allocate new father-buffer patches at FaLv, and new sibling-buffer
//                   and father-buffer patches at FaLv+1
//                6. All MPI lists are NOT reconstructed here
//                7. Several alternative functions are invoked here for better performance
//                   (e.g., LB_AllocateBufferPatch_Sibling --> LB_Refine_AllocateBufferPatch_Sibling)
//                8. All buffer patches at FaLv+1 are deallocated before allocating/deallocating real
//                   patches at FaLv+1
//                   --> make implementation easier
//
// Parameter   :  FaLv           : Targeted refinement level to be refined
//                NNew_Home      : Number of home patches at FaLv to allocate son patches
//                NewPID_Home    : Patch indices of home patches at FaLv to allocate son patches
//                NNew_Away      : Number of away patches at FaLv to allocate son patches
//                NewCr1D_Away   : Padded 1D corner of away patches at FaLv to allocate son patches
//                NewCData_Away  : Coarse-grid data of away patches at FaLv to allocate son patches
//                NDel_Home      : Number of home patches at FaLv to deallocate son patches
//                DelPID_Home    : Patch indices of home patches at FaLv to deallocate son patches
//                NDel_Away      : Number of away patches at FaLv to deallocate son patches
//                DelCr1D_Away   : Padded 1D corner of away patches at FaLv to deallocate son patches
//-------------------------------------------------------------------------------------------------------
void LB_Refine_AllocateNewPatch( const int FaLv, int NNew_Home, int *NewPID_Home, int NNew_Away, 
                                 ulong *NewCr1D_Away, real *NewCData_Away, int NDel_Home, int *DelPID_Home, 
                                 int NDel_Away, ulong *DelCr1D_Away )
{

   const int SonLv    = FaLv + 1;
   const int GraLv    = FaLv + 2;
   const int SonNReal = amr->NPatchComma[SonLv][1];
   const int SonNBuff = amr->NPatchComma[SonLv][3] - SonNReal;
   const int FaNPatch = amr->num[FaLv];


// 1. sort the away patches and get the matching lists
// ==========================================================================================
   int *NewCr1D_Away_IdxTable = new int [NNew_Away];
   int *Match_New             = new int [NNew_Away];
   int *Match_Del             = new int [NDel_Away];
   int *DelPID_Away           = new int [NDel_Away];

   Mis_Heapsort( NNew_Away, NewCr1D_Away, NewCr1D_Away_IdxTable );
   Mis_Heapsort( NDel_Away, DelCr1D_Away, NULL                   );

   Mis_Matching_int( FaNPatch, amr->LB->PaddedCr1DList[FaLv], NNew_Away, NewCr1D_Away, Match_New );
   Mis_Matching_int( FaNPatch, amr->LB->PaddedCr1DList[FaLv], NDel_Away, DelCr1D_Away, Match_Del );

   for (int t=0; t<NDel_Away; t++)  
   {
#     ifdef GAMER_DEBUG
      if ( Match_Del[t] == -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, away patch with Cr1D %lu found no matching !!\n", 
                    FaLv, DelCr1D_Away[t] );
#     endif

      DelPID_Away[t] = amr->LB->PaddedCr1DList_IdxTable[FaLv][ Match_Del[t] ];
   }



// 2. remove all buffer patches at SonLv and backup the fluid and pot arrays
// ==========================================================================================
   const int MirSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const int FSg_Flu    = amr->FluSg[SonLv];
#  ifdef GRAVITY
   const int FSg_Pot    = amr->PotSg[SonLv];
#  endif

   int SibPID, NBufBk=0;   // BufBk : backup the data of buffer patches

   ulong *PCr1D_BufBk                  = new ulong [SonNBuff];
   int   *PCr1D_BufBk_IdxTable         = new int   [SonNBuff];
   real (**fluid_BufBk)[PS1][PS1][PS1] = new ( real (*[SonNBuff])[PS1][PS1][PS1] );
#  ifdef GRAVITY
   real (**pot_BufBk)  [PS1][PS1]      = new ( real (*[SonNBuff])[PS1][PS1] );
#  endif

   if ( SonNBuff != 0 )
   {
      for (int SonPID=SonNReal; SonPID<amr->NPatchComma[SonLv][3]; SonPID++)
      {
//       backup fluid and pot arrays
#        ifdef GRAVITY
         if ( amr->patch[0][SonLv][SonPID]->fluid != NULL  ||  amr->patch[0][SonLv][SonPID]->pot != NULL )
#        else
         if ( amr->patch[0][SonLv][SonPID]->fluid != NULL )
#        endif
         {
            fluid_BufBk[NBufBk] = amr->patch[FSg_Flu][SonLv][SonPID]->fluid;
            amr->patch[FSg_Flu][SonLv][SonPID]->fluid = NULL;

#           ifdef GRAVITY
            pot_BufBk  [NBufBk] = amr->patch[FSg_Pot][SonLv][SonPID]->pot;
            amr->patch[FSg_Pot][SonLv][SonPID]->pot   = NULL;
#           endif

            PCr1D_BufBk[NBufBk] = amr->patch[0][SonLv][SonPID]->PaddedCr1D;
            NBufBk ++;
         }

//       reset sibling indices to -1
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][SonLv][SonPID]->sibling[s];
            if ( SibPID >= 0 )   amr->patch[0][SonLv][SibPID]->sibling[ MirSib[s] ] = -1;
         }

//       deallocate patches
         amr->patch[0][SonLv][SonPID]->son = -1;
         amr->pdelete( SonLv, SonPID );

      } // for (int SonPID=SonNReal; SonPID<amr->NPatchComma[SonLv][3]; SonPID++)

//    sort PCr1D_BufBk
      Mis_Heapsort( NBufBk, PCr1D_BufBk, PCr1D_BufBk_IdxTable );

//    reset NPatchComma
      for (int m=2; m<28; m++)   amr->NPatchComma[SonLv][m] = SonNReal;

//    reset calculate PaddedCr1DList
      amr->LB->PaddedCr1DList         [SonLv] = (ulong*)realloc( amr->LB->PaddedCr1DList         [SonLv], 
                                                                 SonNReal*sizeof(ulong) );
      amr->LB->PaddedCr1DList_IdxTable[SonLv] = (int*  )realloc( amr->LB->PaddedCr1DList_IdxTable[SonLv], 
                                                                 SonNReal*sizeof(int ) );

      for (int SonPID=0; SonPID<SonNReal; SonPID++)
         amr->LB->PaddedCr1DList[SonLv][SonPID] = amr->patch[0][SonLv][SonPID]->PaddedCr1D;

      Mis_Heapsort( SonNReal, amr->LB->PaddedCr1DList[SonLv], amr->LB->PaddedCr1DList_IdxTable[SonLv] );

   } // if ( SonNBuff != 0 )



// 3. allocate new real patches at SonLv
// ==========================================================================================
   const int Padded              = 1<<NLEVEL;
   const int BoxNScale_Padded[3] = { amr->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                     amr->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                     amr->BoxScale[2]/PATCH_SIZE + 2*Padded }; //normalized and padded BoxScale
   const int PScale              = PATCH_SIZE*amr->scale[SonLv];  // scale of a single patch at SonLv
   const int NNew_Real0          = NNew_Home + NNew_Away;

   int FaPID, SonPID0, Cr3D[3], *Cr3D_Ptr=NULL, NNoFa=0;
   int *NewSonPID0_NoFa = new int [ NNew_Away ];         // NNew_Away is the maximum number this array can have
   int *NewSonPID0_All  = (int*)malloc( NNew_Real0*sizeof(int) );
   int *NewSonPID0_Real = NewSonPID0_All;
   int *NewSonPID0_Away = NewSonPID0_All + NNew_Home;

// parameters for spatial interpolation
   int NSide_Flu, CGhost_Flu;

   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, CGhost_Flu );

   const int CSize_Flu = PATCH_SIZE + 2*CGhost_Flu;
#  ifdef GRAVITY
   int NSide_Pot, CGhost_Pot;

   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, CGhost_Pot );

   const int CSize_Pot = PATCH_SIZE + 2*CGhost_Pot;
   const int CSize_Tot = NCOMP*CSize_Flu*CSize_Flu*CSize_Flu + CSize_Pot*CSize_Pot*CSize_Pot;
#  else
   const int CSize_Tot = NCOMP*CSize_Flu*CSize_Flu*CSize_Flu;
#  endif

// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   int BC_Face[26], BC_Face_tmp[3], FluVarIdxList[NCOMP];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1 )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1 )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1 )   BC_Face[s] = BC_Face_tmp[0];
   }

   for (int v=0; v<NCOMP; v++)   FluVarIdxList[v] = v;


// 3.1 home patches
   for (int t=0; t<NNew_Home; t++)
   {
      FaPID    = NewPID_Home[t];
      Cr3D_Ptr = amr->patch[0][FaLv][FaPID]->corner;

#     ifdef GRAVITY
      NewSonPID0_All[t] = AllocateSonPatch( FaLv, Cr3D_Ptr, PScale, FaPID, NULL, CGhost_Flu, NSide_Flu, CGhost_Pot, NSide_Pot,
                                            BC_Face, FluVarIdxList );
#     else
      NewSonPID0_All[t] = AllocateSonPatch( FaLv, Cr3D_Ptr, PScale, FaPID, NULL, CGhost_Flu, NSide_Flu, NULL_INT, NULL_INT,
                                            BC_Face, FluVarIdxList );
#     endif
   }


// 3.2 away patches
   for (int t=0; t<NNew_Away; t++)
   {
//    3.2.1 away patches without father patch
      if ( Match_New[t] == -1 )
      {
         FaPID = -1;
         Mis_Idx1D2Idx3D( BoxNScale_Padded, NewCr1D_Away[t], Cr3D );
         for (int d=0; d<3; d++)    Cr3D[d] = ( Cr3D[d] - Padded )*PATCH_SIZE;

#        ifdef GRAVITY
         NewSonPID0_Away[t] = AllocateSonPatch( FaLv, Cr3D, PScale, FaPID, 
                                                NewCData_Away+NewCr1D_Away_IdxTable[t]*CSize_Tot,
                                                CGhost_Flu, NSide_Flu, CGhost_Pot, NSide_Pot,
                                                NULL, NULL );
#        else
         NewSonPID0_Away[t] = AllocateSonPatch( FaLv, Cr3D, PScale, FaPID, 
                                                NewCData_Away+NewCr1D_Away_IdxTable[t]*CSize_Tot,
                                                CGhost_Flu, NSide_Flu, NULL_INT, NULL_INT,
                                                NULL, NULL );
#        endif

//       record the SonPID (with LocalID == 0 ) with no father at home
#        ifdef GAMER_DEBUG
         if ( NNoFa >= NNew_Away )
            Aux_Error( ERROR_INFO, "FaLv %d, NNoFa (%d) exceeds the maximum number (%d) !!\n",
                       FaLv, NNoFa, NNew_Away );
#        endif

         NewSonPID0_NoFa[ NNoFa ++ ] = NewSonPID0_Away[t];
      }

//    3.2.1 away patches with father patch
      else
      {
         FaPID    = amr->LB->PaddedCr1DList_IdxTable[FaLv][ Match_New[t] ];
         Cr3D_Ptr = amr->patch[0][FaLv][FaPID]->corner;

#        ifdef GRAVITY
         NewSonPID0_Away[t] = AllocateSonPatch( FaLv, Cr3D_Ptr, PScale, FaPID, 
                                                NewCData_Away+NewCr1D_Away_IdxTable[t]*CSize_Tot,
                                                CGhost_Flu, NSide_Flu, CGhost_Pot, NSide_Pot,
                                                NULL, NULL );
#        else
         NewSonPID0_Away[t] = AllocateSonPatch( FaLv, Cr3D_Ptr, PScale, FaPID, 
                                                NewCData_Away+NewCr1D_Away_IdxTable[t]*CSize_Tot,
                                                CGhost_Flu, NSide_Flu, NULL_INT, NULL_INT,
                                                NULL, NULL );
#        endif
      } // if ( Match_New[t] == -1 ) ... else ...
   } // for (int t=0; t<NNew_Away; t++)



// 4. allocate new father-buffer patches at FaLv and construct the relation son->father
// ==========================================================================================
   int  NNewFaBuf0;
   int *NewFaBufPID0 = NULL;

   LB_AllocateBufferPatch_Father( SonLv, false, NNew_Real0, NewSonPID0_Real, true, &NNewFaBuf0, &NewFaBufPID0 );

   LB_FindFather( SonLv, false, NNoFa, NewSonPID0_NoFa );



// 5. deallocate unflagged real patches at SonLv
// ==========================================================================================
   int SwitchIdx = NNew_Real0 - 1;  // Element in NewSonPID0_Real to be reset

// 5.1 home patches
   for (int t=0; t<NDel_Home; t++)
   {
      FaPID = DelPID_Home[t];

      DeallocateSonPatch( FaLv, FaPID, NNew_Real0, NewSonPID0_Real, SwitchIdx-- );
   } 

// 5.2 away patches
   for (int t=0; t<NDel_Away; t++)
   {
      FaPID = DelPID_Away[t];

      DeallocateSonPatch( FaLv, FaPID, NNew_Real0, NewSonPID0_Real, SwitchIdx-- );
   }



// 6. record NPatchComma, NPatchTotal and LB_IdxList_Real at SonLv
// ==========================================================================================
   const int SonNReal_New = amr->NPatchComma[SonLv][1];

// 6.1 record the number of real patches at SonLv for each rank
   if ( SonNReal_New != amr->num[SonLv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d)\" !!\n",
                 SonLv, amr->NPatchComma[SonLv][1], SonLv, amr->num[SonLv] );

   for (int m=2; m<28; m++)   amr->NPatchComma[SonLv][m] = SonNReal_New;


// 6.2 record the total number of real patches at SonLv among all ranks
   Mis_GetTotalPatchNumber( SonLv );


// 6.3 record LB_IdxList_Real at SonLv
   if ( amr->LB->IdxList_Real         [SonLv] != NULL )  delete [] amr->LB->IdxList_Real         [SonLv];
   if ( amr->LB->IdxList_Real_IdxTable[SonLv] != NULL )  delete [] amr->LB->IdxList_Real_IdxTable[SonLv];

   amr->LB->IdxList_Real         [SonLv] = new long [SonNReal_New];
   amr->LB->IdxList_Real_IdxTable[SonLv] = new int  [SonNReal_New];

   for (int SonPID=0; SonPID<SonNReal_New; SonPID++)   
      amr->LB->IdxList_Real[SonLv][SonPID] = amr->patch[0][SonLv][SonPID]->LB_Idx;

   Mis_Heapsort( SonNReal_New, amr->LB->IdxList_Real[SonLv], amr->LB->IdxList_Real_IdxTable[SonLv] );

   
// 6.4 check : no duplicate patches at FaLv and SonLv
#  ifdef GAMER_DEBUG
   for (int lv=FaLv; lv<=SonLv; lv++)
   {
      ulong *TempPaddedCr1D          = new ulong [ amr->num[lv] ];
      int   *TempPaddedCr1D_IdxTable = new  int  [ amr->num[lv] ];

      for (int PID=0; PID<amr->num[lv]; PID++)
         TempPaddedCr1D[PID] = amr->patch[0][lv][PID]->PaddedCr1D;

      Mis_Heapsort( amr->num[lv], TempPaddedCr1D, TempPaddedCr1D_IdxTable );

      for (int t=1; t<amr->num[lv]; t++)
      {
         if ( TempPaddedCr1D[t] == TempPaddedCr1D[t-1] )
            Aux_Error( ERROR_INFO, "duplicate patches at lv %d, PaddedCr1D %lu, PID = %d and %d !!\n",
                       lv, TempPaddedCr1D[t], TempPaddedCr1D_IdxTable[t], TempPaddedCr1D_IdxTable[t-1] );
      }

      delete [] TempPaddedCr1D;
      delete [] TempPaddedCr1D_IdxTable;
   }
#  endif



// 7. construct new father->son relation between SonLv and FaLv and new sibling relation at FaLv
// ==========================================================================================
   LB_FindSonNotHome( FaLv, true, NULL_INT, NULL );

   LB_SiblingSearch( FaLv, false, NNewFaBuf0, NewFaBufPID0 );



// 8. allocate buffer patches at SonLv
// ==========================================================================================
// 8.1 sibling-buffer patches (invoke the alternatvie version which should be faster)
   LB_Refine_AllocateBufferPatch_Sibling( SonLv );
// LB_AllocateBufferPatch_Sibling( SonLv );

// 8.2 father-buffer patches
   if ( GraLv < NLEVEL )
   LB_AllocateBufferPatch_Father( GraLv, true, NULL_INT, NULL, false, NULL, NULL );



// 9. construct new father<->son relation between GraLv and SonLv and new sibling relation at SonLv
// ==========================================================================================
   const int NNew_Buf0 = ( amr->NPatchComma[SonLv][3] - amr->NPatchComma[SonLv][1] ) / 8;
   const int NNew_All0 = NNew_Real0 + NNew_Buf0;
   int *NewSonPID0_Buf, *NewSonPID_All, Count;

// 9.1 record the indices of all newly-allocated buffer patches
   NewSonPID_All  = new int [ 8*NNew_All0 ];
   NewSonPID0_All = (int*)realloc( NewSonPID0_All, NNew_All0*sizeof(int) );
   NewSonPID0_Buf = NewSonPID0_All + NNew_Real0;

   for (SonPID0=amr->NPatchComma[SonLv][1], Count=0; SonPID0<amr->NPatchComma[SonLv][3]; SonPID0+=8, Count++)
      NewSonPID0_Buf[Count] = SonPID0;

   for (int t=0, m=0; t<NNew_All0; t++, m+=8)  
   {
      SonPID0 = NewSonPID0_All[t];

      for (int SonPID=SonPID0, n=m; SonPID<SonPID0+8; SonPID++, n++)    NewSonPID_All[n] = SonPID;
   }

// 9.2 sibling relation
   LB_SiblingSearch( SonLv, false, NNew_All0, NewSonPID0_All );

// 9.3 father <-> son relation
   if ( GraLv < NLEVEL )
   {
      LB_FindFather( GraLv, true, NULL_INT, NULL );

      LB_FindSonNotHome( SonLv, false, 8*NNew_All0, NewSonPID_All );
   }



// 10. restore the fluid and pot arrays in the buffer patches
// ==========================================================================================
   const int FSg_Flu2 = 1 - FSg_Flu;
   real (*fluid_ptr)[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE] = NULL;

#  ifdef GRAVITY
   const int FSg_Pot2 = 1 - FSg_Pot;
   real (*pot_ptr)[PATCH_SIZE][PATCH_SIZE] = NULL;
#  endif

   int *Match_BufBk = new int [NBufBk];
   int  MPID;

// 10.1 get the match lists
   Mis_Matching_int( amr->num[SonLv], amr->LB->PaddedCr1DList[SonLv], NBufBk, PCr1D_BufBk, Match_BufBk );

// 10.2 reset array pointers
   for (int t=0; t<NBufBk; t++)
   {
      fluid_ptr = fluid_BufBk[ PCr1D_BufBk_IdxTable[t] ];
#     ifdef GRAVITY
      pot_ptr   = pot_BufBk  [ PCr1D_BufBk_IdxTable[t] ];
#     endif

      if ( Match_BufBk[t] != -1 )
      {
         MPID = amr->LB->PaddedCr1DList_IdxTable[SonLv][ Match_BufBk[t] ];

#        ifdef GAMER_DEBUG
         if ( MPID < amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "Match_PID = %d matches to a real patch (Match[%d] = %d, SonNReal = %d) !!\n",
                       MPID, t, Match_BufBk[t], amr->NPatchComma[SonLv][1] );
#        endif

         if ( fluid_ptr != NULL )
         {
            amr->patch[FSg_Flu ][SonLv][MPID]->fluid   = fluid_ptr;
#           if ( NPASSIVE > 0 )
            amr->patch[FSg_Flu ][SonLv][MPID]->passive = fluid_ptr + NCOMP;
#           endif
            amr->patch[FSg_Flu2][SonLv][MPID]->hnew();
         }

#        ifdef GRAVITY
         if ( pot_ptr != NULL )
         {
            amr->patch[FSg_Pot ][SonLv][MPID]->pot     = pot_ptr;
            amr->patch[FSg_Pot2][SonLv][MPID]->gnew();
         }
#        endif // #ifdef GARVITY

      } // if ( Match_BufBk[t] != -1 )

      else
      {
         if ( fluid_ptr != NULL )   delete [] fluid_ptr;
#        ifdef GRAVITY
         if ( pot_ptr   != NULL )   delete [] pot_ptr;
#        endif
      } // if ( Match_BufBk[t] != -1 ) ... else ...
   } // for (int t=0; t<NBufBk; t++)



// free memory
   free( NewSonPID0_All );
   delete [] NewCr1D_Away_IdxTable;
   delete [] Match_New;
   delete [] Match_Del;
   delete [] Match_BufBk;
   delete [] DelPID_Away;
   delete [] NewSonPID0_NoFa;
   delete [] NewSonPID_All;
   if ( NewFaBufPID0 != NULL )   delete [] NewFaBufPID0;
   delete [] PCr1D_BufBk;
   delete [] PCr1D_BufBk_IdxTable;
   delete [] fluid_BufBk;
#  ifdef GRAVITY
   delete [] pot_BufBk;
#  endif

} // FUNCTION : LB_Refine_AllocateNewPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  AllocateSonPatch  
// Description :  Allocate eight son patches at FaLv+1
//
// Note        :  Just to avoid duplicate code segment 
//
// Parameter   :  FaLv           : Targeted refinement level to be refined
//                Cr             : Corner coordinates of the son patch with LocalID == 0
//                PScale         : Scale of one patch at SonLv
//                FaPID          : Father patch index (can be -1 for the away patches)
//                CData          : Coarse-grid data for assigning data to son patches by spatial interpolation
//                                 (initialize as NULL if father patch is home --> prepare CData here)
//                BC_Face        : Corresponding boundary faces (0~5) along 26 sibling directions ->for non-periodic B.C. only
//                FluVarIdxList  : List of target fluid variable indices                          ->for non-periodic B.C. only
//
// Return      :  SonPID with LocalID == 0 
//-------------------------------------------------------------------------------------------------------
int AllocateSonPatch( const int FaLv, const int *Cr, const int PScale, const int FaPID, real *CData,
                      const int CGhost_Flu, const int NSide_Flu, const int CGhost_Pot, const int NSide_Pot,
                      const int BC_Face[], const int FluVarIdxList[] )
{

   const int SonLv   = FaLv + 1;
   const int SonPID0 = amr->num[SonLv];
   bool FaIsHome = false;

// 0. check : targeted father patch has no son
#  ifdef GAMER_DEBUG
   if ( FaPID != -1  &&  amr->patch[0][FaLv][FaPID]->son != -1 )
      Aux_Error( ERROR_INFO, "FaLv %d, FaPID (%d) already has sons (SonPID = %d, duplicate SonPID = %d) !!\n",
                 FaLv, FaPID, amr->patch[0][FaLv][FaPID]->son, SonPID0 );
#  endif


// 1. construct relation : father -> child
   if ( FaPID != -1 )   amr->patch[0][FaLv][FaPID]->son = SonPID0;


// 2. allocate child patches and construct relation : child -> father
   amr->pnew( SonLv, Cr[0],        Cr[1],        Cr[2],        FaPID, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1],        Cr[2],        FaPID, true, true );
   amr->pnew( SonLv, Cr[0],        Cr[1]+PScale, Cr[2],        FaPID, true, true );
   amr->pnew( SonLv, Cr[0],        Cr[1],        Cr[2]+PScale, FaPID, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1]+PScale, Cr[2],        FaPID, true, true );
   amr->pnew( SonLv, Cr[0],        Cr[1]+PScale, Cr[2]+PScale, FaPID, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1],        Cr[2]+PScale, FaPID, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1]+PScale, Cr[2]+PScale, FaPID, true, true );

   amr->NPatchComma[SonLv][1] += 8;


// 3. assign data to child patches by spatial interpolation
   const int CRange[3]     = { PS1, PS1, PS1 };
   const int FSize         = PS2;
   const int FStart[3]     = { 0, 0, 0 };

// fluid
   const int CSg_Flu       = amr->FluSg[ FaLv];
   const int FSg_Flu       = amr->FluSg[SonLv];
   const int CSize_Flu     = PS1 + 2*CGhost_Flu;
   const int CStart_Flu[3] = { CGhost_Flu, CGhost_Flu, CGhost_Flu }; 
   real (*FData_Flu)[FSize][FSize][FSize] = new real [NCOMP][FSize][FSize][FSize];

// potential
#  ifdef GRAVITY
   const int CSg_Pot       = amr->PotSg[ FaLv];
   const int FSg_Pot       = amr->PotSg[SonLv];
   const int CSize_Pot     = PS1 + 2*CGhost_Pot;
   const int CStart_Pot[3] = { CGhost_Pot, CGhost_Pot, CGhost_Pot }; 
   real (*FData_Pot)[FSize][FSize] = new real [FSize][FSize][FSize];

   const int CSize_Tot = NCOMP*CSize_Flu*CSize_Flu*CSize_Flu + CSize_Pot*CSize_Pot*CSize_Pot;
#  else
   const int CSize_Tot = NCOMP*CSize_Flu*CSize_Flu*CSize_Flu;
#  endif


// 3.1 prepare the coarse-grid data
   if ( CData == NULL )
   {
#     ifdef GAMER_DEBUG
      if ( FaPID >= amr->NPatchComma[FaLv][1]  ||  FaPID == -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID (%d) is NOT home (FaNReal = %d) !!\n",
                    FaLv, FaPID, amr->NPatchComma[FaLv][1] );
#     endif

      FaIsHome = true;
      CData    = new real [CSize_Tot];

#     ifdef GRAVITY
      PrepareCData( FaLv, FaPID, CData, CSg_Flu, CGhost_Flu, NSide_Flu, CSg_Pot, CGhost_Pot, NSide_Pot, 
                    BC_Face, FluVarIdxList );
#     else
      PrepareCData( FaLv, FaPID, CData, CSg_Flu, CGhost_Flu, NSide_Flu, NULL_INT, NULL_INT, NULL_INT,
                    BC_Face, FluVarIdxList );
#     endif
   }


// 3.2 perform spatial interpolation
   const int   CSize_Flu_Temp[3]      = { CSize_Flu, CSize_Flu, CSize_Flu };
   const int   FSize_Temp    [3]      = { PS2, PS2, PS2 };
   const int   CSize_Flu1v            = CSize_Flu*CSize_Flu*CSize_Flu;
   const bool  PhaseUnwrapping_Yes    = true;
   const bool  PhaseUnwrapping_No     = false;
   const bool  EnsureMonotonicity_Yes = true;
   const bool  EnsureMonotonicity_No  = false;
   real *const CData_Flu              = CData;
#  if ( MODEL == ELBDM )
   real *const CData_Dens             = CData_Flu + DENS*CSize_Flu1v;
   real *const CData_Real             = CData_Flu + REAL*CSize_Flu1v;
   real *const CData_Imag             = CData_Flu + IMAG*CSize_Flu1v;
#  endif

// 3.2.1 determine the variables which must be positive
   bool Monotonicity[NCOMP];

   for (int v=0; v<NCOMP; v++)
   {
#     if ( MODEL == HYDRO )
#     if ( NPASSIVE > 0 )
      if ( v == DENS  ||  v == ENGY  ||  v >= NCOMP )  
#     else
      if ( v == DENS  ||  v == ENGY )
#     endif
                                       Monotonicity[v] = EnsureMonotonicity_Yes;
      else                             Monotonicity[v] = EnsureMonotonicity_No;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      if ( v == DENS )                 Monotonicity[v] = EnsureMonotonicity_Yes;
      else                             Monotonicity[v] = EnsureMonotonicity_No;

#     else
#     warning : WARNING : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION ??          
#     endif // MODEL
   }

// 3.2.2 interpolation
#  if ( MODEL == ELBDM )
   if ( OPT__INT_PHASE )
   {
//    get the wrapped phase (store in the REAL component)
      for (int t=0; t<CSize_Flu1v; t++)   CData_Real[t] = ATAN2( CData_Imag[t], CData_Real[t] );

//    interpolate density 
      Interpolate( CData_Dens, CSize_Flu_Temp, CStart_Flu, CRange, &FData_Flu[DENS][0][0][0], 
                   FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, EnsureMonotonicity_Yes );

//    interpolate phase
      Interpolate( CData_Real, CSize_Flu_Temp, CStart_Flu, CRange, &FData_Flu[REAL][0][0][0], 
                   FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_Yes, EnsureMonotonicity_No );
   }      

   else // if ( OPT__INT_PHASE )
   {
      for (int v=0; v<NCOMP; v++)
      Interpolate( CData_Flu+v*CSize_Flu1v, CSize_Flu_Temp, CStart_Flu, CRange, &FData_Flu[v][0][0][0], 
                   FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, Monotonicity[v] );
   }

   if ( OPT__INT_PHASE )
   {
//    retrieve real and imaginary parts
      real Amp, Phase, Rho;

      for (int k=0; k<FSize; k++)
      for (int j=0; j<FSize; j++)
      for (int i=0; i<FSize; i++)
      {
         Phase = FData_Flu[REAL][k][j][i];
         Rho   = FData_Flu[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho < (real)0.0 )
         {
            if (  FABS( Rho ) < TINY_VALUE  )
            {
               FData_Flu[DENS][k][j][i] = (real)0.0;
               Rho                      = (real)0.0;
            }
            else
               Aux_Error( ERROR_INFO, "negative density (%14.7e) is obtained in %s !!\n", Rho, __FUNCTION__ );
         }

         Amp                      = SQRT( Rho );
         FData_Flu[REAL][k][j][i] = Amp*COS( Phase );
         FData_Flu[IMAG][k][j][i] = Amp*SIN( Phase );
      }
   }

#  else // #if ( MODEL == ELBDM )

   for (int v=0; v<NCOMP; v++)
   Interpolate( CData_Flu+v*CSize_Flu1v, CSize_Flu_Temp, CStart_Flu, CRange, &FData_Flu[v][0][0][0], 
                FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, Monotonicity[v] );

#  endif // #if ( MODEL == ELBDM ) ... else 

#  ifdef GRAVITY
   const int CSize_Pot_Temp[3] = { CSize_Pot, CSize_Pot, CSize_Pot };
   real *const CData_Pot = CData + NCOMP*CSize_Flu*CSize_Flu*CSize_Flu;

   Interpolate( CData_Pot, CSize_Pot_Temp, CStart_Pot, CRange, &FData_Pot[0][0][0],
                FSize_Temp, FStart,     1, OPT__REF_POT_INT_SCHEME, PhaseUnwrapping_No, EnsureMonotonicity_No );
#  endif


// 3.3 copy data from FData_XXX to patch pointers
   int SonPID, Disp_i, Disp_j, Disp_k, I, J, K;

   for (int LocalID=0; LocalID<8; LocalID++)
   {
      SonPID = SonPID0 + LocalID;
      Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE ); 
      Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE ); 
      Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE ); 
         
//    fluid data
      for (int v=0; v<NCOMP; v++)         {
      for (int k=0; k<PATCH_SIZE; k++)    {  K = k + Disp_k;
      for (int j=0; j<PATCH_SIZE; j++)    {  J = j + Disp_j;
      for (int i=0; i<PATCH_SIZE; i++)    {  I = i + Disp_i;

         amr->patch[FSg_Flu][SonLv][SonPID]->fluid[v][k][j][i] = FData_Flu[v][K][J][I];

      }}}}

//    potential data
#     ifdef GRAVITY
      for (int k=0; k<PATCH_SIZE; k++)    {  K = k + Disp_k;
      for (int j=0; j<PATCH_SIZE; j++)    {  J = j + Disp_j;
      for (int i=0; i<PATCH_SIZE; i++)    {  I = i + Disp_i;

         amr->patch[FSg_Pot][SonLv][SonPID]->pot[k][j][i] = FData_Pot[K][J][I];

      }}}
#     endif

//    rescale real and imaginary parts to get the correct density in ELBDM if OPT__INT_PHASE is off
#     if ( MODEL == ELBDM )
      real Real, Imag, Rho_Corr, Rho_Wrong, Rescale;

      if ( !OPT__INT_PHASE )
      for (int k=0; k<PATCH_SIZE; k++)
      for (int j=0; j<PATCH_SIZE; j++)
      for (int i=0; i<PATCH_SIZE; i++)
      {
         Real      = amr->patch[FSg_Flu][SonLv][SonPID]->fluid[REAL][k][j][i];
         Imag      = amr->patch[FSg_Flu][SonLv][SonPID]->fluid[IMAG][k][j][i];
         Rho_Wrong = Real*Real + Imag*Imag;
         Rho_Corr  = amr->patch[FSg_Flu][SonLv][SonPID]->fluid[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
         {
            amr->patch[FSg_Flu][SonLv][SonPID]->fluid[DENS][k][j][i] = (real)0.0;
            Rescale = (real)0.0;
         }
         else
            Rescale = SQRT( Rho_Corr/Rho_Wrong );

         amr->patch[FSg_Flu][SonLv][SonPID]->fluid[REAL][k][j][i] *= Rescale;
         amr->patch[FSg_Flu][SonLv][SonPID]->fluid[IMAG][k][j][i] *= Rescale;
      }
#     endif
   } // for (int LocalID=0; LocalID<8; LocalID++)


// 4. pass particles from father to son
#  ifdef PARTICLE
   Par_PassParticle2Son( FaLv, FaPID );
#  endif


// free memory
   if ( FaIsHome )   delete [] CData;
   delete [] FData_Flu;
#  ifdef GRAVITY
   delete [] FData_Pot;
#  endif


   return SonPID0;

} // FUNCTION : AllocateSonPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  DeallocateSonPatch  
// Description :  Deallocate eight son patches at FaLv+1
//
// Note        :  1. Just to avoid duplicate code segment 
//                2. List "NewSonPID0_Real" will be reset
//
// Parameter   :  FaLv              : Targeted refinement level to be refined
//                FaPID             : Father patch index to remove sons
//                NNew_Real0        : Number of newly-allocated real patches with LocalID==0
//                NewSonPID0_Real   : List recording the indices of all newly-allocated real patches 
//                                    with LocalID==0
//                SwitchIdx         : Element in NewSonPID0_Real to be reset
//-------------------------------------------------------------------------------------------------------
void DeallocateSonPatch( const int FaLv, const int FaPID, const int NNew_Real0, int NewSonPID0_Real[], 
                         int SwitchIdx )
{

   const int MirSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const int SonLv      = FaLv + 1;
   const int GraLv      = FaLv + 2;
   const int SonPID0    = amr->patch[0][FaLv][FaPID]->son;
   const int NewPID0    = SonPID0;
   const int OldPID0    = amr->num[SonLv] - 8;

   int NewPID, OldPID, OldGraPID0, OldFaPID, SibPID, MaxSonPID0, MaxIdx, TempSonPID0;


// check : son patch must exist and be real patches
#  ifdef GAMER_DEBUG
   if ( SonPID0 < 0  ||  SonPID0 >= amr->NPatchComma[SonLv][1] )
      Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, SonNReal %d, incorrect SonPID0 = %d !!\n", 
                 FaLv, FaPID, amr->NPatchComma[SonLv][1], SonPID0 );
#  endif


// pass particles from sons to father
#  ifdef PARTICLE
   Par_PassParticle2Father( FaLv, FaPID );
#  endif


// deallocate the unflagged child patches and reset sibling indices to -1
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
   {
      for (int s=0; s<26; s++)   
      {
         SibPID = amr->patch[0][SonLv][SonPID]->sibling[s];
         if ( SibPID >= 0 )   amr->patch[0][SonLv][SibPID]->sibling[ MirSib[s] ] = -1;
      }

      amr->pdelete( SonLv, SonPID );
   }

// record NPatchComma   
   amr->NPatchComma[SonLv][1] -= 8;

// construct relation : father -> son
   amr->patch[0][FaLv][FaPID]->son = -1;

// relink the child patch pointers so that no patch indices are skipped
   if ( NewPID0 != OldPID0 )  
   {
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         NewPID = NewPID0 + LocalID;
         OldPID = OldPID0 + LocalID;

         for (int Sg=0; Sg<2; Sg++)
         {
//          relink pointers
            amr->patch[Sg][SonLv][NewPID] = amr->patch[Sg][SonLv][OldPID];

//          set redundant patch pointers as NULL
            amr->patch[Sg][SonLv][OldPID] = NULL; 
         }

//       reconstruct relation : grandson -> son
         OldGraPID0 = amr->patch[0][SonLv][NewPID]->son;
         if ( OldGraPID0 >= 0 )
         {
#           ifdef GAMER_DEBUG
            if ( GraLv >= NLEVEL )
               Aux_Error( ERROR_INFO, "SonLv %d, NewPID %d, OldGraPID0 %d, GraLv %d >= NLEVEL (%d) !!\n",
                          SonLv, NewPID, OldGraPID0, GraLv, NLEVEL );
#           endif

            for (int OldGraPID=OldGraPID0; OldGraPID<OldGraPID0+8; OldGraPID++)
               amr->patch[0][GraLv][OldGraPID]->father = NewPID;
         }
         
//       reconstruct relation : sibling
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][SonLv][NewPID]->sibling[s];
            if ( SibPID >= 0 )   amr->patch[0][SonLv][SibPID]->sibling[ MirSib[s] ] = NewPID;
         }  
      } // for (int LocalID=0; LocalID<8; LocalID++)

//    reconstruct relation : father -> son
      OldFaPID = amr->patch[0][SonLv][NewPID0]->father;
#     ifdef GAMER_DEBUG
      if ( OldFaPID == -1 )   
         Aux_Error( ERROR_INFO, "SonLv %d, NewPID0 %d -> father == -1 !!\n", SonLv, NewPID0 );
#     endif
      amr->patch[0][FaLv][OldFaPID]->son = NewPID0;

//    reset the list "NewSonPID0_Real"
      if ( NNew_Real0 > 0 )
      {
         if ( SwitchIdx >= 0 )   // not all elements in NewSonPID0_Real have been reset
         {
#           ifdef GAMER_DEBUG
            if ( NewSonPID0_Real[SwitchIdx] != OldPID0 )
               Aux_Error( ERROR_INFO, "SonLv %d, NewSonPID0_Real[%d] (%d) != OldPID0 (%d) !!\n", 
                          SonLv, SwitchIdx, NewSonPID0_Real[SwitchIdx], OldPID0 );
#           endif

            NewSonPID0_Real[SwitchIdx] = NewPID0;

//          put the largest element on "NewSonPID0_Real[NNew_Real0-1]"
            if ( SwitchIdx-1 < 0 )  
            {
               MaxSonPID0 = -__INT_MAX__;
               MaxIdx     = -__INT_MAX__;

               for (int t=0; t<NNew_Real0; t++)
               {
                  if ( NewSonPID0_Real[t] > MaxSonPID0 )
                  {
                     MaxSonPID0 = NewSonPID0_Real[t];
                     MaxIdx     = t;
                  }
               }

               TempSonPID0                   = NewSonPID0_Real[NNew_Real0-1];
               NewSonPID0_Real[NNew_Real0-1] = MaxSonPID0;
               NewSonPID0_Real[MaxIdx      ] = TempSonPID0;
            }
         } // if ( SwitchIdx >= 0 )

         else  // all elements in NewSonPID0_Real have been reset --> need to resort the list
         {
            SwitchIdx = NNew_Real0 - 1;

            if ( NewSonPID0_Real[SwitchIdx] == OldPID0 )
            {
               NewSonPID0_Real[SwitchIdx] = NewPID0;

//             put the largest element on "NewSonPID0_Real[NNew_Real0-1]"
               MaxSonPID0 = -__INT_MAX__;
               MaxIdx     = -__INT_MAX__;

               for (int t=0; t<NNew_Real0; t++)
               {
                  if ( NewSonPID0_Real[t] > MaxSonPID0 )
                  {
                     MaxSonPID0 = NewSonPID0_Real[t];
                     MaxIdx     = t;
                  }
               }

               TempSonPID0                   = NewSonPID0_Real[NNew_Real0-1];
               NewSonPID0_Real[NNew_Real0-1] = MaxSonPID0;
               NewSonPID0_Real[MaxIdx      ] = TempSonPID0;
            }
         } // if ( SwitchIdx >= 0 ) ... else ...
      } // if ( NNew_Real0 > 0 )
   } // if ( NewPID0 != OldPID0 )

} // FUNCTION : DeallocateSonPatch



#endif // #ifdef LOAD_BALANCE
