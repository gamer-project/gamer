#ifdef LOAD_BALANCE

#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllocateBufferPatch_Sibling_Base
// Description :  Allocate sibling-buffer patches for all real patches at the base level
//
// Note        :  1. This function should be invoked BEFORE LB_AllocateBufferPatch_Father()
//                2. This function is invoked by LB_AllocateBufferPatch_Sibling()
//                3. This function should give exactly the same results as LB_AllocateBufferPatch_Sibling( 0 )
//                   but it is faster
//                4. All buffer patches should be removed in advance
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void LB_AllocateBufferPatch_Sibling_Base()
{

// check the NPatchComma recording
   if ( amr->NPatchComma[0][1] != amr->num[0] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[0][1] (%d) != amr->num[0] (%d)\" !!\n",
                 amr->NPatchComma[0][1], amr->num[0] );


// 1. construct the displacement array of padded 1D corner coordinates
// ==========================================================================================
   const int  Padded              = 1<<NLEVEL;
   const int  BoxNScale_Padded[3] = { amr->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                      amr->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                      amr->BoxScale[2]/PATCH_SIZE + 2*Padded };  //normalized and padded BoxScale
   const int  Scale2              = 2*amr->scale[0];
   const long dr[3]               = { (long)Scale2,
                                      (long)Scale2*BoxNScale_Padded[0],
                                      (long)Scale2*BoxNScale_Padded[0]*BoxNScale_Padded[1] };

   long Cr1D_Disp[3][3][3];

   for (int k=-1; k<=1; k++)
   for (int j=-1; j<=1; j++)
   for (int i=-1; i<=1; i++)
      Cr1D_Disp[k+1][j+1][i+1] = (long)i*dr[0] + (long)j*dr[1] + (long)k*dr[2];



// 2. prepare the allocate list
// ==========================================================================================
   const int MemUnit = amr->NPatchComma[0][1];     // set arbitrarily
   const int PScale  = PATCH_SIZE*amr->scale[0];
   const int PGScale = 2*PScale;

   int  Cr[3], TRank, *Cr0;
   long LBIdx;
   bool Internal, Allocate;

   int    MemSize = MemUnit, NAlloc_Temp = 0;
   ulong *FaPaddedCr1D = (ulong*)malloc( MemSize*sizeof(ulong) );


   for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)
   {
      Cr0 = amr->patch[0][0][PID0]->corner;

      for (int k=-1; k<=1; k++)   {  Cr[2] = Cr0[2] + k*PGScale;
      for (int j=-1; j<=1; j++)   {  Cr[1] = Cr0[1] + j*PGScale;
      for (int i=-1; i<=1; i++)   {  Cr[0] = Cr0[0] + i*PGScale;

//       determine the host rank of the target patch
         LBIdx = LB_Corner2Index( 0, Cr, CHECK_OFF );
         TRank = LB_Index2Rank( 0, LBIdx, CHECK_ON );


//       criteria to allocate a buffer patch:
//       --> internal patch: residing in another rank
//           external patch: the external direction is periodic
//       --> internal patches are defined as those with corner[] within the simulation domain
         Internal = true;
         Allocate = true;

         for (int d=0; d<3; d++)
         {
            if ( Cr[d] < 0  ||  Cr[d] >= amr->BoxScale[d] )
            {
               Internal = false;

               if ( OPT__BC_FLU[2*d] != BC_FLU_PERIODIC )
               {
                  Allocate = false;
                  break;
               }
            }
         }

         if ( Internal )   Allocate = ( TRank != MPI_Rank );


//       record the padded 1D corner coordinates of the sibling-buffer patch to be allocated
         if ( Allocate )
         {
            if ( NAlloc_Temp >= MemSize )
            {
               MemSize      += MemUnit;
               FaPaddedCr1D  = (ulong*)realloc( FaPaddedCr1D, MemSize*sizeof(ulong) );
            }

//###NOTE: Disp = i*dr[0] + j*dr[1] + k*dr[2] can be negative! But it's OK to conduct PaddedCr1D + (ulong)Disp
//         as long as we guarantee "PaddedCr1D + Disp >= 0"
//         --> ulong(Disp) = Disp + UINT_MAX + 1 (if Disp < 0; ==> reduced modulo)
//         --> PaddedCr1D + (ulong)Disp = PaddedCr1D + Disp + UINT_MAX + 1 = PaddedCr1D + Disp + UINT_MAX + 1 - (UINT_MAX + 1)
//                                      = PaddedCr1D + Disp
//             (because PaddedCr1D + Disp >= 0; ==> reduced modulo again)
            FaPaddedCr1D[ NAlloc_Temp ++ ] = amr->patch[0][0][PID0]->PaddedCr1D + (ulong)Cr1D_Disp[k+1][j+1][i+1];
         }
      }}} // k,j,i
   } // for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)



// 3. sort the allocate list and remove duplicates
// ==========================================================================================
   int NAlloc = ( NAlloc_Temp > 0 ) ? 1 : 0;

   Mis_Heapsort( NAlloc_Temp, FaPaddedCr1D, NULL );

   for (int t=1; t<NAlloc_Temp; t++)
      if ( FaPaddedCr1D[t] != FaPaddedCr1D[t-1] )  FaPaddedCr1D[ NAlloc ++ ] = FaPaddedCr1D[t];



// 4. allocate sibling-buffer patches
// ==========================================================================================
   int FaCr3D[3];

   for (int t=0; t<NAlloc; t++)
   {
//    get the 3D corner coordinates
      Mis_Idx1D2Idx3D( BoxNScale_Padded, FaPaddedCr1D[t], FaCr3D );
      for (int d=0; d<3; d++)    FaCr3D[d] = ( FaCr3D[d] - Padded )*PATCH_SIZE;

//    data array is not allocated here
      amr->pnew( 0, FaCr3D[0],        FaCr3D[1],        FaCr3D[2],        -1, false, false );
      amr->pnew( 0, FaCr3D[0]+PScale, FaCr3D[1],        FaCr3D[2],        -1, false, false );
      amr->pnew( 0, FaCr3D[0],        FaCr3D[1]+PScale, FaCr3D[2],        -1, false, false );
      amr->pnew( 0, FaCr3D[0],        FaCr3D[1],        FaCr3D[2]+PScale, -1, false, false );
      amr->pnew( 0, FaCr3D[0]+PScale, FaCr3D[1]+PScale, FaCr3D[2],        -1, false, false );
      amr->pnew( 0, FaCr3D[0],        FaCr3D[1]+PScale, FaCr3D[2]+PScale, -1, false, false );
      amr->pnew( 0, FaCr3D[0]+PScale, FaCr3D[1],        FaCr3D[2]+PScale, -1, false, false );
      amr->pnew( 0, FaCr3D[0]+PScale, FaCr3D[1]+PScale, FaCr3D[2]+PScale, -1, false, false );

      amr->NPatchComma[0][2] += 8;

   } // for (int t=0; t<NAlloc; t++)


// record NPatchComma
   for (int m=3; m<28; m++)   amr->NPatchComma[0][m] = amr->NPatchComma[0][2];

   if ( amr->NPatchComma[0][2] != amr->num[0] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[0][2] (%d) != amr->num[0] (%d)\" !!\n",
                 amr->NPatchComma[0][2], amr->num[0] );



// 5. record the padded 1D corner coordinates
//    --> can be overwritten by LB_AllocateBufferPatch_Father()
// ==========================================================================================
   const int NPatch = amr->NPatchComma[0][2];

   amr->LB->PaddedCr1DList         [0] = (ulong*)realloc( amr->LB->PaddedCr1DList         [0],
                                                          NPatch*sizeof(ulong) );
   amr->LB->PaddedCr1DList_IdxTable[0] = (int*  )realloc( amr->LB->PaddedCr1DList_IdxTable[0],
                                                          NPatch*sizeof(int ) );

   for (int PID=0; PID<NPatch; PID++)
      amr->LB->PaddedCr1DList[0][PID] = amr->patch[0][0][PID]->PaddedCr1D;

   Mis_Heapsort( NPatch, amr->LB->PaddedCr1DList[0], amr->LB->PaddedCr1DList_IdxTable[0] );

// check : no duplicate patches
#  ifdef GAMER_DEBUG
   for (int t=1; t<amr->num[0]; t++)
   {
      if ( amr->LB->PaddedCr1DList[0][t] == amr->LB->PaddedCr1DList[0][t-1] )
         Aux_Error( ERROR_INFO, "duplicate patches at lv 0, PaddedCr1D %lu, PID = %d and %d !!\n",
                    amr->LB->PaddedCr1DList[0][t], amr->LB->PaddedCr1DList_IdxTable[0][t],
                    amr->LB->PaddedCr1DList_IdxTable[0][t-1] );
   }
#  endif


// free memory
   free( FaPaddedCr1D );

} // FUNCTION : LB_AllocateBufferPatch_Sibling_Base



#endif // #ifdef LOAD_BALANCE
