#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_AllocateBufferPatch_Base
// Description :  Allocate buffer patches for the base level
//
// Note        :  1. corner[] recorded in the buffer patches are already non-periodic
//                2. Invoked by Buf_AllocateBufferPatch()
//                3. No buffer patches will be allocated for patches lying outside the simulation domain if
//                   non-periodic B.C. is adopted
//
// Parameter   :  Tamr : Target AMR_t pointer
//-------------------------------------------------------------------------------------------------------
void Buf_AllocateBufferPatch_Base( AMR_t *Tamr )
{

   const int scale0 = Tamr->scale[0];
   const int Width  = PATCH_SIZE*scale0;

   int Cr0[3], Cr[3], loop_size[3];
   bool AllocData[8];      // allocate data or not


   for (int s=0; s<26; s++)
   {
//    no external buffer patches will be allocated for the non-periodic B.C.
      if ( MPI_SibRank[s] <= SIB_OFFSET_NONPERIODIC )    continue;


//    determine the buffer patches that must store physical data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         if (  TABLE_01( s, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
               TABLE_01( s, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
               TABLE_01( s, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )
            AllocData[LocalID] = false;

         else
            AllocData[LocalID] = true;
      }


//    allocate the buffer patches
      for (int d=0; d<3; d++)
      {
         Cr0      [d] = TABLE_01( s, 'x'+d, -2*PATCH_SIZE*scale0, 0, NX0[d]*scale0 ) + MPI_Rank_X[d]*NX0[d]*scale0;
         loop_size[d] = TABLE_01( s, 'x'+d, 1, NX0[d]/PS2, 1 );
      }

      for (int k=0; k<loop_size[2]; k++)     {  Cr[2] = Cr0[2] + k*2*Width;
      for (int j=0; j<loop_size[1]; j++)     {  Cr[1] = Cr0[1] + j*2*Width;
      for (int i=0; i<loop_size[0]; i++)     {  Cr[0] = Cr0[0] + i*2*Width;

         Tamr->pnew( 0, Cr[0],       Cr[1],       Cr[2],       -1, AllocData[0], AllocData[0] );
         Tamr->pnew( 0, Cr[0]+Width, Cr[1],       Cr[2],       -1, AllocData[1], AllocData[1] );
         Tamr->pnew( 0, Cr[0],       Cr[1]+Width, Cr[2],       -1, AllocData[2], AllocData[2] );
         Tamr->pnew( 0, Cr[0],       Cr[1],       Cr[2]+Width, -1, AllocData[3], AllocData[3] );
         Tamr->pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2],       -1, AllocData[4], AllocData[4] );
         Tamr->pnew( 0, Cr[0],       Cr[1]+Width, Cr[2]+Width, -1, AllocData[5], AllocData[5] );
         Tamr->pnew( 0, Cr[0]+Width, Cr[1],       Cr[2]+Width, -1, AllocData[6], AllocData[6] );
         Tamr->pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, -1, AllocData[7], AllocData[7] );

         Tamr->NPatchComma[0][s+2] += 8;

      }}}

      for (int n=s+3; n<28; n++)    Tamr->NPatchComma[0][n] = Tamr->num[0];

   } // for (int s=0; s<26; s++)

   if ( Tamr->NPatchComma[0][27] != Tamr->num[0] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[0][27] (%d) != amr->num[0] (%d) !!\n",
                 Tamr->NPatchComma[0][27], Tamr->num[0] );

} // FUNCTION : Buf_AllocateBufferPatch_Base



#endif // #ifndef SERIAL
