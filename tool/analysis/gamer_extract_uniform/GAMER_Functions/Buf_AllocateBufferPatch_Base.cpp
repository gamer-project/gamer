#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_AllocateBufferPatch_Base
// Description :  Allocate buffer patches for the base level
//
// Note        :  a. The "corner" recorded in the buffer patches are not correct (NOT periodic). They are
//                   useful only for the function "Buf_RecordBasePatchID"
//                b. Invoked by the function "Buf_AllocateBufferPatch"
//-------------------------------------------------------------------------------------------------------
void Buf_AllocateBufferPatch_Base()
{

   const int scale = amr.scale[0];
   const int Width = PATCH_SIZE*scale;

   int Cr0[3], Cr[3], loop_size[3];
   bool AllocateData[8];      // allocate data or not


   bool AllocateData2[8];     // check whether of not the patch is within the candidate box


   for (int s=0; s<26; s++)
   {
//    no external buffer patches will be allocated for the non-periodic B.C.
      if ( SibRank[s] <= SIB_OFFSET_NONPERIODIC )  continue;


//    determine the buffer patches that must store physical data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         if (  TABLE_01( s, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
               TABLE_01( s, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
               TABLE_01( s, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )

            AllocateData[LocalID] = false;

         else

            AllocateData[LocalID] = true;
      }


      for (int d=0; d<3; d++)
      {
         Cr0      [d] = TABLE_01( s, 'x'+d, -2*PATCH_SIZE*scale, 0, NX0[d]*scale ) + MyRank_X[d]*NX0[d]*scale;
         loop_size[d] = TABLE_01( s, 'x'+d, 1, NX0[d]/PS2, 1 );
      }

      for (int k=0; k<loop_size[2]; k++)  {  Cr[2] = Cr0[2] + k*2*Width;
      for (int j=0; j<loop_size[1]; j++)  {  Cr[1] = Cr0[1] + j*2*Width;
      for (int i=0; i<loop_size[0]; i++)  {  Cr[0] = Cr0[0] + i*2*Width;

         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int PS0 = PATCH_SIZE*amr.scale[0];
            int Corner[3];

            for (int d=0; d<3; d++)    Corner[d] = Cr[d] + TABLE_02( LocalID, 'x'+d, 0, Width );

            AllocateData2[LocalID] = ( AllocateData[LocalID] && WithinCandidateBox(Corner,PS0,CanBuf) );
         }

         amr.pnew( 0, Cr[0],       Cr[1],       Cr[2],       -1, AllocateData2[0] );
         amr.pnew( 0, Cr[0]+Width, Cr[1],       Cr[2],       -1, AllocateData2[1] );
         amr.pnew( 0, Cr[0],       Cr[1]+Width, Cr[2],       -1, AllocateData2[2] );
         amr.pnew( 0, Cr[0],       Cr[1],       Cr[2]+Width, -1, AllocateData2[3] );
         amr.pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2],       -1, AllocateData2[4] );
         amr.pnew( 0, Cr[0],       Cr[1]+Width, Cr[2]+Width, -1, AllocateData2[5] );
         amr.pnew( 0, Cr[0]+Width, Cr[1],       Cr[2]+Width, -1, AllocateData2[6] );
         amr.pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, -1, AllocateData2[7] );

         NPatchComma[0][s+2] += 8;

      }}}

      for (int n=s+3; n<28; n++)    NPatchComma[0][n] = amr.num[0];
   }

   if ( NPatchComma[0][27] != amr.num[0] )
   {
      fprintf( stderr, "ERROR : \"NPatchComma[0][27] (%d) != amr.num[0] (%d)\" !!\n",
               NPatchComma[0][27], amr.num[0] );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      MPI_Exit();
   }

}
