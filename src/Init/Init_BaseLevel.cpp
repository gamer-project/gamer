#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_BaseLevel
// Description :  Construct the base level during initialization
//-------------------------------------------------------------------------------------------------------
void Init_BaseLevel()
{

   const int NPatch[3] = { NX0[0]/PATCH_SIZE, NX0[1]/PATCH_SIZE, NX0[2]/PATCH_SIZE };
   const int scale0    = amr->scale[0];
   const int Width     = PATCH_SIZE*scale0;

   int Cr[3];


// allocate the real base-level patches
   for (int Pz=0; Pz<NPatch[2]; Pz+=2)    {  Cr[2] = MPI_Rank_X[2]*NX0[2]*scale0 + Pz*PATCH_SIZE*scale0;
   for (int Py=0; Py<NPatch[1]; Py+=2)    {  Cr[1] = MPI_Rank_X[1]*NX0[1]*scale0 + Py*PATCH_SIZE*scale0;
   for (int Px=0; Px<NPatch[0]; Px+=2)    {  Cr[0] = MPI_Rank_X[0]*NX0[0]*scale0 + Px*PATCH_SIZE*scale0;

      amr->pnew( 0, Cr[0],       Cr[1],       Cr[2],       -1, true, true, true );
      amr->pnew( 0, Cr[0]+Width, Cr[1],       Cr[2],       -1, true, true, true );
      amr->pnew( 0, Cr[0],       Cr[1]+Width, Cr[2],       -1, true, true, true );
      amr->pnew( 0, Cr[0],       Cr[1],       Cr[2]+Width, -1, true, true, true );
      amr->pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2],       -1, true, true, true );
      amr->pnew( 0, Cr[0],       Cr[1]+Width, Cr[2]+Width, -1, true, true, true );
      amr->pnew( 0, Cr[0]+Width, Cr[1],       Cr[2]+Width, -1, true, true, true );
      amr->pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, -1, true, true, true );

   }}}

   for (int m=1; m<28; m++)   amr->NPatchComma[0][m] = amr->num[0];


// allocate the base-level buffer patches
   Buf_AllocateBufferPatch( amr, 0 );

// set up the BaseP List
   Init_RecordBasePatch();

// set up the BounP_IDMap[0]
   Buf_RecordBoundaryPatch( 0 );

// construct the sibling relation for the base level (including the buffer patches)
   SiblingSearch( 0 );

// get the IDs of patches for sending and receiving data between neighbor ranks
   Buf_RecordExchangeDataPatchID( 0 );

// find the base-level home patches of all particles
#  ifdef PARTICLE
   const bool OldParOnly_Yes = true;
   Par_FindHomePatch_UniformGrid( 0, OldParOnly_Yes, NULL_INT, NULL );
#  endif

} // FUNCTION : Init_BaseLevel
