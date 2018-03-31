#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_RecordBasePatch
// Description :  Record the IDs of base-level patches in the BaseP array
//
// Note        :  BaseP array is used to
//                (1) construct the sibling relation at the base level --> SiblingSearch_Base
//                (2) find the father patch                            --> FindFather
//                (3) find the particle's home patch at the base level --> Par_FindHomePatch_Base
//-------------------------------------------------------------------------------------------------------
void Init_RecordBasePatch()
{

   const int NBuf        = 2;
   const int scale0      = amr->scale[0];
   const int NPatch1D[3] = { NX0[0]/PATCH_SIZE+2*NBuf, NX0[1]/PATCH_SIZE+2*NBuf, NX0[2]/PATCH_SIZE+2*NBuf };
   const int NotSetYet   = NULL_INT;

   int order[3], P2[8];    // order[3] : (i,j,k)th patch in the (x,y,z) direction


// 0. initialize BaseP[] so that we can check whether it's properly set later
   for (int t=0; t<NPatch1D[0]*NPatch1D[1]*NPatch1D[2]; t++ )  BaseP[t] = NotSetYet;


// 1. record all patches (real+buffer)
   for (int P=0; P<amr->num[0]; P+=8)
   {
      for (int d=0; d<3; d++)
         order[d] = (amr->patch[0][0][P]->corner[d] - MPI_Rank_X[d]*NX0[d]*scale0) / (PATCH_SIZE*scale0) + NBuf;

      P2[0] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+0);
      P2[1] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+1);
      P2[2] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+0);
      P2[3] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+0);
      P2[4] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);
      P2[5] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+0);
      P2[6] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+1);
      P2[7] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);

//    record the patch ID in the BaseP array
      for (int m=0; m<8; m++)    BaseP[ P2[m] ] = P + m;
   }


// 2. deal with the non-periodic B.C. (and periodicity in the serial code since it does not allocate any buffer patch)
   int Width[3], Disp1[3], NP, ID1;
#  ifdef SERIAL
   int i2, j2, k2, Disp2[3], ID2;
#  endif

   for (int Sib=0; Sib<26; Sib++)
   {
      for (int d=0; d<3; d++)
      {
         NP       = NX0[d]/PATCH_SIZE;
         Width[d] = TABLE_01( Sib, 'x'+d, NBuf, NP,   NBuf    );
         Disp1[d] = TABLE_01( Sib, 'x'+d, 0,    NBuf, NBuf+NP );
#        ifdef SERIAL
         Disp2[d] = TABLE_01( Sib, 'x'+d, NP,   0,    -NP     );
#        endif
      }

//    non-periodic B.C.
      if ( MPI_SibRank[Sib] <= SIB_OFFSET_NONPERIODIC )
      {
         for (int k1=Disp1[2]; k1<Disp1[2]+Width[2]; k1++)  {
         for (int j1=Disp1[1]; j1<Disp1[1]+Width[1]; j1++)  {
         for (int i1=Disp1[0]; i1<Disp1[0]+Width[0]; i1++)  {

            ID1 = ( k1*NPatch1D[1] + j1 )*NPatch1D[0] + i1;

//          ensure it has not been set already
            if ( BaseP[ID1] != NotSetYet )
               Aux_Error( ERROR_INFO, "BaseP[%d] = %d has been set already (Sib %d) !!\n", ID1, BaseP[ID1], Sib );

//          assuming that MPI_SibRank[] has been set to "SIB_OFFSET_NONPERIODIC-sib" for the non-periodic B.C.
//          --> set by Init_Parallelization()
            BaseP[ID1] = MPI_SibRank[Sib];
         }}}
      }

#     ifdef SERIAL
//    periodic B.C.
//    --> necessary only for the serial mode since it does not allocate any buffer patch
      else
      {
         for (int k1=Disp1[2]; k1<Disp1[2]+Width[2]; k1++)  {  k2 = k1 + Disp2[2];
         for (int j1=Disp1[1]; j1<Disp1[1]+Width[1]; j1++)  {  j2 = j1 + Disp2[1];
         for (int i1=Disp1[0]; i1<Disp1[0]+Width[0]; i1++)  {  i2 = i1 + Disp2[0];

            ID1 = ( k1*NPatch1D[1] + j1 )*NPatch1D[0] + i1;
            ID2 = ( k2*NPatch1D[1] + j2 )*NPatch1D[0] + i2;

            BaseP[ID1] = BaseP[ID2];
         }}}
      }
#     endif // #ifdef SERIAL
   } // for (int Sib=0; Sib<26; Sib++)


// 3. check if all elements of BaseP[] have been set properly
   for (int t=0; t<NPatch1D[0]*NPatch1D[1]*NPatch1D[2]; t++ ) {
      if ( BaseP[t] == NotSetYet )
         Aux_Error( ERROR_INFO, "BaseP[%d] has not been set !!\n", t );
   }

} // FUNCTION : Init_RecordBasePatch
