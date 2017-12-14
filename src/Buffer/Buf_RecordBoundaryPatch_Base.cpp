#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_RecordBoundaryPatch_Base
// Description :  Record the base-level patches near the sub-domain boundaries in
//                "amr->ParaVar->BounP_IDList[0]"
//
// Note        :  Invoked by the function "Buf_RecordBoundaryPatch"
//-------------------------------------------------------------------------------------------------------
void Buf_RecordBoundaryPatch_Base()
{

// BaseP_Len : the number of elements of the BaseP array in each direction
   const int BaseP_Len[3] = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };

   int ID_BaseP, ID_BounP, Width[3], Disp[3], ii, jj, kk;

   for (int s=0; s<26; s++)
   {
//    set up the "Width and Disp"
      for (int d=0; d<3; d++)
      {
         Width[d] = TABLE_01( s, 'x'+d, 1, BaseP_Len[d], 1 );
         Disp [d] = TABLE_01( s, 'x'+d, 2, 0, BaseP_Len[d]-3 );
      }


//    initialize the BounP_NList[0][s] as 0
      amr->ParaVar->BounP_NList[0][s] = 0;


//    deallocate memory (only for safety reason)
      if ( amr->ParaVar->BounP_IDList[0][s] != NULL )
      {
         delete [] amr->ParaVar->BounP_IDList[0][s];
         amr->ParaVar->BounP_IDList[0][s] = NULL;
      }

      if ( amr->ParaVar->BounP_PosList[0][s] != NULL )
      {
         delete [] amr->ParaVar->BounP_PosList[0][s];
         amr->ParaVar->BounP_PosList[0][s] = NULL;
      }


//    allocate memory (only necessary during the initialization)
      amr->ParaVar->BounP_IDList [0][s] = new int [ Width[0]*Width[1]*Width[2] ];
      amr->ParaVar->BounP_PosList[0][s] = new int [ Width[0]*Width[1]*Width[2] ];


//    fill up the arrays "BounP_IDList[0][s] and BounP_PosList[0][s]"
      for (int k=0; k<Width[2]; k++)   {  kk = k + Disp[2];
      for (int j=0; j<Width[1]; j++)   {  jj = j + Disp[1];
      for (int i=0; i<Width[0]; i++)   {  ii = i + Disp[0];

         ID_BaseP = (kk*BaseP_Len[1] + jj)*BaseP_Len[0] + ii;
         ID_BounP = (k *Width    [1] + j )*Width    [0] + i;

         if ( BaseP[ ID_BaseP ] >= 0 )    // there are no external buffer patches for non-periodic B.C.
         {
            amr->ParaVar->BounP_IDList [0][s][ amr->ParaVar->BounP_NList[0][s] ] = BaseP[ ID_BaseP ];
            amr->ParaVar->BounP_PosList[0][s][ amr->ParaVar->BounP_NList[0][s] ] = ID_BounP;
            amr->ParaVar->BounP_NList  [0][s] ++;
         }

      }}}

   } // for (int s=0; s<26; s++)

} // FUNCTION : Buf_RecordBoundaryPatch_Base



#endif // #ifndef SERIAL
