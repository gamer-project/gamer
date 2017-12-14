#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_BoundaryPatch
// Description :  Record the IDs and positions of the patches near the sub-domain boundaries in the 
//                BounP_IDList and BounP_PosList, respectively
//
// Note        :  This function is only an auxiliary function, and it is kept only for verifying the results
//                of the functions "RecordBoundaryPatch_Base" and "RecordBoundaryPatch" 
//-------------------------------------------------------------------------------------------------------
void Aux_Record_BoundaryPatch( const int lv, int *NList, int **IDList, int **PosList )
{

   const int NP0[3] = { NX0[0]/PATCH_SIZE, NX0[1]/PATCH_SIZE, NX0[2]/PATCH_SIZE }; 
   const int NP [3] = { NP0[0]*(1<<lv)+4, NP0[1]*(1<<lv)+4, NP0[2]*(1<<lv)+4 };

   const int scale0 = amr->scale[ 0];      
   const int scale  = amr->scale[lv];        // Pos      : Patch position in (x,y,z) direction
   int Pos[3], Pos2, LCorner[3];             // Pos2     : 1-D Patch position
                                             // LCorner  : local coordinate of patch corner in (x,y,z) direction


// initialize counter as zero
   for (int s=0; s<26; s++)   NList[s] = 0;


// scan over all patches at level "lv" to record the IDS of boundary patches
   for (int P=0; P<amr->num[lv]; P++) 
   {
      LCorner[0]  = amr->patch[0][lv][P]->corner[0] - MPI_Rank_X[0]*NX0[0]*scale0;
      LCorner[1]  = amr->patch[0][lv][P]->corner[1] - MPI_Rank_X[1]*NX0[1]*scale0;
      LCorner[2]  = amr->patch[0][lv][P]->corner[2] - MPI_Rank_X[2]*NX0[2]*scale0;

      Pos[0]      = LCorner[0]/(PATCH_SIZE*scale) + 2;
      Pos[1]      = LCorner[1]/(PATCH_SIZE*scale) + 2;
      Pos[2]      = LCorner[2]/(PATCH_SIZE*scale) + 2;


      if ( Pos[0] == 2 )
      {
//             face 0
               Pos2                          = Pos[2]*NP[1] + Pos[1];
               IDList [0][ NList[0] ]        = P;
               PosList[0][ NList[0] ]        = Pos2;
               NList  [0] ++;

         if ( Pos[1] == 2 )
         {
//             face 6
               Pos2                          = Pos[2];
               IDList [6][ NList[6] ]        = P;
               PosList[6][ NList[6] ]        = Pos2;
               NList  [6] ++;

            if ( Pos[2] == 2 )
            {
//             face 18
               Pos2                          = 0;
               IDList [18][ NList[18] ]      = P;
               PosList[18][ NList[18] ]      = Pos2;
               NList  [18] ++;
            }
            else if ( Pos[2] == NP[2]-3 )
            {
//             face 22
               Pos2                          = 0;
               IDList [22][ NList[22] ]      = P;
               PosList[22][ NList[22] ]      = Pos2;
               NList  [22] ++;
            }

         } // if ( Pos[1] == 0 )

         else if ( Pos[1] == NP[1]-3 )
         {
//             face 8
               Pos2                          = Pos[2];
               IDList [8][ NList[8] ]        = P;
               PosList[8][ NList[8] ]        = Pos2;
               NList  [8] ++;

            if ( Pos[2] == 2 )
            {
//             face 20
               Pos2                          = 0;
               IDList [20][ NList[20] ]      = P;
               PosList[20][ NList[20] ]      = Pos2;
               NList  [20] ++;
            }
            else if ( Pos[2] == NP[2]-3 )
            {
//             face 24
               Pos2                          = 0;
               IDList [24][ NList[24] ]      = P;
               PosList[24][ NList[24] ]      = Pos2;
               NList  [24] ++;
            }

         } // else if ( Pos[1] == NP[1]-3 )
      } // if ( Pos[0] == 2 )


      else if ( Pos[0] == NP[0]-3 )
      {
//             face 1
               Pos2                          = Pos[2]*NP[1] + Pos[1];
               IDList [1][ NList[1] ]        = P;
               PosList[1][ NList[1] ]        = Pos2;
               NList  [1] ++;

         if ( Pos[1] == 2 )
         {
//             face 7
               Pos2                          = Pos[2];
               IDList [7][ NList[7] ]        = P;
               PosList[7][ NList[7] ]        = Pos2;
               NList  [7] ++;

            if ( Pos[2] == 2 )
            {
//             face 19
               Pos2                          = 0;
               IDList [19][ NList[19] ]      = P;
               PosList[19][ NList[19] ]      = Pos2;
               NList  [19] ++;

            }
            else if ( Pos[2] == NP[2]-3 )
            {
//             face 23
               Pos2                          = 0;
               IDList [23][ NList[23] ]      = P;
               PosList[23][ NList[23] ]      = Pos2;
               NList  [23] ++;
            }

         } // if ( Pos[1] == 2 )

         else if ( Pos[1] == NP[1]-3 )
         {
//             face 9
               Pos2                          = Pos[2];
               IDList [9][ NList[9] ]        = P;
               PosList[9][ NList[9] ]        = Pos2;
               NList  [9] ++;

            if ( Pos[2] == 2 )
            {
//             face 21
               Pos2                          = 0;
               IDList [21][ NList[21] ]      = P;
               PosList[21][ NList[21] ]      = Pos2;
               NList  [21] ++;
            }
            else if ( Pos[2] == NP[2]-3 )
            {
//             face 25
               Pos2                          = 0;
               IDList [25][ NList[25] ]      = P;
               PosList[25][ NList[25] ]      = Pos2;
               NList  [25] ++;
            }

         } // else if ( Pos[1] == NP[1]-3 )
      } // else if ( Pos[0] == NP[0]-3 )


      if ( Pos[1] == 2 )
      {
//             face 2
               Pos2                          = Pos[2]*NP[0] + Pos[0];
               IDList [2][ NList[2] ]        = P;
               PosList[2][ NList[2] ]        = Pos2;
               NList  [2] ++;

         if ( Pos[2] == 2 )
         {
//             face 10
               Pos2                          = Pos[0];
               IDList [10][ NList[10] ]      = P;
               PosList[10][ NList[10] ]      = Pos2;
               NList  [10] ++;
         }
         else if ( Pos[2] == NP[2]-3 )
         {
//             face 12
               Pos2                          = Pos[0];
               IDList [12][ NList[12] ]      = P;
               PosList[12][ NList[12] ]      = Pos2;
               NList  [12] ++;
         }

      } // if ( Pos[1] == 2 )

      else if ( Pos[1] == NP[1]-3 )
      {
//             face 3
               Pos2                          = Pos[2]*NP[0] + Pos[0];
               IDList [3][ NList[3] ]        = P;
               PosList[3][ NList[3] ]        = Pos2;
               NList  [3] ++;

         if ( Pos[2] == 2 )
         {
//             face 11
               Pos2                          = Pos[0];
               IDList [11][ NList[11] ]      = P;
               PosList[11][ NList[11] ]      = Pos2;
               NList  [11] ++;
         }
         else if ( Pos[2] == NP[2]-3 )
         {
//             face 13
               Pos2                          = Pos[0];
               IDList [13][ NList[13] ]      = P;
               PosList[13][ NList[13] ]      = Pos2;
               NList  [13] ++;
         }

      } // else if ( Pos[1] == NP[1]-3 )


      if ( Pos[2] == 2 )
      {
//             face 4
               Pos2                          = Pos[1]*NP[0] + Pos[0];
               IDList [4][ NList[4] ]        = P;
               PosList[4][ NList[4] ]        = Pos2;
               NList  [4] ++;

         if ( Pos[0] == 2 )
         {
//             face 14
               Pos2                          = Pos[1];
               IDList [14][ NList[14] ]      = P;
               PosList[14][ NList[14] ]      = Pos2;
               NList  [14] ++;
         }
         else if ( Pos[0] == NP[0]-3 )
         {
//             face 16
               Pos2                          = Pos[1];
               IDList [16][ NList[16] ]      = P;
               PosList[16][ NList[16] ]      = Pos2;
               NList  [16] ++;
         }

      } // if ( Pos[2] == 2 )

      else if ( Pos[2] == NP[2]-3 )
      {
//             face 5
               Pos2                          = Pos[1]*NP[0] + Pos[0];
               IDList [5][ NList[5] ]        = P;
               PosList[5][ NList[5] ]        = Pos2;
               NList  [5] ++;

         if ( Pos[0] == 2 )
         {
//             face 15
               Pos2                          = Pos[1];
               IDList [15][ NList[15] ]      = P;
               PosList[15][ NList[15] ]      = Pos2;
               NList  [15] ++;
         }
         else if ( Pos[0] == NP[0]-3 )
         {
//             face 17
               Pos2                          = Pos[1];
               IDList [17][ NList[17] ]      = P;
               PosList[17][ NList[17] ]      = Pos2;
               NList  [17] ++;
         }

      } // else if ( Pos[2] == NP[2]-3 )

   } // for (int P=0; P<amr->num[lv]; P++)


// sort
   for (int s=0; s<26; s++)   Buf_SortBoundaryPatch( NList[s], IDList[s], PosList[s] );

} // FUNCTION : Aux_Record_BoundaryPatch



#endif // #ifndef SERIAL
