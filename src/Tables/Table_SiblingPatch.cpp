#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_GetSibPID_Delta
// Description :  Get the differences between all sibling patch indices (SibPID) and the based index
//                SibPID_Based returned by TABLE_GetSibPID_Based()
//
// Note        :  1. Faster than using TABLE_03() and TABLE_04()
//                2. SibPID_Delta[] must be deallocated manually
//                3. All sibling patches SibPID of a given patch PID0 on level Lv can be iterated over as follows:
//
//                   int NSibPatch[26], *SibPID_Delta[26], SibPID_Based[26];
//                   TABLE_GetSibPID_Delta( NSibPatch, SibPID_Delta );
//
//                   TABLE_GetSibPID_Based( Lv, PID0, SibPID_Based );
//
//                   for (int s=0; s<26; s++) {
//                      const int SibPID0 = SibPID_Based[s];
//                      if ( SibPID0 >= 0 ) {
//                         for (int c=0; c<NSibPatch[s]; c++) {
//                            const int SibPID = SibPID0 + SibPID_Delta[s][c];
//                         }
//                      }
//                   }
//
// Parameter   :  NSibPatch    : Number of sibling patches along different directions
//                SibPID_Delta : SibPID - SibPID_Based
//
// Return      :  NSibPatch, SibPID_Delta
//-------------------------------------------------------------------------------------------------------
void TABLE_GetSibPID_Delta( int NSibPatch[], int *SibPID_Delta[] )
{

   for (int s= 0; s< 6; s++)  NSibPatch[s] = 4;
   for (int s= 6; s<18; s++)  NSibPatch[s] = 2;
   for (int s=18; s<26; s++)  NSibPatch[s] = 1;

   for (int s=0; s<26; s++)   SibPID_Delta[s] = new int [ NSibPatch[s] ];

   SibPID_Delta[ 0][0] = 0;
   SibPID_Delta[ 0][1] = 3;
   SibPID_Delta[ 0][2] = 5;
   SibPID_Delta[ 0][3] = 6;

   SibPID_Delta[ 1][0] = 0;
   SibPID_Delta[ 1][1] = 2;
   SibPID_Delta[ 1][2] = 3;
   SibPID_Delta[ 1][3] = 5;

   SibPID_Delta[ 2][0] = 0;
   SibPID_Delta[ 2][1] = 2;
   SibPID_Delta[ 2][2] = 3;
   SibPID_Delta[ 2][3] = 5;

   SibPID_Delta[ 3][0] = 0;
   SibPID_Delta[ 3][1] = 1;
   SibPID_Delta[ 3][2] = 3;
   SibPID_Delta[ 3][3] = 6;

   SibPID_Delta[ 4][0] = 0;
   SibPID_Delta[ 4][1] = 2;
   SibPID_Delta[ 4][2] = 3;
   SibPID_Delta[ 4][3] = 4;

   SibPID_Delta[ 5][0] = 0;
   SibPID_Delta[ 5][1] = 1;
   SibPID_Delta[ 5][2] = 2;
   SibPID_Delta[ 5][3] = 4;

   SibPID_Delta[ 6][0] = 0;
   SibPID_Delta[ 6][1] = 3;

   SibPID_Delta[ 7][0] = 0;
   SibPID_Delta[ 7][1] = 3;

   SibPID_Delta[ 8][0] = 0;
   SibPID_Delta[ 8][1] = 5;

   SibPID_Delta[ 9][0] = 0;
   SibPID_Delta[ 9][1] = 3;

   SibPID_Delta[10][0] = 0;
   SibPID_Delta[10][1] = 2;

   SibPID_Delta[11][0] = 0;
   SibPID_Delta[11][1] = 3;

   SibPID_Delta[12][0] = 0;
   SibPID_Delta[12][1] = 2;

   SibPID_Delta[13][0] = 0;
   SibPID_Delta[13][1] = 1;

   SibPID_Delta[14][0] = 0;
   SibPID_Delta[14][1] = 1;

   SibPID_Delta[15][0] = 0;
   SibPID_Delta[15][1] = 3;

   SibPID_Delta[16][0] = 0;
   SibPID_Delta[16][1] = 2;

   SibPID_Delta[17][0] = 0;
   SibPID_Delta[17][1] = 2;

   SibPID_Delta[18][0] = 0;
   SibPID_Delta[19][0] = 0;
   SibPID_Delta[20][0] = 0;
   SibPID_Delta[21][0] = 0;
   SibPID_Delta[22][0] = 0;
   SibPID_Delta[23][0] = 0;
   SibPID_Delta[24][0] = 0;
   SibPID_Delta[25][0] = 0;

} // FUNCTION : TABLE_GetSibPID_Delta



//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_GetSibPID_Based
// Description :  Get the starting (based) sibling patch indices of the input patch group along
//                different directions
//
// Note        :  1. Can replace Table_02() in Prepare_PatchData()
//                2. Used together with TABLE_GetSibPID_Delta()
//
// Parameter   :  lv           : Target refinement level
//                PID0         : Patch index with LocalID==0 in the target patch group
//                SibPID_Based : Array to store the output based sibling patch indices
//
// Return      :  SibPID_Based
//-------------------------------------------------------------------------------------------------------
void TABLE_GetSibPID_Based( const int lv, const int PID0, int SibPID_Based[] )
{

   SibPID_Based[ 0] = amr->patch[0][lv][PID0  ]->sibling[ 0];
   SibPID_Based[ 1] = amr->patch[0][lv][PID0+1]->sibling[ 1];
   SibPID_Based[ 2] = amr->patch[0][lv][PID0  ]->sibling[ 2];
   SibPID_Based[ 3] = amr->patch[0][lv][PID0+2]->sibling[ 3];
   SibPID_Based[ 4] = amr->patch[0][lv][PID0  ]->sibling[ 4];
   SibPID_Based[ 5] = amr->patch[0][lv][PID0+3]->sibling[ 5];
   SibPID_Based[ 6] = amr->patch[0][lv][PID0  ]->sibling[ 6];
   SibPID_Based[ 7] = amr->patch[0][lv][PID0+1]->sibling[ 7];
   SibPID_Based[ 8] = amr->patch[0][lv][PID0+2]->sibling[ 8];
   SibPID_Based[ 9] = amr->patch[0][lv][PID0+4]->sibling[ 9];
   SibPID_Based[10] = amr->patch[0][lv][PID0  ]->sibling[10];
   SibPID_Based[11] = amr->patch[0][lv][PID0+2]->sibling[11];
   SibPID_Based[12] = amr->patch[0][lv][PID0+3]->sibling[12];
   SibPID_Based[13] = amr->patch[0][lv][PID0+5]->sibling[13];
   SibPID_Based[14] = amr->patch[0][lv][PID0  ]->sibling[14];
   SibPID_Based[15] = amr->patch[0][lv][PID0+3]->sibling[15];
   SibPID_Based[16] = amr->patch[0][lv][PID0+1]->sibling[16];
   SibPID_Based[17] = amr->patch[0][lv][PID0+6]->sibling[17];
   SibPID_Based[18] = amr->patch[0][lv][PID0  ]->sibling[18];
   SibPID_Based[19] = amr->patch[0][lv][PID0+1]->sibling[19];
   SibPID_Based[20] = amr->patch[0][lv][PID0+2]->sibling[20];
   SibPID_Based[21] = amr->patch[0][lv][PID0+4]->sibling[21];
   SibPID_Based[22] = amr->patch[0][lv][PID0+3]->sibling[22];
   SibPID_Based[23] = amr->patch[0][lv][PID0+6]->sibling[23];
   SibPID_Based[24] = amr->patch[0][lv][PID0+5]->sibling[24];
   SibPID_Based[25] = amr->patch[0][lv][PID0+7]->sibling[25];

} // FUNCTION : TABLE_GetSibPID_Based
