#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  SiblingSearch
// Description :  Construct the sibling relation for level "lv"
//
// Parameter   :  Target refinement level
//-------------------------------------------------------------------------------------------------------
void SiblingSearch( const int lv )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


// lv == 0 : invoke the function "SiblingSearch_Base" for the root level
   if ( lv == 0 )
   {
      SiblingSearch_Base();
      return;
   }


// lv > 0 :
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->num[lv]; PID++)
   {
      int FaSib, FaSibSon;
      patch_t *fa = amr->patch[0][lv-1][ amr->patch[0][lv][PID]->father ];

//    initialize all siblings as -1
      for (int s=0; s<26; s++)      amr->patch[0][lv][PID]->sibling[s] = -1;

      switch ( PID%8 )
      {
         case 0:
            amr->patch[0][lv][PID]->sibling[ 1] = PID+1;
            amr->patch[0][lv][PID]->sibling[ 3] = PID+2;
            amr->patch[0][lv][PID]->sibling[ 5] = PID+3;
            amr->patch[0][lv][PID]->sibling[ 9] = PID+4;
            amr->patch[0][lv][PID]->sibling[13] = PID+5;
            amr->patch[0][lv][PID]->sibling[17] = PID+6;
            amr->patch[0][lv][PID]->sibling[25] = PID+7;

            if (  ( FaSib = fa->sibling[0] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[15] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSib;
                  amr->patch[0][lv][PID]->sibling[15] = FaSib;
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;
            }

            if (  ( FaSib = fa->sibling[2] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[12] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSib;
                  amr->patch[0][lv][PID]->sibling[12] = FaSib;
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;
            }

            if (  ( FaSib = fa->sibling[4] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[11] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[16] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSib;
                  amr->patch[0][lv][PID]->sibling[11] = FaSib;
                  amr->patch[0][lv][PID]->sibling[16] = FaSib;
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;
            }

            if (  ( FaSib = fa->sibling[6] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSib;
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;
            }

            if (  ( FaSib = fa->sibling[10] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[10] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[10] = FaSib;
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;
            }

            if (  ( FaSib = fa->sibling[14] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[14] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[14] = FaSib;
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;
            }

            if (  ( FaSib = fa->sibling[18] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 7;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;

            break;


         case 1:
            amr->patch[0][lv][PID]->sibling[ 8] = PID+1;
            amr->patch[0][lv][PID]->sibling[15] = PID+2;
            amr->patch[0][lv][PID]->sibling[ 3] = PID+3;
            amr->patch[0][lv][PID]->sibling[24] = PID+4;
            amr->patch[0][lv][PID]->sibling[ 5] = PID+5;
            amr->patch[0][lv][PID]->sibling[13] = PID+6;
            amr->patch[0][lv][PID]->sibling[ 0] = PID-1;

            if (  ( FaSib = fa->sibling[1] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[17] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSib;
                  amr->patch[0][lv][PID]->sibling[17] = FaSib;
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;
            }

            if (  ( FaSib = fa->sibling[2] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[12] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSib;
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;
                  amr->patch[0][lv][PID]->sibling[12] = FaSib;
            }

            if (  ( FaSib = fa->sibling[4] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[14] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[11] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[14] = FaSib;
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSib;
                  amr->patch[0][lv][PID]->sibling[11] = FaSib;
            }

            if (  ( FaSib = fa->sibling[7] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSib;
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;
            }

            if (  ( FaSib = fa->sibling[10] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[10] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;
                  amr->patch[0][lv][PID]->sibling[10] = FaSib;
            }

            if (  ( FaSib = fa->sibling[16] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[16] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[16] = FaSib;
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;
            }

            if (  ( FaSib = fa->sibling[19] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon + 5;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;

           break;


         case 2:
            amr->patch[0][lv][PID]->sibling[12] = PID+1;
            amr->patch[0][lv][PID]->sibling[ 1] = PID+2;
            amr->patch[0][lv][PID]->sibling[ 5] = PID+3;
            amr->patch[0][lv][PID]->sibling[23] = PID+4;
            amr->patch[0][lv][PID]->sibling[17] = PID+5;
            amr->patch[0][lv][PID]->sibling[ 7] = PID-1;
            amr->patch[0][lv][PID]->sibling[ 2] = PID-2;

            if (  ( FaSib = fa->sibling[0] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[15] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSib;
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;
                  amr->patch[0][lv][PID]->sibling[15] = FaSib;
            }

            if (  ( FaSib = fa->sibling[3] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[13] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSib;
                  amr->patch[0][lv][PID]->sibling[13] = FaSib;
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;
            }

            if (  ( FaSib = fa->sibling[4] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[10] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[16] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[10] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSib;
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;
                  amr->patch[0][lv][PID]->sibling[16] = FaSib;
            }

            if (  ( FaSib = fa->sibling[8] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSib;
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;
            }

            if (  ( FaSib = fa->sibling[11] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[11] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[11] = FaSib;
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;
            }

            if (  ( FaSib = fa->sibling[14] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[14] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;
                  amr->patch[0][lv][PID]->sibling[14] = FaSib;
            }

            if (  ( FaSib = fa->sibling[20] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon + 6;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;

           break;


         case 3:
            amr->patch[0][lv][PID]->sibling[21] = PID+1;
            amr->patch[0][lv][PID]->sibling[ 3] = PID+2;
            amr->patch[0][lv][PID]->sibling[ 1] = PID+3;
            amr->patch[0][lv][PID]->sibling[ 9] = PID+4;
            amr->patch[0][lv][PID]->sibling[11] = PID-1;
            amr->patch[0][lv][PID]->sibling[16] = PID-2;
            amr->patch[0][lv][PID]->sibling[ 4] = PID-3;

            if (  ( FaSib = fa->sibling[0] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[14] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[14] = FaSib;
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSib;
            }

            if (  ( FaSib = fa->sibling[2] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[10] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[10] = FaSib;
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSib;
            }

            if (  ( FaSib = fa->sibling[5] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[17] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[13] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSib;
                  amr->patch[0][lv][PID]->sibling[17] = FaSib;
                  amr->patch[0][lv][PID]->sibling[13] = FaSib;
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;
            }

            if (  ( FaSib = fa->sibling[6] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSib;
            }

            if (  ( FaSib = fa->sibling[12] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[12] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[12] = FaSib;
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;
            }

            if (  ( FaSib = fa->sibling[15] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[15] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[15] = FaSib;
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;
            }

            if (  ( FaSib = fa->sibling[22] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon + 4;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;

           break;


         case 4:
            amr->patch[0][lv][PID]->sibling[15] = PID+1;
            amr->patch[0][lv][PID]->sibling[12] = PID+2;
            amr->patch[0][lv][PID]->sibling[ 5] = PID+3;
            amr->patch[0][lv][PID]->sibling[22] = PID-1;
            amr->patch[0][lv][PID]->sibling[ 0] = PID-2;
            amr->patch[0][lv][PID]->sibling[ 2] = PID-3;
            amr->patch[0][lv][PID]->sibling[ 6] = PID-4;

            if (  ( FaSib = fa->sibling[1] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[17] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSib;
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;
                  amr->patch[0][lv][PID]->sibling[17] = FaSib;
            }

            if (  ( FaSib = fa->sibling[3] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[13] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSib;
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;
                  amr->patch[0][lv][PID]->sibling[13] = FaSib;
            }

            if (  ( FaSib = fa->sibling[4] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[14] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[10] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;
                  amr->patch[0][lv][PID]->sibling[14] = FaSib;
                  amr->patch[0][lv][PID]->sibling[10] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 4] = FaSib;
            }

            if (  ( FaSib = fa->sibling[9] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon + 3;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSib;
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;
            }

            if (  ( FaSib = fa->sibling[11] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[11] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;
                  amr->patch[0][lv][PID]->sibling[11] = FaSib;
            }

            if (  ( FaSib = fa->sibling[16] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[16] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;
                  amr->patch[0][lv][PID]->sibling[16] = FaSib;
            }

            if (  ( FaSib = fa->sibling[21] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon + 3;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;

           break;


         case 5:
            amr->patch[0][lv][PID]->sibling[ 7] = PID+1;
            amr->patch[0][lv][PID]->sibling[ 1] = PID+2;
            amr->patch[0][lv][PID]->sibling[16] = PID-1;
            amr->patch[0][lv][PID]->sibling[ 2] = PID-2;
            amr->patch[0][lv][PID]->sibling[ 4] = PID-3;
            amr->patch[0][lv][PID]->sibling[19] = PID-4;
            amr->patch[0][lv][PID]->sibling[10] = PID-5;

            if (  ( FaSib = fa->sibling[0] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[14] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSibSon + 6;
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;
                  amr->patch[0][lv][PID]->sibling[14] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 0] = FaSib;
            }

            if (  ( FaSib = fa->sibling[3] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[11] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[11] = FaSib;
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSib;
            }

            if (  ( FaSib = fa->sibling[5] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[12] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[17] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[12] = FaSib;
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSib;
                  amr->patch[0][lv][PID]->sibling[17] = FaSib;
            }

            if (  ( FaSib = fa->sibling[8] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSib;
            }

            if (  ( FaSib = fa->sibling[13] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[13] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon + 1;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[13] = FaSib;
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;
            }

            if (  ( FaSib = fa->sibling[15] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[15] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;
                  amr->patch[0][lv][PID]->sibling[15] = FaSib;
            }

            if (  ( FaSib = fa->sibling[24] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon + 1;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;

           break;


         case 6:
            amr->patch[0][lv][PID]->sibling[ 3] = PID+1;
            amr->patch[0][lv][PID]->sibling[ 8] = PID-1;
            amr->patch[0][lv][PID]->sibling[11] = PID-2;
            amr->patch[0][lv][PID]->sibling[ 0] = PID-3;
            amr->patch[0][lv][PID]->sibling[20] = PID-4;
            amr->patch[0][lv][PID]->sibling[ 4] = PID-5;
            amr->patch[0][lv][PID]->sibling[14] = PID-6;

            if (  ( FaSib = fa->sibling[1] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[16] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[16] = FaSib;
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSib;
            }

            if (  ( FaSib = fa->sibling[2] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[18] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[10] = FaSibSon + 4;
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSibSon + 5;
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSibSon + 7;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[18] = FaSib;
                  amr->patch[0][lv][PID]->sibling[10] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 6] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 2] = FaSib;
            }

            if (  ( FaSib = fa->sibling[5] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[15] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[13] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[15] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSib;
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;
                  amr->patch[0][lv][PID]->sibling[13] = FaSib;
            }

            if (  ( FaSib = fa->sibling[7] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSib;
            }

            if (  ( FaSib = fa->sibling[12] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[12] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;
                  amr->patch[0][lv][PID]->sibling[12] = FaSib;
            }

            if (  ( FaSib = fa->sibling[17] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[17] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon + 2;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[17] = FaSib;
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;
            }

            if (  ( FaSib = fa->sibling[23] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon + 2;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;

           break;


         case 7:
            amr->patch[0][lv][PID]->sibling[ 2] = PID-1;
            amr->patch[0][lv][PID]->sibling[ 0] = PID-2;
            amr->patch[0][lv][PID]->sibling[ 4] = PID-3;
            amr->patch[0][lv][PID]->sibling[ 6] = PID-4;
            amr->patch[0][lv][PID]->sibling[14] = PID-5;
            amr->patch[0][lv][PID]->sibling[10] = PID-6;
            amr->patch[0][lv][PID]->sibling[18] = PID-7;

            if (  ( FaSib = fa->sibling[1] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[19] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[16] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSibSon + 5;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[19] = FaSib;
                  amr->patch[0][lv][PID]->sibling[16] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 7] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 1] = FaSib;
            }

            if (  ( FaSib = fa->sibling[3] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[20] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[11] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSibSon + 3;
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSibSon + 6;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[20] = FaSib;
                  amr->patch[0][lv][PID]->sibling[11] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 8] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 3] = FaSib;
            }

            if (  ( FaSib = fa->sibling[5] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[22] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[12] = FaSibSon + 1;
                  amr->patch[0][lv][PID]->sibling[15] = FaSibSon + 2;
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSibSon + 4;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[22] = FaSib;
                  amr->patch[0][lv][PID]->sibling[12] = FaSib;
                  amr->patch[0][lv][PID]->sibling[15] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 5] = FaSib;
            }

            if (  ( FaSib = fa->sibling[9] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[21] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSibSon + 3;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[21] = FaSib;
                  amr->patch[0][lv][PID]->sibling[ 9] = FaSib;
            }

            if (  ( FaSib = fa->sibling[13] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[24] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[13] = FaSibSon + 1;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[24] = FaSib;
                  amr->patch[0][lv][PID]->sibling[13] = FaSib;
            }

            if (  ( FaSib = fa->sibling[17] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
               {
                  amr->patch[0][lv][PID]->sibling[23] = FaSibSon;
                  amr->patch[0][lv][PID]->sibling[17] = FaSibSon + 2;
               }
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
            {
                  amr->patch[0][lv][PID]->sibling[23] = FaSib;
                  amr->patch[0][lv][PID]->sibling[17] = FaSib;
            }

            if (  ( FaSib = fa->sibling[25] ) >= 0  )
            {
               if (  ( FaSibSon = amr->patch[0][lv-1][FaSib]->son ) != -1  )
                  amr->patch[0][lv][PID]->sibling[25] = FaSibSon;
            }
            else if ( FaSib <= SIB_OFFSET_NONPERIODIC )
                  amr->patch[0][lv][PID]->sibling[25] = FaSib;

           break;

      }  // switch ( PID%8 )
   }  // for (int PID=0; PID<amr->num[lv]; PID++)

} // FUNCTION : SiblingSearch
