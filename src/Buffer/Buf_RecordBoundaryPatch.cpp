#include "GAMER.h"

#ifndef SERIAL




//###OPTIMIZATION : use heap sort or quick sort
//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_RecordBoundaryPatch
// Description :  Record the patches near the sub-domain boundaries in "amr->ParaVar->BounP_IDList"
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Buf_RecordBoundaryPatch( const int lv )
{

// invoke the function "Buf_RecordBoundaryPatch_Base" for the base level
   if ( lv == 0 )
   {
      Buf_RecordBoundaryPatch_Base();

      return;
   }


// check
   if ( lv < 0  ||  lv > NLEVEL-1 )    Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


   const int scale0 = amr->scale[ 0];
   const int scale  = amr->scale[lv];

   int  NSide, Table[4], ListLength[3], ListLength1D, Disp[3], ijk[3];
   int  Pos0, Pos, FaPID, PID0;
   int *Corner;


// begin the main loop of Buf_RecordBoundaryPatch
   for (int s=0; s<26; s++)
   {

//    initialize the BounP_NList[lv][s] as 0
      amr->ParaVar->BounP_NList[lv][s] = 0;


//    deallocate memory
      if ( amr->ParaVar->BounP_IDList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->BounP_IDList[lv][s];
         amr->ParaVar->BounP_IDList[lv][s] = NULL;
      }

      if ( amr->ParaVar->BounP_PosList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->BounP_PosList[lv][s];
         amr->ParaVar->BounP_PosList[lv][s] = NULL;
      }


//    nothing to do if there are no patches at the target level and direction
      if ( amr->ParaVar->BounP_NList[lv-1][s] == 0  ||  amr->num[lv] == 0 )   continue;


//    set up "NSide and Table"
      NSide = TABLE_04(s);
      for (int t=0; t<NSide; t++)   Table[t] = TABLE_07( s, t );


//    allocate the maximum necessary memory
      amr->ParaVar->BounP_IDList [lv][s] = new int [ amr->ParaVar->BounP_NList[lv-1][s]*NSide ];
      amr->ParaVar->BounP_PosList[lv][s] = new int [ amr->ParaVar->BounP_NList[lv-1][s]*NSide ];


//    set up "ListLength and ListLength1D"
      for (int d=0; d<3; d++)
      {
         ListLength[d] = TABLE_01( s, 'x'+d, 1, NX0[d]/PATCH_SIZE*(1<<lv)+4, 1 );
         Disp      [d] = TABLE_01( s, 'x'+d, 0, 2, -NX0[d]/PATCH_SIZE*(1<<lv)+2 );
      }

      ListLength1D = ( s < 2 ) ? ListLength[1] : ListLength[0];


//    fill up the arrays "amr->ParaVar->BounP_IDList[lv][s] and amr->ParaVar->BounP_PosList[lv][s]"
      for (int FaID=0; FaID<amr->ParaVar->BounP_NList[lv-1][s]; FaID++)
      {
         FaPID = amr->ParaVar->BounP_IDList[lv-1][s][FaID];

#        ifdef GAMER_DEBUG
         if ( FaPID < 0  ||  FaPID >= amr->num[lv-1] )   Aux_Error( ERROR_INFO, "incorrect FaPID = %d !!\n", FaPID );
#        endif

         PID0  = amr->patch[0][lv-1][FaPID]->son;

//       we have assumed that the outermost buffer patches will NOT have child patches
         if ( PID0 != -1 )
         {
            Corner = amr->patch[0][lv][PID0]->corner;

            for (int d=0; d<3; d++)    ijk[d] = ( Corner[d]-MPI_Rank_X[d]*NX0[d]*scale0 ) / (PATCH_SIZE*scale) + Disp[d];

            Pos0 = ( ijk[2]*ListLength[1] + ijk[1] )*ListLength[0] + ijk[0];

            for (int Side=0; Side<NSide; Side++)
            {
               Pos = Pos0 + Side%2 + (Side/2)*ListLength1D;

               amr->ParaVar->BounP_IDList [lv][s][ amr->ParaVar->BounP_NList[lv][s] ] = PID0 + Table[Side];
               amr->ParaVar->BounP_PosList[lv][s][ amr->ParaVar->BounP_NList[lv][s] ] = Pos;
               amr->ParaVar->BounP_NList  [lv][s] ++;
            }

         } // if ( PID0 != -1 )
      } // for (int FaID=0; FaID<amr->ParaVar->NList[lv-1][s]; FaID++)


//    sort the boundary patches according to their positions
      Buf_SortBoundaryPatch( amr->ParaVar->BounP_NList[lv][s], amr->ParaVar->BounP_IDList[lv][s],
                             amr->ParaVar->BounP_PosList[lv][s] );


//    deallocate memory if no boundary patches are found
      if ( amr->ParaVar->BounP_NList[lv][s] == 0 )
      {
         delete [] amr->ParaVar->BounP_IDList [lv][s];
         delete [] amr->ParaVar->BounP_PosList[lv][s];

         amr->ParaVar->BounP_IDList [lv][s] = NULL;
         amr->ParaVar->BounP_PosList[lv][s] = NULL;
      }

   } // for (int s=0; s<26; s++)

} // FUNCTION : Buf_RecordBoundaryPatch



#endif // #ifndef SERIAL
