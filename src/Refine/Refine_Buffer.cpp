#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Refine_Buffer
// Description :  Construct buffer patches at level "lv+1" according to the flagging result of buffer patches
//                at level "lv" 
//
// Parameter   :  lv          : Target level to be refined
//                SonTable    : Table recording the linking index of each buffer father patch to GrandTable
//                GrandTable  : Table recording the patch IDs of grandson buffer patches
//-------------------------------------------------------------------------------------------------------
void Refine_Buffer( const int lv, const int *SonTable, const int *GrandTable )
{

   const int Width = PATCH_SIZE * amr->scale[lv+1];   // scale of a single patch at level "lv+1"
   bool AllocData[8];                                 // allocate data or not
   int *Cr;                                           // corner coordinates


// a. allocate buffer patches
// ------------------------------------------------------------------------------------------------
   for (int s=0; s<26; s++)   
   {
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


      for (int PID=amr->NPatchComma[lv][s+1]; PID<amr->NPatchComma[lv][s+2]; PID++)
      {
#        ifdef GAMER_DEBUG
         if ( MPI_SibRank[s] < 0 )  Aux_Error( ERROR_INFO, "why there are buffer patches !!\n" );
#        endif

         if ( amr->patch[0][lv][PID]->flag )
         {
//          construct relation : father -> child   
            amr->patch[0][lv][PID]->son = amr->num[lv+1];


//          allocate child patches and construct relation : child -> father
            Cr = amr->patch[0][lv][PID]->corner;

            amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2],       PID, AllocData[0], AllocData[0] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2],       PID, AllocData[1], AllocData[1] );
            amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2],       PID, AllocData[2], AllocData[2] );
            amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2]+Width, PID, AllocData[3], AllocData[3] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2],       PID, AllocData[4], AllocData[4] );
            amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2]+Width, PID, AllocData[5], AllocData[5] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2]+Width, PID, AllocData[6], AllocData[6] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, PID, AllocData[7], AllocData[7] );


//          record the number of buffer patches in each sibling direction
            amr->NPatchComma[lv+1][s+2] += 8;

         } // if ( amr->patch[0][lv][PID]->flag )
      } // for (int PID=amr->NPatchComma[lv][s+1]; PID<amr->NPatchComma[lv][s+2]; PID++)

      for (int n=s+3; n<28; n++)    amr->NPatchComma[lv+1][n] = amr->num[lv+1];

   } // for (int s=0; s<26; s++)


// check the amr->NPatchComma recording
   if ( amr->NPatchComma[lv+1][27] != amr->num[lv+1] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][27] (%d) != amr->num[%d] (%d) !!\n",
                 lv+1, amr->NPatchComma[lv+1][27], lv+1, amr->num[lv+1] );


// b. construct relations : son <-> grandson
// ------------------------------------------------------------------------------------------------
   if ( lv < NLEVEL-2 )
   {
      /*
      for (int s=0; s<26; s++)   
      {
         if ( amr->NPatchComma[lv+2][s+1] == amr->NPatchComma[lv+2][s+2] )    continue;

         for (int SonPID0=amr->NPatchComma[lv+1][s+1]; SonPID0<amr->NPatchComma[lv+1][s+2]; SonPID0+=8)
         */
         for (int SonPID0=amr->NPatchComma[lv+1][1]; SonPID0<amr->NPatchComma[lv+1][27]; SonPID0+=8)
         {
            const int FaID      = amr->patch[0][lv+1][SonPID0]->father - amr->NPatchComma[lv][1];
            const int TempSonID = SonTable[FaID];

            if ( TempSonID != -1 )
            {
               for (int t=0; t<8; t++)
               {
                  const int SonPID    = SonPID0 + t;
                  const int GrandPID0 = GrandTable[ TempSonID + t ];
                 
                  if ( GrandPID0 != -1 )
                  {
                     amr->patch[0][lv+1][SonPID]->son = GrandPID0;
                 
                     for (int GrandPID=GrandPID0; GrandPID<GrandPID0+8; GrandPID++) 
                        amr->patch[0][lv+2][GrandPID]->father = SonPID;
                  }
               }
            }

         } // for (int SonPID0=amr->NPatchComma[lv+1][s+1]; SonPID0<amr->NPatchComma[lv+1][s+2]; SonPID0+=8)
//      } // for (int s=0; s<26; s++)
   } // if ( lv < NLEVEL-2 )

} // FUNCTION : Refine_Buffer



#endif // #ifndef SERIAL
