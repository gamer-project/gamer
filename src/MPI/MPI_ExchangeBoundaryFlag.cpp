#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_ExchangeBoundaryFlag
// Description :  Get the "BuffFlag_NList" and "BuffFlag_PosList" from 26 neighbor ranks
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void MPI_ExchangeBoundaryFlag( const int lv )
{

   const int v[26] = { 0,1,2,3,4,5,6,9,7,8,10,13,11,12,14,17,16,15,18,25,19,24,20,23,21,22 };
   int SendTarget[2], RecvTarget[2];
   MPI_Request Req[4];


// a. get the BuffFlag_NList
   for (int s=0; s<26; s+=2)
   {
//    properly deal with the non-periodic B.C.
      for (int t=0; t<2; t++)
      {
         SendTarget[t] = ( MPI_SibRank[ v[s+t] ] < 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];
         RecvTarget[t] = ( MPI_SibRank[ v[s+t] ] < 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];
      }

      MPI_Isend( &amr->ParaVar->BounFlag_NList[lv][ v[s  ] ], 1, MPI_INT, SendTarget[0], 0, MPI_COMM_WORLD, &Req[0] );
      MPI_Isend( &amr->ParaVar->BounFlag_NList[lv][ v[s+1] ], 1, MPI_INT, SendTarget[1], 1, MPI_COMM_WORLD, &Req[1] );

      MPI_Irecv( &amr->ParaVar->BuffFlag_NList[lv][ v[s  ] ], 1, MPI_INT, RecvTarget[0], 1, MPI_COMM_WORLD, &Req[2] );
      MPI_Irecv( &amr->ParaVar->BuffFlag_NList[lv][ v[s+1] ], 1, MPI_INT, RecvTarget[1], 0, MPI_COMM_WORLD, &Req[3] );

      MPI_Waitall( 4, Req, MPI_STATUSES_IGNORE );

//    properly deal with the non-periodic B.C.
      for (int t=0; t<2; t++)
         if ( MPI_SibRank[ v[s+t] ] < 0 )    amr->ParaVar->BuffFlag_NList[lv][ v[s+t] ] = 0;

   } // for (int s=0; s<26; s+=2)


// b. allocate memory for "BuffFlag_PosList" (which will be deallocated in the function "Flag_Buffer")
   for (int s=0; s<26; s++)
   {
      if ( amr->ParaVar->BuffFlag_PosList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->BuffFlag_PosList[lv][s];
         amr->ParaVar->BuffFlag_PosList[lv][s] = NULL;
      }

      amr->ParaVar->BuffFlag_PosList[lv][s] = new int [ amr->ParaVar->BuffFlag_NList[lv][s] ];
   }


// c. get the BuffFlag_PosList
   for (int s=0; s<26; s+=2)
   {
      for (int t=0; t<2; t++)
      {
         SendTarget[t] = ( amr->ParaVar->BounFlag_NList[lv][ v[s+t] ] == 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];
         RecvTarget[t] = ( amr->ParaVar->BuffFlag_NList[lv][ v[s+t] ] == 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];

#        ifdef GAMER_DEBUG
         if (  amr->ParaVar->BounFlag_NList[lv][ v[s+t] ] != 0  &&  ( SendTarget[t] < 0 || SendTarget[t] >= MPI_NRank )  )
            Aux_Error( ERROR_INFO, "incorrect SendTarget[%d] = %d !!\n", t, SendTarget[t] );
         if (  amr->ParaVar->BuffFlag_NList[lv][ v[s+t] ] != 0  &&  ( RecvTarget[t] < 0 || RecvTarget[t] >= MPI_NRank )  )
            Aux_Error( ERROR_INFO, "incorrect RecvTarget[%d] = %d !!\n", t, RecvTarget[t] );
#        endif
      }

      MPI_Isend( amr->ParaVar->BounFlag_PosList[lv][ v[s  ] ], amr->ParaVar->BounFlag_NList[lv][ v[s  ] ],
                 MPI_INT, SendTarget[0], 2, MPI_COMM_WORLD, &Req[0] );
      MPI_Isend( amr->ParaVar->BounFlag_PosList[lv][ v[s+1] ], amr->ParaVar->BounFlag_NList[lv][ v[s+1] ],
                 MPI_INT, SendTarget[1], 3, MPI_COMM_WORLD, &Req[1] );

      MPI_Irecv( amr->ParaVar->BuffFlag_PosList[lv][ v[s  ] ], amr->ParaVar->BuffFlag_NList[lv][ v[s  ] ],
                 MPI_INT, RecvTarget[0], 3, MPI_COMM_WORLD, &Req[2] );
      MPI_Irecv( amr->ParaVar->BuffFlag_PosList[lv][ v[s+1] ], amr->ParaVar->BuffFlag_NList[lv][ v[s+1] ],
                 MPI_INT, RecvTarget[1], 2, MPI_COMM_WORLD, &Req[3] );

      MPI_Waitall( 4, Req, MPI_STATUSES_IGNORE );

   } // for (int s=0; s<26; s+=2)


// d. deallocate memory for "BounFlag_PosList"
   for (int s=0; s<26; s++)
   {
      if ( amr->ParaVar->BounFlag_PosList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->BounFlag_PosList[lv][s];
         amr->ParaVar->BounFlag_PosList[lv][s] = NULL;
      }

      amr->ParaVar->BounFlag_NList[lv][s] = 0;
   }

} // FUNCTION : MPI_ExchangeBoundaryFlag



#endif // #ifndef SERIAL
