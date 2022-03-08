#include "ExtractUniform.h"

static int Table_01( const int SibID, const int Count );




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_RecordExchangeDataPatchID
// Description :  Record the information of patches for sending and receiving data between neighbor ranks
//                in the variable "ParaVar"
//
// Parameter   :  lv : The targeted refinement level
//-------------------------------------------------------------------------------------------------------
void Buf_RecordExchangeDataPatchID( const int lv )
{

// check
   if ( lv < 0  ||  lv > NLEVEL-1 )
   {
      fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "lv", lv );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      MPI_Exit();
   }


// MirrorSib : the mirror-symmetric sibling index
   const int MirrorSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const int scale0        = amr.scale[ 0];
   const int scale         = amr.scale[lv];

   int Sib, NSide=0, NFace, ListLength[3], Disp[3], ijk[3], FaceList[4], SibList_Send[9], SibList_Recv[9];
   int NBuf, PID, SibPID, Pos;
   int *Corner=NULL;
   int *RecvP_PosList=NULL;


// begin the main loop of Buf_RecordExchangePatchID
   for (int s=0; s<26; s++)
   {

//    initialize counters as zero
      ParaVar.SendP_NList[lv][s] = 0;
      ParaVar.RecvP_NList[lv][s] = 0;


//    deallocate memory
      if ( ParaVar.SendP_IDList[lv][s] != NULL )
      {
         delete [] ParaVar.SendP_IDList[lv][s];
         ParaVar.SendP_IDList[lv][s] = NULL;
      }

      if ( ParaVar.RecvP_IDList[lv][s] != NULL )
      {
         delete [] ParaVar.RecvP_IDList[lv][s];
         ParaVar.RecvP_IDList[lv][s] = NULL;
      }


//    nothing to do if there are no boundary patches at the targeted direction
      if ( ParaVar.BounP_NList[lv][s] == 0 )    continue;


//    allocate the maximum necessary memory
      NBuf = ( NPatchComma[lv][s+2] - NPatchComma[lv][s+1] ) / TABLE_05( s );

      ParaVar.SendP_IDList [lv][s] = new int [ ParaVar.BounP_NList[lv][s] ];
      ParaVar.RecvP_IDList [lv][s] = new int [ NBuf ];
      RecvP_PosList                = new int [ NBuf ];


//    set up the NSide
      switch ( s )
      {
         case 0: case 1: case 2: case 3: case 4: case 5:
            NSide = 9;
            break;

         case 6: case 7: case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17:
            NSide = 3;
            break;

         case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25:
            NSide = 1;
            break;
      }


//    set up the SibList_Send and SibList_Recv
      for (int Side=0; Side<NSide; Side++)
      {
         SibList_Send[Side] = Table_01(            s, Side );
         SibList_Recv[Side] = Table_01( MirrorSib[s], Side );
      }


//    set up the NFace and FaceList
      NFace = TABLE_04(s);
      for (int f=0; f<NFace; f++)   FaceList[f] = TABLE_07( MirrorSib[s], f );


//    set up the length and displacement of loops in different directions
      for (int d=0; d<3; d++)
      {
         ListLength[d] = TABLE_01( s, 'x'+d, 1, NX0[d]/PATCH_SIZE*(1<<lv), 1 );
         Disp      [d] = TABLE_01( s, 'x'+d, 1, 0, -NX0[d]/PATCH_SIZE*(1<<lv) );
      }


//    a. set up the SendP_IDList
//    --------------------------
      for (int ID=0; ID<ParaVar.BounP_NList[lv][s]; ID++)
      {
         PID = ParaVar.BounP_IDList[lv][s][ID];


         if ( amr.patch[lv][PID]->fluid == NULL )  continue;


         if ( PID < NPatchComma[lv][1] )
         {
            for (int Side=0; Side<NSide; Side++)
            {
               Sib    = SibList_Send[Side];
               SibPID = amr.patch[lv][PID]->sibling[Sib];

               if ( SibPID >= 0  &&  SibPID >= NPatchComma[lv][s+1]  &&  SibPID < NPatchComma[lv][s+2] )
               {
                  ParaVar.SendP_IDList[lv][s][ ParaVar.SendP_NList[lv][s] ] = PID;
                  ParaVar.SendP_NList [lv][s] ++;

                  break;
               }
            }
         }
      } // for (int ID=0; ID<ParaVar.BounP_NList[lv][s]; ID++)


//    b. set up the RecvP_IDList
//    --------------------------
      for (int PID0=NPatchComma[lv][s+1]; PID0<NPatchComma[lv][s+2]; PID0+=8)
      for (int Face=0; Face<NFace; Face++)
      {
         PID = PID0 + FaceList[Face];


         if ( amr.patch[lv][PID]->fluid == NULL )  continue;


         for (int Side=0; Side<NSide; Side++)
         {
            Sib    = SibList_Recv[Side];
            SibPID = amr.patch[lv][PID]->sibling[Sib];

            if ( SibPID >= 0  &&  SibPID < NPatchComma[lv][1] )
            {
               Corner = amr.patch[lv][PID]->corner;

               for (int d=0; d<3; d++)    ijk[d] = ( Corner[d]-MyRank_X[d]*NX0[d]*scale0 ) / (PATCH_SIZE*scale) + Disp[d];

               Pos = ijk[2]*ListLength[1]*ListLength[0] + ijk[1]*ListLength[0] + ijk[0];

               ParaVar.RecvP_IDList[lv][s][ ParaVar.RecvP_NList[lv][s] ] = PID;
               RecvP_PosList              [ ParaVar.RecvP_NList[lv][s] ] = Pos;
               ParaVar.RecvP_NList [lv][s] ++;

               break;
            }
         }
      }


//    c. sort the RecvP_IDList
//    ------------------------
      Buf_SortBoundaryPatch( ParaVar.RecvP_NList[lv][s], ParaVar.RecvP_IDList[lv][s], RecvP_PosList );


//    deallocate memory if no data needed to be sent or received
      if ( ParaVar.SendP_NList[lv][s] == 0 )
      {
         delete [] ParaVar.SendP_IDList[lv][s];
         ParaVar.SendP_IDList [lv][s] = NULL;
      }

      if ( ParaVar.RecvP_NList[lv][s] == 0 )
      {
         delete [] ParaVar.RecvP_IDList[lv][s];
         ParaVar.RecvP_IDList [lv][s] = NULL;
      }

      delete [] RecvP_PosList;

   } // for (int s=0; s<26; s++)

}



// ============
// |  Tables  |
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  Return the sibling index for the function "Buf_RecordExchangePatchID"
//
// Parameter   :  SibID : sibling index (0~25)
//                Count : patch counter (0~3)
//-------------------------------------------------------------------------------------------------------
int Table_01( const int SibID, const int Count )
{

   switch ( SibID )
   {
      case 0:
      {
         switch( Count )
         {
            case 0:     return  0;
            case 1:     return  6;
            case 2:     return  8;
            case 3:     return 14;
            case 4:     return 15;
            case 5:     return 18;
            case 6:     return 20;
            case 7:     return 22;
            case 8:     return 24;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 1:
      {
         switch ( Count )
         {
            case 0:     return  1;
            case 1:     return  7;
            case 2:     return  9;
            case 3:     return 16;
            case 4:     return 17;
            case 5:     return 19;
            case 6:     return 21;
            case 7:     return 23;
            case 8:     return 25;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 2:
      {
         switch ( Count )
         {
            case 0:     return  2;
            case 1:     return  6;
            case 2:     return  7;
            case 3:     return 10;
            case 4:     return 12;
            case 5:     return 18;
            case 6:     return 19;
            case 7:     return 22;
            case 8:     return 23;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 3:
      {
         switch ( Count )
         {
            case 0:     return  3;
            case 1:     return  8;
            case 2:     return  9;
            case 3:     return 11;
            case 4:     return 13;
            case 5:     return 20;
            case 6:     return 21;
            case 7:     return 24;
            case 8:     return 25;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 4:
      {
         switch ( Count )
         {
            case 0:     return  4;
            case 1:     return 10;
            case 2:     return 11;
            case 3:     return 14;
            case 4:     return 16;
            case 5:     return 18;
            case 6:     return 19;
            case 7:     return 20;
            case 8:     return 21;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 5:
      {
         switch ( Count )
         {
            case 0:     return  5;
            case 1:     return 12;
            case 2:     return 13;
            case 3:     return 15;
            case 4:     return 17;
            case 5:     return 22;
            case 6:     return 23;
            case 7:     return 24;
            case 8:     return 25;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 6:
      {
         switch ( Count )
         {
            case 0:     return  6;
            case 1:     return 18;
            case 2:     return 22;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 7:
      {
         switch ( Count )
         {
            case 0:     return  7;
            case 1:     return 19;
            case 2:     return 23;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 8:
      {
         switch ( Count )
         {
            case 0:     return  8;
            case 1:     return 20;
            case 2:     return 24;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 9:
      {
         switch ( Count )
         {
            case 0:     return  9;
            case 1:     return 21;
            case 2:     return 25;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 10:
      {
         switch ( Count )
         {
            case 0:     return 10;
            case 1:     return 18;
            case 2:     return 19;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 11:
      {
         switch ( Count )
         {
            case 0:     return 11;
            case 1:     return 20;
            case 2:     return 21;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 12:
      {
         switch ( Count )
         {
            case 0:     return 12;
            case 1:     return 22;
            case 2:     return 23;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 13:
      {
         switch ( Count )
         {
            case 0:     return 13;
            case 1:     return 24;
            case 2:     return 25;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 14:
      {
         switch ( Count )
         {
            case 0:     return 14;
            case 1:     return 18;
            case 2:     return 20;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 15:
      {
         switch ( Count )
         {
            case 0:     return 15;
            case 1:     return 22;
            case 2:     return 24;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 16:
      {
         switch ( Count )
         {
            case 0:     return 16;
            case 1:     return 19;
            case 2:     return 21;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 17:
      {
         switch ( Count )
         {
            case 0:     return 17;
            case 1:     return 23;
            case 2:     return 25;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }
      }

      case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25:
         return SibID;

      default:
      {
         fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n",  "SibID", SibID );
         fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         MPI_Exit();
         exit(-1);
      }

   } // switch ( SibID )

}
