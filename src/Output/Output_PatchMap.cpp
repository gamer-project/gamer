#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_PatchMap
// Description :  Output the map of one component of a single patch
//
// Parameter   :  lv      : Target refinement level
//                PID     : Target patch index
//                TSg     : Target Sandglass
//                Comp    : Target component to output
//                          --> no   gravity : [0 ... NCOMP_TOTAL-1]
//                              with gravity : [0 ... NCOMP_TOTAL  ] --> output potential if Comp == NCOMP_TOTAL
//                comment : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_PatchMap( const int lv, const int PID, const int TSg, const int Comp, const char *comment )
{

// check
#  ifdef GRAVITY
   if ( Comp > NCOMP_TOTAL   ||  Comp < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Comp", Comp );
#  else
   if ( Comp > NCOMP_TOTAL-1 ||  Comp < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Comp", Comp );
#  endif

   if ( amr->patch[0][lv][PID] == NULL )
   {
      Aux_Message( stderr, "WARNING : lv %d, PID %d does NOT exist !!\n", lv, PID );
      return;
   }

   if ( Comp < NCOMP_TOTAL  &&  amr->patch[TSg][lv][PID]->fluid == NULL )
   {
      Aux_Message( stderr, "WARNING : lv %d, PID %d, Sg %d does NOT have fluid data !!\n", lv, PID, TSg );
      return;
   }

#  ifdef GRAVITY
   if ( Comp == NCOMP_TOTAL  &&  amr->patch[TSg][lv][PID]->pot == NULL )
   {
      Aux_Message( stderr, "WARNING : lv %d, PID %d, Sg %d does NOT have potential data !!\n", lv, PID, TSg );
      return;
   }
#  endif


   patch_t *Relation = amr->patch[ 0][lv][PID];
   patch_t *Data     = amr->patch[TSg][lv][PID];

   char FileName[100];
   sprintf( FileName, "PatchMap_r%d_lv%d_p%d_v%d", MPI_Rank, lv, PID, Comp );
   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// output header
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Rank %d  Lv %d  PID %d  Local ID %d  TSg %d  Time %13.7e  Step %ld  Counter %ld\n",
            MPI_Rank, lv, PID, PID%8, TSg, Time[lv], Step, AdvanceCounter[lv] );

   fprintf( File, "Father  %d    Son  %d    Corner  (%10d,%10d,%10d)\n\n", Relation->father, Relation->son,
            Relation->corner[0], Relation->corner[1], Relation->corner[2] );

   fprintf( File, "Sibling, Sibling->Son, and Father->Sibling Lists :\n" );

   int Fa, Sib, FaSib, SibSon;
   for (int S=0; S<26; S++)
   {
      Fa     = Relation->father;
      Sib    = Relation->sibling[S];
      FaSib  = ( Fa == -1 ) ? -1 : ( amr->patch[0][lv-1][Fa] != NULL ) ?
                                     amr->patch[0][lv-1][Fa]->sibling[S] : -999;
      SibSon = ( Sib < 0 )  ? Sib : amr->patch[0][lv][Sib]->son;

      fprintf( File, "Sib[%2d] = %6d     Sib_Son = %6d     Fa_Sib[%2d] = %6d\n",
               S, Sib, SibSon, S, FaSib );
   }
   fprintf( File, "\n" );


// output data
   if ( Comp < NCOMP_TOTAL ) // output fluid data
   {
      for (int k=0; k<PATCH_SIZE; k++)
      {
         fprintf( File, "k = %d\n\n", k );

         for (int j=PATCH_SIZE-1; j>=0; j--)
         {
            for (int i=0; i<PATCH_SIZE; i++)    fprintf( File, "%12.5e  ", Data->fluid[Comp][k][j][i] );

            fprintf( File, "\n" );
         }

         fprintf( File, "\n\n\n" );
      }
   }

#  ifdef GRAVITY
   else // Comp == NCOMP_TOTAL --> output potential
   {
      for (int k=0; k<PATCH_SIZE; k++)
      {
         fprintf( File, "k = %d\n\n", k );

         for (int j=PATCH_SIZE-1; j>=0; j--)
         {
            for (int i=0; i<PATCH_SIZE; i++)    fprintf( File, "%12.5e  ", Data->pot[k][j][i] );

            fprintf( File, "\n" );
         }

         fprintf( File, "\n\n\n" );
      }
   } // if ( Comp < NCOMP_TOTAL ) ... else ...
#  endif


   fclose( File );

} // FUNCTION : Output_PatchMap
