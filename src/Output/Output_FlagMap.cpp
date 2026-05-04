#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_FlagMap
// Description :  Output the flag result of entire domain at level "lv"
//
// Note        :  ( X  .  O )  <-->  (not flagged  no patch  flagged)
//
// xyz         :  ( 0, 1, 2 )  -->  output (yz, xz, xy) slice
//-------------------------------------------------------------------------------------------------------
void Output_FlagMap( const int lv, const int xyz, const char *comment )
{

// set up parameters for outputing different slice
   const int NPatch1D[3] = { NX0[0]/PATCH_SIZE*(1<<lv)+4, NX0[1]/PATCH_SIZE*(1<<lv)+4,
                             NX0[2]/PATCH_SIZE*(1<<lv)+4 };
   const char *Dim[2]    = { NULL, NULL };
   int *ii = NULL, *jj = NULL, *kk = NULL;
   int i, j, k, i_end, j_end, k_end;

   i_end = j_end = k_end = 0;

   switch ( xyz )
   {
      case 0:
         Dim[0] = "X";  Dim[1] = "YZ";
         ii = &k;    jj = &i;    kk = &j;
         i_end = NPatch1D[1];
         j_end = NPatch1D[2];
         k_end = NPatch1D[0];
         break;

      case 1:
         Dim[0] = "Y";  Dim[1] = "XZ";
         ii = &i;    jj = &k;    kk = &j;
         i_end = NPatch1D[0];
         j_end = NPatch1D[2];
         k_end = NPatch1D[1];
         break;

      case 2:
         Dim[0] = "Z";  Dim[1] = "XY";
         ii = &i;    jj = &j;    kk = &k;
         i_end = NPatch1D[0];
         j_end = NPatch1D[1];
         k_end = NPatch1D[2];
         break;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "xyz", xyz );
   }


   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/FlagMap_%d_%d_%2s", OUTPUT_DIR, MPI_Rank, lv, Dim[1] );
   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// allocate the FlagMap
   char *FlagMap = new char [ NPatch1D[0]*NPatch1D[1]*NPatch1D[2] ];


// initialize the FlagMap as ''
   for (int P=0; P<NPatch1D[0]*NPatch1D[1]*NPatch1D[2]; P++)      FlagMap[P] = '.';


// record the flag result in the FlagMap
   const int scale0 = amr->scale[ 0];
   const int scale  = amr->scale[lv];
   int ip, jp, kp, ID;

   for (int P=0; P<amr->num[lv]; P++)
   {
      ip = ( amr->patch[0][lv][P]->corner[0]-MPI_Rank_X[0]*NX0[0]*scale0 ) / ( PATCH_SIZE*scale ) + 2;
      jp = ( amr->patch[0][lv][P]->corner[1]-MPI_Rank_X[1]*NX0[1]*scale0 ) / ( PATCH_SIZE*scale ) + 2;
      kp = ( amr->patch[0][lv][P]->corner[2]-MPI_Rank_X[2]*NX0[2]*scale0 ) / ( PATCH_SIZE*scale ) + 2;

      ID = kp*NPatch1D[1]*NPatch1D[0] + jp*NPatch1D[0] + ip;

//    the ID should never repeat
      if ( FlagMap[ID] != '.' )
         Aux_Message( stderr, "WARNING : repeated ID (Rank = %d, P = %d, (ip, jp, kp) = (%d,%d,%d)) !!\n",
                      MPI_Rank, P, ip, jp, kp );

      if ( amr->patch[0][lv][P]->son != -1 )    FlagMap[ID] = 'O';
      else                                      FlagMap[ID] = 'X';
   }


// output the FlagMap slice by slice
   FILE *File = fopen( FileName, "w" );
   fprintf( File, "Time = %13.7e  Step = %ld  Rank = %d  Level = %d\n\n", Time[0], Step, MPI_Rank, lv );

   for (k=0; k<k_end; k++)
   {
      fprintf( File, "%s = %d\n\n", Dim[0], k );

      for (j=j_end-1; j>=0; j--)
      {
         for (i=0; i<i_end; i++)
         {
            ID = (*kk)*NPatch1D[1]*NPatch1D[0] + (*jj)*NPatch1D[0] + (*ii);

            fprintf( File, "  %c", FlagMap[ID] );
         }

         fprintf( File, "\n" );
      }

      fprintf( File, "\n\n\n" );
   }

   fclose( File );


   delete [] FlagMap;

} // FUNCTION : Output_FlagMap
