#include "ExtractProfile.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_Patch
// Description :  Output the data of a single patch
//
// Example     :  char comment[10];
//                sprintf( comment, "Step%d", AdvanceCounter[6] );
//                Output_Patch( 6, 5560, comment );
//
// Parameter   :  lv       : Targeted refinement level
//                PID      : Targeted patch index
//                FluSg    : Sandglass of the fluid data
//                PotSg    : Sandglass of the potential data
//                comment  : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_Patch( const int lv, const int PID, const char *comment )
{

   const real dh_min = amr.dh[NLEVEL-1];

   patch_t *Relation = amr.patch[lv][PID];
   patch_t *FluData  = amr.patch[lv][PID];
   patch_t *PotData  = amr.patch[lv][PID];

   char FileName[100];
   sprintf( FileName, "Patch_r%d_lv%d_p%d", 0, lv, PID );
   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }


// output patch information
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Rank %d  Lv %d  PID %d  Local ID %d  FluSg %d  PotSg %d  Time %13.7e  Step %ld\n",
            0, lv, PID, PID%8, 0, 0, Time[lv], Step );

   fprintf( File, "Father %d  Son %d  Corner (%4d,%4d,%4d)  Size %13.7e", Relation->father, Relation->son,
            Relation->corner[0], Relation->corner[1], Relation->corner[2], PATCH_SIZE*amr.dh[lv] );
   fprintf( File, "\n" );
   fprintf( File, "EdgeL = (%20.14e, %20.14e, %20.14e)\n", Relation->corner[0]*dh_min,
                                                           Relation->corner[1]*dh_min,
                                                           Relation->corner[2]*dh_min );
   fprintf( File, "EdgeR = (%20.14e, %20.14e, %20.14e)\n", ( Relation->corner[0] + PATCH_SIZE*amr.scale[lv] )*dh_min,
                                                           ( Relation->corner[1] + PATCH_SIZE*amr.scale[lv] )*dh_min,
                                                           ( Relation->corner[2] + PATCH_SIZE*amr.scale[lv] )*dh_min );


   fprintf( File, "\nSibling, Sibling->Son, and Father->Sibling Lists :\n" );

   int Sib, FaSib, SibSon, Fa;
   for (int S=0; S<26; S++)
   {
      Fa     = Relation->father;
      Sib    = Relation->sibling[S];
      FaSib  = ( Fa == -1 ) ? -1 : ( amr.patch[lv-1][Fa] != NULL ) ?
                                     amr.patch[lv-1][Fa]->sibling[S] : -999;
      SibSon = ( Sib < 0 )  ? Sib : amr.patch[lv][Sib]->son;

      fprintf( File, "Sib[%2d] = %6d     Sib_Son = %6d     Fa_Sib[%2d] = %6d\n",
               S, Sib, SibSon, S, FaSib );
   }
   fprintf( File, "\n" );


// output header
   fprintf( File, "(%2s,%2s,%2s)", "i", "j", "k" );

#  if   ( MODEL == HYDRO )
   fprintf( File, "%14s%14s%14s%14s%14s%14s", "Density", "Px", "Py", "Pz", "Energy", "Pressure" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   fprintf( File, "%14s%14s%14s", "Density", "Real", "Imag" );

#  else
#  warning : WARNING : DO YOU WANT TO ADD the FILE HEADER HERE FOR THE NEW MODEL ??
#  endif // MODEL

   for (int v=0; v<NCOMP_PASSIVE; v++)
   fprintf( File, "%12s%02d", "Passive", v );

   fprintf( File, "%14s", "Potential" );

   fprintf( File, "\n" );


// output data
   real u[NCOMP_FLUID];

   for (int k=0; k<PATCH_SIZE; k++)
   for (int j=0; j<PATCH_SIZE; j++)
   for (int i=0; i<PATCH_SIZE; i++)
   {
//    output cell indices
      fprintf( File, "(%2d,%2d,%2d)", i, j, k );

      if ( FluData->fluid != NULL )
      {
//       output all variables in the fluid array
         for (int v=0; v<NCOMP_FLUID; v++)
         {
            u[v] = FluData->fluid[v][k][j][i];
            fprintf( File, " %13.6e", u[v] );
         }

//       output pressure in HYDRO
#        if   ( MODEL == HYDRO )
         fprintf( File, " %13.6e", ( u[ENGY]-0.5*(u[MOMX]*u[MOMX]+u[MOMY]*u[MOMY]+u[MOMZ]*u[MOMZ])/u[DENS] )*
                                   (GAMMA-1.0) );
#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!
#        endif // MODEL

//       output the passive variables
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         fprintf( File, " %13.6e", FluData->fluid[v][k][j][i] );
      } // if ( FluData->fluid != NULL )

      else
      {
//       output empty strings if the fluid array is not allocated
         for (int v=0; v<NCOMP_FLUID; v++)   fprintf( File, " %13s", "" );

#        if   ( MODEL == HYDRO )
         fprintf( File, " %13s", "" );
#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!
#        endif // MODEL

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         fprintf( File, " %13s", "" );
      } // if ( FluData->fluid != NULL ) ... else ...

//    output potential
      if ( PotData->pot != NULL )   fprintf( File, " %13.6e",  PotData->pot[k][j][i]);

      fprintf( File, "\n" );
   } // i,j,k

   fclose( File );

} // FUNCTION : Output_Patch

