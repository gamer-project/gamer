#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_Patch
// Description :  Output the data of a single patch
//
// Example     :  char comment[10];
//                sprintf( comment, "Step%d", AdvanceCounter[6] );
//                Output_Patch( 6, 5560, comment );
//
// Parameter   :  lv      : Target refinement level
//                PID     : Target patch index
//                FluSg   : Sandglass of the fluid data
//                PotSg   : Sandglass of the potential data
//                comment : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_Patch( const int lv, const int PID, const int FluSg, const int PotSg, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( PID < 0  ||  PID >= MAX_PATCH )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (MAX_PATCH = %d) !!\n", "PID", PID, MAX_PATCH );

   if ( FluSg < 0  ||  FluSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FluSg", FluSg );

#  ifdef GRAVITY
   if ( PotSg < 0  ||  PotSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PotSg", PotSg );
#  endif

   if ( amr->patch[0][lv][PID] == NULL )
   {
      Aux_Message( stderr, "WARNING : level = %d, PID = %d does NOT exist !!\n", lv, PID );
      return;
   }


   patch_t *Relation = amr->patch[    0][lv][PID];
   patch_t *FluData  = amr->patch[FluSg][lv][PID];
#  ifdef GRAVITY
   patch_t *PotData  = amr->patch[PotSg][lv][PID];
#  endif

   char FileName[100];
   sprintf( FileName, "Patch_r%d_lv%d_p%d", MPI_Rank, lv, PID );
   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// output patch information
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Rank %d  Lv %d  PID %d  Local ID %d  FluSg %d  PotSg %d  Time %13.7e  Step %ld  Counter %ld\n",
            MPI_Rank, lv, PID, PID%8, FluSg, PotSg, Time[lv], Step, AdvanceCounter[lv] );

   fprintf( File, "Father %d  Son %d  Corner (%10d,%10d,%10d)  Size %13.7e", Relation->father, Relation->son,
            Relation->corner[0], Relation->corner[1], Relation->corner[2], PS1*amr->dh[lv] );
#  ifdef LOAD_BALANCE
   fprintf( File, "  LB_Idx %ld  PaddedCr1D %lu", Relation->LB_Idx, Relation->PaddedCr1D );
#  endif
   fprintf( File, "\n" );
   fprintf( File, "EdgeL = (%20.14e, %20.14e, %20.14e)\n", Relation->EdgeL[0], Relation->EdgeL[1], Relation->EdgeL[2] );
   fprintf( File, "EdgeR = (%20.14e, %20.14e, %20.14e)\n", Relation->EdgeR[0], Relation->EdgeR[1], Relation->EdgeR[2] );

#  ifdef PARTICLE
   fprintf( File, "NPar = %5d  ParListSize = %5d\n", Relation->NPar, Relation->ParListSize );
#  endif


   fprintf( File, "\nSibling, Sibling->Son, and Father->Sibling Lists :\n" );

   int Sib, FaSib, SibSon, Fa;
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



// check whether or not the target patch stores physical data
   if ( FluData->fluid == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store fluid data !!\n", lv, PID );
#  ifdef GRAVITY
   if ( PotData->pot == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store potential data !!\n", lv, PID );
#  endif


// output header
   fprintf( File, "(%2s,%2s,%2s)", "i", "j", "k" );

   for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, "%14s", FieldLabel[v] );

#  if   ( MODEL == HYDRO )
   fprintf( File, "%14s", "Pressure" );
#  ifdef DUAL_ENERGY
   fprintf( File, "%14s", "DE-status" );
#  endif

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )

#  else
#  warning : WARNING : DO YOU WANT TO ADD THE FILE HEADER HERE FOR THE NEW MODEL ??
#  endif // MODEL

#  ifdef GRAVITY
   fprintf( File, "%14s", PotLabel );
#  endif

   fprintf( File, "\n" );


// output data
   real u[NCOMP_TOTAL];

   for (int k=0; k<PATCH_SIZE; k++)
   for (int j=0; j<PATCH_SIZE; j++)
   for (int i=0; i<PATCH_SIZE; i++)
   {
//    output cell indices
      fprintf( File, "(%2d,%2d,%2d)", i, j, k );

      if ( FluData->fluid != NULL )
      {
//       output all variables in the fluid array
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            u[v] = FluData->fluid[v][k][j][i];
            fprintf( File, " %13.6e", u[v] );
         }

//       output pressure and dual-energy status in HYDRO
#        if   ( MODEL == HYDRO )
         const bool CheckMinPres_No = false;
         fprintf( File, " %13.6e", Hydro_GetPressure(u[DENS],u[MOMX],u[MOMY],u[MOMZ],u[ENGY],GAMMA-1.0,CheckMinPres_No,NULL_REAL) );
#        ifdef DUAL_ENERGY
         fprintf( File, " %13c", Relation->de_status[k][j][i] );
#        endif

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )

#        else
#        warning : WARNING : DO YOU WANT TO ADD THE MODEL-SPECIFIC FIELD HERE FOR THE NEW MODEL ??
#        endif // MODEL
      } // if ( FluData->fluid != NULL )

      else
      {
//       output empty strings if the fluid array is not allocated
         for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %13s", "" );

#        if   ( MODEL == HYDRO )
         fprintf( File, " %13s", "" );
#        ifdef DUAL_ENERGY
         fprintf( File, " %13s", "" );
#        endif

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )

#        else
#        warning : WARNING : DO YOU WANT TO ADD THE MODEL-SPECIFIC FIELD HERE FOR THE NEW MODEL ??
#        endif // MODEL
      } // if ( FluData->fluid != NULL ) ... else ...

//    output potential
#     ifdef GRAVITY
      if ( PotData->pot != NULL )   fprintf( File, " %13.6e",  PotData->pot[k][j][i]);
#     endif

      fprintf( File, "\n" );
   } // i,j,k


// output the particle data
#  ifdef PARTICLE
   long ParID;
   fprintf( File, "\n" );
   fprintf( File, "===================\n" );
   fprintf( File, "== PARTICLE DATA == \n" );
   fprintf( File, "===================\n" );
   fprintf( File, "\n" );
   fprintf( File, "%5s  %10s", "No.", "ParID" );
   for (int v=0; v<PAR_NATT_TOTAL; v++)   fprintf( File, "  %13s", ParAttLabel[v] );
   fprintf( File, "\n" );

   for (int p=0; p<Relation->NPar; p++)
   {
      ParID = Relation->ParList[p];

      fprintf( File, "%5d  %10ld", p, ParID );
      for (int v=0; v<PAR_NATT_TOTAL; v++)   fprintf( File, "  %13.6e", amr->Par->Attribute[v][ParID] );

      fprintf( File, "\n" );
   }
#  endif // #ifdef PARTICLE


   fclose( File );

} // FUNCTION : Output_Patch

