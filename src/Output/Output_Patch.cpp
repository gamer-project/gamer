#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_Patch
// Description :  Output the data of a single patch
//
// Example     :  const int lv  = 2;
//                const int PID = 20;
//                char comment[10];
//                sprintf( comment, "Step%ld", AdvanceCounter[lv] );
//                Output_Patch( lv, PID, amr->FluSg[lv], amr->MagSg[lv], amr->PotSg[lv], comment );
//
// Parameter   :  lv      : Target refinement level
//                PID     : Target patch index
//                FluSg   : Sandglass of the fluid data
//                MagSg   : Sandglass of the magnetic field data
//                PotSg   : Sandglass of the potential data
//                comment : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_Patch( const int lv, const int PID, const int FluSg, const int MagSg, const int PotSg, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( PID < 0  ||  PID >= MAX_PATCH )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (MAX_PATCH = %d) !!\n", "PID", PID, MAX_PATCH );

   if ( FluSg < 0  ||  FluSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FluSg", FluSg );

#  if ( MODEL == MHD )
   if ( MagSg < 0  ||  MagSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "MagSg", MagSg );
#  endif

#  ifdef GRAVITY
   if ( PotSg < 0  ||  PotSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PotSg", PotSg );
#  endif

   if ( amr->patch[0][lv][PID] == NULL )
   {
      Aux_Message( stderr, "WARNING : level = %d, PID = %d does NOT exist !!\n", lv, PID );
      return;
   }


   patch_t *Relation                   = amr->patch[    0][lv][PID];
   real    (*fluid)[PS1][PS1][PS1]     = amr->patch[FluSg][lv][PID]->fluid;
#  if ( MODEL == MHD )
   real    (*magnetic)[PS1_P1*PS1*PS1] = amr->patch[MagSg][lv][PID]->magnetic;
#  endif
#  ifdef GRAVITY
   real    (*pot)[PS1][PS1]            = amr->patch[PotSg][lv][PID]->pot;
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


// 1. output patch information
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
   if ( fluid == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store fluid data !!\n", lv, PID );
#  if ( MODEL == MHD )
   if ( magnetic == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store magnetic field data !!\n", lv, PID );
#  endif
#  ifdef GRAVITY
   if ( pot == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store potential data !!\n", lv, PID );
#  endif



// 2. output physical data
// header
   fprintf( File, "(%2s,%2s,%2s)", "i", "j", "k" );

#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   fprintf( File, "%14s%14s%14s%14s%14s%14s", "Density", "Momentum X", "Momentum Y", "Momentum Z", "Energy", "Pressure" );
#  if ( MODEL == MHD )
   fprintf( File, "%14s%14s%14s%14s", "B_X", "B_Y", "B_Z", "0.5*B^2" );
#  endif
#  ifdef DUAL_ENERGY
   fprintf( File, "%14s", "DE-status" );
#  endif

#  elif ( MODEL == ELBDM )
   fprintf( File, "%14s%14s%14s", "Density", "Real", "Imag" );

#  else
#  warning : WARNING : DO YOU WANT TO ADD the FILE HEADER HERE FOR THE NEW MODEL ??
#  endif // MODEL

   for (int v=0; v<NCOMP_PASSIVE; v++)
   fprintf( File, "%14s", PassiveFieldName_Grid[v] );

#  ifdef GRAVITY
   fprintf( File, "%14s", "Potential" );
#  endif

   fprintf( File, "\n" );


// physical data
   real u[NCOMP_FLUID];

   for (int k=0; k<PATCH_SIZE; k++)
   for (int j=0; j<PATCH_SIZE; j++)
   for (int i=0; i<PATCH_SIZE; i++)
   {
//    cell indices
      fprintf( File, "(%2d,%2d,%2d)", i, j, k );

      if ( fluid != NULL )
      {
//       all variables in the fluid array
         for (int v=0; v<NCOMP_FLUID; v++)
         {
            u[v] = fluid[v][k][j][i];
            fprintf( File, " %13.6e", u[v] );
         }

//       pressure
#        if ( MODEL == HYDRO  ||  MODEL == MHD )

#        if   ( MODEL == HYDRO )
         const real Pres  = CPU_GetPressure( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY],
                                             GAMMA-1.0, false, NULL_REAL, NULL_REAL );
#        elif ( MODEL == MHD )
//       set pressure to NULL_REAL if somehow magnetic[] is not allocated (likely due to a bug)
         const real EngyB = ( magnetic == NULL ) ? NULL_REAL :
                                                   MHD_GetCellCenteredBEnergy( lv, PID, i, j, k );
         const real Pres  = ( magnetic == NULL ) ? NULL_REAL :
                                                   CPU_GetPressure( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY],
                                                                    GAMMA-1.0, false, NULL_REAL, EngyB );
#        endif
         fprintf( File, " %13.6e", Pres );

//       magnetic field
#        if ( MODEL == MHD )
         real B[3] = { NULL_REAL, NULL_REAL, NULL_REAL };
         if ( magnetic != NULL )    MHD_GetCellCenteredBField( B, lv, PID, i, j, k );
         fprintf( File, " %13.6e %13.6e %13.6e %13.6e", B[MAGX], B[MAGY], B[MAGZ], EngyB );
#        endif

//       dual-energy variable
#        ifdef DUAL_ENERGY
         fprintf( File, " %13c", Relation->de_status[k][j][i] );
#        endif
#        endif // if ( MODEL == HYDRO  ||  MODEL == MHD )

//       passive variables
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         fprintf( File, " %13.6e", fluid[v][k][j][i] );
      } // if ( fluid != NULL )

      else
      {
//       output empty strings if the fluid array is not allocated
         for (int v=0; v<NCOMP_FLUID; v++)   fprintf( File, " %13s", "" );

#        if ( MODEL == HYDRO  ||  MODEL == MHD )
         fprintf( File, " %13s", "" );

#        if ( MODEL == MHD )
         fprintf( File, " %13s", "" );
#        endif

#        ifdef DUAL_ENERGY
         fprintf( File, " %13s", "" );
#        endif

#        endif // #if ( MODEL == HYDRO  ||  MODEL == MHD )

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         fprintf( File, " %13s", "" );
      } // if ( fluid != NULL ) ... else ...

//    potential
#     ifdef GRAVITY
      if ( pot != NULL )   fprintf( File, " %13.6e", pot[k][j][i] );
      else                 fprintf( File, " %13s", "" );
#     endif

      fprintf( File, "\n" );
   } // i,j,k



// 3. output face-centered magnetic field
#  if ( MODEL == MHD )
   fprintf( File, "\n" );
   fprintf( File, "====================================\n" );
   fprintf( File, "== MAGNETIC FIELD (face-centered) == \n" );
   fprintf( File, "====================================\n" );
   fprintf( File, "\n" );

   if ( magnetic != NULL )
   {
//    header
      fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s\n", "i", "j", "k", "B_X", "B_Y", "B_Z" );

      for (int k=0; k<PS1_P1; k++)
      for (int j=0; j<PS1_P1; j++)
      for (int i=0; i<PS1_P1; i++)
      {
//       cell indices
         fprintf( File, "(%2d,%2d,%2d)", i, j, k );

//       B_X
         if ( j < PS1  &&  k < PS1 )   fprintf( File, " %13.6e", magnetic[MAGX][ IDX321_BX(i,j,k) ] );
         else                          fprintf( File, " %13s", "" );

//       B_Y
         if ( i < PS1  &&  k < PS1 )   fprintf( File, " %13.6e", magnetic[MAGY][ IDX321_BY(i,j,k) ] );
         else                          fprintf( File, " %13s", "" );

//       B_Z
         if ( i < PS1  &&  j < PS1 )   fprintf( File, " %13.6e", magnetic[MAGZ][ IDX321_BZ(i,j,k) ] );
         else                          fprintf( File, " %13s", "" );

         fprintf( File, "\n" );
      } // i,j,k
   } // if ( magnetic != NULL )
#  endif // #if ( MODEL == MHD )



// 4. output particles
#  ifdef PARTICLE
   long ParID;
   fprintf( File, "\n" );
   fprintf( File, "===================\n" );
   fprintf( File, "== PARTICLE DATA == \n" );
   fprintf( File, "===================\n" );
   fprintf( File, "\n" );
   fprintf( File, "%5s  %10s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s",
            "No.", "ParID", "Mass", "X", "Y", "Z", "Vx", "Vy", "Vz", "Time" );
#  ifdef STORE_PAR_ACC
   fprintf( File, "  %13s  %13s  %13s", "AccX", "AccY", "AccZ" );
#  endif
   for (int v=0; v<PAR_NPASSIVE; v++)
   fprintf( File, "  %13s", PassiveFieldName_Par[v] );
   fprintf( File, "\n" );

   for (int p=0; p<Relation->NPar; p++)
   {
      ParID = Relation->ParList[p];

      fprintf( File, "%5d  %10ld  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e",
               p, ParID, amr->Par->Mass[ParID],
               amr->Par->PosX[ParID], amr->Par->PosY[ParID], amr->Par->PosZ[ParID],
               amr->Par->VelX[ParID], amr->Par->VelY[ParID], amr->Par->VelZ[ParID],
               amr->Par->Time[ParID] );
#     ifdef STORE_PAR_ACC
      fprintf( File, "  %13.6e  %13.6e  %13.6e",
               amr->Par->AccX[ParID], amr->Par->AccY[ParID], amr->Par->AccZ[ParID] );
#     endif
      for (int v=0; v<PAR_NPASSIVE; v++)
      fprintf( File, "  %13.6e", amr->Par->Passive[v][ParID] );

      fprintf( File, "\n" );
   }
#  endif // #ifdef PARTICLE


   fclose( File );

} // FUNCTION : Output_Patch

