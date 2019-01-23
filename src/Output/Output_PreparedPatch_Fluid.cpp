#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_PreparedPatch_Fluid
// Description :  Output the fluid data of a "single patch" plus its ghost zones prepared by Flu_Prepare()
//
// Note        :  This function should be placed after invoking Prepare_PatchData() in Flu_Prepare()
//
//                Example:
//                Prepare_PatchData( ... );
//
//                const int TLv  = 1;
//                const int TPID = 12;
//                char comment[10];
//                sprintf( comment, "Step%ld", AdvanceCounter[TLv] );
//                Output_PreparedPatch_Fluid( TLv, TPID, ( real (*)[FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT] )h_Flu_Array_F_In,
//                                            NPG, PID0_List, lv, comment );
//
// Paremeter   :  TLv         : Level you want to output the prepared patch
//                TPID        : Target patch ID
//                h_Flu_Array : Input fluid array for the fluid solver
//                NPG         : Number of patch groups to be prepared at a time
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//                CLv         : Level of data currently stored in h_Flu_Array
//                comment     : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_PreparedPatch_Fluid( const int TLv, const int TPID,
                                 const real h_Flu_Array[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                                 const int NPG, const int *PID0_List, const int CLv, const char *comment )
{

// nothing to do if the current level is not the target level
   if ( TLv != CLv )    return;


   const int LocalID = TPID % 8;
   const int TPID0   = TPID - LocalID;

   for (int TID=0; TID<NPG; TID++)
   {
//    check if the target patch is within the current patch group
      if ( TPID0 != PID0_List[TID] )   continue;


//    begin to output the prepared data
      patch_t *Relation = amr->patch[0][TLv][TPID];

      char FileName[100];
      sprintf( FileName, "PrePatch_Fluid_r%d_lv%d_p%d", MPI_Rank, TLv, TPID );
      if ( comment != NULL )
      {
         strcat( FileName, "_" );
         strcat( FileName, comment );
      }

      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


//    output header
      FILE *File = fopen( FileName, "w" );

      fprintf( File, "Rank %d  Lv %d  PID %d  Local ID %d  Time %13.7e  Step %ld  Counter %ld\n",
               MPI_Rank, TLv, TPID, LocalID, Time[TLv], Step, AdvanceCounter[TLv] );

      fprintf( File, "Father  %d     Son  %d     Corner  (%10d,%10d,%10d)\n\n", Relation->father, Relation->son,
               Relation->corner[0], Relation->corner[1], Relation->corner[2] );

      fprintf( File, "Sibling, Sibling->Son, and Father->Sibling Lists :\n" );

      int Sib, FaSib, SibSon;
      for (int S=0; S<26; S++)
      {
         Sib    = Relation->sibling[S];
         FaSib  = ( TLv ==  0 ) ? -1 : amr->patch[0][TLv-1][Relation->father]->sibling[S];
         SibSon = ( Sib < 0 )   ? Sib : amr->patch[0][TLv][Sib]->son;

         fprintf( File, "Sib[%2d] = %6d     Sib_Son = %6d     Fa_Sib[%2d] = %6d\n",
                  S, Sib, SibSon, S, FaSib );
      }
      fprintf( File, "\n" );


      fprintf( File, "(%3s,%3s,%3s )", "i", "j", "k" );

#     if ( MODEL == ELBDM )
      fprintf( File, "%16s%16s", "Real", "Imag" );

#     else
      for (int v=0; v<FLU_NIN; v++)    fprintf( File, "%16s", FieldLabel[v] );

#     if ( MODEL == HYDRO )
      fprintf( File, "%16s", "Pressure" );
#     endif
#     endif // MODEL

      fprintf( File, "\n" );


      const int Disp_i  = TABLE_02( LocalID, 'x', FLU_GHOST_SIZE, FLU_GHOST_SIZE+PATCH_SIZE );
      const int Disp_j  = TABLE_02( LocalID, 'y', FLU_GHOST_SIZE, FLU_GHOST_SIZE+PATCH_SIZE );
      const int Disp_k  = TABLE_02( LocalID, 'z', FLU_GHOST_SIZE, FLU_GHOST_SIZE+PATCH_SIZE );
      int  K, J, I, Idx;
      real u[FLU_NIN];

      for (int k=-FLU_GHOST_SIZE; k<FLU_GHOST_SIZE+PATCH_SIZE; k++)   {  K = k + Disp_k;
      for (int j=-FLU_GHOST_SIZE; j<FLU_GHOST_SIZE+PATCH_SIZE; j++)   {  J = j + Disp_j;
      for (int i=-FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PATCH_SIZE; i++)   {  I = i + Disp_i;

         Idx = K*FLU_NXT*FLU_NXT + J*FLU_NXT + I;

         for (int v=0; v<FLU_NIN; v++)    u[v] = h_Flu_Array[TID][v][Idx];

//       output cell indices
         fprintf( File, "(%3d,%3d,%3d )", i, j, k );

//       output all variables in the prepared fluid array
         for (int v=0; v<FLU_NIN; v++)    fprintf( File, "  %14.7e", u[v] );

//       output pressure in HYDRO
#        if   ( MODEL == HYDRO )
         const bool CheckMinPres_No = false;
         fprintf( File, "  %14.7e", Hydro_GetPressure(u[DENS],u[MOMX],u[MOMY],u[MOMZ],u[ENGY],GAMMA-1.0,CheckMinPres_No,NULL_REAL) );

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!
#        endif // MODEL

         fprintf( File, "\n" );
      }}}

      fclose( File );
      break;

   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Output_PreparedPatch_Fluid
