#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_PreparedPatch_Poisson
// Description :  Output the density/potential of a single patch plus its ghost zones prepared by
//                Poi_Prepare_Rho() or Poi_Prepare_Pot()
//
// Note        :  For "TComp == 0", this function should be placed after invoking Prepare_PatchData() in
//                Poi_Prepare_Rho().
//                For "TComp == 1", this function should be placed in the end of Poi_Prepare_Pot().
//
//                Example:
//                const int TLv   = 1;
//                const int TPID  = 12;
//                char comment[10];
//                sprintf( comment, "Step%ld", AdvanceCounter[TLv] );
//
//                const int TComp = 0;
//                Output_PreparedPatch_Poisson( TLv, TPID, TComp, h_Input_Array, NULL, NPG, PID0_List, lv,
//                                              comment );
//
//                const int TComp = 1;
//                Output_PreparedPatch_Poisson( TLv, TPID, TComp, NULL, h_Pot_Array_P_In, NPG, PID0_List, lv,
//                                              comment);
//
// Paremeter   :  TLv              : Level you want to output the prepared patch
//                TPID             : Target patch ID
//                TComp            : Target component
//                                    --> 0 : Output the prepared density   stored in Rho_Array_P
//                                        1 : Output the prepared potential stored in Pot_Array_P_In
//                h_Rho_Array_P    : Input host density array for the Poisson solver
//                h_Pot_Array_P_In : Input host coarse-grid potential array for the Poisson solver
//                NPG              : Number of patch groups to be prepared at a time
//                PID0_List        : List recording the patch indices with LocalID==0 to be udpated
//                CLv              : Level of data currently stored in the input arrays
//                comment          : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_PreparedPatch_Poisson( const int TLv, const int TPID, const int TComp,
                                   const real h_Rho_Array_P   [][RHO_NXT][RHO_NXT][RHO_NXT],
                                   const real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT],
                                   const int NPG, const int *PID0_List, const int CLv, const char *comment )
{

// check
   if ( TComp != 0  &&  TComp != 1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TComp", TComp );


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
      sprintf( FileName, "PrePatch_Poisson_r%d_lv%d_p%d", MPI_Rank, TLv, TPID );
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


//    output data
//###REVISE: support interpolation schemes requiring 2 ghost cells on each side
      const int IntPotGhost = 1;    // assuming interpolation ghost-zone of potential == 1
      const int PotCGhost   = (POT_GHOST_SIZE+1)/2 + IntPotGhost;
      const int N           = 8*TID + LocalID;
      int I, J, K;

      switch ( TComp )
      {
         case 0:
            fprintf( File, "(%3s,%3s,%3s )%16s\n", "i", "j", "k", "Density" );

            for (int k=0; k<RHO_NXT; k++)    {  K = k - RHO_GHOST_SIZE;
            for (int j=0; j<RHO_NXT; j++)    {  J = j - RHO_GHOST_SIZE;
            for (int i=0; i<RHO_NXT; i++)    {  I = i - RHO_GHOST_SIZE;

               fprintf( File, "(%3d,%3d,%3d )  %14.7e\n", I, J, K, h_Rho_Array_P[N][k][j][i] );

            }}}
            break;


         case 1:
            fprintf( File, "(%3s,%3s,%3s )%16s\n", "i", "j", "k", "Potential_In" );

            for (int k=0; k<POT_NXT; k++)    {  K = k - PotCGhost;
            for (int j=0; j<POT_NXT; j++)    {  J = j - PotCGhost;
            for (int i=0; i<POT_NXT; i++)    {  I = i - PotCGhost;

               fprintf( File, "(%3d,%3d,%3d )  %14.7e\n", I, J, K, h_Pot_Array_P_In[N][k][j][i] );

            }}}
            break;

         default :
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TComp", TComp );

      } // switch ( TComp )

      fclose( File );
      break;

   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Output_PreparedPatch_Poisson



#endif // #ifdef GRAVITY
