#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_Pot
// Description :  Fill up the h_Pot_Array_P_In array with the coarse-grid (lv-1) potential for the Poisson solver
//
// Note        :  1. We don't call "Prepare_PatchData" since it always prepare data for a **patch group** but
//                   here we may want to prepare data for a single **patch**
//                2. Use PrepTime to determine the physical time to prepare data
//                   --> Temporal interpolation/extrapolation will be conducted automatically if PrepTime
//                       is NOT equal to the time of data stored previously (e.g., FluSgTime[0/1])
//
// Parameter   :  lv                : Target refinement level
//                PrepTime          : Target physical time to prepare the coarse-grid data
//                h_Pot_Array_P_In  : Host array to store the prepared coarse-grid potential
//                NPG               : Number of patch groups prepared at a time
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Poi_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT],
                      const int NPG, const int *PID0_List )
{

// nothing to do if there is no target patch group
// --> note that this check is necessary for levels without any patch, where Time may be ill-defined and thus
//     the following check code "if ( TimeMin < 0.0 )" will fail
   if ( NPG == 0 )   return;


// check
#  ifdef GAMER_DEBUG
   if ( lv == 0 )    Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

// temporal interpolation is unnecessary for the shared time-step integration
   if ( OPT__DT_LEVEL == DT_LEVEL_SHARED  &&  OPT__INT_TIME )
      Aux_Error( ERROR_INFO, "OPT__INT_TIME should be disabled when \"OPT__DT_LEVEL == DT_LEVEL_SHARED\" !!\n" );
#  endif


//###REVISE: support interpolation schemes requiring 2 ghost cells on each side
   const int IntGhost = 1;    // assuming interpolation ghost-zone == 1
   const int CGhost   = (POT_GHOST_SIZE+1)/2 + IntGhost;
   const int CWidth   = PATCH_SIZE + 2*CGhost;
   const int FaLv     = lv - 1;


// temporal interpolation parameters
   bool PotIntTime;
   int  PotSg, PotSg_IntT;
   real PotWeighting, PotWeighting_IntT;

   if      (  Mis_CompareRealValue( PrepTime, amr->PotSgTime[FaLv][   amr->PotSg[FaLv] ], NULL, false )  )
   {
      PotIntTime        = false;
      PotSg             = amr->PotSg[FaLv];
      PotSg_IntT        = NULL_INT;
      PotWeighting      = NULL_REAL;
      PotWeighting_IntT = NULL_REAL;
   }

   else if (  Mis_CompareRealValue( PrepTime, amr->PotSgTime[FaLv][ 1-amr->PotSg[FaLv] ], NULL, false )  )
   {
      PotIntTime        = false;
      PotSg             = 1 - amr->PotSg[FaLv];
      PotSg_IntT        = NULL_INT;
      PotWeighting      = NULL_REAL;
      PotWeighting_IntT = NULL_REAL;
   }

   else
   {
//    check
      if ( OPT__DT_LEVEL == DT_LEVEL_SHARED )
      Aux_Error( ERROR_INFO, "cannot determine PotSg for OPT__DT_LEVEL == DT_LEVEL_SHARED (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                 FaLv, PrepTime, amr->PotSgTime[FaLv][0], amr->PotSgTime[FaLv][1] );

//    print warning messages if temporal extrapolation is required
      const double TimeMin = MIN( amr->PotSgTime[FaLv][0], amr->PotSgTime[FaLv][1] );
      const double TimeMax = MAX( amr->PotSgTime[FaLv][0], amr->PotSgTime[FaLv][1] );

      if ( TimeMin < 0.0 )
         Aux_Error( ERROR_INFO, "TimeMin (%21.14e) < 0.0 ==> one of the potential arrays has not been initialized !!\n", TimeMin );

      if ( PrepTime < TimeMin  ||  PrepTime-TimeMax >= 1.0e-12*TimeMax )
         Aux_Message( stderr, "WARNING : temporal extrapolation (lv %d, T_Prep %20.14e, T_Min %20.14e, T_Max %20.14e)\n",
                      FaLv, PrepTime, TimeMin, TimeMax );

      if ( OPT__INT_TIME )
      {
         PotIntTime        = true;
         PotSg             = 0;
         PotSg_IntT        = 1;
         PotWeighting      =   ( +amr->PotSgTime[FaLv][PotSg_IntT] - PrepTime )
                             / (  amr->PotSgTime[FaLv][PotSg_IntT] - amr->PotSgTime[FaLv][PotSg] );
         PotWeighting_IntT =   ( -amr->PotSgTime[FaLv][PotSg     ] + PrepTime )
                             / (  amr->PotSgTime[FaLv][PotSg_IntT] - amr->PotSgTime[FaLv][PotSg] );
      }

      else
      {
         PotIntTime        = false;
         PotSg             = amr->PotSg[FaLv]; // set to the current Sg
         PotSg_IntT        = NULL_INT;
         PotWeighting      = NULL_REAL;
         PotWeighting_IntT = NULL_REAL;
      }
   } // Mis_CompareRealValue


#  pragma omp parallel
   {
//    CPot : store the coarse-grid potential for interpoaltion
      real (*CPot)[CWidth][CWidth] = new real [CWidth][CWidth][CWidth];

      int FaPID, FaSibPID, PID0;
      int N, I, J, K, I2, J2, K2, Disp_i, Disp_j, Disp_k, Disp_i2, Disp_j2, Disp_k2, Loop_i, Loop_j, Loop_k;


//    prepare the coarse-grid potential for eight patches (one patch group) at a time
#     pragma omp for schedule( static )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0  = PID0_List[TID];
         FaPID = amr->patch[0][lv][PID0]->father;

#        ifdef GAMER_DEBUG
         if ( FaPID < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPID", FaPID );
#        endif

//       a. fill up the central region (ghost zone is not filled up yet) of the CPot array
// ------------------------------------------------------------------------------------------------------------
         for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost;

            CPot[K][J][I] = amr->patch[PotSg][FaLv][FaPID]->pot[k][j][i];

            if ( PotIntTime ) // temporal interpolation
            CPot[K][J][I] =   PotWeighting     *CPot[K][J][I]
                            + PotWeighting_IntT*amr->patch[PotSg_IntT][FaLv][FaPID]->pot[k][j][i];
         }}}


//       b. fill up the ghost zone (no spatial interpolation is required) of the CPot array
// ------------------------------------------------------------------------------------------------------------
         for (int sib=0; sib<26; sib++)   // we always need coarse-grid ghost zones in 26 sibling directions
         {
            FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

#           ifdef GAMER_DEBUG
            if ( FaSibPID < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaSibPID", FaSibPID );
#           endif

            Loop_i   = TABLE_01( sib, 'x', CGhost, PATCH_SIZE, CGhost );
            Loop_j   = TABLE_01( sib, 'y', CGhost, PATCH_SIZE, CGhost );
            Loop_k   = TABLE_01( sib, 'z', CGhost, PATCH_SIZE, CGhost );
            Disp_i   = TABLE_01( sib, 'x', 0, CGhost, CGhost+PATCH_SIZE );
            Disp_j   = TABLE_01( sib, 'y', 0, CGhost, CGhost+PATCH_SIZE );
            Disp_k   = TABLE_01( sib, 'z', 0, CGhost, CGhost+PATCH_SIZE );
            Disp_i2  = TABLE_01( sib, 'x', PATCH_SIZE-CGhost, 0, 0 );
            Disp_j2  = TABLE_01( sib, 'y', PATCH_SIZE-CGhost, 0, 0 );
            Disp_k2  = TABLE_01( sib, 'z', PATCH_SIZE-CGhost, 0, 0 );

            for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
            for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
            for (int i=0; i<Loop_i; i++)  {  I = i + Disp_i;   I2 = i + Disp_i2;

               CPot[K][J][I] = amr->patch[PotSg][FaLv][FaSibPID]->pot[K2][J2][I2];

               if ( PotIntTime ) // temporal interpolation
               CPot[K][J][I] =   PotWeighting     *CPot[K][J][I]
                               + PotWeighting_IntT*amr->patch[PotSg_IntT][FaLv][FaSibPID]->pot[K2][J2][I2];
            }}}
         } // for (int sib=0; sib<26; sib++)


//       c. copy data from the CPot array to the h_Pot_Array_P_In array
// ------------------------------------------------------------------------------------------------------------
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            N      = 8*TID + LocalID;
            Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE/2 );
            Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE/2 );
            Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE/2 );

            for (int k=0; k<POT_NXT; k++)    {  K = k + Disp_k;
            for (int j=0; j<POT_NXT; j++)    {  J = j + Disp_j;
            for (int i=0; i<POT_NXT; i++)    {  I = i + Disp_i;

               h_Pot_Array_P_In[N][k][j][i] = CPot[K][J][I];

            }}}
         }
      } // for (int TID=0; TID<NPG; TID++)


      delete [] CPot;

   } // OpenMP parallel region

} // FUNCTION : Poi_Prepare_Pot



#endif // #ifdef GRAVITY
