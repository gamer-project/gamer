#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_Pot
// Description :  Fill up the h_Pot_Array_P_In array with the coarse-grid (lv-1) potential for the Poisson solver
//
// Note        :  1. We don't call Prepare_PatchData() since it always prepare data for a **patch group** but
//                   here we may want to prepare data for a single **patch**
//                2. Use PrepTime to determine the physical time to prepare data
//                   --> Temporal interpolation/extrapolation will be conducted automatically if PrepTime
//                       is NOT equal to the time of data stored previously (i.e., PotSgTime[0/1])
//                3. Use extrapolation to obtain data outside the non-periodic boundaries
//
// Parameter   :  lv               : Target refinement level
//                PrepTime         : Target physical time to prepare the coarse-grid data
//                h_Pot_Array_P_In : Host array to store the prepared coarse-grid potential
//                NPG              : Number of patch groups prepared at a time
//                PID0_List        : List recording the patch indicies with LocalID==0 to be udpated
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
// --> it's ok for now as the interpolation schemes adopted in the Poisson solver (CQUAD/QUAD) only require 1 ghost cell
   const int IntGhost = 1;    // assuming interpolation ghost-zone == 1
   const int CGhost   = (POT_GHOST_SIZE+1)/2 + IntGhost;
   const int CWidth   = PS1 + 2*CGhost;
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


// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   int BC_Face[26], BC_Face_tmp[3];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1  &&  OPT__BC_FLU[BC_Face_tmp[2]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1  &&  OPT__BC_FLU[BC_Face_tmp[1]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1  &&  OPT__BC_FLU[BC_Face_tmp[0]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[0];
      else                                                                                   BC_Face[s] = NULL_INT;
   }


#  pragma omp parallel
   {
//    CPot : store the coarse-grid potential for interpoaltion
      real (*CPot)[CWidth][CWidth] = new real [CWidth][CWidth][CWidth];

      int FaPID, FaSibPID, PID0, Idx_Start[3], Idx_End[3], Idx_Start2[3], BC_Sibling;

//    prepare the coarse-grid potential for eight patches (one patch group) at a time
#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0  = PID0_List[TID];
         FaPID = amr->patch[0][lv][PID0]->father;

#        ifdef GAMER_DEBUG
         if ( FaPID < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPID", FaPID );
#        endif

//       a. fill up the central region of CPot[] (excluding ghost zones)
// ------------------------------------------------------------------------------------------------------------
         for (int k=0, K=CGhost; k<PS1; k++, K++)
         for (int j=0, J=CGhost; j<PS1; j++, J++)
         for (int i=0, I=CGhost; i<PS1; i++, I++)
         {
            CPot[K][J][I] = amr->patch[PotSg][FaLv][FaPID]->pot[k][j][i];

            if ( PotIntTime ) // temporal interpolation
            CPot[K][J][I] =   PotWeighting     *CPot[K][J][I]
                            + PotWeighting_IntT*amr->patch[PotSg_IntT][FaLv][FaPID]->pot[k][j][i];
         }


//       b. fill up the ghost zones of CPot[] (do not require spatial interpolation from FaLv-1)
// ------------------------------------------------------------------------------------------------------------
         for (int sib=0; sib<26; sib++)   // we always need coarse-grid ghost zones in 26 sibling directions
         {
            FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[sib];

            for (int d=0; d<3; d++)
            {
               Idx_Start[d] = TABLE_01( sib, 'x'+d, 0, CGhost, CGhost+PS1 );
               Idx_End  [d] = TABLE_01( sib, 'x'+d, CGhost, PS1, CGhost ) + Idx_Start[d] - 1;
            }


//          b1. if the target sibling patch exists --> copy data directly
            if ( FaSibPID >= 0 )
            {
               for (int d=0; d<3; d++)    Idx_Start2[d] = TABLE_01( sib, 'x'+d, PS1-CGhost, 0, 0 );

               for (int k=Idx_Start[2], K=Idx_Start2[2]; k<=Idx_End[2]; k++, K++)
               for (int j=Idx_Start[1], J=Idx_Start2[1]; j<=Idx_End[1]; j++, J++)
               for (int i=Idx_Start[0], I=Idx_Start2[0]; i<=Idx_End[0]; i++, I++)
               {
                  CPot[k][j][i] = amr->patch[PotSg][FaLv][FaSibPID]->pot[K][J][I];

                  if ( PotIntTime ) // temporal interpolation
                  CPot[k][j][i] =   PotWeighting     *CPot[k][j][i]
                                  + PotWeighting_IntT*amr->patch[PotSg_IntT][FaLv][FaSibPID]->pot[K][J][I];
               }
            }


//          b2. if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( FaSibPID <= SIB_OFFSET_NONPERIODIC )
            {
               BC_Sibling = SIB_OFFSET_NONPERIODIC - FaSibPID;

#              ifdef GAMER_DEBUG
               if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
                  Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

               if ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
                  Aux_Error( ERROR_INFO, "FluBC == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, FaPID %d, FaSiBPID %d, PID0 %d, sib %d) !!\n",
                             BC_Sibling, BC_Face[BC_Sibling], FaPID, FaSibPID, PID0, sib );
#              endif

//             extrapolate potential (currently it's the only supported non-periodic BC for potential)
               Poi_BoundaryCondition_Extrapolation( CPot[0][0], BC_Face[BC_Sibling], 1, CGhost,
                                                    CWidth, CWidth, CWidth, Idx_Start, Idx_End );
            }


//          b3. if FaSibPID == -1, something is wrong as it violates the proper-nesting constraint
            else
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaSibPID", FaSibPID );
         } // for (int sib=0; sib<26; sib++)


//       c. copy data from CPot[] to h_Pot_Array_P_In[]
// ------------------------------------------------------------------------------------------------------------
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int N = 8*TID + LocalID;

            for (int d=0; d<3; d++)    Idx_Start2[d] = TABLE_02( LocalID, 'x'+d, 0, PS1/2 );

            for (int k=0, K=Idx_Start2[2]; k<POT_NXT; k++, K++)
            for (int j=0, J=Idx_Start2[1]; j<POT_NXT; j++, J++)
            for (int i=0, I=Idx_Start2[0]; i<POT_NXT; i++, I++)
               h_Pot_Array_P_In[N][k][j][i] = CPot[K][J][I];
         }
      } // for (int TID=0; TID<NPG; TID++)

      delete [] CPot;
   } // OpenMP parallel region

} // FUNCTION : Poi_Prepare_Pot



#endif // #ifdef GRAVITY
