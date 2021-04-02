#include "GAMER.h"

#ifdef FEEDBACK



// user-specified feedback to be set by a test problem initializer
void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                     const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                     real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[] ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_AdvanceDt
// Description :  Feedback from particles to grids
//
// Note        :  1. Invoked by EvolveLevel()
//                2. FB_LEVEL must equal MAX_LEVEL for now
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                SaveSg_Flu : Sandglass to store the updated fluid data
//                SaveSg_Mag : Sandglass to store the updated B field
//
// Return      :  Update both grids and particles
//-------------------------------------------------------------------------------------------------------
void FB_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const int SaveSg_Flu, const int SaveSg_Mag )
{

// only work on FB_LEVEL for now
   if ( lv != FB_LEVEL )   return;


// 1. collect particles for the sibling buffer patches
//    --> exclude father-sibling buffer patches (FaSibBufPatch_No) since currently we assume FB_LEVEL == MAX_LEVEL
//    --> for simplicity, we disable particle postition prediction (PredictPos_No) currently even though particles
//        just crossing from coarse (lv-1) to fine (lv) grids may have time greater than other particles at lv
//        --> to fix it, we will need to correct other particle attributes such as velocity too
   const bool TimingSendPar_Yes = true;
   const bool JustCountNPar_No  = false;
   const bool PredictPos_No     = false;
   const bool SibBufPatch_Yes   = true;
   const bool FaSibBufPatch_No  = false;
//###OPTIMIZATION: only collect necessary particle attributes
   const long ParAttBitIdx      = _PAR_TOTAL;

   Par_CollectParticle2OneLevel( lv, ParAttBitIdx, PredictPos_No, TimeNew, SibBufPatch_Yes, FaSibBufPatch_No,
                                 JustCountNPar_No, TimingSendPar_Yes );


// get the sibling index differences along different directions
   int NSibPID_Delta[26], *SibPID_Delta[26];

   TABLE_GetSibPID_Delta( NSibPID_Delta, SibPID_Delta );

// start of OpenMP parallel region
#  pragma omp parallel
   {

// per-thread variables
   const int NNearbyPatchMax = 64;  // maximum number of neaby patches of a patch group (including 8 local patches)
   int Nearby_PID_List[NNearbyPatchMax], NNearbyPatch, SibPID0_List[26];

   real (*fluid_PG)[PS2][PS2][PS2] = new real [NCOMP_TOTAL][PS2][PS2][PS2];


// iterate over all real patches
#  pragma omp for schedule( runtime )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
      const double PGCenter[3] = { amr->patch[0][lv][PID0+7]->EdgeL[0],
                                   amr->patch[0][lv][PID0+7]->EdgeL[1],
                                   amr->patch[0][lv][PID0+7]->EdgeL[2] };

//    2. prepare the fluid data to be updated
//    --> exclude magnetic field for now
//    --> use patch group as the basic unit
      const int  GhostZone_No        = 0;
      const int  NPG                 = 1;
#     ifndef MHD
      const int  OPT__MAG_INT_SCHEME = INT_NONE;
#     endif
      const bool IntPhase_No         = false;
      const real MinDens_No          = -1.0;
      const real MinPres_No          = -1.0;
      const real MinTemp_No          = -1.0;
      const bool DE_Consistency_No   = false;

      Prepare_PatchData( lv, TimeNew, fluid_PG[0][0][0], NULL, GhostZone_No, NPG, &PID0, _TOTAL, _NONE,
                         OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );



//    3. collect patch information
//    3-1. get nearby patches
//       --> this is not really necessary since currently we use patches instead of patch groups as the basic unit
//           when calling the feedback routines
//       --> but it will be necessary if we switch to patch groups in the future
      NNearbyPatch = 0;

//    local patches
      for (int PID=PID0; PID<PID0+8; PID++)  Nearby_PID_List[ NNearbyPatch ++ ] = PID;

//    sibling patches
      TABLE_GetSibPID_Based( lv, PID0, SibPID0_List );

//###OPTIMIZATION: skip sibling patches if the maximum feedback radius is zero
      for (int s=0; s<26; s++)
      {
         const int SibPID0 = SibPID0_List[s];   // first target patch in the sibling patch group

//       only consider leaf patches on FB_LEVEL (including both real and buffer patches)
         if ( SibPID0 >= 0 )
         for (int c=0; c<NSibPID_Delta[s]; c++)
         {
            const int SibPID = SibPID0 + SibPID_Delta[s][c];
            Nearby_PID_List[ NNearbyPatch ++ ] = SibPID;
         }
      }


//    3-2. record the coarse-fine boundaries
//         --> regard non-periodic boundaries as coarse-fine boundaries too
      bool CoarseFine[26];
      for (int s=0; s<26; s++)   CoarseFine[s] = ( SibPID0_List[s] < 0 ) ? true : false;



//    4. allocate arrays to store the local particle data
//       --> **local** means particles in the leaf **real** patches of this MPI rank
//       --> exclude **buffer** patches since their particle data have been stored in ParAtt_Copy[] already
//       --> allocate the **maximum** required size among all nearby patches of a given patch group **just once**
//           for better performance
      real *ParAtt_Local[PAR_NATT_TOTAL];
      int  *ParSortID = NULL;
      int   NParMax = -1;

      for (int t=0; t<NNearbyPatch; t++)
      {
         const int PID = Nearby_PID_List[t];

//       check both NPar and NPar_Copy (NPar_Copy may be -1, which is fine)
         NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar      );
         NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar_Copy );
      }

      if ( NParMax > 0 )
      {
         for (int v=0; v<PAR_NATT_TOTAL; v++)
         {
            if ( ParAttBitIdx & BIDX(v) )    ParAtt_Local[v] = new real [NParMax];
            else                             ParAtt_Local[v] = NULL;
         }

         ParSortID = new int [NParMax];   // it will fail if "long" is actually required for NParMax
      }


//    iterate over all nearby patches to apply feedback
      for (int t=0; t<NNearbyPatch; t++)
      {
         const int PID = Nearby_PID_List[t];

//       5. prepare the input particle data
//       5-1. get the particle list of the target patch
         long  *ParList = NULL;
         int    NPar;
         bool   UseParAttCopy;

//       the check "son == -1" is actually useless for now since we fix FB_LEVEL == MAX_LEVEL
         if ( amr->patch[0][lv][PID]->son == -1  &&  PID < amr->NPatchComma[lv][1] )
         {
            NPar          = amr->patch[0][lv][PID]->NPar;
            ParList       = amr->patch[0][lv][PID]->ParList;
            UseParAttCopy = false;

#           ifdef DEBUG_PARTICLE
            if ( amr->patch[0][lv][PID]->NPar_Copy != -1 )
               Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                          lv, PID, amr->patch[0][lv][PID]->NPar_Copy );
#           endif
         }

         else
         {
//          note that amr->patch[0][lv][PID]->NPar>0 is still possible
            NPar          = amr->patch[0][lv][PID]->NPar_Copy;
#           ifdef LOAD_BALANCE
            ParList       = NULL;
            UseParAttCopy = true;
#           else
            ParList       = amr->patch[0][lv][PID]->ParList_Copy;
            UseParAttCopy = false;
#           endif
         } // if ( amr->patch[0][lv][PID]->son == -1  &&  PID < amr->NPatchComma[lv][1] ) ... else ...


//       5-2. prepare the particle data
         if ( NPar <= 0 )  continue;   // skip patches without any particle

         real *ParAtt_Ptr[PAR_NATT_TOTAL];

#        ifdef LOAD_BALANCE
         if ( UseParAttCopy )
         {
            for (int v=0; v<PAR_NATT_TOTAL; v++)
            {
               ParAtt_Ptr[v] = amr->patch[0][lv][PID]->ParAtt_Copy[v];

#              ifdef DEBUG_PARTICLE
               if (  ( ParAttBitIdx & BIDX(v) )  &&  ParAtt_Ptr[v] == NULL )
                  Aux_Error( ERROR_INFO, "ParAtt_Ptr[%d] == NULL (NPar %d, lv %d, PID %d) !!\n",
                             v, NPar, lv, PID );
#              endif
            }
         }

         else
#        endif // #ifdef LOAD_BALANCE
         {
#           ifdef DEBUG_PARTICLE
            if ( ParList == NULL )
               Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                          NPar, lv, PID );
#           endif

            for (int v=0; v<PAR_NATT_TOTAL; v++)
            {
               if ( ParAttBitIdx & BIDX(v) )
               {
                  for (int p=0; p<NPar; p++)
                  {
                     const long ParID = ParList[p];
                     ParAtt_Local[v][p] = amr->Par->Attribute[v][ParID];
                  }
               }

               ParAtt_Ptr[v] = ParAtt_Local[v];
            }
         } // if ( UseParAttCopy ) ... else ...


//       5-3. sort particles by positions to fix their order
//            --> necessary when feedback involves random numbers
//            --> because otherwise the same particle in real and buffer patches may have different random numbers
         Par_SortByPos( NPar, ParAtt_Ptr[PAR_POSX], ParAtt_Ptr[PAR_POSY], ParAtt_Ptr[PAR_POSZ], ParSortID );



//       5-4. periodic boundary conditions
//            --> assuming there are at least TWO patch groups along each spatial direction (i.e., NX0_TOT[*]/PS2>=2)
//            --> so each particle won't affect the target patch group more than once
         const bool   Periodic [3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                       OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                       OPT__BC_FLU[4] == BC_FLU_PERIODIC };
         const double HalfBox  [3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
         const int    ParPosIdx[3] = { PAR_POSX, PAR_POSY, PAR_POSZ };

         if ( Periodic[0]  ||  Periodic[1]  ||  Periodic[2] )
         for (int p=0; p<NPar; p++)
         {
            for (int d=0; d<3; d++)
            {
               if ( Periodic[d] )
               {
                  real *ParPos = ParAtt_Ptr[ ParPosIdx[d] ] + p;
                  const double dr = *ParPos - PGCenter[d];

                  if      ( dr > +HalfBox[d] )  *ParPos -= amr->BoxSize[d];
                  else if ( dr < -HalfBox[d] )  *ParPos += amr->BoxSize[d];
               }
            }
         } // for (int p=0; p<NPar; p++)


//       6. invoke all feedback routines
         if ( FB_USER )    FB_User_Ptr( lv, TimeNew, TimeOld, dt, NPar, ParSortID, ParAtt_Ptr, fluid_PG,
                                        amr->patch[0][lv][PID0]->EdgeL, amr->dh[lv], CoarseFine );

      } // for (int t=0; t<NNearbyPatch; t++)



//    7. store the updated fluid data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID    = PID0 + LocalID;
         const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1 );
         const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1 );
         const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1 );

         for (int v=0; v<NCOMP_TOTAL; v++)   {
         for (int k_o=0; k_o<PS1; k_o++)     {  const int k_i = Disp_k + k_o;
         for (int j_o=0; j_o<PS1; j_o++)     {  const int j_i = Disp_j + j_o;
         for (int i_o=0; i_o<PS1; i_o++)     {  const int i_i = Disp_i + i_o;

            amr->patch[SaveSg_Flu][lv][PID]->fluid[v][k_o][j_o][i_o] = fluid_PG[v][k_i][j_i][i_i];

         }}}}
      } // for (int LocalID=0; LocalID<8; LocalID++)

//    free memory
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParAtt_Local[v];
      delete [] ParSortID;

   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

// free memory
   delete [] fluid_PG;

   } // end of OpenMP parallel region


// free particles
   Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_Yes, FaSibBufPatch_No );

} // FUNCTION : FB_AdvanceDt



#endif // #ifdef FEEDBACK
