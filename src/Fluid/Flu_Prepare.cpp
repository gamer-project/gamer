#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Prepare
// Description :  Prepare input arrays for the fluid solver
//
// Note        :  Invoke Prepare_PatchData()
//
// Parameter   :  lv                    : Target refinement level
//                PrepTime              : Target physical time to prepare the coarse-grid data
//                h_Flu_Array_F_In      : Host array to store the prepared fluid data
//                h_Mag_Array_F_In      : Host array to store the prepared B field (for MHD onlhy)
//                h_Pot_Array_USG_F     : Host array to store the prepared potential data (for UNSPLIT_GRAVITY only)
//                h_Corner_Array_USG_F  : Host array to store the prepared corner data (for UNSPLIT_GRAVITY only)
//                h_IsCompletelyRefined : Host array storing which patch groups are completely refined ( ELBDM only )
//                h_HasWaveCounterpart  : Host array storing which cells have wave counterpart ( ELBDM_HYBRID only )
//                NPG                   : Number of patch groups to be prepared at a time
//                PID0_List             : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Flu_Prepare( const int lv, const double PrepTime,
                  real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                  real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                  real h_Pot_Array_USG_F[][ CUBE(USG_NXT_F) ],
                  double h_Corner_Array_F[][3],
                  bool h_IsCompletelyRefined[],
                  bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                  const int NPG, const int *PID0_List,
                  LB_GlobalTree* GlobalTree )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef UNSPLIT_GRAVITY
   if (  ( OPT__SELF_GRAVITY || OPT__EXT_POT )  &&  h_Pot_Array_USG_F == NULL  )
      Aux_Error( ERROR_INFO, "h_Pot_Array_USG_F == NULL !!\n" );

   if ( OPT__EXT_ACC  &&  h_Corner_Array_F == NULL )
      Aux_Error( ERROR_INFO, "h_Corner_Array_F == NULL !!\n" );
#  endif
#  endif


#  ifndef MHD
   const int  OPT__MAG_INT_SCHEME = INT_NONE;
#  endif
   const bool IntPhase_No         = false;
   const real MinDens_No          = -1.0;
   const real MinPres_No          = -1.0;
   const real MinTemp_No          = -1.0;
   const real MinEntr_No          = -1.0;
   const bool DE_Consistency_Yes  = true;
   const bool DE_Consistency_No   = false;
   const bool DE_Consistency      = ( OPT__OPTIMIZE_AGGRESSIVE ) ? DE_Consistency_No : DE_Consistency_Yes;
   const real MinDens             = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;


// prepare the fluid array
// --> exclude passive scalars for ELBDM for now since it is not supported yet
// --> consistent with FLU_NIN == NCOMP_FLUID - 1 in Macro.h
#  if ( MODEL == ELBDM )

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
      Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], NULL,
//                      FLU_GHOST_SIZE, NPG, PID0_List, _REAL|_IMAG|_PASSIVE, _NONE,
                        FLU_GHOST_SIZE, NPG, PID0_List, _REAL|_IMAG, _NONE,
                        OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, OPT__INT_PHASE,
                        OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
      Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], NULL,
//                      HYB_GHOST_SIZE, NPG, PID0_List, _DENS|_PHAS|_PASSIVE, _NONE,
                        HYB_GHOST_SIZE, NPG, PID0_List, _DENS|_PHAS, _NONE,
                        OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, OPT__INT_PHASE,
                        OPT__BC_FLU, BC_POT_NONE, MinDens,    MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
   }
#  endif

#  else // #if ( MODEL == ELBDM )

#  ifdef MHD
   real *Mag_Array = h_Mag_Array_F_In[0][0];
#  else
   real *Mag_Array = NULL;
#  endif
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], Mag_Array,
                      FLU_GHOST_SIZE, NPG, PID0_List, _TOTAL, _MAG,
                      OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, BC_POT_NONE, MinDens, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency );

#  endif // #if ( MODEL == ELBDM ) ... else ...

#  ifdef UNSPLIT_GRAVITY
// prepare the potential array
   if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
   Prepare_PatchData( lv, PrepTime, h_Pot_Array_USG_F[0], NULL,
                      USG_GHOST_SIZE_F, NPG, PID0_List, _POTE, _NONE,
                      OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


// prepare the corner array
   if ( OPT__EXT_ACC )
   {
      const double dh_half = 0.5*amr->dh[lv];

      int PID0;

//#     pragma omp parallel for private( PID0 ) schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

//       not considering ghost zones
         for (int d=0; d<3; d++)    h_Corner_Array_F[TID][d] = amr->patch[0][lv][PID0]->EdgeL[d] + dh_half;
      } // for (int TID=0; TID<NPG; TID++)
   }
#  endif // #ifdef UNSPLIT_GRAVITY


// prepare boolean array that indicates whether patch group is fully refined (in other words has 8*8=64 children)
#  if ( MODEL == ELBDM )
#  pragma omp parallel for schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      bool PGIsCompletelyRefined = true;

      const int PID0 = PID0_List[TID];
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;
         if ( amr->patch[0][lv][PID]->son == -1 )
         {
            PGIsCompletelyRefined = false;
            break;
         }
      }

      h_IsCompletelyRefined[TID] = PGIsCompletelyRefined;
   }
#  endif


// prepare h_HasWaveCounterpart with information which cells have wave counterparts
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( !amr->use_wave_flag[lv] )
      Prepare_PatchData_HasWaveCounterpart( lv, h_HasWaveCounterpart, HYB_GHOST_SIZE, NPG, PID0_List, NSIDE_26, GlobalTree );
#  endif


// validate input arrays for debugging purposes
   if ( OPT__CK_INPUT_FLUID )
   {
      bool CheckFailed = false;

#     pragma omp parallel for reduction ( ||:CheckFailed ) schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         real fluid[FLU_NIN];

         int flu_nxt = FLU_NXT;

//       distinguish between FLU_NXT and HYB_NXT
#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         if ( !amr->use_wave_flag[lv] )   flu_nxt = HYB_NXT;
#        endif

//       a. fluid
         for (int k=0; k<flu_nxt; k++)
         for (int j=0; j<flu_nxt; j++)
         for (int i=0; i<flu_nxt; i++)
         {
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            if ( amr->use_wave_flag[lv] ) {
#           endif

            const int t = IDX321( i, j, k, FLU_NXT, FLU_NXT );
            for (int v=0; v<FLU_NIN; v++)    fluid[v] = h_Flu_Array_F_In[TID][v][t];


#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            } else {
//          convert to smaller array for HYB_GHOST_SIZE ghost zones
            real (*smaller_h_Flu_Array_F_In)[FLU_NIN][ CUBE(HYB_NXT) ] = (real (*)[FLU_NIN][ CUBE(HYB_NXT) ]) h_Flu_Array_F_In;
            const int t = IDX321( i, j, k, HYB_NXT, HYB_NXT );

            for (int v=0; v<FLU_NIN; v++)    fluid[v] = smaller_h_Flu_Array_F_In[TID][v][t];
            }
#           endif

//          HYDRO
#           if ( MODEL == HYDRO )
            real Emag=NULL_REAL;

#           if ( FLU_NIN != NCOMP_TOTAL )
#           error : ERROR : FLU_NIN != NCOMP_TOTAL for HYDRO !!
#           endif

#           ifdef MHD
            Emag = MHD_GetCellCenteredBEnergy( h_Mag_Array_F_In[TID][MAGX],
                                               h_Mag_Array_F_In[TID][MAGY],
                                               h_Mag_Array_F_In[TID][MAGZ],
                                               FLU_NXT, FLU_NXT, FLU_NXT, i, j, k );
#           endif

            if (  Hydro_IsUnphysical( UNPHY_MODE_CONS, fluid, NULL, NULL_REAL, NULL_REAL, Emag,
                                      EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                      EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                      ERROR_INFO, UNPHY_VERBOSE )  )
               CheckFailed = true;

//          generic
#           else // #if ( MODEL == HYDRO )

            for (int v=0; v<FLU_NIN; v++)
            {
               bool isDataInvalid = !Aux_IsFinite( fluid[v] );

//             check for negative densities on fluid levels
#              if ( MODEL == ELBDM )
               if ( !amr->use_wave_flag[lv] )   isDataInvalid |= ( (v == DENS) && fluid[DENS] < (real)0.0 );
#              endif

               if ( isDataInvalid )
               {
                  Aux_Message( stderr, "Invalid input fluid data on level %d:\n", lv );
                  Aux_Message( stderr, "Fluid: " );
                  for (int v=0; v<FLU_NIN; v++)    Aux_Message( stderr, " [%d]=%14.7e", v, fluid[v] );
                  Aux_Message( stderr, "\n" );

                  CheckFailed = true;
               }
            }
#           endif // MODEL
         } // k, j, i


//       b. magnetic field
#        if ( MODEL == HYDRO  &&  defined MHD )
         for (int v=0; v<NCOMP_MAG; v++)
         for (int t=0; t<FLU_NXT_P1*SQR(FLU_NXT); t++)
         {
            const real B = h_Mag_Array_F_In[TID][v][t];

            if ( !Aux_IsFinite(B) )
            {
               Aux_Message( stderr, "Invalid input magnetic field: B[%d] = %14.7e\n", v, B );
               CheckFailed = true;
            }
         }
#        endif


//       c. gravitational potential
#        ifdef UNSPLIT_GRAVITY
         if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
         for (int t=0; t<CUBE(USG_NXT_F); t++)
         {
            const real pot = h_Pot_Array_USG_F[TID][t];

            if ( !Aux_IsFinite(pot) )
            {
               Aux_Message( stderr, "Invalid input gravitational potential: %14.7e\n", pot );
               CheckFailed = true;
            }
         }
#        endif
      } // for (int TID=0; TID<NPG; TID++)

      if ( CheckFailed )
         Aux_Error( ERROR_INFO, "OPT__CK_INPUT_FLUID failed (lv %d, Time %20.14e) !!\n", lv, PrepTime  );

   } // if ( OPT__CK_INPUT_FLUID )

} // FUNCTION : Flu_Prepare
