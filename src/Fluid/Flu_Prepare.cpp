#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Prepare
// Description :  Prepare input arrays for the fluid solver
//
// Note        :  Invoke Prepare_PatchData()
//
// Parameter   :  lv                   : Target refinement level
//                PrepTime             : Target physical time to prepare the coarse-grid data
//                h_Flu_Array_F_In     : Host array to store the prepared fluid data
//                h_Mag_Array_F_In     : Host array to store the prepared B field (for MHD onlhy)
//                h_Pot_Array_USG_F    : Host array to store the prepared potential data (for UNSPLIT_GRAVITY only)
//                h_Corner_Array_USG_F : Host array to store the prepared corner data (for UNSPLIT_GRAVITY only)
//                NPG                  : Number of patch groups to be prepared at a time
//                PID0_List            : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Flu_Prepare( const int lv, const double PrepTime,
                  real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                  real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                  real h_Pot_Array_USG_F[][ CUBE(USG_NXT_F) ],
                  double h_Corner_Array_F[][3], const int NPG, const int *PID0_List )
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


#  if ( MODEL != HYDRO )
   const double MIN_DENS            = -1.0;  // set to an arbitrarily negative value to disable it
#  endif
#  ifndef MHD
   const int    OPT__MAG_INT_SCHEME = INT_NONE;
#  endif
   const bool   IntPhase_No         = false;
   const real   MinDens_No          = -1.0;
   const real   MinPres_No          = -1.0;
   const real   MinTemp_No          = -1.0;
   const real   MinEntr_No          = -1.0;
   const bool   DE_Consistency_Yes  = true;
   const bool   DE_Consistency_No   = false;
   const bool   DE_Consistency      = ( OPT__OPTIMIZE_AGGRESSIVE ) ? DE_Consistency_No : DE_Consistency_Yes;
   const real   MinDens             = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;


// prepare the fluid array
#  if ( MODEL == ELBDM )
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], NULL,
                      FLU_GHOST_SIZE, NPG, PID0_List, _REAL|_IMAG|_PASSIVE, _NONE,
                      OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, OPT__INT_PHASE,
                      OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
#  else
#  ifdef MHD
   real *Mag_Array = h_Mag_Array_F_In[0][0];
#  else
   real *Mag_Array = NULL;
#  endif
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], Mag_Array,
                      FLU_GHOST_SIZE, NPG, PID0_List, _TOTAL, _MAG,
                      OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, BC_POT_NONE, MinDens,    MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency );
#  endif

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


// validate input arrays for debugging purposes
   if ( OPT__CK_INPUT_FLUID )
   {
      bool CheckFailed = false;

#     pragma omp parallel for reduction ( ||:CheckFailed ) schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         real fluid[FLU_NIN];

//       a. fluid
         for (int k=0; k<FLU_NXT; k++)
         for (int j=0; j<FLU_NXT; j++)
         for (int i=0; i<FLU_NXT; i++)
         {
            const int t = IDX321( i, j, k, FLU_NXT, FLU_NXT );

            for (int v=0; v<FLU_NIN; v++)    fluid[v] = h_Flu_Array_F_In[TID][v][t];

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
               if (  !Aux_IsFinite( fluid[v] )  )
               {
                  Aux_Message( stderr, "Invalid input fluid data:\n" );
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
