#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Close
// Description :  Copy the momentum and energy density stored in the h_Flu_Array_G array back into the
//                patch pointers
//
// Note        :  1. Use SaveSg to determine where to store the data
//                   --> Currently it's set to the same Sg as the fluid data when calling
//                       Gra_AdvanceDt() in EvolveLevel()
//
// Parameter   :  lv             : Target refinement level
//                SaveSg         : Sandglass to store the updated data
//                h_Flu_Array_G  : Host array storing the updated fluid variables
//                h_DE_Array_G   : Host array storing the dual-energy status
//                h_Emag_Array_G : Host array storing the cell-centered magnetic energy (MHD with DUAL_ENERGY only)
//                NPG            : Number of patch groups to store the updated data
//                PID0_List      : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Close( const int lv, const int SaveSg, const real h_Flu_Array_G[][GRA_NIN][PS1][PS1][PS1],
                const char h_DE_Array_G[][PS1][PS1][PS1], const real h_Emag_Array_G[][PS1][PS1][PS1],
                const int NPG, const int *PID0_List )
{

   int N, PID, PID0;


#  pragma omp parallel for private( N, PID, PID0 ) schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

#        if ( MODEL == HYDRO )
//       density field is sent in and out but NOT updated in the hydro gravity solver
         for (int v=1; v<GRA_NIN; v++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
            amr->patch[SaveSg][lv][PID]->fluid[v][k][j][i] = h_Flu_Array_G[N][v][k][j][i];

//       for the dual-energy formalism only
#        ifdef DUAL_ENERGY
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
         {
//          update the dual-energy status (which is always stored in Sg=0)
            amr->patch[0][lv][PID]->de_status[k][j][i] = h_DE_Array_G[N][k][j][i];

//          correct the dual-energy variable to be consistent with the updated internal energy
//          --> only necessary for cells with the dual-energy status labelled as DE_UPDATED_BY_ETOT_GRA
//              since for all other cases we fix the internal energy in the gravity solver
#           ifdef UNSPLIT_GRAVITY
            if ( h_DE_Array_G[N][k][j][i] == DE_UPDATED_BY_ETOT_GRA )
            {
#              ifdef MHD
               const real Emag = h_Emag_Array_G[N][k][j][i];
#              else
               const real Emag = NULL_REAL;
#              endif

#              ifdef DUAL_ENERGY
               amr->patch[SaveSg][lv][PID]->fluid[DUAL][k][j][i]
                  = Hydro_Con2Dual( amr->patch[SaveSg][lv][PID]->fluid[DENS][k][j][i],
                                    amr->patch[SaveSg][lv][PID]->fluid[MOMX][k][j][i],
                                    amr->patch[SaveSg][lv][PID]->fluid[MOMY][k][j][i],
                                    amr->patch[SaveSg][lv][PID]->fluid[MOMZ][k][j][i],
                                    amr->patch[SaveSg][lv][PID]->fluid[ENGY][k][j][i],
                                    Emag, EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#              endif
            }
#           endif // #ifdef UNSPLIT_GRAVITY
         } // i,j,k
#        endif // #ifdef DUAL_ENERGY

#        elif ( MODEL == ELBDM )
         for (int v=0; v<GRA_NIN; v++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
         {
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            if ( amr->use_wave_flag[lv] ) {
#           endif
//          density field is NOT sent in and out in the ELBDM gravity solver for the wave scheme --> v+1
            amr->patch[SaveSg][lv][PID]->fluid[v+1][k][j][i] = h_Flu_Array_G[N][v][k][j][i];
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            } else {
//          in fluid scheme, send both density and phase fields --> v
//###OPTIMIZATION: no need to transfer and update the density field
            amr->patch[SaveSg][lv][PID]->fluid[v][k][j][i] = h_Flu_Array_G[N][v][k][j][i];
            }
#           endif
         }

#        else
#        error : unsupported MODEL !!
#        endif // MODEL
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Gra_Close



#endif // #ifdef GRAVITY
