#include "GAMER.h"

#if ( MODEL == ELBDM )


// global variable to store the ELBDM total mass
double ELBDM_MassPsi = NULL_REAL;
double ELBDM_MassPsi_AErr = NULL_REAL;
double ELBDM_MassPsi_RErr = NULL_REAL;



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_RescaleMassError
// Description :  Remove the mass error created bo floating numerical error, only for base level.
//
// Note        :  1. Work with the option ELBDM_RESCALE_MASS_ERROR
//                   --> Must also enable OPT__CK_CONSERVATION since it relies on Aux_Check_Conservation()
//                       to calculate the total ELBDM mass (ELBDM_MassPsi)
//                2. Invoked by Main()
//
// Parameter   :  None
//
// Return      :  amr->fluid[REAL/IMAG]
//-------------------------------------------------------------------------------------------------------
void ELBDM_RescaleMassError()
{

// check

   if ( ELBDM_MassPsi == NULL_REAL )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi == NULL_REAL !!\n");

   if ( ! Aux_IsFinite(ELBDM_MassPsi) )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi = %14.7e !!\n", ELBDM_MassPsi );

   if ( ELBDM_MassPsi_AErr == NULL_REAL )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi_AErr == NULL_REAL !!\n");

   if ( ! Aux_IsFinite(ELBDM_MassPsi_AErr) )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi_AErr = %14.7e !!\n", ELBDM_MassPsi_AErr );

   if ( ELBDM_MassPsi_RErr == NULL_REAL )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi_RErr == NULL_REAL !!\n");

   if ( ! Aux_IsFinite(ELBDM_MassPsi_RErr) )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi_RErr = %14.7e !!\n", ELBDM_MassPsi_RErr );

// Rescale the total ELBDM mass
   for (int lv=0; lv<NLEVEL; lv++)
   {
      const double dh = amr->dh[lv];

#     pragma omp parallel for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         real (*const fluid)[PS1][PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
         double x, y, z, x0, y0, z0;

         x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

            fluid[REAL][k][j][i] /= (1-ELBDM_MassPsi_RErr);
            fluid[IMAG][k][j][i] /= (1-ELBDM_MassPsi_RErr);

         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

//    update the data on MPI buffer patches
#     if ( amr->use_wave_flag[lv] ) {
#     endif
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _REAL|_IMAG, _NONE, Flu_ParaBuf, USELB_YES );
      } else {
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _PHAS, _NONE, Flu_ParaBuf, USELB_YES );
      }

   } // for (int lv=0; lv<NLEVEL; lv++)


// reset ELBDM_Vcm[] to check whether it is properly recalculated by Aux_Check_Conservation()
ELBDM_MassPsi = NULL_REAL;
ELBDM_MassPsi_AErr = NULL_REAL;
ELBDM_MassPsi_RErr = NULL_REAL;

} // FUNCTION : ELBDM_RescaleMassError



#endif // #if ( MODEL == ELBDM )
