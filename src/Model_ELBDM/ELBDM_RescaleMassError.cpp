#include "GAMER.h"

#if ( MODEL == ELBDM )

// global variable to store the ELBDM total mass
       double ELBDM_MassPsi     = NULL_REAL;
static double ELBDM_InitMassPsi = NULL_REAL;

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
// Return      :  amr->fluid[DENS/REAL/IMAG]
//-------------------------------------------------------------------------------------------------------
void ELBDM_RescaleMassError()
{
   if ( ELBDM_InitMassPsi == NULL_REAL )
   {
      if ( MPI_Rank == 0 )
      {
         ELBDM_InitMassPsi = ConRef[1];
      }
      MPI_Bcast( &ELBDM_InitMassPsi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }
// check

   if ( ELBDM_MassPsi == NULL_REAL )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi == NULL_REAL !!\n");

   if ( ! Aux_IsFinite(ELBDM_MassPsi) )
      Aux_Error( ERROR_INFO, "ELBDM_MassPsi = %14.7e !!\n", ELBDM_MassPsi );

// Rescale the total ELBDM mass
   for (int lv=0; lv<NLEVEL; lv++)
   {
#     pragma omp parallel for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         real (*const fluid)[PS1][PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;

         for (int k=0; k<PS1; k++)  {
         for (int j=0; j<PS1; j++)  {
         for (int i=0; i<PS1; i++)  {

#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            if ( amr->use_wave_flag[lv] ) {
#           endif
    
               fluid[REAL][k][j][i] *= SQRT(ELBDM_InitMassPsi/ELBDM_MassPsi);
               fluid[IMAG][k][j][i] *= SQRT(ELBDM_InitMassPsi/ELBDM_MassPsi);
               fluid[DENS][k][j][i]  = SQR(fluid[REAL][k][j][i]) + SQR(fluid[IMAG][k][j][i]);

#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            } else {
               fluid[DENS][k][j][i] *= (ELBDM_InitMassPsi/ELBDM_MassPsi);
            }
#           endif

         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

//    update the data on MPI buffer patches
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( amr->use_wave_flag[lv] ) {
#     endif
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _REAL|_IMAG|_DENS, _NONE, Flu_ParaBuf, USELB_YES );
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      } else {
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _DENS, _NONE, Flu_ParaBuf, USELB_YES );
      }
#     endif

   } // for (int lv=0; lv<NLEVEL; lv++)


// reset ELBDM_MassPsi[] to check whether it is properly recalculated by Aux_Check_Conservation()
ELBDM_MassPsi = NULL_REAL;

} // FUNCTION : ELBDM_RescaleMassError



#endif // #if ( MODEL == ELBDM )
