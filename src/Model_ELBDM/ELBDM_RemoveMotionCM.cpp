#include "GAMER.h"

#if ( MODEL == ELBDM )


// global variable to store the ELBDM center-of-mass velocity
double ELBDM_Vcm[3] = { NULL_REAL, NULL_REAL, NULL_REAL };




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_RemoveMotionCM
// Description :  Remove the motion of center-of-mass
//
// Note        :  1. Work with the option ELBDM_REMOVE_MOTION_CM
//                   --> Must also enable OPT__CK_CONSERVATION since it relies on Aux_Check_Conservation()
//                       to calculate the CM velocity (ELBDM_Vcm)
//                2. Invoked by Main()
//
// Parameter   :  None
//
// Return      :  amr->fluid[REAL/IMAG]
//-------------------------------------------------------------------------------------------------------
void ELBDM_RemoveMotionCM()
{

// check
   for (int d=0; d<3; d++)
   {
      if ( ELBDM_Vcm[d] == NULL_REAL )
         Aux_Error( ERROR_INFO, "ELBDM_Vcm[%d] == NULL_REAL !!\n", d );

      if ( ! Aux_IsFinite(ELBDM_Vcm[d]) )
         Aux_Error( ERROR_INFO, "ELBDM_Vcm[%d] = %14.7e !!\n", d, ELBDM_Vcm[d] );
   }


// remove the CM velocity
   for (int d=0; d<3; d++)    ELBDM_Vcm[d] *= ELBDM_ETA;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const double dh = amr->dh[lv];

#     pragma omp parallel for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         real (*const fluid)[PS1][PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
         double x, y, z, x0, y0, z0;
         real   R, I, S;

         x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

            S = ELBDM_Vcm[0]*x + ELBDM_Vcm[1]*y + ELBDM_Vcm[2]*z;

#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            if ( amr->use_wave_flag[lv] ) {
#           endif
            R = fluid[REAL][k][j][i];
            I = fluid[IMAG][k][j][i];

            fluid[REAL][k][j][i] = +R*COS(S) + I*SIN(S);
            fluid[IMAG][k][j][i] = -R*SIN(S) + I*COS(S);

#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            } else {
            fluid[PHAS][k][j][i] -= S;
            }
#           endif

         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

//    update the data on MPI buffer patches
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( amr->use_wave_flag[lv] ) {
#     endif
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _REAL|_IMAG, _NONE, Flu_ParaBuf, USELB_YES );
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      } else {
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _PHAS, _NONE, Flu_ParaBuf, USELB_YES );
      }
#     endif

   } // for (int lv=0; lv<NLEVEL; lv++)


// reset ELBDM_Vcm[] to check whether it is properly recalculated by Aux_Check_Conservation()
   for (int d=0; d<3; d++)    ELBDM_Vcm[d] = NULL_REAL;

} // FUNCTION : ELBDM_RemoveMotionCM



#endif // #if ( MODEL == ELBDM )
