#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_SwapFixUpTempArray
// Description :  Swap the flux and flux_tmp pointers (and electric and electric_tmp pointers in MHD)
//                for all patches on the target level
//
// Note        :  1. Work on both real and buffer patches
//                2. Work for the option "AUTO_REDUCE_DT"
//                   --> Flu_Close()->CorrectFlux() will store the updated (corrected) coarse-grid fluxes
//                       in the temporary array flux_tmp[] instead of flux[]
//                   --> For MHD, Flu_Close()->CorrectElectric() will store the updated (corrected) coarse-grid
//                       electric field in the temporary array electric_tmp[] instead of electric[]
//                3. Invoked by Flu_AdvanceDt()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Flu_SwapFixUpTempArray( const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv >= NLEVEL )  Aux_Error( ERROR_INFO, "incorrect lv (%d) !!\n", lv );
#  endif


// nothing to do if all related fix-up operations are disabled
#  ifndef MHD
   const bool OPT__FIXUP_ELECTRIC = false;
#  endif
   if ( !OPT__FIXUP_FLUX  &&  !OPT__FIXUP_ELECTRIC )  return;


// swap pointers
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   {
      if ( OPT__FIXUP_FLUX )
      for (int s=0; s<6; s++)
      {
         if ( amr->patch[0][lv][PID]->flux_tmp[s] != NULL )
         {
            Aux_SwapPointer( (void**)&amr->patch[0][lv][PID]->flux    [s],
                             (void**)&amr->patch[0][lv][PID]->flux_tmp[s] );
         }
      }

#     ifdef MHD
      if ( OPT__FIXUP_ELECTRIC )
      for (int s=0; s<18; s++)
      {
         if ( amr->patch[0][lv][PID]->electric_tmp[s] != NULL )
         {
            Aux_SwapPointer( (void**)&amr->patch[0][lv][PID]->electric    [s],
                             (void**)&amr->patch[0][lv][PID]->electric_tmp[s] );
         }
      }
#     endif
   } // for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)

} // FUNCTION : Flu_SwapFixUpTempArray



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_InitFixUpTempArray
// Description :  Initialize flux_tmp[] (and electric_tmp[] in MHD)
//
// Note        :  1. Copy flux[] to flux_tmp[] (and electric[] to electric_tmp[] in MHD)
//                2. Work for the option "AUTO_REDUCE_DT", for which Flu_SwapFixUpTempArray() will swap the
//                   flux and flux_tmp pointers (and electric and electric_tmp pointers in MHD)
//                3. Invoked by Flu_AdvanceDt()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Flu_InitFixUpTempArray( const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv >= NLEVEL )  Aux_Error( ERROR_INFO, "incorrect lv (%d) !!\n", lv );
#  endif


// nothing to do if all related fix-up operations are disabled
#  ifndef MHD
   const bool OPT__FIXUP_ELECTRIC = false;
#  endif
   if ( !OPT__FIXUP_FLUX  &&  !OPT__FIXUP_ELECTRIC )  return;


// copy data
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   {
      if ( OPT__FIXUP_FLUX )
      for (int s=0; s<6; s++)
      {
         if ( amr->patch[0][lv][PID]->flux_tmp[s] != NULL )
            memcpy( amr->patch[0][lv][PID]->flux_tmp[s], amr->patch[0][lv][PID]->flux[s], SQR(PS1)*NFLUX_TOTAL*sizeof(real) );
      }

#     ifdef MHD
      if ( OPT__FIXUP_ELECTRIC )
      for (int s=0; s<18; s++)
      {
         if ( amr->patch[0][lv][PID]->electric_tmp[s] != NULL )
         {
            const int Size = ( s < 6 ) ? NCOMP_ELE*PS1M1*PS1*sizeof(real) : PS1*sizeof(real);
            memcpy( amr->patch[0][lv][PID]->electric_tmp[s], amr->patch[0][lv][PID]->electric[s], Size );
         }
      }
#     endif
   } // for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)

} // FUNCTION : Flu_InitFixUpTempArray
