#include "GAMER.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Microphysics_Init
// Description :  Initialize the microphysics
//
// Note        :  1. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Microphysics_Init()
{
// check if microphysics has been initialized already
   static bool MicroPhy_Initialized = false;

   if ( MicroPhy_Initialized )  return;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

#  ifdef CR_DIFFUSION
#  ifndef COSMIC_RAY
#    error : ERROR : CR_DIFFUSION must enable COSMIC_RAY !!
#  endif

#  ifndef MHD
#    error : ERROR : CR_DIFFUSION must enable MHD !!
#  endif

   MicroPhy.CR_safety          = DT_DIFFUSION;
   MicroPhy.CR_diff_coeff_para = CR_DIFF_PARA;
   MicroPhy.CR_diff_coeff_perp = CR_DIFF_PERP;
#  endif // #ifdef CR_DIFFUSION

   MicroPhy_Initialized = true;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Microphysics_Init
