#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_LB_ResetBufferElectric
// Description :  Reset all electric field in the buffer patches as zero
//
// Note        :  1. Invoked by Flu_AdvanceDt()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void MHD_LB_ResetBufferElectric( const int lv )
{

// check
   if ( !amr->WithElectric )
   {
      Aux_Message( stderr, "WARNING : invoking %s is useless since no electric field is required !!\n", __FUNCTION__ );
      return;
   }


#  pragma omp parallel for schedule( runtime )
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][27]; PID++)
   {
      for (int s=0; s<18; s++)
      {
         const int EleSize = ( s < 6 ) ? NCOMP_ELE*PS1M1*PS1 : PS1;
         real *ElePtr = amr->patch[0][lv][PID]->electric[s];

         if ( ElePtr != NULL )
            for (int t=0; t<EleSize; t++)    ElePtr[t] = (real)0.0;
      }
   }

} // FUNCTION : MHD_LB_ResetBufferElectric



#endif // #if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )
