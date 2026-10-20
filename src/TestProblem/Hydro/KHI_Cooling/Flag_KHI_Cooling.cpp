#include "GAMER.h"

extern int  KHC_RefineShearMaxLv;
extern int  KHC_PeriodicYFactor;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_KHI_Cooling
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  Users can put their favorite flag criteria in this function
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_KHI_Cooling( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh          = amr->dh[lv];
   const double dy_periodic = amr->BoxSize[1] / KHC_PeriodicYFactor;
   const double y_shear     = 0.5*dy_periodic;
   const double y           = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
   const double y_periodic  = fmod( y, dy_periodic );

// flag all cells below level KHC_RefineShearMaxLv and adjacent to the shear plane
   return (  lv < KHC_RefineShearMaxLv  &&  ( fabs(y_periodic-y_shear)<=dh || y_periodic<=dh || y_periodic>=dy_periodic-dh )  );

} // FUNCTION : Flag_KHI_Cooling

