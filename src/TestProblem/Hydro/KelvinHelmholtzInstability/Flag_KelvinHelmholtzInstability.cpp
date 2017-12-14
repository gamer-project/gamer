#include "GAMER.h"

extern bool KH_AllRankSame;
extern int  KH_RefineShearMaxLv;
extern int  KH_PeriodicZFactor;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_KelvinHelmholtzInstability
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  Users can put their favorite flag criteria in this function
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the
//                              file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_KelvinHelmholtzInstability( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{

   const double dh          = amr->dh[lv];
   const double dz_periodic = ( KH_AllRankSame ) ? amr->BoxSize[2] / KH_PeriodicZFactor / MPI_NRank_X[2]
                                                 : amr->BoxSize[2] / KH_PeriodicZFactor;
   const double z_shear     = 0.5*dz_periodic;
   const double z           = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;
   const double z_periodic  = fmod( z, dz_periodic );

// flag all cells below level KH_RefineShearMaxLv and adjacent to the shear plane
   return (  lv < KH_RefineShearMaxLv  &&  ( fabs(z_periodic-z_shear)<=dh || z_periodic<=dh || z_periodic>=dz_periodic-dh )  );

} // FUNCTION : Flag_KelvinHelmholtzInstability

