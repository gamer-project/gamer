#include "GAMER.h"


#if ( MODEL == SR_HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_Flag_Lorentz
// Description :  Check if the Lorentz factr of the input data at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Flag if Lorentz factor > threshold
//
// Parameter   :  i,j,k       : Target array indices
//                lv          : Target refinement level
//                PID         : Target patch ID
//                Threshold   : Refinement threshold
//
// Return      :  "true"  if the Lorentz factor is greater          than the given threshold
//                "false" if the Lorentz factor is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool SRHydro_Flag_Lorentz( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{
// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0 ||  j >= PS1  ||  k < 0  ||  k >= PS1   )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif


   const real (*Dens)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];
   const real (*MomX)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX];
   const real (*MomY)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY];
   const real (*MomZ)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ];
   const real (*Engy)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY];

   bool Flag = false;

   real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID];                                                                                            
	           
	           
   Cons[DENS] = Dens[k][j][i];
   Cons[MOMX] = MomX[k][j][i];
   Cons[MOMY] = MomY[k][j][i];
   Cons[MOMZ] = MomZ[k][j][i];
   Cons[ENGY] = Engy[k][j][i];


// calculate the Lorentz factor
   double LorentzFactor;

   LorentzFactor = SRHydro_Con2Pri( Cons, Prim, GAMMA, MIN_TEMP );

// flag if the Lorentz factor > threshold
   Flag = ( LorentzFactor >= Threshold );

   return Flag;

} // FUNCTION : SRHydro_Flag_Lorentz



#endif // #if ( MODEL == SR_HYDRO )
