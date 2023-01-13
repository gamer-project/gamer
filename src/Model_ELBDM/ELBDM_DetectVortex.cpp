#include "GAMER.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Flag_DetectVortex
// Description :  Detect vortex by using Lap(Density)/Density
//
// Note        :  1. Flag the input cell if Lap(Density)/(3 * Density) exceeds threshold
//
// Parameter   :  i,j,k                    : Indices of the target cell in the array "Dens"
//                i_size, j_size, k_size   : Sizes of the array "Dens" in all directions
//                Dens_Array               : Density array
//                Threshodl                : Vortex detection threshold
//                Eps         : Soften factor
//
// Return      :  "true"  if vortex is detected
//                "false" if vortex is not detected
//-------------------------------------------------------------------------------------------------------
bool ELBDM_DetectVortex( const int i, const int j, const int k, const int i_size, const int j_size, const int k_size, const real Dens_Array[], const double Threshold )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= i_size  ||  j < 0  ||  j >= j_size  ||  k < 0  ||  k >= k_size  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   const int    ijk[3]    = { i, j, k };
   const int    sizes[3]  = { i_size, j_size, k_size };
   const int    Idx       = k*i_size*j_size + j*i_size + i;
   const int    dIdx[3]   = { 1, i_size, i_size*j_size };
   const real   Dens      = Dens_Array[Idx];

   int  Idx_p, Idx_c, Idx_m;
   bool isVortex;


   real Lap_Rho  = 0; 

// evaluate laplacians
   for (int d=0; d<3; d++)
   {
      if ( ijk[d] == 0 ) 
      {
         Idx_m = Idx              ;   Idx_c = Idx + dIdx[d]; Idx_p = Idx + 2 * dIdx[d];
      } else 
      if ( ijk[d] == sizes[d]-1 )
      {
         Idx_m = Idx - 2 * dIdx[d];   Idx_c = Idx - dIdx[d]; Idx_p = Idx              ;
      } else
      {
         Idx_m = Idx - 1 * dIdx[d];   Idx_c = Idx          ; Idx_p = Idx + 1 * dIdx[d];
      } 

      Lap_Rho += FABS( Dens_Array[Idx_p] - 2 * Dens_Array[Idx_c]  + Dens_Array[Idx_m]);

   } // for (int d=0; d<3; d++)


// evaluate energy density and check the flag criterion
   isVortex        = Lap_Rho / (3.0 * (Dens + 1e-12)) > Threshold;

   return isVortex;

} // FUNCTION : ELBDM_DetectVortex




#endif // #if ( MODEL == ELBDM )
