#include "GAMER.h"


#if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Flag_Vorticity
// Description :  Check if the vorticity of the input data at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Flag if curl(v)*dh/|v| > threshold
//                2. For cells adjacent to the patch boundary, only first-order approximation is adopted
//                   to estimate derivatives. Otherwise, second-order approximation is adopted.
//                   --> Do NOT need to prepare the ghost-zone data for the target patch
//
// Parameter   :  i,j,k       : Target array indices
//                lv          : Target refinement level
//                PID         : Target patch ID
//                Threshold   : Refinement threshold
//
// Return      :  "true"  if the vorticity is greater          than the given threshold
//                "false" if the vorticity is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool Hydro_Flag_Vorticity( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0 ||  j >= PS1  ||  k < 0  ||  k >= PS1   )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif


   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];
   const real (*MomX)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX];
   const real (*MomY)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY];
   const real (*MomZ)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ];

   int  im, ip, jm, jp, km, kp;
   real _dx, _dy, _dz, _Rho, Vx, Vy, Vz, V2, Wx, Wy, Wz, W2;
   bool Flag = false;


// check if the target cell indices are adjacent to the patch boundaries
   if      ( i == 0     ) { im = 0;       ip = 1;       _dx = (real)1.0; }
   else if ( i == PS1-1 ) { im = PS1-2;   ip = PS1-1;   _dx = (real)1.0; }
   else                   { im = i-1;     ip = i+1;     _dx = (real)0.5; }

   if      ( j == 0     ) { jm = 0;       jp = 1;       _dy = (real)1.0; }
   else if ( j == PS1-1 ) { jm = PS1-2;   jp = PS1-1;   _dy = (real)1.0; }
   else                   { jm = j-1;     jp = j+1;     _dy = (real)0.5; }

   if      ( k == 0     ) { km = 0;       kp = 1;       _dz = (real)1.0; }
   else if ( k == PS1-1 ) { km = PS1-2;   kp = PS1-1;   _dz = (real)1.0; }
   else                   { km = k-1;     kp = k+1;     _dz = (real)0.5; }


// calculate velocity magnitude
   _Rho = (real)1.0/Rho[k][j][i];
   Vx   = _Rho*MomX[k][j][i];
   Vy   = _Rho*MomY[k][j][i];
   Vz   = _Rho*MomZ[k][j][i];
   V2   = Vx*Vx + Vy*Vy + Vz*Vz;


// calculate vorticity
   Wx =   _dy*( MomZ[k ][jp][i ] / Rho[k ][jp][i ] - MomZ[k ][jm][i ] / Rho[k ][jm][i ] )
        - _dz*( MomY[kp][j ][i ] / Rho[kp][j ][i ] - MomY[km][j ][i ] / Rho[km][j ][i ] );
   Wy =   _dz*( MomX[kp][j ][i ] / Rho[kp][j ][i ] - MomX[km][j ][i ] / Rho[km][j ][i ] )
        - _dx*( MomZ[k ][j ][ip] / Rho[k ][j ][ip] - MomZ[k ][j ][im] / Rho[k ][j ][im] );
   Wz =   _dx*( MomY[k ][j ][ip] / Rho[k ][j ][ip] - MomY[k ][j ][im] / Rho[k ][j ][im] )
        - _dy*( MomX[k ][jp][i ] / Rho[k ][jp][i ] - MomX[k ][jm][i ] / Rho[k ][jm][i ] );
   W2   = Wx*Wx + Wy*Wy + Wz*Wz;


// flag if curl(v)*dh/|v| > threshold
   Flag = ( W2/V2 > SQR(Threshold) );

   return Flag;

} // FUNCTION : Hydro_Flag_Vorticity



#endif // #if ( MODEL == HYDRO )
