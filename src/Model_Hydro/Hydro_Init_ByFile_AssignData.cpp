#include "GAMER.h"

#if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Init_ByFile_AssignData
// Description :  Use the input uniform-mesh array to assign data to all patches at level "lv"
//
// Note        :  1. Work in the model HYDRO
//                2. Only load "density". Momentum x/y/z are initialized as zero. Total energy is initialized
//                   by specifying the sound speed parameter "Cs".
//                3. Data format in the UM_IC file: [k][j][i][v]
//
// Parameter   :  lv       : Target refinement level to assign data
//                UM_Data  : Input uniform-mesh array
//                NVar     : Number of variables stored in UM_Data
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_ByFile_AssignData( const int lv, real *UM_Data, const int NVar )
{

// check
   if ( NVar != 1 )
      Aux_Error( ERROR_INFO, "%s only supports OPT__UM_IC_NVAR == 1 (for total gas density) !!\n", __FUNCTION__ );

#  if ( NCOMP_PASSIVE > 0 )
   Aux_Error( ERROR_INFO, "%s does NOT support passive scalars !!\n", __FUNCTION__ );
#  endif


#  if ( defined COMOVING  &&  defined GRAVITY )
   const real  a_z1000    = 1.0/1001.0;
   const real  T_z1000    = 3000.0;
   const real  T          = T_z1000*pow( a_z1000/A_INIT, 2.0 );
   const real  Boltzmann  = 1.381e-29;   // (kg*km*km/s/s/K)
   const real  ProtonMass = 1.673e-27;   // (kg)
   const real  Hubble0    = 100.0;       // (h*km/s/Mpc)
   const real  V_unit     = 1.0*Hubble0; // (h*km/s)

   const real  Cs         = A_INIT * sqrt( GAMMA*Boltzmann*T/ProtonMass ) / V_unit;    // sound speed

   static bool FirstTime  = true;

   if ( MPI_Rank == 0  &&  FirstTime )
   {
      Aux_Message( stdout, "NOTE : sound speed is set to %13.7e in the cosmological simulations\n", Cs );
      FirstTime = false;
   }
#  else
   const real Cs          = 1.0;
#  endif // #if ( defined COMOVING  &&  defined GRAVITY )


   const int NX[3]  = { NX0[0]*(1<<lv), NX0[1]*(1<<lv), NX0[2]*(1<<lv) };
   const int scale0 = amr->scale[ 0];
   const int scale  = amr->scale[lv];

   int *Corner;
   int ii, jj, kk;
   long Idx;


   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      Corner = amr->patch[0][lv][PID]->corner;

      for (int k=0; k<PATCH_SIZE; k++)    { kk = ( Corner[2] - MPI_Rank_X[2]*NX0[2]*scale0 ) / scale + k;
      for (int j=0; j<PATCH_SIZE; j++)    { jj = ( Corner[1] - MPI_Rank_X[1]*NX0[1]*scale0 ) / scale + j;
      for (int i=0; i<PATCH_SIZE; i++)    { ii = ( Corner[0] - MPI_Rank_X[0]*NX0[0]*scale0 ) / scale + i;

         Idx = (long)NVar*( (long)kk*NX[1]*NX[0] + jj*NX[0]+ ii );

//       assuming that UM_Data only stores density, momentum == 0, and sound speed = Cs
         if ( NVar == 1 )
         {
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] = UM_Data[Idx];
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i] = 0.0;
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i] = 0.0;
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i] = 0.0;
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i] = UM_Data[Idx]*Cs*Cs/( GAMMA*(GAMMA-1.0) );
         }

         /* for arbitrary number of input variables (NVar == 3 in the following example)
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] = UM_Data[Idx+0];
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i] = 0.0;
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i] = UM_Data[Idx+1];
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i] = 0.0;
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i] = UM_Data[Idx+2];
         */
      }}}
   }

} // FUNCTION : Hydro_Init_ByFile_AssignData



#endif // #if ( MODEL == HYDRO )
