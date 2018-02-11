#include "GAMER.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Init_ByFile_AssignData
// Description :  Use the input uniform-mesh array to assign data to all patches at level "lv"
//
// Note        :  1. Work in the model ELBDM
//                2. Two modes:
//                   (1) NVar == 1: load density and set real=sqrt(dens), imag=0.0
//                   (2) NVar == 2: load both real and imaginary parts, and then *calculate* the density field
//                3. Data format in the UM_START file : [k][j][i][v]
//
// Parameter   :  lv       : Target refinement level to assign data
//                UM_Data  : Input uniform-mesh array
//                NVar     : Number of variables stored in UM_Data
//-------------------------------------------------------------------------------------------------------
void ELBDM_Init_ByFile_AssignData( const int lv, real *UM_Data, const int NVar )
{

// check
   if ( NVar != 1  &&  NVar != 2 )
      Aux_Error( ERROR_INFO, "%s only supports OPT__UM_START_NVAR == 1 or 2 !!\n", __FUNCTION__ );

#  if ( NCOMP_PASSIVE > 0 )
   Aux_Error( ERROR_INFO, "%s does NOT support passive scalars !!\n", __FUNCTION__ );
#  endif


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

         if ( NVar == 1 )
         {
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] = UM_Data[Idx];
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i] = SQRT( UM_Data[Idx] );
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i] = (real)0.0;
         }

         else // NVar == 2
         {
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i] = UM_Data[Idx+0];
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i] = UM_Data[Idx+1];
            amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] = SQR( UM_Data[Idx+0] ) + SQR( UM_Data[Idx+1] );
         }
      }}}
   }

} // FUNCTION : ELBDM_Init_ByFile_AssignData



#endif // #if ( MODEL == ELBDM )
