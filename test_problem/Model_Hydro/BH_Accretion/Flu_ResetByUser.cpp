#include "GAMER.h"

extern double BH_InBC_Rho;
extern double BH_InBC_R;
extern double BH_InBC_E;

extern double BH_SinkMass;
extern double BH_SinkMomX;
extern double BH_SinkMomY;
extern double BH_SinkMomZ;
extern double BH_SinkMomXAbs;
extern double BH_SinkMomYAbs;
extern double BH_SinkMomZAbs;
extern double BH_SinkEk;
extern double BH_SinkEt;
extern int    BH_SinkNCell;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser 
// Description :  Reset the fluid array
//
// Note        :  1. Work for the option "OPT__RESET_FLUID == true"
//                2. Invokded by both "Init_StartOver" and "EvolveLevel"
//                3. Will not be applied to reset the input uniform array
//                   --> Init_UM does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//
// Parameter   :  lv    :  Target refinement level 
//                FluSg :  Target fluid sandglass
//                TTime :  Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser( const int lv, const int FluSg, const double TTime )
{

   const double dh = amr->dh[lv];
   const real   dv = CUBE(dh);

   real   fluid[NCOMP], Ek, Et;
   double x, y, z, x0, y0, z0;
   int    SinkNCell = 0;

#  pragma omp parallel for private( fluid, Ek, Et, x, y, z, x0, y0, z0 ) schedule( runtime ) \
   reduction(+:BH_SinkMass, BH_SinkMomX, BH_SinkMomY, BH_SinkMomZ, BH_SinkMomXAbs, BH_SinkMomYAbs, BH_SinkMomZAbs, \
               BH_SinkEk, BH_SinkEt, SinkNCell)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

         for (int v=0; v<NCOMP; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//       Flu_ResetByUser_Func( fluid, x, y, z, TTime );

         const double dr[3] = { x-0.5*amr->BoxSize[0], y-0.5*amr->BoxSize[1], z-0.5*amr->BoxSize[2] };
         const double r     = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

         if ( r <= BH_InBC_R )
         {
//          record the amount of sunk variables removed at the maximum level
            if ( lv == MAX_LEVEL )
            {
               Ek = (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
               Et = fluid[ENGY] - Ek;

               BH_SinkMass    += dv*fluid[DENS];
               BH_SinkMomX    += dv*fluid[MOMX];
               BH_SinkMomY    += dv*fluid[MOMY];
               BH_SinkMomZ    += dv*fluid[MOMZ];
               BH_SinkMomXAbs += dv*FABS( fluid[MOMX] );
               BH_SinkMomYAbs += dv*FABS( fluid[MOMY] );
               BH_SinkMomZAbs += dv*FABS( fluid[MOMZ] );
               BH_SinkEk      += dv*Ek;
               BH_SinkEt      += dv*Et;
            }

            fluid[DENS] = BH_InBC_Rho;
            fluid[MOMX] = 0.0;
            fluid[MOMY] = 0.0;
            fluid[MOMZ] = 0.0;
            fluid[ENGY] = BH_InBC_E;

            SinkNCell ++;
         }

         for (int v=0; v<NCOMP; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   if ( lv == MAX_LEVEL )     BH_SinkNCell = SinkNCell; 

} // FUNCTION : Flu_ResetByUser
