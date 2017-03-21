#include "GAMER.h"

static void Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func 
// Description :  Function to reset the fluid field 
//
// Note        :  1. Invoked by "Flu_ResetByUser"
//                2. Input "fluid" array stores the original values
//
// Parameter   :  fluid : Fluid array storing both the input (origial) and reset values
//                        --> Including both active and passive variables
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time )
{

// Example : reset fluid variables to extremely small values if the cell is within a specific sphere
   /*
   const real dr[3]   = { x-0.5*amr->BoxSize[0], y-0.5*amr->BoxSize[1], z-0.5*amr->BoxSize[2] };
   const real r       = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

   const real TRad    = 0.3;
   const real MaxDens = 1.0e15;
   const real MaxPres = 1.0e15;

   if ( r <= TRad )
   {
//    set active scalars
      fluid[DENS] = MaxDens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = MaxPres / ( GAMMA-(real)1.0 );

//    set passive scalars
   }
   */

} // FUNCTION : Flu_ResetByUser_Func



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

   real   fluid[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


#  pragma omp parallel for private( fluid, x, y, z, x0, y0, z0 ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

         Flu_ResetByUser_Func( fluid, x, y, z, TTime );

         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   }

} // FUNCTION : Flu_ResetByUser
