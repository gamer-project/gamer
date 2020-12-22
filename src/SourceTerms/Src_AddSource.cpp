#include "GAMER.h"

// function pointer for the user-defined source terms
extern void (*Src_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                             const int lv, double AuxArray[], const double dt );




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_AddSource
// Description :  Add various local source terms
//
// Note        :  1. Invoked by Evolve()
//                2. Grackle library is treated separately
//
// Parameter   :  lv      : Target refinement level
//                TimeNew : Target physical time to reach
//                TimeOld : Physical time before update
//                          --> This function updates physical time from TimeOld to TimeNew
//                dt      : Time interval to advance solution
//                FluSg   : Target fluid sandglass (for both input and output)
//
// Return      : fluid[] in all patches
//-------------------------------------------------------------------------------------------------------
void Src_AddSource( const int lv, const double TimeNew, const double TimeOld, const double dt, const int FluSg )
{

// check
   if ( SRC_USER  &&  Src_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "Src_User_Ptr == NULL !!\n" );


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

//       get the input array
         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//###REVISE: should we make different source terms "commutative"?
//       add source terms
//       (1) user-defined
         if ( SRC_USER )   Src_User_Ptr( fluid, x, y, z, TimeNew, lv, NULL, dt );

//       store the updated results
         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Src_AddSource
