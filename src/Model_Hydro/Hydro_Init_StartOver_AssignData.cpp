#include "GAMER.h"

#if ( MODEL == HYDRO )

static void Init_Function_User( real fluid[], const double x, const double y, const double z, const double Time );
void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time ) = Init_Function_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User 
// Description :  Function to initialize the fluid field 
//
// Note        :  Invoked by "Hydro_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double Gamma2  = 1.0/GAMMA/(GAMMA-1.0);
   const double C1[3] = { 0.5*amr->BoxSize[0]+100.0, 
                          0.5*amr->BoxSize[1]+200.0, 
                          0.5*amr->BoxSize[2]+300.0 };
   const double C2[3] = { 20.0, 40.0, 10.0 };

   const double Cs      =   1.0;
   const double Height1 = 100.0;
   const double Height2 = 400.0;
   const double Width1  = 640.0;
   const double Width2  = 512.0;


   fluid[DENS] = 1.0 + Height1*exp(  -( SQR(x-C1[0])+ SQR(y-C1[1]) + SQR(z-C1[2]) ) / SQR(Width1)  );
   fluid[DENS] +=      Height2*exp(  -( SQR(x-C2[0])+ SQR(y-C2[1]) + SQR(z-C2[2]) ) / SQR(Width2)  );
   fluid[MOMX] = 1.0;
   fluid[MOMY] = 2.0;
   fluid[MOMZ] = 3.0;
   fluid[ENGY] = Cs*Cs*fluid[DENS]*Gamma2 + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : Init_Function_User



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Init_StartOver_AssignData
// Description :  Construct the initial condition in HYDRO
//
// Note        :  1. Work for the option "OPT__INIT == INIT_STARTOVER"
//                2. The initialization function should be specified in "Init_Function_Ptr", which is a function
//                   pointer pointing to either "Init_Function_User" or the test problem specified function
//                   (e.g., Hydro_TestProbSol_Riemann) set in "Init_TestProb"
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_StartOver_AssignData( const int lv )
{

   const int    NSub   = INIT_SUBSAMPLING_NCELL;
   const double dh     = amr->dh[lv];
   const double dh_sub = dh / NSub;
   const double _NSub3 = 1.0/(NSub*NSub*NSub);

   real   fluid[NCOMP], fluid_sub[NCOMP];
   double x, y, z, x0, y0, z0;


   if ( NSub > 1 )   // with sub-sampling
   {
#     pragma omp parallel for private( fluid, fluid_sub, x, y, z, x0, y0, z0 ) schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)  {  z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + 0.5*dh_sub;
      for (int j=0; j<PS1; j++)  {  y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + 0.5*dh_sub;
      for (int i=0; i<PS1; i++)  {  x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + 0.5*dh_sub;

         for (int v=0; v<NCOMP; v++)   fluid[v] = 0.0;

         for (int kk=0; kk<NSub; kk++)    {  z = z0 + kk*dh_sub;
         for (int jj=0; jj<NSub; jj++)    {  y = y0 + jj*dh_sub;
         for (int ii=0; ii<NSub; ii++)    {  x = x0 + ii*dh_sub;

            Init_Function_Ptr( fluid_sub, x, y, z, Time[lv] );

            for (int v=0; v<NCOMP; v++)   fluid[v] += fluid_sub[v];

         }}}

         for (int v=0; v<NCOMP; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v]*_NSub3;

      }}}
   } // if ( NSub > 1 )

   else // without sub-sampling
   {
#     pragma omp parallel for private( fluid, x, y, z ) schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)  {  z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;
      for (int j=0; j<PS1; j++)  {  y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
      for (int i=0; i<PS1; i++)  {  x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh;

         Init_Function_Ptr( fluid, x, y, z, Time[lv] );

         for (int v=0; v<NCOMP; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   } // if ( NSub > 1 ) ... else ...

} // FUNCTION : Hydro_Init_StartOver_AssignData



#endif // #if ( MODEL == HYDRO )
