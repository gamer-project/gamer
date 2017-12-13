#include "GAMER.h"

#if ( MODEL == ELBDM )

static void Init_Function_User( real fluid[], const real x, const real y, const real z, const double Time );
void (*Init_Function_Ptr)( real fluid[], const real x, const real y, const real z, const double Time ) = Init_Function_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User 
// Description :  Function to initialize the scalar field
//
// Note        :  Invoked by "ELBDM_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Scalar field array to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User( real fluid[], const real x, const real y, const real z, const double Time )
{

   const real C1[3]   = { 0.5*amr->BoxSize[0]+100.0, 
                          0.5*amr->BoxSize[1]+200.0, 
                          0.5*amr->BoxSize[2]+300.0 };
   const real C2[3]   = { 20.0, 40.0, 10.0 };
   const real Height1 =   5.0;
   const real Height2 =   8.0;
   const real Width1  =  64.0;
   const real Width2  = 512.0;

   fluid[REAL] = 1.0 + Height1*exp(  -( SQR(x-C1[0]) + SQR(y-C1[1]) + SQR(z-C1[2]) ) /SQR(Width1)  );
   fluid[IMAG] = 1.0 + Height2*exp(  -( SQR(x-C2[0]) + SQR(y-C2[1]) + SQR(z-C2[2]) ) /SQR(Width2)  );
   fluid[DENS] = fluid[REAL]*fluid[REAL] + fluid[IMAG]*fluid[IMAG];

} // FUNCTION : Init_Function_User



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Init_StartOver_AssignData
// Description :  Construct the initial condition in ELBDM
//
// Note        :  1. Work for the option "OPT__INIT == INIT_STARTOVER"
//                2. The initialization function should be specified in "Init_Function_Ptr", which is a function
//                   pointer pointing to either "Init_Function_User" or the test problem specified function
//                   (e.g., ELBDM_TestProbSol_JeansInstability_Comoving) set in "Init_TestProb"
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void ELBDM_Init_StartOver_AssignData( const int lv )
{

   const real scale  = (real)amr->scale[lv];
   const real dh     = amr->dh[lv];
   const int  NSub   = INIT_SUBSAMPLING_NCELL;
   const real dh_sub = dh / NSub;
   const real _NSub3 = 1.0/(NSub*NSub*NSub);

   real fluid[NCOMP+1], fluid_sub[NCOMP+1];
   real x, y, z, x0, y0, z0;


   if ( NSub > 1 )   // with sub-sampling
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)  {  z0 = ( amr->patch[0][lv][PID]->corner[2]/scale + k )*dh;
      for (int j=0; j<PS1; j++)  {  y0 = ( amr->patch[0][lv][PID]->corner[1]/scale + j )*dh;
      for (int i=0; i<PS1; i++)  {  x0 = ( amr->patch[0][lv][PID]->corner[0]/scale + i )*dh;

         for (int v=0; v<NCOMP; v++)   fluid[v] = 0.0;

         for (int kk=0; kk<NSub; kk++)    {  z = z0 + (kk+0.5)*dh_sub;
         for (int jj=0; jj<NSub; jj++)    {  y = y0 + (jj+0.5)*dh_sub;
         for (int ii=0; ii<NSub; ii++)    {  x = x0 + (ii+0.5)*dh_sub;

            Init_Function_Ptr( fluid_sub, x, y, z, Time[lv] );

            for (int v=0; v<NCOMP; v++)   fluid[v] += fluid_sub[v];

         }}}

         for (int v=0; v<NCOMP; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v]*_NSub3;

      }}}
   } // if ( NSub > 1 )

   else // without sub-sampling
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)  {  z = ( amr->patch[0][lv][PID]->corner[2]/scale + k + 0.5 )*dh;
      for (int j=0; j<PS1; j++)  {  y = ( amr->patch[0][lv][PID]->corner[1]/scale + j + 0.5 )*dh;
      for (int i=0; i<PS1; i++)  {  x = ( amr->patch[0][lv][PID]->corner[0]/scale + i + 0.5 )*dh;

         Init_Function_Ptr( fluid, x, y, z, Time[lv] );

//       for (int v=0; v<NCOMP; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];
         for (int v=0; v<NCOMP+1; v++) amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   } // if ( NSub > 1 ) ... else ...

} // FUNCTION : ELBDM_Init_StartOver_AssignData



#endif // #if ( MODEL == ELBDM )
