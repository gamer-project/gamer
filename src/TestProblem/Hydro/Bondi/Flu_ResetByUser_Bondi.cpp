#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )



extern double Bondi_InBC_Rho;
extern double Bondi_InBC_R;
extern double Bondi_InBC_E;

extern double Bondi_SinkMass;
extern double Bondi_SinkMomX;
extern double Bondi_SinkMomY;
extern double Bondi_SinkMomZ;
extern double Bondi_SinkMomXAbs;
extern double Bondi_SinkMomYAbs;
extern double Bondi_SinkMomZAbs;
extern double Bondi_SinkEk;
extern double Bondi_SinkEt;
extern int    Bondi_SinkNCell;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func_Bondi
// Description :  Function to reset the fluid field in the Bondi accretion problem
//
// Note        :  1. Invoked by "Flu_ResetByUser_API_Bondi()" and "Model_Init_StartOver_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                   --> This function pointer is reset by "Init_TestProb_Bondi()"
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_StartOver_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API_Bondi()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API_Bondi()" and "Model_Init_StartOver_AssignData()"
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser_Func_Bondi( real fluid[], const double x, const double y, const double z, const double Time,
                                 const int lv, double AuxArray[] )
{

   const double Pos[3]  = { x, y, z };
   const double InBC_R2 = SQR( Bondi_InBC_R );
   double dr2[3], r2;

   for (int d=0; d<3; d++)
   {
      dr2[d] = Pos[d] - 0.5*amr->BoxSize[d];
      dr2[d] = SQR( dr2[d] );

//    skip cells far away from the void region
      if ( dr2[d] > InBC_R2 )    return false;
   }

   r2 = dr2[0] + dr2[1] + dr2[2];

   if ( r2 <= InBC_R2 )
   {
      fluid[DENS] = Bondi_InBC_Rho;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = Bondi_InBC_E;

      return true;
   }

   else
      return false;

} // FUNCTION : Flu_ResetByUser_Func_Bondi



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API_Bondi
// Description :  API for resetting the fluid array in the Bondi accretion problem
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked by either "Flu_AdvanceDt()" or "Gra_AdvanceDt()" using the function pointer
//                   "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer is reset by "Init_TestProb_Bondi()"
//                3. Currently does not work with "OPT__OVERLAP_MPI"
//                4. Invoke "Flu_ResetByUser_Func_Bondi" directly
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Target fluid sandglass
//                TTime : Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_Bondi( const int lv, const int FluSg, const double TTime )
{

   const double dh       = amr->dh[lv];
   const real   dv       = CUBE(dh);
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   const real   Gamma_m1 = GAMMA - (real)1.0;
   const real  _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif

   bool   Reset;
   real   fluid[NCOMP_TOTAL], fluid_bk[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;

// reset to 0 since we only want to record the number of void cells **for one sub-step**
   Bondi_SinkNCell = 0;


#  pragma omp parallel for private( Reset, fluid, fluid_bk, x, y, z, x0, y0, z0 ) schedule( runtime ) \
   reduction(+:Bondi_SinkMass, Bondi_SinkMomX, Bondi_SinkMomY, Bondi_SinkMomZ, Bondi_SinkMomXAbs, Bondi_SinkMomYAbs, Bondi_SinkMomZAbs, \
               Bondi_SinkEk, Bondi_SinkEt, Bondi_SinkNCell)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            fluid   [v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//          backup the unmodified values since we want to record the amount of sunk variables removed at the maximum level
            fluid_bk[v] = fluid[v];
         }

//       reset this cell
         Reset = Flu_ResetByUser_Func_Bondi( fluid, x, y, z, TTime, lv, NULL );

//       operations necessary only when this cell has been reset
         if ( Reset )
         {
#           if ( MODEL == HYDRO  ||  MODEL == MHD )
//          check minimum density and pressure
            fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
            fluid[ENGY] = CPU_CheckMinPresInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                  Gamma_m1, _Gamma_m1, MIN_PRES );

//          calculate the dual-energy variable (entropy or internal energy)
#           if   ( DUAL_ENERGY == DE_ENPY )
            fluid[ENPY] = CPU_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1 );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               CPU_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif
#           endif // if ( MODEL == HYDRO  ||  MODEL == MHD )

//          store the reset values
            for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];


//          record the amount of sunk variables removed at the maximum level
            if ( lv == MAX_LEVEL )
            {
               real Ek = (real)0.5*( SQR(fluid_bk[MOMX]) + SQR(fluid_bk[MOMY]) + SQR(fluid_bk[MOMZ]) ) / fluid_bk[DENS];
               real Et = fluid_bk[ENGY] - Ek;

               Bondi_SinkMass    += dv*fluid_bk[DENS];
               Bondi_SinkMomX    += dv*fluid_bk[MOMX];
               Bondi_SinkMomY    += dv*fluid_bk[MOMY];
               Bondi_SinkMomZ    += dv*fluid_bk[MOMZ];
               Bondi_SinkMomXAbs += dv*FABS( fluid_bk[MOMX] );
               Bondi_SinkMomYAbs += dv*FABS( fluid_bk[MOMY] );
               Bondi_SinkMomZAbs += dv*FABS( fluid_bk[MOMZ] );
               Bondi_SinkEk      += dv*Ek;
               Bondi_SinkEt      += dv*Et;
               Bondi_SinkNCell   ++;
            }
         } // if ( Reset )

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Flu_ResetByUser_API_Bondi



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
