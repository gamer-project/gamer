#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )

extern double Bondi_MassBH;
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

extern bool   Bondi_void;
extern bool   Bondi_dynBH;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func_Bondi
// Description :  Function to reset the fluid field in the Bondi accretion problem
//
// Note        :  1. Invoked by Flu_ResetByUser_API_Bondi() and Hydro_Init_ByFunction_AssignData() using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                   --> This function pointer is reset by Init_TestProb_Hydro_Bondi()
//                   --> Hydro_Init_ByFunction_AssignData(): constructing initial condition
//                       Flu_ResetByUser_API_Bondi()       : after each update
//                2. Input fluid[] stores the original values
//                3. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically
//                4. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                Emag     : Magnetic energy (MHD only)
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                dt       : Time interval to advance solution
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
int Flu_ResetByUser_Func_Bondi( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                const double dt, const int lv, double AuxArray[] )
{

   if ( !Bondi_void )   return false;


   const double Pos[3]  = { x, y, z };
   const double InBC_R2 = SQR( Bondi_InBC_R );
   double dr2[3], r2;

   for (int d=0; d<3; d++)
   {
      dr2[d] = SQR( Pos[d] - amr->BoxCenter[d] );

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
//                2. Invoked using the function pointer "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer is reset by Init_TestProb_Hydro_Bondi()
//                3. Currently does not work with "OPT__OVERLAP_MPI"
//                4. Invoke Flu_ResetByUser_Func_Bondi() directly
//
// Parameter   :  lv      : Target refinement level
//                FluSg   : Target fluid sandglass
//                MagSg   : Target B field sandglass
//                TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_Bondi( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt )
{

   const double dh       = amr->dh[lv];
   const real   dv       = CUBE(dh);
#  if ( MODEL == HYDRO )
   const real   Gamma_m1 = GAMMA - (real)1.0;
   const real  _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif

   int    Reset;
   real   fluid[NCOMP_TOTAL], fluid_bk[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;
   double SinkMass_OneSubStep_ThisRank = 0.0;   // variables to record sink mass at every time step
   double SinkMass_OneSubStep_AllRank;

// reset to 0 since we only want to record the number of void cells **for one sub-step**
   Bondi_SinkNCell = 0;


#  pragma omp parallel for private( Reset, fluid, fluid_bk, x, y, z, x0, y0, z0 ) schedule( runtime ) \
   reduction(+:Bondi_SinkMass, Bondi_SinkMomX, Bondi_SinkMomY, Bondi_SinkMomZ, Bondi_SinkMomXAbs, Bondi_SinkMomYAbs, Bondi_SinkMomZAbs, \
               Bondi_SinkEk, Bondi_SinkEt, Bondi_SinkNCell, SinkMass_OneSubStep_ThisRank )
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

#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#        else
         const real Emag = (real)0.0;
#        endif

//       reset this cell
         Reset = Flu_ResetByUser_Func_Bondi( fluid, Emag, x, y, z, TimeNew, dt, lv, NULL );

//       operations necessary only when this cell has been reset
         if ( Reset )
         {
//          apply density and energy floors
#           if ( MODEL == HYDRO )

            fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
            fluid[ENGY] = Hydro_CheckMinEintInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                    (real)MIN_EINT, Emag );

//          calculate the dual-energy variable
#           ifdef DUAL_ENERGY
            fluid[DUAL] = Hydro_Con2Dual( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Emag,
                                          EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#           endif

//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif
#           endif // if ( MODEL == HYDRO )

//          store the reset values
            for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];


//          record the amount of sunk variables removed at the maximum level
            if ( lv == MAX_LEVEL )
            {
               real Ek = (real)0.5*( SQR(fluid_bk[MOMX]) + SQR(fluid_bk[MOMY]) + SQR(fluid_bk[MOMZ]) ) / fluid_bk[DENS];
               real Et = fluid_bk[ENGY] - Ek - Emag;

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

               SinkMass_OneSubStep_ThisRank += dv*fluid_bk[DENS];
            }

            else if ( amr->patch[0][lv][PID]->son == -1 )
            {
//             void region must be completely refined to the max level
               Aux_Error( ERROR_INFO, "void region lies outside the max-level region !!\n" );
            }
         } // if ( Reset )
      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   if ( Bondi_dynBH )
   {
      MPI_Allreduce( &SinkMass_OneSubStep_ThisRank, &SinkMass_OneSubStep_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      Bondi_MassBH += SinkMass_OneSubStep_AllRank;
   }

} // FUNCTION : Flu_ResetByUser_API_Bondi



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
