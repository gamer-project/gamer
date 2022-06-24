#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser()" and "Model_Init_StartOver_AssignData()"
//                2. Input "fluid" array stores the original values
//                3. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser()" and "Model_Init_StartOver_AssignData()"
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
bool Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time,
                           const int lv, double AuxArray[] )
{

   fluid[DENS] = 1.0e-20;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = 1.0e-22;

   return true;

} // FUNCTION : Flu_ResetByUser_Func



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser
// Description :  Reset the fluid array
//
// Note        :  1. Work for the option "OPT__RESET_FLUID == true"
//                2. Invokded by either "Flu_AdvanceDt()" or "Gra_AdvanceDt()"
//                3. Currently NOT applied to reset the input uniform array
//                   --> Init_UM() does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//
// Parameter   :  lv    :  Target refinement level
//                FluSg :  Target fluid sandglass
//                TTime :  Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser( const int lv, const int FluSg, const double TTime )
{

   const double dh       = amr->dh[lv];
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   const real   Gamma_m1 = GAMMA - (real)1.0;
   const real  _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif

   bool   Reset;
   real   fluid[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


#  pragma omp parallel for private( Reset, fluid, x, y, z, x0, y0, z0 ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//       reset this cell
         Reset = Flu_ResetByUser_Func( fluid, x, y, z, TTime, lv, NULL );

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
            fluid[DUAL] = CPU_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1 );
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
         } // if ( Reset )

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Flu_ResetByUser
