#include "GAMER.h"

// declare as static so that other functions cannot invoke them directly and must use the function pointers
static bool Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time,
                                  const int lv, double AuxArray[] );
static void Flu_ResetByUser_API( const int lv, const int FluSg, const double TTime );

// these function pointers may be overwritten by various test problem initializers
bool (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                  const int lv, double AuxArray[] ) = Flu_ResetByUser_Func;
void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const double TTime ) = Flu_ResetByUser_API;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_ByFunction_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()"
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
bool Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time,
                           const int lv, double AuxArray[] )
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

      return true;
   }

   else
      return false;
   */

   return false;

} // FUNCTION : Flu_ResetByUser_Func



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API
// Description :  API for resetting the fluid array
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked by either "Flu_AdvanceDt()" or "Gra_AdvanceDt()" using the function pointer
//                   "Flu_ResetByUser_API_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                3. Currently NOT applied to the input uniform array
//                   --> Init_ByFile() does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//                5. The function pointer "Flu_ResetByUser_Func_Ptr" points to "Flu_ResetByUser_Func()" by default
//                   but may be overwritten by various test problem initializers
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Target fluid sandglass
//                TTime : Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API( const int lv, const int FluSg, const double TTime )
{

// check
   if ( Flu_ResetByUser_Func_Ptr == NULL )
   {
      OPT__RESET_FLUID = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : OPT__RESET_FLUID has been disabled since Flu_ResetByUser_Func_Ptr == NULL !!\n" );

      return;
   }


   const double dh       = amr->dh[lv];
#  if ( MODEL == HYDRO )
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
         Reset = Flu_ResetByUser_Func_Ptr( fluid, x, y, z, TTime, lv, NULL );

//       operations necessary only when this cell has been reset
         if ( Reset )
         {
#           if ( MODEL == HYDRO )
#           ifdef MHD
            const real EngyB = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#           else
            const real EngyB = NULL_REAL;
#           endif

//          check minimum density and pressure
            fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
            fluid[ENGY] = Hydro_CheckMinPresInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                    Gamma_m1, _Gamma_m1, MIN_PRES, EngyB );

//          calculate the dual-energy variable (entropy or internal energy)
#           if   ( DUAL_ENERGY == DE_ENPY )
            fluid[ENPY] = Hydro_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1, EngyB );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
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
         } // if ( Reset )

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Flu_ResetByUser_API
