#include "GAMER.h"

// declare as static so that other functions cannot invoke them directly and must use the function pointers
static int Flu_ResetByUser_Func_Template( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                          const double dt, const int lv, double AuxArray[] );
static void Flu_ResetByUser_API_Default( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );

// this function pointer **must** be set by a test problem initializer
int (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                 const double dt, const int lv, double AuxArray[] ) = NULL;
// this function pointer **may** be overwritten by a test problem initializer
void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt ) = Flu_ResetByUser_API_Default;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func_Template
// Description :  Function template to reset the fluid field
//
// Note        :  1. Invoked by Flu_ResetByUser_API_Default() and Model_Init_ByFunction_AssignData() using the
//                   function pointer "Flu_ResetByUser_Func_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked when constructing the initial condition in
//                   Model_Init_ByFunction_AssignData() and after each update in Flu_ResetByUser_API_Default()
//                3. Input fluid[] stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in Flu_ResetByUser_API_Default() and
//                       Model_Init_ByFunction_AssignData()
//                5. Enabled by the runtime options "OPT__RESET_FLUID" or "OPT__RESET_FLUID_INIT"
//
// Parameter   :  fluid    : Fluid array storing both the input (original) and reset values
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
int Flu_ResetByUser_Func_Template( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                   const double dt, const int lv, double AuxArray[] )
{

// Example : reset fluid variables to extremely small values if the cell is within a specific sphere
   /*
   const real dr[3]   = { x-amr->BoxCenter[0], y-amr->BoxCenter[1], z-amr->BoxCenter[2] };
   const real r       = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

   const real TRad                      = 0.3;
   const real MinDens                   = 1.0e-10;
#  if ( NCOMP_PASSIVE > 0 )
   const real MinPassive[NCOMP_PASSIVE] = { ... };
#  else
   const real *MinPassive = NULL;
#  endif
   const real MinPres                   = 1.0e-10;
   const real MinEint                   = EoS_DensPres2Eint_CPUPtr( MinDens, MinPres, MinPassive,
                                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   const real MinEmag                   = 0.0;  // assuming MHD is not adopted

   if ( r <= TRad )
   {
//    set active scalars
      fluid[DENS] = MinDens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], MinEint, MinEmag );
#     ifdef MHD
      fluid[ENGY] += Emag;
#     endif

//    set passive scalars
#     if ( NCOMP_PASSIVE > 0 )
//    fluid[XXXX] = ...;
#     endif

      return true;
   }

   else
      return false;
   */

   return false;

} // FUNCTION : Flu_ResetByUser_Func_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API_Default
// Description :  Default API for resetting the fluid array
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked by EvolveLevel() using the function pointer "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer may be reset by a test problem initializer, in which case
//                       this funtion will become useless
//                3. Currently NOT applied to the input uniform array
//                   --> Init_ByFile() does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//
// Parameter   :  lv      : Target refinement level
//                FluSg   : Target fluid sandglass
//                MagSg   : Target B field sandglass
//                TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_Default( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt )
{

// check
#  ifdef MHD
   if ( Flu_ResetByUser_Func_Ptr == NULL  &&  MHD_ResetByUser_BField_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL and MHD_ResetByUser_BField_Ptr == NULL for OPT__RESET_FLUID !!\n" );
#  else
   if ( Flu_ResetByUser_Func_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL for OPT__RESET_FLUID !!\n" );
#  endif


   const bool   ResetFlu  = ( Flu_ResetByUser_Func_Ptr   != NULL );
#  ifdef MHD
   const bool   ResetMag  = ( MHD_ResetByUser_BField_Ptr != NULL );
   const bool   UseVecPot = ( MHD_ResetByUser_VecPot_Ptr != NULL );
#  endif
   const double dh   = amr->dh[lv];
   const double dh_2 = 0.5*dh;
#  ifdef OPENMP
   const int    NT   = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int    NT   = 1;
#  endif


// 1. reset the magnetic field
#  ifdef MHD
   if ( ResetMag ) {

// 1-1. allocate memory
   real (*Ax)[ CUBE(PS1+1) ] = NULL;
   real (*Ay)[ CUBE(PS1+1) ] = NULL;
   real (*Az)[ CUBE(PS1+1) ] = NULL;
   real (*Emag_old)[PS1][PS1][PS1] = new real [NT][PS1][PS1][PS1];

   if ( UseVecPot )
   {
      Ax = new real [NT][ CUBE(PS1+1) ];
      Ay = new real [NT][ CUBE(PS1+1) ];
      Az = new real [NT][ CUBE(PS1+1) ];
   }


#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    1-2. store the magnetic energy density before resetting
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
         Emag_old[TID][k][j][i] = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );


//    1-3. compute vector potential
//         --> compute the entire patch at once to avoid redundant calculations
      real *AxTID = ( UseVecPot ) ? Ax[TID] : NULL;
      real *AyTID = ( UseVecPot ) ? Ay[TID] : NULL;
      real *AzTID = ( UseVecPot ) ? Az[TID] : NULL;

      if ( UseVecPot )
      {
         int idx = 0;

         for (int k=0; k<PS1+1; k++) {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + k*dh;
         for (int j=0; j<PS1+1; j++) {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + j*dh;
         for (int i=0; i<PS1+1; i++) {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + i*dh;

            AxTID[idx] = (real)0.0;
            AyTID[idx] = (real)0.0;
            AzTID[idx] = (real)0.0;

            if ( i != PS1 )   AxTID[idx] = MHD_ResetByUser_VecPot_Ptr( x+dh_2, y,      z,      TimeNew, dt, lv, 'x', NULL );
            if ( j != PS1 )   AyTID[idx] = MHD_ResetByUser_VecPot_Ptr( x,      y+dh_2, z,      TimeNew, dt, lv, 'y', NULL );
            if ( k != PS1 )   AzTID[idx] = MHD_ResetByUser_VecPot_Ptr( x,      y,      z+dh_2, TimeNew, dt, lv, 'z', NULL );

            idx ++;
         }}} // i,j,k
      } // if ( UseVecPot )


//    1-4. reset B field
//         --> set one component at a time since different components are defined at different cell faces
      for (int v=0; v<NCOMP_MAG; v++)
      {
         int    ijk_end[3], idx=0;
         double dxyz0[3];

         for (int d=0; d<3; d++)
         {
            ijk_end[d] = ( d == v ) ? PS1+1 : PS1;
            dxyz0  [d] = ( d == v ) ? 0.0   : dh_2;
         }

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + dxyz0[0];
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + dxyz0[1];
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + dxyz0[2];

         for (int k=0; k<ijk_end[2]; k++)    {  const double z = z0 + k*dh;
         for (int j=0; j<ijk_end[1]; j++)    {  const double y = y0 + j*dh;
         for (int i=0; i<ijk_end[0]; i++)    {  const double x = x0 + i*dh;

            real B_in, B_out;

            B_in  = amr->patch[MagSg][lv][PID]->magnetic[v][idx];
            B_out = MHD_ResetByUser_BField_Ptr( x, y, z, TimeNew, dt, lv, 'x'+v, NULL, B_in,
                                                UseVecPot, AxTID, AyTID, AzTID, i, j, k );

            amr->patch[MagSg][lv][PID]->magnetic[v][ idx ++ ] = B_out;
         }}} // i,j,k
      } // for (int v=0; v<NCOMP_MAG; v++)


//    1-5. update the total energy density
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const real Emag_new = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
         amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i] += Emag_new - Emag_old[TID][k][j][i];
      }
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// 1-6. free memory
   if ( UseVecPot )
   {
      delete [] Ax;
      delete [] Ay;
      delete [] Az;
   }
   delete [] Emag_old;
   } // if ( ResetMag )
#  endif // #ifdef MHD



// 2. reset the fluid fields
   if ( ResetFlu ) {
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + dh_2;
      const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + dh_2;
      const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + dh_2;

      real fluid[NCOMP_TOTAL];

      for (int k=0; k<PS1; k++)  {  const double z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  const double y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  const double x = x0 + i*dh;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#        else
         const real Emag = (real)0.0;
#        endif

//       2-1. reset this cell
         const bool Reset = Flu_ResetByUser_Func_Ptr( fluid, Emag, x, y, z, TimeNew, dt, lv, NULL );

//       2-2. operations necessary only when this cell has been reset
         if ( Reset )
         {
#           if ( MODEL == HYDRO )
//          apply density and internal energy floors
            fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
#           ifndef SRHD
            fluid[ENGY] = Hydro_CheckMinEintInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                    MIN_EINT, PassiveFloorMask, Emag );
#           endif

//          calculate the dual-energy variable (entropy or internal energy)
#           ifdef DUAL_ENERGY
            fluid[DUAL] = Hydro_Con2Dual( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Emag,
                                          EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                          PassiveFloorMask );
#           endif

//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
               if ( PassiveFloorMask & BIDX(v) )  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif
#           endif // if ( MODEL == HYDRO )

//          store the reset values
            for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
         } // if ( Reset )

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // if ( ResetFlu )

} // FUNCTION : Flu_ResetByUser_API_Default
