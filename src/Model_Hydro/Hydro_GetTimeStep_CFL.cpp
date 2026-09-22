#include "GAMER.h"

#if ( MODEL == HYDRO )

static real GetMaxCFL( const int lv, const bool isCRDiffusion );




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetTimeStep_CFL
// Description :  Estimate the evolution time-step from the CFL condition of the hydro/MHD solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Time-step is estimated by the stability criterion from the von Neumann stability analysis
//                   --> CFL condition
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Hydro_GetTimeStep_CFL( const int lv )
{

// get the maximum CFL
   const bool   CRDiffusion_No     = false;
   const real   MaxCFL             = GetMaxCFL( lv, CRDiffusion_No  );
#  ifdef CR_DIFFUSION
   const bool   CRDiffusion_Yes    = true;
   const real   MaxCFL_CRDiffusion = GetMaxCFL( lv, CRDiffusion_Yes );
#  endif

// get the time-step
   const double dh     = amr->dh[lv];
   const double Safety = (Step==0) ? DT__FLUID_INIT : DT__FLUID;
         double dt     = Safety*dh/MaxCFL;

#  ifdef CR_DIFFUSION
   const double dt_CR  = MicroPhy.CR_safety*0.5*dh*dh/MaxCFL_CRDiffusion;
   dt = fmin( dt, dt_CR );
#  endif

   return dt;

} // FUNCTION : Hydro_GetTimeStep_CFL



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxCFL
// Description :  Get the maximum CFL at the target level among all MPI ranks
//
// Note        :  1. Invoked by Hydro_GetTimeStep_CFL()
//
// Parameter   :  lv            : Target refinement level
//                isCRDiffusion : Whether it is the CFL for CR diffusion
//
// Return      :  MaxCFL
//-------------------------------------------------------------------------------------------------------
real GetMaxCFL( const int lv, const bool isCRDiffusion )
{

   const int FluSg = amr->FluSg[lv];
#  ifdef MHD
   const int MagSg = amr->MagSg[lv];
#  endif


// maximum CFL in this rank
   real MaxCFL  = (real)0.0;
   bool AnyCell = false;


// get the maximum CFL in this rank
#  pragma omp parallel for reduction( max:MaxCFL ) reduction( ||:AnyCell ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      bool Skip = false;

//    if OPT__DT_FLUID_LEAF, skip all non-leaf patches
//    because their data are later overwritten by the refined patches
//    note that this leads to the fluid timestep being "inf" when a level is completely refined
      if ( OPT__DT_FLUID_LEAF  &&  amr->patch[0][lv][PID]->son != -1 )
      {
         Skip = true;
      }

      if ( Skip )    continue;
      else           AnyCell = true;


//    calculate the CFL
      for (int k=0; k<PATCH_SIZE; k++) {
      for (int j=0; j<PATCH_SIZE; j++) {
      for (int i=0; i<PATCH_SIZE; i++) {

         real fluid[FLU_NIN_T];

         for (int v=0; v<FLU_NIN_T; v++)  fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

         real B[3] = { (real)0.0, (real)0.0, (real)0.0 };
#        ifdef MHD
         MHD_GetCellCenteredBField( B,
                                    amr->patch[MagSg][lv][PID]->magnetic[MAGX],
                                    amr->patch[MagSg][lv][PID]->magnetic[MAGY],
                                    amr->patch[MagSg][lv][PID]->magnetic[MAGZ],
                                    PS1, PS1, PS1, i, j, k );
#        endif // #ifdef MHD

         const real CFL =
#                         ifdef CR_DIFFUSION
                          ( isCRDiffusion ) ?
                          Hydro_GetCFL_CRDiffusion( MicroPhy ) :
#                         endif
                          Hydro_GetCFL( fluid, B, (real)MIN_PRES, PassiveFloorMask,
                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                        EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

         MaxCFL = FMAX( MaxCFL, CFL );

      }}} // k,j,i
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// get the maximum CFL in all ranks
   real MaxCFL_AllRank;
   bool AnyCell_AllRank;

   MPI_Allreduce(  &MaxCFL,  &MaxCFL_AllRank, 1, MPI_GAMER_REAL, MPI_MAX,    MPI_COMM_WORLD );
   MPI_Reduce   ( &AnyCell, &AnyCell_AllRank, 1,     MPI_C_BOOL, MPI_LOR, 0, MPI_COMM_WORLD );


// check
   if ( MaxCFL_AllRank == 0.0  &&  AnyCell_AllRank  &&  MPI_Rank == 0 )
      Aux_Error( ERROR_INFO, "MaxCFL == 0.0 at lv %d !!\n", lv );


   return MaxCFL_AllRank;

} // FUNCTION : GetMaxCFL



#endif // #if ( MODEL == HYDRO )
