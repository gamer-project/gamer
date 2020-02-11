#include "GAMER.h"

#ifdef GRAVITY

extern real (*Poi_AddExtraMassForGravity_Ptr)( const double x, const double y, const double z, const double Time,
                                               const int lv, double AuxArray[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_Rho
// Description :  Prepare h_Rho_Array_P[] for the Poisson solver
//
// Note        :  1. Invoke Prepare_PatchData()
//                2. Minimum density threshold (MIN_DENS) is applied
//
// Parameter   :  lv            : Target refinement level
//                PrepTime      : Target physical time to prepare the coarse-grid data
//                h_Rho_Array_P : Host array to store the prepared data
//                NPG           : Number of patch groups to be prepared at a time
//                PID0_List     : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Poi_Prepare_Rho( const int lv, const double PrepTime, real h_Rho_Array_P[][RHO_NXT][RHO_NXT][RHO_NXT],
                      const int NPG, const int *PID0_List )
{

   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinPres_No        = -1.0;

   Prepare_PatchData( lv, PrepTime, &h_Rho_Array_P[0][0][0][0], NULL, RHO_GHOST_SIZE, NPG, PID0_List, _TOTAL_DENS, _NONE,
                      OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                      MIN_DENS, MinPres_No, DE_Consistency_No );


// add extra mass source for gravity if required
   if ( OPT__GRAVITY_EXTRA_MASS )
   {
      const double dh          = amr->dh[lv];
      const double L[3]        = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
      const bool   Periodic[3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                   OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                   OPT__BC_FLU[4] == BC_FLU_PERIODIC };

      for (int TID=0; TID<NPG; TID++)
      {
         const int PID0 = PID0_List[TID];

         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int    PID = PID0 + LocalID;
            const int    N   = 8*TID + LocalID;
            const double x0  = amr->patch[0][lv][PID]->EdgeL[0] + (0.5-RHO_GHOST_SIZE)*dh;
            const double y0  = amr->patch[0][lv][PID]->EdgeL[1] + (0.5-RHO_GHOST_SIZE)*dh;
            const double z0  = amr->patch[0][lv][PID]->EdgeL[2] + (0.5-RHO_GHOST_SIZE)*dh;

            double x, y, z;

            for (int k=0; k<RHO_NXT; k++)  {  z = z0 + k*dh;  if ( Periodic[2] )  z = fmod( z+L[2], L[2] );
            for (int j=0; j<RHO_NXT; j++)  {  y = y0 + j*dh;  if ( Periodic[1] )  y = fmod( y+L[1], L[1] );
            for (int i=0; i<RHO_NXT; i++)  {  x = x0 + i*dh;  if ( Periodic[0] )  x = fmod( x+L[0], L[0] );

               h_Rho_Array_P[N][k][j][i] += Poi_AddExtraMassForGravity_Ptr( x, y, z, Time[lv], lv, NULL );

            }}}
         } // for (int LocalID=0; LocalID<8; LocalID++)
      } // for (int TID=0; TID<NPG; TID++)
   } // if ( OPT__GRAVITY_EXTRA_MASS )


// subtract the background density to be consistent with the periodic (and/or comoving) base-level FFT solver
#  ifdef GAMER_DEBUG
   if ( OPT__BC_POT == BC_POT_PERIODIC  &&  AveDensity_Init <= 0.0 )
      Aux_Error( ERROR_INFO, "AveDensity_Init (%14.7e) hasn't been set properly !!\n", AveDensity_Init );
#  endif

#  ifdef COMOVING
   const bool Comoving = true;
#  else
   const bool Comoving = false;
#  endif

   real RhoSubtract = NULL_REAL;
   if      ( OPT__BC_POT == BC_POT_PERIODIC )   RhoSubtract = AveDensity_Init;
   else if ( Comoving )                         RhoSubtract = (real)1.0;

   if ( OPT__BC_POT == BC_POT_PERIODIC  ||  Comoving )
   {
#     pragma omp parallel for schedule( static )
      for (int TID=0; TID<8*NPG; TID++)
      for (int k=0; k<RHO_NXT; k++)
      for (int j=0; j<RHO_NXT; j++)
      for (int i=0; i<RHO_NXT; i++)
         h_Rho_Array_P[TID][k][j][i] -= RhoSubtract;
   }

} // FUNCTION : Poi_Prepare_Rho



#endif // #ifdef GRAVITY
