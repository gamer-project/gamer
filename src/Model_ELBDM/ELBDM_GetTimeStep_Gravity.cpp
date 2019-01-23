#include "GAMER.h"

#if ( MODEL == ELBDM  &&  defined GRAVITY )

#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];

static real GetMaxPot( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Gravity
// Description :  Estimate the evolution time-step from the ELBDM potential energy solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Time-step is set to restrict the real-space rotation angle to be "DT__GRAVITY*2*pi"
//                   --> dt = DT__GRAVITY*2*pi/Eta/Max(Pot)
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Gravity( const int lv )
{

   real   MaxPot;
   double dt;

// get the maximum potential
   MaxPot = GetMaxPot( lv );

// get the time-step
   dt = DT__GRAVITY*2.0*M_PI/(ELBDM_ETA*MaxPot);

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Gravity



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxPot
// Description :  Get the maximum potential at the target level among all MPI ranks
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Gravity()
//                2. Include gravitational, self-interaction, and external potentials
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MaxPot
//-------------------------------------------------------------------------------------------------------
real GetMaxPot( const int lv )
{

   real   PotG, PotS, PotE;      // PotG/S/E: gravitational/self-interaction/external potential
   real   Pot, MaxPot=0.0;       // Pot = PotG + PotS + PotE
   double x0, y0, z0, x, y, z;
   int    SibPID;
   bool   Skip;


// get the maximum potential in this rank
#  pragma omp parallel for private( PotG, PotS, PotE, Pot, x0, y0, z0, x, y, z, SibPID, Skip ) reduction( max:MaxPot ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip all non-leaf patches not adjacent to any coarse-fine boundaries (since their data are useless)
      Skip = true;

      if ( amr->patch[0][lv][PID]->son == -1 )  Skip = false;
      else
      {
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][lv][PID]->sibling[s];

//          proper-nesting check
#           ifdef GAMER_DEBUG
            if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 (Lv %d, PID %d) !!\n", lv, PID );
#           endif

//          non-periodic BC.
            if ( SibPID < -1 )   continue;

//          check whether this patch is adjacent to a coarse-fine boundary
            if ( amr->patch[0][lv][SibPID]->son == -1 )
            {
               Skip = false;
               break;
            }
         } // for (int s=0; s<26; s++)
      } // if ( amr->patch[0][lv][PID]->son == -1 ) ... else ...

      if ( Skip )    continue;


//    calculate the potential
      const double dh = amr->dh[lv];

      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PATCH_SIZE; k++) {  z = z0 + (double)k*dh;
      for (int j=0; j<PATCH_SIZE; j++) {  y = y0 + (double)j*dh;
      for (int i=0; i<PATCH_SIZE; i++) {  x = x0 + (double)i*dh;

         PotG   = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];
#        ifdef QUARTIC_SELF_INTERACTION
         PotS   = ELBDM_LAMBDA*amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
#        else
         PotS   = (real)0.0;
#        endif
         PotE   = ( OPT__EXTERNAL_POT ) ? ExternalPot( x, y, z, Time[lv], ExtPot_AuxArray ) : (real)0.0;

         Pot    = FABS( PotG + PotS + PotE );  // remember to take the absolute value
         MaxPot = MAX( MaxPot, Pot );
      }}} // k,j,i
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// get the maximum potential in all ranks
   real MaxPot_AllRank;
#  ifdef FLOAT8
   MPI_Allreduce( &MaxPot, &MaxPot_AllRank, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#  else
   MPI_Allreduce( &MaxPot, &MaxPot_AllRank, 1, MPI_FLOAT,  MPI_MAX, MPI_COMM_WORLD );
#  endif


// check
   if ( MaxPot_AllRank == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MaxPot == 0.0 at lv %d !!\n", lv );


   return MaxPot_AllRank;

} // FUNCTION : GetMaxPot



#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )
