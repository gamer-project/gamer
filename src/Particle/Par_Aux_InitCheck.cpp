#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Aux_InitCheck
// Description :  Check the initial condition of particles 
//
// Note        :  1. Invoked by "Init_GAMER"
//                2. Check if all particles lie within the simulation box
//                3. Remove particles outside the active region for non-periodic B.C.
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Par_Aux_InitCheck()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };


// 1. all active particles should lie within the simulation domain
// (periodicity should be taken care of in the initial condition, not here)
   for (long ParID=0; ParID<amr->Par->NPar; ParID++)
   {
//    skip inactive particles
      if ( Mass[ParID] < 0.0 )   continue;

      for (int d=0; d<3; d++)
      {
         if ( Pos[d][ParID] < (real)0.0  ||  Pos[d][ParID] >= amr->BoxSize[d] )
            Aux_Error( ERROR_INFO, "Pos[%d][%ld] = %14.7e lies outside the simulation domain (0.0 ... %13.7e) !!\n",
                       d, ParID, Pos[d][ParID], amr->BoxSize[d] );
      }
   }


// 2. remove particles outside the active region for non-periodic B.C.
   if ( OPT__BC_POT != BC_POT_PERIODIC )
   {
      for (long ParID=0; ParID<amr->Par->NPar; ParID++)
      {
//       skip inactive particles
         if ( Mass[ParID] < 0.0 )   continue;

         if (  !Par_WithinActiveRegion( Pos[0][ParID], Pos[1][ParID], Pos[2][ParID] )  )
         {
            amr->Par->RemoveOneParticle( ParID, PAR_INACTIVE_OUTSIDE, NULL_INT, NULL, NULL_REAL );

            if ( OPT__VERBOSE )
               Aux_Message( stderr, "\nWARNING : removing particle %10d (Pos = [%14.7e, %14.7e, %14.7e], Time = %13.7e)\n",
                            ParID, Pos[0][ParID], Pos[1][ParID], Pos[2][ParID], Time[0] );
         }
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Aux_InitCheck



#endif // #ifdef PARTICLE
