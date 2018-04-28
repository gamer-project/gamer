#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Aux_InitCheck
// Description :  Check the initial condition of particles
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Check if all particles lie within the simulation box
//                3. Remove particles outside the active region for non-periodic B.C.
//                4. There should be no inactive particles before calling this function
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Par_Aux_InitCheck()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };


// 1. all active particles should lie within the simulation domain
//    --> periodicity should be taken care of in the initial condition, not here
   for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
//    there should be no inactive particles initially
      if ( Mass[ParID] < 0.0 )   Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, Mass[ParID] );

      for (int d=0; d<3; d++)
      {
         if ( Pos[d][ParID] < (real)0.0  ||  Pos[d][ParID] >= amr->BoxSize[d] )
         {
//          periodicity should be taken care of in advance
            if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "Pos[%d][%ld] = %14.7e lies outside the simulation domain (0.0 ... %13.7e) !!\n",
                       d, ParID, Pos[d][ParID], amr->BoxSize[d] );

//          for non-periodic BC., particles lying outside the box will be removed in the next check
            else
            Aux_Message( stderr, "WARNING : Pos[%d][%ld] = %14.7e lies outside the simulation domain (0.0 ... %13.7e) !!\n",
                        d, ParID, Pos[d][ParID], amr->BoxSize[d] );
         }
      }
   }


// 2. remove particles outside the active region for non-periodic B.C.
   for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
      if (  !Par_WithinActiveRegion( Pos[0][ParID], Pos[1][ParID], Pos[2][ParID] )  )
      {
         amr->Par->RemoveOneParticle( ParID, PAR_INACTIVE_OUTSIDE );

         Aux_Message( stderr, "WARNING : removing particle %10d (Pos = [%14.7e, %14.7e, %14.7e], Time = %13.7e)\n",
                      ParID, Pos[0][ParID], Pos[1][ParID], Pos[2][ParID], Time[0] );
      }
   }


// 3. get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Aux_InitCheck



#endif // #ifdef PARTICLE
