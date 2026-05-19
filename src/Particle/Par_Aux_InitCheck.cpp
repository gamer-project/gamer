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
//                5. Check particle types
//                6. Check particle UID (which might not be assigned yet when calling this function)
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Par_Aux_InitCheck()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const real_par *Mass   =   amr->Par->Mass;
   const real_par *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const long_par *Type   =   amr->Par->Type;
   const long_par *PUID   =   amr->Par->PUID;
   const long_par *Flag   =   amr->Par->Flag;


// 1. all active particles should lie within the simulation domain
//    --> periodicity should be taken care of in the initial condition, not here
//    --> also check particle types and UID here
   for (long ParIdx=0; ParIdx<amr->Par->NPar_AcPlusInac; ParIdx++)
   {
//    there should be no inactive particles initially
      if ( Mass[ParIdx] < 0.0 )   Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParIdx, Mass[ParIdx] );

//    check particle types
      if ( Type[ParIdx] < (real_par)0  ||  Type[ParIdx] >= (real_par)PAR_NTYPE )
         Aux_Error( ERROR_INFO, "Type[%ld] = %d (accepted range: 0<=index<%d) !!\n", ParIdx, (int)Type[ParIdx], PAR_NTYPE );

//    check particle UID
//    a particle here can be initialized via PAR_INIT_BY_FUNCTION or PAR_INIT_BY_FILE
//    --> if the particle UID is loaded by Par_Init_ByFile(),
//        it must lie within [1, NextPUID-1]
//    --> if it is not loaded from a file,
//        the initialized particle should not be assigned a new UID beforehand
//        and must have particle UID == PUID_TBA,
//        so that Par_SetParUID() can be invoked to assign it properly
      if ( PUID[ParIdx] != PUID_TBA  &&  ( PUID[ParIdx] <= 0  ||  PUID[ParIdx] >= amr->Par->NextPUID ) )
         Aux_Error( ERROR_INFO, "PUID[%ld] = %ld (accepted range: 0<index<%ld) !!\n", ParIdx, (long)PUID[ParIdx], amr->Par->NextPUID );

//    check particle flags
      if ( Flag[ParID] == PFLAG_TBA )
         Aux_Error( ERROR_INFO, "Flag[%ld] = %d (PFLAG_TBA) !!\n", ParID, (int)Flag[ParID] );

//    only support tracer particles when disabling GRAVITY
#     ifndef GRAVITY
      if ( Type[ParIdx] != PTYPE_TRACER )
         Aux_Error( ERROR_INFO, "Type[%ld] = %d != PTYPE_TRACER (%d) when disabling GRAVITY !!\n", ParIdx, (int)Type[ParIdx], (int)PTYPE_TRACER );
#     endif

//    must enable TRACER for tracer particles
#     ifndef TRACER
      if ( Type[ParIdx] == PTYPE_TRACER )
         Aux_Error( ERROR_INFO, "Type[%ld] = %d (PTYPE_TRACER) when disabling TRACER !!\n", ParIdx, (int)Type[ParIdx] );
#     endif

//    tracer particles must be massless
#     ifdef TRACER
      if ( Type[ParIdx] == PTYPE_TRACER  &&  Mass[ParIdx] != (real_par)0.0 )
         Aux_Error( ERROR_INFO, "Tracer[%ld] has non-zero mass (%13.7e) !!\n", ParIdx, Mass[ParIdx] );
#     endif

      for (int d=0; d<3; d++)
      {
         if ( Pos[d][ParIdx] < (real_par)0.0  ||  Pos[d][ParIdx] >= amr->BoxSize[d] )
         {
//          periodicity should be taken care of in advance
            if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "Pos[%d][%ld] = %14.7e lies outside the simulation domain (0.0 ... %13.7e) !!\n",
                       d, ParIdx, Pos[d][ParIdx], amr->BoxSize[d] );

//          for non-periodic BC., particles lying outside the box will be removed in the next check
            else
            Aux_Message( stderr, "WARNING : Pos[%d][%ld] = %14.7e lies outside the simulation domain (0.0 ... %13.7e) !!\n",
                        d, ParIdx, Pos[d][ParIdx], amr->BoxSize[d] );
         }
      }
   } // for (long ParIdx=0; ParIdx<amr->Par->NPar_AcPlusInac; ParIdx++)


// 2. remove particles outside the active region for non-periodic B.C.
   for (long ParIdx=0; ParIdx<amr->Par->NPar_AcPlusInac; ParIdx++)
   {
      if (  !Par_WithinActiveRegion( Pos[0][ParIdx], Pos[1][ParIdx], Pos[2][ParIdx] )  )
      {
         amr->Par->RemoveOneParticle( ParIdx, PAR_INACTIVE_OUTSIDE );

         Aux_Message( stderr, "WARNING : removing particle %10d (Pos = [%14.7e, %14.7e, %14.7e], Time = %13.7e)\n",
                      ParIdx, Pos[0][ParIdx], Pos[1][ParIdx], Pos[2][ParIdx], Time[0] );
      }
   }


// 3. get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Aux_InitCheck



#endif // #ifdef PARTICLE
