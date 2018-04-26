#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PredictPos
// Description :  Predict particle position to a target time
//
// Note        :  1. Skip particles with time < 0.0
//                   --> These are the particles waiting for velocity correction in the KDK scheme
//                   --> We assume that they have already been synchronized with TargetTime
//                2. Called by Par_MassAssignment(), Par_LB_CollectParticle2OneLevel(), and
//                   Par_LB_CollectParticleFromRealPatch()
//                3. This function does NOT take care of periodicity
//                   --> Particle may lie outside the simulation domain after prediction
//
// Parmaeter   :  NPar        : Number of target particles
//                ParList     : List of target particle IDs
//                ParPosX/Y/Z : x/y/z particle position arrays
//                TargetTime  : Target physical time
//
// Return      :  ParPosX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_PredictPos( const long NPar, const long *ParList, real *ParPosX, real *ParPosY, real *ParPosZ,
                     const double TargetTime )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( NPar < 0 )            Aux_Error( ERROR_INFO, "NPar = %ld < 0 !!\n", NPar );
   else if ( NPar > 0 ) {
   if ( ParList == NULL )     Aux_Error( ERROR_INFO, "ParList == NULL for NPar = %ld !!\n", NPar );
   if ( ParPosX == NULL )     Aux_Error( ERROR_INFO, "ParPosX == NULL for NPar = %ld !!\n", NPar );
   if ( ParPosY == NULL )     Aux_Error( ERROR_INFO, "ParPosY == NULL for NPar = %ld !!\n", NPar );
   if ( ParPosZ == NULL )     Aux_Error( ERROR_INFO, "ParPosZ == NULL for NPar = %ld !!\n", NPar ); }
   if ( TargetTime < 0.0 )    Aux_Error( ERROR_INFO, "TargetTime = %14.7e < 0.0 !!\n", TargetTime );
#  endif


   long ParID;
   real dt;

   for (long p=0; p<NPar; p++)
   {
      ParID = ParList[p];

#     ifdef DEBUG_PARTICLE
      if ( ParID < 0  ||  ParID >= amr->Par->NPar_AcPlusInac )
         Aux_Error( ERROR_INFO, "ParID (%ld) lies outside the accepted range (0 <= ParID < %ld) !!\n",
                    ParID, amr->Par->NPar_AcPlusInac );

      if ( amr->Par->Mass[ParID] < (real)0.0 )
         Aux_Error( ERROR_INFO, "Found inactive particle (ParID %ld, Mass %14.7e) !!\n", ParID, amr->Par->Mass[ParID] );
#     endif

//    skip particles waiting for velocity correction (they should already be synchronized with TargetTime)
      if ( amr->Par->Time[ParID] < (real)0.0 )  continue;

      dt = (real)TargetTime - amr->Par->Time[ParID];

//    note that we don't have to worry about the periodic BC here (in other words, ParPos can lie outside the simulation box)
      ParPosX[p] += amr->Par->VelX[ParID]*dt;
      ParPosY[p] += amr->Par->VelY[ParID]*dt;
      ParPosZ[p] += amr->Par->VelZ[ParID]*dt;
   } // for (long p=0; p<NPar; p++)

} // FUNCTION : Par_PredictPos



#endif // #ifdef PARTICLE
