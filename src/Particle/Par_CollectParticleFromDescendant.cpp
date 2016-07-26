#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE

# ifndef LOAD_BALANCE
static void CollectParticle( const int FaLv, const int FaPID, int &NPar_SoFar, long *ParList );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_CollectParticleFromDescendant
// Description :  Collect particles from all descendants (sons, grandsons, ...) to father
//                --> Will call Par_LB_CollectParticleFromDescendant for LOAD_BALANCE
//
// Note        :  1. ParList in all descendants will NOT be changeed after calling this function
//                2. This function assumes that father and all descendants are always in the same rank
//                   --> Does NOT work with LOAD_BALANCE
//                3. The array "ParList_Desc" will be allocated for the target patch
//                   --> Must be deleted afterward by calling Par_CollectParticleFromDescendant_FreeMemory
//                       (e.g., in the function "Gra_AdvanceDt")
//                4. Do not take into account it's own particles (particles at FaLv)
//                   --> Do nothing if this patch is a leaf patch (ParList_Desc will NOT be allocated)
//                5. This function assumes that all particles live in the leaf patch
//                6. When using OpenMP, one must ensure that different threads do NOT invoke this function
//                   for the same patch at the same time !!!
//                   --> Because this function will modify "NPar_Desc & ParList_Desc" for the target patch
//                7. Invoked by Gra_AdvanceDt (and Main when DEBUG is on)
//
// Parameter   :  FaLv        : Father's refinement leve
//                PredictPos  : true --> predict particle position to TargetTime (for LOAD_BALANCE only)
//                TargetTime  : Target time for predicting the particle position (for LOAD_BALANCE only)
//
// Return      :  NPar_Desc and ParList_Desc for all non-leaf patches at FaLv
//-------------------------------------------------------------------------------------------------------
void Par_CollectParticleFromDescendant( const int FaLv, const bool PredictPos, const double TargetTime )
{

// nothing to do for the max level
   if ( FaLv >= MAX_LEVEL )   return;


// call the parallel version instead
#  ifdef LOAD_BALANCE
   Par_LB_CollectParticleFromDescendant( FaLv, PredictPos, TargetTime );

   return;


#  else


//###NOTE: OpenMP may not improve performance here
#  pragma omp parallel for schedule( runtime )
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
//    check
#     ifdef DEBUG_PARTICLE
      if ( amr->patch[0][FaLv][FaPID]->NPar_Desc != -1  ||  amr->patch[0][FaLv][FaPID]->ParList_Desc != NULL )
         Aux_Error( ERROR_INFO, "Desc variables have been initialized already (FaLv %d, FaPID %d, NPar_Dens %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar_Desc );
#     endif

//    nothing to do if father has no son
      if ( amr->patch[0][FaLv][FaPID]->son == -1 )    continue;


//    1. get the total number of particles in all descendants
      amr->patch[0][FaLv][FaPID]->NPar_Desc = Par_CountParticleInDescendant( FaLv, FaPID );


//    nothing to do if there are no particles in descendants
      if ( amr->patch[0][FaLv][FaPID]->NPar_Desc == 0 )  continue;


//    2. allocate the array ParList_Desc
      amr->patch[0][FaLv][FaPID]->ParList_Desc = new long [ amr->patch[0][FaLv][FaPID]->NPar_Desc ];


//    3. gather the particle lists from all descendants
      int NPar_SoFar = 0;
      CollectParticle( FaLv, FaPID, NPar_SoFar, amr->patch[0][FaLv][FaPID]->ParList_Desc );


//    check the particle count
#     ifdef DEBUG_PARTICLE
      if ( NPar_SoFar != amr->patch[0][FaLv][FaPID]->NPar_Desc )
         Aux_Error( ERROR_INFO, "NPar_SoFar (%d) != NPar_Desc (%d) !!\n", 
                    NPar_SoFar, amr->patch[0][FaLv][FaPID]->NPar_Desc );
#     endif
   } // for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)


#  endif // #ifdef LOAD_BALANCE ... else ...


} // FUNCTION : Par_CollectParticleFromDescendant



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_CollectParticleFromDescendant_FreeMemory
// Description :  Release the memory allocated by Par_CollectParticleFromDescendant
//
// Note        :  1. Invoded by Gra_AdvanceDt (and Main when DEBUG is on)
//
// Parameter   :  FaLv  : Father's refinement level
//-------------------------------------------------------------------------------------------------------
void Par_CollectParticleFromDescendant_FreeMemory( const int FaLv )
{

#  ifdef LOAD_BALANCE

#  else

   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->ParList_Desc != NULL )
      {
         delete [] amr->patch[0][FaLv][FaPID]->ParList_Desc;
         amr->patch[0][FaLv][FaPID]->ParList_Desc = NULL;
      }

      amr->patch[0][FaLv][FaPID]->NPar_Desc = -1;
   }

#  endif // #ifdef LOAD_BALANCE ... else ...

} // FUNCTION : Par_CollectParticleFromDescendant_FreeMemory



# ifndef LOAD_BALANCE
//-------------------------------------------------------------------------------------------------------
// Function    :  CollectParticle 
// Description :  Collect particles from all descendants (sons, grandsons, ...) of the target patch
//
// Note        :  This function will search over all descendants recursively
//
// Parameter   :  FaLv        : Father patch level
//                FaPID       : Father patch ID
//                NPar_SoFar  : Number of particles counted so for 
//                              --> used as the starting array index in ParList
//                ParList     : Array to store the particle IDs
//
// Return      :  NPar_SoFar, ParList
//-------------------------------------------------------------------------------------------------------
void CollectParticle( const int FaLv, const int FaPID, int &NPar_SoFar, long *ParList )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( FaLv < 0  ||  FaLv > TOP_LEVEL )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );
#  endif


   const int SonPID0 = amr->patch[0][FaLv][FaPID]->son;
   const int SonLv   = FaLv + 1;

   if ( SonPID0 == -1 )    return;  // nothing to do if the target patch has no son
   else
   {
      for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      {
//       search over all grandsons recursively
         if ( amr->patch[0][SonLv][SonPID]->son != -1 )  CollectParticle( SonLv, SonPID, NPar_SoFar, ParList ); 

         else
            for (int p=0; p<amr->patch[0][SonLv][SonPID]->NPar; p++)    
               ParList[ NPar_SoFar ++ ] = amr->patch[0][SonLv][SonPID]->ParList[p];
      }
   }

} // FUNCTION : CollectParticle
#endif // #ifndef LOAD_BALANCE



#endif // #ifdef PARTICLE
