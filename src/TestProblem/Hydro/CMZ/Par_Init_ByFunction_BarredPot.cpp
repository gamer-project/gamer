#include "GAMER.h"

#ifdef PARTICLE


//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_BarredPot
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank   : Number of particles to be set by this MPI rank
//                NPar_AllRank    : Total Number of particles in all MPI ranks
//                ParMass         : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z     : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z     : Particle velocity array with the size of NPar_ThisRank
//                ParTime         : Particle time     array with the size of NPar_ThisRank
//                ParType         : Particle type     array with the size of NPar_ThisRank
//                AllAttributeFlt : Pointer array for all particle floating-point attributes
//                                  --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                      to access the data
//                AllAttributeInt : Pointer array for all particle integer attributes
//                                  --> Dimension = [PAR_NATT_INT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttributeFlt, AllAttributeInt
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_BarredPot( const long NPar_ThisRank, const long NPar_AllRank,
                                  real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                  real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                  long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                  long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing problem-specific particle attributes ... " );

   const real_par Useless = -1.0;

#  ifdef STAR_FORMATION
   for (int p=0; p<NPar_ThisRank; p++)    AllAttributeFlt[Idx_ParCreTime  ][p] = Useless;
#  endif

   for (int p=0; p<NPar_ThisRank; p++)    AllAttributeFlt[Idx_ParMetalFrac][p] = Useless;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


}


#endif // #ifdef PARTICLE
