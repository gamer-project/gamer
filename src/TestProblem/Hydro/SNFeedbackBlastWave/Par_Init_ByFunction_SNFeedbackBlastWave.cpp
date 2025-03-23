#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

extern double SNBlast_ParMass;
extern double SNBlast_ParCenter[3];
extern double SNBlast_ParVelocity[3];
extern double SNBlast_ParMetalMassFrac;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_SNFeedbackBlastWave
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
void Par_Init_ByFunction_SNFeedbackBlastWave( const long NPar_ThisRank, const long NPar_AllRank,
                                              real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                              real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                              long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                              long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// assume there is only one particle in the box for this test problem
   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParMass[p] = (real_par)SNBlast_ParMass;
      ParPosX[p] = (real_par)SNBlast_ParCenter[0];
      ParPosY[p] = (real_par)SNBlast_ParCenter[1];
      ParPosZ[p] = (real_par)SNBlast_ParCenter[2];
      ParVelX[p] = (real_par)SNBlast_ParVelocity[0];
      ParVelY[p] = (real_par)SNBlast_ParVelocity[1];
      ParVelZ[p] = (real_par)SNBlast_ParVelocity[2];
      ParTime[p] = (real_par)Time[0];
      ParType[p] = PTYPE_STAR;

#     ifdef FEEDBACK
      AllAttributeFlt[Idx_ParSNIITime][p]  = (real_par)Time[0] + (real_par)FB_RESOLVED_SNEII_DELAY_TIME;
#     endif

      AllAttributeFlt[Idx_ParMetalFrac][p] = (real_par)SNBlast_ParMetalMassFrac;

   } // for (long p=0; p<NPar_AllRank; p++)


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_SNFeedbackBlastWave



#endif // #ifdef MASSIVE_PARTICLES
