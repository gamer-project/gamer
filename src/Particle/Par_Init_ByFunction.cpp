#include "GAMER.h"

#ifdef PARTICLE

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Par_Init_ByFunction_Template( const long NPar_ThisRank, const long NPar_AllRank,
                                          real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                          real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                          long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                          long_par *AllAttributeInt[PAR_NATT_INT_TOTAL]);

// this function pointer must be set by a test problem initializer
void (*Par_Init_ByFunction_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                 real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                 real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                 long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                 long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Template
// Description :  Template of user-specified particle initializer
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr",
//                   which must be set by a test problem initializer
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
void Par_Init_ByFunction_Template( const long NPar_ThisRank, const long NPar_AllRank,
                                   real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                   real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                   long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                   long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// synchronize all particles to the physical time on the base level
// and assign particle type
   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParTime[p] = Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }


// initialize the particle creation time by an arbitrary negative value since it is
// only used for star particles created during evolution and is useless during initialization
#  ifdef STAR_FORMATION
   for (int p=0; p<NPar_ThisRank; p++)    AllAttributeFlt[Idx_ParCreTime][p] = -1.0;
#  endif


// set other particle attributes
// ============================================================================================================
   real_par *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real_par *ParVel[3] = { ParVelX, ParVelY, ParVelZ };

// example : randomly initialize
   /*
   const uint     RSeed     = 2;                                             // random seed
   const real_par MassMin   = 1.0e-2;                                        // minimum value of particle mass
   const real_par MassMax   = 1.0;                                           // maximum value of particle mass
   const real_par PosMin[3] = { 0.0, 0.0, 0.0 };                             // minimum value of particle position
   const real_par PosMax[3] = { real_par( amr->BoxSize[0]*(1.0-1.0e-5) ),    // maximum value of particle position
                                real_par( amr->BoxSize[1]*(1.0-1.0e-5) ),
                                real_par( amr->BoxSize[2]*(1.0-1.0e-5) ) };
   const real_par VelMin[3] = { -1.0, -1.0, -1.0 };                          // minimum value of particle velocity
   const real_par VelMax[3] = { +1.0, +1.0, +1.0 };                          // maximum value of particle velocity

   srand( RSeed );

   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParMass[p] = ( (real_par)rand()/RAND_MAX )*( MassMax - MassMin ) + MassMin;

      for (int d=0; d<3; d++)
      {
         ParPos[d][p] = ( (real_par)rand()/RAND_MAX )*( PosMax[d] - PosMin[d] ) + PosMin[d];
         ParVel[d][p] = ( (real_par)rand()/RAND_MAX )*( VelMax[d] - VelMin[d] ) + VelMin[d];
      }
   }
   */
// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Template



#endif // #ifdef PARTICLE
