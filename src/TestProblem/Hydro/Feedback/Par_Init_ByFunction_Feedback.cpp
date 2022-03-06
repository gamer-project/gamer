#include "GAMER.h"
#include <stdlib.h>

#ifdef PARTICLE

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
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Feedback( const long NPar_ThisRank, const long NPar_AllRank,
 	                           real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
          	                   real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                  	           real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   real *ParData_AllRank[PAR_NATT_TOTAL];
   for (int v=0; v<PAR_NATT_TOTAL; v++)   ParData_AllRank[v] = NULL;

// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)   ParTime[p] = Time[0];


// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
   ParData_AllRank[PAR_MASS] = new real [NPar_AllRank];
   ParData_AllRank[PAR_POSX] = new real [NPar_AllRank];
   ParData_AllRank[PAR_POSY] = new real [NPar_AllRank];
   ParData_AllRank[PAR_POSZ] = new real [NPar_AllRank];
   ParData_AllRank[PAR_VELX] = new real [NPar_AllRank];
   ParData_AllRank[PAR_VELY] = new real [NPar_AllRank];
   ParData_AllRank[PAR_VELZ] = new real [NPar_AllRank];


// initialize the particle creation time by an arbitrary negative value since it is
// only used for star particles created during evolution and is useless during initialization
// #  ifdef STAR_FORMATION
//    for (int p=0; p<NPar_AllRank; p++)    AllAttribute[Idx_ParCreTime][p] = 0.000000;
// #  endif


// set other particle attributes
// ============================================================================================================
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real *ParVel[3] = { ParVelX, ParVelY, ParVelZ };
   
   const uint RSeed     = 2;                                 // random seed
   const real MassMin   = 10.0;                              // minimum value of particle mass
   const real MassMax   = 100.0;                             // maximum value of particle mass

   srand( RSeed );

   for (long p=0; p<NPar_AllRank; p++)
   {
      ParData_AllRank[PAR_MASS][p] = MassMax * Const_Msun / UNIT_M;
   }


// example : randomly initialize
   /*
   const uint RSeed     = 2;                                         // random seed
   const real MassMin   = 1.0e-2;                                    // minimum value of particle mass
   const real MassMax   = 1.0;                                       // maximum value of particle mass
   const real PosMin[3] = { 0.0, 0.0, 0.0 };                         // minimum value of particle position
   const real PosMax[3] = { real( amr->BoxSize[0]*(1.0-1.0e-5) ),    // maximum value of particle position
                            real( amr->BoxSize[1]*(1.0-1.0e-5) ),
                            real( amr->BoxSize[2]*(1.0-1.0e-5) ) };
   const real VelMin[3] = { -1.0, -1.0, -1.0 };                      // minimum value of particle velocity
   const real VelMax[3] = { +1.0, +1.0, +1.0 };                      // maximum value of particle velocity

   srand( RSeed );
*/
   for (long p=0; p<NPar_AllRank; p++)
   {
      for (int d=0; d<3; d++)
      {
//         ParPos[d][p] = 0.5;
	 ParData_AllRank[PAR_POSX+d][p] = (double) rand() / (RAND_MAX + 1.0 );
         ParData_AllRank[PAR_VELX+d][p] = 0.0;
      }
   }

// ============================================================================================================
   }// if ( MPI_Rank == 0 )


// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, ParData_AllRank, AllAttribute );

// free resource
   if ( MPI_Rank == 0 )
      {
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParData_AllRank[v];
      }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Template



#endif // #ifdef PARTICLE
