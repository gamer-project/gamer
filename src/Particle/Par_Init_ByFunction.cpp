#include "GAMER.h"

#ifdef PARTICLE

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Par_Init_ByFunction();

// this function pointer may be overwritten by various test problem initializers
void (*Par_Init_ByFunction_Ptr)() = Par_Init_ByFunction;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize particle attributes by the user-specified function
//
// Note        :  1. Invoked by "Init_GAMER" using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       in this function. They will be removed automatically later when calling Par_Aux_InitCheck()
//                       in Init_GAMER()
//                       --> So as particles lying outside the "active" region
//                4. Particles set here are only temporarily stored in this rank
//                   --> They will be redistributed when calling "Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()"
//                5. For LOAD_BALANCE, the number of particles in each rank must be set in advance
//                   --> Currently it's set by "Init_Parallelization()" and stored in "amr->Par->NPar_AcPlusInac"
//
// Parameter   :  None
//
// Return      :  amr->Par->Time,Mass,PosX/Y/Z,VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// synchronize all particles to the physical time at the base level
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   amr->Par->Time[p] = Time[0];


// set other particle attributes
// ============================================================================================================
   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

// exmaple : randomly initialize
   /*
   const uint RSeed     = 2;                                         // random seed
   const real MassMin   = 1.0e-2;                                    // minimum value of particle mass
   const real MassMax   = 1.0;                                       // maximum value of particle mass
   const real PosMin[3] = { 0.0, 0.0, 0.0 };                         // minimum value of particle position
   const real PosMax[3] = { amr->BoxSize[0]*(1.0-1.0e-5),            // maximum value of particle position
                            amr->BoxSize[1]*(1.0-1.0e-5),
                            amr->BoxSize[2]*(1.0-1.0e-5) };
   const real VelMin[3] = { -1.0, -1.0, -1.0 };                      // minimum value of particle velocity
   const real VelMax[3] = { +1.0, +1.0, +1.0 };                      // maximum value of particle velocity


   srand( RSeed );

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
      Mass[p] = ( (real)rand()/RAND_MAX )*( MassMax - MassMin ) + MassMin;

      for (int d=0; d<3; d++)
      {
         Pos[d][p] = ( (real)rand()/RAND_MAX )*( PosMax[d] - PosMin[d] ) + PosMin[d];
         Vel[d][p] = ( (real)rand()/RAND_MAX )*( VelMax[d] - VelMin[d] ) + VelMin[d];
      }
   }
   */
// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
