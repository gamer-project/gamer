#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Function
// Description :  Initialize the particle position and velocity
//
// Note        :  1. Invoked by "Init_GAMER"
//                2. Periodicity should be taken care of in this function
//                3. Particles lying outside the active region will be removed by "Par_Aux_InitCheck"
//                   if non-periodic B.C. is adopted
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass,PosX/Y/Z,VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_Function()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// synchronize all particles to the physical time at the base level
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   amr->Par->Time[p] = Time[0];


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };


// ============================================================================================================
// exmaple : randomly initialize
   /*
   const uint RSeed     = 2;                                         // random seed
   const real MassMin   = TINY_VALUE;                                // minimum value of particle mass
   const real MassMax   = 1.0;                                       // maximum value of particle mass
   const real PosMin[3] = { 0.0, 0.0, 0.0 };                         // minimum value of particle position
   const real PosMax[3] = { amr->BoxSize[0]*(1.0-TINY_VALUE),        // maximum value of particle position
                            amr->BoxSize[1]*(1.0-TINY_VALUE),
                            amr->BoxSize[2]*(1.0-TINY_VALUE) };
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

} // FUNCTION : Par_Init_Function



#endif // #ifdef PARTICLE
