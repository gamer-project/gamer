#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFile
// Description :  Initialize particle attributes from a file
//
// Note        :  1. Invoked by "Init_GAMER"
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box even for the periodic BC
//                3. Particles lying outside the active region will be removed by "Par_Aux_InitCheck"
//                   if non-periodic B.C. is adopted
//                4. Particles loaded here are only temporarily stored in this rank
//                   --> They will be redistributed when calling "Par_LB_Init_RedistributeByRectangular
//                       and LB_Init_LoadBalance"
//
// Parameter   :  None
//
// Return      :  amr->Par->Time,Mass,PosX/Y/Z,VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// verify the file size


// synchronize all particles to the physical time at the base level
// for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   amr->Par->Time[p] = Time[0];


// set other particle attributes


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFile



#endif // #ifdef PARTICLE
