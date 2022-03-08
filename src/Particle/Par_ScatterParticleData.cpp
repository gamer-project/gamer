#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_ScatterParticleData
// Description :  Scatter particle data from the root rank to all ranks
//
// Note        :  1. Typically called by a particle initial condition routine like Par_Init_ByFunction_*()
//                2. Do not put this file in "Particle/LoadBalance" since we still need to copy data from
//                   Data_Send[] to Data_Recv[] in the serial mode
//
// Parameter   :  NPar_ThisRank : Number of particles to be received by this MPI rank
//                NPar_AllRank  : Total number of particles in all MPI ranks
//                AttBitIdx     : Bitwise indices of the target particle attributes (e.g., _PAR_MASS | _PAR_VELX)
//                                --> A user-defined attribute with an integer index AttIntIdx returned by
//                                    AddParticleAttribute() can be converted to a bitwise index by BIDX(AttIntIdx)
//                Data_Send     : Pointer array for all particle attributes to be sent
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_AllRank]
//                                --> Target particle attributes are set by "AttBitIdx"
//                Data_Recv     : Pointer array for all particle attributes to be received
//
// Return      :  Data_Recv
//-------------------------------------------------------------------------------------------------------
void Par_ScatterParticleData( const long NPar_ThisRank, const long NPar_AllRank, const long AttBitIdx,
                              real *Data_Send[PAR_NATT_TOTAL], real *Data_Recv[PAR_NATT_TOTAL] )
{

// check integer overflow in MPI
   if ( NPar_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "Total number of particles (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                 NPar_AllRank, (long)__INT_MAX__ );


// get the number of particles in each rank and set the corresponding offsets
   int NSend[MPI_NRank], SendDisp[MPI_NRank];
   int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                             // --> (ii) is necessary for OpenMPI version < 1.7

   MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particle attributes (one at a time) from the root rank to all ranks
   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
      if ( AttBitIdx & (1L<<v) )    MPI_Scatterv( Data_Send[v], NSend, SendDisp, MPI_GAMER_REAL,
                                                  Data_Recv[v], NPar_ThisRank,   MPI_GAMER_REAL,
                                                  0, MPI_COMM_WORLD );
   }

} // FUNCTION : Par_ScatterParticleData



#endif // #ifdef PARTICLE
