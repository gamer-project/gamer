#include "GAMER.h"

#ifdef OPENMP



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_OpenMP
// Description :  Initialize OpenMP
//
// Note        :  1. The schedule clause specified here only applies to those using the "runtime" schedule
//                   --> We adopt the "runtime" schedule if the workload of each iteration may be different
//                   --> We adopt the "static" schedule if the workload of each iteration is the same
//                   --> The schedule of most particle routines is however controlled by the symbolic constant
//                       "PAR_OMP_SCHED" and "PAR_OMP_SCHED_CHUNK" defined in "Macro.h"
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_OpenMP()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// numbef of OMP threads
   omp_set_num_threads( OMP_NTHREAD );

// enable/disable nested parallelization
   omp_set_nested( false );

// schedule
   const int chunk_size = 1;
// omp_set_schedule( omp_sched_static,  chunk_size );
   omp_set_schedule( omp_sched_dynamic, chunk_size );
// omp_set_schedule( omp_sched_guided,  chunk_size );
// omp_set_schedule( omp_sched_auto,    chunk_size );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_OpenMP



#endif // #ifdef OPENMP
