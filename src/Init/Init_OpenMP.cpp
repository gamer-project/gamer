#include "Copyright.h"
#include "GAMER.h"

#ifdef OPENMP



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_OpenMP
// Description :  Initialize OpenMP
//
// Note        :  The schedule clause specified here only applies to all solver-related routines
//                (e.g., XXX_Prepare_XXX, XXX_Close, CPU_XXXSolver_XXX, Flu_FixUp, Flu_Restrict, Flag_Real,
//                 Hydro_GetTimeStep_Gravity, ...)
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
// omp_set_schedule( omp_sched_dynamic, chunk_size );
   omp_set_schedule( omp_sched_guided,  chunk_size );
// omp_set_schedule( omp_sched_auto,    chunk_size );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_OpenMP



#endif // #ifdef OPENMP
