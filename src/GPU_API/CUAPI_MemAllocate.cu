#include "CUAPI.h"

#ifdef GPU



void CUAPI_SetMemSize( int &GPU_NStream, int &Flu_GPU_NPGroup, int &Pot_GPU_NPGroup, int &Che_GPU_NPGroup,
                       int &Src_GPU_NPGroup );
int CUAPI_MemAllocate_Fluid( const int Flu_NPG, const int Pot_NPG, const int Src_NPG, const int GPU_NStream );
void CUAPI_MemFree_Fluid( const int GPU_NStream );
#ifdef GRAVITY
int CUAPI_MemAllocate_PoissonGravity( const int Pot_NPatchGroup );
void CUAPI_MemFree_PoissonGravity();
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate
// Description :  Allocate GPU memory for several global arrays
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
#  ifndef SUPPORT_GRACKLE
   int CHE_GPU_NPGROUP = NULL_INT;
#  endif

   const bool AutoNPG_Flu = ( FLU_GPU_NPGROUP <= 0 );
   const bool AutoNPG_Pot = ( POT_GPU_NPGROUP <= 0 );
   const bool AutoNPG_Che = ( CHE_GPU_NPGROUP <= 0 );
   const bool AutoNPG_Src = ( SRC_GPU_NPGROUP <= 0 );

   while ( true )
   {
      int Status = GAMER_SUCCESS;

//    set parameters related to GPU memory allocation
      CUAPI_SetMemSize( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP, CHE_GPU_NPGROUP, SRC_GPU_NPGROUP );


//    allocate memory for all GPU solvers
      if ( Status == GAMER_SUCCESS )   Status = CUAPI_MemAllocate_Fluid( FLU_GPU_NPGROUP, POT_GPU_NPGROUP, SRC_GPU_NPGROUP, GPU_NSTREAM );
#     ifdef GRAVITY
      if ( Status == GAMER_SUCCESS )   Status = CUAPI_MemAllocate_PoissonGravity( POT_GPU_NPGROUP );
#     endif


//    check if there is **any** MPI rank with Status != GAMER_SUCCESS
//    --> necessary since multiple ranks may share the same GPU (in which case the ranks that fail are nondeterministic)
//        or different nodes may have different GPUs
//    --> use MPI_BAND instead of MPI_BOR since GAMER_SUCCESS=1
#     ifndef SERIAL
      MPI_Allreduce( MPI_IN_PLACE, &Status, 1, MPI_INT, MPI_BAND, MPI_COMM_WORLD );
#     endif


//    reduce GPU_NSTREAM when running out of GPU global memory
      if ( Status == GAMER_SUCCESS )   break;
      else
      {
         if ( GPU_NSTREAM > 1 )
         {
//          free memory to avoid leakage
            CUAPI_MemFree_Fluid( GPU_NSTREAM );
#           ifdef GRAVITY
            CUAPI_MemFree_PoissonGravity();
#           endif

            const int NStreamOld = GPU_NSTREAM;
            GPU_NSTREAM /= 2;

//          set *_GPU_NPGROUP<=0 if they are not fixed manually so that they can be reset by CUAPI_SetMemSize() later
            if ( AutoNPG_Flu )   FLU_GPU_NPGROUP = -1;
            if ( AutoNPG_Pot )   POT_GPU_NPGROUP = -1;
            if ( AutoNPG_Che )   CHE_GPU_NPGROUP = -1;
            if ( AutoNPG_Src )   SRC_GPU_NPGROUP = -1;


            if ( MPI_Rank == 0 )
               Aux_Message( stderr, "WARNING : CUDA out of memory --> reducing GPU_NSTREAM from %d to %d automatically\n",
                            NStreamOld, GPU_NSTREAM );
         }

         else
            Aux_Error( ERROR_INFO, "CUDA out of memory even with GPU_NSTREAM = 1\n"
                       "        --> Try reducing *_GPU_NPGROUP in Input__Parameter manually (or set to -1) !!\n" );
      }
   } // while ( true )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_MemAllocate



#endif // #ifdef GPU
