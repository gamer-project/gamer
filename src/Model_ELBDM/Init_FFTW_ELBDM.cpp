#include "GAMER.h"

#if ( MODEL == ELBDM )



#ifdef SERIAL
fftwnd_plan      FFTW_Plan_Psi, FFTW_Plan_Psi_Inv;
#else
fftwnd_mpi_plan  FFTW_Plan_Psi, FFTW_Plan_Psi_Inv;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW_ELBDM
// Description :  Create the FFTW plans for ELBDM solver
//-------------------------------------------------------------------------------------------------------
void Init_FFTW_ELBDM()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// determine the FFT size
   int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };

// check
   if ( MPI_Rank == 0 )
   for (int d=0; d<3; d++)
   {
      if ( FFT_Size[d] <= 0 )    Aux_Error( ERROR_INFO, "FFT_Size[%d] = %d < 0 !!\n", d, FFT_Size[d] );
   }


// create plans for the ELBDM solver
#  ifdef SERIAL
   FFTW_Plan_Psi     = fftw3d_create_plan( NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], FFTW_FORWARD,
                                        FFTW_ESTIMATE | FFTW_IN_PLACE );

   FFTW_Plan_Psi_Inv = fftw3d_create_plan( NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], FFTW_BACKWARD,
                                        FFTW_ESTIMATE | FFTW_IN_PLACE );

#  else

   FFTW_Plan_Psi     = fftw3d_mpi_create_plan( MPI_COMM_WORLD, NX0_TOT[2], NX0_TOT[1], NX0_TOT[0],
                                            FFTW_FORWARD, FFTW_ESTIMATE );

   FFTW_Plan_Psi_Inv = fftw3d_mpi_create_plan( MPI_COMM_WORLD, NX0_TOT[2], NX0_TOT[1], NX0_TOT[0],
                                            FFTW_BACKWARD, FFTW_ESTIMATE );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_FFTW_ELBDM



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW_ELBDM
// Description :  Delete the FFTW plans for ELBDM solver
//-------------------------------------------------------------------------------------------------------
void End_FFTW_ELBDM()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );

#  ifdef SERIAL
   fftwnd_destroy_plan     ( FFTW_Plan_Psi     );
   fftwnd_destroy_plan     ( FFTW_Plan_Psi_Inv );
#  else
   fftwnd_mpi_destroy_plan ( FFTW_Plan_Psi     );
   fftwnd_mpi_destroy_plan ( FFTW_Plan_Psi_Inv );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_FFTW_ELBDM



#endif // #if ( MODEL == ELBDM )
