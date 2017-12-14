#include "GAMER.h"

#ifdef GRAVITY



#ifdef SERIAL
rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv, FFTW_Plan_PS;    // PS : plan for calculating the power spectrum
#else
rfftwnd_mpi_plan FFTW_Plan, FFTW_Plan_Inv, FFTW_Plan_PS;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW
// Description :  Create the FFTW plans 
//-------------------------------------------------------------------------------------------------------
void Init_FFTW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// determine the FFT size
   int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };

// the zero-padding method is adopted for the isolated BC.
   if ( OPT__BC_POT == BC_POT_ISOLATED )  
      for (int d=0; d<3; d++)    FFT_Size[d] *= 2;

// check
   if ( MPI_Rank == 0 )
   for (int d=0; d<3; d++)
   {
      if ( FFT_Size[d] <= 0 )    Aux_Error( ERROR_INFO, "FFT_Size[%d] = %d < 0 !!\n", d, FFT_Size[d] );
   }


// create plans for the self-gravity solver
#  ifdef SERIAL
   FFTW_Plan     = rfftw3d_create_plan( FFT_Size[2], FFT_Size[1], FFT_Size[0], FFTW_REAL_TO_COMPLEX, 
                                        FFTW_ESTIMATE | FFTW_IN_PLACE );

   FFTW_Plan_Inv = rfftw3d_create_plan( FFT_Size[2], FFT_Size[1], FFT_Size[0], FFTW_COMPLEX_TO_REAL, 
                                        FFTW_ESTIMATE | FFTW_IN_PLACE );

#  else

   FFTW_Plan     = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, FFT_Size[2], FFT_Size[1], FFT_Size[0], 
                                            FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );

   FFTW_Plan_Inv = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, FFT_Size[2], FFT_Size[1], FFT_Size[0], 
                                            FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE );
#  endif


// create plans for calculating the power spectrum
   if ( OPT__BC_POT == BC_POT_ISOLATED )
   {
#     ifdef SERIAL
      FFTW_Plan_PS = rfftw3d_create_plan( NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], FFTW_REAL_TO_COMPLEX, 
                                          FFTW_ESTIMATE | FFTW_IN_PLACE );
#     else
      FFTW_Plan_PS = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], 
                                              FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );
#     endif
   }

   else
      FFTW_Plan_PS = FFTW_Plan;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" ); 

} // FUNCTION : Init_FFTW



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW
// Description :  Delete the FFTW plans 
//-------------------------------------------------------------------------------------------------------
void End_FFTW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );

#  ifdef SERIAL
   if ( FFTW_Plan_PS != FFTW_Plan ) 
   rfftwnd_destroy_plan    ( FFTW_Plan_PS  );

   rfftwnd_destroy_plan    ( FFTW_Plan     );
   rfftwnd_destroy_plan    ( FFTW_Plan_Inv );
#  else
   if ( FFTW_Plan_PS != FFTW_Plan ) 
   rfftwnd_mpi_destroy_plan( FFTW_Plan_PS  );

   rfftwnd_mpi_destroy_plan( FFTW_Plan     );
   rfftwnd_mpi_destroy_plan( FFTW_Plan_Inv );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_FFTW



#endif // #ifdef GRAVITY
