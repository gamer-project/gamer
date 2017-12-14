#include "GAMER.h"

#ifdef GRAVITY

#ifdef SERIAL
extern rfftwnd_plan     FFTW_Plan;
#else
extern rfftwnd_mpi_plan FFTW_Plan;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GreenFuncK
// Description :  Evaluate the k-space Green's function for the Possin solver with the isolated BC.
//
// Note        :  1. We only need to calculate it once during the initialization stage
//                2. The zero-padding method is implemented
//                3. Slab decomposition is assumed in FFTW
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_GreenFuncK()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Message( stderr, "OPT__BC_POT != BC_POT_ISOLATED, why do you need to calculate the Green's function !?\n" );


// 1. get the array indices used by FFTW
   const int FFT_Size[3] = { 2*NX0_TOT[0], 2*NX0_TOT[1], 2*NX0_TOT[2] };
   int local_nx, local_ny, local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

// note: total_local_size is NOT necessary to be equal to local_nx*local_ny*local_nz
   local_nx = 2*( FFT_Size[0]/2 + 1 );
   local_ny = FFT_Size[1];

#  ifdef SERIAL
   local_nz                      = FFT_Size[2];
   local_z_start                 = 0;
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = local_nx*local_ny*local_nz;
#  else
   rfftwnd_mpi_local_sizes( FFTW_Plan, &local_nz, &local_z_start, &local_ny_after_transpose,
                            &local_y_start_after_transpose, &total_local_size );

#  endif


// 2. calculate the Green's function in the real space
   const double dh0   = amr->dh[0];
   const double Coeff = -NEWTON_G*CUBE(dh0)/( (double)FFT_Size[0]*FFT_Size[1]*FFT_Size[2] );
   double x, y, z, r;
   int    kk;
   long   idx;

   GreenFuncK = new real [ total_local_size ];

   for (int k=0; k<local_nz; k++)   {  kk = k + local_z_start;
                                       z  = ( kk <= NX0_TOT[2] ) ? kk*dh0 : (FFT_Size[2]-kk)*dh0;
   for (int j=0; j<local_ny; j++)   {  y  = ( j  <= NX0_TOT[1] ) ? j *dh0 : (FFT_Size[1]-j )*dh0;
   for (int i=0; i<local_nx; i++)   {  x  = ( i  <= NX0_TOT[0] ) ? i *dh0 : (FFT_Size[0]-i )*dh0;

      r   = sqrt( x*x + y*y + z*z );
      idx = ( (long)k*local_ny + j )*local_nx + i;

      GreenFuncK[idx] = real( Coeff / r );

   }}}


// 3. reset the Green's function at the origin
// ***by setting it equal to zero, we ignore the contribution from the mass within the same cell***
   if ( MPI_Rank == 0 )    GreenFuncK[0] = GFUNC_COEFF0*Coeff/dh0;


// 4. convert the Green's function to the k space
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan, GreenFuncK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan, 1, GreenFuncK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_GreenFuncK



#endif // #ifdef GRAVITY
