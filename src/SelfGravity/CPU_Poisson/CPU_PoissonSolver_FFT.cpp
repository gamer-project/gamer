#include "GAMER.h"

#if ( defined GRAVITY  &&  defined SUPPORT_FFTW )

static void FFT_Periodic( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const long RhoK_Size );
static void FFT_Isolated( real *RhoK, const real *gFuncK, const real Poi_Coeff, const long RhoK_Size );

extern root_fftw::real_plan_nd FFTW_Plan_Poi, FFTW_Plan_Poi_Inv;




//-------------------------------------------------------------------------------------------------------
// Function    :  FFT_Periodic
// Description :  Evaluate the gravitational potential by FFT for the periodic BC
//
// Note        :  Effect from the homogenerous background density (DC) will be ignored by setting the k=0 mode
//                equal to zero
//
// Parameter   :  RhoK      : Array storing the input density and output potential
//                Poi_Coeff : Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)
//                j_start   : Starting j index
//                dj        : Size of array in the j (y) direction after the forward FFT
//                RhoK_Size : Size of the array "RhoK"
//-------------------------------------------------------------------------------------------------------
void FFT_Periodic( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const long RhoK_Size )
{

   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const int Nx_Padded = Nx/2 + 1;
   const real dh       = amr->dh[0];
   real Deno;
   gamer_fftw::fft_complex *cdata;


// forward FFT
   root_fftw_r2c( FFTW_Plan_Poi, RhoK );

// the data are now complex, so typecast a pointer
   cdata = (gamer_fftw::fft_complex*) RhoK;


// set up the dimensionless wave number and the corresponding sin(k)^2 function
   real kx[Nx_Padded], ky[Ny], kz[Nz];
   real sinkx2[Nx_Padded], sinky2[Ny], sinkz2[Nz];

   for (int i=0; i<Nx_Padded; i++) {   kx    [i] = 2.0*M_PI/Nx*i;
                                       sinkx2[i] = SQR(  SIN( (real)0.5*kx[i] )  );    }
   for (int j=0; j<Ny;        j++) {   ky    [j] = ( j <= Ny/2 ) ? 2.0*M_PI/Ny*j : 2.0*M_PI/Ny*(j-Ny);
                                       sinky2[j] = SQR(  SIN( (real)0.5*ky[j] )  );    }
   for (int k=0; k<Nz;        k++) {   kz    [k] = ( k <= Nz/2 ) ? 2.0*M_PI/Nz*k : 2.0*M_PI/Nz*(k-Nz);
                                       sinkz2[k] = SQR(  SIN( (real)0.5*kz[k] )  );    }


// divide the Rho_K by -k^2
   long ID;
#  ifdef SERIAL // serial mode

   for (int k=0; k<Nz; k++)
   {
      for (int j=0; j<Ny; j++)
      for (int i=0; i<Nx_Padded; i++)
      {
         ID = ((long)k*Ny + j)*Nx_Padded + i;

#  else // parallel mode
   int j;

   for (int jj=0; jj<dj; jj++)
   {
      j = j_start + jj;

      for (int k=0; k<Nz; k++)
      for (int i=0; i<Nx_Padded; i++)
      {
         ID = ((long)jj*Nz + k)*Nx_Padded + i;

#  endif // #ifdef SERIAL ... else ...


//       this form is more consistent with the "second-order discrete" Laplacian operator
         Deno = -4.0 * ( sinkx2[i] + sinky2[j] + sinkz2[k] );
//       Deno = -( kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k] );


//       remove the DC mode
         if ( Deno == 0.0 )
         {
            c_re(cdata[ID]) = 0.0;
            c_im(cdata[ID]) = 0.0;
         }

         else
         {
            c_re(cdata[ID]) =  c_re(cdata[ID]) * Poi_Coeff / Deno;
            c_im(cdata[ID]) =  c_im(cdata[ID]) * Poi_Coeff / Deno;
         }
      } // i,j,k
   } // i,j,k


// backward FFT
   root_fftw_c2r( FFTW_Plan_Poi_Inv, RhoK );

// normalization
   const real norm = dh*dh / ( (real)Nx*Ny*Nz );

   for (long t=0; t<RhoK_Size; t++)  RhoK[t] *= norm;

} // FUNCTION : FFT_Periodic



//-------------------------------------------------------------------------------------------------------
// Function    :  FFT_Isolated
// Description :  Evaluate the gravitational potential by FFT for the isolated BC
//
// Note        :  1. Green's function in the k space has been set by Init_GreenFuncK()
//                2. 4*PI*NEWTON_G and FFT normalization coefficient has been included in gFuncK
//                   --> The only coefficient that hasn't been taken into account is the scale factor in the comoving frame
//
// Parameter   :  RhoK      : Array storing the input density and output potential
//                Poi_Coeff : Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)
//                RhoK_Size : Size of the array "RhoK"
//-------------------------------------------------------------------------------------------------------
void FFT_Isolated( real *RhoK, const real *gFuncK, const real Poi_Coeff, const long RhoK_Size )
{

   gamer_fftw::fft_complex *RhoK_cplx   = (gamer_fftw::fft_complex *)RhoK;
   gamer_fftw::fft_complex *gFuncK_cplx = (gamer_fftw::fft_complex *)gFuncK;
   gamer_fftw::fft_complex  Temp_cplx;


// forward FFT
   root_fftw_r2c( FFTW_Plan_Poi, RhoK );


// multiply density and Green's function in the k space
   const long RhoK_Size_cplx = RhoK_Size/2;

   for (long t=0; t<RhoK_Size_cplx; t++)
   {
      c_re(Temp_cplx) = c_re(RhoK_cplx[t]);
      c_im(Temp_cplx) = c_im(RhoK_cplx[t]);

      c_re(RhoK_cplx[t]) = c_re(Temp_cplx)*c_re(gFuncK_cplx[t]) - c_im(Temp_cplx)*c_im(gFuncK_cplx[t]);
      c_im(RhoK_cplx[t]) = c_re(Temp_cplx)*c_im(gFuncK_cplx[t]) + c_im(Temp_cplx)*c_re(gFuncK_cplx[t]);
   }


// backward FFT
   root_fftw_c2r( FFTW_Plan_Poi_Inv, RhoK );

// effect of "4*PI*NEWTON_G" has been included in gFuncK, but the scale factor in the comoving frame hasn't
#  ifdef COMOVING
   const real Coeff = Poi_Coeff / ( 4.0*M_PI*NEWTON_G );    // == Time[0] == scale factor at the base level

   for (long t=0; t<RhoK_Size; t++)  RhoK[t] *= Coeff;
#  endif

} // FUNCTION : FFT_Isolated



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonSolver_FFT
// Description :  Evaluate the base-level potential by FFT
//
// Note        :  Work with both periodic and isolated BC's
//
// Parameter   :  Poi_Coeff : Coefficient in front of the RHS in the Poisson eq.
//                SaveSg    : Sandglass to store the updated data
//                PrepTime  : Physical time for preparing the density field
//-------------------------------------------------------------------------------------------------------
void CPU_PoissonSolver_FFT( const real Poi_Coeff, const int SaveSg, const double PrepTime )
{

// determine the FFT size (the zero-padding method is adopted for the isolated BC)
   int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };

   if ( OPT__BC_POT == BC_POT_ISOLATED )
      for (int d=0; d<3; d++)    FFT_Size[d] *= 2;


// get the array indices using by FFTW
   mpi_index_int local_nx, local_ny, local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

// note: total_local_size is NOT necessarily equal to local_nx*local_ny*local_nz
   local_nx = 2*( FFT_Size[0]/2 + 1 );
   local_ny = FFT_Size[1];

#  ifdef SERIAL
   local_nz                      = FFT_Size[2];
   local_z_start                 = 0;
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = local_nx*local_ny*local_nz;
#  else // #ifdef SERIAL
#  if ( SUPPORT_FFTW == FFTW3 )
   total_local_size = fftw_mpi_local_size_3d_transposed( FFT_Size[2], local_ny, local_nx, MPI_COMM_WORLD,
                                                         &local_nz, &local_z_start, &local_ny_after_transpose,
                                                         &local_y_start_after_transpose );
#  else
   rfftwnd_mpi_local_sizes( FFTW_Plan_Poi, &local_nz, &local_z_start, &local_ny_after_transpose,
                            &local_y_start_after_transpose, &total_local_size );
#  endif
#  endif // #ifdef SERIAL ... else ...

// check integer overflow (assuming local_nx*local_ny*local_nz ~ total_local_size)
   const long local_nxyz = (long)local_nx*(long)local_ny*(long)local_nz;

   if ( local_nx < 0  ||  local_ny < 0  ||  local_nz < 0 )
      Aux_Error( ERROR_INFO, "local_nx/y/z (%ld, %ld, %ld) < 0 for FFT !!\n", local_nx, local_ny, local_nz );

   if (  ( sizeof(mpi_index_int) == sizeof(int) && local_nxyz > __INT_MAX__ )  ||  total_local_size < 0  )
      Aux_Error( ERROR_INFO, "local_nx*local_ny*local_nz = %d*%d*%d = %ld > __INT_MAX__ (%d)\n"
                     "        and/or total_local_size (%ld) < 0 for FFT, suggesting integer overflow !!\n"
                     "        --> Try using more MPI processes or switching to FFTW3\n",
                 local_nx, local_ny, local_nz, local_nxyz, __INT_MAX__, total_local_size );


// collect "local_nz" from all ranks and set the corresponding list "List_z_start"
   int List_nz     [MPI_NRank  ];   // slab thickness of each rank in the FFTW slab decomposition
   int List_z_start[MPI_NRank+1];   // starting z coordinate of each rank in the FFTW slab decomposition

   const int local_nz_int = local_nz;  // necessary since "mpi_index_int" maps to "long int" for FFTW3
   MPI_Allgather( &local_nz_int, 1, MPI_INT, List_nz, 1, MPI_INT, MPI_COMM_WORLD );

   List_z_start[0] = 0;
   for (int r=0; r<MPI_NRank; r++)  List_z_start[r+1] = List_z_start[r] + List_nz[r];

   if ( List_z_start[MPI_NRank] != FFT_Size[2] )
      Aux_Error( ERROR_INFO, "List_z_start[%d] (%d) != expectation (%d) !!\n",
                 MPI_NRank, List_z_start[MPI_NRank], FFT_Size[2] );


// allocate memory (properly taking into account the zero-padding regions, where no data need to be exchanged)
   const int NRecvSlice = MIN( List_z_start[MPI_Rank]+local_nz, NX0_TOT[2] ) - MIN( List_z_start[MPI_Rank], NX0_TOT[2] );

   real *RhoK         = (real*)root_fftw::fft_malloc( sizeof(real)*total_local_size ); // array storing both density and potential
   real *SendBuf      = new real [ (long)amr->NPatchComma[0][1]*CUBE(PS1) ];           // MPI send buffer for density and potential
   real *RecvBuf      = new real [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice ];           // MPI recv buffer for density and potentia
   long *SendBuf_SIdx = new long [ (long)amr->NPatchComma[0][1]*PS1 ];                 // MPI send buffer for 1D coordinate in slab
   long *RecvBuf_SIdx = new long [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/SQR(PS1) ];  // MPI recv buffer for 1D coordinate in slab

   int  *List_PID    [MPI_NRank];   // PID of each patch slice sent to each rank
   int  *List_k      [MPI_NRank];   // local z coordinate of each patch slice sent to each rank
   long  List_NSend  [MPI_NRank];   // size of data (density/potential) sent to each rank
   long  List_NRecv  [MPI_NRank];   // size of data (density/potential) received from each rank
   const bool ForPoisson  = true;   // preparing the density field for the Poisson solver
   const bool InPlacePad  = true;   // pad the array for in-place real-to-complex FFT


// initialize RhoK as zeros for the isolated BC where the zero-padding method is adopted
   if ( OPT__BC_POT == BC_POT_ISOLATED )
      for (long t=0; t<total_local_size; t++)   RhoK[t] = (real)0.0;


// rearrange data from patch to slab
   Patch2Slab( RhoK, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
               local_nz, FFT_Size, NRecvSlice, PrepTime, _TOTAL_DENS, InPlacePad, ForPoisson, OPT__GRAVITY_EXTRA_MASS );


// evaluate potential by FFT
   if      ( OPT__BC_POT == BC_POT_PERIODIC )
      FFT_Periodic( RhoK, Poi_Coeff, local_y_start_after_transpose, local_ny_after_transpose, total_local_size );

   else if ( OPT__BC_POT == BC_POT_ISOLATED )
      FFT_Isolated( RhoK, GreenFuncK, Poi_Coeff, total_local_size );

   else
      Aux_Error( ERROR_INFO, "unsupported paramter %s = %d !!\n", "OPT__BC_POT", OPT__BC_POT );


// rearrange data from slab back to patch
   Slab2Patch( RhoK, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID, List_k, List_NRecv, List_NSend,
               local_nz, FFT_Size, NRecvSlice, _POTE, InPlacePad );


   root_fftw::fft_free( RhoK );
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;

} // FUNCTION : CPU_PoissonSolver_FFT



#endif // #if ( defined GRAVITY  &&  defined SUPPORT_FFTW )
