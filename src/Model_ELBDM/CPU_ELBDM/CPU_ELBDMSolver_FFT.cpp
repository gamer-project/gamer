#include "GAMER.h"

#if ( MODEL == ELBDM  &&  defined SUPPORT_FFTW )

static void Psi_Advance_FFT( real *PsiR, real *PsiI, const int j_start, const int dj, const long PsiK_Size, const real dt );

extern root_fftw::complex_plan_nd FFTW_Plan_Psi, FFTW_Plan_Psi_Inv;  // Psi : plan for the ELBDM spectral solver




//-------------------------------------------------------------------------------------------------------
// Function    :  Psi_Advance_FFT
// Description :  Use FFT to advance the wave function by the kinetic operator (the pseudo-spectral method)
//
// Note        :  1. Invoked by CPU_ELBDMSolver_FFT()
//                2. Advance wave function by exp( -i*dt*k^2/(2*ELBDM_ETA) ) in the k-space
//
// Parameter   :  PsiR      : Array storing the real part of wave function (input and output)
//                PsiI      : Array storing the imag part of wave function (input and output)
//                j_start   : Starting j index
//                dj        : Size of array in the j (y) direction after the forward FFT
//                PsiK_Size : Size of the array "PsiK"
//                dt        : Time interval to advance solution
//-------------------------------------------------------------------------------------------------------
void Psi_Advance_FFT( real *PsiR, real *PsiI, const int j_start, const int dj, const long PsiK_Size, const real dt )
{

   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const real dh       = amr->dh[0];
   const real Dt_2Eta  = (real)0.5*dt/ELBDM_ETA;

   real PsiKR, PsiKI, DtKK_2Eta;
   gamer_fftw::fft_complex *PsiK;
   PsiK = (gamer_fftw::fft_complex*)root_fftw::fft_malloc( PsiK_Size*sizeof(gamer_fftw::fft_complex) );

   for (long t=0; t<PsiK_Size; t++)
   {
      c_re(PsiK[t]) = PsiR[t];
      c_im(PsiK[t]) = PsiI[t];
   }


// forward FFT
   root_fftw_c2c( FFTW_Plan_Psi, PsiK );


// set up the dimensional wave number
// and the k-space evolution phase (K_i*K_i)*dt/(2*ELBDM_ETA), where i=x,y,z
   real Kx[Nx], Ky[Ny], Kz[Nz];
   real DtKxKx_2Eta[Nx], DtKyKy_2Eta[Ny], DtKzKz_2Eta[Nz];

   for (int i=0; i<Nx; i++)
   {
      Kx[i]          = ( i <= Nx/2 ) ? 2.0*M_PI/(Nx*dh)*i : 2.0*M_PI/(Nx*dh)*(i-Nx);
      DtKxKx_2Eta[i] = SQR( Kx[i] )*Dt_2Eta;
   }
   for (int j=0; j<Ny; j++)
   {
      Ky[j]          = ( j <= Ny/2 ) ? 2.0*M_PI/(Ny*dh)*j : 2.0*M_PI/(Ny*dh)*(j-Ny);
      DtKyKy_2Eta[j] = SQR( Ky[j] )*Dt_2Eta;
   }
   for (int k=0; k<Nz; k++)
   {
      Kz[k]          = ( k <= Nz/2 ) ? 2.0*M_PI/(Nz*dh)*k : 2.0*M_PI/(Nz*dh)*(k-Nz);
      DtKzKz_2Eta[k] = SQR( Kz[k] )*Dt_2Eta;
   }


// multiply the wave function by exp( -i*dt*k^2/(2*ELBDM_ETA) ) in the k-space
#  ifdef SERIAL // serial mode

#  pragma omp parallel for schedule( runtime )
   for (int k=0; k<Nz; k++)
   {
      for (int j=0; j<Ny; j++)
      for (int i=0; i<Nx; i++)
      {
         const long ID = ((long)k*Ny + j)*Nx + i;

#  else // parallel mode

#  pragma omp parallel for schedule( runtime )
   for (int jj=0; jj<dj; jj++)
   {
      const int j = j_start + jj;

      for (int k=0; k<Nz; k++)
      for (int i=0; i<Nx; i++)
      {
         const long ID = ((long)jj*Nz + k)*Nx + i;

#  endif // #ifdef SERIAL ... else ...

         const real PsiKR = c_re(PsiK[ID]);
         const real PsiKI = c_im(PsiK[ID]);
         const real DtKK_2Eta = DtKxKx_2Eta[i] + DtKyKy_2Eta[j] + DtKzKz_2Eta[k];

         c_re(PsiK[ID]) =  PsiKR*COS( DtKK_2Eta ) + PsiKI*SIN( DtKK_2Eta );
         c_im(PsiK[ID]) =  PsiKI*COS( DtKK_2Eta ) - PsiKR*SIN( DtKK_2Eta );
      } // i,j,k
   } // i,j,k


// backward FFT
   root_fftw_c2c( FFTW_Plan_Psi_Inv, PsiK );

// normalization
   const real norm = 1.0 / ( (real)Nx*Ny*Nz );

   for (long t=0; t<PsiK_Size; t++)
   {
      PsiR[t] = c_re(PsiK[t]) * norm;
      PsiI[t] = c_im(PsiK[t]) * norm;
   }

   root_fftw::fft_free( PsiK );

} // FUNCTION : Psi_Advance_FFT



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ELBDMSolver_FFT
// Description :  CPU ELBDM kinetic solver of base-level wave function using FFT (the pseudo-spectral method)
//
// Note        :  1. Work with the option ELBDM_BASE_SPECTRAL
//                2. Invoked by Flu_AdvanceDt()
//
// Parameter   :  dt       : Time interval to advance solution
//                PrepTime : Physical time for preparing the wave function
//                SaveSg   : Sandglass to store the updated data
//-------------------------------------------------------------------------------------------------------
void CPU_ELBDMSolver_FFT( const real dt, const double PrepTime, const int SaveSg )
{

// determine the FFT size
   const int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };


// get the array indices using by FFTW
   mpi_index_int local_nx, local_ny, local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

// note: total_local_size is NOT necessarily equal to local_nx*local_ny*local_nz
   local_nx = FFT_Size[0];
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
   fftwnd_mpi_local_sizes( FFTW_Plan_Psi, &local_nz, &local_z_start, &local_ny_after_transpose,
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


// allocate memory
   const bool InPlacePad_No = false;   // not pad the array for in-place real-to-complex FFT
   const bool ForPoisson_No = false;   // not for the Poisson solver
   const int  NRecvSlice    = MIN( List_z_start[MPI_Rank]+local_nz, NX0_TOT[2] ) - MIN( List_z_start[MPI_Rank], NX0_TOT[2] );

   real *PsiR         = (real*)root_fftw::fft_malloc( sizeof(real)*total_local_size ); // array storing real and imaginary parts of wave function
   real *PsiI         = (real*)root_fftw::fft_malloc( sizeof(real)*total_local_size );
   real *SendBuf      = new real [ (long)amr->NPatchComma[0][1]*CUBE(PS1) ];           // MPI send buffer
   real *RecvBuf      = new real [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice ];           // MPI recv buffer
   long *SendBuf_SIdx = new long [ (long)amr->NPatchComma[0][1]*PS1 ];                 // MPI send buffer for 1D coordinate in slab
   long *RecvBuf_SIdx = new long [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/SQR(PS1) ];  // MPI recv buffer for 1D coordinate in slab

   int  *List_PID_R  [MPI_NRank];   // PID of each patch slice sent to each rank for the real part
   int  *List_k_R    [MPI_NRank];   // local z coordinate of each patch slice sent to each rank for the real part
   int  *List_PID_I  [MPI_NRank];   // PID of each patch slice sent to each rank for the imag part
   int  *List_k_I    [MPI_NRank];   // local z coordinate of each patch slice sent to each rank for the imag part
   long  List_NSend  [MPI_NRank];   // size of data sent to each rank
   long  List_NRecv  [MPI_NRank];   // size of data received from each rank


// rearrange data from patch to slab
   Patch2Slab( PsiR, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID_R, List_k_R, List_NSend, List_NRecv, List_z_start,
               local_nz, FFT_Size, NRecvSlice, PrepTime, _REAL, InPlacePad_No, ForPoisson_No, false );
   Patch2Slab( PsiI, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID_I, List_k_I, List_NSend, List_NRecv, List_z_start,
               local_nz, FFT_Size, NRecvSlice, PrepTime, _IMAG, InPlacePad_No, ForPoisson_No, false );


// advance wave function by exp( -i*dt*k^2/(2*ELBDM_ETA) ) in the k-space using FFT
   Psi_Advance_FFT( PsiR, PsiI, local_y_start_after_transpose, local_ny_after_transpose, total_local_size, dt );


// rearrange data from slab back to patch
   Slab2Patch( PsiR, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID_R, List_k_R, List_NRecv, List_NSend,
               local_nz, FFT_Size, NRecvSlice, _REAL, InPlacePad_No );
   Slab2Patch( PsiI, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID_I, List_k_I, List_NRecv, List_NSend,
               local_nz, FFT_Size, NRecvSlice, _IMAG, InPlacePad_No );


// update density according to the updated wave function
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[0][1]; PID++)
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
      const real NewReal = amr->patch[SaveSg][0][PID]->fluid[REAL][k][j][i];
      const real NewImag = amr->patch[SaveSg][0][PID]->fluid[IMAG][k][j][i];
      const real NewDens = SQR( NewReal ) + SQR( NewImag );

      amr->patch[SaveSg][0][PID]->fluid[DENS][k][j][i] = NewDens;
   } // PID,i,j,k


// free memory
   root_fftw::fft_free( PsiR );
   root_fftw::fft_free( PsiI );
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;

} // FUNCTION : CPU_ELBDMSolver_FFT



#endif // #if ( MODEL == ELBDM  &&  defined SUPPORT_FFTW )
