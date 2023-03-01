#include "GAMER.h"

#if ( defined GRAVITY  &&  defined SUPPORT_FFTW )



static void FFT_Periodic( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const int RhoK_Size );
static void FFT_Isolated( real *RhoK, const real *gFuncK, const real Poi_Coeff, const int RhoK_Size );

#ifdef SERIAL
extern rfftwnd_plan     FFTW_Plan_Poi, FFTW_Plan_Poi_Inv;
#else
extern rfftwnd_mpi_plan FFTW_Plan_Poi, FFTW_Plan_Poi_Inv;
#endif




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
void FFT_Periodic( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const int RhoK_Size )
{

   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const int Nx_Padded = Nx/2 + 1;
   const real dh       = amr->dh[0];
   real Deno;
   fftw_complex *cdata;


// forward FFT
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan_Poi, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Poi, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// the data are now complex, so typecast a pointer
   cdata = (fftw_complex*) RhoK;


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
            cdata[ID].re = 0.0;
            cdata[ID].im = 0.0;
         }

         else
         {
            cdata[ID].re =  cdata[ID].re * Poi_Coeff / Deno;
            cdata[ID].im =  cdata[ID].im * Poi_Coeff / Deno;
         }
      } // i,j,k
   } // i,j,k


// backward FFT
#  ifdef SERIAL
   rfftwnd_one_complex_to_real( FFTW_Plan_Poi_Inv, cdata, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Poi_Inv, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// normalization
   const real norm = dh*dh / ( (real)Nx*Ny*Nz );

   for (int t=0; t<RhoK_Size; t++)  RhoK[t] *= norm;

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
void FFT_Isolated( real *RhoK, const real *gFuncK, const real Poi_Coeff, const int RhoK_Size )
{

   fftw_complex *RhoK_cplx   = (fftw_complex *)RhoK;
   fftw_complex *gFuncK_cplx = (fftw_complex *)gFuncK;
   fftw_complex  Temp_cplx;


// forward FFT
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan_Poi, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Poi, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// multiply density and Green's function in the k space
   const int RhoK_Size_cplx = RhoK_Size/2;

   for (int t=0; t<RhoK_Size_cplx; t++)
   {
      Temp_cplx = RhoK_cplx[t];

      RhoK_cplx[t].re = Temp_cplx.re*gFuncK_cplx[t].re - Temp_cplx.im*gFuncK_cplx[t].im;
      RhoK_cplx[t].im = Temp_cplx.re*gFuncK_cplx[t].im + Temp_cplx.im*gFuncK_cplx[t].re;
   }


// backward FFT
#  ifdef SERIAL
   rfftwnd_one_complex_to_real( FFTW_Plan_Poi_Inv, RhoK_cplx, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Poi_Inv, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// effect of "4*PI*NEWTON_G" has been included in gFuncK, but the scale factor in the comoving frame hasn't
#  ifdef COMOVING
   const real Coeff = Poi_Coeff / ( 4.0*M_PI*NEWTON_G );    // == Time[0] == scale factor at the base level

   for (int t=0; t<RhoK_Size; t++)  RhoK[t] *= Coeff;
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
   int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

#  ifdef SERIAL
   local_nz                      = FFT_Size[2];
   local_z_start                 = 0;
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = 2*(FFT_Size[0]/2+1)*FFT_Size[1]*FFT_Size[2];
#  else
   rfftwnd_mpi_local_sizes( FFTW_Plan_Poi, &local_nz, &local_z_start, &local_ny_after_transpose,
                            &local_y_start_after_transpose, &total_local_size );
#  endif


// collect "local_nz" from all ranks and set the corresponding list "List_z_start"
   int List_nz     [MPI_NRank  ];   // slab thickness of each rank in the FFTW slab decomposition
   int List_z_start[MPI_NRank+1];   // starting z coordinate of each rank in the FFTW slab decomposition

   MPI_Allgather( &local_nz, 1, MPI_INT, List_nz, 1, MPI_INT, MPI_COMM_WORLD );

   List_z_start[0] = 0;
   for (int r=0; r<MPI_NRank; r++)  List_z_start[r+1] = List_z_start[r] + List_nz[r];

   if ( List_z_start[MPI_NRank] != FFT_Size[2] )
      Aux_Error( ERROR_INFO, "List_z_start[%d] (%d) != expectation (%d) !!\n",
                 MPI_NRank, List_z_start[MPI_NRank], FFT_Size[2] );


// allocate memory (properly taking into account the zero-padding regions, where no data need to be exchanged)
   const int NRecvSlice = MIN( List_z_start[MPI_Rank]+local_nz, NX0_TOT[2] ) - MIN( List_z_start[MPI_Rank], NX0_TOT[2] );

   real *RhoK         = new real [ total_local_size ];                           // array storing both density and potential
   real *SendBuf      = new real [ (long)amr->NPatchComma[0][1]*CUBE(PS1) ];     // MPI send buffer for density and potential
   real *RecvBuf      = new real [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice ];     // MPI recv buffer for density and potentia
   long *SendBuf_SIdx = new long [ amr->NPatchComma[0][1]*PS1 ];                 // MPI send buffer for 1D coordinate in slab
   long *RecvBuf_SIdx = new long [ NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/SQR(PS1) ];  // MPI recv buffer for 1D coordinate in slab

   int  *List_PID    [MPI_NRank];   // PID of each patch slice sent to each rank
   int  *List_k      [MPI_NRank];   // local z coordinate of each patch slice sent to each rank
   int   List_NSend  [MPI_NRank];   // size of data (density/potential) sent to each rank
   int   List_NRecv  [MPI_NRank];   // size of data (density/potential) received from each rank
   const bool ForPoisson  = true;   // preparing the density field for the Poisson solver
   const bool InPlacePad  = true;   // pad the array for in-place real-to-complex FFT


// initialize RhoK as zeros for the isolated BC where the zero-padding method is adopted
   if ( OPT__BC_POT == BC_POT_ISOLATED )
      for (int t=0; t<total_local_size; t++)    RhoK[t] = (real)0.0;


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


   delete [] RhoK;
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;

} // FUNCTION : CPU_PoissonSolver_FFT



#endif // #if ( defined GRAVITY  &&  defined SUPPORT_FFTW )
