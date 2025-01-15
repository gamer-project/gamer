#include "GAMER.h"

#ifdef SUPPORT_FFTW

//output the dimensionless power spectrum
//#define DIMENSIONLESS_FORM

static void GetBasePowerSpectrum( real *VarK, const int j_start, const int dj, double *PS_total, double *NormDC );

extern root_fftw::real_plan_nd FFTW_Plan_PS;




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_BasePowerSpectrum
// Description :  Evaluate and output the base-level power spectrum of the target variable by FFT
//
// Parameter   :  FileName : Name of the output file
//                TVar     : Target variable
//-------------------------------------------------------------------------------------------------------
void Output_BasePowerSpectrum( const char *FileName, const long TVar )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );

// check
// check only single field
   if ( TVar == 0  ||  TVar & (TVar-1) )
      Aux_Error( ERROR_INFO, "number of target variables is not one !!\n" );
// check cubic box
   if ( NX0_TOT[0] != NX0_TOT[1]  ||  NX0_TOT[0] != NX0_TOT[2] )
      Aux_Error( ERROR_INFO, "%s only works with CUBIC domain !!\n", __FUNCTION__ );


// 1. determine the FFT size
   const int Nx_Padded   = NX0_TOT[0]/2+1;
   const int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };

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
   rfftwnd_mpi_local_sizes( FFTW_Plan_PS, &local_nz, &local_z_start, &local_ny_after_transpose,
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


// 2. allocate memory
   const int NRecvSlice = MIN( List_z_start[MPI_Rank]+local_nz, NX0_TOT[2] ) - MIN( List_z_start[MPI_Rank], NX0_TOT[2] );

   double *PS_total     = NULL;
   real   *VarK         = (real*)root_fftw::fft_malloc( sizeof(real)*total_local_size );  // array storing data
   real   *SendBuf      = new real [ (long)amr->NPatchComma[0][1]*CUBE(PS1) ];            // MPI send buffer for data
   real   *RecvBuf      = new real [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice ];            // MPI recv buffer for data
   long   *SendBuf_SIdx = new long [ (long)amr->NPatchComma[0][1]*PS1 ];                  // MPI send buffer for 1D coordinate in slab
   long   *RecvBuf_SIdx = new long [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/SQR(PS1) ];   // MPI recv buffer for 1D coordinate in slab

   int  *List_PID    [MPI_NRank];   // PID of each patch slice sent to each rank
   int  *List_k      [MPI_NRank];   // local z coordinate of each patch slice sent to each rank
   long  List_NSend  [MPI_NRank];   // size of data sent to each rank
   long  List_NRecv  [MPI_NRank];   // size of data received from each rank
   const bool ForPoisson  = false;  // preparing the density field for the Poisson solver
   const bool InPlacePad  = true;   // pad the array for in-place real-to-complex FFT

   if ( MPI_Rank == 0 )    PS_total = new double [Nx_Padded];


// 3. initialize the particle density array (rho_ext) and collect particles to the target level
#  ifdef MASSIVE_PARTICLES
   const bool TimingSendPar_No = false;
   const bool JustCountNPar_No = false;
#  ifdef LOAD_BALANCE
   const bool PredictPos       = amr->Par->PredictPos;
   const bool SibBufPatch      = true;
   const bool FaSibBufPatch    = true;
#  else
   const bool PredictPos       = false;
   const bool SibBufPatch      = NULL_BOOL;
   const bool FaSibBufPatch    = NULL_BOOL;
#  endif

   if ( TVar == _TOTAL_DENS ) {
      Par_CollectParticle2OneLevel( 0, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ, _PAR_TYPE, PredictPos, Time[0],
                                    SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

      Prepare_PatchData_InitParticleDensityArray( 0, Time[0] );
   } // if ( TVar == _TOTAL_DENS )
#  endif // #ifdef MASSIVE_PARTICLES


// 4. rearrange data from patch to slab
   Patch2Slab( VarK, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
               local_nz, FFT_Size, NRecvSlice, Time[0], TVar, InPlacePad, ForPoisson, false );


// 5. evaluate the base-level power spectrum by FFT
   double NormDC;  // to record the FFT DC value used for normalization
   GetBasePowerSpectrum( VarK, local_y_start_after_transpose, local_ny_after_transpose, PS_total, &NormDC );


// 6. output the power spectrum
   if ( MPI_Rank == 0 )
   {
//    check if the target file already exists
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

//    output the power spectrum
      const double WaveK0 = 2.0*M_PI/amr->BoxSize[0];
      FILE *File = fopen( FileName, "w" );

      fprintf( File, "# average value (DC) used for normalization = %20.14e\n", NormDC );
      fprintf( File, "\n" );
      fprintf( File, "#%*s %*s\n", StrLen_Flt, "k", StrLen_Flt, "Power" );

//    DC mode is not output
      for (int b=1; b<Nx_Padded; b++) {
         fprintf( File, BlankPlusFormat_Flt, WaveK0*b );
         fprintf( File, BlankPlusFormat_Flt, PS_total[b] );
         fprintf( File, "\n");
      }

      fclose( File );
   } // if ( MPI_Rank == 0 )


// 7. free memory
   root_fftw::fft_free( VarK );
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;

   for (int r=0; r<MPI_NRank; r++)
   {
      free( List_PID[r] );
      free( List_k  [r] );
   }

   if ( MPI_Rank == 0 )    delete [] PS_total;

// free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#  ifdef MASSIVE_PARTICLES
   if ( TVar == _TOTAL_DENS ) {
      Par_CollectParticle2OneLevel_FreeMemory( 0, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( 0 );
   }
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_BasePowerSpectrum



//-------------------------------------------------------------------------------------------------------
// Function    :  GetBasePowerSpectrum
// Description :  Evaluate and base-level power spectrum by FFT
//
// Note        :  Invoked by the function "Output_BasePowerSpectrum"
//
// Parameter   :  VarK        : Array storing the input data
//                j_start     : Starting j index
//                dj          : Size of array in the j (y) direction after the forward FFT
//                PS_total    : Power spectrum summed over all MPI ranks
//                NormDC      : Record of the average (DC) value used for normalization of power spectrum
//
// Return      :  PS_total, NormDC
//-------------------------------------------------------------------------------------------------------
void GetBasePowerSpectrum( real *VarK, const int j_start, const int dj, double *PS_total, double *NormDC )
{

// check
   if ( MPI_Rank == 0  &&  PS_total == NULL )   Aux_Error( ERROR_INFO, "PS_total == NULL at the root rank !!\n" );


   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const int Nx_Padded = Nx/2 + 1;

   gamer_fftw::fft_complex *cdata=NULL;
   double PS_local[Nx_Padded];
   long   Count_local[Nx_Padded], Count_total[Nx_Padded];
   int    bin, bin_i[Nx_Padded], bin_j[Ny], bin_k[Nz];

   root_fftw_r2c( FFTW_Plan_PS, VarK );

// the data are now complex, so typecast a pointer
   cdata = (gamer_fftw::fft_complex*) VarK;


// set up the dimensionless wave number coefficients according to the FFTW data format
   for (int i=0; i<Nx_Padded; i++)     bin_i[i] = i;
   for (int j=0; j<Ny;        j++)     bin_j[j] = ( j <= Ny/2 ) ? j : j-Ny;
   for (int k=0; k<Nz;        k++)     bin_k[k] = ( k <= Nz/2 ) ? k : k-Nz;


// estimate the power spectrum
   long Idx;

   for (int b=0; b<Nx_Padded; b++)
   {
      PS_local   [b] = 0.0;
      Count_local[b] = 0;
   }

#  ifdef SERIAL // serial mode

   for (int k=0; k<Nz; k++)
   {
      for (int j=0; j<Ny; j++)
      for (int i=0; i<Nx_Padded; i++)
      {
         Idx = ((long)k*Ny + j)*Nx_Padded + i;

#  else // parallel mode
   int j;

   for (int jj=0; jj<dj; jj++)
   {
      j = j_start + jj;

      for (int k=0; k<Nz; k++)
      for (int i=0; i<Nx_Padded; i++)
      {
         Idx = ((long)jj*Nz + k)*Nx_Padded + i;

#  endif // #ifdef SERIAL ... else ...

//       round to nearest bin
//       bin =     int(   SQRT(   real( SQR(bin_i[i]) + SQR(bin_j[j]) + SQR(bin_k[k]) )  )   );
#        ifdef FLOAT8
         bin = lround (   SQRT(   real( SQR(bin_i[i]) + SQR(bin_j[j]) + SQR(bin_k[k]) )  )   );
#        else
         bin = lroundf(   SQRT(   real( SQR(bin_i[i]) + SQR(bin_j[j]) + SQR(bin_k[k]) )  )   );
#        endif

         if ( bin < Nx_Padded )
         {
            PS_local   [bin] += double(  SQR( c_re(cdata[Idx]) ) + SQR( c_im(cdata[Idx])  ) );
            Count_local[bin] ++;
         }
      } // i,j,k
   } // i,j,k


// sum over all ranks
   MPI_Reduce( PS_local,    PS_total,    Nx_Padded, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( Count_local, Count_total, Nx_Padded, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );


// normalization
   double Coeff;
   double AveVar;
   double Norm;

#  ifdef DIMENSIONLESS_FORM
   const double k0 = 2.0*M_PI/amr->BoxSize[0];     // assuming cubic box
   double WaveK;
#  endif

   if ( MPI_Rank == 0 )
   {
//    normalization: SQR(AveVar) accounts for Delta=Var/AveVar
      AveVar = SQRT(PS_total[0]/(double)Count_total[0]) / ( (double)Nx*(double)Ny*(double)Nz );  // from DC mode of FFT
      Coeff  = amr->BoxSize[0]*amr->BoxSize[1]*amr->BoxSize[2] / SQR( (double)Nx*(double)Ny*(double)Nz*AveVar );

      for (int b=0; b<Nx_Padded; b++)
      {
//       average
         PS_total[b] /= (double)Count_total[b];

//       normalization
#        ifdef DIMENSIONLESS_FORM
         WaveK        = b*k0;
         Norm         = Coeff*CUBE(WaveK)/(2.0*M_PI*M_PI);     // dimensionless power spectrum
#        else
         Norm         = Coeff;                                 // dimensional power spectrum [Mpc^3/h^3]
#        endif
         PS_total[b] *= Norm;
      }

      *NormDC = AveVar;
   }

} // FUNCTION : GetBasePowerSpectrum



#endif // #ifdef SUPPORT_FFTW
