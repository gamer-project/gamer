#include "GAMER.h"

#ifdef GRAVITY

//output the dimensionless power spectrum
//#define DIMENSIONLESS_FORM


static void GetBasePowerSpectrum( real *RhoK, const int j_start, const int dj, double *PS_total );

#ifdef SERIAL
extern rfftwnd_plan     FFTW_Plan_PS;
#else
extern rfftwnd_mpi_plan FFTW_Plan_PS;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_BasePowerSpectrum
// Description :  Evaluate and output the base-level power spectrum by FFT
//
// Parameter   :  FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Output_BasePowerSpectrum( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// 1. determine the FFT size
   const int Nx_Padded   = NX0_TOT[0]/2+1;
   const int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };

// get the array indices using by FFTW
   int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

#  ifdef SERIAL
   local_nz                      = FFT_Size[2];
   local_z_start                 = 0;
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = 2*Nx_Padded*FFT_Size[1]*FFT_Size[2];
#  else
   rfftwnd_mpi_local_sizes( FFTW_Plan_PS, &local_nz, &local_z_start, &local_ny_after_transpose,
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


// 2. allocate memory
   const int NRecvSlice = MIN( List_z_start[MPI_Rank]+local_nz, NX0_TOT[2] ) - MIN( List_z_start[MPI_Rank], NX0_TOT[2] );

   double *PS_total     = NULL;
   real   *RhoK         = new real [ total_local_size ];                         // array storing both density and potential
   real   *SendBuf      = new real [ amr->NPatchComma[0][1]*CUBE(PS1) ];         // MPI send buffer for density and potential
   real   *RecvBuf      = new real [ NX0_TOT[0]*NX0_TOT[1]*NRecvSlice ];         // MPI recv buffer for density and potentia
   long   *SendBuf_SIdx = new long [ amr->NPatchComma[0][1]*PS1 ];               // MPI send buffer for 1D coordinate in slab
   long   *RecvBuf_SIdx = new long [ NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/SQR(PS1) ];// MPI recv buffer for 1D coordinate in slab

   int  *List_PID    [MPI_NRank];   // PID of each patch slice sent to each rank
   int  *List_k      [MPI_NRank];   // local z coordinate of each patch slice sent to each rank
   int   List_NSend  [MPI_NRank];   // size of data (density/potential) sent to each rank
   int   List_NRecv  [MPI_NRank];   // size of data (density/potential) received from each rank

   if ( MPI_Rank == 0 )    PS_total = new double [Nx_Padded];


// 3. initialize the particle density array (rho_ext) and collect particles to the target level
#  ifdef PARTICLE
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

   Prepare_PatchData_InitParticleDensityArray( 0 );

   Par_CollectParticle2OneLevel( 0, PredictPos, Time[0], SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                 TimingSendPar_No );
#  endif // #ifdef PARTICLE


// 4. rearrange data from patch to slab
   Patch2Slab( RhoK, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
               local_nz, FFT_Size, NRecvSlice, Time[0] );


// 5. evaluate the base-level power spectrum by FFT
   GetBasePowerSpectrum( RhoK, local_y_start_after_transpose, local_ny_after_transpose, PS_total );


// 6. output the power spectrum
   if ( MPI_Rank == 0 )
   {
//    check if the target file already exists
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

//    output the power spectrum
      const double WaveK0 = 2.0*M_PI/amr->BoxSize[0];
      FILE *File = fopen( FileName, "w" );

      fprintf( File, "%13s%4s%13s\n", "k", "", "Power" );

//    DC mode is not output
      for (int b=1; b<Nx_Padded; b++)     fprintf( File, "%13.6e%4s%13.6e\n", WaveK0*b, "", PS_total[b] );

      fclose( File );
   } // if ( MPI_Rank == 0 )


// 7. free memory
   delete [] RhoK;
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;
   if ( MPI_Rank == 0 )    delete [] PS_total;

// free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#  ifdef PARTICLE
   Par_CollectParticle2OneLevel_FreeMemory( 0, SibBufPatch, FaSibBufPatch );

   Prepare_PatchData_FreeParticleDensityArray( 0 );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_BasePowerSpectrum



//-------------------------------------------------------------------------------------------------------
// Function    :  GetBasePowerSpectrum
// Description :  Evaluate and base-level power spectrum by FFT
//
// Note        :  Invoked by the function "Output_BasePowerSpectrum"
//
// Parameter   :  RhoK        : Array storing the input density and output potential
//                j_start     : Starting j index
//                dj          : Size of array in the j (y) direction after the forward FFT
//                PS_total    : Power spectrum summed over all MPI ranks
//
// Return      :  PS_total
//-------------------------------------------------------------------------------------------------------
void GetBasePowerSpectrum( real *RhoK, const int j_start, const int dj, double *PS_total )
{

// check
   if ( MPI_Rank == 0  &&  PS_total == NULL )   Aux_Error( ERROR_INFO, "PS_total == NULL at the root rank !!\n" );


   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const int Nx_Padded = Nx/2 + 1;

   fftw_complex *cdata=NULL;
   double PS_local[Nx_Padded];
   long   Count_local[Nx_Padded], Count_total[Nx_Padded];
   int    bin, bin_i[Nx_Padded], bin_j[Ny], bin_k[Nz];


// forward FFT
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan_PS, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_PS, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// the data are now complex, so typecast a pointer
   cdata = (fftw_complex*) RhoK;


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
            PS_local   [bin] += double(  SQR( cdata[Idx].re ) + SQR( cdata[Idx].im )  );
            Count_local[bin] ++;
         }
      } // i,j,k
   } // i,j,k


// sum over all ranks
   MPI_Reduce( PS_local,    PS_total,    Nx_Padded, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( Count_local, Count_total, Nx_Padded, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );


// normalization: SQR(AveRho) accounts for Delta=Rho/AveRho
// --> we have assumed that the total mass in the simulation is conserved (since we don't recalculate it here)
   const double Coeff = amr->BoxSize[0]*amr->BoxSize[1]*amr->BoxSize[2] / SQR( (double)Nx*(double)Ny*(double)Nz*AveDensity_Init );
   double Norm;

#  ifdef DIMENSIONLESS_FORM
   const double k0 = 2.0*M_PI/amr->BoxSize[0];     // assuming cubic box
   double WaveK;
#  endif

   if ( MPI_Rank == 0 )
   {
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
   }

} // FUNCTION : GetBasePowerSpectrum



#endif // #ifdef GRAVITY
