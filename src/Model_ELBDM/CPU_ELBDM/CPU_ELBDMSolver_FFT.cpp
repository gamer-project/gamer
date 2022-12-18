#include "GAMER.h"

#if ( MODEL == ELBDM )



static int ZIndex2Rank_Psi( const int IndexZ, const int *List_z_start, const int TRank_Guess );

#ifdef SERIAL
extern fftwnd_plan     FFTW_Plan_Psi, FFTW_Plan_Psi_Inv;
#else
extern fftwnd_mpi_plan FFTW_Plan_Psi, FFTW_Plan_Psi_Inv;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Patch2Slab_Psi
// Description :  Patch-based data --> slab domain decomposition (for wavefunction)
//
// Parameter   :  Psi            : Array storing wavefunction
//                SendBuf_Psi    : Sending MPI buffer of wavefunction
//                RecvBuf_Psi    : Receiving MPI buffer of wavefunction
//                SendBuf_SIdx   : Sending MPI buffer of 1D coordinate in slab
//                RecvBuf_SIdx   : Receiving MPI buffer of 1D coordinate in slab
//                List_PID       : PID of each patch slice sent to each rank
//                List_k         : Local z coordinate of each patch slice sent to each rank
//                List_NSend_Psi : Size of wavefunction data sent to each rank
//                List_NRecv_Psi : Size of wavefunction data received from each rank
//                List_z_start   : Starting z coordinate of each rank in the FFTW slab decomposition
//                local_nz       : Slab thickness of this MPI rank
//                FFT_Size       : Size of the FFT operation
//                NRecvSlice     : Total number of z slices received from other ranks
//                Time           : Physical time for preparing the wavefunction
//                Target         : Target variable to be prepared (REAL or IMAG or DENS)
//-------------------------------------------------------------------------------------------------------
void Patch2Slab_Psi( real *Psi, real *SendBuf_Psi, real *RecvBuf_Psi, long *SendBuf_SIdx, long *RecvBuf_SIdx,
                 int **List_PID, int **List_k, int *List_NSend_Psi, int *List_NRecv_Psi,
                 const int *List_z_start, const int local_nz, const int FFT_Size[], const int NRecvSlice,
                 const double Time, const int Target )
{

// check
#  ifdef GAMER_DEBUG
   if ( List_z_start[MPI_Rank+1] - List_z_start[MPI_Rank] != local_nz )
      Aux_Error( ERROR_INFO, "local_nz (%d) != expectation (%d) !!\n",
                 local_nz, List_z_start[MPI_Rank+1] - List_z_start[MPI_Rank] );
#  endif // GAMER_DEBUG


   const int SSize[2]   = { FFT_Size[0], FFT_Size[1] };           // slab size in the x and y directions
   const int PSSize     = PS1*PS1;                                // patch slice size
   const int MemUnit    = amr->NPatchComma[0][1]*PS1;             // set arbitrarily
   const int AveNz      = FFT_Size[2]/MPI_NRank + ( (FFT_Size[2]%MPI_NRank == 0 ) ? 0 : 1 );    // average slab thickness
   const int Scale0     = amr->scale[0];

   int   Cr[3];                        // corner coordinates of each patch normalized to the base-level grid size
   int   BPos_z;                       // z coordinate of each patch slice in the simulation box
   int   SPos_z;                       // z coordinate of each patch slice in the slab
   long  SIdx;                         // 1D coordinate of each patch slice in the slab
   int   List_NSend_SIdx[MPI_NRank];   // number of patch slices sent to each rank
   int   List_NRecv_SIdx[MPI_NRank];   // number of patch slices received from each rank
   long *TempBuf_SIdx   [MPI_NRank];   // 1D slab coordinate of each patch slice sent to each rank
   real *TempBuf_Psi    [MPI_NRank];   // wavefunction of each patch slice sent to each rank
   real *TempBuf_Psi_Ptr = NULL;

   int TRank, TRank_Guess, MemSize[MPI_NRank], idx;


// 1. set memory allocation unit
   for (int r=0; r<MPI_NRank; r++)
   {
      MemSize        [r] = MemUnit;
      if ( Target == REAL  ) // REAL is the first one to use them
      {
         List_PID    [r] = (int* )malloc( MemSize[r]*sizeof(int)         );
         List_k      [r] = (int* )malloc( MemSize[r]*sizeof(int)         );
      }
      TempBuf_SIdx   [r] = (long*)malloc( MemSize[r]*sizeof(long)        );
      TempBuf_Psi    [r] = (real*)malloc( MemSize[r]*sizeof(real)*PSSize );
      List_NSend_SIdx[r] = 0;
   }


// 2. prepare the temporary send buffer and record lists
   const OptPotBC_t  PotBC_None        = BC_POT_NONE;
   const IntScheme_t IntScheme         = INT_NONE;
   const NSide_t     NSide_None        = NSIDE_00;
   const bool        IntPhase_No       = false;
   const bool        DE_Consistency_No = false;
   const real        MinDens_No        = -1.0;
   const real        MinPres_No        = -1.0;
   const real        MinTemp_No        = -1.0;
   const real        MinEntr_No        = -1.0;
   const int         GhostSize         = 0;
   const int         NPG               = 1;

   real (*Psi_prep)[PS1][PS1][PS1] = new real [8*NPG][PS1][PS1][PS1];

   for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)
   {
      if ( Target ==  REAL  )
         Prepare_PatchData( 0, Time, Psi_prep[0][0][0], NULL, GhostSize, NPG, &PID0, _REAL, _NONE,
                            IntScheme, INT_NONE, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
      else if ( Target == IMAG  )
         Prepare_PatchData( 0, Time, Psi_prep[0][0][0], NULL, GhostSize, NPG, &PID0, _IMAG, _NONE,
                            IntScheme, INT_NONE, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
      else if ( Target == DENS  )
         Prepare_PatchData( 0, Time, Psi_prep[0][0][0], NULL, GhostSize, NPG, &PID0, _DENS, _NONE,
                            IntScheme, INT_NONE, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


//    copy data to the send buffer
      for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
      {
         for (int d=0; d<3; d++)    Cr[d] = amr->patch[0][0][PID]->corner[d] / Scale0;

         for (int k=0; k<PS1; k++)
         {
            BPos_z      = Cr[2] + k;
            TRank_Guess = BPos_z / AveNz;
            TRank       = ZIndex2Rank_Psi( BPos_z, List_z_start, TRank_Guess );
            SPos_z      = BPos_z - List_z_start[TRank];
            SIdx        = ( (long)SPos_z*SSize[1] + Cr[1] )*SSize[0] + Cr[0];

#           ifdef GAMER_DEBUG
            if ( SPos_z < 0  ||  SPos_z >= List_z_start[TRank+1] - List_z_start[TRank] )
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SPos_z", SPos_z );
#           endif

//          allocate enough memory
            if ( List_NSend_SIdx[TRank] >= MemSize[TRank] )
            {
               MemSize     [TRank] += MemUnit;
               if ( Target == REAL  ) // REAL is the first one to use them
               {
                  List_PID [TRank]  = (int* )realloc( List_PID    [TRank], MemSize[TRank]*sizeof(int)         );
                  List_k   [TRank]  = (int* )realloc( List_k      [TRank], MemSize[TRank]*sizeof(int)         );
               }
               TempBuf_SIdx[TRank]  = (long*)realloc( TempBuf_SIdx[TRank], MemSize[TRank]*sizeof(long)        );
               TempBuf_Psi [TRank]  = (real*)realloc( TempBuf_Psi [TRank], MemSize[TRank]*sizeof(real)*PSSize );
            }

//          record list
            if ( Target == REAL  ) // REAL is the first one to use them
            {
               List_PID [TRank][ List_NSend_SIdx[TRank] ] = PID;
               List_k   [TRank][ List_NSend_SIdx[TRank] ] = k;
            }
            TempBuf_SIdx[TRank][ List_NSend_SIdx[TRank] ] = SIdx;

//          store data
            TempBuf_Psi_Ptr = TempBuf_Psi[TRank] + List_NSend_SIdx[TRank]*PSSize;

            idx = 0;
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
               TempBuf_Psi_Ptr[ idx ++ ] = Psi_prep[LocalID][k][j][i];


            List_NSend_SIdx[TRank] ++;
         } // for (int k=0; k<PS1; k++)
      } // for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
   } // for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)

   delete [] Psi_prep;


// 3. prepare the send buffer
   int   Send_Disp_Psi[MPI_NRank], Recv_Disp_Psi[MPI_NRank], Send_Disp_SIdx[MPI_NRank], Recv_Disp_SIdx[MPI_NRank];
   long *SendPtr_SIdx = NULL;
   real *SendPtr_Psi  = NULL;

// 3.1 broadcast the number of elements sending to different ranks
   MPI_Alltoall( List_NSend_SIdx, 1, MPI_INT, List_NRecv_SIdx, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   {
      List_NSend_Psi[r] = List_NSend_SIdx[r]*PSSize;
      List_NRecv_Psi[r] = List_NRecv_SIdx[r]*PSSize;
   }

// 3.2 calculate the displacement
   Send_Disp_SIdx[0] = 0;
   Recv_Disp_SIdx[0] = 0;
   Send_Disp_Psi [0] = 0;
   Recv_Disp_Psi [0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_SIdx[r] = Send_Disp_SIdx[r-1] + List_NSend_SIdx[r-1];
      Recv_Disp_SIdx[r] = Recv_Disp_SIdx[r-1] + List_NRecv_SIdx[r-1];
      Send_Disp_Psi [r] = Send_Disp_Psi [r-1] + List_NSend_Psi [r-1];
      Recv_Disp_Psi [r] = Recv_Disp_Psi [r-1] + List_NRecv_Psi [r-1];
   }

// check
#  ifdef GAMER_DEBUG
   const int NSend_Total  = Send_Disp_Psi[MPI_NRank-1] + List_NSend_Psi[MPI_NRank-1];
   const int NRecv_Total  = Recv_Disp_Psi[MPI_NRank-1] + List_NRecv_Psi[MPI_NRank-1];
   const int NSend_Expect = amr->NPatchComma[0][1]*CUBE(PS1);
   const int NRecv_Expect = NX0_TOT[0]*NX0_TOT[1]*NRecvSlice;

   if ( NSend_Total != NSend_Expect )  Aux_Error( ERROR_INFO, "NSend_Total = %d != expected value = %d !!\n",
                                                  NSend_Total, NSend_Expect );

   if ( NRecv_Total != NRecv_Expect )  Aux_Error( ERROR_INFO, "NRecv_Total = %d != expected value = %d !!\n",
                                                  NRecv_Total, NRecv_Expect );
#  endif

// 3.3 prepare the send buffer of SIdx
   SendPtr_SIdx = SendBuf_SIdx;
   for (int r=0; r<MPI_NRank; r++)
   {
      memcpy( SendPtr_SIdx, TempBuf_SIdx[r], List_NSend_SIdx[r]*sizeof(long) );
      SendPtr_SIdx += List_NSend_SIdx[r];
   }

// 3.4 prepare the send buffer of data
   SendPtr_Psi  = SendBuf_Psi;
   for (int r=0; r<MPI_NRank; r++)
   {
      memcpy( SendPtr_Psi, TempBuf_Psi[r], List_NSend_Psi[r]*sizeof(real) );
      SendPtr_Psi += List_NSend_Psi[r];
   }


// 4. exchange data by MPI
   MPI_Alltoallv( SendBuf_SIdx, List_NSend_SIdx, Send_Disp_SIdx, MPI_LONG,
                  RecvBuf_SIdx, List_NRecv_SIdx, Recv_Disp_SIdx, MPI_LONG,   MPI_COMM_WORLD );

#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_Psi ,  List_NSend_Psi ,  Send_Disp_Psi ,  MPI_DOUBLE,
                  RecvBuf_Psi ,  List_NRecv_Psi ,  Recv_Disp_Psi ,  MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_Psi ,  List_NSend_Psi ,  Send_Disp_Psi ,  MPI_FLOAT,
                  RecvBuf_Psi ,  List_NRecv_Psi ,  Recv_Disp_Psi ,  MPI_FLOAT,  MPI_COMM_WORLD );
#  endif


// 5. store the received data to array "Psi" for FFTW
   const long NPSlice = (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/PSSize;  // total number of received patch slices
   long  dSIdx, Counter = 0;
   real *Psi_Ptr = NULL;

   for (long t=0; t<NPSlice; t++)
   {
      SIdx     = RecvBuf_SIdx[t];
      Psi_Ptr = Psi + SIdx;

      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         dSIdx           = j*SSize[0] + i;
         Psi_Ptr[dSIdx] = RecvBuf_Psi[ Counter ++ ];
      }
   }


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( TempBuf_SIdx[r] );
      free( TempBuf_Psi [r] );
   }

} // FUNCTION : Patch2Slab_Psi



//-------------------------------------------------------------------------------------------------------
// Function    :  ZIndex2Rank_Psi
// Description :  Return the MPI rank which the input z coordinates belongs to in the FFTW slab decomposition
//
// Note        :  1. "List_z_start[r] <= IndexZ < List_z_start[r+1]" belongs to rank r
//                2. List_z_start[MPI_NRank] can be set to any value >= FFT_Size[2]
//
// Parameter   :  IndexZ       : Input z coordinate
//                List_z_start : Starting z coordinate of each rank in the FFTW slab decomposition
//                TRank_Guess  : First guess of the targeting MPI rank
//
// Return      :  MPI rank
//-------------------------------------------------------------------------------------------------------
int ZIndex2Rank_Psi( const int IndexZ, const int *List_z_start, const int TRank_Guess )
{

// check
#  ifdef GAMER_DEBUG
   const int FFT_SizeZ = NX0_TOT[2];

   if ( IndexZ < 0  ||  IndexZ >= NX0_TOT[2] )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IndexZ", IndexZ );

   if ( List_z_start[MPI_NRank] < FFT_SizeZ )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "List_z_start[MPI_NRank]", List_z_start[MPI_NRank] );

   if ( TRank_Guess < 0  ||  TRank_Guess >= MPI_NRank )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TRank_Guess", TRank_Guess );
#  endif


   int TRank = TRank_Guess;   // have a first guess to improve the performance

   while ( true )
   {
#     ifdef GAMER_DEBUG
      if ( TRank < 0  ||  TRank >= MPI_NRank )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TRank", TRank );
#     endif

      if ( IndexZ < List_z_start[TRank] )    TRank --;
      else
      {
         if ( IndexZ < List_z_start[TRank+1] )  return TRank;
         else                                   TRank ++;
      }
   }

} // FUNCTION : ZIndex2Rank_Psi



//-------------------------------------------------------------------------------------------------------
// Function    :  Slab2Patch_Psi
// Description :  Slab domain decomposition --> patch-based data (for wavefunction)
//
// Parameter   :  Psi        : Array storing wavefunction
//                SendBuf    : Sending MPI buffer of wavefunction
//                RecvBuf    : Receiving MPI buffer of wavefunction
//                SaveSg     : Sandglass to store the updated data
//                List_SIdx  : 1D coordinate in slab
//                List_PID   : PID of each patch slice sent to each rank
//                List_k     : Local z coordinate of each patch slice sent to each rank
//                List_NSend : Size of wavefunction data sent to each rank
//                List_NRecv : Size of wavefunction data received from each rank
//                local_nz   : Slab thickness of this MPI rank
//                FFT_Size   : Size of the FFT operation
//                NSendSlice : Total number of z slices need to be sent to other ranks
//                Target     : Target variable to be prepared (REAL or IMAG or DENS)
//-------------------------------------------------------------------------------------------------------
void Slab2Patch_Psi( const real *Psi, real *SendBuf, real *RecvBuf, const int SaveSg, const long *List_SIdx,
                 int **List_PID, int **List_k, int *List_NSend, int *List_NRecv, const int local_nz, const int FFT_Size[],
                 const int NSendSlice, const int Target )
{

// 1. store the evaluated wavefunction to the send buffer
   const int   SSize[2]   = { FFT_Size[0], FFT_Size[1] };                     // slab size in the x and y directions
   const int   PSSize     = PS1*PS1;                                          // patch slice size
   const long  NPSlice    = (long)NX0_TOT[0]*NX0_TOT[1]*NSendSlice/PSSize;    // total number of patch slices to be sent
   const real *Psi_Ptr  = NULL;

   long SIdx, dSIdx, Counter = 0;

   for (long t=0; t<NPSlice; t++)
   {
      SIdx     = List_SIdx[t];
      Psi_Ptr = Psi + SIdx;

      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         dSIdx                 = j*SSize[0] + i;
         SendBuf[ Counter ++ ] = Psi_Ptr[dSIdx];
      }
   }


// 2. calculate the displacement and exchange data by MPI
   int Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];

   Send_Disp[0] = 0;
   Recv_Disp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + List_NSend[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + List_NRecv[r-1];
   }

#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf, List_NSend, Send_Disp, MPI_DOUBLE,
                  RecvBuf, List_NRecv, Recv_Disp, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf, List_NSend, Send_Disp, MPI_FLOAT,
                  RecvBuf, List_NRecv, Recv_Disp, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif


// 3. store the received data to different patch objects
   int   PID, k, NRecvSlice;
   real *RecvPtr = RecvBuf;

   for (int r=0; r<MPI_NRank; r++)
   {
      NRecvSlice = List_NRecv[r]/PSSize;

      for (int t=0; t<NRecvSlice; t++)
      {
         PID = List_PID[r][t];
         k   = List_k  [r][t];

         if ( Target == REAL  )
            memcpy( amr->patch[SaveSg][0][PID]->fluid[REAL][k], RecvPtr, PSSize*sizeof(real) );
         else if ( Target == IMAG  )
            memcpy( amr->patch[SaveSg][0][PID]->fluid[IMAG][k], RecvPtr, PSSize*sizeof(real) );
         else if ( Target == DENS  )
            memcpy( amr->patch[SaveSg][0][PID]->fluid[DENS][k], RecvPtr, PSSize*sizeof(real) );

         RecvPtr += PSSize;
      }
   }


// free memory
   if ( Target == DENS  ) // DENS is the last one to use them
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         free( List_PID[r] );
         free( List_k  [r] );
      }
   }

} // FUNCTION : Slab2Patch_Psi



//-------------------------------------------------------------------------------------------------------
// Function    :  Psi_Advance_FFT
// Description :  Evolve the wavefunction by FFT (the pseudo-spectral method)
//
// Note        :
//
// Parameter   :  PsiR     : Array storing the real part of wavefunction (input and output)
//                PsiI     : Array storing the imag part of wavefunction (input and output)
//                PsiD     : Array storing the density (output)
//                j_start  : Starting j index
//                dj       : Size of array in the j (y) direction after the forward FFT
//                PsiK_Size: Size of the array "PsiK"
//                dt       : Time interval to advance solution
//-------------------------------------------------------------------------------------------------------
void Psi_Advance_FFT( real *PsiR, real *PsiI, real *PsiD, const int j_start, const int dj, const int PsiK_Size, const real dt )
{

   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const real dh       = amr->dh[0];
   real PsiKR, PsiKI, DtKK_2Eta;
   fftw_complex *PsiK;
   PsiK = (fftw_complex*) malloc( PsiK_Size * sizeof(fftw_complex));

   for (int t=0; t<PsiK_Size; t++)
   {
      PsiK[t].re = PsiR[t];
      PsiK[t].im = PsiI[t];
   }


// forward FFT
#  ifdef SERIAL
   fftwnd_one( FFTW_Plan_Psi, PsiK, NULL );
#  else
   fftwnd_mpi( FFTW_Plan_Psi, 1, PsiK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// set up the dimensional wave number
   real kx[Nx], ky[Ny], kz[Nz];

   for (int i=0; i<Nx; i++) { kx[i] = ( i <= Nx/2 ) ? 2.0*M_PI/(Nx*dh)*i : 2.0*M_PI/(Nx*dh)*(i-Nx);}
   for (int j=0; j<Ny; j++) { ky[j] = ( j <= Ny/2 ) ? 2.0*M_PI/(Ny*dh)*j : 2.0*M_PI/(Ny*dh)*(j-Ny);}
   for (int k=0; k<Nz; k++) { kz[k] = ( k <= Nz/2 ) ? 2.0*M_PI/(Nz*dh)*k : 2.0*M_PI/(Nz*dh)*(k-Nz);}


// multiply wavefunction and exp(-i*dt*k^2/(2*ELBDM_ETA)) in the k space
   long ID;
#  ifdef SERIAL // serial mode

   for (int k=0; k<Nz; k++)
   {
      for (int j=0; j<Ny; j++)
      for (int i=0; i<Nx; i++)
      {
         ID = ((long)k*Ny + j)*Nx + i;

#  else // parallel mode
   int j;

   for (int jj=0; jj<dj; jj++)
   {
      j = j_start + jj;

      for (int k=0; k<Nz; k++)
      for (int i=0; i<Nx; i++)
      {
         ID = ((long)jj*Nz + k)*Nx + i;

#  endif // #ifdef SERIAL ... else ...

         PsiKR = PsiK[ID].re;
         PsiKI = PsiK[ID].im;
         DtKK_2Eta = (real)0.5*dt*( kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k] )/ELBDM_ETA;

         PsiK[ID].re =  PsiKR * COS(DtKK_2Eta) + PsiKI * SIN(DtKK_2Eta);
         PsiK[ID].im =  PsiKI * COS(DtKK_2Eta) - PsiKR * SIN(DtKK_2Eta);

      } // i,j,k
   } // i,j,k


// backward FFT
#  ifdef SERIAL
   fftwnd_one( FFTW_Plan_Psi_Inv, PsiK, NULL );
#  else
   fftwnd_mpi( FFTW_Plan_Psi_Inv, 1, PsiK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif

// normalization
   const real norm = 1.0 / ( (real)Nx*Ny*Nz );

   for (int t=0; t<PsiK_Size; t++)
   {
      PsiR[t] = PsiK[t].re * norm;
      PsiI[t] = PsiK[t].im * norm;
      PsiD[t] = SQR( PsiR[t] ) + SQR( PsiI[t] );
   }

   free( PsiK );

} // FUNCTION : Psi_Advance_FFT



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonSolver_FFT
// Description :  Evolve base-level wavefunction by FFT
//
// Note        :  1. Work with the option ELBDM_BASE_SPECTRAL
//                2. Invoked by Flu_AdvanceDt()
//
// Parameter   :  dt       : Time interval to advance solution
//                Time     : Physical time for preparing the wavefunction
//                SaveSg   : Sandglass to store the updated data
//-------------------------------------------------------------------------------------------------------
void CPU_ELBDMSolver_FFT( const real dt, const double Time, const int SaveSg )
{
   Aux_Message( stdout, "ELBDM_BASE_SPECTRAL is working.\n" );
// determine the FFT size
   int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };


// get the array indices using by FFTW
   int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

#  ifdef SERIAL
   local_nz                      = FFT_Size[2];
   local_z_start                 = 0;
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = FFT_Size[0]*FFT_Size[1]*FFT_Size[2];
#  else
   fftwnd_mpi_local_sizes( FFTW_Plan_Psi, &local_nz, &local_z_start, &local_ny_after_transpose,
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


// allocate memory
   const int NRecvSlice = MIN( List_z_start[MPI_Rank]+local_nz, NX0_TOT[2] ) - MIN( List_z_start[MPI_Rank], NX0_TOT[2] );

   real *PsiR         = new real [ total_local_size ];                           // array storing the real part of wavefunction
   real *PsiI         = new real [ total_local_size ];                           // array storing the imag part of wavefunction
   real *PsiD         = new real [ total_local_size ];                           // array storing the density
   real *SendBuf      = new real [ (long)amr->NPatchComma[0][1]*CUBE(PS1) ];     // MPI send buffer
   real *RecvBuf      = new real [ (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice ];     // MPI recv buffer
   long *SendBuf_SIdx = new long [ amr->NPatchComma[0][1]*PS1 ];                 // MPI send buffer for 1D coordinate in slab
   long *RecvBuf_SIdx = new long [ NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/SQR(PS1) ];  // MPI recv buffer for 1D coordinate in slab

   int  *List_PID    [MPI_NRank];   // PID of each patch slice sent to each rank
   int  *List_k      [MPI_NRank];   // local z coordinate of each patch slice sent to each rank
   int   List_NSend  [MPI_NRank];   // size of data sent to each rank
   int   List_NRecv  [MPI_NRank];   // size of data received from each rank


// rearrange data from patch to slab
   Patch2Slab_Psi( PsiR, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
                   local_nz, FFT_Size, NRecvSlice, Time, REAL );
   Patch2Slab_Psi( PsiI, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
                   local_nz, FFT_Size, NRecvSlice, Time, IMAG );
   Patch2Slab_Psi( PsiD, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
                   local_nz, FFT_Size, NRecvSlice, Time, DENS );


// evolve wavefunction by FFT
   Psi_Advance_FFT( PsiR, PsiI, PsiD, local_y_start_after_transpose, local_ny_after_transpose, total_local_size , dt );


// rearrange data from slab back to patch
   Slab2Patch_Psi( PsiR, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID, List_k, List_NRecv, List_NSend,
                   local_nz, FFT_Size, NRecvSlice, REAL );
   Slab2Patch_Psi( PsiI, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID, List_k, List_NRecv, List_NSend,
                   local_nz, FFT_Size, NRecvSlice, IMAG );
   Slab2Patch_Psi( PsiD, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID, List_k, List_NRecv, List_NSend,
                   local_nz, FFT_Size, NRecvSlice, DENS );


   delete [] PsiR;
   delete [] PsiI;
   delete [] PsiD;
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;

} // FUNCTION : CPU_ELBDMSolver_FFT



#endif // #if ( MODEL == ELBDM )
