#include "GAMER.h"

#ifdef GRAVITY



static void FFT_Periodic( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const int RhoK_Size );
static void FFT_Isolated( real *RhoK, const real *gFuncK, const real Poi_Coeff, const int RhoK_Size );
static int ZIndex2Rank( const int IndexZ, const int *List_z_start, const int TRank_Guess );

#ifdef SERIAL
extern rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv;
#else
extern rfftwnd_mpi_plan FFTW_Plan, FFTW_Plan_Inv;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Patch2Slab
// Description :  Patch-based data --> slab domain decomposition (for density)
//
// Parameter   :  RhoK           : In-place FFT array
//                SendBuf_Rho    : Sending MPI buffer of density
//                RecvBuf_Rho    : Receiving MPI buffer of density
//                SendBuf_SIdx   : Sending MPI buffer of 1D coordinate in slab
//                RecvBuf_SIdx   : Receiving MPI buffer of 1D coordinate in slab
//                List_PID       : PID of each patch slice sent to each rank
//                List_k         : Local z coordinate of each patch slice sent to each rank
//                List_NSend_Rho : Size of density data sent to each rank
//                List_NRecv_Rho : Size of density data received from each rank
//                List_z_start   : Starting z coordinate of each rank in the FFTW slab decomposition
//                local_nz       : Slab thickness of this MPI rank
//                FFT_Size       : Size of the FFT operation including the zero-padding regions
//                NRecvSlice     : Total number of z slices received from other ranks (could be zero in the isolated BC)
//                PrepTime       : Physical time for preparing the density field
//-------------------------------------------------------------------------------------------------------
void Patch2Slab( real *RhoK, real *SendBuf_Rho, real *RecvBuf_Rho, long *SendBuf_SIdx, long *RecvBuf_SIdx,
                 int **List_PID, int **List_k, int *List_NSend_Rho, int *List_NRecv_Rho,
                 const int *List_z_start, const int local_nz, const int FFT_Size[], const int NRecvSlice,
                 const double PrepTime )
{

// check
#  ifdef GAMER_DEBUG
   if ( List_z_start[MPI_Rank+1] - List_z_start[MPI_Rank] != local_nz )
      Aux_Error( ERROR_INFO, "local_nz (%d) != expectation (%d) !!\n",
                 local_nz, List_z_start[MPI_Rank+1] - List_z_start[MPI_Rank] );
#  endif // GAMER_DEBUG


   const int SSize[2]   = { 2*(FFT_Size[0]/2+1), FFT_Size[1] };   // padded slab size in the x and y directions
   const int PSSize     = PS1*PS1;                                // patch slice size
// const int MemUnit    = amr->NPatchComma[0][1]*PS1/MPI_NRank;   // set arbitrarily
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
   real *TempBuf_Rho    [MPI_NRank];   // density of each patch slice sent to each rank
   real *TempBuf_Rho_Ptr = NULL;

   int TRank, TRank_Guess, MemSize[MPI_NRank], idx;


// 1. set memory allocation unit
   for (int r=0; r<MPI_NRank; r++)
   {
      MemSize        [r] = MemUnit;
      List_PID       [r] = (int* )malloc( MemSize[r]*sizeof(int)         );
      List_k         [r] = (int* )malloc( MemSize[r]*sizeof(int)         );
      TempBuf_SIdx   [r] = (long*)malloc( MemSize[r]*sizeof(long)        );
      TempBuf_Rho    [r] = (real*)malloc( MemSize[r]*sizeof(real)*PSSize );
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
   const int         GhostSize         = 0;
   const int         NPG               = 1;

   real (*Dens)[PS1][PS1][PS1] = new real [8*NPG][PS1][PS1][PS1];

   for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)
   {
//    even with NSIDE_00 and GhostSize=0, we still need OPT__BC_FLU to determine whether periodic BC is adopted
//    also note that we do not check minimum density here since no ghost zones are required
      Prepare_PatchData( 0, PrepTime, Dens[0][0][0], GhostSize, NPG, &PID0, _TOTAL_DENS,
                         IntScheme, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                         MinDens_No, MinPres_No, DE_Consistency_No );

      for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
      {
         for (int d=0; d<3; d++)    Cr[d] = amr->patch[0][0][PID]->corner[d] / Scale0;

         for (int k=0; k<PS1; k++)
         {
            BPos_z      = Cr[2] + k;
            TRank_Guess = BPos_z / AveNz;
            TRank       = ZIndex2Rank( BPos_z, List_z_start, TRank_Guess );
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
               List_PID    [TRank]  = (int* )realloc( List_PID    [TRank], MemSize[TRank]*sizeof(int)         );
               List_k      [TRank]  = (int* )realloc( List_k      [TRank], MemSize[TRank]*sizeof(int)         );
               TempBuf_SIdx[TRank]  = (long*)realloc( TempBuf_SIdx[TRank], MemSize[TRank]*sizeof(long)        );
               TempBuf_Rho [TRank]  = (real*)realloc( TempBuf_Rho [TRank], MemSize[TRank]*sizeof(real)*PSSize );
            }

//          record list
            List_PID    [TRank][ List_NSend_SIdx[TRank] ] = PID;
            List_k      [TRank][ List_NSend_SIdx[TRank] ] = k;
            TempBuf_SIdx[TRank][ List_NSend_SIdx[TRank] ] = SIdx;

//          store density
            TempBuf_Rho_Ptr = TempBuf_Rho[TRank] + List_NSend_SIdx[TRank]*PSSize;

            idx = 0;
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
               TempBuf_Rho_Ptr[ idx ++ ] = Dens[LocalID][k][j][i];
//             TempBuf_Rho_Ptr[ idx ++ ] = amr->patch[Sg][0][PID]->fluid[DENS][k][j][i];

//          subtract the background density (which is assumed to be UNITY) for the isolated BC in the comoving frame
//          --> to be consistent with the comoving-frame Poisson eq.
#           ifdef COMOVING
            if ( OPT__BC_POT == BC_POT_ISOLATED )
            {
               for (int t=0; t<PSSize; t++)  TempBuf_Rho_Ptr[t] -= (real)1.0;
            }
#           endif

            List_NSend_SIdx[TRank] ++;
         } // for (int k=0; k<PS1; k++)
      } // for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
   } // for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)

   delete [] Dens;


// 3. prepare the send buffer
   int   Send_Disp_Rho[MPI_NRank], Recv_Disp_Rho[MPI_NRank], Send_Disp_SIdx[MPI_NRank], Recv_Disp_SIdx[MPI_NRank];
   long *SendPtr_SIdx = NULL;
   real *SendPtr_Rho  = NULL;

// 3.1 broadcast the number of elements sending to different ranks
   MPI_Alltoall( List_NSend_SIdx, 1, MPI_INT, List_NRecv_SIdx, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   {
      List_NSend_Rho[r] = List_NSend_SIdx[r]*PSSize;
      List_NRecv_Rho[r] = List_NRecv_SIdx[r]*PSSize;
   }

// 3.2 calculate the displacement
   Send_Disp_SIdx[0] = 0;
   Recv_Disp_SIdx[0] = 0;
   Send_Disp_Rho [0] = 0;
   Recv_Disp_Rho [0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_SIdx[r] = Send_Disp_SIdx[r-1] + List_NSend_SIdx[r-1];
      Recv_Disp_SIdx[r] = Recv_Disp_SIdx[r-1] + List_NRecv_SIdx[r-1];
      Send_Disp_Rho [r] = Send_Disp_Rho [r-1] + List_NSend_Rho [r-1];
      Recv_Disp_Rho [r] = Recv_Disp_Rho [r-1] + List_NRecv_Rho [r-1];
   }

// check
#  ifdef GAMER_DEBUG
   const int NSend_Total  = Send_Disp_Rho[MPI_NRank-1] + List_NSend_Rho[MPI_NRank-1];
   const int NRecv_Total  = Recv_Disp_Rho[MPI_NRank-1] + List_NRecv_Rho[MPI_NRank-1];
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

// 3.4 prepare the send buffer of density
   SendPtr_Rho = SendBuf_Rho;
   for (int r=0; r<MPI_NRank; r++)
   {
      memcpy( SendPtr_Rho, TempBuf_Rho[r], List_NSend_Rho[r]*sizeof(real) );
      SendPtr_Rho += List_NSend_Rho[r];
   }


// 4. exchange data by MPI
   MPI_Alltoallv( SendBuf_SIdx, List_NSend_SIdx, Send_Disp_SIdx, MPI_LONG,
                  RecvBuf_SIdx, List_NRecv_SIdx, Recv_Disp_SIdx, MPI_LONG,   MPI_COMM_WORLD );

#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_Rho,  List_NSend_Rho,  Send_Disp_Rho,  MPI_DOUBLE,
                  RecvBuf_Rho,  List_NRecv_Rho,  Recv_Disp_Rho,  MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_Rho,  List_NSend_Rho,  Send_Disp_Rho,  MPI_FLOAT,
                  RecvBuf_Rho,  List_NRecv_Rho,  Recv_Disp_Rho,  MPI_FLOAT,  MPI_COMM_WORLD );
#  endif


// 5. store the received density to the padded array "RhoK" for FFTW
   const long NPSlice = (long)NX0_TOT[0]*NX0_TOT[1]*NRecvSlice/PSSize;  // total number of received patch slices
   long  dSIdx, Counter = 0;
   real *RhoK_Ptr = NULL;

   for (long t=0; t<NPSlice; t++)
   {
      SIdx     = RecvBuf_SIdx[t];
      RhoK_Ptr = RhoK + SIdx;

      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         dSIdx           = j*SSize[0] + i;
         RhoK_Ptr[dSIdx] = RecvBuf_Rho[ Counter ++ ];
      }
   }


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( TempBuf_SIdx[r] );
      free( TempBuf_Rho [r] );
   }

} // FUNCTION : Patch2Slab



//-------------------------------------------------------------------------------------------------------
// Function    :  ZIndex2Rank
// Description :  Return the MPI rank which the input z coordinates belongs to in the FFTW slab decomposition
//
// Note        :  1. "List_z_start[r] <= IndexZ < List_z_start[r+1]" belongs to rank r
//                2. List_z_start[MPI_NRank] can be set to any value >= FFT_Size[2]
//
// Parameter   :  IndexZ         : Input z coordinate
//                List_z_start   : Starting z coordinate of each rank in the FFTW slab decomposition
//                TRank_Guess    : First guess of the targeting MPI rank
//
// Return      :  MPI rank
//-------------------------------------------------------------------------------------------------------
int ZIndex2Rank( const int IndexZ, const int *List_z_start, const int TRank_Guess )
{

// check
#  ifdef GAMER_DEBUG
   const int FFT_SizeZ = ( OPT__BC_POT == BC_POT_ISOLATED ) ? 2*NX0_TOT[2] : NX0_TOT[2];

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

} // FUNCTION : ZIndex2Rank



//-------------------------------------------------------------------------------------------------------
// Function    :  Slab2Patch
// Description :  Slab domain decomposition --> patch-based data (for potential)
//
// Parameter   :  RhoK        : In-place FFT array
//                SendBuf     : Sending MPI buffer of potential
//                RecvBuf     : Receiving MPI buffer of potential
//                SaveSg      : Sandglass to store the updated data
//                List_SIdx   : 1D coordinate in slab
//                List_PID    : PID of each patch slice sent to each rank
//                List_k      : Local z coordinate of each patch slice sent to each rank
//                List_NSend  : Size of potential data sent to each rank
//                List_NRecv  : Size of potential data received from each rank
//                local_nz    : Slab thickness of this MPI rank
//                FFT_Size    : Size of the FFT operation including the zero-padding regions
//                NSendSlice  : Total number of z slices need to be sent to other ranks (could be zero in the isolated BC)
//-------------------------------------------------------------------------------------------------------
void Slab2Patch( const real *RhoK, real *SendBuf, real *RecvBuf, const int SaveSg, const long *List_SIdx,
                 int **List_PID, int **List_k, int *List_NSend, int *List_NRecv, const int local_nz, const int FFT_Size[],
                 const int NSendSlice )
{

// 1. store the evaluated potential to the send buffer
   const int   SSize[2]   = { 2*(FFT_Size[0]/2+1), FFT_Size[1] };             // padded slab size in the x and y directions
   const int   PSSize     = PS1*PS1;                                          // patch slice size
   const long  NPSlice    = (long)NX0_TOT[0]*NX0_TOT[1]*NSendSlice/PSSize;    // total number of patch slices to be sent
   const real *RhoK_Ptr   = NULL;

   long SIdx, dSIdx, Counter = 0;

   for (long t=0; t<NPSlice; t++)
   {
      SIdx     = List_SIdx[t];
      RhoK_Ptr = RhoK + SIdx;

      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         dSIdx                 = j*SSize[0] + i;
         SendBuf[ Counter ++ ] = RhoK_Ptr[dSIdx];
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


// 3. store the received potential data to different patch objects
   int   PID, k, NRecvSlice;
   real *RecvPtr = RecvBuf;

   for (int r=0; r<MPI_NRank; r++)
   {
      NRecvSlice = List_NRecv[r]/PSSize;

      for (int t=0; t<NRecvSlice; t++)
      {
         PID = List_PID[r][t];
         k   = List_k  [r][t];

         memcpy( amr->patch[SaveSg][0][PID]->pot[k], RecvPtr, PSSize*sizeof(real) );

         RecvPtr += PSSize;
      }
   }


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( List_PID[r] );
      free( List_k  [r] );
   }

} // FUNCTION : Slab2Patch



//-------------------------------------------------------------------------------------------------------
// Function    :  FFT_Periodic
// Description :  Evaluate the gravitational potential by FFT for the periodic BC
//
// Note        :  Effect from the homogenerous background density (DC) will be ignored by setting the k=0 mode
//                equal to zero
//
// Parameter   :  RhoK        : Array storing the input density and output potential
//                Poi_Coeff   : Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)
//                j_start     : Starting j index
//                dj          : Size of array in the j (y) direction after the forward FFT
//                RhoK_Size   : Size of the array "RhoK"
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
   rfftwnd_one_real_to_complex( FFTW_Plan, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
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
   rfftwnd_one_complex_to_real( FFTW_Plan_Inv, cdata, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Inv, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// normalization
   const real norm = dh*dh / ( (real)Nx*Ny*Nz );

   for (int t=0; t<RhoK_Size; t++)  RhoK[t] *= norm;

} // FUNCTION : FFT_Periodic



//-------------------------------------------------------------------------------------------------------
// Function    :  FFT_Isolated
// Description :  Evaluate the gravitational potential by FFT for the isolated BC
//
// Note        :  1. The Green's function in the k space has been initialized by the function "Init_GreenFuncK"
//                2. 4*PI*NEWTON_G and FFT normalization coefficient has been included in gFuncK
//                   --> The only coefficient that hasn't been taken into account is the scale factor in the comoving frame
//
// Parameter   :  RhoK        : Array storing the input density and output potential
//                Poi_Coeff   : Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)
//                RhoK_Size   : Size of the array "RhoK"
//-------------------------------------------------------------------------------------------------------
void FFT_Isolated( real *RhoK, const real *gFuncK, const real Poi_Coeff, const int RhoK_Size )
{

   fftw_complex *RhoK_cplx   = (fftw_complex *)RhoK;
   fftw_complex *gFuncK_cplx = (fftw_complex *)gFuncK;
   fftw_complex  Temp_cplx;


// forward FFT
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
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
   rfftwnd_one_complex_to_real( FFTW_Plan_Inv, RhoK_cplx, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Inv, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
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
// Parameter   :  Poi_Coeff   : Coefficient in front of the RHS in the Poisson eq.
//                SaveSg      : Sandglass to store the updated data
//                PrepTime    : Physical time for preparing the density field
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
   rfftwnd_mpi_local_sizes( FFTW_Plan, &local_nz, &local_z_start, &local_ny_after_transpose,
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


// initialize RhoK as zeros for the isolated BC where the zero-padding method is adopted
   if ( OPT__BC_POT == BC_POT_ISOLATED )
      for (int t=0; t<total_local_size; t++)    RhoK[t] = (real)0.0;


// rearrange data from patch to slab
   Patch2Slab( RhoK, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv, List_z_start,
               local_nz, FFT_Size, NRecvSlice, PrepTime );


// evaluate potential by FFT
   if      ( OPT__BC_POT == BC_POT_PERIODIC )
      FFT_Periodic( RhoK, Poi_Coeff, local_y_start_after_transpose, local_ny_after_transpose, total_local_size );

   else if ( OPT__BC_POT == BC_POT_ISOLATED )
      FFT_Isolated( RhoK, GreenFuncK, Poi_Coeff, total_local_size );

   else
      Aux_Error( ERROR_INFO, "unsupported paramter %s = %d !!\n", "OPT__BC_POT", OPT__BC_POT );


// rearrange data from slab back to patch
   Slab2Patch( RhoK, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID, List_k, List_NRecv, List_NSend,
               local_nz, FFT_Size, NRecvSlice );


   delete [] RhoK;
   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;

} // FUNCTION : CPU_PoissonSolver_FFT



#endif // #ifdef GRAVITY
