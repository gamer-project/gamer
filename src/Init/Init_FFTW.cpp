#include "GAMER.h"

#ifdef SUPPORT_FFTW



static int ZIndex2Rank( const int IndexZ, const int *List_z_start, const int TRank_Guess );

#ifdef SERIAL
rfftwnd_plan     FFTW_Plan_PS;                        // PS  : plan for calculating the power spectrum

#ifdef GRAVITY
rfftwnd_plan     FFTW_Plan_Poi, FFTW_Plan_Poi_Inv;    // Poi : plan for the self-gravity Poisson solver
#endif

#else
rfftwnd_mpi_plan FFTW_Plan_PS;

#ifdef GRAVITY
rfftwnd_mpi_plan FFTW_Plan_Poi, FFTW_Plan_Poi_Inv;
#endif

#endif

#ifdef GRAVITY
extern real (*Poi_AddExtraMassForGravity_Ptr)( const double x, const double y, const double z, const double Time,
                                               const int lv, double AuxArray[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW
// Description :  Create the FFTW plans
//-------------------------------------------------------------------------------------------------------
void Init_FFTW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );

// create plans for calculating the power spectrum
#  ifdef SERIAL
   FFTW_Plan_PS = rfftw3d_create_plan( NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], FFTW_REAL_TO_COMPLEX,
                                       FFTW_ESTIMATE | FFTW_IN_PLACE );
#  else
   FFTW_Plan_PS = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, NX0_TOT[2], NX0_TOT[1], NX0_TOT[0],
                                           FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );
#  endif


#  ifdef GRAVITY
// determine the FFT size for the self-gravity solver
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
   FFTW_Plan_Poi     = rfftw3d_create_plan( FFT_Size[2], FFT_Size[1], FFT_Size[0], FFTW_REAL_TO_COMPLEX,
                                            FFTW_ESTIMATE | FFTW_IN_PLACE );

   FFTW_Plan_Poi_Inv = rfftw3d_create_plan( FFT_Size[2], FFT_Size[1], FFT_Size[0], FFTW_COMPLEX_TO_REAL,
                                            FFTW_ESTIMATE | FFTW_IN_PLACE );

#  else

   FFTW_Plan_Poi     = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, FFT_Size[2], FFT_Size[1], FFT_Size[0],
                                                FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );

   FFTW_Plan_Poi_Inv = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, FFT_Size[2], FFT_Size[1], FFT_Size[0],
                                                FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE );
#  endif
#  endif // #ifdef GRAVITY


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
   rfftwnd_destroy_plan    ( FFTW_Plan_PS      );

#  ifdef GRAVITY
   rfftwnd_destroy_plan    ( FFTW_Plan_Poi     );
   rfftwnd_destroy_plan    ( FFTW_Plan_Poi_Inv );
#  endif

#  else // #ifdef SERIAL
   rfftwnd_mpi_destroy_plan( FFTW_Plan_PS      );

#  ifdef GRAVITY
   rfftwnd_mpi_destroy_plan( FFTW_Plan_Poi     );
   rfftwnd_mpi_destroy_plan( FFTW_Plan_Poi_Inv );
#  endif

#  endif // #ifdef SERIAL ... else ...

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_FFTW



//-------------------------------------------------------------------------------------------------------
// Function    :  Patch2Slab_Rho
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
//                AddExtraMass   : Adding an extra density field for computing gravitational potential only
//-------------------------------------------------------------------------------------------------------
void Patch2Slab_Rho( real *RhoK, real *SendBuf_Rho, real *RecvBuf_Rho, long *SendBuf_SIdx, long *RecvBuf_SIdx,
                     int **List_PID, int **List_k, int *List_NSend_Rho, int *List_NRecv_Rho,
                     const int *List_z_start, const int local_nz, const int FFT_Size[], const int NRecvSlice,
                     const double PrepTime, const bool AddExtraMass )
{

#  ifdef GRAVITY
// check
   if ( AddExtraMass  &&  Poi_AddExtraMassForGravity_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Poi_AddExtraMassForGravity_Ptr == NULL for AddExtraMass !!\n" );
#  endif // GRAVITY

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
   const real        MinTemp_No        = -1.0;
   const real        MinEntr_No        = -1.0;
   const int         GhostSize         = 0;
   const int         NPG               = 1;

   real (*Dens)[PS1][PS1][PS1] = new real [8*NPG][PS1][PS1][PS1];

   for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)
   {
//    even with NSIDE_00 and GhostSize=0, we still need OPT__BC_FLU to determine whether periodic BC is adopted
//    for depositing particle mass onto grids.
//    also note that we do not check minimum density here since no ghost zones are required
      Prepare_PatchData( 0, PrepTime, Dens[0][0][0], NULL, GhostSize, NPG, &PID0, _TOTAL_DENS, _NONE,
                         IntScheme, INT_NONE, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


#     ifdef GRAVITY
//    add extra mass source for gravity if required
      if ( AddExtraMass )
      {
         const double dh = amr->dh[0];

         for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
         {
            const double x0 = amr->patch[0][0][PID]->EdgeL[0] + 0.5*dh;
            const double y0 = amr->patch[0][0][PID]->EdgeL[1] + 0.5*dh;
            const double z0 = amr->patch[0][0][PID]->EdgeL[2] + 0.5*dh;

            double x, y, z;

            for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
            for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
            for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;
               Dens[LocalID][k][j][i] += Poi_AddExtraMassForGravity_Ptr( x, y, z, Time[0], 0, NULL );
            }}}
         }
      }
#     endif // #ifdef GRAVITY


//    copy data to the send buffer
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

#           ifdef GRAVITY
//          subtract the background density (which is assumed to be UNITY) for the isolated BC in the comoving frame
//          --> to be consistent with the comoving-frame Poisson eq.
#           ifdef COMOVING
            if ( OPT__BC_POT == BC_POT_ISOLATED )
            {
               for (int t=0; t<PSSize; t++)  TempBuf_Rho_Ptr[t] -= (real)1.0;
            }
#           endif
#           endif // #ifdef GRAVITY

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

} // FUNCTION : Patch2Slab_Rho



//-------------------------------------------------------------------------------------------------------
// Function    :  ZIndex2Rank
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
int ZIndex2Rank( const int IndexZ, const int *List_z_start, const int TRank_Guess )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef GRAVITY
   const int FFT_SizeZ = ( OPT__BC_POT == BC_POT_ISOLATED ) ? 2*NX0_TOT[2] : NX0_TOT[2];

   if ( List_z_start[MPI_NRank] < FFT_SizeZ )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "List_z_start[MPI_NRank]", List_z_start[MPI_NRank] );
#  endif // #ifdef GRAVITY

   if ( IndexZ < 0  ||  IndexZ >= NX0_TOT[2] )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IndexZ", IndexZ );

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



#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Slab2Patch_Pot
// Description :  Slab domain decomposition --> patch-based data (for potential)
//
// Parameter   :  RhoK       : In-place FFT array
//                SendBuf    : Sending MPI buffer of potential
//                RecvBuf    : Receiving MPI buffer of potential
//                SaveSg     : Sandglass to store the updated data
//                List_SIdx  : 1D coordinate in slab
//                List_PID   : PID of each patch slice sent to each rank
//                List_k     : Local z coordinate of each patch slice sent to each rank
//                List_NSend : Size of potential data sent to each rank
//                List_NRecv : Size of potential data received from each rank
//                local_nz   : Slab thickness of this MPI rank
//                FFT_Size   : Size of the FFT operation including the zero-padding regions
//                NSendSlice : Total number of z slices need to be sent to other ranks (could be zero in the isolated BC)
//-------------------------------------------------------------------------------------------------------
void Slab2Patch_Pot( const real *RhoK, real *SendBuf, real *RecvBuf, const int SaveSg, const long *List_SIdx,
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

} // FUNCTION : Slab2Patch_Pot
#endif // #ifdef GRAVITY



#endif // #ifdef SUPPORT_FFTW
