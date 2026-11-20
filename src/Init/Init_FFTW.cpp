#include "GAMER.h"

#ifdef SUPPORT_FFTW

static int ZIndex2Rank( const int IndexZ, const int *List_z_start, const int TRank_Guess );

// PS : plan for calculating the power spectrum
root_fftw::real_plan_nd FFTW_Plan_PS;
static void Init_FFTW_PowerSpectrum( const int StartupFlag );
static void End_FFTW_PowerSpectrum();

// Poi : plan for the self-gravity Poisson solver
#ifdef GRAVITY
root_fftw::real_plan_nd FFTW_Plan_Poi[NLEVEL], FFTW_Plan_Poi_Inv[NLEVEL];
static void Init_FFTW_Poisson( const int StartupFlag, const int lv );
static void End_FFTW_Poisson();
#endif // #ifdef GRAVITY

// Psi : plan for the ELBDM spectral solver
#if ( MODEL == ELBDM )
root_fftw::complex_plan_nd FFTW_Plan_Psi, FFTW_Plan_Psi_Inv;
static void Init_FFTW_ELBDMSpectral( const int StartupFlag );
static void End_FFTW_ELBDMSpectral();

// ExtPsi : plan for the Gram Fourier extension solver
#if ( WAVE_SCHEME == WAVE_GRAMFE )
gramfe_fftw::complex_plan_1d FFTW_Plan_ExtPsi, FFTW_Plan_ExtPsi_Inv;
static void Init_FFTW_GramFE( const int StartupFlag );
static void End_FFTW_GramFE();
#endif // #if ( WAVE_SCHEME == WAVE_GRAMFE )
#endif // #if ( MODEL == ELBDM )




#if ( SUPPORT_FFTW == FFTW3 )
//-------------------------------------------------------------------------------------------------------
// Function    :  ComputePaddedTotalSize
// Description :  Return padded total size for complex-to-real and real-to-complex 3D FFTW transforms
//
// Parameter   :  size : 3D array with size of FFT block
//
// Return      :  length of array that is large enough to store FFT input and output
//-------------------------------------------------------------------------------------------------------
size_t ComputePaddedTotalSize( int* size )
{

#  ifdef SERIAL
   return 2*( (size_t)size[0]/2 + 1 )*size[1]*size[2];
#  else
   mpi_index_int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose;

   return (size_t)fftw_mpi_local_size_3d_transposed( size[2], size[1], 2*(size[0]/2+1),
                                                     MPI_COMM_WORLD, &local_nz, &local_z_start, &local_ny_after_transpose,
                                                     &local_y_start_after_transpose);
#  endif // #if SERIAL

} // FUNCTION : ComputePaddedTotalSize



//-------------------------------------------------------------------------------------------------------
// Function    :  ComputeTotalSize
// Description :  Return total size for complex-to-complex 3D FFTW transforms
//
// Parameter   :  size : 3D array with size of FFT block
//
// Return      :  length of array that is large enough to store FFT input and output
//-------------------------------------------------------------------------------------------------------
size_t ComputeTotalSize( int* size )
{

#  ifdef SERIAL
   return (size_t)size[0]*size[1]*size[2];
#  else
   mpi_index_int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose;

   return (size_t)fftw_mpi_local_size_3d_transposed( size[2], size[1], size[0],
                                                     MPI_COMM_WORLD, &local_nz, &local_z_start, &local_ny_after_transpose,
                                                     &local_y_start_after_transpose);
#  endif // #if SERIAL

} // FUNCTION : ComputeTotalSize
#endif // #if ( SUPPORT_FFTW == FFTW3 )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW
// Description :  Create the FFTW plans
//
// Parameter   :  lv : Target level
//-------------------------------------------------------------------------------------------------------
void Init_FFTW( const int lv )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );

   if ( FFTW_Inited[lv] )   return;

   static bool FirstTime = true;

// determine how to initialise fftw plans
   int StartupFlag;

   switch ( OPT__FFTW_STARTUP )
   {
      case FFTW_STARTUP_ESTIMATE:    StartupFlag = FFTW_ESTIMATE;    break;
      case FFTW_STARTUP_MEASURE:     StartupFlag = FFTW_MEASURE;     break;
#     if ( SUPPORT_FFTW == FFTW3 )
      case FFTW_STARTUP_PATIENT:     StartupFlag = FFTW_PATIENT;     break;
#     endif // # if ( SUPPORT_FFTW == FFTW3 )

      default:                      Aux_Error( ERROR_INFO, "unrecognised FFTW startup option %d  !!\n", OPT__FFTW_STARTUP );
   } // switch ( OPT__FFTW_STARTUP )


   if ( FirstTime )
   {
#     if ( SUPPORT_FFTW == FFTW3 )
      FFTW3_Double_OMP_Enabled = false;
      FFTW3_Single_OMP_Enabled = false;


#     ifdef OPENMP

#     ifndef SERIAL
//    check the level of MPI thread support
      int MPI_Thread_Status;
      MPI_Query_thread( &MPI_Thread_Status );

//    enable multithreading if possible
      FFTW3_Double_OMP_Enabled = MPI_Thread_Status >= MPI_THREAD_FUNNELED;
      FFTW3_Single_OMP_Enabled = FFTW3_Double_OMP_Enabled;
#     else // # ifndef SERIAL

//    always enable multithreading in serial mode with openmp
      FFTW3_Double_OMP_Enabled = true;
      FFTW3_Single_OMP_Enabled = true;
#     endif // # ifndef SERIAL ... # else

//    initialise fftw multithreading
      if ( FFTW3_Double_OMP_Enabled )
      {
         FFTW3_Double_OMP_Enabled = fftw_init_threads();
         if ( !FFTW3_Double_OMP_Enabled )    Aux_Error( ERROR_INFO, "fftw_init_threads() failed !!\n" );
      }
      if ( FFTW3_Single_OMP_Enabled )
      {
         FFTW3_Single_OMP_Enabled = fftwf_init_threads();
         if ( !FFTW3_Single_OMP_Enabled )    Aux_Error( ERROR_INFO, "fftwf_init_threads() failed !!\n" );
      }
#     endif // # ifdef OPENMP

//    initialise fftw mpi support
#     ifndef SERIAL
      fftw_mpi_init();
      fftwf_mpi_init();
#     endif // # ifndef SERIAL

//    tell all subsequent fftw3 planners to use OMP_NTHREAD threads
#     ifdef OPENMP
      if ( FFTW3_Double_OMP_Enabled ) fftw_plan_with_nthreads ( OMP_NTHREAD );
      if ( FFTW3_Single_OMP_Enabled ) fftwf_plan_with_nthreads( OMP_NTHREAD );
#     endif // # ifdef OPENMP
#     endif // # if ( SUPPORT_FFTW == FFTW3 )

      FirstTime = false;
   } // if ( FirstTime )

   if ( lv == 0 )   Init_FFTW_PowerSpectrum( StartupFlag );

#  ifdef GRAVITY
   Init_FFTW_Poisson( StartupFlag, lv );
#  endif

#  if ( MODEL == ELBDM )
   if ( lv == 0 )   Init_FFTW_ELBDMSpectral( StartupFlag );

#  if ( WAVE_SCHEME == WAVE_GRAMFE )
   if ( lv == 0 )   Init_FFTW_GramFE( StartupFlag );
#  endif // #if ( WAVE_SCHEME == WAVE_GRAMFE )
#  endif // #if ( MODEL == ELBDM )

   FFTW_Inited[lv] = true;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_FFTW



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW
// Description :  Delete the FFTW plans
//-------------------------------------------------------------------------------------------------------
void End_FFTW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


   End_FFTW_PowerSpectrum();

#  ifdef GRAVITY
   End_FFTW_Poisson();
#  endif

#  if ( MODEL == ELBDM )
   End_FFTW_ELBDMSpectral();

#  if ( WAVE_SCHEME == WAVE_GRAMFE )
   End_FFTW_GramFE();
#  endif // # if ( WAVE_SCHEME == WAVE_GRAMFE )
#  endif // #if ( MODEL == ELBDM )


#  if ( SUPPORT_FFTW == FFTW3 )
#  ifdef OPENMP
   if ( FFTW3_Double_OMP_Enabled )  fftw_cleanup_threads();
   if ( FFTW3_Single_OMP_Enabled ) fftwf_cleanup_threads();
#  endif // # ifdef OPENMP

#  ifdef SERIAL
   fftw_cleanup();
   fftwf_cleanup();
#  else // # ifdef SERIAL
// fftwx_mpi_cleanup also calls fftwx_cleanup
   fftw_mpi_cleanup();
   fftwf_mpi_cleanup();
#  endif // # ifdef SERIAL ... # else
#  endif // # if ( SUPPORT_FFTW == FFTW3 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_FFTW



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW_PowerSpectrum
// Description :  Create the FFTW plan for the power spectrum
//
// Parameter   :  StartupFlag : Initialize FFTW plan method
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void Init_FFTW_PowerSpectrum( const int StartupFlag )
{

// determine the FFT size
   int PS_FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };

   real* PS = NULL;

// allocate memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   PS = (real*)root_fftw::fft_malloc( ComputePaddedTotalSize( PS_FFT_Size ) * sizeof(real) );
#  endif

// create plans
   FFTW_Plan_PS = root_fftw_create_3d_r2c_plan( PS_FFT_Size, PS, StartupFlag );

// free memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   root_fftw::fft_free( PS );
#  endif

} // FUNCITON : Init_FFTW_PowerSpectrum



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW_PowerSpectrum
// Description :  Delete the FFTW plan for the power spectrum
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void End_FFTW_PowerSpectrum()
{

   root_fftw::destroy_real_plan_nd( FFTW_Plan_PS );

} // FUNCITON : End_FFTW_PowerSpectrum



#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW_Poisson
// Description :  Create the FFTW plan for the self-gravity Poisson solver
//
// Parameter   :  StartupFlag : Initialize FFTW plan method
//                lv          : Target level
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void Init_FFTW_Poisson( const int StartupFlag, const int lv )
{

   const int CellFactor = (int)(1L<<lv);

// determine the FFT size
   int Gravity_FFT_Size[3] = { NX0_TOT[0]*CellFactor, NX0_TOT[1]*CellFactor, NX0_TOT[2]*CellFactor };

// the zero-padding method is adopted for the isolated BC
   if ( OPT__BC_POT == BC_POT_ISOLATED )
      for (int d=0; d<3; d++)   Gravity_FFT_Size[d] *= 2;

// check
   if ( MPI_Rank == 0 )
   for (int d=0; d<3; d++)
   {
      if ( Gravity_FFT_Size[d] <= 0 )   Aux_Error( ERROR_INFO, "Gravity_FFT_Size[%d] = %d < 0 !!\n", d, Gravity_FFT_Size[d] );
   }

   real* RhoK = NULL;

// allocate memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   RhoK = (real*)root_fftw::fft_malloc( ComputePaddedTotalSize( Gravity_FFT_Size ) * sizeof(real) );
#  endif

// create plans
   FFTW_Plan_Poi[lv]     = root_fftw_create_3d_r2c_plan( Gravity_FFT_Size, RhoK, StartupFlag );
   FFTW_Plan_Poi_Inv[lv] = root_fftw_create_3d_c2r_plan( Gravity_FFT_Size, RhoK, StartupFlag );

// free memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   root_fftw::fft_free( RhoK );
#  endif

} // FUNCTION : Init_FFTW_Poisson



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW_Poisson
// Description :  Delete the FFTW plan for the self-gravity Poisson solver
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void End_FFTW_Poisson()
{

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( ! FFTW_Inited[lv] )   continue;

      root_fftw::destroy_real_plan_nd( FFTW_Plan_Poi[lv]     );
      root_fftw::destroy_real_plan_nd( FFTW_Plan_Poi_Inv[lv] );
   }

} // FUNCITON : End_FFTW_Poisson
#endif // #ifdef GRAVITY



#if ( MODEL == ELBDM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW_ELBDMSpectral
// Description :  Create the FFTW plan for the ELBDM spectral solver
//
// Parameter   :  StartupFlag : Initialize FFTW plan method
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void Init_FFTW_ELBDMSpectral( const int StartupFlag )
{

// determine the FFT size
   int Psi_FFT_Size[3]    = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };
#  if ( defined( SERIAL )  ||  SUPPORT_FFTW == FFTW3 )
   int InvPsi_FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] };
#  else // # ifdef SERIAL || FFTW3
// Note that the dimensions of the inverse transform in FFTW2,
// which are given by the dimensions of the output of the forward transform,
// are Ny*Nz*Nx because we are using "FFTW_TRANSPOSED_ORDER" in fftwnd_mpi().
   int InvPsi_FFT_Size[3] = { NX0_TOT[0], NX0_TOT[2], NX0_TOT[1] };
#  endif // # ifdef SERIAL || FFTW3 ... else

   real* PsiK = NULL;

// allocate memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   PsiK = (real*)root_fftw::fft_malloc( ComputeTotalSize( Psi_FFT_Size ) * sizeof(real) * 2 );  // 2 * real for size of complex number
#  endif

// create plans
   FFTW_Plan_Psi     = root_fftw_create_3d_forward_c2c_plan ( Psi_FFT_Size,    PsiK, StartupFlag );
   FFTW_Plan_Psi_Inv = root_fftw_create_3d_backward_c2c_plan( InvPsi_FFT_Size, PsiK, StartupFlag );

// free memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   root_fftw::fft_free( PsiK );
#  endif

} // FUNCTION : Init_FFTW_ELBDMSpectral



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW_ELBDMSpectral
// Description :  Delete the FFTW plan for the ELBDM spectral solver
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void End_FFTW_ELBDMSpectral()
{

   root_fftw::destroy_complex_plan_nd( FFTW_Plan_Psi     );
   root_fftw::destroy_complex_plan_nd( FFTW_Plan_Psi_Inv );

} // FUNCITON : End_FFTW_ELBDMSpectral



#if ( WAVE_SCHEME == WAVE_GRAMFE )
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW_GramFE
// Description :  Create the FFTW plan for the Gram Fourier extension solver
//
// Note        :  1. The Gram-Fourier extension planners only use one thread because OMP parallelisation
//                   evolves different patches parallely.
//                   From the FFTW3 documentation: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
//                   --> "You can call fftw_plan_with_nthreads, create some plans, call fftw_plan_with_nthreads
//                        again with a different argument, and create some more plans for a new number of threads."
//
// Parameter   :  StartupFlag : Initialize FFTW plan method
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void Init_FFTW_GramFE( const int StartupFlag )
{

// determine the FFT size
   int ExtPsi_FFT_Size = GRAMFE_FLU_NXT;

   gramfe_fftw::fft_complex* ExtPsiK = NULL;

// allocate memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   ExtPsiK = (gramfe_fftw::fft_complex*)gramfe_fftw::fft_malloc( ExtPsi_FFT_Size * sizeof(gramfe_fftw::fft_complex) );
#  endif

// See note 1.
#  if ( defined( SUPPORT_FFTW3 )  &&  defined( OPENMP ) )
   if ( FFTW3_Double_OMP_Enabled )  fftw_plan_with_nthreads( 1 );
   if ( FFTW3_Single_OMP_Enabled ) fftwf_plan_with_nthreads( 1 );
#  endif

// create plans
   FFTW_Plan_ExtPsi      = gramfe_fftw_create_1d_forward_c2c_plan ( ExtPsi_FFT_Size, ExtPsiK, StartupFlag );
   FFTW_Plan_ExtPsi_Inv  = gramfe_fftw_create_1d_backward_c2c_plan( ExtPsi_FFT_Size, ExtPsiK, StartupFlag );

// restore regular settings
#  if ( defined( SUPPORT_FFTW3 )  &&  defined( OPENMP ) )
   if ( FFTW3_Double_OMP_Enabled )  fftw_plan_with_nthreads( OMP_NTHREAD );
   if ( FFTW3_Single_OMP_Enabled ) fftwf_plan_with_nthreads( OMP_NTHREAD );
#  endif

// free memory for arrays in fftw3
#  if ( SUPPORT_FFTW == FFTW3 )
   gramfe_fftw::fft_free( ExtPsiK );
#  endif

} // FUNCTION : Init_FFTW_Poisson



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW_GramFE
// Description :  Delete the FFTW plan for the Gram Fourier extension solver
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void End_FFTW_GramFE()
{

   gramfe_fftw::destroy_complex_plan_1d( FFTW_Plan_ExtPsi     );
   gramfe_fftw::destroy_complex_plan_1d( FFTW_Plan_ExtPsi_Inv );

} // FUNCITON : End_FFTW_Poisson
#endif // #if ( WAVE_SCHEME == WAVE_GRAMFE )
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Patch2Slab
// Description :  Patch-based data --> slab domain decomposition
//
// Note        :  1. List_PID[] and List_k[] will be allocated here; user needs to either call Slab2Patch()
//                   to free the memory or do manual deallocation
//
// Parameter   :  VarS           : Slab array of target variable for FFT
//                SendBuf_Var    : Sending MPI buffer of the target field
//                RecvBuf_Var    : Receiving MPI buffer of the target field
//                SendBuf_SIdx   : Sending MPI buffer of 1D coordinate in slab
//                RecvBuf_SIdx   : Receiving MPI buffer of 1D coordinate in slab
//                List_PID       : PID of each patch slice sent to each rank
//                List_k         : Local z coordinate of each patch slice sent to each rank
//                List_NSend_Var : Size of data sent to each rank
//                List_NRecv_Var : Size of data received from each rank
//                List_z_start   : Starting z coordinate of each rank in the FFTW slab decomposition
//                local_nz       : Slab thickness of this MPI rank
//                FFT_Size       : Size of the FFT operation including the zero-padding regions
//                NRecvSlice     : Total number of z slices received from other ranks (could be zero in the isolated BC)
//                PrepTime       : Physical time for preparing the target variable field
//                TVar           : Target variable to be prepared
//                InPlacePad     : Whether or not to pad the array size for in-place real-to-complex FFT
//                ForPoisson     : Preparing the density field for the Poisson solver
//                AddExtraMass   : Adding an extra density field for computing gravitational potential (only works with ForPoisson)
//                lv             : Target level
//-------------------------------------------------------------------------------------------------------
void Patch2Slab( real *VarS, real *SendBuf_Var, real *RecvBuf_Var, long *SendBuf_SIdx, long *RecvBuf_SIdx,
                 int **List_PID, int **List_k, long *List_NSend_Var, long *List_NRecv_Var,
                 const int *List_z_start, const int local_nz, const int FFT_Size[], const int NRecvSlice,
                 const double PrepTime, const long TVar, const bool InPlacePad, const bool ForPoisson,
                 const bool AddExtraMass, const int lv )
{

// check
// check only single field
   if ( TVar == 0  ||  TVar & (TVar-1) )
      Aux_Error( ERROR_INFO, "number of target variables is not one !!\n" );

#  ifdef GRAVITY
// check
   if ( ForPoisson  &&  TVar != _TOTAL_DENS )
      Aux_Error( ERROR_INFO, "TVar != _TOTAL_DENS for Poisson solver !!\n" );

   if ( ForPoisson  &&  AddExtraMass  &&  Poi_AddExtraMassForGravity_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Poi_AddExtraMassForGravity_Ptr == NULL for AddExtraMass !!\n" );
#  endif // GRAVITY

#  ifdef GAMER_DEBUG
   if ( List_z_start[MPI_Rank+1] - List_z_start[MPI_Rank] != local_nz )
      Aux_Error( ERROR_INFO, "local_nz (%d) != expectation (%d) !!\n",
                 local_nz, List_z_start[MPI_Rank+1] - List_z_start[MPI_Rank] );
#  endif // GAMER_DEBUG


   const int SSize[2]   = { ( InPlacePad ? 2*(FFT_Size[0]/2+1) : FFT_Size[0] ), FFT_Size[1] }; // padded slab size in the x and y directions
   const int PSSize     = PS1*PS1;                                                             // patch slice size
   const int MemUnit    = MAX( 1, amr->NPatchComma[lv][1]/MPI_NRank );                         // set arbitrarily; divided by MPI_NRank so that the
                                                                                               // memory consumption per rank for arrays like
                                                                                               // TempBuf_Var[MPI_NRank] reduces when launching
                                                                                               // more ranks
   const int AveNz      = FFT_Size[2]/MPI_NRank + ( ( FFT_Size[2]%MPI_NRank == 0 ) ? 0 : 1 );  // average slab thickness
   const int Scale      = amr->scale[lv];

   int   Cr[3];                        // corner coordinates of each patch normalized to the base-level grid size
   int   BPos_z;                       // z coordinate of each patch slice in the simulation box
   int   SPos_z;                       // z coordinate of each patch slice in the slab
   long  SIdx;                         // 1D coordinate of each patch slice in the slab
   int   List_NSend_SIdx[MPI_NRank];   // number of patch slices sent to each rank
   int   List_NRecv_SIdx[MPI_NRank];   // number of patch slices received from each rank
   long *TempBuf_SIdx   [MPI_NRank];   // 1D slab coordinate of each patch slice sent to each rank
   real *TempBuf_Var    [MPI_NRank];   // data of each patch slice sent to each rank
   real *TempBuf_Var_Ptr = NULL;

   int TRank, TRank_Guess, MemSize[MPI_NRank], idx;


// 1. set memory allocation unit
   for (int r=0; r<MPI_NRank; r++)
   {
      MemSize        [r] = MemUnit;
      List_PID       [r] = (int* )malloc( MemSize[r]*sizeof(int)         );
      List_k         [r] = (int* )malloc( MemSize[r]*sizeof(int)         );
      TempBuf_SIdx   [r] = (long*)malloc( MemSize[r]*sizeof(long)        );
      TempBuf_Var    [r] = (real*)malloc( MemSize[r]*sizeof(real)*PSSize );
      List_NSend_SIdx[r] = 0;

      if ( List_PID    [r] == NULL )
          Aux_Error( ERROR_INFO, "List_PID[%d] is NULL on Rank %d !!\n",     r, MPI_Rank );
      if ( List_k      [r] == NULL )
          Aux_Error( ERROR_INFO, "List_k[%d] is NULL on Rank %d !!\n",       r, MPI_Rank );
      if ( TempBuf_SIdx[r] == NULL )
          Aux_Error( ERROR_INFO, "TempBuf_SIdx[%d] is NULL on Rank %d !!\n", r, MPI_Rank );
      if ( TempBuf_Var [r] == NULL )
          Aux_Error( ERROR_INFO, "TempBuf_Var[%d] is NULL on Rank %d !!\n",  r, MPI_Rank );
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
   const long        CellFactor        = (long)(1L<<lv);

   real (*VarPatch)[PS1][PS1][PS1] = new real [8*NPG][PS1][PS1][PS1];

   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    even with NSIDE_00 and GhostSize=0, we still need OPT__BC_FLU to determine whether periodic BC is adopted
//    for depositing particle mass onto grids.
//    also note that we do not check minimum density here since no ghost zones are required
      Prepare_PatchData( lv, PrepTime, VarPatch[0][0][0], NULL, GhostSize, NPG, &PID0, TVar, _NONE,
                         IntScheme, INT_NONE, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


#     ifdef GRAVITY
//    add extra mass source for gravity if required
      if ( ForPoisson  &&  AddExtraMass )
      {
         const double dh = amr->dh[lv];

         for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
         {
            const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
            const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
            const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

            double x, y, z;

            for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
            for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
            for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;
               VarPatch[LocalID][k][j][i] += Poi_AddExtraMassForGravity_Ptr( x, y, z, Time[lv], lv, NULL );
            }}}
         }
      }
#     endif // #ifdef GRAVITY


//    copy data to the send buffer
      for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
      {
         for (int d=0; d<3; d++)    Cr[d] = amr->patch[0][lv][PID]->corner[d] / Scale;

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
               TempBuf_Var [TRank]  = (real*)realloc( TempBuf_Var [TRank], MemSize[TRank]*sizeof(real)*PSSize );

               if ( List_PID    [TRank] == NULL )
                  Aux_Error( ERROR_INFO, "List_PID[%d] is NULL on Rank %d !!\n",     TRank, MPI_Rank );
               if ( List_k      [TRank] == NULL )
                  Aux_Error( ERROR_INFO, "List_k[%d] is NULL on Rank %d !!\n",       TRank, MPI_Rank );
               if ( TempBuf_SIdx[TRank] == NULL )
                  Aux_Error( ERROR_INFO, "TempBuf_SIdx[%d] is NULL on Rank %d !!\n", TRank, MPI_Rank );
               if ( TempBuf_Var [TRank] == NULL )
                  Aux_Error( ERROR_INFO, "TempBuf_Var[%d] is NULL on Rank %d !!\n",  TRank, MPI_Rank );
            }

//          record list
            List_PID    [TRank][ List_NSend_SIdx[TRank] ] = PID;
            List_k      [TRank][ List_NSend_SIdx[TRank] ] = k;
            TempBuf_SIdx[TRank][ List_NSend_SIdx[TRank] ] = SIdx;

//          store data
            TempBuf_Var_Ptr = TempBuf_Var[TRank] + List_NSend_SIdx[TRank]*PSSize;

            idx = 0;
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
               TempBuf_Var_Ptr[ idx ++ ] = VarPatch[LocalID][k][j][i];

#           ifdef GRAVITY
//          subtract the background density (which is assumed to be UNITY) for the isolated BC in the comoving frame
//          --> to be consistent with the comoving-frame Poisson eq.
#           ifdef COMOVING
            if ( ForPoisson  &&  OPT__BC_POT == BC_POT_ISOLATED )
            {
               for (int t=0; t<PSSize; t++)  TempBuf_Var_Ptr[t] -= (real)1.0;
            }
#           endif
#           endif // #ifdef GRAVITY

            List_NSend_SIdx[TRank] ++;
         } // for (int k=0; k<PS1; k++)
      } // for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++)
   } // for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8)

   delete [] VarPatch;


// 3. prepare the send buffer
   int   Send_Disp_SIdx[MPI_NRank], Recv_Disp_SIdx[MPI_NRank];
   long  Send_Disp_Var[MPI_NRank], Recv_Disp_Var[MPI_NRank];
   long *SendPtr_SIdx = NULL;
   real *SendPtr_Var  = NULL;

// 3.1 broadcast the number of elements sending to different ranks
   MPI_Alltoall( List_NSend_SIdx, 1, MPI_INT, List_NRecv_SIdx, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   {
      List_NSend_Var[r] = (long)List_NSend_SIdx[r]*(long)PSSize;
      List_NRecv_Var[r] = (long)List_NRecv_SIdx[r]*(long)PSSize;
   }

// 3.2 calculate the displacement
   Send_Disp_SIdx[0] = 0;
   Recv_Disp_SIdx[0] = 0;
   Send_Disp_Var [0] = 0L;
   Recv_Disp_Var [0] = 0L;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_SIdx[r] = Send_Disp_SIdx[r-1] + List_NSend_SIdx[r-1];
      Recv_Disp_SIdx[r] = Recv_Disp_SIdx[r-1] + List_NRecv_SIdx[r-1];
      Send_Disp_Var [r] = Send_Disp_Var [r-1] + List_NSend_Var [r-1];
      Recv_Disp_Var [r] = Recv_Disp_Var [r-1] + List_NRecv_Var [r-1];
   }

// check
#  ifdef GAMER_DEBUG
   const long NSend_Total  = Send_Disp_Var[MPI_NRank-1] + List_NSend_Var[MPI_NRank-1];
   const long NRecv_Total  = Recv_Disp_Var[MPI_NRank-1] + List_NRecv_Var[MPI_NRank-1];
   const long NSend_Expect = (long)amr->NPatchComma[lv][1]*(long)CUBE(PS1);
   const long NRecv_Expect = (long)NX0_TOT[0]*CellFactor*(long)NX0_TOT[1]*CellFactor*(long)NRecvSlice;

   if ( NSend_Total != NSend_Expect )  Aux_Error( ERROR_INFO, "NSend_Total = %ld != expected value = %ld !!\n",
                                                  NSend_Total, NSend_Expect );

   if ( NRecv_Total != NRecv_Expect )  Aux_Error( ERROR_INFO, "NRecv_Total = %ld != expected value = %ld !!\n",
                                                  NRecv_Total, NRecv_Expect );
#  endif

// 3.3 prepare the send buffer of SIdx
   SendPtr_SIdx = SendBuf_SIdx;
   for (int r=0; r<MPI_NRank; r++)
   {
      memcpy( SendPtr_SIdx, TempBuf_SIdx[r], List_NSend_SIdx[r]*sizeof(long) );
      SendPtr_SIdx += List_NSend_SIdx[r];
   }

// 3.4 prepare the send buffer of the target field
   SendPtr_Var = SendBuf_Var;
   for (int r=0; r<MPI_NRank; r++)
   {
      memcpy( SendPtr_Var, TempBuf_Var[r], List_NSend_Var[r]*sizeof(real) );
      SendPtr_Var += List_NSend_Var[r];
   }


// 4. exchange data by MPI
   MPI_Alltoallv      ( SendBuf_SIdx, List_NSend_SIdx, Send_Disp_SIdx, MPI_LONG,
                        RecvBuf_SIdx, List_NRecv_SIdx, Recv_Disp_SIdx, MPI_LONG,       MPI_COMM_WORLD );

   MPI_Alltoallv_GAMER( SendBuf_Var,  List_NSend_Var,  Send_Disp_Var,  MPI_GAMER_REAL,
                        RecvBuf_Var,  List_NRecv_Var,  Recv_Disp_Var,  MPI_GAMER_REAL, MPI_COMM_WORLD );


// 5. store the received data to the padded array "VarS" for FFTW
   const long NPSlice = (long)NX0_TOT[0]*CellFactor*NX0_TOT[1]*CellFactor*NRecvSlice/PSSize;  // total number of received patch slices
   long  dSIdx, Counter = 0;
   real *VarS_Ptr = NULL;

   for (long t=0; t<NPSlice; t++)
   {
      SIdx     = RecvBuf_SIdx[t];
      VarS_Ptr = VarS + SIdx;

      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         dSIdx           = j*SSize[0] + i;
         VarS_Ptr[dSIdx] = RecvBuf_Var[ Counter ++ ];
      }
   }


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( TempBuf_SIdx[r] );
      free( TempBuf_Var [r] );
   }

} // FUNCTION : Patch2Slab



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
// Description :  Slab domain decomposition --> patch-based data
//
// Parameter   :  VarS       : Slab array of target variable after FFT
//                SendBuf    : Sending MPI buffer of the target field
//                RecvBuf    : Receiving MPI buffer of the target field
//                SaveSg     : Sandglass to store the updated data
//                List_SIdx  : 1D coordinate in slab
//                List_PID   : PID of each patch slice sent to each rank
//                List_k     : Local z coordinate of each patch slice sent to each rank
//                List_NSend : Size of data sent to each rank
//                List_NRecv : Size of data received from each rank
//                local_nz   : Slab thickness of this MPI rank
//                FFT_Size   : Size of the FFT operation including the zero-padding regions
//                NSendSlice : Total number of z slices need to be sent to other ranks (could be zero in the isolated BC)
//                TVar       : Target variable to be prepared
//                InPlacePad : Whether or not to pad the array size for in-place real-to-complex FFT
//                lv         : Target level
//-------------------------------------------------------------------------------------------------------
void Slab2Patch( const real *VarS, real *SendBuf, real *RecvBuf, const int SaveSg, const long *List_SIdx,
                 int **List_PID, int **List_k, long *List_NSend, long *List_NRecv, const int local_nz, const int FFT_Size[],
                 const int NSendSlice, const long TVar, const bool InPlacePad, const int lv )
{

// check
// check only single field
   if ( TVar == 0  ||  TVar & (TVar-1) )
      Aux_Error( ERROR_INFO, "number of target variables is not one !!\n" );

// check TVar is one of the fields in NCOMP_TOTAL or POTE
   long SupportedVar = _TOTAL;
#  ifdef GRAVITY
   SupportedVar |= _POTE;
#  endif

   if ( TVar & ~SupportedVar )
      Aux_Error( ERROR_INFO, "unsupported variable %s = %d !!\n", "TVar", TVar );

// find the variable index
   int TVarIdx = -1;   // target variable index
   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVar & (1L<<v) )  TVarIdx = v;  // assuming only single field

#  ifdef GRAVITY
   if ( TVar == _POTE )  TVarIdx = NCOMP_TOTAL+NDERIVE;  // assuming only single field
#  endif

   if ( TVarIdx < 0 )
      Aux_Error( ERROR_INFO, "TVarIdx is not found !!\n" );


// 1. store the evaluated data to the send buffer
   const long  CellFactor = (long)(1L<<lv);
   const int   SSize[2]   = { ( InPlacePad ? 2*(FFT_Size[0]/2+1) : FFT_Size[0] ), FFT_Size[1] };  // padded slab size in the x and y directions
   const int   PSSize     = PS1*PS1;                                                              // patch slice size
   const long  NPSlice    = (long)NX0_TOT[0]*CellFactor*NX0_TOT[1]*CellFactor*NSendSlice/PSSize;  // total number of patch slices to be sent
   const real *VarS_Ptr   = NULL;

   long SIdx, dSIdx, Counter = 0;

   for (long t=0; t<NPSlice; t++)
   {
      SIdx     = List_SIdx[t];
      VarS_Ptr = VarS + SIdx;

      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         dSIdx                 = j*SSize[0] + i;
         SendBuf[ Counter ++ ] = VarS_Ptr[dSIdx];
      }
   }


// 2. calculate the displacement and exchange data by MPI
   long Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];

   Send_Disp[0] = 0L;
   Recv_Disp[0] = 0L;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + List_NSend[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + List_NRecv[r-1];
   }

   MPI_Alltoallv_GAMER( SendBuf, List_NSend, Send_Disp, MPI_GAMER_REAL,
                        RecvBuf, List_NRecv, Recv_Disp, MPI_GAMER_REAL, MPI_COMM_WORLD );


// 3. store the received data to different patch objects
   int   PID, k, NRecvSlice;
   real *RecvPtr = RecvBuf;

   for (int r=0; r<MPI_NRank; r++)
   {
      NRecvSlice = int(List_NRecv[r]/(long)PSSize);

      for (int t=0; t<NRecvSlice; t++)
      {
         PID = List_PID[r][t];
         k   = List_k  [r][t];

         if ( TVarIdx < NCOMP_TOTAL )
            memcpy( amr->patch[SaveSg][lv][PID]->fluid[TVarIdx][k], RecvPtr, PSSize*sizeof(real) );
#        ifdef GRAVITY
         else if ( TVarIdx == NCOMP_TOTAL+NDERIVE ) // TVar == _POTE
            memcpy( amr->patch[SaveSg][lv][PID]->pot[k], RecvPtr, PSSize*sizeof(real) );
#        endif
         else
            Aux_Error( ERROR_INFO, "incorrect target variable index %s = %d !!\n", "TVarIdx", TVarIdx );

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



#endif // #ifdef SUPPORT_FFTW
