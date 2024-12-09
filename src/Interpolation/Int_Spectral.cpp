#include "GAMER.h"

#ifdef SUPPORT_SPECTRAL_INT


/*******************************************************************************************************/
/*******************************************************************************************************/
/**  How to recreate the tables used in the below interpolation routines?                             **/
/**                                                                                                   **/
/**  They were created using the Python script                                                        **/
/**  tool/table_maker/GramFE/compute_interpolation_tables.py. The script can benefit from up to 42    **/
/**  MPI ranks in parallel and requires the mpi4py library.                                           **/
/**  The Python code implements the Gram Fourier extension algorithm using arbitrary precision        **/
/**  arithmetic and outputs interpolation and extension tables. The idea is to periodically extend    **/
/**  the interpolant and then interpolate it using a DFT - either precomputed as matrix               **/
/**  (PrecomputedInterpolationContext) or via the FFT algorithm (GramFEInterpolationContext)          **/
/**                                                                                                   **/
/**  The tables for the PrecomputedInterpolationContext are computed using the GramFEInterpolation    **/
/**  class and put in the folder "interpolation_tables". The tables for the                           **/
/**  GramFEInterpolationContext are computed using the GramFEFixedSizeExtension class and put in the  **/
/**  folder "boundary2extension_tables". The respective constructors take the parameters              **/
/**  Gamma (number of sampling points), g (number of Fourier modes), nd (size of extension region),   **/
/**  and m (the order of boundary polynomials).                                                       **/
/**  The tables use the Gram Fourier extension algorithm to find a periodic extension of a given      **/
/**  non-periodic function - the interpolant. This periodic function can then be interpolated using   **/
/**  a DFT. This interpolation is not monotonic and not conservative, but highly accurate.            **/
/**  Its accuracy is determined by the quality of the periodic extension.                             **/
/**  Gamma = 150, g = 63, nd = 32 and m = 8 have been found to be a good compromise between stability **/
/**  and accuracy for interpolation. Higher values of m led to interpolation artifacts in some tests. **/
/**  Lower values of m should in general be less accurate, but lead to better-behaved interpolation.  **/
/**                                                                                                   **/
/**  interpolation_tables:                                                                            **/
/**     The tables are named N.bin where N is the input interpolation size.                           **/
/**     They contain double-precision 2D matrices that represent the complete interpolation           **/
/**     GramFE extension * DFT * shift in k space * IDFT as a matrix.                                 **/
/**     It is an efficient way to implement the Gram FE interpolation for small N                     **/
/**     with a cost of N^2 operations for a single matrix multiplication.                             **/
/**     For m = N, this interpolation is equivalent to interpolating with a polynomial of             **/
/**     order N (as observed in numerical tests).                                                     **/
/**                                                                                                   **/
/**  boundary2extension_tables:                                                                       **/
/**     The tables are named %ld_%ld_%d_%d.bin with the arguments nd, m, Gamma, g.                    **/
/**     They contain double-precision 2D matrices that represent GramFE extension as a matrix.        **/
/**     The interpolation itself is carried out using the FFTW library which is faster for larger     **/
/**     interpolants with an asymptotic cost of 2 * mxm matrix multiplication + N log N.              **/
/**                                                                                                   **/
/*******************************************************************************************************/
/*******************************************************************************************************/

#include <complex>
#include <memory>
#include "GramFE_Interpolation.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Int_Spectral
// Description :  Perform spatial interpolation based on the Gram-Fourier extension method
//
// Note        :  1. The interpolation is spectrally accurate
//                2. The interpolation result is neither conservative nor monotonic
//                3. 3D interpolation is achieved by performing interpolation along x, y, and z directions
//                   in order
//
// Parameter   :  CData           : Input coarse-grid array
//                CSize           : Size of the CData array
//                CStart          : (x,y,z) starting indices to perform interpolation on the CData array
//                CRange          : Number of grids in each direction to perform interpolation
//                FData           : Output fine-grid array
//                FStart          : (x,y,z) starting indcies to store the interpolation results
//                NComp           : Number of components in the CData and FData array
//                UnwrapPhase     : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//                                  This option requires that the input field includes only the density and the phase field in this order
//                Monotonic       : Ensure that all interpolation results are monotonic
//                MonoCoeff       : Slope limiter coefficient for the option "Monotonic"
//                OppSign0thOrder : See Int_MinMod1D()
//-------------------------------------------------------------------------------------------------------
void Int_Spectral( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                   real FData[], const int FSize[3], const int FStart[3], const int NComp,
                   const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder )
{

// interpolation-scheme-dependent parameters
// ===============================================================================
// number of coarse-grid ghost zone
   const int CGhost = SPEC_INT_GHOST_BOUNDARY;
// ===============================================================================


// allocate memory to store largest input and output column
   const size_t MaxSize    = MAX( MAX(CRange[0],CRange[1]), CRange[2] );
   const size_t InputDisp  = MaxSize + 2*CGhost;
   const size_t OutputDisp = 2*MaxSize;

// add interpolation contexts (make sure that everything is thread-safe in AddInterpolationContext)
   for (size_t dim=0; dim<3; ++dim) {
      Int_InterpolationHandler.AddInterpolationContext( CRange[dim] + 2*CGhost, CGhost );
   }

// determine workspace size
   size_t WorkspaceSize = 0;
   for (size_t dim=0; dim<3; ++dim) {
      WorkspaceSize = MAX(  WorkspaceSize, Int_InterpolationHandler.GetWorkspaceSize( CRange[dim] + 2*CGhost, CGhost )  );
   }

// allocate workspace using fftw_malloc for contiguous memory used in FFT
   char* Workspace = (char*) gamer_fftw::fft_malloc( WorkspaceSize );


// index stride of the coarse-grid input array
   const size_t Cdx     = 1;
   const size_t Cdy     = Cdx*CSize[0];
   const size_t Cdz     = Cdy*CSize[1];

// index stride of the temporary arrays storing the data after x and y interpolations
   const size_t Tdx     = 1;
   const size_t Tdy     = Tdx* CRange[0]*2;
   const size_t TdzX    = Tdy*(CRange[1]+2*CGhost);   // array after x interpolation
   const size_t TdzY    = Tdy* CRange[1]*2;           // array after y interpolation

// index stride of the fine-grid output array
   const size_t Fdx     = 1;
   const size_t Fdy     = Fdx*FSize[0];
   const size_t Fdz     = Fdy*FSize[1];

// index stride of different components
   const int CDisp      = CSize[0]*CSize[1]*CSize[2];
   const int FDisp      = FSize[0]*FSize[1]*FSize[2];
   const int TDataXDisp = (CRange[2]+2*CGhost)*TdzX;
   const int TDataYDisp = (CRange[2]+2*CGhost)*TdzY;

// x = 0, y = 1, z = 2
// first interpolate x, then y, then z
   const size_t P[3][3] = { {2, 1, 0},    // z, y, x
                            {2, 0, 1},    // z, x, y
                            {1, 0, 2} };  // x, y, z

// starting indices for 3D loops [x, y, z][1st iteration, 2nd iteration, 3rd iteration]
   const int StartIndex[3][3] = { {CStart[0]-CGhost, 0, 0},
                                  {CStart[1]-CGhost, 0, 0},
                                  {CStart[2]-CGhost, 0, 0} };

// end indices for 3D loops [x, y, z][1st iteration, 2nd iteration, 3rd iteration]
   const int EndIndex[3][3] = { {CStart[0]+CRange[0]+CGhost, CRange[0]+CRange[0], CRange[0]+CRange[0]},
                                {CStart[1]+CRange[1]+CGhost, CRange[1]+2*CGhost,  CRange[1]+CRange[1]},
                                {CStart[2]+CRange[2]+CGhost, CRange[2]+2*CGhost,  CRange[2]+2*CGhost } };

// index strides for input arrays [x, y, z][1st iteration, 2nd iteration, 3rd iteration]
   const size_t InStride[3][3] = { {Cdx, Tdx,  Tdx },
                                   {Cdy, Tdy,  Tdy },
                                   {Cdz, TdzX, TdzY} };

// index strides for output arrays [x, y, z][1st iteration, 2nd iteration, 3rd iteration]
   const size_t OutStride[3][3] = { {Tdx,  Tdx,  Fdx},
                                    {Tdy,  Tdy,  Fdy},
                                    {TdzX, TdzY, Fdz} };

// input sizes [x, y, z] (coarse array + ghost zone)
   const int InSize[3] = { CRange[0]+2*CGhost, CRange[1]+2*CGhost, CRange[2]+2*CGhost };

// output sizes [x, y, z] (2*coarse array)
   const int OutSize[3] = { CRange[0]+CRange[0], CRange[1]+CRange[1], CRange[2]+CRange[2] };


   real *Input    = new real [ NComp*InputDisp  ]; // hold one column of input data
   real *Output   = new real [ NComp*OutputDisp ]; // hold one column of output data
   real *TDataX   = new real [ NComp*TDataXDisp ]; // temporary array after x interpolation
   real *TDataY   = new real [ NComp*TDataYDisp ]; // temporary array after y interpolation
   real *CPtr     = CData;
   real *FPtr     = FData;
   real *InPtr3D  = NULL;
   real *OutPtr3D = NULL;
   real *InPtr1D  = NULL;
   real *OutPtr1D = NULL;

#  if ( MODEL == ELBDM )
#  ifdef GAMER_DEBUG
   if ( UnwrapPhase  &&  NComp != 2 ) {
      Aux_Error( ERROR_INFO, "NComp = %d != 2 for OPT__INT_PHASE and INT_SPECTRAL !!\n", NComp );
   }
#  endif

   real *DensInput  = NULL;
   real *DensOutput = NULL;
#  endif // #if ( MODEL == ELBDM )

   for (size_t XYZ=0; XYZ<3; ++XYZ)
   {
//    map indices i, j, k to dimensions x, y, z
      const size_t Pi = P[XYZ][0];
      const size_t Pj = P[XYZ][1];
      const size_t Pk = P[XYZ][2];

      for (int ii=StartIndex[Pi][XYZ], io=0; ii<EndIndex[Pi][XYZ]; ii++, io++)
      for (int ji=StartIndex[Pj][XYZ], jo=0; ji<EndIndex[Pj][XYZ]; ji++, jo++)
      {
         for (int v=0; v<NComp; v++)
         {
//          set input 3D arrays during each iteration
            switch ( XYZ )
            {
               case 0: InPtr3D = CPtr   + v*CDisp;       break;
               case 1: InPtr3D = TDataX + v*TDataXDisp;  break;
               case 2: InPtr3D = TDataY + v*TDataYDisp;  break;
            }

            InPtr1D = Input + v*InputDisp;

            for (int ki=StartIndex[Pk][XYZ], ko=0; ki<EndIndex[Pk][XYZ]; ki++, ko++)
            {
               const size_t IdxIn = ii*InStride[Pi][XYZ] + ji*InStride[Pj][XYZ] + ki*InStride[Pk][XYZ];

               InPtr1D[ko] = InPtr3D[IdxIn];
            } // ki
         } // v


#        if ( MODEL == ELBDM )
         bool InterpolateReIm = false;

         if ( UnwrapPhase )
         {
            real* Real = Input + 0*InputDisp;
            real* Imag = Input + 1*InputDisp;

//          unwrap phase
            for (int k=1; k<InSize[XYZ]; k++)
            {
               Imag[k] = ELBDM_UnwrapPhase( Imag[k-1], Imag[k] );
            }

//          convert density and phase to real and imaginary part if vortex is detected
//          NOTE:
//          Two other strategies that were thought to improve the interpolation around vortices have been thoroughly tested
//          Strategy 1: Interpolate rho and rho * cos(S/N)
//          This strategy performs poorly because rho's derivatives have a kink at the vortex.
//          Taking powers of rho to shift the kink to higher derivatives helps, but increases the error when taking a high root to obtain rho again
//
//          Strategy 2:
//          Perform a singular gauge transformation rho -> -rho and S -> S +- pi at the vortex
//          This leads to smooth fields at the vortex, but fails around the vortex where S might jump by less then pi and the density is non-zero
//
//          The safest strategy is to simply interpolate the real and imaginary part around the vortex with a generous vortex detection threshold
            if ( SPEC_INT_XY_INSTEAD_DEPHA )
            {
//             detect phase discontinuities
               for (int k=1; k<InSize[XYZ]-1; k++) {
//                assuming Lap(S) * dx**2 > threshold implies a significant phase jump
                  if (  FABS( Imag[k+1] - (real)2.0*Imag[k] + Imag[k-1] ) > SPEC_INT_VORTEX_THRESHOLD  ) {
                     InterpolateReIm = true;
                     break;
                  }
               }

//             convert back to real & imaginary part
               if ( InterpolateReIm ) {
                  for (int k=0; k<InSize[XYZ]; k++) {
                     const real f  = SQRT( Real[k] );
                     const real Re = f*COS( Imag[k] );
                     const real Im = f*SIN( Imag[k] );
                     Real[k] = Re;
                     Imag[k] = Im;
                  }
               }
            } // if ( SPEC_INT_XY_INSTEAD_DEPHA )
         } // if ( UnwrapPhase )
#        endif // #if ( MODEL == ELBDM )

//       interpolate data
         for (int v=0; v<NComp; v++)
         {
            InPtr1D  = Input  + v*InputDisp;
            OutPtr1D = Output + v*OutputDisp;
            Int_InterpolationHandler.InterpolateReal( InPtr1D, OutPtr1D, InSize[XYZ], CGhost, Workspace, Monotonic[v], MonoCoeff, OppSign0thOrder );
         }

#        if ( MODEL == ELBDM )
         if ( UnwrapPhase )
         {
            real* Real = Output + 0*OutputDisp;
            real* Imag = Output + 1*OutputDisp;

            if ( SPEC_INT_XY_INSTEAD_DEPHA ) {
//             convert back to density and phase from real/imag
               if ( InterpolateReIm ) {
                  for (int k=0; k<OutSize[XYZ]; k++) {
                     const real De  = SQR( Real[k] ) + SQR( Imag[k] );
                     const real Pha = SATAN2( Imag[k], Real[k] );
                     Real[k] = De;
                     Imag[k] = Pha;
                  }

//                unwrap phase to be consistent with the results of interpolating density and phase
                  for (int k=1; k<OutSize[XYZ]; k++) {
                     Imag[k] = ELBDM_UnwrapPhase( Imag[k-1], Imag[k] );
                  }
               }
            }

//          density floor
            for (int k=0; k<OutSize[XYZ]; k++) {
               Real[k] = MAX( TINY_NUMBER, Real[k] );
            } // if ( SPEC_INT_XY_INSTEAD_DEPHA )
         } // if ( UnwrapPhase )
#        endif // #if ( MODEL == ELBDM )

//       write result of interpolation (excluding ghost zones) to temporary array
         for (int v=0; v<NComp; v++)
         {
//          set output 3D arrays during each iteration
            switch ( XYZ )
            {
               case 0: OutPtr3D = TDataX + v*TDataXDisp;    break;
               case 1: OutPtr3D = TDataY + v*TDataYDisp;    break;
               case 2: OutPtr3D = FPtr   + v*FDisp;         break;
            }

            OutPtr1D = Output + v*OutputDisp;

            for (int ki=0, ko=0; ki<CRange[XYZ]; ki++, ko+=2)
            {
               const size_t Idx_Out = io*OutStride[Pi][XYZ] + jo*OutStride[Pj][XYZ] + ko*OutStride[Pk][XYZ];
               OutPtr3D[ Idx_Out                      ] = OutPtr1D[ 2*ki     ];
               OutPtr3D[ Idx_Out + OutStride[Pk][XYZ] ] = OutPtr1D[ 2*ki + 1 ];
            } // ki
         } // v
      } // ii, ji
   } // XYZ

   delete [] Input;
   delete [] Output;
   delete [] TDataX;
   delete [] TDataY;

   gamer_fftw::fft_free( Workspace );

} // FUNCTION : Int_Spectral



//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::InterpolationContext
// Description :  Constructor of InterpolationContext; perform parameter checks
//
// Parameter   :  nInput         : Size of input array
//                nGhostBoundary : Size of ghost boundary
//-------------------------------------------------------------------------------------------------------
InterpolationContext::InterpolationContext( size_t nInput, size_t nGhostBoundary ) :
   nInput         ( nInput ),
   nGhostBoundary ( nGhostBoundary ),
   nInterpolated  ( 2*(nInput - 2*nGhostBoundary) )
{

#  ifdef GAMER_DEBUG
   if ( nInput < 3 )
      Aux_Error( ERROR_INFO, "Input size %ld is smaller than minimum supported input size %d !!\n", nInput, 3 );

   if ( nGhostBoundary < 1 )
      Aux_Error( ERROR_INFO, "Ghost boundary must be greater than 0 !!\n" );

   if ( nInput <= 2*nGhostBoundary )
      Aux_Error( ERROR_INFO, "Input array of size %ld smaller than left and right ghost boundary of size %ld !!\n", nInput, nGhostBoundary );
#  endif

} // CONSTRUCTOR : InterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::ReadBinaryFile
// Description :  Read "size" doubles from binary file "filename" into memory at pointer "array"
//
// Note        :  Array must be preallocated, e.g. via array = (double*) malloc(sizeof(double) * size);
//
// Parameter   :  filename : Input filename
//                array    : Pointer to preallocated memory with sizeof(double) * size bytes
//                size     : Number of double-precision values to load
//
// Return      :  array
//-------------------------------------------------------------------------------------------------------
void InterpolationContext::ReadBinaryFile( const char* filename, double* array, const int size ) const
{

   FILE* file = fopen( filename, "rb" );
   if ( file == NULL ) {
      Aux_Error( ERROR_INFO, "Error opening interpolation file %s !!\n"
                 "        --> See https://github.com/gamer-project/gamer/wiki/ELBDM-Spectral-Interpolation#obtaining-interpolation-tables\n",
                 filename );
      return;
   }

// read the array from the binary file
   int num = fread( array, sizeof(double), size, file );

   if ( num != size ) {
//    fread failed
      if      ( ferror(file) )   // possibility 1
        Aux_Error( ERROR_INFO, "Error reading interpolation file %s (num %d, size %d) !!\n", filename, num, size );
      else if ( feof(file) )     // possibility 2
        Aux_Error( ERROR_INFO, "EOF found in interpolation file %s (num %d, size %d) !!\n", filename, num, size );
   }

// close the file
   fclose( file );

} // FUNCTION : ReadBinaryFile



//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::Preprocess
// Description :  Unwrap phase in ELBDM
//
// Parameter   :  input       : Input array
//                UnwrapPhase : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//-------------------------------------------------------------------------------------------------------
void InterpolationContext::Preprocess( real* input, const bool UnwrapPhase ) const
{

} // FUNCTION : Preprocess



//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::Postprocess
// Description :  Ensure monotonicity of interpolation result assuming nonconservative quartic interpolation
//
// Parameter   :  input           : Input array
//                output          : Output array
//                Monotonic       : Ensure that all interpolation results are monotonic
//                MonoCoeff       : Slope limiter coefficient for the option "Monotonic"
//                OppSign0thOrder : See Int_MinMod1D()
//-------------------------------------------------------------------------------------------------------
void InterpolationContext::Postprocess( const real* input, real* output, const bool Monotonic, const real MonoCoeff,
                                        const bool OppSign0thOrder) const
{

// ensure monotonicity
   if ( Monotonic )
   {
//    MonoCoeff/4
      const real MonoCoeff_4 = (real)0.25*MonoCoeff;
      real LSlopeDh_4, RSlopeDh_4, SlopeDh_4, Sign, CDataMax, CDataMin;

      size_t Idx_InC, Idx_InL1, Idx_InR1, Idx_Out, Tdy;

      for (size_t i=nGhostBoundary; i<nInput-nGhostBoundary; i++) {
         Idx_InC    = i;
         Idx_InL1   = i - 1;
         Idx_InR1   = i + 1;
         Idx_Out    = 2*( i - nGhostBoundary );
         Tdy        = 1;

         LSlopeDh_4 = input[Idx_InC ] - input[Idx_InL1];
         RSlopeDh_4 = input[Idx_InR1] - input[Idx_InC ];

         if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )
         {
            if ( LSlopeDh_4 > (real)0.0 )
            {
               CDataMax = input[Idx_InR1];
               CDataMin = input[Idx_InL1];
            }
            else
            {
               CDataMax = input[Idx_InL1];
               CDataMin = input[Idx_InR1];
            }

            if (  FMAX( output[Idx_Out], output[Idx_Out+Tdy] ) > CDataMax  ||
                  FMIN( output[Idx_Out], output[Idx_Out+Tdy] ) < CDataMin  ||
                  ( output[Idx_Out+Tdy] - output[Idx_Out] )*LSlopeDh_4 < (real)0.0  )
            {
               SlopeDh_4   = (real)0.125*( input[Idx_InR1] - input[Idx_InL1] );
               LSlopeDh_4 *= MonoCoeff_4;
               RSlopeDh_4 *= MonoCoeff_4;
               Sign        = SIGN( LSlopeDh_4 );

               SlopeDh_4  *= Sign;
               SlopeDh_4   = FMIN( Sign*LSlopeDh_4, SlopeDh_4 );
               SlopeDh_4   = FMIN( Sign*RSlopeDh_4, SlopeDh_4 );
               SlopeDh_4  *= Sign;

               output[ Idx_Out       ] = input[Idx_InC] - SlopeDh_4;
               output[ Idx_Out + Tdy ] = input[Idx_InC] + SlopeDh_4;
            }
         } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )

         else
         {
            output[ Idx_Out       ] = input[Idx_InC];
            output[ Idx_Out + Tdy ] = input[Idx_InC];
         } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 ) ... else
      } // if ( Monotonic[v] )

      if ( OppSign0thOrder  &&  input[Idx_InL1]*input[Idx_InR1] < (real)0.0 )
      {
         output[ Idx_Out       ] = input[Idx_InC];
         output[ Idx_Out + Tdy ] = input[Idx_InC];
      }
   } // if ( Monotonic )

} // FUNCTION : Postprocess



//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext
// Description :  Constructor of GramFEInterpolationContext, sets up FFTs and translation arrays in k-space
//                for fast interpolation, load Gram-Fourier extension table from disk
//
// Parameter   :  nInput         : Size of input function to be interpolated
//                nGhostBoundary : Size of ghost boundary
//                nExtension     : Size of extension domain
//                nDelta         : Size of boundary domain used for Gram-Fourier extension
//-------------------------------------------------------------------------------------------------------
GramFEInterpolationContext::GramFEInterpolationContext( size_t nInput, size_t nGhostBoundary, size_t nExtension, size_t nDelta )
   :  InterpolationContext( nInput, nGhostBoundary ),
      nExtension          ( nExtension ),
      nExtended           ( nInput + nExtension ),
      nExtendedPadded     ( nExtended/2 + 1 ),
      nDelta              ( nDelta )
{

// sanity checks
#  ifdef GAMER_DEBUG
   if ( nInput < nDelta )
      Aux_Error( ERROR_INFO, "Input array of size %ld smaller than Gram polynomial boundary of size %ld !!\n", nInput, nDelta );
#  endif

// load Gram-Fourier extension tables from file
   char filename[2*MAX_STRING];
   sprintf( filename, "%s/boundary2extension_tables/%ld_%ld_%d_%d.bin", SPEC_INT_TABLE_PATH, nExtension, nDelta, GRAMFE_GAMMA, GRAMFE_G );

   size_t  size  = nExtension*2*nDelta;
   double* table = (double*) malloc( size*sizeof(double) );
   ReadBinaryFile( filename, table, size );

// load extension table into GSL matrix
   extensionMatrix = gfei_gsl::matrix_alloc( nExtension, 2*nDelta );
   for (size_t i=0; i<nExtension; ++i) {
      for (size_t j=0; j<2*nDelta; ++j) {
            gfei_gsl::matrix_set( extensionMatrix, i, j, (gfei_gsl::gsl_real) table[ 2*nDelta*i + j ] );
      }
   }

// free table
   free( table );

// allocate memory for planning FFTs
   gfei_fftw::fft_complex* planner_array = (gfei_fftw::fft_complex*) gfei_fftw::fft_malloc( nExtended*sizeof(gfei_fftw::fft_complex) );

// plan FFT of function to be interpolated
   r2cPlan = gfei_fftw_create_1d_r2c_plan( nExtended, planner_array, FFTW_STARTUP_MEASURE );
   c2rPlan = gfei_fftw_create_1d_c2r_plan( nExtended, planner_array, FFTW_STARTUP_MEASURE );

// free memory
   gfei_fftw::fft_free( planner_array );

// set up translation arrays
   translationCoeffL = (gfei_fftw::fft_complex*) malloc( nExtendedPadded*sizeof(gfei_fftw::fft_complex) );
   translationCoeffR = (gfei_fftw::fft_complex*) malloc( nExtendedPadded*sizeof(gfei_fftw::fft_complex) );

   const gfei_fftw::fft_real Norm = 1.0 / ( (gfei_fftw::fft_real) nExtended );
   gfei_fftw::fft_real K;
   for (size_t i=0; i<nExtendedPadded; i++)
   {
      K                          = 2.0*M_PI/nExtended*i;
      c_re(translationCoeffL[i]) = cos(-K/4.0)*Norm;
      c_im(translationCoeffL[i]) = sin(-K/4.0)*Norm;
      c_re(translationCoeffR[i]) = cos(+K/4.0)*Norm;
      c_im(translationCoeffR[i]) = sin(+K/4.0)*Norm;
   }

} // CONSTRUCTOR : GramFEInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext:~GramFEInterpolationContext
// Description :  Destructor of GramFEInterpolationContext, frees memory
//-------------------------------------------------------------------------------------------------------
GramFEInterpolationContext::~GramFEInterpolationContext()
{

// free interpolation matrix
   gfei_gsl::matrix_free( extensionMatrix );

// free FFTW plans
   gfei_fftw::destroy_real_plan_1d( r2cPlan );
   gfei_fftw::destroy_real_plan_1d( c2rPlan );

// free arrays with translation coefficients
   free( translationCoeffL );
   free( translationCoeffR );

} // DESTRUCTOR : GramFEInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext:GetWorkspaceSize
// Description :  Return size of workspace in bytes
//
// Parameter   :  nInput         : Size of input array
//                nGhostBoundary : Size of ghost boundary
//
// Return      :  Workspace size in bytes as unsigned integer
//-------------------------------------------------------------------------------------------------------
size_t GramFEInterpolationContext::GetWorkspaceSize() const
{

   size_t total_size = 0;
   total_size += 2*nDelta;          // left and right boundary for Gram polynomials
   total_size += 2*nExtendedPadded; // extension region
   total_size += 2*nExtendedPadded; // left interpolated values
   total_size += 2*nExtendedPadded; // right interpolated values

   return total_size*sizeof( gfei_fftw::fft_real );

} // FUNCTION : GetWorkspaceSize



//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext:InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2*(nInput - nGhostBoundary) in output array
//
// Note        :  This method might suffer from large round-off errors because the extension matrix has a high-condition number.
//                Tests of using the Gram-Fourier extension method with precomputed extension matrices indicate that the matrix has a condition number of
//                roughly 1e11. This implies that the extension should be computed with quadruple precision in order to be accurate.
//                However, tests of this interpolation function in single precision indicate that the roundoff error may be acceptable
//                for interpolation because it does not accumulate over time. Use with care.
//
// Parameter   :  input     : Real input  array of size nInput
//                output    : Real output array of size 2*(nInput - nGhostBoundary)
//                workspace : Workspace with size GetWorkspaceSize in bytes that needs to be allocated by user using fftw_malloc
//-------------------------------------------------------------------------------------------------------
void GramFEInterpolationContext::InterpolateReal( const real* input, real *output, char* workspace ) const
{

// define arrays in workspace
   gfei_fftw::fft_real* boundary      = (gfei_fftw::fft_real*) workspace;
   gfei_fftw::fft_real* inputExtended = boundary      + 2*nDelta;          // extension region
   gfei_fftw::fft_real* outputL       = inputExtended + 2*nExtendedPadded; // left output
   gfei_fftw::fft_real* outputR       = outputL       + 2*nExtendedPadded; // right output

// convert input to FFT precision
   for (size_t i=0; i<nInput; ++i) {
      inputExtended[i] = (gfei_fftw::fft_real) input[i];
   }

// write boundaries of input array to boundary array
   for (size_t i=0; i<nDelta; ++i) {
      boundary[ i          ] = (gfei_fftw::fft_real) input[                   i ];
      boundary[ i + nDelta ] = (gfei_fftw::fft_real) input[ nInput - nDelta + i ];
   }

// compute extension using GSL BLAS real matrix-vector multiplication
   gfei_gsl::vector_const_view Boundary_view = gfei_gsl::vector_const_view_array( boundary, 2*nDelta );
   gfei_gsl::vector_view Extension_view      = gfei_gsl::vector_view_array      ( inputExtended + nInput, nExtension );
   gfei_gsl::blas_sgemv( CblasNoTrans, 1.0, extensionMatrix, &Boundary_view.vector, 0.0, &Extension_view.vector );

// compute forward FFT
   gfei_fftw_r2c( r2cPlan, inputExtended );

// set up arrays with complex multiplication
   std::complex<gfei_fftw::fft_real>* InputK   = ( std::complex<gfei_fftw::fft_real>* ) inputExtended;
   std::complex<gfei_fftw::fft_real>* OutputLK = ( std::complex<gfei_fftw::fft_real>* ) outputL;
   std::complex<gfei_fftw::fft_real>* OutputRK = ( std::complex<gfei_fftw::fft_real>* ) outputR;
   std::complex<gfei_fftw::fft_real>* TransL   = ( std::complex<gfei_fftw::fft_real>* ) translationCoeffL;
   std::complex<gfei_fftw::fft_real>* TransR   = ( std::complex<gfei_fftw::fft_real>* ) translationCoeffR;

// interpolate by shifting samples to the left by 0.25 dh
   for (size_t i=0; i<nExtendedPadded; i++)
   {
      OutputLK[i] = InputK[i]*TransL[i];
      OutputRK[i] = InputK[i]*TransR[i];
   }

// compute backward FFT
   gfei_fftw_c2r( c2rPlan, OutputLK );

// compute backward FFT
   gfei_fftw_c2r( c2rPlan, OutputRK );

// store output data
   for (size_t i=nGhostBoundary; i<nInput-nGhostBoundary; i++)
   {
      output[ 2*(i - nGhostBoundary)     ] = outputL[i];
      output[ 2*(i - nGhostBoundary) + 1 ] = outputR[i];
   }

} // FUNCTION : InterpolateReal



//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::PrecomputedInterpolationContext
// Description :  Constructor of PrecomputedInterpolationContext, loads interpolation table into memory
//-------------------------------------------------------------------------------------------------------
PrecomputedInterpolationContext::PrecomputedInterpolationContext( size_t nInput, size_t nGhostBoundary )
   :  InterpolationContext( nInput, nGhostBoundary )
{

   char filename[2*MAX_STRING];
   sprintf( filename, "%s/interpolation_tables/N=%ld.bin", SPEC_INT_TABLE_PATH, nInput );

   size_t size   = 2*(nInput - 2)*nInput;
   double* table = (double*) malloc( size*sizeof(double) );
   ReadBinaryFile( filename, table, size );

// write table to GSL matrix in desired precision
   interpolationMatrix = pi_gsl::matrix_alloc( nInterpolated, nInput );
   for (size_t i=0; i<nInterpolated; ++i) {
      for (size_t j=0; j<nInput; ++j) {
         size_t index = ( i + 2*(nGhostBoundary - 1) )*nInput + j;
         pi_gsl::matrix_set( interpolationMatrix, i, j, (pi_gsl::gsl_real) table[index] );
      }
   }

// free table array
   free( table );

} // CONSTRUCTOR : PrecomputedInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::~PrecomputedInterpolationContext
// Description :  Destructor of PrecomputedInterpolationContext, frees memory
//-------------------------------------------------------------------------------------------------------
PrecomputedInterpolationContext::~PrecomputedInterpolationContext()
{

// free interpolation matrix
   pi_gsl::matrix_free( interpolationMatrix );

} // DESTRUCTOR : PrecomputedInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::GetWorkspaceSize
// Description :  Return size of input and output arrays in pi_gsl::gsl_real precision for matrix multiplication
//
// Return      :  Returns the workspace size in bytes
//-------------------------------------------------------------------------------------------------------
size_t PrecomputedInterpolationContext::GetWorkspaceSize() const
{

   size_t total_size = 0;
   total_size += nInput;        // input
   total_size += nInterpolated; // output

   return total_size*sizeof( pi_gsl::gsl_real );

} // FUNCTION : GetWorkspaceSize



//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
// Parameter   :  input     : Real input  array of size nInput
//                output    : Real output array of size 2 * (nInput - nGhostBoundary)
//                workspace : Workspace containing input and output in double precision pi_gsl::gsl_real
//-------------------------------------------------------------------------------------------------------
void PrecomputedInterpolationContext::InterpolateReal( const real *input, real *output, char* workspace ) const
{

// define arrays in workspace
   pi_gsl::gsl_real* pi_gsl_input  = (pi_gsl::gsl_real*) workspace;
   pi_gsl::gsl_real* pi_gsl_output = pi_gsl_input + nInput;

// convert input to matrix multiplication precision
   for (size_t i=0; i<nInput; ++i) {
      pi_gsl_input[i] = (pi_gsl::gsl_real) input[i];
   }

   pi_gsl::vector_const_view in = pi_gsl::vector_const_view_array( pi_gsl_input, nInput );
   pi_gsl::vector_view out      = pi_gsl::vector_view_array      ( pi_gsl_output, nInterpolated );

   pi_gsl::blas_sgemv( CblasNoTrans, 1.0, interpolationMatrix, &in.vector, 0.0, &out.vector );


   for (size_t i=0; i<nInterpolated; ++i) {
      output[i] = (real) pi_gsl_output[i];
   }

} // FUNCTION : InterpolateReal



// fixed interpolation constants
const real QuarticInterpolationContext::QuarticR[5] = { +35.0/2048.0, -252.0/2048.0, +1890.0/2048.0, +420.0/2048.0, -45.0/2048.0 };
const real QuarticInterpolationContext::QuarticL[5] = { -45.0/2048.0, +420.0/2048.0, +1890.0/2048.0, -252.0/2048.0, +35.0/2048.0 };

//-------------------------------------------------------------------------------------------------------
// Function    :  QuarticInterpolationContext::QuarticInterpolationContext
// Description :  Constructor of QuarticInterpolationContext
//-------------------------------------------------------------------------------------------------------
QuarticInterpolationContext::QuarticInterpolationContext( size_t nInput, size_t nGhostBoundary )
   : InterpolationContext( nInput, nGhostBoundary )
{

#  ifdef GAMER_DEBUG
   if ( nGhostBoundary != 2 )
      Aux_Error( ERROR_INFO, "QuarticInterpolationContext requires nGhostBoundary = 2 !!\n" );
#  endif

} // CONSTRUCTOR : QuarticInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  QuarticInterpolationContext::GetWorkspaceSize
// Description :  Return 0 since QuarticInterpolationContext does not require a workspace
//
// Return      :  Returns the unsigned integer 0
//-------------------------------------------------------------------------------------------------------
size_t QuarticInterpolationContext::GetWorkspaceSize() const
{

   return 0;

} // FUNCTION : GetWorkspaceSize



//-------------------------------------------------------------------------------------------------------
// Function    :  QuadraticInterpolationContext::InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2*(nInput - nGhostBoundary) in output array
//
// Parameter   :  input     : Real input  array of size nInput
//                output    : Real output array of size 2*(nInput - nGhostBoundary)
//                workspace : Useless for QuadraticInterpolationContext
//-------------------------------------------------------------------------------------------------------
void QuarticInterpolationContext::InterpolateReal( const real *input, real *output, char* workspace ) const
{

   for (size_t i=nGhostBoundary; i<nInput-nGhostBoundary; ++i) {
      real r=0, l=0;
      for (size_t j=0; j<5; ++j) {
         l += QuarticInterpolationContext::QuarticL[j]*input[ i - 2 + j ];
         r += QuarticInterpolationContext::QuarticR[j]*input[ i - 2 + j ];
      }
      output[ (i-nGhostBoundary)*2     ] = l;
      output[ (i-nGhostBoundary)*2 + 1 ] = r;
   }

} // FUNCTION : InterpolateReal



const real CQuarticInterpolationContext::CQuartic[5] = { +3.0/128.0, -22.0/128.0, +128.0/128.0, +22.0/128.0, -3.0/128.0 };

//-------------------------------------------------------------------------------------------------------
// Function    :  CQuarticInterpolationContext::CQuarticInterpolationContext
// Description :  Constructor of CQuarticInterpolationContext
//-------------------------------------------------------------------------------------------------------
CQuarticInterpolationContext::CQuarticInterpolationContext( size_t nInput, size_t nGhostBoundary )
   : InterpolationContext( nInput, nGhostBoundary )
{

#  ifdef GAMER_DEBUG
   if ( nGhostBoundary != 2 )
      Aux_Error( ERROR_INFO, "CQuarticInterpolationContext requires nGhostBoundary = 2 !!\n" );
#  endif

} // CONSTRUCTOR : CQuarticInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Function    :  CQuarticInterpolationContext::GetWorkspaceSize
// Description :  Return 0 since CQuarticInterpolationContext does not require a workspace
//
// Return      :  Returns the unsigned integer 0
//-------------------------------------------------------------------------------------------------------
size_t CQuarticInterpolationContext::GetWorkspaceSize() const
{

   return 0;

} // FUNCTION : GetWorkspaceSize



//-------------------------------------------------------------------------------------------------------
// Function    :  QuadraticInterpolationContext::InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2*(nInput - nGhostBoundary) in output array
//
// Parameter   :  input     : Real input  array of size nInput
//                output    : Real output array of size 2*(nInput - nGhostBoundary)
//                workspace : Useless for QuadraticInterpolationContext
//-------------------------------------------------------------------------------------------------------
void CQuarticInterpolationContext::InterpolateReal( const real *input, real *output, char* workspace ) const
{

   for (size_t i=nGhostBoundary; i<nInput-nGhostBoundary; ++i) {
      const real slope = CQuarticInterpolationContext::CQuartic[0]*input[ i - 2 ]
                       + CQuarticInterpolationContext::CQuartic[1]*input[ i - 1 ]
                       + CQuarticInterpolationContext::CQuartic[3]*input[ i + 1 ]
                       + CQuarticInterpolationContext::CQuartic[4]*input[ i + 2 ];

      output[ (i-nGhostBoundary)*2     ] = input[i] - slope;
      output[ (i-nGhostBoundary)*2 + 1 ] = input[i] + slope;
   }

} // FUNCTION : InterpolateReal



const size_t NFast = 24;
const size_t fastNExtended[NFast] = { 16, 18, 24, 32, 36, 40, 48, 54, 60, 64, 72, 80, 90, 96, 100, 108, 112, 120, 128, 135, 256, 288, 512, 540 };

void InterpolationHandler::AddInterpolationContext( size_t nInput, size_t nGhostBoundary )
{

// ensure thread safety when adding new interpolation contexts
// only one thread in OMP enviroment may add a new interpolation context
#  pragma omp critical
   if ( contexts.find(nInput) == contexts.end() ) {
//    for small N <= 32 pick precomputed interpolation with cost of N^2
      if ( nInput <= 32 ) {
            contexts.emplace( nInput, new PrecomputedInterpolationContext(nInput, nGhostBoundary) );
//    for large N >  32 use Gram-Fourier extension scheme with cost of N log(N)
      } else {
         size_t nExtension = 32, nDelta = 14;
//
         const size_t minimumExtensionSize = 24;
         const size_t maximumExtensionSize = 36;
         for (size_t i=0; i<NFast; ++i) {
            if ( nInput + minimumExtensionSize < fastNExtended[i]  &&  nInput + maximumExtensionSize >= fastNExtended[i] ) {
               nExtension = fastNExtended[i] - nInput;
               break;
            }
         }
         contexts.emplace( nInput, new GramFEInterpolationContext(nInput, nGhostBoundary, nExtension, nDelta) );
      }
   }

} // FUNCTION : AddInterpolationContext



size_t InterpolationHandler::GetWorkspaceSize( size_t nInput, size_t nGhostBoundary ) const
{

   return contexts.at( nInput )->GetWorkspaceSize();

} // FUNCTION : GetWorkspaceSize



void InterpolationHandler::InterpolateReal( real* input, real* output, size_t nInput, size_t nGhostBoundary, char* workspace,
                                            const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder ) const
{

   auto context = contexts.at(nInput);

// unwrap phase
// context->Preprocess( input, UnwrapPhase );

// interpolate
   context->InterpolateReal( input, output, workspace );

// ensure monotonicity
// context->Postprocess( input, output, Monotonic, MonoCoeff, OppSign0thOrder );

} // FUNCTION : InterpolateReal



#endif // #endif SUPPORT_SPECTRAL_INT
