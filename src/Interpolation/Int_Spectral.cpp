#include "GAMER.h"


#ifdef SUPPORT_SPECTRAL_INTERPOLATION

#include <complex>
#include <memory>
#include "GramFE_Interpolation.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Int_Spectral
// Description :  Perform spatial interpolation based on the Gram-Fourier extension method
//
// Note        :  1. The interpolation is spectrally accurate
//                2. The interpolation result is neither  nor monotonic
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
//                Monotonic       : Ensure that all interpolation results are monotonic
//                MonoCoeff       : Slope limiter coefficient for the option "Monotonic"
//                OppSign0thOrder : See Int_MinMod1D()
//-------------------------------------------------------------------------------------------------------
void Int_Spectral(  real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                    real FData[], const int FSize[3], const int FStart[3], const int NComp,
                    const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder )
{

// interpolation-scheme-dependent parameters
// ===============================================================================
// number of coarse-grid ghost zone
   const int CGhost    = 2;
// ===============================================================================


// allocate memory to store largest input and output column
   size_t maxSize       = MAX(MAX(CRange[0], CRange[1]), CRange[2]);
   real* Input, *Output;
   Input  = (real*) malloc( (maxSize +  2 * CGhost) * sizeof(real) );
   Output = (real*) malloc(  2 * maxSize            * sizeof(real) );

// add interpolation contexts (make sure that everything is thread-safe in AddInterpolationContext)
   for (size_t i = 0; i < 3; ++i) {
      INTERPOLATION_HANDLER.AddInterpolationContext( CRange[i] +  2 * CGhost, CGhost );
   }

// determine workspace size
   size_t workspaceSize = 0;
   for (size_t i = 0; i < 3; ++i) {
      workspaceSize = MAX(workspaceSize, INTERPOLATION_HANDLER.GetWorkspaceSize( CRange[i] +  2 * CGhost, CGhost ));
   }

// allocate workspace using fftw_malloc for contiguous memory used in FFT
   char* workspace = (char*) gamer_fftw::fft_malloc(workspaceSize);


// index stride of the coarse-grid input array
   const int Cdx    = 1;
   const int Cdy    = Cdx*CSize[0];
   const int Cdz    = Cdy*CSize[1];

// index stride of the temporary arrays storing the data after x and y interpolations
   const int Tdx    = 1;
   const int Tdy    = Tdx* CRange[0]*2;
   const int TdzX   = Tdy*(CRange[1]+2*CGhost);    // array after x interpolation
   const int TdzY   = Tdy* CRange[1]*2;            // array after y interpolation

// index stride of the fine-grid output array
   const int Fdx    = 1;
   const int Fdy    = Fdx*FSize[0];
   const int Fdz    = Fdy*FSize[1];

// index stride of different components
   const int CDisp  = CSize[0]*CSize[1]*CSize[2];
   const int FDisp  = FSize[0]*FSize[1]*FSize[2];


   real *CPtr   = CData;
   real *FPtr   = FData;
   real *TDataX = new real [ (CRange[2]+2*CGhost)*TdzX ];   // temporary array after x interpolation
   real *TDataY = new real [ (CRange[2]+2*CGhost)*TdzY ];   // temporary array after y interpolation

   int Idx_InL, Idx_InC, Idx_InR, Idx_Out;
   int nInput;

   for (int v=0; v<NComp; v++)
   {
//    interpolation along x direction
      for (int In_z=CStart[2]-CGhost, Out_z=0;  In_z<CStart[2]+CRange[2]+CGhost;  In_z++, Out_z++)
      for (int In_y=CStart[1]-CGhost, Out_y=0;  In_y<CStart[1]+CRange[1]+CGhost;  In_y++, Out_y++)
      {
         nInput = CRange[0] + 2 * CGhost;

//       fill FFT array with one column of input data (including ghost zones) in x-direction
         for (int In_x=CStart[0]-CGhost, Out_x=0;  In_x<CStart[0]-CGhost+nInput;  In_x++, Out_x++)
         {
            Idx_InC                 = In_z*Cdz  +  In_y*Cdy +  In_x*Cdx;
            Input[Out_x]            = CPtr[Idx_InC];
         } // i

//       interpolate data
         INTERPOLATION_HANDLER.InterpolateReal(Input, Output, nInput, CGhost, workspace, UnwrapPhase, Monotonic[v], MonoCoeff, OppSign0thOrder);

//       write result of Fourier interpolation (excluding ghost zones) to temporary array
         for (int In_x=0, Out_x=0; In_x<CRange[0]; In_x++, Out_x+=2)
         {
            Idx_Out                 = Out_z*TdzX + Out_y*Tdy + Out_x*Tdx;
            TDataX[ Idx_Out       ] = Output[2 * In_x    ];
            TDataX[ Idx_Out + Tdx ] = Output[2 * In_x + 1];
         } // i
      } // k,j

//    interpolation along y direction
      for (int InOut_z=0;             InOut_z<CRange[2]+2*CGhost;  InOut_z++)
      for (int InOut_x=0;             InOut_x<2*CRange[0];         InOut_x++)
      {

         nInput = CRange[1] + 2 * CGhost;
//       fill FFT array with one column of input data (including ghost boundary) in y-direction
         for (int In_y=0, Out_y=0;    In_y < nInput;   In_y++, Out_y++) {
            Idx_InC                 = InOut_z*TdzX +  In_y*Tdy + InOut_x*Tdx;
            Input[Out_y]            = TDataX[Idx_InC];
         } // j

//       interpolate data
         INTERPOLATION_HANDLER.InterpolateReal(Input, Output, nInput, CGhost, workspace, UnwrapPhase, Monotonic[v], MonoCoeff, OppSign0thOrder);


//       write result of Fourier interpolation to temporary array
         for (int In_y=0, Out_y=0;  In_y<CRange[1];   In_y++, Out_y+=2) {
            Idx_Out = InOut_z*TdzY + Out_y*Tdy + InOut_x*Tdx;
            TDataY[ Idx_Out       ] = Output[2 * In_y    ];
            TDataY[ Idx_Out + Tdy ] = Output[2 * In_y + 1];
         } // j
      } // k,i

//    interpolation along z direction
      for (int In_y=0,      Out_y=FStart[1];  In_y<2*CRange[1];       In_y++, Out_y++)
      for (int In_x=0,      Out_x=FStart[0];  In_x<2*CRange[0];       In_x++, Out_x++)
      {

         nInput = CRange[2] + 2 * CGhost;
//       fill FFT array with one column of input data in z-direction
         for (int In_z=0, Out_z=0;  In_z<nInput;  In_z++, Out_z++)
         {
            Idx_InC                 = In_z*TdzY +  In_y*Tdy +  In_x*Tdx;
            Input[Out_z]            = TDataY[Idx_InC];
         }

//       interpolate data
         INTERPOLATION_HANDLER.InterpolateReal(Input, Output, nInput, CGhost, workspace, UnwrapPhase, Monotonic[v], MonoCoeff, OppSign0thOrder);

         for (int In_z=0, Out_z=FStart[2];  In_z<CRange[2];  In_z++, Out_z+=2)
         {
            Idx_Out                 = Out_z*Fdz  + Out_y*Fdy + Out_x*Fdx;
            FPtr[ Idx_Out       ]   = Output[ 2 * In_z     ];
            FPtr[ Idx_Out + Fdz ]   = Output[ 2 * In_z + 1 ];
         }
      } // k,j,i

      CPtr += CDisp;
      FPtr += FDisp;
   } // for (int v=0; v<NComp; v++)

   delete [] TDataX;
   delete [] TDataY;

   free( Input  );
   free( Output );
   gamer_fftw::fft_free(workspace);

} // FUNCTION : Int_Spectral


//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::InterpolationContext
// Description :  Constructor of InterpolationContext; perform parameter checks
//
// Parameter   :  nInput             : Size of input array
//                nGhostBoundary     : Size of ghost boundary
//
//-------------------------------------------------------------------------------------------------------
InterpolationContext::InterpolationContext(size_t nInput, size_t nGhostBoundary) :
   nInput         (nInput),
   nGhostBoundary (nGhostBoundary),
   nInterpolated  (2 * (nInput - 2 * nGhostBoundary))
{
#  ifdef GAMER_DEBUG
   if ( nInput <  3 )
   {
      Aux_Error(ERROR_INFO, "Input size %ld is smaller than minimum supported input size %d!!\n", nInput, 3);
   }
   if ( nGhostBoundary <  1 )
   {
      Aux_Error(ERROR_INFO, "Ghost boundary must be greater than 0!!\n");
   }
   if ( nInput <= 2 * nGhostBoundary )
   {
      Aux_Error(ERROR_INFO, "Input array of size %ld smaller than left and right ghost boundary of size %ld!!\n", nInput, nGhostBoundary);
   }
#  endif
} // CONSTRUCTOR : InterpolationContext

//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::ReadBinaryFile
// Description :  Read "size" doubles from binary file "filename" into memory at pointer "array"
//
// Note        :  Array must be preallocated, e.g. via array = (double*) malloc(sizeof(double) * size);
//
// Parameter   :  filename    : Input filename
//                array       : Pointer to preallocated memory with sizeof(double) * size bytes
//                size        : Number of double precision values to load
//
//-------------------------------------------------------------------------------------------------------
void InterpolationContext::ReadBinaryFile(const char* filename, double* array, int size) const
{
   FILE* file = fopen(filename, "rb");
   if (file == NULL) {
      Aux_Error(ERROR_INFO, "Error opening interpolation file %s.\n", filename);
      return;
   }

// read the array from the binary file
   int num = fread(array, sizeof(double), size, file);

   if ( num != size ) {
       /* fread failed */
      if ( ferror(file) )    /* possibility 1 */
        Aux_Error(ERROR_INFO, "Error reading interpolation file %s.\n", filename);
      else if ( feof(file))  /* possibility 2 */
        Aux_Error(ERROR_INFO, "EOF found in interpolation file %s.\n", filename);
   }

// close the file
   fclose(file);
} // FUNCTION : ReadBinaryFile


//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolationContext::Preprocess
// Description :  Unwrap phase in ELBDM
//
// Parameter   :  input           : Input array
//                UnwrapPhase     : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//
//-------------------------------------------------------------------------------------------------------
void InterpolationContext::Preprocess(real* input, const bool UnwrapPhase) const {
//    unwrap phase
#     if ( MODEL == ELBDM )
      if ( UnwrapPhase )
      {
         size_t Idx_InC, Idx_InL1;
         real* CPtr = input;
         for (size_t i = 1;  i<nInput;  i++)
         {
            Idx_InC       = i;
            Idx_InL1      = i - 1;
            CPtr[Idx_InC] = ELBDM_UnwrapPhase( CPtr[Idx_InL1], CPtr[Idx_InC] );
         }
      }
#     endif
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
//
//-------------------------------------------------------------------------------------------------------
void InterpolationContext::Postprocess(const real* input, real* output, const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder) const {
// ensure monotonicity
   if ( Monotonic )
   {
//    MonoCoeff/4
      const real MonoCoeff_4 = (real)0.25*MonoCoeff;
      real LSlopeDh_4, RSlopeDh_4, SlopeDh_4, Sign, CDataMax, CDataMin;

      size_t Idx_InC, Idx_InL1, Idx_InR1, Idx_Out, Tdy;

      for (size_t i = nGhostBoundary;  i < nInput - nGhostBoundary;  i++) {
         Idx_InC    = i;
         Idx_InL1   = i - 1;
         Idx_InR1   = i + 1;
         Idx_Out    = 2 * (i - nGhostBoundary);
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

               SlopeDh_4 *= Sign;
               SlopeDh_4  = FMIN( Sign*LSlopeDh_4, SlopeDh_4 );
               SlopeDh_4  = FMIN( Sign*RSlopeDh_4, SlopeDh_4 );
               SlopeDh_4 *= Sign;

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
   }
} // FUNCTION : Postprocess

//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext
// Description :  Constructor of GramFEInterpolationContext, sets up FFTs and translation arrays in k-space for fast interpolation, load Gram-Fourier extension table from disk
//
// Parameter   :  nInput         : Size of input function to be interpolated
//                nGhostBoundary : Size of ghost boundary
//                nExtension     : Size of extension domain
//                nDelta         : Size of boundary domain used for Gram-Fourier extension
//
//-------------------------------------------------------------------------------------------------------
GramFEInterpolationContext::GramFEInterpolationContext(size_t nInput, size_t nGhostBoundary, size_t nExtension, size_t nDelta)
   :  InterpolationContext(nInput, nGhostBoundary),
      nExtension          (nExtension),
      nExtended           (nInput + nExtension),
      nExtendedPadded     (nExtended / 2 + 1),
      nDelta              (nDelta)
{
// sanity checks
#  ifdef GAMER_DEBUG
   if ( nInput < nDelta )
   {
      Aux_Error(ERROR_INFO, "Input array of size %ld smaller than Gram polynomial boundary of size %ld!!\n", nInput, nDelta);
   }
#  endif // # ifdef GAMER_DEBUG

// load Gram-Fourier extension tables from file
   char filename[2 * MAX_STRING];
   sprintf(filename, "%s/boundary2extension_tables/%ld_%ld_%d_%d.bin", INT_TABLE_PATH, nExtension, nDelta, 150, 63);

   size_t  size  = nExtension * 2 * nDelta;
   double* table = (double*) malloc(size * sizeof(double));
   ReadBinaryFile(filename, table, size);


// load extension table into GSL matrix
   extensionMatrix = gfei_gsl::matrix_alloc(nExtension, 2 * nDelta);
   for (size_t i = 0; i < nExtension; ++i) {
      for (size_t j = 0; j < 2 * nDelta; ++j) {
            gfei_gsl::matrix_set(extensionMatrix, i, j, (gfei_gsl::gsl_real) table[2 * nDelta * i + j]);
      }
   }

// free table
   free(table);

// allocate memory for planning FFTs
   gfei_fftw::fft_complex* planner_array = (gfei_fftw::fft_complex*) gfei_fftw::fft_malloc(nExtended * sizeof(gfei_fftw::fft_complex));

// plan FFT of function to be interpolated
   r2cPlan  = gfei_fftw_create_1d_r2c_plan(nExtended, planner_array, FFTW_STARTUP_MEASURE);
   c2rPlan  = gfei_fftw_create_1d_c2r_plan(nExtended, planner_array, FFTW_STARTUP_MEASURE);

// free memory
   gfei_fftw::fft_free(planner_array);

// set up translation arrays
   translationCoeffL = ( gfei_fftw::fft_complex* ) malloc( nExtendedPadded * sizeof(gfei_fftw::fft_complex) );
   translationCoeffR = ( gfei_fftw::fft_complex* ) malloc( nExtendedPadded * sizeof(gfei_fftw::fft_complex) );

   const gfei_fftw::fft_real Norm = 1.0 / ( (gfei_fftw::fft_real) nExtended );
   gfei_fftw::fft_real  K;
   for (size_t i = 0; i < nExtendedPadded; i++)
   {
      K                    = 2.0 * M_PI / nExtended * i;
      c_re(translationCoeffL[i]) = cos(-K/4.0)* Norm;
      c_im(translationCoeffL[i]) = sin(-K/4.0)* Norm;
      c_re(translationCoeffR[i]) = cos(+K/4.0)* Norm;
      c_im(translationCoeffR[i]) = sin(+K/4.0)* Norm;
   }

} // CONSTRUCTOR : GramFEInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext:~GramFEInterpolationContext
// Description :  Destructor of GramFEInterpolationContext, frees memory
//
//-------------------------------------------------------------------------------------------------------
GramFEInterpolationContext::~GramFEInterpolationContext()
{
// free interpolation matrix
   gfei_gsl::matrix_free(extensionMatrix);
// free FFTW plans
   gfei_fftw::destroy_real_plan_1d(r2cPlan);
   gfei_fftw::destroy_real_plan_1d(c2rPlan);
// free arrays with translation coefficients
   free(translationCoeffL);
   free(translationCoeffR);
} // DESTRUCTOR : GramFEInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext:GetWorkspaceSize
// Description :  Allocate block of memory using fftw_malloc
//
// Parameter   :  nInput         : Size of input array
//                nGhostBoundary : Size of ghost boundary
//
//-------------------------------------------------------------------------------------------------------
size_t GramFEInterpolationContext::GetWorkspaceSize() const {
   size_t total_size = 0;
   total_size += 2  * nDelta;          // left and right boundary for Gram polynomials
   total_size += 2  * nExtendedPadded; // extension region
   total_size += 2  * nExtendedPadded; // left interpolated values
   total_size += 2  * nExtendedPadded; // right interpolated values
   return total_size * sizeof(gfei_fftw::fft_real);
} // FUNCTION : GetWorkspaceSize


//-------------------------------------------------------------------------------------------------------
// Function    :  GramFEInterpolationContext:InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
// Parameter   :  input          : Real input  array of size nInput
//                output         : Real output array of size 2 * (nInput - nGhostBoundary)
//                workspace      : Workspace with size GetWorkspaceSize in bytes that needs to be allocated by user using fftw_malloc
//
//-------------------------------------------------------------------------------------------------------
void GramFEInterpolationContext::InterpolateReal(const real* input, real *output, char* workspace) const
{
// define arrays in workspace
   gfei_fftw::fft_real* boundary      = (gfei_fftw::fft_real*) workspace;
   gfei_fftw::fft_real* inputExtended = boundary      + 2 * nDelta;          // extension region
   gfei_fftw::fft_real* outputL       = inputExtended + 2 * nExtendedPadded; // left output
   gfei_fftw::fft_real* outputR       = outputL       + 2 * nExtendedPadded; // right output


// convert input to FFT precision
   for (size_t i = 0; i < nInput; ++i) {
      inputExtended[i]     = (gfei_fftw::fft_real) input[ i];
   }

// write boundaries of input array to boundary array
   for (size_t i = 0; i < nDelta; ++i) {
      boundary[i         ] = (gfei_fftw::fft_real) input[                  i];
      boundary[i + nDelta] = (gfei_fftw::fft_real) input[nInput - nDelta + i];
   }

// compute extension using GSL BLAS real matrix-vector multiplication
   gfei_gsl::vector_const_view Boundary_view = gfei_gsl::vector_const_view_array(boundary, 2 * nDelta);
   gfei_gsl::vector_view Extension_view      = gfei_gsl::vector_view_array      (inputExtended + nInput, nExtension);
   gfei_gsl::blas_sgemv(CblasNoTrans, 1.0, extensionMatrix, &Boundary_view.vector, 0.0, &Extension_view.vector);

// compute forward FFT
   gfei_fftw_r2c   (r2cPlan, inputExtended);

// set up arrays with complex multiplication
   std::complex<gfei_fftw::fft_real>* InputK   = (std::complex<gfei_fftw::fft_real>*) inputExtended;
   std::complex<gfei_fftw::fft_real>* OutputLK = (std::complex<gfei_fftw::fft_real>*) outputL;
   std::complex<gfei_fftw::fft_real>* OutputRK = (std::complex<gfei_fftw::fft_real>*) outputR;
   std::complex<gfei_fftw::fft_real>* TransL   = (std::complex<gfei_fftw::fft_real>*) translationCoeffL;
   std::complex<gfei_fftw::fft_real>* TransR   = (std::complex<gfei_fftw::fft_real>*) translationCoeffR;

// interpolate by shifting samples to the left by 0.25 dh
   for (size_t i = 0; i < nExtendedPadded; i++)
   {
      OutputLK[i] = InputK[i] * TransL[i];
      OutputRK[i] = InputK[i] * TransR[i];
   }

// compute backward FFT
   gfei_fftw_c2r( c2rPlan, OutputLK );

// compute backward FFT
   gfei_fftw_c2r( c2rPlan, OutputRK );

// store output data
   for (size_t i = nGhostBoundary; i < nInput - nGhostBoundary; i++)
   {
      output[2 * (i - nGhostBoundary)    ] = outputL[i];
      output[2 * (i - nGhostBoundary) + 1] = outputR[i];
   }
} // FUNCTION : InterpolateReal


//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::PrecomputedInterpolationContext
// Description :  Constructor of PrecomputedInterpolationContext, loads interpolation table into memory
//
//-------------------------------------------------------------------------------------------------------
PrecomputedInterpolationContext::PrecomputedInterpolationContext(size_t nInput, size_t nGhostBoundary)
   :  InterpolationContext(nInput, nGhostBoundary)
{
   char filename[2 * MAX_STRING];
   sprintf(filename, "%s/interpolation_tables/N=%ld.bin", INT_TABLE_PATH, nInput);

   size_t size   = 2 * (nInput - 2) * nInput;
   double* table = (double*) malloc(size * sizeof(double));
   ReadBinaryFile(filename, table, size);

// write table to GSL matrix in desired precision
   interpolationMatrix = pi_gsl::matrix_alloc(nInterpolated, nInput);
   for (size_t i = 0; i < nInterpolated; ++i) {
      for (size_t j = 0; j < nInput; ++j) {
         size_t index = (i + 2 * (nGhostBoundary - 1)) * nInput + j;
         pi_gsl::matrix_set(interpolationMatrix, i, j, (pi_gsl::gsl_real) table[index]);
      }
   }

// free table array
   free(table);
} // CONSTRUCTOR : PrecomputedInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::~PrecomputedInterpolationContext
// Description :  Destructor of PrecomputedInterpolationContext, frees memory
//
//-------------------------------------------------------------------------------------------------------
PrecomputedInterpolationContext::~PrecomputedInterpolationContext()
{
// free interpolation matrix
   pi_gsl::matrix_free(interpolationMatrix);
} // DESTRUCTOR : PrecomputedInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::GetWorkspaceSize
// Description :  Return NULL since PrecomputedInterpolationContext does not require a workspace
//
//-------------------------------------------------------------------------------------------------------
size_t PrecomputedInterpolationContext::GetWorkspaceSize() const {
   return 0;
} // FUNCTION : GetWorkspaceSize

//-------------------------------------------------------------------------------------------------------
// Function    :  PrecomputedInterpolationContext::InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
// Parameter   :  input          : Real input  array of size nInput
//                output         : Real output array of size 2 * (nInput - nGhostBoundary)
//                workspace      : Useless for PrecomputedInterpolationContext
//
//-------------------------------------------------------------------------------------------------------
void PrecomputedInterpolationContext::InterpolateReal(const real *input, real *output, char* workspace) const
{
   pi_gsl::vector_const_view in = pi_gsl::vector_const_view_array( input, nInput       );
   pi_gsl::vector_view out      = pi_gsl::vector_view_array      (output, nInterpolated);

   pi_gsl::blas_sgemv(CblasNoTrans, 1.0, interpolationMatrix, &in.vector, 0.0, &out.vector);
} // FUNCTION : InterpolateReal


//fixed interpolation constants
const real QuarticInterpolationContext::QuarticR[5] = { +35.0/2048.0, -252.0/2048.0, +1890.0/2048.0, +420.0/2048.0, -45.0/2048.0 };
const real QuarticInterpolationContext::QuarticL[5] = { -45.0/2048.0, +420.0/2048.0, +1890.0/2048.0, -252.0/2048.0, +35.0/2048.0 };

//-------------------------------------------------------------------------------------------------------
// Function    :  QuarticInterpolationContext::QuarticInterpolationContext
// Description :  Constructor of QuarticInterpolationContext
//
//-------------------------------------------------------------------------------------------------------
QuarticInterpolationContext::QuarticInterpolationContext(size_t nInput, size_t nGhostBoundary)
   :  InterpolationContext(nInput, nGhostBoundary)
{
#  ifdef GAMER_DEBUG
   if (nGhostBoundary != 2)
   {
      Aux_Error(ERROR_INFO, "QuarticInterpolationContext requires nGhostBoundary = 2!!\n");
   }
#  endif
}

//-------------------------------------------------------------------------------------------------------
// Function    :  QuarticInterpolationContext::GetWorkspaceSize
// Description :  Return 0 since QuarticInterpolationContext does not require a workspace
//
//-------------------------------------------------------------------------------------------------------
size_t QuarticInterpolationContext::GetWorkspaceSize() const
{
   return 0;
} // FUNCTION : GetWorkspaceSize

//-------------------------------------------------------------------------------------------------------
// Function    :  QuadraticInterpolationContext::InterpolateReal
// Description :  Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
// Parameter   :  input          : Real input  array of size nInput
//                output         : Real output array of size 2 * (nInput - nGhostBoundary)
//                workspace      : Useless for QuadraticInterpolationContext
//
//-------------------------------------------------------------------------------------------------------
void QuarticInterpolationContext::InterpolateReal(const real *input, real *output, char* workspace) const
{
   for (size_t i = nGhostBoundary; i < nInput - nGhostBoundary; ++i) {
      real r = 0, l = 0;
      for ( size_t j = 0; j < 5; ++j ) {
         l += QuarticInterpolationContext::QuarticL[j] * input[i - 2 + j];
         r += QuarticInterpolationContext::QuarticR[j] * input[i - 2 + j];
      }
      output[(i - nGhostBoundary) * 2    ] = l;
      output[(i - nGhostBoundary) * 2 + 1] = r;
   }
} // FUNCTION : InterpolateReal


const size_t NFast = 24;
const size_t fastNExtended[NFast] = {16, 18, 24, 32, 36, 40, 48, 54, 60, 64, 72, 80, 90, 96, 100, 108, 112, 120, 128, 135, 256, 288, 512, 540};

void InterpolationHandler::AddInterpolationContext(size_t nInput, size_t nGhostBoundary)
{

   if (contexts.find(nInput) == contexts.end()) {
//    ensure thread safety when adding new interpolation contexts
//    only one thread in OMP enviroment may add a new interpolation context
      #pragma omp critical
      {
//       for small N <= 32 pick precomputed interpolation with cost of N^2
         if ( nInput <= 32 ) {
               contexts.emplace(nInput, new QuarticInterpolationContext(nInput, nGhostBoundary));
               //contexts.emplace(nInput, new PrecomputedInterpolationContext(nInput, nGhostBoundary));
//       for large N >  32 use Gram-Fourier extension scheme with cost of N log(N)
         } else {
               size_t nExtension = 32, nDelta = 14;

               const size_t minimumExtensionSize = 24;
               const size_t maximumExtensionSize = 36;
               for (size_t i = 0; i < NFast; ++i) {
                  if (nInput + minimumExtensionSize < fastNExtended[i] && nInput + maximumExtensionSize >= fastNExtended[i]) {
                     nExtension = fastNExtended[i] - nInput;
                     break;
                  }
               }
               contexts.emplace(nInput, new GramFEInterpolationContext(nInput, nGhostBoundary, nExtension, nDelta));

         }
      }
   }
} // FUNCTION : AddContext

size_t InterpolationHandler::GetWorkspaceSize(size_t nInput, size_t nGhostBoundary) const
{
   return contexts.at(nInput)->GetWorkspaceSize();
} // FUNCTION : GetWorkspaceSize

void InterpolationHandler::InterpolateReal(real* input, real* output, size_t nInput, size_t nGhostBoundary, char* workspace, const bool UnwrapPhase, const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder) const
{
   auto context = contexts.at(nInput);

// unwrap phase
   context->Preprocess(input, UnwrapPhase);
// interpolate
   context->InterpolateReal(input, output, workspace);
// ensure monotonicity
   context->Postprocess(input, output, Monotonic, MonoCoeff, OppSign0thOrder);
} // FUNCTION : InterpolateReal


#endif // #endif SUPPORT_SPECTRAL_INTERPOLATION