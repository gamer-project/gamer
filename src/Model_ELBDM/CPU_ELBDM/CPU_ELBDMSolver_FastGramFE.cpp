#include "CUFLU.h"
#include "GAMER.h"

#if ( ( !defined(__CUDACC__) && defined(SUPPORT_FFTW) ) || ( defined(__CUDACC__) && defined(GRAMFE_ENABLE_GPU) ) )

#if ( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE && defined(ENABLE_FAST_GRAMFE) )

#include "GramFE_ExtensionTables.h"
#include "GSL.h"

#define GRAMFE_MM

#ifdef GRAMFE_MM
namespace gramfe_mm_gsl = gsl_double_precision;
#define gramfe_mm_real double
#else // #ifdef GRAMFE_MM
namespace gramfe_mm_gsl = gsl_single_precision;
#define gramfe_mm_real float
#endif // #ifdef GRAMFE_MM ... # else


#ifdef GRAMFE_EV
namespace gramfe_ev_gsl = gsl_double_precision;
#define gramfe_ev_real double
#else // #ifdef GRAMFE_FLOAT8
namespace gramfe_ev_gsl = gsl_single_precision;
#define gramfe_ev_real float
#endif // #ifdef GRAMFE_FLOAT8 ... # else



#ifdef __CUDACC__
#include "cuda_complex.h"
using gramfe_mm_complex_type = complex<gramfe_mm_real>;
using gramfe_ev_complex_type = complex<gramfe_ev_real>;
#else
#include <complex.h>
using gramfe_mm_complex_type = std::complex<gramfe_mm_real>;
using gramfe_ev_complex_type = std::complex<gramfe_ev_real>;
#endif

#ifndef __CUDACC__
//#define GRAMFE_EV

//-------------------------------------------------------------------------------------------------------
// Function    :  Factorial
// Description :  Compute the factorial of n.
// Parameter   :  n  : integer
// Return      :  n!
//-------------------------------------------------------------------------------------------------------
int Factorial(int n) {
     return (n==0) || (n==1) ? 1 : n* Factorial(n-1);
} // FUNCTION : Factorial

//-------------------------------------------------------------------------------------------------------
// Function    :  CosineTaylorExpansion
// Description :  Compute the taylor expansion of the cosine function at the point x with Nterms terms
// Parameter   :  x       : Point at which to evaluate Taylor expansion
//                NTerms  : Number of terms to retain in expansion
// Return      :  Value of expansion
//-------------------------------------------------------------------------------------------------------
gramfe_mm_real CosineTaylorExpansion(gramfe_mm_real x, int Nterms) {
   gramfe_mm_real result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * (1 / ((gramfe_mm_real) Factorial(2 * i))) * pow(x, 2 * i   );
   }

   return result;
} // FUNCTION : CosineTaylorExpansion

//-------------------------------------------------------------------------------------------------------
// Function    :  SineTaylorExpansion
// Description :  Compute the taylor expansion of the sine function at the point x with Nterms terms
// Parameter   :  x       : Point at which to evaluate Taylor expansion
//                NTerms  : Number of terms to retain in expansion
// Return      :  Value of expansion
//-------------------------------------------------------------------------------------------------------
gramfe_mm_real SineTaylorExpansion(gramfe_mm_real x, int Nterms) {
   gramfe_mm_real result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * ( 1 / ((gramfe_mm_real) Factorial(2 * i + 1)) ) * pow(x, 2 * i + 1);
   }

   return result;
} // FUNCTION : SineTaylorExpansion

void GramFE_SetupTimeEvolutionMatrix(real (*output)[2 * FLU_NXT], real dt, real dh, real Eta) {

// set up time evolution operator and filter
   gramfe_mm_real K, Filter, Coeff;
   gramfe_mm_complex_type ExpCoeff;
   gramfe_ev_complex_type (* out)    [FLU_NXT]        = (gramfe_ev_complex_type (*)[       FLU_NXT]) output;
   gramfe_mm_complex_type (*ifft2)   [GRAMFE_FLU_NXT] = (gramfe_mm_complex_type (*)[GRAMFE_FLU_NXT]) GramFE_IFFT;
   gramfe_mm_complex_type (* fft)    [GRAMFE_FLU_NXT] = (gramfe_mm_complex_type (*)[GRAMFE_FLU_NXT]) GramFE_FFT;
   gramfe_mm_complex_type GramFE_DFFT     [GRAMFE_FLU_NXT][GRAMFE_FLU_NXT];
   gramfe_mm_complex_type GramFE_IFFTDFFT [           PS2][GRAMFE_FLU_NXT];  // exp(- 1j * dt/(2*ELBDM_ETA) * k^2)
   gramfe_mm_complex_type GramFE_Evolution[           PS2][       FLU_NXT];
   const gramfe_mm_real filterDecay  = (gramfe_mm_real) 32.0 * (gramfe_mm_real) 2.302585092994046; // decay of k-space filter ( 32 * log(10) )
   const gramfe_mm_real filterDegree = (gramfe_mm_real) 100;                                     // degree of k-space filter
   const gramfe_mm_real kmax         = (gramfe_mm_real) M_PI / dh;                              // maximum value of k
   const gramfe_mm_real dk           = (gramfe_mm_real) + 2.0 * kmax / GRAMFE_FLU_NXT;           // k steps in k-space
   const gramfe_mm_real dT           = (gramfe_mm_real) - 0.5 * dt / Eta;                        // coefficient in time evolution operator

// naively "exp(1j * Coeff)" should give the exact time evolution of the free Schrödinger equation
// however, the time-evolution operator depends on higher-order derivatives
// since these are unavailable for a finite ghost boundary size, the series needs to be truncated
// the ideal Taylor expansion order using all available derivatives is FLU_GHOST_SIZE - 1
// for FLU_GHOST_SIZE == 8, 4 terms in the cosine series and 3 terms in the sine series are retained
   const int cosineNTerms = FLU_GHOST_SIZE / 2;
   const int sineNTerms   = cosineNTerms - (int) ((FLU_GHOST_SIZE % 2) == 0);

// set up momentum, filter and time evolution array

   for (int i=0; i<GRAMFE_FLU_NXT; i++)
   {
      K           = ( i <= GRAMFE_FLU_NXT/2 ) ? dk*i : dk*(i-GRAMFE_FLU_NXT);
      Filter      = exp(-filterDecay * pow(fabs(K/kmax), 2*filterDegree));
      Coeff       = SQR(K)*dT;
      ExpCoeff    = gramfe_mm_complex_type(CosineTaylorExpansion(Coeff, cosineNTerms), SineTaylorExpansion(Coeff, sineNTerms)) * Filter;
      for (int k=0; k<GRAMFE_FLU_NXT; k++)
      {
         GramFE_DFFT[i][k] = ExpCoeff * fft[i][k];
      }
   }
   gramfe_mm_gsl::matrix_complex_const_view extend     = gramfe_mm_gsl::matrix_complex_const_view_array(GramFE_Extend                     , GRAMFE_FLU_NXT, FLU_NXT);
   gramfe_mm_gsl::matrix_complex_const_view ifft       = gramfe_mm_gsl::matrix_complex_const_view_array(GramFE_IFFT                       , PS2           , GRAMFE_FLU_NXT);
   gramfe_mm_gsl::matrix_complex_const_view dfft       = gramfe_mm_gsl::matrix_complex_const_view_array((gramfe_mm_real*) GramFE_DFFT     , GRAMFE_FLU_NXT, GRAMFE_FLU_NXT);
   gramfe_mm_gsl::matrix_complex_view       ifftdfft   = gramfe_mm_gsl::matrix_complex_view_array      ((gramfe_mm_real*) GramFE_IFFTDFFT , PS2,            GRAMFE_FLU_NXT);
   gramfe_mm_gsl::matrix_complex_view       evolution  = gramfe_mm_gsl::matrix_complex_view_array      ((gramfe_mm_real*) GramFE_Evolution, PS2           , FLU_NXT);

   gramfe_mm_gsl::blas_cgemm(CblasNoTrans, CblasNoTrans, {1.0,  0.0}, &ifft.matrix,     &dfft.matrix,   {0.0, 0.0}, &ifftdfft.matrix);
   gramfe_mm_gsl::blas_cgemm(CblasNoTrans, CblasNoTrans, {1.0,  0.0}, &ifftdfft.matrix, &extend.matrix, {0.0, 0.0}, &evolution.matrix);

   for (size_t i = 0; i < PS2; ++i) {
      for (size_t j = 0; j < FLU_NXT; ++j) {
         out[i][j].real((real) GramFE_Evolution[i][j].real());
         out[i][j].imag((real) GramFE_Evolution[i][j].imag());
      }
   }
}
#endif // #ifndef __CUDACC__


// useful macros

// convert to 1D index with ghost boundary
# define to1D1(z,y,x) (  (z)                 * FLU_NXT * FLU_NXT +  (y)                 * FLU_NXT +  (x)                  )
// convert to 1D index without ghost boundary
# define to1D2(z,y,x) ( ((z)-FLU_GHOST_SIZE) * PS2     * PS2     + ((y)-FLU_GHOST_SIZE) * PS2     + ((x)-FLU_GHOST_SIZE)  )


#ifdef __CUDACC__
# define CGPU_FLU_BLOCK_SIZE_X FLU_BLOCK_SIZE_X
# define CGPU_FLU_BLOCK_SIZE_Y FLU_BLOCK_SIZE_Y
#else
# define CGPU_FLU_BLOCK_SIZE_X 1
# define CGPU_FLU_BLOCK_SIZE_Y 1
#endif



GPU_DEVICE
static uint get1D1(uint k, uint j, uint i, int XYZ) {
   switch ( XYZ )
   {
      case 0:  return to1D1( k, j, i );
      case 3:  return to1D1( k, i, j );
      case 6:  return to1D1( i, k, j );
   }
   return 0;
}

GPU_DEVICE
static uint get1D2(uint k, uint j, uint i, int XYZ) {
   switch ( XYZ )
   {
      case 0:  return to1D2( k, j, i );
      case 3:  return to1D2( k, i, j );
      case 6:  return to1D2( i, k, j );
   }
   return 0;
}

// multithreaded CPU/GPU loop over array respecting left and right ghost zones
# define CELL_LOOP( NCell, leftGhost, rightGhost )    for ( (Idx   = tid, \
                                                            (NStep = (NCell) - (leftGhost) - (rightGhost), \
                                                            (si    = Idx % NStep + (leftGhost), \
                                                             sj    = Idx / NStep))); \
                                                             Idx   < NColumnOnce * NStep; \
                                                            (Idx  += NThread, (si = Idx % NStep + (leftGhost), sj = Idx / NStep)) )


GPU_DEVICE
static void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                           real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                           int NPatchGroup,
                           const uint j_gap, const uint k_gap,
                           gramfe_ev_complex_type  s_In       [][FLU_NXT],
                           gramfe_ev_complex_type  s_Out      [][PS2],
                           gramfe_ev_complex_type* s_Evolve,
                           const bool FinalOut, const int XYZ, const real MinDens
                           );


//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver_GramFE
// Description :  CPU and GPU ELBDM kinematic solver based on computing Gram (FE) extension and evolving wave function using pseudo-spectral method on extended domain
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input pamameter "XYZ" is actually meaningless (if CONSERVE_MASS is off)
//                   --> Nevertheless, the symmetry in different directions will be broken if CONSERVE_MASS is on
//                2. The implementation is very similar to the function "CPU_FluidSolver_RTVD"
//
// Parameter   :  Flu_Array_In   : Array storing the input variables (only REAL/IMAG)
//                Flu_Array_Out  : Array to store the output variables (DENS/REAL/IMAG)
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta            : Particle mass / Planck constant
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                                 --> Meaningless since the operators along different directions
//                                     commute
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_ELBDMSolver_GramFE(   real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                 real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                 real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                 real g_Evolve   [][FLU_NXT * 2],
                                 const real dt, const real dh, const real Eta, const bool StoreFlux,
                                 const bool XYZ, const real MinDens )
#else
void CPU_ELBDMSolver_GramFE(     real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                 real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                 real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                 real g_Evolve   [][FLU_NXT * 2],
                                 const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                                 const bool XYZ, const real MinDens )
#endif
{

#  ifdef __CUDACC__
// create memories for columns of input field
   __shared__ gramfe_ev_complex_type s_In    [FLU_BLOCK_SIZE_Y][FLU_NXT];
   __shared__ gramfe_ev_complex_type s_Out   [FLU_BLOCK_SIZE_Y][PS2];

   const int NPatchGroup                           = NULL_INT;

#  else // #  ifdef __CUDACC__
// allocate memory on stack within loop for CPU run
   gramfe_ev_complex_type (*s_In)        [FLU_NXT] = NULL;
   gramfe_ev_complex_type (*s_Out)       [PS2]     = NULL;
#  endif // #  ifdef __CUDACC__ ... else

// time evolution matrix
#  ifdef __CUDACC__
   __shared__ gramfe_ev_complex_type s_Evolve[PS2 * FLU_NXT];
#  else // #  ifdef __CUDACC__
              gramfe_ev_complex_type* s_Evolve = (gramfe_ev_complex_type *) g_Evolve;
#  endif // #  ifdef __CUDACC__ ... else

#  ifdef __CUDACC__
   gramfe_ev_complex_type* s_LinEvolve = (gramfe_ev_complex_type *) s_Evolve;
   gramfe_ev_complex_type* g_LinEvolve = (gramfe_ev_complex_type *) g_Evolve;


   const uint tx           = threadIdx.x;
   const uint ty           = threadIdx.y;
   const uint tid          = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
   const uint NThread      = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block
   uint row, col;

// tranpose input evolution matrix from PS2 x FLU_NXT to FLU_NXT x PS2 for it to be in row-major order
   for (uint i = tid; i < PS2 * FLU_NXT; i += NThread) {
      row = i / FLU_NXT;
      col = i % FLU_NXT;
      s_LinEvolve[col * PS2 + row] = g_LinEvolve[i];
   }
#  endif // #  ifdef __CUDACC__


   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                                  0,              0, s_In, s_Out, s_Evolve, false, 0, MinDens);
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                     FLU_GHOST_SIZE,              0, s_In, s_Out, s_Evolve, false, 3, MinDens);
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Out, s_Evolve, true,  6, MinDens);
   } else  {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                                  0,              0, s_In, s_Out, s_Evolve, false, 6, MinDens);
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                                  0, FLU_GHOST_SIZE, s_In, s_Out, s_Evolve, false, 3, MinDens);
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Out, s_Evolve, true,  0, MinDens);
   }

} // FUNCTION : CUFLU_ELBDMSolver_GramFE



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use CPU/GPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on Gram-Fourier extension with pseudo-spectral solver on extended domain
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input variables
//                g_Fluid_Out    : Global memory array to store the output variables
//                NPatchGroup    : Number of patch groups to advance in parallel
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Eta            : Particle mass / Planck constant
//                j_gap          : Number of useless grids on each side in the j direction (j may not be equal to y)
//                k_gap          : Number of useless grids on each side in the k direction (k mya not be equal to z)
//                s_In           : Shared memory array to store the input data
//                s_Ae           : Shared memory array to store the even Gram polynomial coefficients
//                s_Ao           : Shared memory array to store the odd Gram polynomial coefficients
//                ExpCoeff       : Array to store the values of the time evolution operator with filter
//                FinalOut       : true --> store the updated data to g_Fluid_Out
//                XYZ            : 0 : Update the solution in the x direction
//                                 3 : Update the solution in the y direction
//                                 6 : Update the solution in the z direction
//                                 --> This parameter is also used to determine the place to store the output fluxes
//                MinDens        : Minimum allowed density
//                Workspace      : Workspace for forward GPU FFT (useless in CPU mode)
//                WorkspaceInv   : Workspace for inverse GPU FFT (useless in CPU mode)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CUFLU_Advance(  real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                     real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                     int NPatchGroup,
                     const uint j_gap, const uint k_gap,
                     gramfe_ev_complex_type  s_In       [][FLU_NXT],
                     gramfe_ev_complex_type  s_Out      [][PS2],
                     gramfe_ev_complex_type* s_Evolve,
                     const bool FinalOut,
                     const int XYZ, const real MinDens
                  )
{
   const uint size_j       = FLU_NXT -  j_gap * 2; // number of y-columns to be updated
   const uint size_k       = FLU_NXT -  k_gap * 2; // number of z-columns to be updated
   const uint NColumnTotal = size_j * size_k;      // total number of data columns to be updated

// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
#     ifdef __CUDACC__
      const int bx = blockIdx.x;
#     else

//    create arrays for columns of various intermediate fields on the stack
      gramfe_ev_complex_type* s_In_1PG  = (gramfe_ev_complex_type* ) malloc( FLU_NXT * sizeof(gramfe_ev_complex_type) ); // allocate memory for fourier transform
      gramfe_ev_complex_type* s_Out_1PG = (gramfe_ev_complex_type* ) malloc(     PS2 * sizeof(gramfe_ev_complex_type) ); // allocate memory for fourier transform

      gramfe_ev_gsl::vector_complex_const_view Input_view  = gramfe_ev_gsl::vector_complex_const_view_array ((gramfe_ev_gsl::gsl_real*) s_In_1PG ,      FLU_NXT);
      gramfe_ev_gsl::vector_complex_view       Output_view = gramfe_ev_gsl::vector_complex_view_array       ((gramfe_ev_gsl::gsl_real*) s_Out_1PG, PS2         );
      gramfe_ev_gsl::matrix_complex_const_view Evo_view    = gramfe_ev_gsl::matrix_complex_const_view_array ((gramfe_ev_gsl::gsl_real*) s_Evolve , PS2, FLU_NXT);

//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_Out )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif
      {

#        ifdef __CUDACC__
//       use two-dimensional thread blocks in GPU mode
         const uint tx            = threadIdx.x;
         const uint ty            = threadIdx.y;
#        else  // # ifdef __CUDACC__
//       every block just has a single thread with temporary memory on the stack in CPU mode
         const uint tx            = 0;
         const uint ty            = 0;


         s_In     = (gramfe_ev_complex_type (*)[FLU_NXT]) s_In_1PG;
         s_Out    = (gramfe_ev_complex_type (*)[PS2])    s_Out_1PG;

#        endif // # ifdef __CUDACC__ ... # else

         const uint tid          = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
         const uint NThread      = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block

         uint j, k;             // (i,j,k): array indices used in g_Fluid_In

         uint Column0 = 0;      // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;  // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out

         uint si, sj;           // array indices used in the shared memory array
         uint NStep;            // number of iterations for updating each column

         uint NColumnOnce = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration

         gramfe_ev_complex_type Psi_New;             // buffer for left and right Gram coefficients
         gramfe_ev_real   Amp_New, Re_New, Im_New;  // store density, real and imaginary part to apply minimum density check

//       loop over all data columns
         while ( Column0 < NColumnTotal )
         {

//          1. load data into shared memory
            CELL_LOOP(FLU_NXT, 0, 0) {
               j = j_gap + ( sj + Column0 ) % size_j;
               k = k_gap + ( sj + Column0 ) / size_j;

//             1.1 determine the array indices for loading global memory data along different directions
               Idx1 = get1D1( k, j, si, XYZ );

//             1.2 load the interior data into shared memory
               s_In[sj][si].real(g_Fluid_In[bx][0][Idx1]);
               s_In[sj][si].imag(g_Fluid_In[bx][1][Idx1]);
            } // CELL_LOOP(FLU_NXT, 0, 0)


//          1.4 sync data read into S_In
#           ifdef __CUDACC__
            __syncthreads();

//          2. evolve wave function
            CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            {
               Psi_New = {0, 0};

               for (int t=0; t < FLU_NXT; t++) {
                  Psi_New += s_Evolve[(si - FLU_GHOST_SIZE) + t * PS2] * s_In[sj][t];
               } // for t

               s_Out[sj][si - FLU_GHOST_SIZE] = Psi_New;
            }

            __syncthreads();
#           else
            gramfe_ev_gsl::blas_cgemv(CblasNoTrans, {1.0, 0.0}, &Evo_view.matrix, &Input_view.vector, {0.0, 0.0}, &Output_view.vector);
#           endif


            CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            {
               //printf("si %d %3.7e + i %3.7e\n", si, s_Out[sj][si - FLU_GHOST_SIZE].real(), s_Out[sj][si - FLU_GHOST_SIZE].imag());
            }

            if ( FinalOut )
            {
               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  Re_New = s_Out[sj][si - FLU_GHOST_SIZE].real();
                  Im_New = s_Out[sj][si - FLU_GHOST_SIZE].imag();

   //             4.1 write FFT array back to output array
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx2 = get1D2( k, j, si, XYZ );

                  Amp_New =  SQR(Re_New) + SQR(Im_New);

   //                apply the the minimum density check
                  if ( Amp_New < MinDens )
                  {
                     const real Rescale = SQRT( MinDens / (real)Amp_New );

                     Re_New *= Rescale;
                     Im_New *= Rescale;
                     Amp_New = MinDens;
                  }
                  g_Fluid_Out[bx][DENS][Idx2] = Amp_New;
                  g_Fluid_Out[bx][REAL][Idx2] = Re_New;
                  g_Fluid_Out[bx][IMAG][Idx2] = Im_New;
               }
            } else { // if ( FinalOut )

               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  Re_New = s_Out[sj][si - FLU_GHOST_SIZE].real();
                  Im_New = s_Out[sj][si - FLU_GHOST_SIZE].imag();
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx1 = get1D1( k, j,si, XYZ );
                  g_Fluid_In[bx][0][Idx1] = Re_New;
                  g_Fluid_In[bx][1][Idx1] = Im_New;
               } // CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            } // if ( FinalOut ) ... else

#           ifdef  __CUDACC__
            __syncthreads();
#           endif

//          4.2 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      } // # pragma  for (int bx=0; bx<NPatchGroup; bx++)
#     ifndef __CUDACC__
      free(s_In_1PG);
      free(s_Out_1PG);
#     endif
   } // # pragma omp parallel
} // FUNCTION : CUFLU_Advance


#endif // #if ( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE)
#endif // #if ( ( !defined(__CUDACC__) && defined(SUPPORT_FFTW) ) || ( defined(__CUDACC__) && defined(GRAMFE_ENABLE_GPU) ) )