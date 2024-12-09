#include "GAMER.h"
#include "CUFLU.h"

#if (  ( !defined(__CUDACC__) && defined(SUPPORT_GSL) )  ||  defined(__CUDACC__)  )

#if ( GRAMFE_SCHEME == GRAMFE_MATMUL )


#ifdef __CUDACC__

// implement a complex type with the required operations for matrix multiplication
template <typename T>
class complex {
public:
   T re;
   T im;

// unfortunately, we cannot defined custom constructors here
// otherwise, CUDA will complain about unsupported dynamic initialisation of shared arrays

   __device__ complex<T> operator+(const complex<T>& other) const {
      return complex<T>(re + other.re, im + other.im);
   }

   __device__ complex<T> operator*(const complex<T>& other) const {
      complex<T> out;
      out.real(re * other.re - im * other.im);
      out.imag(re * other.im + im * other.re);
      return out;
   }

   __device__ complex<T>& operator+=(const complex<T>& other) {
      re += other.re;
      im += other.im;
      return *this;
   }

   __device__ T real() const {return re;}
   __device__ T imag() const {return im;}
   __device__ void real(T r) {re = r;}
   __device__ void imag(T i) {im = i;}
};

using gramfe_matmul_complex_type = complex<gramfe_matmul_float>;
// matrix-multiplication is carried out manually
// therefore the types of the input vector (gramfe_input_complex_type) and the matrix (gramfe_matmul_complex_type) can be different
using gramfe_input_complex_type  = complex<real>;

#else // #ifdef __CUDACC__

#include <complex.h>
#include "GSL.h"


// for matrix-multiplication via GSL using BLAS the types of the input vector (gramfe_input_complex_type)
// and the matrix (gramfe_matmul_complex_type) must be the same
using gramfe_matmul_complex_type = std::complex<gramfe_matmul_float>;
using gramfe_input_complex_type  = std::complex<gramfe_matmul_float>;

// precision of matrix multiplication
#ifdef GRAMFE_MATMUL_FLOAT8
namespace gramfe_matmul_gsl = gsl_double_precision;
#else
namespace gramfe_matmul_gsl = gsl_single_precision;
#endif

#endif // #ifdef __CUDACC__ ... else ...

#ifdef __CUDACC__
# define CGPU_FLU_BLOCK_SIZE_X FLU_BLOCK_SIZE_X
# define CGPU_FLU_BLOCK_SIZE_Y FLU_BLOCK_SIZE_Y
#else
# define CGPU_FLU_BLOCK_SIZE_X 1
# define CGPU_FLU_BLOCK_SIZE_Y 1
#endif

// useful macros
// convert to 1D index with ghost boundary
# define to1D1(z,y,x) (   (z)                 * FLU_NXT * FLU_NXT +  (y)                 * FLU_NXT +  (x)                  )
// convert to 1D index without ghost boundary
# define to1D2(z,y,x) (  ((z)-FLU_GHOST_SIZE) * PS2     * PS2     + ((y)-FLU_GHOST_SIZE) * PS2     + ((x)-FLU_GHOST_SIZE)  )


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
#define CELL_LOOP( NCell, leftGhost, rightGhost )  for ( (Idx   = tid, \
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
                           gramfe_input_complex_type s_In       [][FLU_NXT],
                           gramfe_input_complex_type s_Out      [][PS2],
                           gramfe_matmul_complex_type* s_TimeEvo,
                           const bool FinalOut, const int XYZ, const real MinDens );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver_GramFE_MAMTUL
// Description :  CPU and GPU ELBDM kinematic solver based on computing Gram (FE) extension and evolving wave function
//                using pseudo-spectral method on extended domain using matrix multiplication
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//
// Parameter   :  g_Fluid_In  : Array storing the input variables (only REAL/IMAG)
//                g_Fluid_Out : Array to store the output variables (DENS/REAL/IMAG)
//                g_Flux      : Useless
//                g_TimeEvo   : Array to store complex PS2 * FLU_NXT matrix that evolves wave function in time
//                NPatchGroup : Number of patch groups to be evaluated
//                dt          : Time interval to advance solution
//                dh          : Grid size
//                Eta         : Particle mass / Planck constant
//                XYZ         : true  : x->y->z ( forward sweep)
//                              false : z->y->x (backward sweep)
//                              --> Meaningless since the operators along different directions
//                                  commute
//                MinDens     : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_ELBDMSolver_GramFE_MATMUL( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      gramfe_matmul_float g_TimeEvo[][ FLU_NXT*2 ],
                                      const real dt, const real dh, const real Eta, const bool StoreFlux,
                                      const bool XYZ, const real MinDens )
#else
void CPU_ELBDMSolver_GramFE_MATMUL(   real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      gramfe_matmul_float g_TimeEvo[][ FLU_NXT*2 ],
                                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                                      const bool XYZ, const real MinDens )
#endif
{

#  ifdef __CUDACC__
// create memories for columns of input field
   __shared__ gramfe_input_complex_type s_In [FLU_BLOCK_SIZE_Y][FLU_NXT];
   __shared__ gramfe_input_complex_type s_Out[FLU_BLOCK_SIZE_Y][PS2];

   const int NPatchGroup = NULL_INT;

#  else
// allocate memory on heap within loop for CPU run
//###OPTIMIZATION: allocate memory only once right here
   gramfe_input_complex_type (*s_In) [FLU_NXT] = NULL;
   gramfe_input_complex_type (*s_Out)[PS2]     = NULL;
#  endif

// time evolution matrix
// GPU: transpose input evolution matrix from PS2 x FLU_NXT to FLU_NXT x PS2 for it to be in row-major order
// CPU: no transposition since matrix already is in column-major order
#  ifdef __CUDACC__
   __shared__ gramfe_matmul_complex_type s_TimeEvo[PS2 * FLU_NXT];

   gramfe_matmul_complex_type* s_LinEvolve = (gramfe_matmul_complex_type *) s_TimeEvo;
   gramfe_matmul_complex_type* g_LinEvolve = (gramfe_matmul_complex_type *) g_TimeEvo;
   uint row, col;

   const uint tx      = threadIdx.x;
   const uint ty      = threadIdx.y;
   const uint tid     = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
   const uint NThread = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block

   for (uint i=tid; i<PS2*FLU_NXT; i+=NThread) {
      row = i / FLU_NXT;
      col = i % FLU_NXT;
      s_LinEvolve[ col*PS2 + row ] = g_LinEvolve[i];
   }

   __syncthreads();

#  else // #ifdef __CUDACC__
   gramfe_matmul_complex_type* s_TimeEvo = (gramfe_matmul_complex_type *) g_TimeEvo;
#  endif // #ifdef __CUDACC__ ... else ...


   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                                  0,              0, s_In, s_Out, s_TimeEvo, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                     FLU_GHOST_SIZE,              0, s_In, s_Out, s_TimeEvo, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Out, s_TimeEvo, true,  6, MinDens );
   } else {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                                  0,              0, s_In, s_Out, s_TimeEvo, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                                  0, FLU_GHOST_SIZE, s_In, s_Out, s_TimeEvo, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Out, s_TimeEvo, true,  0, MinDens );
   }

} // FUNCTION : CUFLU_ELBDMSolver_GramFE_MATMUL



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use CPU/GPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on Gram-Fourier extension with pseudo-spectral solver on extended domain
//
// Parameter   :  g_Fluid_In  : Global memory array storing the input variables
//                g_Fluid_Out : Global memory array to store the output variables
//                NPatchGroup : Number of patch groups to advance in parallel
//                dt          : Time interval to advance solution
//                _dh         : 1 / grid size
//                Eta         : Particle mass / Planck constant
//                j_gap       : Number of useless grids on each side in the j direction (j may not be equal to y)
//                k_gap       : Number of useless grids on each side in the k direction (k mya not be equal to z)
//                s_In        : Shared memory array to store the input data
//                s_Out       : Shared memory array to store the output data
//                s_TimeEvo   : Shared memory array that stores time evolution matrix
//                ExpCoeff    : Array to store the values of the time evolution operator with filter
//                FinalOut    : true --> store the updated data to g_Fluid_Out
//                XYZ         : 0 : Update the solution in the x direction
//                              3 : Update the solution in the y direction
//                              6 : Update the solution in the z direction
//                              --> This parameter is also used to determine the place to store the output fluxes
//                MinDens     : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                    real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                    int NPatchGroup,
                    const uint j_gap, const uint k_gap,
                    gramfe_input_complex_type s_In [][FLU_NXT],
                    gramfe_input_complex_type s_Out[][PS2],
                    gramfe_matmul_complex_type* s_TimeEvo,
                    const bool FinalOut,
                    const int XYZ, const real MinDens )
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
//    create arrays for columns of various intermediate fields on the heap
      gramfe_matmul_complex_type* s_In_1PG  = (gramfe_matmul_complex_type*) malloc( FLU_NXT * sizeof(gramfe_matmul_complex_type) );
      gramfe_matmul_complex_type* s_Out_1PG = (gramfe_matmul_complex_type*) malloc(     PS2 * sizeof(gramfe_matmul_complex_type) );

      gramfe_matmul_gsl::vector_complex_const_view Input_view  = gramfe_matmul_gsl::vector_complex_const_view_array ( (gramfe_matmul_gsl::gsl_real*) s_In_1PG ,      FLU_NXT );
      gramfe_matmul_gsl::vector_complex_view       Output_view = gramfe_matmul_gsl::vector_complex_view_array       ( (gramfe_matmul_gsl::gsl_real*) s_Out_1PG, PS2          );
      gramfe_matmul_gsl::matrix_complex_const_view Evo_view    = gramfe_matmul_gsl::matrix_complex_const_view_array ( (gramfe_matmul_gsl::gsl_real*) s_TimeEvo, PS2, FLU_NXT );

//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_Out )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif
      {

#        ifdef __CUDACC__
//       use two-dimensional thread blocks in GPU mode
         const uint tx = threadIdx.x;
         const uint ty = threadIdx.y;

         // define register variables for conversion of complex types
         gramfe_matmul_complex_type Psi_In, Psi_New;

#        else  // # ifdef __CUDACC__
//       every block just has a single thread with temporary memory on the stack in CPU mode
         const uint tx = 0;
         const uint ty = 0;


         s_In  = (gramfe_input_complex_type (*)[FLU_NXT]) s_In_1PG;
         s_Out = (gramfe_input_complex_type (*)[PS2])    s_Out_1PG;

#        endif // # ifdef __CUDACC__ ... else ...

         const uint tid     = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
         const uint NThread = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block

         uint j, k;             // (i,j,k): array indices used in g_Fluid_In

         uint Column0 = 0;      // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;  // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out

         uint si, sj;           // array indices used in the shared memory array
         uint NStep;            // number of iterations for updating each column

         uint NColumnOnce = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration

         real Amp_New, Re_New, Im_New;  // store density, real and imaginary part to apply minimum density check

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
            }


//          2. evolve wave function via matrix multiplication
#           ifdef __CUDACC__
            __syncthreads();

            CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            {
               Psi_New = {0, 0};

               for (int t=0; t < FLU_NXT; t++) {
                  Psi_In  = {(gramfe_matmul_float) s_In[sj][t].real(), (gramfe_matmul_float) s_In[sj][t].imag()};
                  Psi_New += s_TimeEvo[(si - FLU_GHOST_SIZE) + t * PS2] * Psi_In;
               }

               s_Out[sj][si - FLU_GHOST_SIZE] = {(real) Psi_New.real(), (real) Psi_New.imag()};
            }

            __syncthreads();
#           else
            gramfe_matmul_gsl::blas_cgemv(CblasNoTrans, {1.0, 0.0}, &Evo_view.matrix, &Input_view.vector, {0.0, 0.0}, &Output_view.vector);
#           endif


//          3. store output back to global memory
            if ( FinalOut )
            {
               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  Re_New = s_Out[sj][si - FLU_GHOST_SIZE].real();
                  Im_New = s_Out[sj][si - FLU_GHOST_SIZE].imag();

//                3.1 write FFT array back to output array
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
               } // CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            } else { // if ( FinalOut )

               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  Re_New = s_Out[sj][si - FLU_GHOST_SIZE].real();
                  Im_New = s_Out[sj][si - FLU_GHOST_SIZE].imag();
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx1 = get1D1( k, j, si, XYZ );
                  g_Fluid_In[bx][0][Idx1] = Re_New;
                  g_Fluid_In[bx][1][Idx1] = Im_New;
               } // CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            } // if ( FinalOut ) ... else ...

#           ifdef  __CUDACC__
            __syncthreads();
#           endif

//          3.2 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      } // #pragma for (int bx=0; bx<NPatchGroup; bx++)
#     ifndef __CUDACC__
      free(s_In_1PG);
      free(s_Out_1PG);
#     endif
   } // #pragma omp parallel
} // FUNCTION : CUFLU_Advance



#endif // #if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
#endif // #if (  ( !defined(__CUDACC__) && defined(SUPPORT_GSL) )  ||  defined(__CUDACC__)  )
