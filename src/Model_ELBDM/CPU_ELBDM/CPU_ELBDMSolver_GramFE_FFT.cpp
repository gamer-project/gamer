#include "GAMER.h"
#include "CUFLU.h"

#if (  ( !defined(__CUDACC__) && defined(SUPPORT_FFTW) )  ||  defined(__CUDACC__)  )

#if ( GRAMFE_SCHEME == GRAMFE_FFT )


#include "GramFE_ExtensionTables.h"

// useful macros

// convert to 1D index with ghost boundary
# define to1D1(z,y,x) (   (z)                 * FLU_NXT * FLU_NXT +  (y)                 * FLU_NXT +  (x)                  )
// convert to 1D index without ghost boundary
# define to1D2(z,y,x) (  ((z)-FLU_GHOST_SIZE) * PS2     * PS2     + ((y)-FLU_GHOST_SIZE) * PS2     + ((x)-FLU_GHOST_SIZE)  )

// use cufftdx library for FFTs on GPU
#ifdef __CUDACC__

using forward_workspace_type = typename FFT::workspace_type;
using inverse_workspace_type = typename IFFT::workspace_type;

// overload cufftdx's complex_type to allow for multiplication, addition and subtraction of complex and real numbers

// multiplication of complex and real
template<class OtherType>
__device__ __forceinline__ complex_type operator*(const complex_type& a, const OtherType& other) {
      complex_type result(a);
      result *= other;
      return result;
}

// multiplication of real and complex
template<class OtherType>
__device__ __forceinline__ complex_type operator*(const OtherType& other, const complex_type& a) {
      return a * other;
}

// multiplication of complex and complex
__device__ __forceinline__ complex_type operator*(const complex_type& a, const complex_type& b) {
      complex_type result(a);
      result *= b;
      return result;
}

// addition of complex and complex
__device__ __forceinline__ complex_type operator+(const complex_type& a, const complex_type& b) {
      complex_type result(a);
      result += b;
      return result;
}

// subtraction of complex and complex
__device__ __forceinline__ complex_type operator-(const complex_type& a, const complex_type& b) {
      complex_type result(a);
      result -= b;
      return result;
}

#else   // #ifdef __CUDACC__

extern gramfe_fftw::complex_plan_1d FFTW_Plan_ExtPsi, FFTW_Plan_ExtPsi_Inv;

#if ( SUPPORT_FFTW == FFTW3 )

#include <complex.h>

using complex_type = std::complex<gramfe_fft_float>;

#else // #if ( SUPPORT_FFTW == FFTW3 )

// derive from gramfe_fftw::fft_complex which is an alias for FFTW2's complex_type in gramfe_fft_float precision
// to allow for complex arithmetic operations
struct complex_type : public gramfe_fftw::fft_complex {
   complex_type() {
   }

   complex_type(gramfe_fft_float re, gramfe_fft_float im) {
      real(re);
      imag(im);
   }

   complex_type& operator=(const complex_type& other) {
      real(other.real());
      imag(other.imag());
      return *this;
   }

   gramfe_fft_float real() const {
      return c_re(*this);
   }

   gramfe_fft_float imag() const {
      return c_im(*this);
   }

   void real(gramfe_fft_float re) {
      c_re(*this) = re;
   }

   void imag(gramfe_fft_float im) {
      c_im(*this) = im;
   }

   complex_type operator*(const complex_type& other) {
         return complex_type(this->real() * other.real() - this->imag() * other.imag(), this->real() * other.imag() + this->imag() * other.real());
   }

   complex_type operator+(const complex_type& other) {
         return complex_type(this->real() + other.real(), this->imag() + other.imag());
   }

   complex_type operator-(const complex_type& other) {
         return complex_type(this->real() - other.real(), this->imag() - other.imag());
   }

   complex_type& operator*=(const complex_type& other) {
      (*this) = (*this) * other;
      return *this;
   }

   complex_type& operator+=(const complex_type& other) {
      (*this) = (*this) + other;
      return *this;
   }
}; // struct complex_type : public gramfe_fftw::fft_complex


template<class OtherType>
complex_type operator*(const complex_type& a, const OtherType& other) {
      return complex_type (a.real() * other, a.imag() * other);
}

template<class OtherType>
complex_type operator*(const OtherType& other, const complex_type& a) {
      return a * other;
}
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... else ...


// no workspaces required in CPU solver
using forward_workspace_type = bool;
using inverse_workspace_type = bool;
#endif // #ifdef __CUDACC__



#ifdef __CUDACC__
# define CGPU_FLU_BLOCK_SIZE_X FFT::block_dim.x
# define CGPU_FLU_BLOCK_SIZE_Y FFT::block_dim.y
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
#define CELL_LOOP( NCell, leftGhost, rightGhost )  for ( (Idx   = tid, \
                                                         (NStep = (NCell) - (leftGhost) - (rightGhost), \
                                                         (si    = Idx % NStep + (leftGhost), \
                                                          sj    = Idx / NStep))); \
                                                          Idx   < NColumnOnce * NStep; \
                                                         (Idx  += NThread, (si = Idx % NStep + (leftGhost), sj = Idx / NStep)) )


GPU_DEVICE
static void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                           real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                           int NPatchGroup,
                           const gramfe_fft_float dt, const gramfe_fft_float _dh, const gramfe_fft_float Eta,
                           const uint j_gap, const uint k_gap,
                           complex_type s_In    [][GRAMFE_FLU_NXT],
                           complex_type s_Ae    [][GRAMFE_NDELTA],
                           complex_type s_Ao    [][GRAMFE_NDELTA],
                           complex_type ExpCoeff[],
                           const bool FinalOut, const int XYZ, const gramfe_fft_float MinDens,
                           forward_workspace_type Workspace,
                           inverse_workspace_type WorkspaceInv );




//-------------------------------------------------------------------------------------------------------
// Function    :  Factorial
// Description :  Compute the factorial of n.
// Parameter   :  n  : integer
// Return      :  n!
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
int Factorial(int n) {
     return (n==0) || (n==1) ? 1 : n* Factorial(n-1);
} // FUNCTION : Factorial



//-------------------------------------------------------------------------------------------------------
// Function    :  CosineTaylorExpansion
// Description :  Compute the Taylor expansion of the cosine function at the point x with Nterms terms
// Parameter   :  x       : Point at which to evaluate Taylor expansion
//                NTerms  : Number of terms to retain in expansion
// Return      :  Value of expansion
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
gramfe_fft_float CosineTaylorExpansion(gramfe_fft_float x, int Nterms) {
   gramfe_fft_float result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * (1 / ((gramfe_fft_float) Factorial(2 * i))) * pow(x, 2 * i   );
   }

   return result;
} // FUNCTION : CosineTaylorExpansion



//-------------------------------------------------------------------------------------------------------
// Function    :  SineTaylorExpansion
// Description :  Compute the Taylor expansion of the sine function at the point x with Nterms terms
// Parameter   :  x       : Point at which to evaluate Taylor expansion
//                NTerms  : Number of terms to retain in expansion
// Return      :  Value of expansion
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
gramfe_fft_float SineTaylorExpansion(gramfe_fft_float x, int Nterms) {
   gramfe_fft_float result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * ( 1 / ((gramfe_fft_float) Factorial(2 * i + 1)) ) * pow(x, 2 * i + 1);
   }

   return result;
} // FUNCTION : SineTaylorExpansion



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver_GramFE_FFT
// Description :  CPU and GPU ELBDM kinematic solver based on computing Gram (FE) extension and evolving wave function
//                using pseudo-spectral method on extended domain
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input pamameter "XYZ" is actually meaningless (if CONSERVE_MASS is off)
//                   --> Nevertheless, the symmetry in different directions will be broken if CONSERVE_MASS is on
//                2. The implementation is very similar to the function "CPU_FluidSolver_RTVD"
//
// Parameter   :  Flu_Array_In  : Array storing the input variables (only REAL/IMAG)
//                Flu_Array_Out : Array to store the output variables (DENS/REAL/IMAG)
//                NPatchGroup   : Number of patch groups to be evaluated
//                dt            : Time interval to advance solution
//                dh            : Grid size
//                Eta           : Particle mass / Planck constant
//                XYZ           : true  : x->y->z ( forward sweep)
//                                false : z->y->x (backward sweep)
//                                --> Meaningless since the operators along different directions
//                                    commute
//                MinDens       : Minimum allowed density
//                Only in GPU mode:
//                Workspace     : Workspace for forward GPU FFT
//                WorkspaceInv  : Workspace for inverse GPU FFT
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__launch_bounds__(FFT::max_threads_per_block)
__global__
void CUFLU_ELBDMSolver_GramFE_FFT( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                   real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                   const bool XYZ, const real MinDens,
                                   typename FFT::workspace_type  Workspace,
                                   typename IFFT::workspace_type WorkspaceInv )
#else
void CPU_ELBDMSolver_GramFE_FFT(   real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                   real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                   const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                                   const bool XYZ, const real MinDens )
#endif
{

#  ifdef __CUDACC__
// shared memory array for cufftdx
   extern __shared__ complex_type shared_mem[];

// create memories for columns of various intermediate fields in shared GPU memory
   complex_type (*s_In)[GRAMFE_FLU_NXT]    = (complex_type (*)[GRAMFE_FLU_NXT]) (shared_mem);
// even extension coefficients
   complex_type (*s_Ae)[GRAMFE_NDELTA]     = (complex_type (*)[GRAMFE_NDELTA])  (shared_mem + CGPU_FLU_BLOCK_SIZE_Y * (GRAMFE_FLU_NXT                ));
// odd extension coefficients
   complex_type (*s_Ao)[GRAMFE_NDELTA]     = (complex_type (*)[GRAMFE_NDELTA])  (shared_mem + CGPU_FLU_BLOCK_SIZE_Y * (GRAMFE_FLU_NXT + GRAMFE_NDELTA));
   const int NPatchGroup                   = NULL_INT;

#  else // #ifdef __CUDACC__
// allocate memory on heap within loop for CPU run
   complex_type (*s_In)   [GRAMFE_FLU_NXT] = NULL;
   complex_type (*s_Ae)   [GRAMFE_NDELTA]  = NULL;
   complex_type (*s_Ao)   [GRAMFE_NDELTA]  = NULL;
   const gramfe_fft_float _dh              = gramfe_fft_float(1.0)/dh;
   bool Workspace                          = NULL_BOOL;
   bool WorkspaceInv                       = NULL_BOOL;
#  endif // #ifdef __CUDACC__ ... else ...


// set up time evolution operator and filter
   gramfe_fft_float K;
   gramfe_fft_float Filter;               // exp(-filterDecay * (k/kMax)**(2*filterDegree))
   gramfe_fft_float Coeff;                // dT * k^2
   complex_type ExpCoeff[GRAMFE_FLU_NXT]; // exp(- 1j * dt/(2*ELBDM_ETA) * k^2)

   const gramfe_fft_float filterDecay  = (gramfe_fft_float) 32.0 * (gramfe_fft_float) 2.302585092994046; // decay of k-space filter (32*log(10))
   const gramfe_fft_float filterDegree = (gramfe_fft_float) 100;                                         // degree of k-space filter
   const gramfe_fft_float kmax         = (gramfe_fft_float) M_PI * _dh;                                  // maximum value of k
   const gramfe_fft_float dk           = (gramfe_fft_float) + 2.0 * kmax / GRAMFE_FLU_NXT;               // k steps in k-space
   const gramfe_fft_float dT           = (gramfe_fft_float) - 0.5 * dt / Eta;                            // coefficient in time evolution operator
   const gramfe_fft_float Norm         = (gramfe_fft_float) + 1.0 / GRAMFE_FLU_NXT;                      // norm for inverse Fourier transform

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
      ExpCoeff[i] = complex_type(CosineTaylorExpansion(Coeff, cosineNTerms), SineTaylorExpansion(Coeff, sineNTerms)) * Norm * Filter;
   }


   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup, dt, _dh, Eta,
                                  0,              0, s_In, s_Ae, s_Ao, ExpCoeff, false, 0, MinDens, Workspace, WorkspaceInv );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup, dt, _dh, Eta,
                     FLU_GHOST_SIZE,              0, s_In, s_Ae, s_Ao, ExpCoeff, false, 3, MinDens, Workspace, WorkspaceInv );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup, dt, _dh, Eta,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Ae, s_Ao, ExpCoeff, true,  6, MinDens, Workspace, WorkspaceInv );
   } else  {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup, dt, _dh, Eta,
                                  0,              0, s_In, s_Ae, s_Ao, ExpCoeff, false, 6, MinDens, Workspace, WorkspaceInv );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup, dt, _dh, Eta,
                                  0, FLU_GHOST_SIZE, s_In, s_Ae, s_Ao, ExpCoeff, false, 3, MinDens, Workspace, WorkspaceInv );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, NPatchGroup, dt, _dh, Eta,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Ae, s_Ao, ExpCoeff, true,  0, MinDens, Workspace, WorkspaceInv );
   }

} // FUNCTION : CUFLU_ELBDMSolver_GramFE_FFT



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use CPU/GPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on Gram-Fourier extension with pseudo-spectral solver on extended domain
//
// Parameter   :  g_Fluid_In   : Global memory array storing the input variables
//                g_Fluid_Out  : Global memory array to store the output variables
//                NPatchGroup  : Number of patch groups to advance in parallel
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Eta          : Particle mass / Planck constant
//                j_gap        : Number of useless grids on each side in the j direction (j may not be equal to y)
//                k_gap        : Number of useless grids on each side in the k direction (k mya not be equal to z)
//                s_In         : Shared memory array to store the input data
//                s_Ae         : Shared memory array to store the even Gram polynomial coefficients
//                s_Ao         : Shared memory array to store the odd Gram polynomial coefficients
//                ExpCoeff     : Array to store the values of the time evolution operator with filter
//                FinalOut     : true --> store the updated data to g_Fluid_Out
//                XYZ          : 0 : Update the solution in the x direction
//                               3 : Update the solution in the y direction
//                               6 : Update the solution in the z direction
//                               --> This parameter is also used to determine the place to store the output fluxes
//                MinDens      : Minimum allowed density
//                Workspace    : Workspace for forward GPU FFT (useless in CPU mode)
//                WorkspaceInv : Workspace for inverse GPU FFT (useless in CPU mode)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                    real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                    int NPatchGroup,
                    const gramfe_fft_float dt, const gramfe_fft_float _dh, const gramfe_fft_float Eta,
                    const uint j_gap, const uint k_gap,
                    complex_type s_In    [][GRAMFE_FLU_NXT],
                    complex_type s_Ae    [][GRAMFE_NDELTA],
                    complex_type s_Ao    [][GRAMFE_NDELTA],
                    complex_type ExpCoeff[],
                    const bool FinalOut,
                    const int XYZ, const gramfe_fft_float MinDens,
                    forward_workspace_type Workspace,
                    inverse_workspace_type WorkspaceInv )
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
      complex_type* s_In_1PG = (complex_type*) gramfe_fftw::fft_malloc( GRAMFE_FLU_NXT*sizeof(complex_type) ); // allocate memory for fourier transform

      complex_type  s_Ae_1PG [CGPU_FLU_BLOCK_SIZE_Y][GRAMFE_NDELTA]; // left projection polynomials
      complex_type  s_Ao_1PG [CGPU_FLU_BLOCK_SIZE_Y][GRAMFE_NDELTA]; // right projection polynomials

//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_Ae, s_Ao )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif
      {

#        ifdef __CUDACC__
//       use two-dimensional thread blocks in GPU mode
         const uint tx = threadIdx.x;
         const uint ty = threadIdx.y;
#        else
//       every block just has a single thread with temporary memory on the stack in CPU mode
         const uint tx = 0;
         const uint ty = 0;

         s_In = (complex_type (*)[GRAMFE_FLU_NXT]) s_In_1PG;
         s_Ae = s_Ae_1PG;
         s_Ao = s_Ao_1PG;
#        endif

         const uint tid     = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
         const uint NThread = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block

         uint j, k;             // (i,j,k): array indices used in g_Fluid_In

         uint Column0 = 0;      // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;  // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out

         uint si, sj;           // array indices used in the shared memory array
         uint NStep;            // number of iterations for updating each column

         uint NColumnOnce = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration

//       use register variables Al, Ar and Psi_Ext to speed up summation
         complex_type Al, Ar;             // buffer for left and right Gram coefficients
         complex_type Psi_Ext;            // buffer for wave function in extension region

         real   Amp_New, Re_New, Im_New;  // store density, real and imaginary part to apply minimum density check

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

//          1.4 sync data read into S_In
#           ifdef __CUDACC__
            __syncthreads();
#           endif

//          2.1 compute Gram-Polynomial expansion coefficients via semi-discrete scalar products on boundary
            CELL_LOOP(GRAMFE_ORDER, 0, 0)
            {
               Al.real(0);
               Al.imag(0);
               Ar.real(0);
               Ar.imag(0);
               for (int t=0; t < GRAMFE_NDELTA; t++) {
                  Al += Pl[si][t] * s_In[sj][t];                           // left boundary
                  Ar += Pr[si][t] * s_In[sj][FLU_NXT - GRAMFE_NDELTA + t]; // right boundary
               } // for t

               s_Ae[sj][si] = (gramfe_fft_float) 0.5 * (Ar + Al);
               s_Ao[sj][si] = (gramfe_fft_float) 0.5 * (Ar - Al);
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif

//          2.2 function values in extension domain given as linear combinations of extended Gram polynomials
            CELL_LOOP(GRAMFE_FLU_NXT, FLU_NXT, 0)
            {
               Psi_Ext.real(0);
               Psi_Ext.imag(0);

               for (int order=0; order < GRAMFE_ORDER; order++) {
                  Psi_Ext += s_Ae[sj][order] *  Fe[order][si - FLU_NXT];
                  Psi_Ext += s_Ao[sj][order] *  Fo[order][si - FLU_NXT];
               } // for (int order=0; order < GRAMFE_ORDER; order++)

               s_In[sj][si] = Psi_Ext;
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif

//          3.1 forward FFT
#           ifdef __CUDACC__
            FFT().execute(reinterpret_cast<void*>(s_In), Workspace);
#           else
            gramfe_fftw_c2c( FFTW_Plan_ExtPsi, s_In );
#           endif

#           ifdef __CUDACC__
            __syncthreads();
#           endif

//          3.2 evolve wave function via time evolution operator with filter
            CELL_LOOP(GRAMFE_FLU_NXT, 0, 0)
            {
               s_In[sj][si] *= ExpCoeff[si];
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif


//          3.3 backward FFT
#           ifdef __CUDACC__
            IFFT().execute(reinterpret_cast<void*>(s_In), WorkspaceInv);
#           else
            gramfe_fftw_c2c( FFTW_Plan_ExtPsi_Inv, s_In );
#           endif

#           ifdef  __CUDACC__
            __syncthreads();
#           endif

//          4.1 write FFT array back to output array
            if ( FinalOut )
            {
               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx2 = get1D2( k, j, si, XYZ );

                  Amp_New = SQR(s_In[sj][si].real()) + SQR(s_In[sj][si].imag());
                  Re_New  = s_In[sj][si].real();
                  Im_New  = s_In[sj][si].imag();

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
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx1 = get1D1( k, j,si, XYZ );
                  g_Fluid_In[bx][0][Idx1]     = s_In[sj][si].real();
                  g_Fluid_In[bx][1][Idx1]     = s_In[sj][si].imag();
               }  // CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            } // if ( FinalOut ) ... else

#           ifdef  __CUDACC__
            __syncthreads();
#           endif

//          4.2 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      } // #pragma for (int bx=0; bx<NPatchGroup; bx++)
#     ifndef __CUDACC__
      gramfe_fftw::fft_free(s_In_1PG);
#     endif
   } // #pragma omp parallel

} // FUNCTION : CUFLU_Advance



#endif // #if ( GRAMFE_SCHEME == GRAMFE_FFT )
#endif // #if (  ( !defined(__CUDACC__) && defined(SUPPORT_FFTW) )  ||  defined(__CUDACC__)  )
