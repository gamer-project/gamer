#include "CUFLU.h"
#include "GAMER.h"

#if ( ( !defined(__CUDACC__) && defined(SUPPORT_FFTW) ) || ( defined(__CUDACC__) && defined(GRAMFE_ENABLE_GPU) ) )

#if ( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE )
#include "GramExtensionTables.h"


// useful macros
# define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
# define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )


#ifndef __umul24
#  define __umul24( a, b )   ( (a)*(b) )
#endif
#ifndef __mul24
#  define  __mul24( a, b )   ( (a)*(b) )
#endif


// use cufftdx library for FFTs on GPU
#if __CUDACC__

// include cufftdx and define FFTs
#include <cufftdx.hpp>

using namespace cufftdx;

using fft_incomplete = decltype(Block() + Size<GRAMFE_FLU_NXT>() + Type<fft_type::c2c>() + Precision<gramfe_float>() + SM<750>());
using fft_base       = decltype(fft_incomplete() + Direction<fft_direction::forward>());
using ifft_base      = decltype(fft_incomplete() + Direction<fft_direction::inverse>());
static constexpr unsigned int elements_per_thread = use_suggested ? fft_base::elements_per_thread : custom_elements_per_thread;
static constexpr unsigned int ffts_per_block      = use_suggested ? fft_base::suggested_ffts_per_block : custom_ffts_per_block;


using FFT          = decltype(fft_base()  + ElementsPerThread<elements_per_thread>() + FFTsPerBlock<ffts_per_block>());
using IFFT         = decltype(ifft_base() + ElementsPerThread<elements_per_thread>() + FFTsPerBlock<ffts_per_block>());
using complex_type          = typename FFT::value_type;
using forward_workspace_type = typename FFT::workspace_type;
using inverse_workspace_type = typename IFFT::workspace_type;


//overload cufftdx's complex_type to allow for all the arithmetic operations required
template<class OtherType>
__device__ __forceinline__ complex_type operator*(const complex_type& a, const OtherType& other) {
      complex_type result(a);
      result *= other;
      return result;
}
template<class OtherType>
__device__ __forceinline__ complex_type operator*(const OtherType& other, const complex_type& a) {
      return a * other;
}

__device__ __forceinline__ complex_type operator*(const complex_type& a, const complex_type& b) {
      complex_type result(a);
      result *= b;
      return result;
}
__device__ __forceinline__ complex_type operator+(const complex_type& a, const complex_type& b) {
      complex_type result(a);
      result += b;
      return result;
}
__device__ __forceinline__ complex_type operator-(const complex_type& a, const complex_type& b) {
      complex_type result(a);
      result -= b;
      return result;
}

#else   // #if ( defined(__CUDACC__) && defined(GRAMFE_ENABLE_GPU) )

extern gramfe_complex_fftw_plan FFTW_Plan_ExtPsi, FFTW_Plan_ExtPsi_Inv;

#if ( SUPPORT_FFTW == FFTW3 )
#include <complex.h>

using complex_type = std::complex<gramfe_float>;

#else // #if ( SUPPORT_FFTW == FFTW3 )
//derive from gramfe_float_complex which is an alias for FFTW2's and FFTW3's complex_types in gramfe_float precision to allow for complex arithmetic operations
struct complex_type : public gramfe_float_complex {
   complex_type() {
      real(0);
      imag(0);
   }

   complex_type(gramfe_float re, gramfe_float im) {
      real(re);
      imag(im);
   }

   complex_type& operator=(const complex_type& other) {
      real(other.real());
      imag(other.imag());
      return *this;
   }

   gramfe_float real() const {
      return c_re(*this);
   }
   gramfe_float imag() const {
      return c_im(*this);
   }

   void real(gramfe_float re) {
      c_re(*this) = re;
   }
   void imag(gramfe_float im) {
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
};


template<class OtherType>
complex_type operator*(const complex_type& a, const OtherType& other) {
      return complex_type (a.real() * other, a.imag() * other);
}

template<class OtherType>
complex_type operator*(const OtherType& other, const complex_type& a) {
      return a * other;
}
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... # else


// no workspaces required in CPU solver
using forward_workspace_type = bool;
using inverse_workspace_type = bool;
#endif // #if ( defined(__CUDACC__) && defined(GRAMFE_ENABLE_GPU) )



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


# define CELL_LOOP( NCell, leftGhost, rightGhost )    for ( (Idx   = tid, (NStep = (NCell) - (leftGhost) - (rightGhost), (si = Idx % NStep + (leftGhost), sj = Idx / NStep))); \
                                                             Idx   < NColumnOnce * NStep; \
                                                            (Idx  += NThread, (si = Idx % NStep + (leftGhost), sj = Idx / NStep)) )


GPU_DEVICE
static void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                           real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                           real Flux_Array [][9][NFLUX_TOTAL][ SQR(PS2) ],
                           int NPatchGroup,
                           const gramfe_float dt, const gramfe_float _dh, const gramfe_float Eta, const bool StoreFlux,
                           const uint j_gap, const uint k_gap,
                           complex_type s_In  [][GRAMFE_FLU_NXT],
                           complex_type s_Al  [][GRAMFE_NDELTA],
                           complex_type s_Ar  [][GRAMFE_NDELTA],
                           complex_type ExpCoeff [],
                           const bool FinalOut, const int XYZ, const gramfe_float MinDens,
                           forward_workspace_type workspace,
                           inverse_workspace_type workspace_inverse );

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
//                Flux_Array     : Array to store the output flux
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                                 --> Meaningless if CONSERVE_MASS is off since the operators along different directions
//                                     commute
//                                 --> Meaningful if CONSERVE_MASS is on, in which the symmetry along different directions
//                                     are broken ...
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------

GPU_DEVICE
int fact(int n) {
     return (n==0) || (n==1) ? 1 : n* fact(n-1);
}
GPU_DEVICE
gramfe_float taylor_fact(int n) {
   return (1 / ((gramfe_float) fact(n)));
}

GPU_DEVICE
gramfe_float cosine_series(gramfe_float x, int Nterms) {
   gramfe_float result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * taylor_fact(2 * i   ) * pow(x, 2 * i   );
   }

   return result;
}

GPU_DEVICE
gramfe_float sine_series(gramfe_float x, int Nterms) {
   gramfe_float result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * taylor_fact(2 * i + 1) * pow(x, 2 * i + 1);
   }

   return result;
}

#ifdef __CUDACC__
__launch_bounds__(FFT::max_threads_per_block)
__global__
void CUFLU_ELBDMSolver_GramFE(    real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                  real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                  real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                  const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                  const bool XYZ, const real MinDens,
                                  typename FFT::workspace_type workspace,
                                  typename IFFT::workspace_type workspace_inverse  )
#else
void CPU_ELBDMSolver_GramFE(      real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                  real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                  real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                  const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                                  const bool XYZ, const real MinDens )
#endif
{

#  ifdef __CUDACC__
    // Execute FFT
    extern __shared__ complex_type shared_mem[];

// create memories for columns of various intermediate fields in shared GPU memory
   complex_type (*s_In)[GRAMFE_FLU_NXT]     = (complex_type (*)[GRAMFE_FLU_NXT]) (shared_mem);
   complex_type (*s_Ae)[GRAMFE_NDELTA]      = (complex_type (*)[GRAMFE_NDELTA])  (shared_mem + CGPU_FLU_BLOCK_SIZE_Y * (GRAMFE_FLU_NXT                   ));    // 0.5 * log(rho)
   complex_type (*s_Ao)[GRAMFE_NDELTA]      = (complex_type (*)[GRAMFE_NDELTA])  (shared_mem + CGPU_FLU_BLOCK_SIZE_Y * (GRAMFE_FLU_NXT + GRAMFE_NDELTA));    // the fluxes for every thread block

#  ifdef CONSERVE_MASS
   __shared__ real s_Flux    [CGPU_FLU_BLOCK_SIZE_Y][GRAMFE_FLU_NXT];
#  else  // #  ifdef CONSERVE_MASS
              real (*s_Flux)[GRAMFE_FLU_NXT] = NULL;  // useless if CONSERVE_MASS is off
#  endif // #  ifdef CONSERVE_MASS ... # else

   const int NPatchGroup = 0;
   const gramfe_float dh                            = gramfe_float(1.0)/_dh;

#  else // #  ifdef __CUDACC__
// allocate memory on stack within loop for CPU run
   complex_type (*s_In)   [GRAMFE_FLU_NXT]          = NULL;
   complex_type (*s_Ae)   [GRAMFE_NDELTA]           = NULL;
   complex_type (*s_Ao)   [GRAMFE_NDELTA]           = NULL;
   const gramfe_float _dh                           = gramfe_float(1.0)/dh;
   bool workspace, workspace_inverse;
#  endif // #  ifdef __CUDACC__ ... else


   const gramfe_float filterDecay      = 32.0 * 2.302585f;   // decay of k-space filter
   const gramfe_float filterDegree     = 100;                // degree of k-space filter

// set up wave numbers
// and the k-space evolution phase k**2*dt/(2*ELBDM_ETA)
   gramfe_float K;
   gramfe_float Filter;                               // exp(-filterDecay * (k/kMax)**(2*filterDegree))
   gramfe_float Coeff;                                // dT * k^2
   complex_type ExpCoeff[GRAMFE_FLU_NXT];             // exp(1j * dT * k^2)

   const gramfe_float kmax     = M_PI*_dh;

   const gramfe_float dT       = -(gramfe_float)0.5*dt/Eta;                     // coefficient in time evolution operator

   const gramfe_float Norm     = 1.0 / ( (gramfe_float)GRAMFE_FLU_NXT );

   for (int i=0; i<GRAMFE_FLU_NXT; i++)
   {
//    set up momentum array
      K           = ( i <= GRAMFE_FLU_NXT/2 ) ? 2.0*M_PI/(GRAMFE_FLU_NXT*dh)*i : 2.0*M_PI/(GRAMFE_FLU_NXT*dh)*(i-GRAMFE_FLU_NXT);
//    set up filter array
      Filter      = EXP(-filterDecay * POW(FABS(K/kmax), 2*filterDegree));
//    set up time evolution array
      Coeff       = SQR(K)*dT;

      //ExpCoeff[i] = complex_type(COS(Coeff), SIN(Coeff)) * Norm * Filter;

//    ideal Taylor expansion order == FLU_GHOST_SIZE - 1
//    for FLU_GHOST_SIZE == 8, an order of 7 provides optimal results
      ExpCoeff[i] = complex_type(cosine_series(Coeff, 4), sine_series(Coeff, 3)) * Norm * Filter;
   }


   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Ae, s_Ao, ExpCoeff, false, 0, MinDens, workspace, workspace_inverse );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE,              0, s_In, s_Ae, s_Ao, ExpCoeff, false, 3, MinDens, workspace, workspace_inverse );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Ae, s_Ao, ExpCoeff, true, 6, MinDens, workspace, workspace_inverse );
   } else  {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Ae, s_Ao, ExpCoeff, false, 6, MinDens, workspace, workspace_inverse );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0, FLU_GHOST_SIZE, s_In, s_Ae, s_Ao, ExpCoeff, false, 3, MinDens, workspace, workspace_inverse );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Ae, s_Ao, ExpCoeff, true, 0, MinDens, workspace, workspace_inverse );
   }

} // FUNCTION : CUFLU_ELBDMSolver_GramFE



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_AdvanceX
// Description :  Use CPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on Gram-Fourier extension with pseudo-spectral solver on extended domain
//
// Parameter   :  u              : Array storing the input variables (only REAL/IMAG)
//                Flux_Array     : Array to store the output flux (only density)
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion
//                j_gap          : Number of cells to be skipped on each side in the y direction
//                k_gap          : Number of cells to be skipped on each side in the z direction
//                Flux_XYZ       : Parameter used to determine the place to store the output fluxes
//                                 --> (0,3,6) <-> (x/y/z) fluxes
//                                 --> useful only if CONSERVE_MASS is defined
//-------------------------------------------------------------------------------------------------------

GPU_DEVICE
void CUFLU_Advance(  real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                     real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                     real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                     int NPatchGroup,
                     const gramfe_float dt, const gramfe_float _dh, const gramfe_float Eta, const bool StoreFlux,
                     const uint j_gap, const uint k_gap,
                     complex_type s_In    [][GRAMFE_FLU_NXT],
                     complex_type s_Ae    [][GRAMFE_NDELTA],
                     complex_type s_Ao    [][GRAMFE_NDELTA],
                     complex_type ExpCoeff [],
                     const bool FinalOut,
                     const int XYZ, const gramfe_float MinDens,
                     forward_workspace_type workspace,
                     inverse_workspace_type workspace_inverse  )

{
   const uint j_end        = FLU_NXT -  j_gap    ;          // last y-column to be updated
   const uint size_j       = FLU_NXT -  j_gap * 2;          // number of y-columns to be updated
   const uint size_k       = FLU_NXT -  k_gap * 2;          // number of z-columns to be updated
   const uint NColumnTotal = __umul24( size_j, size_k );    // total number of data columns to be updated
   int Idx;

// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
#     ifdef __CUDACC__
      const int bx = blockIdx.x;
#     else

//    create arrays for columns of various intermediate fields on the stack
      complex_type* s_In_1PG = (complex_type* ) gramfe_fftw_malloc( GRAMFE_FLU_NXT * sizeof(complex_type) ); // allocate memory for fourier transform

      complex_type s_Ae_1PG     [CGPU_FLU_BLOCK_SIZE_Y][GRAMFE_NDELTA]; // left projection polynomials
      complex_type s_Ao_1PG     [CGPU_FLU_BLOCK_SIZE_Y][GRAMFE_NDELTA]; // right projection polynomials

//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_Ae, s_Ao )
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


         s_In     = (complex_type (*)[GRAMFE_FLU_NXT]) s_In_1PG;
         s_Ae     = s_Ae_1PG;
         s_Ao     = s_Ao_1PG;

#        endif // # ifdef __CUDACC__ ... # else

         const uint tid          = __umul24( ty , CGPU_FLU_BLOCK_SIZE_X ) + tx;    // thread ID within block
         const uint NThread      = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block

         uint j, k;                         // (i,j,k): array indices used in g_Fluid_In

         uint Column0 = 0;                // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;            // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out

         uint si, sj;                     // array indices used in the shared memory array
         uint NStep;                      // number of iterations for updating each column

         uint NColumnOnce        = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration

         complex_type Al, Ar;

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
            } // CELL_LOOP(FLU_NXT, 0, 0) {

//          1.4 sync data read into S_In
#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          2.1 compute Gram-Polynomial expansion coefficients via semi-discrete scalar products on boundary
            CELL_LOOP(GRAMFE_ORDER, 0, 0)
            {
               Al.real(0);
               Al.imag(0);
               Ar.real(0);
               Ar.imag(0);
               for (int t=0; t < GRAMFE_NDELTA; t++) {
                  Al += Pl[si][t] * s_In[sj][t];                              // left boundary
                  Ar += Pr[si][t] * s_In[sj][FLU_NXT - GRAMFE_NDELTA + t]; // right boundary
               } // for t

               s_Ae[sj][si] = (gramfe_float) 0.5 * (Ar + Al);
               s_Ao[sj][si] = (gramfe_float) 0.5 * (Ar - Al);
            } // CELL_LOOP(GRAMFE_ORDER, 0, 0)

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          2.2 function values in extension domain given as linear combinations of extended Gram polynomials
            CELL_LOOP(GRAMFE_FLU_NXT, FLU_NXT, 0)
            {
               s_In[sj][si].real(0);
               s_In[sj][si].imag(0);


               for (int order=0; order < GRAMFE_ORDER; order++) {
                  s_In[sj][si] += s_Ae[sj][order] *  Fe[order][si - FLU_NXT + GRAMFE_NDELTA];
                  s_In[sj][si] += s_Ao[sj][order] *  Fo[order][si - FLU_NXT + GRAMFE_NDELTA];
               } // for (int order=0; order < GRAMFE_ORDER; order++)
            } // CELL_LOOP(GRAMFE_FLU_NXT, FLU_NXT, 0)

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          3.1 forward FFT
#           ifdef __CUDACC__
            FFT().execute(reinterpret_cast<void*>(s_In), workspace);//, shared_memory, workspace);
#           else // # ifdef __CUDACC_
            gramfe_fftw_c2c( FFTW_Plan_ExtPsi, s_In );
#           endif // # ifdef __CUDACC_ ... # else

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          3.2 evolve wave function via time evolution operator with filter
            CELL_LOOP(GRAMFE_FLU_NXT, 0, 0)
            {
               s_In[sj][si] *= ExpCoeff[si];
            } // CELL_LOOP(GRAMFE_FLU_NXT, 0, 0)

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_


//          3.3 backward FFT
#           ifdef __CUDACC__
            IFFT().execute(reinterpret_cast<void*>(s_In), workspace_inverse);//, shared_memory, workspace);
#           else // # ifdef __CUDACC_
            gramfe_fftw_c2c( FFTW_Plan_ExtPsi_Inv, s_In );
#           endif // # ifdef __CUDACC_ ... # else

#           ifdef  __CUDACC__
            __syncthreads();
#           endif

//          3.4 write FFT array back to output array
            if ( FinalOut )
            {
               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx2 = get1D2( k, j, si, XYZ);
                  g_Fluid_Out[bx][DENS][Idx2] = SQR(s_In[sj][si].real()) + SQR(s_In[sj][si].imag());
                  g_Fluid_Out[bx][REAL][Idx2] = s_In[sj][si].real();
                  g_Fluid_Out[bx][IMAG][Idx2] = s_In[sj][si].imag();
               } // CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            } else { // if ( FinalOut )
               CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
               {
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx1 = get1D1( k, j,si, XYZ);
                  g_Fluid_In[bx][0][Idx1]     = s_In[sj][si].real();
                  g_Fluid_In[bx][1][Idx1]     = s_In[sj][si].imag();
               }  // CELL_LOOP(FLU_NXT, FLU_GHOST_SIZE, FLU_GHOST_SIZE)
            } // if ( FinalOut ) ... else#

#           ifdef  __CUDACC__
            __syncthreads();
#           endif


//          3.7 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      } // # pragma  for (int bx=0; bx<NPatchGroup; bx++)
#     ifndef __CUDACC__
      gramfe_fftw_free(s_In_1PG);
#     endif
   } // # pragma omp parallel
} // FUNCTION : CUFLU_Advance


#endif // #if ( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE)
#endif // #if ( ( !defined(__CUDACC__) && defined(SUPPORT_FFTW) ) || ( defined(__CUDACC__) && defined(GRAMFE_ENABLE_GPU) ) )