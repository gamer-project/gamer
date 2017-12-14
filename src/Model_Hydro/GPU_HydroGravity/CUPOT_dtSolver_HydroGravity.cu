#include "Macro.h"
#include "CUPOT.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )



#include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"


// parallel reduction routine
#define RED_NTHREAD  ( PS1*PS1*DT_GRA_BLOCK_SIZE_Z )
#define RED_MAX

#ifdef DT_GRA_USE_SHUFFLE
#  include "../../GPU_Utility/CUUTI_BlockReduction_Shuffle.cu"
#else
#  include "../../GPU_Utility/CUUTI_BlockReduction_WarpSync.cu"
#endif


// variables reside in constant memory
__constant__ double ExtAcc_AuxArray_d[EXT_ACC_NAUX_MAX];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_dtSolver_HydroGravity_SetConstMem
// Description :  Set the constant memory used by CUPOT_dtSolver_HydroGravity
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUPOT_dtSolver_HydroGravity_SetConstMem( double ExtAcc_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtAcc_AuxArray_d, ExtAcc_AuxArray_h, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUPOT_dtSolver_HydroGravity_SetConstMem



//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_dtSolver_HydroGravity
// Description :  Estimate the evolution time-step (dt) required for the hydro gravity solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the free-fall time of the maximum gravitational acceleration
//
// Parameter   :  g_dt_Array     : Global memory array to store the minimum dt in each target patch
//                g_Pot_Array    : Global memory array storing the prepared potential data of each target patch
//                g_Corner_Array : Global memory array storing the physical corner coordinates of each patch
//                dh             : Grid size
//                Safety         : dt safety factor
//                P5_Gradient    : Use 5-points stencil to evaluate the potential gradient
//                GravityType    : Types of gravity --> self-gravity, external gravity, both
//                ExtAcc_Time    : Physical time for adding the external acceleration
//
// Return      :  dt_Array
//-----------------------------------------------------------------------------------------
__global__ void CUPOT_dtSolver_HydroGravity( real g_dt_Array[],
                                             const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                             const double g_Corner_Array[][3],
                                             const real dh, const real Safety, const bool P5_Gradient,
                                             const OptGravityType_t GravityType, const double ExtAcc_Time )
{

   const uint bx        = blockIdx.x;
   const uint tx        = threadIdx.x;
   const uint ty        = threadIdx.y;
   const uint tz        = threadIdx.z;
   const uint ID        = __umul24( tz, PS1*PS1 ) + __umul24( ty, PS1 ) + tx;
   const uint NSlice    = DT_GRA_BLOCK_SIZE_Z;
   const real Gra_Const = ( P5_Gradient ) ? (real)-1.0/((real)12.0*dh) : (real)-1.0/((real)2.0*dh);

   uint s_idx           =   __umul24( GRA_GHOST_SIZE+tz, GRA_NXT*GRA_NXT )
                          + __umul24( GRA_GHOST_SIZE+ty, GRA_NXT   ) + (GRA_GHOST_SIZE+tx);

   uint   ip1, jp1, kp1, im1, jm1, km1, t;
   uint   ip2, jp2, kp2, im2, jm2, km2;
   real   Acc[3], AccMax=(real)0.0;
   double x, y, z;

   __shared__ real s_Pot[ CUBE(GRA_NXT) ];


// set the physical coordinates of each cell for the external gravity solver
   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   {
      x = g_Corner_Array[bx][0] + (double)tx*dh;
      y = g_Corner_Array[bx][1] + (double)ty*dh;
      z = g_Corner_Array[bx][2] + (double)tz*dh;
   }


// load the potential from the global memory to the shared memory
   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   {
      t = ID;

      while ( t < CUBE(GRA_NXT) )
      {
         s_Pot[t] = g_Pot_Array[bx][t];
         t       += RED_NTHREAD;
      }
   }

   __syncthreads();


// loop over all z slices
   for (uint Slice=tz; Slice<PS1; Slice+=NSlice)
   {
      ip1 = s_idx + 1;
      jp1 = s_idx + GRA_NXT;
      kp1 = s_idx + GRA_NXT*GRA_NXT;
      im1 = s_idx - 1;
      jm1 = s_idx - GRA_NXT;
      km1 = s_idx - GRA_NXT*GRA_NXT;

      if ( P5_Gradient )
      {
         ip2 = s_idx + 2;
         jp2 = s_idx + 2*GRA_NXT;
         kp2 = s_idx + 2*GRA_NXT*GRA_NXT;
         im2 = s_idx - 2;
         jm2 = s_idx - 2*GRA_NXT;
         km2 = s_idx - 2*GRA_NXT*GRA_NXT;
      } // if ( P5_Gradient )


//    1. evaluate the gravitational acceleration
      Acc[0] = (real)0.0;
      Acc[1] = (real)0.0;
      Acc[2] = (real)0.0;

//    1.1 external gravity
      if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
         CUPOT_ExternalAcc( Acc, x, y, z, ExtAcc_Time, ExtAcc_AuxArray_d );

//    1.2 self-gravity
      if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
      {
         if ( P5_Gradient )   // 5-point gradient
         {
            Acc[0] += Gra_Const*( - s_Pot[ip2] + (real)8.0*s_Pot[ip1] - (real)8.0*s_Pot[im1] + s_Pot[im2] );
            Acc[1] += Gra_Const*( - s_Pot[jp2] + (real)8.0*s_Pot[jp1] - (real)8.0*s_Pot[jm1] + s_Pot[jm2] );
            Acc[2] += Gra_Const*( - s_Pot[kp2] + (real)8.0*s_Pot[kp1] - (real)8.0*s_Pot[km1] + s_Pot[km2] );
         }

         else                 // 3-point gradient
         {
            Acc[0] += Gra_Const*( s_Pot[ip1] - s_Pot[im1] );
            Acc[1] += Gra_Const*( s_Pot[jp1] - s_Pot[jm1] );
            Acc[2] += Gra_Const*( s_Pot[kp1] - s_Pot[km1] );
         }
      } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )


//    2. get the maximum acceleration
      AccMax = FMAX( AccMax, FABS(Acc[0]) );
      AccMax = FMAX( AccMax, FABS(Acc[1]) );
      AccMax = FMAX( AccMax, FABS(Acc[2]) );


//    3. update target cell indices
      s_idx += NSlice*SQR(GRA_NXT);
      z     += NSlice*dh;
   } // for (uint Slice=tz; Slice<PS1; Slice+=NSlice)


// 4. perform parallel reduction to get the maximum acceleration in each thread block
#  ifdef DT_GRA_USE_SHUFFLE
   AccMax = BlockReduction_Shuffle ( AccMax );
#  else
   AccMax = BlockReduction_WarpSync( AccMax );
#  endif


// 5. store the minimum dt in each patch back to the global memory
   if ( ID == 0 )    g_dt_Array[bx] = Safety*SQRT( (real)2.0*dh/AccMax );

} // FUNCTION : CUPOT_dtSolver_HydroGravity



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )
