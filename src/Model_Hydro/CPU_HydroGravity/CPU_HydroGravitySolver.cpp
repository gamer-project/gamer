#include "CUPOT.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"


// variables reside in constant memory
__constant__ double c_ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_SetConstMem_HydroGravitySolver
// Description :  Set the constant memory used by CUPOT_HydroGravitySolver()
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Init_ExternalAccPot()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUPOT_SetConstMem_HydroGravitySolver( double h_ExtAcc_AuxArray[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( c_ExtAcc_AuxArray, h_ExtAcc_AuxArray, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUPOT_SetConstMem_HydroGravitySolver

#endif // ifdef __CUDACC__




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUPOT_HydroGravitySolver
// Description :  Advances the momentum and energy density of a group of patches by gravitational acceleration
//                (including external gravity)
//
// Note        :  1. Currently this function does NOT ensure the consistency between internal energy and
//                   dual-energy variable (e.g., entropy)
//                   --> This consistency breaks only for cells with the dual-energy status labelled
//                       as DE_UPDATED_BY_ETOT_GRA
//                   --> We restore this consistency in Gra_Close()
//                2. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_Flu_Array_New   : Array to store the input and output fluid variables
//                g_Pot_Array_New   : Array storing the input potential (at the current step)
//                                    --> _New: to be distinguishable from g_Pot_Array_USG[], which is defined at the previous step
//                g_Corner_Array    : Array storing the physical corner coordinates of each patch
//                g_Pot_Array_USG   : Array storing the input potential          for UNSPLIT_GRAVITY (at the previous step)
//                g_Flu_Array_USG   : Array storing the input density + momentum for UNSPLIT_GRAVITY (at the previous step)
//                g_DE_Array        : Array storing the dual-energy status (for both input and output)
//                g_EngyB_Array     : Array storing the cell-centered magnetic energy
//                                    --> Only for checking minimum internal energy in MHD
//                NPatchGroup       : Number of input patch groups (for CPU only)
//                dt                : Time interval to advance solution
//                dh                : Cell size
//                P5_Gradient       : Use 5-points stencil to evaluate the potential gradient
//                GravityType       : Types of gravity --> self-gravity, external gravity, both
//                c_ExtAcc_AuxArray : Auxiliary array for adding external acceleration (for CPU only)
//                                    --> When using GPU, this array is stored in the constant memory and does
//                                        not need to be passed as a function argument
//                                        --> Declared on top of this file with the prefix "c_" to
//                                            highlight that this is a constant variable on GPU
//                TimeNew           : Physical time at the current  step (for the external gravity solver)
//                TimeOld           : Physical time at the previous step (for the external gravity solver in UNSPLIT_GRAVITY)
//                MinEint           : Minimum allowed internal energy (== MIN_PRES / (GAMMA-1))
//
// Return      :  g_Flu_Array_New, g_DE_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUPOT_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real   g_EngyB_Array  [][ CUBE(PS1) ],
   const real dt, const real dh, const bool P5_Gradient,
   const OptGravityType_t GravityType,
   const double TimeNew, const double TimeOld, const real MinEint )
#else
void CPU_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real   g_EngyB_Array  [][ CUBE(PS1) ],
   const int NPatchGroup,
   const real dt, const real dh, const bool P5_Gradient,
   const OptGravityType_t GravityType, const double c_ExtAcc_AuxArray[],
   const double TimeNew, const double TimeOld, const real MinEint )
#endif
{

// check
#  ifdef GAMER_DEBUG
   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   if ( TimeNew < 0.0 )
      printf( "ERROR : incorrect TimeNew (%14.7e) !!\n", TimeNew );

#  ifdef UNSPLIT_GRAVITY
   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   if ( g_Pot_Array_USG == NULL  ||  g_Flu_Array_USG == NULL )
      printf( "ERROR : g_Pot_Array_USG == NULL  ||  g_Flu_Array_USG == NULL !!\n" );

   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   if ( TimeOld >= TimeNew  ||  TimeOld < 0.0 )
      printf( "ERROR : incorrect time (TimeOld %14.7e, TimeNew = %14.7e) !!\n", TimeOld, TimeNew );
#  endif

#  ifdef DUAL_ENERGY
   if ( g_DE_Array == NULL )
      printf( "ERROR : g_DE_Array == NULL !!\n" );
#  endif

#  ifdef MHD
   if ( g_EngyB_Array == NULL )
      printf( "ERROR : g_EngyB_Array == NULL !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG


   const real Gra_Const   = ( P5_Gradient ) ? -dt/(12.0*dh) : -dt/(2.0*dh);
   const int  PS1_sqr     = SQR(PS1);
   const int  didx_new[3] = { 1, GRA_NXT,   SQR(GRA_NXT)   };
#  ifdef UNSPLIT_GRAVITY
   const int  didx_old[3] = { 1, USG_NXT_G, SQR(USG_NXT_G) };
#  endif


// load potential from global to shared memory to improve the GPU performance
#  ifdef __CUDACC__
   __shared__ real s_pot_new[ CUBE(GRA_NXT) ];
#  ifdef UNSPLIT_GRAVITY
   __shared__ real s_pot_old[ CUBE(USG_NXT_G) ];
#  endif

   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   {
      for (int t=threadIdx.x; t<CUBE(GRA_NXT); t+=GRA_BLOCK_SIZE)
         s_pot_new[t] = g_Pot_Array_New[blockIdx.x][t];

#     ifdef UNSPLIT_GRAVITY
      for (int t=threadIdx.x; t<CUBE(USG_NXT_G); t+=GRA_BLOCK_SIZE)
         s_pot_old[t] = g_Pot_Array_USG[blockIdx.x][t];
#     endif
   }

   __syncthreads();
#  endif // #ifdef __CUDACC__


// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//     to work on different patches
#  ifdef __CUDACC__
   const int P = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int P=0; P<NPatchGroup*8; P++)
#  endif
   {
//    point to the potential array of the target patch
#     ifdef __CUDACC__
      const real *const pot_new = s_pot_new;
#     ifdef UNSPLIT_GRAVITY
      const real *const pot_old = s_pot_old;
#     endif
#     else // #ifdef __CUDACC__
      const real *const pot_new = g_Pot_Array_New[P];
#     ifdef UNSPLIT_GRAVITY
      const real *const pot_old = g_Pot_Array_USG[P];
#     endif
#     endif // #ifdef __CUDACC__ ... else ...


//    loop over all cells of the target patch
//    _g0: indices for the arrays without any ghost zone
      CGPU_LOOP( idx_g0, CUBE(PS1) )
      {
//       Enk = Etot - Ek = Eint in hydro and Eint+Eb in MHD
         real acc_new[3]={0.0, 0.0, 0.0}, px_new, py_new, pz_new, rho_new, Enk_in, Ek_out, Etot_in, Etot_out, _rho2;
#        ifdef UNSPLIT_GRAVITY
         real acc_old[3]={0.0, 0.0, 0.0}, px_old, py_old, pz_old, rho_old;
#        ifdef MHD
         real Eb_in;
#        endif
#        endif

         const int i_g0    = idx_g0 % PS1;
         const int j_g0    = idx_g0 % PS1_sqr / PS1;
         const int k_g0    = idx_g0 / PS1_sqr;

         const int i_new   = i_g0 + GRA_GHOST_SIZE;
         const int j_new   = j_g0 + GRA_GHOST_SIZE;
         const int k_new   = k_g0 + GRA_GHOST_SIZE;
         const int idx_new = IDX321( i_new, j_new, k_new, GRA_NXT, GRA_NXT );

#        ifdef UNSPLIT_GRAVITY
         const int i_old   = i_g0 + USG_GHOST_SIZE_G;
         const int j_old   = j_g0 + USG_GHOST_SIZE_G;
         const int k_old   = k_g0 + USG_GHOST_SIZE_G;
         const int idx_old = IDX321( i_old, j_old, k_old, USG_NXT_G, USG_NXT_G );
#        endif


//       external gravity
         if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
         {
            double x, y, z;

            x = g_Corner_Array[P][0] + (double)(i_g0*dh);
            y = g_Corner_Array[P][1] + (double)(j_g0*dh);
            z = g_Corner_Array[P][2] + (double)(k_g0*dh);

            ExternalAcc( acc_new, x, y, z, TimeNew, c_ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    acc_new[d] *= dt;

#           ifdef UNSPLIT_GRAVITY
            ExternalAcc( acc_old, x, y, z, TimeOld, c_ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    acc_old[d] *= dt;
#           endif
         }


//       self-gravity
         if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
         {
            const int ip1_new = idx_new + didx_new[0];
            const int jp1_new = idx_new + didx_new[1];
            const int kp1_new = idx_new + didx_new[2];
            const int im1_new = idx_new - didx_new[0];
            const int jm1_new = idx_new - didx_new[1];
            const int km1_new = idx_new - didx_new[2];

#           ifdef UNSPLIT_GRAVITY
            const int ip1_old = idx_old + didx_old[0];
            const int jp1_old = idx_old + didx_old[1];
            const int kp1_old = idx_old + didx_old[2];
            const int im1_old = idx_old - didx_old[0];
            const int jm1_old = idx_old - didx_old[1];
            const int km1_old = idx_old - didx_old[2];
#           endif

            if ( P5_Gradient )
            {
               const real Const_8 = (real)8.0;

               const int ip2_new = ip1_new + didx_new[0];
               const int jp2_new = jp1_new + didx_new[1];
               const int kp2_new = kp1_new + didx_new[2];
               const int im2_new = im1_new - didx_new[0];
               const int jm2_new = jm1_new - didx_new[1];
               const int km2_new = km1_new - didx_new[2];

#              ifdef UNSPLIT_GRAVITY
               const int ip2_old = ip1_old + didx_old[0];
               const int jp2_old = jp1_old + didx_old[1];
               const int kp2_old = kp1_old + didx_old[2];
               const int im2_old = im1_old - didx_old[0];
               const int jm2_old = jm1_old - didx_old[1];
               const int km2_old = km1_old - didx_old[2];
#              endif

               acc_new[0] += Gra_Const*( - pot_new[ip2_new] + Const_8*pot_new[ip1_new] - Const_8*pot_new[im1_new] + pot_new[im2_new] );
               acc_new[1] += Gra_Const*( - pot_new[jp2_new] + Const_8*pot_new[jp1_new] - Const_8*pot_new[jm1_new] + pot_new[jm2_new] );
               acc_new[2] += Gra_Const*( - pot_new[kp2_new] + Const_8*pot_new[kp1_new] - Const_8*pot_new[km1_new] + pot_new[km2_new] );

#              ifdef UNSPLIT_GRAVITY
               acc_old[0] += Gra_Const*( - pot_old[ip2_old] + Const_8*pot_old[ip1_old] - Const_8*pot_old[im1_old] + pot_old[im2_old] );
               acc_old[1] += Gra_Const*( - pot_old[jp2_old] + Const_8*pot_old[jp1_old] - Const_8*pot_old[jm1_old] + pot_old[jm2_old] );
               acc_old[2] += Gra_Const*( - pot_old[kp2_old] + Const_8*pot_old[kp1_old] - Const_8*pot_old[km1_old] + pot_old[km2_old] );
#              endif
            } // if ( P5_Gradient )

            else
            {
               acc_new[0] += Gra_Const*( pot_new[ip1_new] - pot_new[im1_new] );
               acc_new[1] += Gra_Const*( pot_new[jp1_new] - pot_new[jm1_new] );
               acc_new[2] += Gra_Const*( pot_new[kp1_new] - pot_new[km1_new] );

#              ifdef UNSPLIT_GRAVITY
               acc_old[0] += Gra_Const*( pot_old[ip1_old] - pot_old[im1_old] );
               acc_old[1] += Gra_Const*( pot_old[jp1_old] - pot_old[jm1_old] );
               acc_old[2] += Gra_Const*( pot_old[kp1_old] - pot_old[km1_old] );
#              endif
            } // if ( P5_Gradient ) ... else ...
         } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )


//       advance fluid
#        ifdef UNSPLIT_GRAVITY

         rho_new = g_Flu_Array_New[P][DENS][idx_g0];
         rho_old = g_Flu_Array_USG[P][DENS][idx_g0];
         px_new  = g_Flu_Array_New[P][MOMX][idx_g0];
         px_old  = g_Flu_Array_USG[P][MOMX][idx_g0];
         py_new  = g_Flu_Array_New[P][MOMY][idx_g0];
         py_old  = g_Flu_Array_USG[P][MOMY][idx_g0];
         pz_new  = g_Flu_Array_New[P][MOMZ][idx_g0];
         pz_old  = g_Flu_Array_USG[P][MOMZ][idx_g0];

//       backup the original internal energy so that we can restore it later if necessary
         _rho2   = (real)0.5/rho_new;
         Etot_in = g_Flu_Array_New[P][ENGY][idx_g0];
         Enk_in  = Etot_in - _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
#        ifdef MHD
         Eb_in   = g_EngyB_Array[P][idx_g0];
#        endif

//       update the momentum density
         px_new += (real)0.5*( rho_old*acc_old[0] + rho_new*acc_new[0] );
         py_new += (real)0.5*( rho_old*acc_old[1] + rho_new*acc_new[1] );
         pz_new += (real)0.5*( rho_old*acc_old[2] + rho_new*acc_new[2] );

         g_Flu_Array_New[P][MOMX][idx_g0] = px_new;
         g_Flu_Array_New[P][MOMY][idx_g0] = py_new;
         g_Flu_Array_New[P][MOMZ][idx_g0] = pz_new;

//       record the updated kinematic energy density
         Ek_out = _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );

//       update the total energy density
#        ifdef DUAL_ENERGY

//       for the unsplitting method with the dual-energy formalism, we correct the **total energy density**
//       only if the dual-energy status != DE_UPDATED_BY_DUAL
//       --> for (a) DE_UPDATED_BY_DUAL     --> Eint has been updated by the dual-energy variable
//               (b) DE_UPDATED_BY_MIN_PRES --> Eint has been set to the minimum threshold
//       --> currently for (b) we still update the total energy density

         if ( g_DE_Array[P][idx_g0] == DE_UPDATED_BY_DUAL )
         {
//          fix the internal energy and the dual-energy variable
            Etot_out = Enk_in + Ek_out;
         }

         else
         {
//          update the total energy, where internal energy and dual-energy variable may change as well
            Etot_out = Etot_in + (real)0.5*( px_old*acc_old[0] + py_old*acc_old[1] + pz_old*acc_old[2] +
                                             px_new*acc_new[0] + py_new*acc_new[1] + pz_new*acc_new[2] );

//          check the minimum internal energy
//          (a) if the updated internal energy is greater than the threshold, set the dual-energy status == DE_UPDATED_BY_ETOT_GRA
#           ifdef MHD
            if ( Etot_out - Ek_out - Eb_in >= MinEint )
#           else
            if ( Etot_out - Ek_out         >= MinEint )
#           endif
               g_DE_Array[P][idx_g0] = DE_UPDATED_BY_ETOT_GRA;

//          (b) otherwise restore the original internal energy and keep the original dual-energy status
            else
               Etot_out = Enk_in + Ek_out;
         }

#        else // # ifdef DUAL_ENERGY

//       for the unsplitting method without the dual-energy formalism, we always correct the total energy density
//       instead of the kinematic energy density
//       --> internal energy may change
//       --> we must check the minimum internal energy after this update
         Etot_out = Etot_in + (real)0.5*( px_old*acc_old[0] + py_old*acc_old[1] + pz_old*acc_old[2] +
                                          px_new*acc_new[0] + py_new*acc_new[1] + pz_new*acc_new[2] );

//       check the minimum internal energy
//       --> restore the original internal energy if the updated value becomes smaller than the threshold
#        ifdef MHD
         if ( Etot_out - Ek_out - Eb_in < MinEint )
#        else
         if ( Etot_out - Ek_out         < MinEint )
#        endif
            Etot_out = Enk_in + Ek_out;

#        endif // #ifdef DUAL_ENERGY ... else ...


#        else  // #ifdef UNSPLIT_GRAVITY


         rho_new = g_Flu_Array_New[P][DENS][idx_g0];
         px_new  = g_Flu_Array_New[P][MOMX][idx_g0];
         py_new  = g_Flu_Array_New[P][MOMY][idx_g0];
         pz_new  = g_Flu_Array_New[P][MOMZ][idx_g0];

//       backup the original internal energy so that we can restore it later
         _rho2   = (real)0.5/rho_new;
         Etot_in = g_Flu_Array_New[P][ENGY][idx_g0];
         Enk_in  = Etot_in - _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );

//       update the momentum density
         px_new += rho_new*acc_new[0];
         py_new += rho_new*acc_new[1];
         pz_new += rho_new*acc_new[2];

         g_Flu_Array_New[P][MOMX][idx_g0] = px_new;
         g_Flu_Array_New[P][MOMY][idx_g0] = py_new;
         g_Flu_Array_New[P][MOMZ][idx_g0] = pz_new;

//       for the splitting method, we ensure that the internal energy is unchanged
         Ek_out   = _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
         Etot_out = Enk_in + Ek_out;

#        endif // #ifdef UNSPLIT_GRAVITY ... else ...


//       store the updated total energy density to the output array
         g_Flu_Array_New[P][ENGY][idx_g0] = Etot_out;

      } // CGPU_LOOP( idx_g0, CUBE(PS1) )
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU/CUPOT_HydroGravitySolver



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
