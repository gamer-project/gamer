#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == WAF )


// check before compiling anything else
#if ( NCOMP_PASSIVE != 0 )
#  error : WAF scheme does NOT support passive scalars !!
#endif


#define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
#define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )

#include "CUFLU_Shared_FluUtility.cu"
#if(    RSOLVER == EXACT )
#include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( RSOLVER == ROE )
static __device__ void Solve_StarRoe( real eival[5], real L_star[5], real R_star[5], const real L[5],
                                      const real R[5], const real Gamma, const real MinPres );
#endif
#ifdef WAF_DISSIPATE
static __device__ void Dis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5],
                                 const real R_star[5], const real limit[5], const real theta[5],
                                 const real Gamma );
#else
static __device__ void Undis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5],
                                   const real R_star[5], const real limit[5], const real theta[5],
                                   const real Gamma );
#endif
static __device__ void CUFLU_Advance( real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                      real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                      real g_Flux[][9][5][ PS2*PS2 ],
                                      const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                      const int j_gap, const int k_gap, real s_u[][FLU_NXT][5],
                                      real s_flux[][PS2+1][5], real s_Lstar[][PS2+3][5], real s_Rstar[][PS2+3][5],
                                      const bool FinalOut, const int XYZ, const WAF_Limiter_t WAF_Limiter,
                                      const real MinDens, const real MinPres );
static __device__ void Solve_Flux( real flux[5], const real lL_star[5], const real lR_star[5],
                                   const real cL_star[5], const real cR_star[5], const real rL_star[5],
                                   const real rR_star[5], const real eival[5] ,const real L_2[5],
                                   const real L_1[5], const real R_1[5],const real R_2[5],
                                   const real Gamma, const real ratio, const WAF_Limiter_t WAF_Limiter );
static __device__ void set_flux( real flux[5], const real val[5], const real Gamma );
static __device__ real set_limit( const real r, const real c, const WAF_Limiter_t WAF_Limiter );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_WAF
// Description :  GPU fluid solver based on the Weighted-Average-Flux (WAF) scheme
//
// Note        :  a. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                b. The three-dimensional evolution is achieved by using the dimensional-split method
//
// Parameter   :  g_Fluid_In   : Global memory array to store the input fluid variables
//                g_Fluid_Out  : Global memory array to store the output fluid variables
//                g_Flux       : Global memory array to store the output flux
//                g_Corner     : Global memory array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                g_Pot_USG    : Global memory array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in RTVD)
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Gamma        : Ratio of specific heats
//                StoreFlux    : true --> store the coarse-fine fluxes
//                XYZ          : true  : x->y->z ( forward sweep)
//                               false : z->y->x (backward sweep)
//                WAF_Limiter  : Selection of the limit function
//                                 0 : superbee
//                                 1 : van-Leer
//                                 2 : van-Albada
//                                 3 : minbee
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_FluidSolver_WAF( real g_Fluid_In []   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Fluid_Out[]   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
                                       real g_Flux     [][9][NCOMP_TOTAL][ PS2*PS2 ],
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[][ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const bool XYZ, const WAF_Limiter_t WAF_Limiter,
                                       const real MinDens, const real MinPres )
{

   __shared__ real s_u     [FLU_BLOCK_SIZE_Y][FLU_NXT][5];
   __shared__ real s_flux  [FLU_BLOCK_SIZE_Y][PS2+1][5];
   __shared__ real s_L_st  [FLU_BLOCK_SIZE_Y][PS2+3][5];
   __shared__ real s_R_st  [FLU_BLOCK_SIZE_Y][PS2+3][5];

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux,              0,              0,
                     s_u, s_flux, s_L_st, s_R_st, false, 0, WAF_Limiter, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux, FLU_GHOST_SIZE,              0,
                     s_u, s_flux, s_L_st, s_R_st, false, 3, WAF_Limiter, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE,
                     s_u, s_flux, s_L_st, s_R_st,  true, 6, WAF_Limiter, MinDens, MinPres );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux,              0,              0,
                     s_u, s_flux, s_L_st, s_R_st, false, 6, WAF_Limiter, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux,              0, FLU_GHOST_SIZE,
                     s_u, s_flux, s_L_st, s_R_st, false, 3, WAF_Limiter, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE,
                     s_u, s_flux, s_L_st, s_R_st,  true, 0, WAF_Limiter, MinDens, MinPres );
   }

} // FUNCTION : CUFLU_FluidSolver_WAF



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  GPU device function, which performs a one-dimensional sweep based on the WAF scheme
//
// Note        :  a. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                b. The direction of the one dimensional sweep is determined by the input parameter "XYZ"
//
// Parameter   :  g_Fluid_In   : Global memory array to store the input fluid variables
//                g_Fluid_Out  : Global memory array to store the output fluid variables
//                g_Flux       : Global memory array to store the output flux
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Gamma        : ratio of specific heats
//                StoreFlux    : true --> store the coarse-fine fluxes
//                j_gap        : Number of useless grids in each side in the j direction (j may not be equal to y)
//                k_gap        : Number of useless grids in each side in the k direction (k mya not be equal to z)
//                s_u          : Shared memory array storing the fluid variables used to compute the intercell flux
//                s_flux       : Shared memory array storing the final flux used to update the fluid variables
//                s_Lstar      : Shared memory array storing the left region in the solution of Riemann problem
//                s_Rstar      : Shared memory array storing the right region in the solution of Riemann problem
//                FinalOut     : true  : output data
//                               false : don't output data
//                XYZ          : 0 : Update the solution in the x direction
//                               3 : Update the solution in the y direction
//                               6 : Update the solution in the z direction
//                               --> This parameter is also used to determine the place to store the output flux
//                WAF_Limiter  : Selection of the limit function
//                                0 : superbee
//                                1 : van-Leer
//                                2 : van-Albada
//                                3 : minbee
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_Advance( real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                               real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                               real g_Flux[][9][5][ PS2*PS2 ],
                               const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                               const int j_gap, const int k_gap, real s_u[][FLU_NXT][5], real s_flux[][PS2+1][5],
                               real s_Lstar[][PS2+3][5], real s_Rstar[][PS2+3][5], const bool FinalOut,
                               const int XYZ, const WAF_Limiter_t WAF_Limiter, const real MinDens, const real MinPres )
{

   const uint bx        = blockIdx.x;
   const uint tx        = threadIdx.x;
   const uint ty        = threadIdx.y;
   const uint dj        = blockDim.y;
   const uint size_j    = FLU_NXT - (j_gap<<1);
   const uint size_k    = FLU_NXT - (k_gap<<1);
   const uint NColumn   = __umul24( size_j, size_k );
   const uint i         = tx;                   // (i,j) the element in shared memory under evaluation
         uint j         = j_gap + ty%size_j;
         uint k         = k_gap + ty/size_j;
         uint Column0   = 0;                    // the total number of columns that have been updated
   const uint j_end     = FLU_NXT - j_gap;
   const uint k_end     = FLU_NXT - k_gap;
   const real ratio     = dt*_dh;               // dt over dx
   const real Gamma_m1  = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
         bool RuleOut   = false;

   real   Fluid[5], eval[5];
   int    ID1, ID2, ID3, ii, delta_k, Comp[5];
   FluVar ConVar;


// set the order of component for update in different directions
   switch ( XYZ )
   {
      case 0:  Comp[0] = 0;   Comp[1] = 1;   Comp[2] = 2;   Comp[3] = 3;   Comp[4] = 4;   break;
      case 3:  Comp[0] = 0;   Comp[1] = 2;   Comp[2] = 1;   Comp[3] = 3;   Comp[4] = 4;   break;
      case 6:  Comp[0] = 0;   Comp[1] = 3;   Comp[2] = 2;   Comp[3] = 1;   Comp[4] = 4;   break;
   }


// start the WAF scheme
   do
   {
//    determine the array indices for updating in different directions
      switch ( XYZ )
      {
         case 0:  ID1 = to1D1( k, j, i );    break;
         case 3:  ID1 = to1D1( k, i, j );    break;
         case 6:  ID1 = to1D1( i, k, j );    break;
      }


//    load data into per-thread registers
      for (int v=0; v<5; v++)    Fluid[v] = g_Fluid_In[bx][ Comp[v] ][ID1];


//    load the primitive variables into shared memory
      s_u[ty][i][0] = Fluid[0];
      s_u[ty][i][1] = Fluid[1] / Fluid[0];
      s_u[ty][i][2] = Fluid[2] / Fluid[0];
      s_u[ty][i][3] = Fluid[3] / Fluid[0];
      s_u[ty][i][4] = Gamma_m1*( Fluid[4] - (real)0.5*( Fluid[1]*Fluid[1] + Fluid[2]*Fluid[2] +
                                                        Fluid[3]*Fluid[3] ) / Fluid[0] );
      s_u[ty][i][4] = CUFLU_CheckMinPres( s_u[ty][i][4], MinPres );

      __syncthreads();


//    solve the Riemann problem
      if ( i >= 1  &&  i <= FLU_GHOST_SIZE + PS2 + 1 )
      {
         ii = i - 1;

#        if ( RSOLVER == EXACT )
         FluVar5 eival_st, L_star_st, R_star_st;
         FluVar  L_st, R_st;

         L_st.Rho = s_u[ty][ii][0];
         L_st.Px  = s_u[ty][ii][1];
         L_st.Py  = s_u[ty][ii][2];
         L_st.Pz  = s_u[ty][ii][3];
         L_st.Egy = s_u[ty][ii][4];

         R_st.Rho = s_u[ty][ i][0];
         R_st.Px  = s_u[ty][ i][1];
         R_st.Py  = s_u[ty][ i][2];
         R_st.Pz  = s_u[ty][ i][3];
         R_st.Egy = s_u[ty][ i][4];

         CUFLU_RiemannSolver_Exact( 0, &eival_st, &L_star_st, &R_star_st, L_st, R_st, Gamma );

         eval[0] = eival_st.Rho;
         eval[1] = eival_st.Px;
         eval[2] = eival_st.Py;
         eval[3] = eival_st.Pz;
         eval[4] = eival_st.Egy;

         s_Lstar[ty][ii][0] = L_star_st.Rho;
         s_Lstar[ty][ii][1] = L_star_st.Px;
         s_Lstar[ty][ii][2] = L_star_st.Py;
         s_Lstar[ty][ii][3] = L_star_st.Pz;
         s_Lstar[ty][ii][4] = L_star_st.Egy;

         s_Rstar[ty][ii][0] = R_star_st.Rho;
         s_Rstar[ty][ii][1] = R_star_st.Px;
         s_Rstar[ty][ii][2] = R_star_st.Py;
         s_Rstar[ty][ii][3] = R_star_st.Pz;
         s_Rstar[ty][ii][4] = R_star_st.Egy;

#        elif ( RSOLVER == ROE )
         Solve_StarRoe( eval, s_Lstar[ty][ii], s_Rstar[ty][ii], s_u[ty][ii], s_u[ty][i], Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
#        endif
      }

      __syncthreads();


//    solve the intercell flux
      if ( i >= FLU_GHOST_SIZE  &&  i <= FLU_GHOST_SIZE+PS2  )
      {
         ii = i - FLU_GHOST_SIZE;
         int ii_p1 = ii + 1;

         Solve_Flux( s_flux[ty][ii], s_Lstar[ty][ii], s_Rstar[ty][ii],
                     s_Lstar[ty][ii_p1], s_Rstar[ty][ii_p1], s_Lstar[ty][i], s_Rstar[ty][i], eval,
                     s_u[ty][i-2], s_u[ty][i-1], s_u[ty][i], s_u[ty][i+1],
                     Gamma, ratio, WAF_Limiter );
      }

      __syncthreads();


//    update the conservative variables
      if ( i >= FLU_GHOST_SIZE  &&  i < FLU_GHOST_SIZE+PS2  &&  RuleOut == false )
      {
         ii = i - FLU_GHOST_SIZE;
         for (int v=0; v<5; v++)    Fluid[v] += ratio*( s_flux[ty][ii][v] - s_flux[ty][ii+1][v] );

//       enforce positive density and pressure
         ConVar.Rho = Fluid[0];
         ConVar.Px  = Fluid[1];
         ConVar.Py  = Fluid[2];
         ConVar.Pz  = Fluid[3];
         ConVar.Egy = Fluid[4];

         ConVar.Rho = FMAX( ConVar.Rho, MinDens );
         Fluid[0]   = ConVar.Rho;
         Fluid[4]   = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


//       check negative density and energy
#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( CUFLU_CheckNegative(Fluid[0]) )
            printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid[0], __FILE__, __LINE__, __FUNCTION__ );

         if ( CUFLU_CheckNegative(Fluid[4]) )
            printf( "ERROR : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid[4], __FILE__, __LINE__, __FUNCTION__ );
#        endif


//       store the updated data back to the global memory
         if ( FinalOut )
         {
            switch ( XYZ )
            {
               case 0:  ID2 = to1D2( k, j, i );    break;
               case 3:  ID2 = to1D2( k, i, j );    break;
               case 6:  ID2 = to1D2( i, k, j );    break;
            }

            for (int v=0; v<5; v++)    g_Fluid_Out[bx][ Comp[v] ][ID2] = Fluid[v];
         }

         else
            for (int v=0; v<5; v++)    g_Fluid_In [bx][ Comp[v] ][ID1] = Fluid[v];
      }


//    paste the s_flux into g_Flux
      if ( StoreFlux )
      if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
      if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
      if ( i == 0 )
      {
         ID3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

         for (int v=0; v<5; v++)
         {
            g_Flux[bx][XYZ+0][v][ID3] = s_flux[ty][     0 ][ Comp[v] ];
            g_Flux[bx][XYZ+1][v][ID3] = s_flux[ty][ PS2/2 ][ Comp[v] ];
            g_Flux[bx][XYZ+2][v][ID3] = s_flux[ty][ PS2   ][ Comp[v] ];
         }
      }


//    reset the target array indices
      j += dj;

      if ( j >= j_end )
      {
         delta_k  = ( j - j_end )/size_j + 1;
         k       += delta_k;
         j       -= __umul24( size_j, delta_k );
      }

      Column0 += dj;

//    if the index k exceeds the maximum allowed value --> reset (j,k) to harmless values and wait for other
//    threads (all threads must exist the while loop "at the same time", otherwise __syncthreads will fail !!)
      if ( k >= k_end )
      {
         j       = 0;
         k       = 0;
         RuleOut = true;
      }

      __syncthreads();
   }
   while ( Column0 < NColumn );

} // CUFLU_Advance



//-------------------------------------------------------------------------------------------------------
// Function    :  Solve_Flux
// Description :  Solve the intercell flux
//
// Parameter   :  flux        : Intercell flux
//                lL_star     : Primitive variables in the left star region of the left region
//                lR_star     : Primitive variables in the right star region of the left region
//                cL_star     : Primitive variables in the left star region of the centor region
//                cR_star     : Primitive variables in the right star region of the centor region
//                rL_star     : Primitive variables in the left star region of the right region
//                rR_star     : Primitive variables in the right star region of the right region
//                eival       : Eigenvalue
//                L_2         : Primitive variables in the region left to the left region
//                L_1         : Primitive variables in the left region
//                R_1         : Primitive variables in the right region
//                R_2         : Primitive variables in the region right to the right region
//                ratio       : dt over dx
//                Gamma       : Ratio of specific heats
//                WAF_Limiter : Selection of te limit function
//                                 0 : superbee
//                                 1 : van-Leer
//                                 2 : van-Albada
//                                 3 : minbee
//-------------------------------------------------------------------------------------------------------
__device__ void Solve_Flux( real flux[5], const real lL_star[5], const real lR_star[5],
                            const real cL_star[5], const real cR_star[5], const real rL_star[5],
                            const real rR_star[5], const real eival[5], const real L_2[5], const real L_1[5],
                            const real R_1[5],const real R_2[5], const real Gamma, const real ratio,
                            const WAF_Limiter_t WAF_Limiter )
{

   real theta[5];            // the sign of speed of waves
   real limit[5];            // limit functions
   real mean [3][5];
   real delta[3][5];

   delta[0][0] = lL_star[0] -     L_2[0];
   delta[0][1] = lR_star[0] - lL_star[0];
   delta[0][2] = lR_star[2] - lL_star[2];
   delta[0][3] = lR_star[3] - lL_star[3];
   delta[0][4] =     L_1[0] - lR_star[0];
   mean[0][0] = (real)0.5*( FABS( lL_star[0] ) + FABS(     L_2[0] ) );
   mean[0][1] = (real)0.5*( FABS( lR_star[0] ) + FABS( lL_star[0] ) );
   mean[0][2] = (real)0.5*( FABS( lR_star[2] ) + FABS( lL_star[2] ) );
   mean[0][3] = (real)0.5*( FABS( lR_star[3] ) + FABS( lL_star[3] ) );
   mean[0][4] = (real)0.5*( FABS(     L_1[0] ) + FABS( lR_star[0] ) );

   delta[1][0] = cL_star[0] -     L_1[0];
   delta[1][1] = cR_star[0] - cL_star[0];
   delta[1][2] = cR_star[2] - cL_star[2];
   delta[1][3] = cR_star[3] - cL_star[3];
   delta[1][4] =     R_1[0] - cR_star[0];
   mean[1][0] = (real)0.5*( FABS( cL_star[0] ) + FABS(     L_1[0] ) );
   mean[1][1] = (real)0.5*( FABS( cR_star[0] ) + FABS( cL_star[0] ) );
   mean[1][2] = (real)0.5*( FABS( cR_star[2] ) + FABS( cL_star[2] ) );
   mean[1][3] = (real)0.5*( FABS( cR_star[3] ) + FABS( cL_star[3] ) );
   mean[1][4] = (real)0.5*( FABS(     R_1[0] ) + FABS( cR_star[0] ) );

   delta[2][0] = rL_star[0] -     R_1[0];
   delta[2][1] = rR_star[0] - rL_star[0];
   delta[2][2] = rR_star[2] - rL_star[2];
   delta[2][3] = rR_star[3] - rL_star[3];
   delta[2][4] =     R_2[0] - rR_star[0];
   mean[2][0] = (real)0.5*( FABS( rL_star[0] ) + FABS(     R_1[0] ) );
   mean[2][1] = (real)0.5*( FABS( rR_star[0] ) + FABS( rL_star[0] ) );
   mean[2][2] = (real)0.5*( FABS( rR_star[2] ) + FABS( rL_star[2] ) );
   mean[2][3] = (real)0.5*( FABS( rR_star[3] ) + FABS( rL_star[3] ) );
   mean[2][4] = (real)0.5*( FABS(     R_2[0] ) + FABS( rR_star[0] ) );


// set limit function
   for (int i=0; i<5; i++)
   {
      if (  FABS( eival[i] ) < MAX_ERROR  )      limit[i] = (real)1.0;
      else
      {
         if ( eival[i] > (real)0.0 )
         {
            if (  mean[0][i] == (real)0.0  ||  mean[1][i] == (real)0.0  )     limit[i] = (real)1.0;
            else
            {
               if (  ( delta[0][i]*delta[1][i] ) / ( mean[0][i]*mean[1][i] ) < MAX_ERROR*MAX_ERROR  )
                  limit[i] = (real)1.0;
               else
               {
                  real r = delta[0][i] / delta[1][i];
                  limit[i] = set_limit( r, eival[i] * ratio, WAF_Limiter );
               }
            }
         }

         else
         {
            if ( mean[2][i] == (real)0.0  ||  mean[1][i] == (real)0.0  )      limit[i] = (real)1.0;
            else
            {
               if (  ( delta[2][i]*delta[1][i] ) / ( mean[2][i]*mean[1][i] ) < MAX_ERROR*MAX_ERROR  )
                  limit[i] = (real)1.0;
               else
               {
                  real r = delta[2][i] / delta[1][i];
                  limit[i] = set_limit( r, eival[i] * ratio, WAF_Limiter );
               }
            }
         }
      }
   } // for (int i=0; i<5; i++)


// solve the sign of waves
   for (int i=0; i<5; i++)
   {
      if (  FABS( eival[i] ) < MAX_ERROR  )     theta[i] =  (real)0.0;
      else if ( eival[i] > (real)0.0 )          theta[i] =  (real)1.0;
      else                                      theta[i] = -(real)1.0;
   }

// solve the intercell flux
#  ifdef WAF_DISSIPATE
   Dis_Stru  ( flux, L_1, R_1, cL_star, cR_star, limit, theta, Gamma );
#  else
   Undis_Stru( flux, L_1, R_1, cL_star, cR_star, limit, theta, Gamma );
#  endif

} // FUNCTION : Solve_Flux



//-----------------------------------------------------------------------------------------------------
// Function    :  set_limit
// Description :  set the limit function
//
// parameter   :  r            : flow variable
//                c            : Courant number
//                WAF_Limiter  : Selection of te limit function
//                                  0 : superbee
//                                  1 : van-Leer
//                                  2 : van-Albada
//                                  3 : minbee
//-------------------------------------------------------------------------------------------------------
__device__ real set_limit( const real r, const real c, const WAF_Limiter_t WAF_Limiter )
{

   real limit;

// choose the limit function
   switch ( WAF_Limiter )
   {
      case WAF_SUPERBEE :
      {
         if ( r > (real)0.0  &&  r <= (real)0.5 )     limit = (real)1.0 - (real)2.0*r*( (real)1.0 - FABS(c) );
         else if ( r <= (real)1.0 )                   limit = FABS(c);
         else if ( r <= (real)2.0 )                   limit = (real)1.0 - r*( (real)1.0 - FABS(c) );
         else                                         limit = (real)2.0*FABS(c) - (real)1.0;
         break;
      }

      case WAF_VANLEER :
      {
         limit = (real)1.0 - (real)2.0*r*( (real)1.0 - FABS(c) ) / ( (real)1.0 + r );
         break;
      }

      case WAF_ALBADA :
      {
         limit = (real)1.0 - r*( (real)1.0 + r )*( (real)1.0 - FABS(c) ) / ( (real)1.0 + r*r );
         break;
      }

      case WAF_MINBEE :
      {
         if ( r > (real)0.0  &&  r <= (real)1.0 )  limit = (real)1.0 - r*( (real)1.0 - FABS(c) );
         else                                      limit = FABS(c);
         break;
      }

      default:
         break;
   }

   return limit;

} // FUNCTION : set_limit



#ifdef WAF_DISSIPATE
//------------------------------------------------------------------------------------------------------
// Function    :  Dis_Stru
// Description :  Set the intercell flux by dissipative wave structure
//
// Parameter   :  flux     : Intercel flux
//                L        : Primitive variables in the left region
//                R        : Primitive variables in the right region
//                L_star   : Primitive variables in the left star region
//                R_star   : Primitive variables in the right star region
//                limit    : Limit functions
//                theta    : Sign of wave speed
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void Dis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5],
                          const real R_star[5], const real limit[5], const real theta[5], const real Gamma )
{

   real iflux[6][5];
   real lim[5];

   for (int i=0; i<5; i++)    lim[i] = limit[i];


// flux function evaluated at the given stat
   set_flux( iflux[0], L,      Gamma );
   set_flux( iflux[1], L_star, Gamma );
   set_flux( iflux[4], R_star, Gamma );
   set_flux( iflux[5], R,      Gamma );


// determine the ghost stats
   real stat[2][5];

   if ( limit[1] <= limit[2] )
   {
      if ( limit[3] <= limit[1] )
      {
         stat[0][0] = L_star[0];
         stat[0][1] = L_star[1];
         stat[0][2] = L_star[2];
         stat[0][3] = R_star[3];
         stat[0][4] = L_star[4];

         stat[1][0] = R_star[0];
         stat[1][1] = L_star[1];
         stat[1][2] = L_star[2];
         stat[1][3] = R_star[3];
         stat[1][4] = L_star[4];
      }

      else if ( limit[3] <= limit[2] )
      {
         stat[0][0] = R_star[0];
         stat[0][1] = L_star[1];
         stat[0][2] = L_star[2];
         stat[0][3] = L_star[3];
         stat[0][4] = L_star[4];

         stat[1][0] = R_star[0];
         stat[1][1] = L_star[1];
         stat[1][2] = L_star[2];
         stat[1][3] = R_star[3];
         stat[1][4] = L_star[4];
      }

      else
      {
         stat[0][0] = R_star[0];
         stat[0][1] = L_star[1];
         stat[0][2] = L_star[2];
         stat[0][3] = L_star[3];
         stat[0][4] = L_star[4];

         stat[1][0] = R_star[0];
         stat[1][1] = L_star[1];
         stat[1][2] = R_star[2];
         stat[1][3] = L_star[3];
         stat[1][4] = L_star[4];
      }
   } // if ( limit[1] <= limit[2] )

   else // limit[1] > limit[2]
   {
      if ( limit[3] <= limit[2] )
      {
         stat[0][0] = L_star[0];
         stat[0][1] = L_star[1];
         stat[0][2] = L_star[2];
         stat[0][3] = R_star[3];
         stat[0][4] = L_star[4];

         stat[1][0] = L_star[0];
         stat[1][1] = L_star[1];
         stat[1][2] = R_star[2];
         stat[1][3] = R_star[3];
         stat[1][4] = L_star[4];
      }

      else if ( limit[3] <= limit[1] )
      {
         stat[0][0] = L_star[0];
         stat[0][1] = L_star[1];
         stat[0][2] = R_star[2];
         stat[0][3] = L_star[3];
         stat[0][4] = L_star[4];

         stat[1][0] = L_star[0];
         stat[1][1] = L_star[1];
         stat[1][2] = R_star[2];
         stat[1][3] = R_star[3];
         stat[1][4] = L_star[4];
      }

      else
      {
         stat[0][0] = L_star[0];
         stat[0][1] = L_star[1];
         stat[0][2] = R_star[2];
         stat[0][3] = L_star[3];
         stat[0][4] = L_star[4];

         stat[1][0] = R_star[0];
         stat[1][1] = L_star[1];
         stat[1][2] = R_star[2];
         stat[1][3] = L_star[3];
         stat[1][4] = L_star[4];
      }
   } // if ( limit[1] <= limit[2] ) ... else ...


// set flux in ghost region
   set_flux( iflux[2], stat[0], Gamma );
   set_flux( iflux[3], stat[1], Gamma );


// reoder the limit
   for (int i=1; i<3; i++)
   {
      if ( lim[i] > lim[i+1] )
      {
         real tmp = lim[i+1];
         lim[i+1] = lim[i  ];
         lim[i  ] = tmp;
      }
   }

   if ( lim[1] > lim[2] )
   {
      real tmp = lim[2];
      lim[2]   = lim[1];
      lim[1]   = tmp;
   }


// set the intercell flux
   for (int i=0; i<5; i++)
    {
       flux[i] =   (real)0.5*( iflux[0][i] + iflux[5][i] )
                 - (real)0.5*(  theta[0]*lim[0]*( iflux[1][i] - iflux[0][i] ) +
                                theta[1]*lim[1]*( iflux[2][i] - iflux[1][i] ) +
                                theta[2]*lim[2]*( iflux[3][i] - iflux[2][i] ) +
                                theta[3]*lim[3]*( iflux[4][i] - iflux[3][i] ) +
                                theta[4]*lim[4]*( iflux[5][i] - iflux[4][i] )   );
    }

} // FUNCTION : Dis_Stru
#endif // #ifdef WAF_DISSIPATE



#ifndef WAF_DISSIPATE
//------------------------------------------------------------------------------------------------------
// Function    :  Undis_Stru
// Description :  Set the intercell flux by non-dissipative wave structure
//
// Parameter   :  flux     : Intercel flux
//                L        : Primitive variables in the left region
//                R        : Primitive variables in the right region
//                L_star   : Primitive variables in the left star region
//                R_star   : Primitive variables in the right star region
//                limit    : Limit functions
//                theta    : Sign of wave speed
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void Undis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5],
                            const real R_star[5], const real limit[5], const real theta[5], const real Gamma )
{

// flux function evaluated at the given stat
   real iflux[4][5];

   set_flux( iflux[0], L,      Gamma );
   set_flux( iflux[1], L_star, Gamma );
   set_flux( iflux[2], R_star, Gamma );
   set_flux( iflux[3], R,      Gamma );


// set the intercell flux
   flux[0] =   (real)0.5*( iflux[0][0] + iflux[3][0] )
             - (real)0.5*(  theta[0]*limit[0]*( iflux[1][0] - iflux[0][0] ) +
                            theta[1]*limit[1]*( iflux[2][0] - iflux[1][0] ) +
                            theta[4]*limit[4]*( iflux[3][0] - iflux[2][0] )   );

   flux[1] =   (real)0.5*( iflux[0][1] + iflux[3][1] )
             - (real)0.5*(  theta[0]*limit[0]*( iflux[1][1] - iflux[0][1] ) +
                            theta[1]*limit[1]*( iflux[2][1] - iflux[1][1] ) +
                            theta[4]*limit[4]*( iflux[3][1] - iflux[2][1] )   );

   flux[4] =   (real)0.5*( iflux[0][4] + iflux[3][4] )
             - (real)0.5*(  theta[0]*limit[0]*( iflux[1][4] - iflux[0][4] ) +
                            theta[1]*limit[1]*( iflux[2][4] - iflux[1][4] ) +
                            theta[4]*limit[4]*( iflux[3][4] - iflux[2][4] )   );

   flux[2] =   (real)0.5*( iflux[0][2] + iflux[3][2] )
             - (real)0.5*(  theta[0]*limit[0]*( iflux[1][2] - iflux[0][2] ) +
                            theta[2]*limit[2]*( iflux[2][2] - iflux[1][2] ) +
                            theta[4]*limit[4]*( iflux[3][2] - iflux[2][2] )   );

   flux[3] =   (real)0.5*( iflux[0][3] + iflux[3][3] )
             - (real)0.5*(  theta[0]*limit[0]*( iflux[1][3] - iflux[0][3] ) +
                            theta[3]*limit[3]*( iflux[2][3] - iflux[1][3] ) +
                            theta[4]*limit[4]*( iflux[3][3] - iflux[2][3] )   );

} // FUNCTION : Undis_Stru
#endif // #ifndef WAF_DISSIPATE



//-------------------------------------------------------------------------------------------------------
// Function    :  set_flux
// Description :  Set the flux function evaluated at the given stat
//
// Parameter   :  flux  : Flux function
//                val   : Primitive variables
//                Gamma : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void set_flux( real flux[5], const real val[5], const real Gamma )
{

  const real Gamma_m1 = Gamma -(real)1.0;

// set flux
   flux[0] = val[0]*val[1];
   flux[1] = val[0]*val[1]*val[1] + val[4];
   flux[2] = val[0]*val[1]*val[2];
   flux[3] = val[0]*val[1]*val[3];
   flux[4] = val[1]*(  (real)0.5*val[0]*( val[1]*val[1] + val[2]*val[2] + val[3]*val[3] )
                      + val[4]/Gamma_m1 + val[4]  );

} // FUNCTION : set_flux



#if ( RSOLVER == ROE )
//-------------------------------------------------------------------------------------------------------
// Function    :  Solve_StarRoe
// Description :  Solve the star region and speed of waves by Roe's method
//
// Parameter   :  eival   : Speed of waves
//                L_star  : Primitive variables in the left star region
//                R_star  : Primitive variables in the right star region
//                L       : Primitive variables in the left rrgion
//                R       : Primitive variables in the right rrgion
//                Gamma   : Ratio of specific heats
//                MinPres : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
__device__ void Solve_StarRoe( real eival[5], real L_star[5], real R_star[5], const real L[5], const real R[5],
                               const real Gamma, const real MinPres )
{

   const real Gamma_m1 = Gamma - (real)1.0; // for evaluating pressure and sound speed

   real u_bar, v_bar, w_bar, h_bar, a_bar, a_bar_inv;   // Roe's average of vx, vy, vz, enthapy, sound speed, and
                                                        // one over a_bar
   real coef[5];                                        // Roe's coefficients
   real TempPres, TempRho, _TempRho;

// solve Roe's average
   {
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(L[0]) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L[0], __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R[0]) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real n_L_sq = SQRT( L[0] ); // rooting number of left density
      real n_R_sq = SQRT( R[0] ); // rooting number of right density
      real h_L = (real)0.5*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] ) + Gamma/Gamma_m1*L[4]/L[0]; // left enthapy
      real h_R = (real)0.5*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] ) + Gamma/Gamma_m1*R[4]/R[0]; // right enthapy
      real n_bar_inv = (real)1.0 / ( n_L_sq + n_R_sq ); // one over (n_L_sq plus n_L_sq)

      u_bar = ( n_L_sq*L[1] + n_R_sq*R[1] )*n_bar_inv;
      v_bar = ( n_L_sq*L[2] + n_R_sq*R[2] )*n_bar_inv;
      w_bar = ( n_L_sq*L[3] + n_R_sq*R[3] )*n_bar_inv;
      h_bar = ( n_L_sq*h_L  + n_R_sq*h_R  )*n_bar_inv;

      real GammaP_Rho = Gamma_m1*(  h_bar - (real)0.5*( u_bar*u_bar + v_bar*v_bar + w_bar*w_bar )  );

      TempRho    = (real)0.5*( L[0] + R[0] );
      _TempRho   = (real)1.0/TempRho;
      TempPres   = GammaP_Rho*TempRho/Gamma;
      TempPres   = CUFLU_CheckMinPres( TempPres, MinPres );
      GammaP_Rho = Gamma*TempPres*_TempRho;

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(GammaP_Rho) )
         printf( "ERROR : negative GammaP_Rho (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 GammaP_Rho, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      a_bar     = SQRT( GammaP_Rho );
      a_bar_inv = (real)1.0 / a_bar;
   }


// solve Roe's coefficients
   {
//    the difference of conservative variables
      real du_1 = R[0] - L[0];
      real du_2 = R[0]*R[1] - L[0]*L[1];
      real du_3 = R[0]*R[2] - L[0]*L[2];
      real du_4 = R[0]*R[3] - L[0]*L[3];
      real du_5 = + (real)0.5*R[0]*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] ) + R[4]/Gamma_m1
                  - (real)0.5*L[0]*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] ) - L[4]/Gamma_m1;

      coef[2] = du_3 - v_bar*du_1;
      coef[3] = du_4 - w_bar*du_1;
      coef[1] = Gamma_m1*a_bar_inv*a_bar_inv*(  du_1*( h_bar - u_bar*u_bar ) + u_bar*du_2 - du_5
                                               + coef[2]*v_bar + coef[3]*w_bar  );
      coef[0] = (real)0.5*a_bar_inv*( du_1*( u_bar + a_bar ) - du_2 - a_bar*coef[1] );
      coef[4] = du_1 - ( coef[0] + coef[1] );
   }


// solve the star region
   {
      L_star[0] = L[0] + coef[0];
      R_star[0] = R[0] - coef[4];
      L_star[1] = (real)0.5*(   (  L[0]*L[1] + coef[0]*( u_bar - a_bar )  ) / L_star[0]
                              + (  R[0]*R[1] - coef[4]*( u_bar + a_bar )  ) / R_star[0]   );
      R_star[1] = L_star[1];
      L_star[2] = L[2];
      R_star[2] = R[2];
      L_star[3] = L[3];
      R_star[3] = R[3];
      real E_L = (real)0.5*L[0]*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] );
      real E_R = (real)0.5*R[0]*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] );
      real e_L_star = (real)0.5*L_star[0]*( L_star[1]*L_star[1] + L_star[2]*L_star[2] + L_star[3]*L_star[3] );
      real e_R_star = (real)0.5*R_star[0]*( R_star[1]*R_star[1] + R_star[2]*R_star[2] + R_star[3]*R_star[3] );
      L_star[4] = (real)0.5*Gamma_m1*(  E_L - e_L_star + L[4]/Gamma_m1 + coef[0]*( h_bar - u_bar*a_bar )
                                      + E_R - e_R_star + R[4]/Gamma_m1 - coef[4]*( h_bar + u_bar*a_bar )  );
      L_star[4] = CUFLU_CheckMinPres( L_star[4], MinPres );
      R_star[4] = L_star[4];
   }

// solve the speed of waves
   {
      real eigen[2];
      eival[1] = L_star[1];
      eival[2] = L_star[1];
      eival[3] = L_star[1];

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(L[4]) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L[4],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(L[0]) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L[0],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(L_star[4]) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L_star[4], __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(L_star[0]) )
         printf( "ERROR : negative density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L_star[0], __FILE__, __LINE__, __FUNCTION__ );

      if ( CUFLU_CheckNegative(R[4]) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R[4],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R[0]) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R[0],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R_star[4]) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R_star[4], __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R_star[0]) )
         printf( "ERROR : negative density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R_star[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif

      eigen[0] = L     [1] - SQRT( Gamma*L     [4]/L     [0] );
      eigen[1] = L_star[1] - SQRT( Gamma*L_star[4]/L_star[0] );
      if ( eigen[0] <= eigen[1] )   eival[0] = eigen[0];
      else                          eival[0] = eigen[1];

      eigen[0] = R     [1] + SQRT( Gamma*R     [4]/R     [0] );
      eigen[1] = R_star[1] + SQRT( Gamma*R_star[4]/R_star[0] );
      if ( eigen[0] <= eigen[1] )   eival[4] = eigen[1];
      else                          eival[4] = eigen[0];
   }

} // FUNCTION : Solve_StarRoe
#endif // #if ( RSOLVER == ROE )



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == WAF )
