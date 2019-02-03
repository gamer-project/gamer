#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == RTVD )


// check before compiling anything else
#if ( NCOMP_PASSIVE != 0 )
#  error : RTVD scheme does NOT support passive scalars !!
#endif


#include "CUFLU_Shared_FluUtility.cu"

#define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
#define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )

static __device__ void CUFLU_Advance( real g_Fluid_In [][5][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][5][ CUBE(PS2) ],
                                      real g_Flux[][9][5][ SQR(PS2) ],
                                      const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                      const int j_gap, const int k_gap, real s_cu[][5][FLU_NXT],
                                      real s_cw[][5][FLU_NXT], real s_flux[][5][FLU_NXT],
                                      real s_RLflux[][5][FLU_NXT], const bool FinalOut, const int XYZ,
                                      const real MinDens, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_RTVD
// Description :  GPU fluid solver based on the relaxing TVD (RTVD) scheme
//
// Note        :  a. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                b. The three-dimensional evolution is achieved by using the dimensional-split method
//
// Parameter   :  g_Fluid_In   : Global memory array to store the input fluid variables
//                g_Fluid_Out  : Global memory array to store the output fluid variables
//                g_Flux       : Global memory array to store the output fluxes
//                g_Corner     : Global memory array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                g_Pot_USG    : Global memory array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in RTVD)
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Gamma        : Ratio of specific heats
//                StoreFlux    : true --> store the coarse-fine fluxes
//                XYZ          : true  : x->y->z ( forward sweep)
//                               false : z->y->x (backward sweep)
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_FluidSolver_RTVD(
   real g_Fluid_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   real g_Fluid_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   real g_Flux     [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double g_Corner[][3],
   const real g_Pot_USG[][ CUBE(USG_NXT_F) ],
   const real dt, const real _dh, const real Gamma, const bool StoreFlux,
   const bool XYZ, const real MinDens, const real MinPres )
{

   __shared__ real s_cu    [FLU_BLOCK_SIZE_Y][5][FLU_NXT];
   __shared__ real s_cw    [FLU_BLOCK_SIZE_Y][5][FLU_NXT];
   __shared__ real s_flux  [FLU_BLOCK_SIZE_Y][5][FLU_NXT];
   __shared__ real s_RLflux[FLU_BLOCK_SIZE_Y][5][FLU_NXT];

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux,              0,              0,
                     s_cu, s_cw, s_flux, s_RLflux, false, 0, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux, FLU_GHOST_SIZE,              0,
                     s_cu, s_cw, s_flux, s_RLflux, false, 3, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE,
                     s_cu, s_cw, s_flux, s_RLflux,  true, 6, MinDens, MinPres );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux,              0,              0,
                     s_cu, s_cw, s_flux, s_RLflux, false, 6, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux,              0, FLU_GHOST_SIZE,
                     s_cu, s_cw, s_flux, s_RLflux, false, 3, MinDens, MinPres );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Gamma, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE,
                     s_cu, s_cw, s_flux, s_RLflux,  true, 0, MinDens, MinPres );
   }

} // FUNCTION : CUFLU_FluidSolver_RTVD



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  GPU device function, which performs a one-dimensional sweep based on the TVD scheme
//
// Note        :  a. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                b. The direction of the one dimensional sweep is determined by the input parameter "XYZ"
//
// Parameter   :  g_Fluid_In   : Global memory array to store the input fluid variables
//                g_Fluid_Out  : Global memory array to store the output fluid variables
//                g_Flux       : Global memory array to store the output fluxes
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Gamma        : Ratio of specific heats
//                StoreFlux    : true --> store the coarse-fine fluxes
//                j_gap        : Number of useless grids in each side in the j direction (j may not be equal to y)
//                k_gap        : Number of useless grids in each side in the k direction (k mya not be equal to z)
//                s_cu         : Shared memory array storing the normal flux
//                s_cw         : Shared memory array storing the auxiliary flux
//                s_flux       : Shared memory array storing the final flux used to update the fluid variables
//                s_RLflux     : Shared memory array storing the left/right-moving flux
//                XYZ          : 0 : Update the solution in the x direction
//                               3 : Update the solution in the y direction
//                               6 : Update the solution in the z direction
//                               --> This parameter is also used to determine the place to store the output fluxes
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_Advance( real g_Fluid_In [][5][ CUBE(FLU_NXT) ],
                               real g_Fluid_Out[][5][ CUBE(PS2) ],
                               real g_Flux[][9][5][ SQR(PS2) ],
                               const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                               const int j_gap, const int k_gap, real s_cu[][5][FLU_NXT],
                               real s_cw[][5][FLU_NXT], real s_flux[][5][FLU_NXT], real s_RLflux[][5][FLU_NXT],
                               const bool FinalOut, const int XYZ, const real MinDens, const real MinPres )
{

   const unsigned int bx       = blockIdx.x;
   const unsigned int tx       = threadIdx.x;
   const unsigned int ty       = threadIdx.y;
   const unsigned int dj       = blockDim.y;
   const unsigned int size_j   = FLU_NXT - (j_gap<<1);
   const unsigned int size_k   = FLU_NXT - (k_gap<<1);
   const unsigned int NColumn  = __umul24( size_j, size_k );
   const unsigned int i        = tx;                        // (i,j) the element in shared memory under evaluation
   const unsigned int ip       = i+1;
   const unsigned int im       = i-1;
         unsigned int j        = j_gap + ty%size_j;
         unsigned int k        = k_gap + ty/size_j;
         unsigned int Column0  = 0;                         // the total number of columns that have been updated
   const unsigned int j_end    = FLU_NXT - j_gap;
   const unsigned int k_end    = FLU_NXT - k_gap;
   const real         Gamma_m1 = Gamma - (real)1.0;
   const real        _Gamma_m1 = (real)1.0 / Gamma_m1;
   const real         dt_half  = (real)0.5*dt;
         bool         RuleOut  = false;

   real   _rho, vx, p, c, Temp, Fluid[5], Fluid_half[5];
   int    ID1, ID2, ID3, Comp[5], delta_k;


// set the order of component for update in different directions
   switch ( XYZ )
   {
      case 0:  Comp[0] = 0;   Comp[1] = 1;   Comp[2] = 2;   Comp[3] = 3;   Comp[4] = 4;   break;
      case 3:  Comp[0] = 0;   Comp[1] = 2;   Comp[2] = 1;   Comp[3] = 3;   Comp[4] = 4;   break;
      case 6:  Comp[0] = 0;   Comp[1] = 3;   Comp[2] = 2;   Comp[3] = 1;   Comp[4] = 4;   break;
   }


// start the TVD scheme
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


//    a. Evaluate the half-step values of fluid variables
//-----------------------------------------------------------------------------

//    (a1). set variables defined in the center of cell
      _rho = (real)1.0 / Fluid[0];
      vx   = _rho * Fluid[1];
      p    = Gamma_m1*(  Fluid[4] - (real)0.5*_rho*( Fluid[1]*Fluid[1]+Fluid[2]*Fluid[2]+Fluid[3]*Fluid[3] )  );
      p    = Hydro_CheckMinPres( p, MinPres );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(p) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 p, __FILE__, __LINE__, __FUNCTION__ );

      if ( Hydro_CheckNegative(Fluid[0]) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Fluid[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif
      c    = FABS( vx ) + SQRT( Gamma*p*_rho );

      s_cw[ty][0][i] = Fluid[1];
      s_cw[ty][1][i] = Fluid[1]*vx + p;
      s_cw[ty][2][i] = Fluid[2]*vx;
      s_cw[ty][3][i] = Fluid[3]*vx;
      s_cw[ty][4][i] = ( Fluid[4]+p )*vx;

      s_cu[ty][0][i] = c*Fluid[0];
      s_cu[ty][1][i] = c*Fluid[1];
      s_cu[ty][2][i] = c*Fluid[2];
      s_cu[ty][3][i] = c*Fluid[3];
      s_cu[ty][4][i] = c*Fluid[4];

      __syncthreads();


//    (a2). set flux defined in the right-hand surface of cell by the upwind scheme
      if ( i < FLU_NXT-1 )
      {
         for (int v=0; v<5; v++)
            s_flux[ty][v][i] = (real)0.5*(  ( s_cu[ty][v][i ]+s_cw[ty][v][i ] ) -
                                            ( s_cu[ty][v][ip]-s_cw[ty][v][ip] )  );
      }

      __syncthreads();


//    (a3). evaluate the intermidiate values (u_half)
//    if ( i > 0 )
      if ( i > 0  &&  i < FLU_NXT-1 )
      {
         for (int v=0; v<5; v++)    Fluid_half[v] = Fluid[v] - _dh*dt_half*( s_flux[ty][v][i] - s_flux[ty][v][im] ) ;

//       enforce positive density and pressure
         Fluid_half[0] = FMAX( Fluid_half[0], MinDens );
         Fluid_half[4] = Hydro_CheckMinPresInEngy( Fluid_half[0], Fluid_half[1], Fluid_half[2], Fluid_half[3], Fluid_half[4],
                                                   Gamma_m1, _Gamma_m1, MinPres );
      }



//    Evaluate the full-step values of fluid variables
//-----------------------------------------------------------------------------

//    (b1). reset variables defined in the center of cell at the intermidate state
      if ( i > 0  &&  i < FLU_NXT-1 )
      {
         _rho = (real)1.0 / Fluid_half[0];
         vx   = _rho * Fluid_half[1];
         p    = Gamma_m1*(  Fluid_half[4]
                           - (real)0.5*_rho*( Fluid_half[1]*Fluid_half[1] + Fluid_half[2]*Fluid_half[2] +
                                              Fluid_half[3]*Fluid_half[3] )  );
         p    = Hydro_CheckMinPres( p, MinPres );

#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( Hydro_CheckNegative(p) )
            printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    p, __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(Fluid_half[0]) )
            printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid_half[0], __FILE__, __LINE__, __FUNCTION__ );
#        endif

         c    = FABS( vx ) + SQRT( Gamma*p*_rho );

         s_cw[ty][0][i] = Fluid_half[1];
         s_cw[ty][1][i] = Fluid_half[1]*vx + p;
         s_cw[ty][2][i] = Fluid_half[2]*vx;
         s_cw[ty][3][i] = Fluid_half[3]*vx;
         s_cw[ty][4][i] = ( Fluid_half[4]+p )*vx;

         s_cu[ty][0][i] = c*Fluid_half[0];
         s_cu[ty][1][i] = c*Fluid_half[1];
         s_cu[ty][2][i] = c*Fluid_half[2];
         s_cu[ty][3][i] = c*Fluid_half[3];
         s_cu[ty][4][i] = c*Fluid_half[4];
      }


//    (b2). set the right-moving flux defined in the right-hand surface by the TVD scheme
      if ( i > 0  &&  i < FLU_NXT-2 )
      {
         for (int v=0; v<5; v++)    s_RLflux[ty][v][i] = (real)0.5*( s_cu[ty][v][i] + s_cw[ty][v][i] );
      }

      __syncthreads();


      if ( i > 1  &&  i < FLU_NXT-3 )
      {
         for (int v=0; v<5; v++)
         {
            s_flux[ty][v][i] = s_RLflux[ty][v][i];

            Temp =   ( s_RLflux[ty][v][ip]-s_RLflux[ty][v][i ] )
                   * ( s_RLflux[ty][v][i ]-s_RLflux[ty][v][im] );

            if ( Temp > (real)0.0 )
               s_flux[ty][v][i] += Temp / ( s_RLflux[ty][v][ip]-s_RLflux[ty][v][im] );
         }
      }

      __syncthreads();


//    (b3). set the left-moving flux defined in the left-hand surface by the TVD scheme, get the total flux
//    if ( i < FLU_NXT-2 )
      if ( i > 0  &&  i < FLU_NXT-2 )
      {
         for (int v=0; v<5; v++)    s_RLflux[ty][v][i] = (real)0.5*( s_cu[ty][v][ip] - s_cw[ty][v][ip] );
      }

      __syncthreads();


      if ( i > 1  &&  i < FLU_NXT-3 )
      {
         for (int v=0; v<5; v++)
         {
            s_flux[ty][v][i] -= s_RLflux[ty][v][i];

            Temp =   ( s_RLflux[ty][v][im]-s_RLflux[ty][v][i ] )
                   * ( s_RLflux[ty][v][i ]-s_RLflux[ty][v][ip] );

            if ( Temp > (real)0.0 )
               s_flux[ty][v][i] -= Temp / ( s_RLflux[ty][v][im]-s_RLflux[ty][v][ip] );
         }
      }

      __syncthreads();


//    (b4). advance fluid by one full time-step
//    if ( i > 2 )
//    if ( i > 2  &&  i < FLU_NXT-3 )
      if ( i > 2  &&  i < FLU_NXT-3  &&  RuleOut == false )
      {
         for (int v=0; v<5; v++)    Fluid[v] -= _dh*dt*( s_flux[ty][v][i] - s_flux[ty][v][im] ) ;

//       enforce positive density and pressure
         Fluid[0] = FMAX( Fluid[0], MinDens );
         Fluid[4] = Hydro_CheckMinPresInEngy( Fluid[0], Fluid[1], Fluid[2], Fluid[3], Fluid[4],
                                              Gamma_m1, _Gamma_m1, MinPres );


//       check negative density and energy
#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( Hydro_CheckNegative(Fluid[0]) )
            printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid[0], __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(Fluid[4]) )
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


//    (b5). save the flux required by the flux-correction operation
      if ( StoreFlux )
      if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
      if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
      if ( i == 0 )
      {
         ID3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

         for (int v=0; v<5; v++)
         {
            g_Flux[bx][XYZ+0][v][ID3] = s_flux[ty][ Comp[v] ][          2];
            g_Flux[bx][XYZ+1][v][ID3] = s_flux[ty][ Comp[v] ][FLU_NXT/2-1];
            g_Flux[bx][XYZ+2][v][ID3] = s_flux[ty][ Comp[v] ][FLU_NXT - 4];
         }
      }


//    reset the target array indices
      j += dj;

      if ( j >= j_end )
      {
         delta_k   = ( j - j_end )/size_j + 1;
         k        += delta_k;
         j        -= __umul24( size_j, delta_k );
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

} // FUNCTION : CUFLU_Advance



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == RTVD )
