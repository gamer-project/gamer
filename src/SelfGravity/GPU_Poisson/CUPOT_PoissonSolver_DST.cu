#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == DST )


#define POT_NXT_F    ( PATCH_SIZE+2*POT_GHOST_SIZE           )
#define POT_PAD      ( WARP_SIZE/2 - (POT_NXT_F*2%WARP_SIZE) )
#define POT_NTHREAD  ( RHO_NXT*RHO_NXT/2    )
#define POT_USELESS  ( POT_GHOST_SIZE%2                      )

#include "SineTransformTables.h"
#include "CUDA_ConstMemory.h"

__shared__  real shared_mem[RHO_NXT*RHO_NXT];

__device__ uint Rhoid_3Dto1D(uint x, uint y, uint z, uint N, uint XYZ){

  if (XYZ==0) return __umul24(z , N*N) + __umul24(y , N) + x;
  if (XYZ==1) return __umul24(z , N*N) + __umul24(x , N) + y;
  if (XYZ==2) return __umul24(x , N*N) + __umul24(z , N) + y;

  return 0;
}

__device__ void DST_Scheme(const uint ID,
			   real *Rho_Array,
			   uint XYZ){

  uint t,Rhoid_x,Rhoid_y,Rhoid_z,Rhoid;
  uint N = RHO_NXT ;
  uint stride = blockDim.x * blockDim.y;
  real c;


  // pass part of Rho to shared memory
  for (int step=0; step<RHO_NXT; step++){

    t = ID ;
    stride = blockDim.x * blockDim.y;
    do 
    {

     Rhoid_z = t / (N*N);
     Rhoid_y = t / N % N;
     Rhoid_x = t % N ;

     Rhoid =  Rhoid_3Dto1D( Rhoid_x, 
			    Rhoid_y, 
			    Rhoid_z + step,
			    N,XYZ);
       
     shared_mem[t ]  =   Rho_Array[Rhoid];
     t += stride;
     
    } while (t < N*N*N/RHO_NXT);
    __syncthreads();
     
  // execute FFT
    t=ID;
    do 
    {

     Rhoid_z = t / (N*N);
     Rhoid_y = t / N % N;
     Rhoid_x = t % N ;
     Rhoid =  Rhoid_3Dto1D( Rhoid_x, 
			    Rhoid_y, 
			    Rhoid_z + step ,
			    N,XYZ);

     c = 0;
     for (int e=0; e<N; e++)
       c += shared_mem[(Rhoid_z * N + Rhoid_y)*N + e] * M[Rhoid_x * N + e];
     
     Rho_Array[Rhoid] = c * 2.0;
     t += stride;

    } while (t < N*N*N/RHO_NXT);
    __syncthreads();

  } //   for (int step=0; step<RHO_NXT; step++){

    __syncthreads();
}

__device__ void Assign_sFPot(   real TempFPot,int FID,
				int bid,
				int FIDxx, int FIDyy,int FIDzz)

{
  // Note: Naive way to assign boundary condition only at the outermost layer using if-else statement.
  
    int  FIDz = FID /(POT_NXT_F*POT_NXT_F);
    int  FIDy = (FID /POT_NXT_F) % POT_NXT_F ;
    int  FIDx = FID % POT_NXT_F;

    if (FIDx==0         and FIDy >= 1 and FIDy <= RHO_NXT and FIDz >= 1 and FIDz <= RHO_NXT) shared_mem[RHO_NXT*RHO_NXT * 0 + (FIDz-1)*RHO_NXT + FIDy-1 ] = TempFPot; 
    if (FIDy==0         and FIDx >= 1 and FIDx <= RHO_NXT and FIDz >= 1 and FIDz <= RHO_NXT) shared_mem[RHO_NXT*RHO_NXT * 1 + (FIDz-1)*RHO_NXT + FIDx-1] = TempFPot;
    if (FIDz==0         and FIDx >= 1 and FIDx <= RHO_NXT and FIDy >= 1 and FIDy <= RHO_NXT) shared_mem[RHO_NXT*RHO_NXT * 2 + (FIDy-1)*RHO_NXT + FIDx-1] = TempFPot;
      
    if (FIDx==RHO_NXT+1 and FIDy >= 1 and FIDy <= RHO_NXT and FIDz >= 1 and FIDz <= RHO_NXT) shared_mem[RHO_NXT*RHO_NXT * 3 + (FIDz-1)*RHO_NXT + FIDy-1] = TempFPot;
    if (FIDy==RHO_NXT+1 and FIDx >= 1 and FIDx <= RHO_NXT and FIDz >= 1 and FIDz <= RHO_NXT) shared_mem[RHO_NXT*RHO_NXT * 4 + (FIDz-1)*RHO_NXT + FIDx-1] = TempFPot;
    if (FIDz==RHO_NXT+1 and FIDx >= 1 and FIDx <= RHO_NXT and FIDy >= 1 and FIDy <= RHO_NXT) shared_mem[RHO_NXT*RHO_NXT * 5 + (FIDy-1)*RHO_NXT + FIDx-1] = TempFPot;

}
//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_PoissonSolver_SOR
// Description :  GPU Poisson solver using the Discret Sine transform
//
// Method (detail) :  1) Assign boundary condition to the outermost layer of density array. 
//                    2) Apply 1D sine transform Type-I to all three directions
//                    3) mulitply by the the eigen value of the 2nd order finte difference matrix
//                    4) apply 1D sine transform Type-I again to all three directions
//
// Parameter   :  g_Rho_Array     : Global memory array to store the input density
//                g_Pot_Array_In  : Global memory array storing the input "coarse-grid" potential for
//                                  interpolation
//                g_Pot_Array_Out : Global memory array to store the output potential
//                Const           : (Coefficient in front of the RHS in the Poisson eq.) / dh^2
//                IntScheme       : Interpolation scheme for potential
//                                  --> currently supported schemes include
//                                      INT_CQUAD : conservative quadratic interpolation
//                                      INT_QUAD  : quadratic interpolation
//---------------------------------------------------------------------------------------------------
//                           
__global__ void CUPOT_PoissonSolver_DST(       real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ], 
                                               real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ], 
                                               real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ], 
					       const real Const, 
					       const IntScheme_t IntScheme)
{


   const uint bid       = blockIdx.x;
   const uint tid_x     = threadIdx.x;
   const uint tid_y     = threadIdx.y;
   const uint tid_z     = threadIdx.z;
   const uint ID        = __umul24( tid_z, __umul24(blockDim.x,blockDim.y) ) + __umul24( tid_y, blockDim.x ) + tid_x;

   uint t, stride;
   uint Rhoid, Rhoid_x,Rhoid_y,Rhoid_z;




// a2. load the coarse-grid potential into the shared memory
// -----------------------------------------------------------------------------------------------------------
#  ifdef SOR_CPOT_SHARED
   t = ID;
   do {  s_CPot[t] = g_Pot_Array_In[bid][t];    t += POT_NTHREAD; }     while ( t < POT_NXT*POT_NXT*POT_NXT );
   __syncthreads();
#  else
   const real *s_CPot = g_Pot_Array_In[bid];
#  endif

// b. evaluate the "fine-grid" potential by interpolation (as the initial guess and the B.C.)
// -----------------------------------------------------------------------------------------------------------
   const int N_CSlice = POT_NTHREAD / ( (POT_NXT-2)*(POT_NXT-2) );

   if ( ID < N_CSlice*(POT_NXT-2)*(POT_NXT-2) )
   {
      const real Const_8   = 1.0/8.0;
      const real Const_64  = 1.0/64.0;
      const real Const_512 = 1.0/512.0;

      const int Cdx  = 1;
      const int Cdy  = POT_NXT;
      const int Cdz  = POT_NXT*POT_NXT;
      const int CIDx = 1 + ID % ( POT_NXT-2 );
      const int CIDy = 1 + (  ID % ( (POT_NXT-2)*(POT_NXT-2) )  ) / ( POT_NXT-2 );
      const int CIDz = 1 + ID / ( (POT_NXT-2)*(POT_NXT-2) );
      int       CID  = __mul24( CIDz, Cdz ) + __mul24( CIDy, Cdy ) + __mul24( CIDx, Cdx );
      const int Fdx  = 1;
      const int Fdy  = POT_NXT_F;
      const int FIDx = ( (CIDx-1)<<1 ) - POT_USELESS;
      const int FIDy = ( (CIDy-1)<<1 ) - POT_USELESS;
      int       FIDz = ( (CIDz-1)<<1 ) - POT_USELESS;
#     ifdef SOR_USE_PADDING
      const int Fpad = ( FIDy < 3 ) ? 0 : POT_PAD*((FIDy-3)/4 + 1);    // padding logic
      const int Fdz  = POT_NXT_F*POT_NXT_F + POT_PAD*4;                // added padding
#     else
      const int Fpad = 0;
      const int Fdz  = POT_NXT_F*POT_NXT_F;
#     endif
      int       FID  = Fpad + __mul24( FIDz, Fdz ) + __mul24( FIDy, Fdy ) + __mul24( FIDx, Fdx );

      real TempFPot1, TempFPot2, TempFPot3, TempFPot4, TempFPot5, TempFPot6, TempFPot7, TempFPot8;
      real Slope_00, Slope_01, Slope_02, Slope_03, Slope_04, Slope_05, Slope_06, Slope_07;
      real Slope_08, Slope_09, Slope_10, Slope_11, Slope_12;
      int  Idx, Idy, Idz, ii, jj, kk;


      for (int z=CIDz; z<POT_NXT-1; z+=N_CSlice)
      {
         switch ( IntScheme )
         {
            /*
            case INT_CENTRAL :
            {
               Slope_00 = (real)0.125 * ( s_CPot[CID+Cdx] - s_CPot[CID-Cdx] );
               Slope_01 = (real)0.125 * ( s_CPot[CID+Cdy] - s_CPot[CID-Cdy] );
               Slope_02 = (real)0.125 * ( s_CPot[CID+Cdz] - s_CPot[CID-Cdz] );

               TempFPot1 = s_CPot[CID] - Slope_00 - Slope_01 - Slope_02;
               TempFPot2 = s_CPot[CID] + Slope_00 - Slope_01 - Slope_02;
               TempFPot3 = s_CPot[CID] - Slope_00 + Slope_01 - Slope_02;
               TempFPot4 = s_CPot[CID] + Slope_00 + Slope_01 - Slope_02;
               TempFPot5 = s_CPot[CID] - Slope_00 - Slope_01 + Slope_02;
               TempFPot6 = s_CPot[CID] + Slope_00 - Slope_01 + Slope_02;
               TempFPot7 = s_CPot[CID] - Slope_00 + Slope_01 + Slope_02;
               TempFPot8 = s_CPot[CID] + Slope_00 + Slope_01 + Slope_02;
            }
            break; // INT_CENTRAL
            */


            case INT_CQUAD :
            {
               Slope_00 = Const_8   * ( s_CPot[CID+Cdx        ] - s_CPot[CID-Cdx        ] );
               Slope_01 = Const_8   * ( s_CPot[CID    +Cdy    ] - s_CPot[CID    -Cdy    ] );
               Slope_02 = Const_8   * ( s_CPot[CID        +Cdz] - s_CPot[CID        -Cdz] );

               Slope_03 = Const_64  * ( s_CPot[CID+Cdx    -Cdz] - s_CPot[CID-Cdx    -Cdz] );
               Slope_04 = Const_64  * ( s_CPot[CID    +Cdy-Cdz] - s_CPot[CID    -Cdy-Cdz] );
               Slope_05 = Const_64  * ( s_CPot[CID+Cdx-Cdy    ] - s_CPot[CID-Cdx-Cdy    ] );
               Slope_06 = Const_64  * ( s_CPot[CID+Cdx+Cdy    ] - s_CPot[CID-Cdx+Cdy    ] );
               Slope_07 = Const_64  * ( s_CPot[CID+Cdx    +Cdz] - s_CPot[CID-Cdx    +Cdz] );
               Slope_08 = Const_64  * ( s_CPot[CID    +Cdy+Cdz] - s_CPot[CID    -Cdy+Cdz] );

               Slope_09 = Const_512 * ( s_CPot[CID+Cdx-Cdy-Cdz] - s_CPot[CID-Cdx-Cdy-Cdz] );
               Slope_10 = Const_512 * ( s_CPot[CID+Cdx+Cdy-Cdz] - s_CPot[CID-Cdx+Cdy-Cdz] );
               Slope_11 = Const_512 * ( s_CPot[CID+Cdx-Cdy+Cdz] - s_CPot[CID-Cdx-Cdy+Cdz] );
               Slope_12 = Const_512 * ( s_CPot[CID+Cdx+Cdy+Cdz] - s_CPot[CID-Cdx+Cdy+Cdz] );

               TempFPot1 = - Slope_00 - Slope_01 - Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                           + Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot2 = + Slope_00 - Slope_01 - Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                           - Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];

               TempFPot3 = - Slope_00 + Slope_01 - Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                           + Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];

               TempFPot4 = + Slope_00 + Slope_01 - Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                           - Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot5 = - Slope_00 - Slope_01 + Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                           - Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];

               TempFPot6 = + Slope_00 - Slope_01 + Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                           + Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot7 = - Slope_00 + Slope_01 + Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                           - Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot8 = + Slope_00 + Slope_01 + Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                           + Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];
            }
            break; // INT_CQUAD

            case INT_QUAD :
            {
               TempFPot1 = TempFPot2 = TempFPot3 = TempFPot4 = (real)0.0;
               TempFPot5 = TempFPot6 = TempFPot7 = TempFPot8 = (real)0.0;

               for (int dk=-1; dk<=1; dk++)  {  Idz = dk+1;    kk = __mul24( dk, Cdz );
               for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
               for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

                  TempFPot1 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mm[Idy] * c_Mm[Idx];
                  TempFPot2 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mm[Idy] * c_Mp[Idx];
                  TempFPot3 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mp[Idy] * c_Mm[Idx];
                  TempFPot4 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mp[Idy] * c_Mp[Idx];
                  TempFPot5 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mm[Idy] * c_Mm[Idx];
                  TempFPot6 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mm[Idy] * c_Mp[Idx];
                  TempFPot7 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mp[Idy] * c_Mm[Idx];
                  TempFPot8 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mp[Idy] * c_Mp[Idx];

               }}}
            }
            break; // INT_QUAD

         } // switch ( IntScheme )

//       save data to the shared-memory array.
//       Currently this part is highly diverged. However, since the interpolation takes much less time than the
//       SOR iteration does, we have not yet tried to optimize this part

         if ( FIDz >= 0 )
         {
            if ( FIDx >= 0            &&  FIDy >= 0           )   Assign_sFPot(TempFPot1, FID,bid,FIDx,FIDy,FIDz);;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy >= 0           )   Assign_sFPot(TempFPot2, FID+Fdx,bid,FIDx,FIDy,FIDz);
            if ( FIDx >= 0            &&  FIDy <= POT_NXT_F-2 )   Assign_sFPot(TempFPot3, FID    +Fdy,bid,FIDx,FIDy,FIDz);
            if ( FIDx <= POT_NXT_F-2  &&  FIDy <= POT_NXT_F-2 )   Assign_sFPot(TempFPot4, FID+Fdx+Fdy  ,bid,FIDx,FIDy,FIDz);
         }
	 
         if ( FIDz <= POT_NXT_F-2 )
         {
            if ( FIDx >= 0            &&  FIDy >= 0           )   Assign_sFPot(TempFPot5, FID        +Fdz,bid,FIDx,FIDy,FIDz);
            if ( FIDx <= POT_NXT_F-2  &&  FIDy >= 0           )   Assign_sFPot(TempFPot6, FID+Fdx    +Fdz,bid,FIDx,FIDy,FIDz);
            if ( FIDx >= 0            &&  FIDy <= POT_NXT_F-2 )   Assign_sFPot(TempFPot7, FID    +Fdy+Fdz,bid,FIDx,FIDy,FIDz);
            if ( FIDx <= POT_NXT_F-2  &&  FIDy <= POT_NXT_F-2 )   Assign_sFPot(TempFPot8, FID+Fdx+Fdy+Fdz,bid,FIDx,FIDy,FIDz);
         }


	 
         CID  += __mul24(   N_CSlice, Cdz );
         FID  += __mul24( 2*N_CSlice, Fdz );
         FIDz += 2*N_CSlice;

      } // for (int z=CIDz; z<POT_NXT-1; z+=N_CSlice)
   } // if ( ID < N_CSlice*(POT_NXT-2)*(POT_NXT-2) )
   __syncthreads();
   

// c. Assign interpolated potential to Rho_array as boundary condition
// -----------------------------------------------------------------------------------------------------------

   float bc_xm,bc_xp,bc_ym,bc_yp,bc_zm,bc_zp;
   unsigned int N = RHO_NXT;
         

   t = ID;
   stride = blockDim.x * blockDim.y;

   do 
   {
     if (t < RHO_NXT*RHO_NXT*RHO_NXT){
     Rhoid_x = t % RHO_NXT;
     Rhoid_y = t/RHO_NXT % RHO_NXT;
     Rhoid_z = t / (RHO_NXT*RHO_NXT);
	     
     // if boundary condition
     if (Rhoid_x==0)         {bc_xm = shared_mem[RHO_NXT*RHO_NXT * 0 + Rhoid_y + RHO_NXT * Rhoid_z];} //[k+2][j+2][im+2];}
     else bc_xm = 0.0;
     if (Rhoid_y==0)         {bc_ym = shared_mem[RHO_NXT*RHO_NXT * 1 + Rhoid_x + RHO_NXT * Rhoid_z];} //[k+2][j+2][im+2];}
     else bc_ym = 0.0;
     if (Rhoid_z==0)         {bc_zm = shared_mem[RHO_NXT*RHO_NXT * 2 + Rhoid_x + RHO_NXT * Rhoid_y];} //[k+2][j+2][im+2];}
     else bc_zm = 0.0;
	 
     if (Rhoid_x==RHO_NXT-1) {bc_xp = shared_mem[RHO_NXT*RHO_NXT * 3 + Rhoid_y + RHO_NXT * Rhoid_z];} //[k+2][j+2][ip+2];}
     else bc_xp = 0.0;
     if (Rhoid_y==RHO_NXT-1) {bc_yp = shared_mem[RHO_NXT*RHO_NXT * 4 + Rhoid_x + RHO_NXT * Rhoid_z];} //[k+2][j+2][ip+2];}
     else bc_yp = 0.0;
     if (Rhoid_z==RHO_NXT-1) {bc_zp = shared_mem[RHO_NXT*RHO_NXT * 5 + Rhoid_x + RHO_NXT * Rhoid_y];} //[k+2][j+2][ip+2];}
     else bc_zp = 0.0;

     g_Rho_Array[bid][t] *= -Const; 
     g_Rho_Array[bid][t] += bc_xm + bc_ym + bc_zm + bc_xp + bc_yp + bc_zp;

     }
     t += stride;    

   } while (t < RHO_NXT*RHO_NXT*RHO_NXT);
   __syncthreads();
   

   //  ----------------------------------------------   Discret Sine Transform Poisson Solver -- START -----------------------------------------------//

   //  Forward FFT
   DST_Scheme(ID, g_Rho_Array[bid],0);
   DST_Scheme(ID, g_Rho_Array[bid],1);
   DST_Scheme(ID, g_Rho_Array[bid],2);

    // Poisson Eigen
  __shared__ real Eigen[RHO_NXT];
  for (int i=0; i<RHO_NXT; i++)
    Eigen[i] = 1.-COS(M_PI*(i+1)/(RHO_NXT+1));

   __syncthreads();
   
   t=ID;
   do
     {
     Rhoid_z = t /(RHO_NXT*RHO_NXT);
     Rhoid_y = (t /RHO_NXT) % RHO_NXT;
     Rhoid_x = t % RHO_NXT;
     g_Rho_Array[bid][t] /=  2. * (Eigen[Rhoid_x] + Eigen[Rhoid_y] + Eigen[Rhoid_z]) * (8*(N+1)*(N+1)*(N+1));

     t +=stride;

     } while (t < RHO_NXT*RHO_NXT*RHO_NXT);
     __syncthreads();
 
   // // reverse FFT  
   DST_Scheme(ID, g_Rho_Array[bid],0);
   DST_Scheme(ID, g_Rho_Array[bid],1);
   DST_Scheme(ID, g_Rho_Array[bid],2);
   
  //  ----------------------------------------------   Discret Sine Transform Poisson Solver -- END  -----------------------------------------------//

   t=ID;
   do
   {
     Rhoid_z = t /(GRA_NXT*GRA_NXT);
     Rhoid_y = (t /GRA_NXT) % GRA_NXT;
     Rhoid_x = t % GRA_NXT;

     Rhoid    =  (Rhoid_z + 2) * (RHO_NXT) * (RHO_NXT)
               + (Rhoid_y + 2) * (RHO_NXT)
               + (Rhoid_x + 2) ;

     if (t< GRA_NXT * GRA_NXT * GRA_NXT){
       g_Pot_Array_Out[bid][t] = g_Rho_Array[bid][Rhoid] ; /// (8*(N+1)*(N+1)*(N+1));
     }
     t += stride;

   } while (t< GRA_NXT * GRA_NXT * GRA_NXT);
    __syncthreads();
} // FUNCTION : CUPOT_PoissonSolver_SOR



#endif // #if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == SOR )
