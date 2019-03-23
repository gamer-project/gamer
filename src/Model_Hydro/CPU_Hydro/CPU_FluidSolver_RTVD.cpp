#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == RTVD )



// check before compiling anything else
#if ( NCOMP_PASSIVE != 0 )
#  error : RTVD scheme does NOT support passive scalars !!
#endif


#define to1D(z,y,x) ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

real Hydro_CheckMinPres( const real InPres, const real MinPres );
real Hydro_CheckMinPresInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                               const real Gamma_m1, const real _Gamma_m1, const real MinPres );

static void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], const real dt, const real dx, const real Gamma,
                          const bool StoreFlux, const int j_gap, const int k_gap, const real MinDens, const real MinPres );
static void TransposeXY( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] );
static void TransposeXZ( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver_RTVD
// Description :  CPU fluid solver based on the relaxing TVD (RTVD) scheme
//
// Note        :  The three-dimensional evolution is achieved by using the dimensional-split method
//                --> Use the input pamameter "XYZ" to control the order of update
//
// Parameter   :  Flu_Array_In   : Array storing the input fluid variables
//                Flu_Array_Out  : Array to store the output fluid variables
//                Flux_Array     : Array to store the output flux
//                Corner_Array   : Array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                Pot_Array_USG  : Array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in RTVD)
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Gamma          : Ratio of specific heats
//                StoreFlux      : true --> store the coarse-fine fluxes
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                MinDens/Pres   : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver_RTVD(
   real Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   real Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   real Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double Corner_Array[][3],
   const real Pot_Array_USG[][ CUBE(USG_NXT_F) ],
   const int NPatchGroup, const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const bool XYZ, const real MinDens, const real MinPres )
{

   if ( XYZ )
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         CPU_AdvanceX( Flu_Array_In[P], dt, dh, Gamma, StoreFlux,              0,              0, MinDens, MinPres );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], dt, dh, Gamma, StoreFlux, FLU_GHOST_SIZE,              0, MinDens, MinPres );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], dt, dh, Gamma, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE, MinDens, MinPres );

         TransposeXZ ( Flu_Array_In[P] );
         TransposeXY ( Flu_Array_In[P] );
      }
   }

   else
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         TransposeXY ( Flu_Array_In[P] );
         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], dt, dh, Gamma, StoreFlux,              0,              0, MinDens, MinPres );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], dt, dh, Gamma, StoreFlux,              0, FLU_GHOST_SIZE, MinDens, MinPres );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], dt, dh, Gamma, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE, MinDens, MinPres );
      }
   }


// copy the updated fluid variables to Flu_Array_Out
   int ID1, ID2, ii, jj, kk;

#  pragma omp parallel for private( ID1, ID2, ii, jj, kk ) schedule( runtime )
   for (int P=0; P<NPatchGroup; P++)   {
   for (int v=0; v<5; v++)             {
   for (int k=0; k<PS2; k++)           {  kk = k + FLU_GHOST_SIZE;
   for (int j=0; j<PS2; j++)           {  jj = j + FLU_GHOST_SIZE;
   for (int i=0; i<PS2; i++)           {  ii = i + FLU_GHOST_SIZE;

      ID1 = k*PS2*PS2 + j*PS2 + i;
      ID2 = to1D(kk,jj,ii);

      Flu_Array_Out[P][v][ID1] = Flu_Array_In[P][v][ID2];

   }}}}}


// copy the coarse-fine fluxes into Flux_Array
   int mm, nn;

   if ( StoreFlux ) {
#  pragma omp parallel for private( ID1, mm, nn ) schedule( runtime )
   for (int P=0; P<NPatchGroup; P++)
   {
      for (int v=0; v<5; v++)    {
      for (int m=0; m<PS2; m++)  {  mm = m + FLU_GHOST_SIZE;
      for (int n=0; n<PS2; n++)  {  nn = n + FLU_GHOST_SIZE;

         ID1 = m*PS2 + n;

//       left yz plane
         Flux_Array[P][0][v][ID1] = Flu_Array_In[P][v][(       mm)*FLU_NXT*FLU_NXT+(       nn)*FLU_NXT+        2];

//       central yz plane
         Flux_Array[P][1][v][ID1] = Flu_Array_In[P][v][(       mm)*FLU_NXT*FLU_NXT+(       nn)*FLU_NXT+        0];

//       right yz plane
         Flux_Array[P][2][v][ID1] = Flu_Array_In[P][v][(       mm)*FLU_NXT*FLU_NXT+(       nn)*FLU_NXT+FLU_NXT-3];

//       left xz plane
         Flux_Array[P][3][v][ID1] = Flu_Array_In[P][v][(       mm)*FLU_NXT*FLU_NXT+(        2)*FLU_NXT+       nn];

//       central xz plane
         Flux_Array[P][4][v][ID1] = Flu_Array_In[P][v][(       mm)*FLU_NXT*FLU_NXT+(        0)*FLU_NXT+       nn];

//       right xz plane
         Flux_Array[P][5][v][ID1] = Flu_Array_In[P][v][(       mm)*FLU_NXT*FLU_NXT+(FLU_NXT-3)*FLU_NXT+       nn];

//       left xy plane
         Flux_Array[P][6][v][ID1] = Flu_Array_In[P][v][(        2)*FLU_NXT*FLU_NXT+(       mm)*FLU_NXT+       nn];

//       central xy plane
         Flux_Array[P][7][v][ID1] = Flu_Array_In[P][v][(        0)*FLU_NXT*FLU_NXT+(       mm)*FLU_NXT+       nn];

//       right xy plane
         Flux_Array[P][8][v][ID1] = Flu_Array_In[P][v][(FLU_NXT-3)*FLU_NXT*FLU_NXT+(       mm)*FLU_NXT+       nn];

      }}}
   } // for (int P=0; P<NPatchGroup; P++)
   } // if ( StoreFlux )

} // FUNCTION : CPU_FluidSolver_RTVD



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_AdvanceX
// Description :  Use CPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on the TVD scheme
//
// Parameter   :  u            : Input fluid array
//                dt           : Time interval to advance solution
//                dx           : Grid size
//                Gamma        : Ratio of specific heats
//                StoreFlux    : true --> store the coarse-fine fluxes
//                j_gap        : Number of cells that can be skipped on each side in the y direction
//                k_gap        : Number of cells that can be skipped on each side in the z direction
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_AdvanceX( real u[][ CUBE(FLU_NXT) ], const real dt, const real dx, const real Gamma,
                   const bool StoreFlux, const int j_gap, const int k_gap, const real MinDens, const real MinPres )
{

   const real  Gamma_m1 = Gamma - (real)1.0;    // for evaluating pressure
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
   const real _dx       = (real)1.0/dx;         // one over dx
   const real dt_half   = (real)0.5*dt;         // for evaluating u_half

   const int j_start    = j_gap;
   const int k_start    = k_gap;
   const int j_end      = FLU_NXT-j_gap;
   const int k_end      = FLU_NXT-k_gap;

// set local variables
   real ux     [5][FLU_NXT];              // one column of u in x direction
   real u_half [5][FLU_NXT];              // u in the midpoint
   real flux   [5][FLU_NXT];              // flux defined in the right-hand surface of cell
   real cu     [5][FLU_NXT];              // freezing speed c * u
   real cw     [5][FLU_NXT];              // freezing speed c * w ( == flux defined in the center of cell )
   real RLflux [5][FLU_NXT];              // right/left-moving flux ( defined in the right-hand surface of cell )
   real c;                                // freezing speed ( == abs(vx) + sound speed )
   real _rho, vx, p;                      // one over rho, velocity in x direction, pressure

   int ip, im;
   real Temp;


   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {

//    copy one column of data from u to ux
      for (int v=0; v<5; v++)    memcpy( ux[v], &u[v][to1D(k,j,0)], FLU_NXT*sizeof(real) );


//    a. Evaluate the half-step values of fluid variables
//-----------------------------------------------------------------------------

//    (a1). set variables defined in the center of cell
      for (int i=0; i<FLU_NXT; i++)
      {
         _rho = (real)1.0 / ux[0][i];
         vx   = _rho * ux[1][i];
         p    = Gamma_m1 * ( ux[4][i] - (real)0.5*_rho*( ux[1][i]*ux[1][i] + ux[2][i]*ux[2][i] +
                                                         ux[3][i]*ux[3][i] ) );
         p    = Hydro_CheckMinPres( p, MinPres );

#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( Hydro_CheckNegative(p) )
            Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         p, __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(ux[0][i]) )
            Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         ux[0][i], __FILE__, __LINE__, __FUNCTION__ );
#        endif

         c    = FABS( vx ) + SQRT( Gamma*p*_rho );

         cw[0][i] = ux[1][i];
         cw[1][i] = ux[1][i] * vx + p;
         cw[2][i] = ux[2][i] * vx;
         cw[3][i] = ux[3][i] * vx;
         cw[4][i] = ( ux[4][i] + p ) * vx;

         cu[0][i] = c*ux[0][i];
         cu[1][i] = c*ux[1][i];
         cu[2][i] = c*ux[2][i];
         cu[3][i] = c*ux[3][i];
         cu[4][i] = c*ux[4][i];
      }


//    (a2). set flux defined in the right-hand surface of cell by the upwind scheme
      for (int v=0; v<5; v++)
      for (int i=0; i<FLU_NXT-1; i++)
      {
         ip = i+1;
         flux[v][i] = (real)0.5*(  ( cu[v][i]+cw[v][i] ) - ( cu[v][ip]-cw[v][ip] )  );
      }


//    (a3). evaluate the intermidiate values (u_half)
      for (int v=0; v<5; v++)
      for (int i=1; i<FLU_NXT-1; i++)
      {
         im = i-1;
         u_half[v][i] = ux[v][i] - _dx*dt_half*( flux[v][i]-flux[v][im] ) ;
      }


//    (a4). enforce positive density and pressure
      for (int i=1; i<FLU_NXT-1; i++)
      {
         u_half[0][i] = FMAX( u_half[0][i], MinDens );
         u_half[4][i] = Hydro_CheckMinPresInEngy( u_half[0][i], u_half[1][i], u_half[2][i], u_half[3][i], u_half[4][i],
                                                  Gamma_m1, _Gamma_m1, MinPres );
      }



//    b. Evaluate the full-step values of fluid variables
//-----------------------------------------------------------------------------

//    (b1). reset variables defined in the center of cell at the intermidate state
      for (int i=1; i<FLU_NXT-1; i++)
      {
         _rho = (real)1.0 / u_half[0][i];
         vx   = _rho * u_half[1][i];
         p    = Gamma_m1 * (  u_half[4][i] - (real)0.5*_rho*(  u_half[1][i]*u_half[1][i] +
                                                               u_half[2][i]*u_half[2][i] +
                                                               u_half[3][i]*u_half[3][i] )  );
         p    = Hydro_CheckMinPres( p, MinPres );

#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( Hydro_CheckNegative(p) )
            Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         p, __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(u_half[0][i]) )
            Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         u_half[0][i], __FILE__, __LINE__, __FUNCTION__ );
#        endif

         c    = FABS( vx ) + SQRT( Gamma*p*_rho );

         cw[0][i] = u_half[1][i];
         cw[1][i] = u_half[1][i] * vx + p;
         cw[2][i] = u_half[2][i] * vx;
         cw[3][i] = u_half[3][i] * vx;
         cw[4][i] = ( u_half[4][i] + p ) * vx;

         cu[0][i] = c*u_half[0][i];
         cu[1][i] = c*u_half[1][i];
         cu[2][i] = c*u_half[2][i];
         cu[3][i] = c*u_half[3][i];
         cu[4][i] = c*u_half[4][i];
      }


//    (b2). set the right-moving flux defined in the right-hand surface by the TVD scheme
      for (int v=0; v<5; v++)
      for (int i=1; i<FLU_NXT-2; i++)
         RLflux[v][i] = (real)0.5*( cu[v][i] + cw[v][i] );


      for (int v=0; v<5; v++)
      for (int i=2; i<FLU_NXT-3; i++)
      {
         im = i-1; ip = i+1;

         flux[v][i] = RLflux[v][i];

         Temp = ( RLflux[v][ip]-RLflux[v][i] ) * ( RLflux[v][i]-RLflux[v][im] );

         if ( Temp > (real)0.0 )    flux[v][i] += Temp / ( RLflux[v][ip]-RLflux[v][im] );
      }


//    (b3). set the left-moving flux defined in the left-hand surface by the TVD scheme, get the total flux
      for (int v=0; v<5; v++)
      for (int i=1; i<FLU_NXT-2; i++)
      {
         ip = i+1;
         RLflux[v][i] = (real)0.5*( cu[v][ip] - cw[v][ip] );
      }

      for (int v=0; v<5; v++)
      for (int i=2; i<FLU_NXT-3; i++)
      {
         im = i-1;
         ip = i+1;

         flux[v][i] -= RLflux[v][i];

         Temp = ( RLflux[v][im]-RLflux[v][i] ) * ( RLflux[v][i]-RLflux[v][ip] );

         if ( Temp > (real)0.0 )    flux[v][i] -= Temp / ( RLflux[v][im]-RLflux[v][ip] );
      }


//    (b4). advance fluid by one full time-step
      for (int v=0; v<5; v++)
      for (int i=3; i<FLU_NXT-3; i++)
      {
         im = i-1;
         ux[v][i] -= _dx*dt*( flux[v][i] - flux[v][im] );
      }


//    (b5). enforce positive density and pressure
      for (int i=3; i<FLU_NXT-3; i++)
      {
         ux[0][i] = FMAX( ux[0][i], MinDens );
         ux[4][i] = Hydro_CheckMinPresInEngy( ux[0][i], ux[1][i], ux[2][i], ux[3][i], ux[4][i], Gamma_m1, _Gamma_m1, MinPres );
      }


//    (b6). check negative density and energy
#     ifdef CHECK_NEGATIVE_IN_FLUID
      for (int i=3; i<FLU_NXT-3; i++)
      {
         if ( Hydro_CheckNegative(ux[0][i]) )
            Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         ux[0][i], __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(ux[4][i]) )
            Aux_Message( stderr, "ERROR : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         ux[4][i], __FILE__, __LINE__, __FUNCTION__ );
      }
#     endif



//    c. Save results
//-----------------------------------------------------------------------------

//    (c1). save the final result back to array u
      for (int v=0; v<5; v++)    memcpy( &u[v][to1D(k,j,3)], ux[v]+3, (FLU_NXT-6)*sizeof(real) );


//    (c2). save the flux required by the flux-correction operation
      if ( StoreFlux )
      if (  ( j>=3 && j<FLU_NXT-3 ) && ( k>=3 && k<FLU_NXT-3 )  )
      {
         for (int v=0; v<5; v++)
         {
            u[v][ to1D(k,j,        2) ] = flux[v][          2];
            u[v][ to1D(k,j,FLU_NXT-3) ] = flux[v][FLU_NXT - 4];
            u[v][ to1D(k,j,        0) ] = flux[v][FLU_NXT/2-1];
         }
      }

   } // for (int j=j_start; j<j_end; j++) ... for (int k=k_start; k<k_end; k++)

} // FUNCTION : CPU_AdvanceX



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXY
// Description :  Transpose the x and y components
//
// Parameter   :  u : Input fluid array
//-------------------------------------------------------------------------------------------------------
void TransposeXY( real u[][ CUBE(FLU_NXT) ] )
{

   real (*u_xy)[FLU_NXT*FLU_NXT] = new real [5][FLU_NXT*FLU_NXT];
   int ID;

   for (int k=0; k<FLU_NXT; k++)
   {
      for (int j=0; j<FLU_NXT; j++)
      for (int i=0; i<FLU_NXT; i++)
      {
         ID = to1D(k,j,i);

         u_xy[0][j+i*FLU_NXT] = u[0][ID];
         u_xy[1][j+i*FLU_NXT] = u[2][ID];
         u_xy[2][j+i*FLU_NXT] = u[1][ID];
         u_xy[3][j+i*FLU_NXT] = u[3][ID];
         u_xy[4][j+i*FLU_NXT] = u[4][ID];
      }

      for (int v=0; v<5; v++)    memcpy( &u[v][to1D(k,0,0)], u_xy[v], FLU_NXT*FLU_NXT*sizeof(real) );
   }

   delete [] u_xy;

} // FUNCTION : TrasposeXY



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXZ
// Description :  Transpose the x and z components
//
// Parameter   :  u : Input fluid array
//-------------------------------------------------------------------------------------------------------
void TransposeXZ( real u[][ CUBE(FLU_NXT) ] )
{

   real u_temp[5];
   int ID1, ID2;

   for (int j=0; j<FLU_NXT; j++)
   for (int k=0; k<FLU_NXT; k++)
   {
      for (int i=0; i<k; i++)
      {
         ID1 = to1D(k,j,i);
         ID2 = to1D(i,j,k);

         u_temp[0] = u[0][ID1];
         u_temp[1] = u[3][ID1];
         u_temp[2] = u[2][ID1];
         u_temp[3] = u[1][ID1];
         u_temp[4] = u[4][ID1];

         u[0][ID1] = u[0][ID2];
         u[1][ID1] = u[3][ID2];
         u[2][ID1] = u[2][ID2];
         u[3][ID1] = u[1][ID2];
         u[4][ID1] = u[4][ID2];

         u[0][ID2] = u_temp[0];
         u[1][ID2] = u_temp[1];
         u[2][ID2] = u_temp[2];
         u[3][ID2] = u_temp[3];
         u[4][ID2] = u_temp[4];
      }

      ID1 = to1D(k,j,k);

      u_temp[0] = u[0][ID1];
      u_temp[1] = u[3][ID1];
      u_temp[2] = u[2][ID1];
      u_temp[3] = u[1][ID1];
      u_temp[4] = u[4][ID1];

      u[0][ID1] = u_temp[0];
      u[1][ID1] = u_temp[1];
      u[2][ID1] = u_temp[2];
      u[3][ID1] = u_temp[3];
      u[4][ID1] = u_temp[4];

   } // j,k

} // FUNCTION : TrasposeXZ



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == RTVD )
