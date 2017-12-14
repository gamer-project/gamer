#include "GAMER.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == WAF )


// check before compiling anything else
#if ( NCOMP_PASSIVE != 0 )
#  error : WAF scheme does NOT support passive scalars !!
#endif


#define to1D(z,y,x) ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

extern real CPU_CheckMinPres( const real InPres, const real MinPres );

static real set_limit( const real r, const real c, const WAF_Limiter_t WAF_Limiter );
static void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], real fc[PS2*PS2][3][5], const real dt,
                          const real dx, const real Gamma, const int j_gap, const int k_gap,
                          const WAF_Limiter_t WAF_Limiter, const real MinDens, const real MinPres );
static void Store_flux( real L_f[5][PS2*PS2], real C_f[5][PS2*PS2], real R_f[5][PS2*PS2],
                        const real fc[PS2*PS2][3][5] );
static void Solve_Flux( real flux[5], const real lL_star[5], const real lR_star[5], const real cL_star[5],
                        const real cR_star[5], const real rL_star[5], const real rR_star[5], const real eival[5],
                        const real L_2[5], const real L_1[5], const real R_1[5], const real R_2[5],
                        const real Gamma, const real ratio, const WAF_Limiter_t WAF_Limiter );
static void set_flux( real flux[5], const real val[5], const real Gamma );
static void Change_variable( real pval[5], const real cval[5], const real Gamma, const real MinPres );
static void TransposeXY( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] );
static void TransposeXZ( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] );
static void TransposeFluxXY( real fc[PS2*PS2][3][5] );
static void TransposeFluxXZ( real fc[PS2*PS2][3][5] );
#ifdef WAF_DISSIPATE
static void Dis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5], const real R_star[5],
                      const real limit[5], const real theta[5], const real Gamma );
#else
static void Undis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5],
                        const real R_star[5], const real limit[5], const real theta[5], const real Gamma );
#endif
#if   ( RSOLVER == EXACT )
extern void CPU_RiemannSolver_Exact( const int XYZ, real eival_out[], real L_star_out[], real R_star_out[],
                                     real Flux_Out[], const real L_In[], const real R_In[], const real Gamma );
#elif ( RSOLVER == ROE )
static void Solve_StarRoe( real eival[5], real L_star[5], real R_star[5], const real L[5], const real R[5],
                           const real Gamma, const real MinPres );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver_WAF
// Description :  CPU fluid solver based on the Weighted-Average-Flux (WAF) scheme
//
// Note        :  The three-dimensional evolution is achieved by using the dimensional-split method
//                --> Use the input pamameter "XYZ" to control the order of update
//
// Parameter   :  Flu_Array_In   : Array to store the input fluid variables
//                Flu_Array_Out  : Array to store the output fluid variables
//                Flux_Array     : Array to store the output flux
//                Corner_Array   : Array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                Pot_Array_USG  : Array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in WAF)
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Gamma          : Ratio of specific heats
//                StoreFlux      : true --> store the coarse-fine fluxes
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                WAF_Limiter    : Selection of the limit function
//                                    0 : superbee
//                                    1 : van-Leer
//                                    2 : van-Albada
//                                    3 : minbee
//                MinDens/Pres   : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver_WAF( real Flu_Array_In [][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                          real Flu_Array_Out[][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                          real Flux_Array[][9][NCOMP_TOTAL][ PS2*PS2 ],
                          const double Corner_Array[][3],
                          const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                          const int NPatchGroup, const real dt, const real dh, const real Gamma,
                          const bool StoreFlux, const bool XYZ, const WAF_Limiter_t WAF_Limiter,
                          const real MinDens, const real MinPres )
{

#  pragma omp parallel
   {
      real (*FC)[3][5]= new real [PS2*PS2][3][5];


      if ( XYZ ) // x->y->z
      {
#        pragma omp for schedule( runtime )
         for (int P=0; P<NPatchGroup; P++)
         {
//          solve the x direction
            CPU_AdvanceX( Flu_Array_In[P], FC, dt, dh, Gamma,              0,              0, WAF_Limiter, MinDens, MinPres );

//          store the intercell flux in x direction
            if ( StoreFlux )
            Store_flux( Flux_Array[P][0], Flux_Array[P][1], Flux_Array[P][2], FC );

//          x-y-z --> y-x-z for conservative variables
            TransposeXY ( Flu_Array_In[P] );

//          solve the y direction
            CPU_AdvanceX( Flu_Array_In[P], FC, dt, dh, Gamma, FLU_GHOST_SIZE,              0, WAF_Limiter, MinDens, MinPres );

//          y-x-z --> x-y-z for flux
            TransposeFluxXY( FC );

//          store the intercell flux in y direction
            if ( StoreFlux )
            Store_flux( Flux_Array[P][3], Flux_Array[P][4], Flux_Array[P][5], FC );

//          y-x-z --> z-x-y for conservative variables
            TransposeXZ ( Flu_Array_In[P] );

//          solve the z direction
            CPU_AdvanceX( Flu_Array_In[P], FC, dt, dh, Gamma, FLU_GHOST_SIZE, FLU_GHOST_SIZE, WAF_Limiter, MinDens, MinPres );

//          z-x-y --> x-y-z for flux
            TransposeFluxXZ( FC );
            TransposeFluxXY( FC );

//          store intercell flux of z direction
            if ( StoreFlux )
            Store_flux( Flux_Array[P][6], Flux_Array[P][7], Flux_Array[P][8], FC );

//          z-x-y --> x-y-z for conservative variables
            TransposeXZ ( Flu_Array_In[P] );
            TransposeXY ( Flu_Array_In[P] );

         } // for (int P=0; P<NPatchGroup; P++)
      } // if ( XYZ )


      else // z->y->x
      {
#        pragma omp for schedule( runtime )
         for (int P=0; P<NPatchGroup; P++)
         {
//          x-y-z --> z-x-y for conservative variables
            TransposeXY ( Flu_Array_In[P] );
            TransposeXZ ( Flu_Array_In[P] );

//          solve the z direction
            CPU_AdvanceX( Flu_Array_In[P], FC, dt, dh, Gamma,              0,              0, WAF_Limiter, MinDens, MinPres );

//          z-x-y --> x-y-z for flux
            TransposeFluxXZ( FC );
            TransposeFluxXY( FC );

//          store the intercell flux of z direction
            if ( StoreFlux )
            Store_flux( Flux_Array[P][6], Flux_Array[P][7], Flux_Array[P][8], FC );

//          z-x-y --> y-x-z for conservative variables
            TransposeXZ ( Flu_Array_In[P] );

//          solve the y direction
            CPU_AdvanceX( Flu_Array_In[P], FC, dt, dh, Gamma,              0, FLU_GHOST_SIZE, WAF_Limiter, MinDens, MinPres );

//          y-x-z --> x-y-z for flux
            TransposeFluxXY( FC );

//          store the intercell flux of y direction
            if ( StoreFlux )
            Store_flux( Flux_Array[P][3], Flux_Array[P][4], Flux_Array[P][5], FC );

//          y-x-z --> x-y-z for conservative variables
            TransposeXY ( Flu_Array_In[P] );

//          solve the x direction
            CPU_AdvanceX( Flu_Array_In[P], FC, dt, dh, Gamma, FLU_GHOST_SIZE, FLU_GHOST_SIZE, WAF_Limiter, MinDens, MinPres );

//          store the intercell flux of x direction
            if ( StoreFlux )
            Store_flux( Flux_Array[P][0], Flux_Array[P][1], Flux_Array[P][2], FC );

         } // for (int P=0; P<NPatchGroup; P++)
      } // if ( XYZ ) ... else ...


//    copy the updated fluid variables into the output array "Flu_Array_Out"
      int ID1, ID2, ii, jj, kk;

#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)   {
      for (int v=0; v<5; v++)             {
      for (int k=0; k<PS2; k++)           {  kk = k + FLU_GHOST_SIZE;
      for (int j=0; j<PS2; j++)           {  jj = j + FLU_GHOST_SIZE;
      for (int i=0; i<PS2; i++)           {  ii = i + FLU_GHOST_SIZE;

         ID1 = k*PS2*PS2 + j*PS2 + i;
         ID2 = to1D(kk,jj,ii);

         Flu_Array_Out[P][v][ID1] = Flu_Array_In[P][v][ID2];

      }}}}}


      delete [] FC;

   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_WAF



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_AdvanceX
// Description :  Use CPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on the WAF scheme
//
// Parameter   :  u            : Input fluid array
//                fc           : Intercell flux
//                dt           : Time interval to advance solution
//                dx           : Grid size
//                Gamma        : Ratio of specific heats
//                j_gap        : Number of cells that can be skipped on each side in the y direction
//                k_gap        : Number of cells that can be skipped on each side in the z direction
//                WAF_Limiter  : Selection of te limit function
//                               0 : superbee
//                               1 : van-Leer
//                               2 : van-Albada
//                               3 : minbee
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], real fc[PS2*PS2][3][5], const real dt, const real dx,
                   const real Gamma, const int j_gap, const int k_gap, const WAF_Limiter_t WAF_Limiter,
                   const real MinDens, const real MinPres )
{

   const real  Gamma_m1 = Gamma - (real)1.0;    // for evaluating pressure
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
   const real ratio     = dt/dx;                // dt over dx
   const int j_start    = j_gap;
   const int k_start    = k_gap;
   const int j_end      = FLU_NXT-j_gap;
   const int k_end      = FLU_NXT-k_gap;

// set local variables
   real ux  [FLU_NXT][5];              // one column of u in x direction
   real flux[PS2+1][5];                // flux defined in the right-hand surface of cell
   real L_st[PS2+3][5];
   real R_st[PS2+3][5];
   real eval[PS2+3][5];

   int ii, ID;


   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {

//    copy one column of data from u to ux
      for (int i=0; i<FLU_NXT; i++)
      {
         ID = to1D(k,j,i);

         for (int v=0; v<5; v++)    ux[i][v] = u[v][ID];
      }

//    solve Riemann problem
      for (int i=0; i<PS2+3; i++)
      {
         real c_L[5], c_R[5];

         Change_variable( c_L, ux[i]  , Gamma, MinPres );
         Change_variable( c_R, ux[i+1], Gamma, MinPres );

#        if   ( RSOLVER == EXACT )
         CPU_RiemannSolver_Exact( 0, eval[i], L_st[i], R_st[i], NULL, c_L, c_R, Gamma );
#        elif ( RSOLVER == ROE )
         Solve_StarRoe( eval[i], L_st[i], R_st[i], c_L, c_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
#        endif
      }

//    solve intecell flux
      for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2+1; i++)
      {
         real c_L2[5], c_L1[5], c_R1[5], c_R2[5];

         Change_variable( c_L2, ux[i-2], Gamma, MinPres );
         Change_variable( c_L1, ux[i-1], Gamma, MinPres );
         Change_variable( c_R1, ux[i]  , Gamma, MinPres );
         Change_variable( c_R2, ux[i+1], Gamma, MinPres );

         int ii    = i - FLU_GHOST_SIZE;
         int ii_p1 = ii + 1;
         Solve_Flux( flux[ii], L_st[ii], R_st[ii], L_st[ii_p1], R_st[ii_p1], L_st[i], R_st[i], eval[ii_p1],
                     c_L2, c_L1, c_R1, c_R2, Gamma, ratio, WAF_Limiter );

//       paste the intercell flux to fc
         switch ( i )
         {
            case FLU_GHOST_SIZE : // the left intercell flux
            {
               if (  j >= FLU_GHOST_SIZE  &&  j < FLU_GHOST_SIZE+PS2  &&  k >= FLU_GHOST_SIZE  &&
                     k < FLU_GHOST_SIZE+PS2 )
               {
                  int ID1 = ( k - FLU_GHOST_SIZE )*PS2 + ( j - FLU_GHOST_SIZE );
                  for (int v=0; v<5; v++)    fc[ID1][0][v] = flux[ i - FLU_GHOST_SIZE ][v];
               }
               break;
            }

            case FLU_GHOST_SIZE + PS2/2 : // the center intercell flux
            {
               if (  j >= FLU_GHOST_SIZE  &&  j < FLU_GHOST_SIZE+PS2  &&  k >= FLU_GHOST_SIZE  &&
                     k < FLU_GHOST_SIZE+PS2  )
               {
                  int ID1 = ( k - FLU_GHOST_SIZE )*PS2 + ( j - FLU_GHOST_SIZE );
                  for (int v=0; v<5 ; v++)   fc[ID1][1][v] = flux[ i - FLU_GHOST_SIZE ][v];
               }
               break;
            }

            case FLU_GHOST_SIZE + PS2 : //the right intercell flux
            {
               if (  j >= FLU_GHOST_SIZE  &&  j < FLU_GHOST_SIZE+PS2  &&  k >= FLU_GHOST_SIZE  &&
                     k < FLU_GHOST_SIZE+PS2  )
               {
                  int ID1 = ( k - FLU_GHOST_SIZE )*PS2 + ( j - FLU_GHOST_SIZE );
                  for (int v=0; v<5; v++)    fc[ID1][2][v] = flux[ i - FLU_GHOST_SIZE ][v];
               }
               break;
            }

            default : /* nothing to do */
               break;
         }
      }

//    update the conservative variables
      for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2; i++)
      {
         ii = i - FLU_GHOST_SIZE;
         for (int v=0; v<5; v++)    ux[i][v] += ratio*( flux[ii][v] - flux[ii+1][v] );

//       enforce positive density and pressure
         ux[i][0] = FMAX( ux[i][0], MinDens );
         ux[i][4] = CPU_CheckMinPresInEngy( ux[i][0], ux[i][1], ux[i][2], ux[i][3], ux[i][4], Gamma_m1, _Gamma_m1, MinPres );

//       check negative density and energy
#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( CPU_CheckNegative(ux[i][0]) )
            Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         ux[i][0], __FILE__, __LINE__, __FUNCTION__ );

         if ( CPU_CheckNegative(ux[i][4]) )
            Aux_Message( stderr, "ERROR : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                         ux[i][4], __FILE__, __LINE__, __FUNCTION__ );
#        endif

      } // for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2; i++)

//    paste to the original space
      for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2; i++)
      {
         ID = to1D(k,j,i);
         for (int v=0; v<5; v++)    u[v][ID] = ux[i][v];
      }
   } // k,j

} // FUNCTION : CPU_AdvanceX



//-------------------------------------------------------------------------------------------------------
// Function    : Store_flux
// Description : Store the intercell flux for AMR
//
// Parameter   : L_f : Left intercell flux
//               C_f : Center intercell flux
//               R_f : Right intercel flux
//               fc  : Intercell flux obtained by function solve_Flux
//-------------------------------------------------------------------------------------------------------
void Store_flux( real L_f[5][PS2*PS2], real C_f[5][PS2*PS2], real R_f[5][PS2*PS2], const real fc[PS2*PS2][3][5] )
{

   int ID1;

   for (int v=0; v<5 ; v++)
   for (int m=0; m<PS2; m++)
   for (int n=0; n<PS2; n++)
   {
      ID1 = m*PS2 + n;

      L_f[v][ID1] = fc[ID1][0][v];
      C_f[v][ID1] = fc[ID1][1][v];
      R_f[v][ID1] = fc[ID1][2][v];
   }

} // FUNCTION : Store_flux



//-------------------------------------------------------------------------------------------------------
// Function    : Solve_Flux
// Description : Solve the intercell flux
//
// Parameter   : flux         : Intercell flux
//               lL_star      : Primitive variables in the left star region of the left region
//               lR_star      : Primitive variables in the right star region of the left region
//               cL_star      : Primitive variables in the left star region of the centor region
//               cR_star      : Primitive variables in the right star region of the centor region
//               rL_star      : Primitive variables in the left star region of the right region
//               rR_star      : Primitive variables in the right star region of the right region
//               eival        : Eigenvalue
//               L_2          : Primitive variables in the region left to the left region
//               L_1          : Primitive variables in the left region
//               R_1          : Primitive variables in the right region
//               R_2          : Primitive variables in the region right to the right region
//               ratio        : dt over dx
//               Gamma        : Ratio of specific heats
//               WAF_Limiter  : Selection of te limit function
//                               0 : superbee
//                               1 : van-Leer
//                               2 : van-Albada
//                               3 : minbee
//-------------------------------------------------------------------------------------------------------
static void Solve_Flux( real flux[5], const real lL_star[5], const real lR_star[5], const real cL_star[5],
                        const real cR_star[5], const real rL_star[5], const real rR_star[5], const real eival[5],
                        const real L_2[5], const real L_1[5], const real R_1[5], const real R_2[5],
                        const real Gamma, const real ratio, const WAF_Limiter_t WAF_Limiter )
{

   real theta[5];            //the sign of sped of waves
   real limit[5];            //limit functions
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
      if (  FABS( eival[i] ) < MAX_ERROR  )  limit[i] = (real)1.0;
      else
      {
         if ( eival[i] > (real)0.0 )
         {
            if (  mean[0][i] == (real)0.0  ||  mean[1][i] == (real)0.0  )  limit[i] = (real)1.0;
            else
            {
               if (  ( delta[0][i]*delta[1][i] ) / ( mean[0][i]*mean[1][i] ) < MAX_ERROR*MAX_ERROR  )
                  limit[i] = (real)1.0;
               else
               {
                  real r   = delta[0][i] / delta[1][i];
                  limit[i] = set_limit( r, eival[i]*ratio, WAF_Limiter );
               }
            }
         }

         else
         {
            if (  mean[2][i] == (real)0.0  ||  mean[1][i] == (real)0.0  )     limit[i] = (real)1.0;
            else
            {
               if (  ( delta[2][i]*delta[1][i] ) / ( mean[2][i]*mean[1][i] ) < MAX_ERROR*MAX_ERROR )
                  limit[i] = (real)1.0;
               else
               {
                  real r   = delta[2][i] / delta[1][i];
                  limit[i] = set_limit( r, eival[i]*ratio, WAF_Limiter );
               }
            }
         }
      }
   } // for (int i=0; i<5; i++)


// solve the sign of waves
   for ( int i=0; i<5 ; i++)
   {
      if (  FABS( eival[i] ) < MAX_ERROR  )  theta[i] =  (real)0.0;
      else if ( eival[i] > (real)0.0 )       theta[i] =  (real)1.0;
      else                                   theta[i] = -(real)1.0;
   }


// solve the intercell flux
#  ifdef WAF_DISSIPATE
   Dis_Stru(   flux, L_1, R_1, cL_star, cR_star, limit, theta, Gamma );
#  else
   Undis_Stru( flux, L_1, R_1, cL_star, cR_star, limit, theta, Gamma );
#  endif

} // FUNCTION : Solve_Flux



//-----------------------------------------------------------------------------------------------------
// Function    : set_limit
// Description : set the limit function
//
// parameter   : r            : Fluid variable
//               c            : Courant number
//               WAF_Limiter  : Selection of te limit function
//                               0 : superbee
//                               1 : van-Leer
//                               2 : van-Albada
//                               3 : minbee
//-------------------------------------------------------------------------------------------------------
real set_limit( const real r, const real c, const WAF_Limiter_t WAF_Limiter )
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
         limit = (real)1.0 - r*( (real)1.0 + r )*((real)1.0 - FABS(c) )/ ( (real)1.0 + r*r );
         break;
      }

      case WAF_MINBEE :
      {
         if ( r > (real)0.0  &&  r <= (real)1.0 )     limit = (real)1.0 - r*( (real)1.0 - FABS(c) );
         else                                         limit = FABS(c);
         break;
      }

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "WAF_Limiter", WAF_Limiter );
   }

   return limit;

} // FUNCTION : set_limit



//-----------------------------------------------------------------------------------------------------
// Function    : Change_variable
// Description : Change conservative variables to primitive variables
//
// Parameter   : pval    : Primitive variables
//               cval    : Conservative variables
//               Gamma   : Ratio of specific heats
//               MinPres : Minimum allowed pressure
//------------------------------------------------------------------------------------------------------
void Change_variable( real pval[5], const real cval[5], const real Gamma, const real MinPres )
{

   const real Gamma_m1 = Gamma - (real)1.0;

   pval[0] = cval[0];
   pval[1] = cval[1] / cval[0];
   pval[2] = cval[2] / cval[0];
   pval[3] = cval[3] / cval[0];
   pval[4] = Gamma_m1*( cval[4] - (real)0.5*( cval[1]*cval[1] + cval[2]*cval[2] + cval[3]*cval[3] ) / cval[0] );

// ensure a positive pressure
   pval[4] = CPU_CheckMinPres( pval[4], MinPres );

} // FUNCTION : Change_variable



#ifdef WAF_DISSIPATE
//------------------------------------------------------------------------------------------------------
// Function    : Dis_Stru
// Description : Set the intercell flux by dissipative wave structure
//
// Parameter   : flux   : Intercel flux
//               L      : Primitive variables in the left region
//               R      : Primitive variables in the right region
//               L_star : Primitive variables in the left star region
//               R_star : Primitive variables in the right star region
//               limit  : Limit functions
//               theta  : Sign of wave speed
//               Gamma  : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void Dis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5], const real R_star[5],
               const real limit[5], const real theta[5], const real Gamma )
{

   real iflux[6][5];
   real lim[5];

   for (int i=0; i<5 ; i++)   lim[i] = limit[i];


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
       flux[i] = (real)0.5*( iflux[0][i] + iflux[5][i] ) -
                 (real)0.5*(  theta[0]*lim[0]*( iflux[1][i] - iflux[0][i] ) +
                              theta[1]*lim[1]*( iflux[2][i] - iflux[1][i] ) +
                              theta[2]*lim[2]*( iflux[3][i] - iflux[2][i] ) +
                              theta[3]*lim[3]*( iflux[4][i] - iflux[3][i] ) +
                              theta[4]*lim[4]*( iflux[5][i] - iflux[4][i] )   );
   }

} // FUNCTION : Dis_Stru
#endif // #ifdef WAF_DISSIPATE



#ifndef WAF_DISSIPATE
//------------------------------------------------------------------------------------------------------
// Function    : Undis_Stru
// Description : Set the intercell flux by non-dissipative wave structure
//
// Parameter   : flux   : Intercel flux
//               L      : Primitive variables in the left region
//               R      : Primitive variables in the right region
//               L_star : Primitive variables in the left star region
//               R_star : Primitive variables in the right star region
//               limit  : Limit functions
//               theta  : Sign of wave speed
//               Gamma  : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void Undis_Stru( real flux[5], const real L[5], const real R[5], const real L_star[5], const real R_star[5],
                 const real limit[5], const real theta[5], const real Gamma )
{

   real iflux[4][5];

// flux function evaluated at the given stat
   set_flux( iflux[0], L,      Gamma );
   set_flux( iflux[1], L_star, Gamma );
   set_flux( iflux[2], R_star, Gamma );
   set_flux( iflux[3], R,      Gamma );

// set the intercell flux
   flux[0] =    (real)0.5*( iflux[0][0] + iflux[3][0] )
              - (real)0.5*(  theta[0]*limit[0]*( iflux[1][0] - iflux[0][0] ) +
                             theta[1]*limit[1]*( iflux[2][0] - iflux[1][0] ) +
                             theta[4]*limit[4]*( iflux[3][0] - iflux[2][0] )   );

   flux[1] =    (real)0.5*( iflux[0][1] + iflux[3][1] )
              - (real)0.5*( theta[0]*limit[0]*( iflux[1][1] - iflux[0][1] ) +
                            theta[1]*limit[1]*( iflux[2][1] - iflux[1][1] ) +
                            theta[4]*limit[4]*( iflux[3][1] - iflux[2][1] )   );

   flux[4] =    (real)0.5*( iflux[0][4] + iflux[3][4] )
              - (real)0.5*( theta[0]*limit[0]*( iflux[1][4] - iflux[0][4] ) +
                            theta[1]*limit[1]*( iflux[2][4] - iflux[1][4] ) +
                            theta[4]*limit[4]*( iflux[3][4] - iflux[2][4] )   );

   flux[2] =    (real)0.5*( iflux[0][2] + iflux[3][2] )
              - (real)0.5*( theta[0]*limit[0]*( iflux[1][2] - iflux[0][2] ) +
                            theta[2]*limit[2]*( iflux[2][2] - iflux[1][2] ) +
                            theta[4]*limit[4]*( iflux[3][2] - iflux[2][2] )   );

   flux[3] =    (real)0.5*( iflux[0][3] + iflux[3][3] )
              - (real)0.5*( theta[0]*limit[0]*( iflux[1][3] - iflux[0][3] ) +
                            theta[3]*limit[3]*( iflux[2][3] - iflux[1][3] ) +
                            theta[4]*limit[4]*( iflux[3][3] - iflux[2][3] )   );

} // FUNCTION : Undis_Stru
#endif // #ifndef WAF_DISSIPATE



//-------------------------------------------------------------------------------------------------------
// Function    : set_flux
// Description : Set the flux function evaluated at the given stat
//
// Parameter   : flux   : Flux function
//               val    : Primitive variables
//               Gamma  : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void set_flux( real flux[5], const real val[5], const real Gamma )
{

   const real Gamma_m1 = Gamma -(real)1.0;

// set flux
   flux[0] = val[0]*val[1];
   flux[1] = val[0]*val[1]*val[1] + val[4];
   flux[2] = val[0]*val[1]*val[2];
   flux[3] = val[0]*val[1]*val[3];
   flux[4] = val[1]*( (real)0.5*val[0]*( val[1]*val[1] + val[2]*val[2] + val[3]*val[3] )
                      + val[4]/Gamma_m1 + val[4] );

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
void Solve_StarRoe( real eival[5], real L_star[5], real R_star[5], const real L[5], const real R[5],
                    const real Gamma, const real MinPres )
{

   const real Gamma_m1 = Gamma - (real)1.0; // for evaluating pressure and sound speed

   real u_bar, v_bar, w_bar, h_bar, a_bar, a_bar_inv; // Roe's average of vx, vy, vz, enthapy, sound speed,
                                                      // one over a_bar
   real coef[5]; //Roe's coefficients
   real TempPres, TempRho, _TempRho;

// solve Roe's average
   {
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CPU_CheckNegative(L[0]) )
         Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      L[0], __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(R[0]) )
         Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      R[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real n_L_sq = SQRT( L[0] ); // rooting number of left density
      real n_R_sq = SQRT( R[0] ); // rooting number of right density
      real h_L = (real)0.5*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] ) + Gamma/Gamma_m1*L[4]/L[0]; // left enthapy
      real h_R = (real)0.5*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] ) + Gamma/Gamma_m1*R[4]/R[0]; // right enthapy
      real n_bar_inv = (real)1.0 / ( n_L_sq + n_R_sq ); // one over (n_L_sq plus n_L_sq)

      u_bar = ( n_L_sq * L[1] + n_R_sq * R[1] ) * n_bar_inv;
      v_bar = ( n_L_sq * L[2] + n_R_sq * R[2] ) * n_bar_inv;
      w_bar = ( n_L_sq * L[3] + n_R_sq * R[3] ) * n_bar_inv;
      h_bar = ( n_L_sq * h_L  + n_R_sq * h_R  ) * n_bar_inv;

      real GammaP_Rho = Gamma_m1*(  h_bar - (real)0.5*( u_bar*u_bar + v_bar*v_bar + w_bar*w_bar )  );

//    enforce the pressure to be positive
      TempRho    = (real)0.5*( L[0] + R[0] );
      _TempRho   = (real)1.0/TempRho;
      TempPres   = GammaP_Rho*TempRho/Gamma;
      TempPres   = CPU_CheckMinPres( TempPres, MinPres );
      GammaP_Rho = Gamma*TempPres*_TempRho;

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CPU_CheckNegative(GammaP_Rho) )
         Aux_Message( stderr, "ERROR : negative GammaP_Rho (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      GammaP_Rho, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      a_bar      = SQRT( GammaP_Rho );
      a_bar_inv  = (real)1.0 / a_bar;
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
      coef[0] = (real)0.5*a_bar_inv*(  du_1*( u_bar + a_bar ) - du_2 - a_bar*coef[1]  );
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
      L_star[4] = (real)0.5*Gamma_m1*(  ( E_L - e_L_star + L[4]/Gamma_m1 + coef[0]*(h_bar - u_bar*a_bar) ) +
                                        ( E_R - e_R_star + R[4]/Gamma_m1 - coef[4]*(h_bar + u_bar*a_bar) )   );
      L_star[4] = CPU_CheckMinPres( L_star[4], MinPres );
      R_star[4] = L_star[4];
   }


// solve the speed of waves
   {
      real eigen[2];
      eival[1] = L_star[1];
      eival[2] = L_star[1];
      eival[3] = L_star[1];

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CPU_CheckNegative(L[4]) )
         Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      L[4],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(L[0]) )
         Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      L[0],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(L_star[4]) )
         Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      L_star[4], __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(L_star[0]) )
         Aux_Message( stderr, "ERROR : negative density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      L_star[0], __FILE__, __LINE__, __FUNCTION__ );

      if ( CPU_CheckNegative(R[4]) )
         Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      R[4],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(R[0]) )
         Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      R[0],      __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(R_star[4]) )
         Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      R_star[4], __FILE__, __LINE__, __FUNCTION__ );
      if ( CPU_CheckNegative(R_star[0]) )
         Aux_Message( stderr, "ERROR : negative density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
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



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXY
// Description :  Transpose the x and y components
//
// Parameter   :  u  : Input fluid array
//-------------------------------------------------------------------------------------------------------
void TransposeXY( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] )
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
// Parameter   :  u  : Input fluid array
//-------------------------------------------------------------------------------------------------------
void TransposeXZ( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] )
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

} // FUNCTION : TransposeXZ



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeFluxXY
// Description :  Transpose the x and y components for intercell flux
//
// Parameter   :  fc : Input flux array
//-------------------------------------------------------------------------------------------------------
void TransposeFluxXY( real fc[PS2*PS2][3][5] )
{

   real tmp;
   int ID;

   for (int i=0; i<3;   i++)
   for (int k=0; k<PS2; k++)
   for (int j=0; j<PS2; j++)
   {
      ID = k*PS2 + j;

      tmp          = fc[ID][i][2];
      fc[ID][i][2] = fc[ID][i][1];
      fc[ID][i][1] = tmp;
   }

} // FUNCTION : TrasposeFluxXY



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeFluxXZ
// Description :  Transpose the x and z components for intercell flux
//
// Parameter   :  fc : Input flux array
//-------------------------------------------------------------------------------------------------------
void TransposeFluxXZ( real fc[PS2*PS2][3][5] )
{

   real tmp;
   int ID;

   for (int i=0; i<3;   i++)
   for (int k=0; k<PS2; k++)
   for (int j=0; j<PS2; j++)
   {
      ID = k*PS2 + j;

      tmp          = fc[ID][i][3];
      fc[ID][i][3] = fc[ID][i][1];
      fc[ID][i][1] = tmp;
   }

} // FUNCTION : TrasposeFluxXZ



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == WAF )
