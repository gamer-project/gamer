#ifndef __CUFLU_RIEMANNSOLVER_EXACT__
#define __CUFLU_RIEMANNSOLVER_EXACT__



#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  \
       ( RSOLVER == EXACT || CHECK_INTERMEDIATE == EXACT )  &&  \
       ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );

#endif // #ifdef __CUDACC__ ... else ...


// internal functions (GPU_DEVICE is defined in CUFLU.h)
GPU_DEVICE static real Solve_f( const real rho,const real p,const real p_star,const real Gamma );
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
GPU_DEVICE static void Set_Flux( real flux[], const real val[], const real Gamma );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemmanSolver_Exact
// Description :  Exact Riemann solver
//
// Note        :  1. Input data should be primitive variables
//                2. This function is shared by MHM, MHM_RP, and CTU schemes
//                3. Currently it does NOT check the minimum density and pressure criteria
//
// Parameter   :  XYZ         : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out    : Output array to store the average flux along t axis
//                L_In        : Input **primitive** variables in the left region
//                              --> But note that the input passive scalars should be mass density instead of mass fraction
//                R_In        : Input **primitive** variables in the right region
//                              --> But note that the input passive scalars should be mass density instead of mass fraction
//                Gamma       : Ratio of specific heats
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_Exact( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], const real Gamma )
{

   const real Gamma_p1 = Gamma + (real)1.0;
   const real Gamma_m1 = Gamma - (real)1.0;
   const real c        = Gamma_m1 / Gamma_p1;

   real eival[5], L_star[5], R_star[5];
   real L[5], R[5], Temp;


// reorder the input variables for different spatial directions
   for (int v=0; v<5; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );


// solution of the tangential velocity
   L_star[2] = L[2];
   L_star[3] = L[3];
   R_star[2] = R[2];
   R_star[3] = R[3];


// solution of pressure
   {
      const real du = R[1] - L[1];

      real f;
      real f_L;
      real f_R;
      real bound[2];
      real compare[2];

      bound[0] = FMIN( L[4], R[4] );
      bound[1] = FMAX( L[4], R[4] );

      for (int i=0; i<2; i++)
          {
             f_L = Solve_f( L[0], L[4], bound[i], Gamma );
             f_R = Solve_f( R[0], R[4], bound[i], Gamma );

             compare[i] = f_L + f_R + du;
          }

      if( compare[0]*compare[1] > (real)0.0 )
      {
        if( compare[0] > (real)0.0 )
        {
           bound[1] = bound[0];
           bound[0] = (real)0.0;
        }
        else if( compare[1] < (real)0.0 )
        {
           bool Continue;

           bound[0] = bound[1];
           bound[1] = (real)2.0*bound[0];

           do
           {
              for (int i=0; i<2; i++)
              {
                 f_L = Solve_f( L[0], L[4], bound[i], Gamma );
                 f_R = Solve_f( R[0], R[4], bound[i], Gamma );

                 compare[i] = f_L + f_R + du;
              }

              Continue = ( compare[0]*compare[1] > (real)0.0 );

              if ( Continue )
              {
                 bound[0] = bound[1];
                 bound[1] = bound[0]*(real)2.0;
              }
           }
           while ( Continue );
        }
      }

//    search p_star
      do
      {
         L_star[4] = (real)0.5 * ( bound[0] + bound[1] );

         if (  ( L_star[4] == bound[0] )  ||  ( L_star[4] == bound[1] )  )    break;
         else
         {
            f_L = Solve_f( L[0], L[4], L_star[4], Gamma );
            f_R = Solve_f( R[0], R[4], L_star[4], Gamma );
            f   = f_L + f_R + du;

            if ( f > (real)0.0 )    bound[1] = L_star[4];
            else                    bound[0] = L_star[4];
         }
      }
      while ( FABS(f) >= MAX_ERROR );
   }

   R_star[4] = L_star[4];


// solution normal velocity
   {
      real f_L = Solve_f( L[0], L[4], L_star[4], Gamma );
      real f_R = Solve_f( R[0], R[4], R_star[4], Gamma );

      L_star[1] = (real)0.5*( L[1] + R[1] ) + (real)0.5*( f_R - f_L );
   }

   R_star[1] = L_star[1];


// left complete solution
   if ( L_star[4] > L[4] ) // left shock
   {
        real r = L_star[4] / L[4];
        L_star[0] = L[0]*(  ( r + c ) / ( c*r + (real)1.0 )  );   // solution of density
   }

   else                   // left rarefaction
   {
      L_star[0] = L[0]*POW( L_star[4]/L[4], (real)1.0/Gamma );    // solution of density

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(L[4]) )
         printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L[4],      __FILE__, __LINE__, __FUNCTION__ );
      if ( Hydro_CheckNegative(L[0]) )
         printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L[0],      __FILE__, __LINE__, __FUNCTION__ );
      if ( Hydro_CheckNegative(L_star[4]) )
         printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L_star[4], __FILE__, __LINE__, __FUNCTION__ );
      if ( Hydro_CheckNegative(L_star[0]) )
         printf( "ERROR : invalid density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L_star[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real a_L    = SQRT( Gamma * L[4] / L[0] );                  // sound speed of left region
      real a_star = SQRT( Gamma * L_star[4] / L_star[0] );        // sound speed of left star region

      if ( L[1] < a_L  &&  L_star[1] > a_star  ) // sonic rarefaction
      {
         L_star[0] = L[0]*POW(  (real)2.0/Gamma_p1 + c*L[1]/a_L, (real)2.0/Gamma_m1  );
         L_star[1] = (real)2.0*( a_L + (real)0.5*Gamma_m1*L[1] ) / Gamma_p1;
         L_star[4] = L[4]*POW(  (real)2.0/Gamma_p1 + c*L[1]/a_L, (real)2.0*Gamma/Gamma_m1  );
      }
   }

// right complete solution
   if ( R_star[4] > R[4] ) // right shock
   {
      real r    = R_star[4] / R[4];
      R_star[0] = R[0]*(  ( r + c ) / ( c * r + (real)1.0 )  );         // solution of density
   }

   else                    // right rarefaction
   {
      R_star[0] = R[0]*POW( R_star[4]/R[4], (real)1.0/Gamma ); // solution of density

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(R[4]) )
         printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R[4],      __FILE__, __LINE__, __FUNCTION__ );
      if ( Hydro_CheckNegative(R[0]) )
         printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R[0],      __FILE__, __LINE__, __FUNCTION__ );
      if ( Hydro_CheckNegative(R_star[4]) )
         printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R_star[4], __FILE__, __LINE__, __FUNCTION__ );
      if ( Hydro_CheckNegative(R_star[0]) )
         printf( "ERROR : invalid density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R_star[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real a_R    = SQRT( Gamma * R[4] / R[0] ); // sound speed of right region
      real a_star = SQRT( Gamma * R_star[4] / R_star[0] ); // sound speed of right star region

      if ( R[1] > -a_R  &&  R_star[1] < -a_star ) // sonic rarefaction
      {
         R_star[0] = R[0]*POW (  (real)2.0/Gamma_p1 - c*R[1]/a_R, (real)2.0/Gamma_m1  );
         R_star[1] = (real)2.0*( -a_R + (real)0.5*Gamma_m1*R[1] ) / Gamma_p1;
         R_star[4] = R[4]*POW (  (real)2.0/Gamma_p1 - c*R[1]/a_R, (real)2.0*Gamma/Gamma_m1 );
      }
   }

// solve speed of waves
   eival[1] = L_star[1];
   eival[2] = L_star[1];
   eival[3] = L_star[1];

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative( R[4]) )
      printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R[4], __FILE__, __LINE__, __FUNCTION__ );
   if ( Hydro_CheckNegative(R[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R[0], __FILE__, __LINE__, __FUNCTION__ );
   if ( Hydro_CheckNegative(L[4]) )
      printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L[4], __FILE__, __LINE__, __FUNCTION__ );
   if ( Hydro_CheckNegative(L[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif

   if ( L[4] < L_star[4] ) // left shock
   {
      Temp = (real)0.5/Gamma*( Gamma_p1*L_star[4]/L[4] + Gamma_m1 );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(Temp) )
         printf( "ERROR : invalid value (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Temp, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      eival[0] = L[1] - SQRT( Gamma*L[4]/L[0] )*SQRT( Temp );
   }
   else                    // left rarefaction
      eival[0] = L[1] - SQRT ( Gamma*L[4]/L[0] );

   if ( R[4] < R_star[4] ) // right shock
   {
      Temp = (real)0.5/Gamma*( Gamma_p1*R_star[4]/R[4] + Gamma_m1 );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(Temp) )
         printf( "ERROR : invalid value (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Temp, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      eival[4] = R[1] + SQRT( Gamma*R[4]/R[0] )*SQRT( Temp );
   }
   else                    // right rarefaction
      eival[4] = R[1] + SQRT ( Gamma*R[4]/R[0] );


// evaluate the average fluxes along the t axis
   if (  FABS( eival[1] ) < MAX_ERROR  ) // contact wave is zero
   {
      Flux_Out[0] = (real)0.0;
      Flux_Out[1] = L_star[4];
      Flux_Out[2] = (real)0.0;
      Flux_Out[3] = (real)0.0;
      Flux_Out[4] = (real)0.0;
   }

   else
   {
      if ( eival[1] > (real)0.0 )
      {
         if ( eival[0] > (real)0.0 )   Set_Flux( Flux_Out, L,      Gamma );
         else                          Set_Flux( Flux_Out, L_star, Gamma );
      }

      else
      {
         if ( eival[4] < (real)0.0 )   Set_Flux( Flux_Out, R,      Gamma );
         else                          Set_Flux( Flux_Out, R_star, Gamma );
      }
   }


// evaluate the fluxes for passive scalars
// --> note that L_In and R_In are mass density instead of mass fraction for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out[FLUX_DENS] >= (real)0.0 )
   {
      const real vx = Flux_Out[FLUX_DENS] / L_In[DENS];

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L_In[v]*vx;
   }

   else
   {
      const real vx = Flux_Out[FLUX_DENS] / R_In[DENS];

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R_In[v]*vx;
   }
#  endif


// restore the correct order
   Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_RiemannSolve_Exact



//-------------------------------------------------------------------------------------------------------
// Function    :  Solve_f
// Description :  Solve the parameter f in Godunov's method
//
// paremater   :  rho      : Density
//                p        : Pressure
//                p_star   : Pressure in star region
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Solve_f( const real rho, const real p, const real p_star, const real Gamma )
{

   const real Gamma_m1 = Gamma - (real)1.0;
   const real Gamma_p1 = Gamma + (real)1.0;

   real f, Temp;

   if ( p_star > p )
   {
      real A = (real)2.0/( rho*Gamma_p1 );
      real B = p*Gamma_m1/Gamma_p1;
      Temp   = A/(p_star+B);

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(Temp) )
         printf( "ERROR : invalid value (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Temp, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      f = (p_star-p)*SQRT( Temp );
   }

   else
   {
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(p) )
         printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 p, __FILE__, __LINE__, __FUNCTION__ );

      if ( Hydro_CheckNegative(rho) )
         printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 rho, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real a = SQRT( Gamma*p/rho );
      real c = p_star/p;
      f = (real)2.0*a*(  POW( c, (real)0.5*Gamma_m1/Gamma ) - (real)1.0  ) / Gamma_m1;
   }

   return f;

} // FUNCTION : Solve_f



#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Flux
// Description :  Set the flux function evaluated at the given state
//
// Parameter   :  flux  : Output flux
//                val   : Input primitive variables
//                Gamma : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Set_Flux( real flux[], const real val[], const real Gamma )
{

  const real Gamma_m1 = Gamma -(real)1.0;

// set flux
  flux[0] = val[0]*val[1];
  flux[1] = val[0]*val[1]*val[1] + val[4];
  flux[2] = val[0]*val[1]*val[2];
  flux[3] = val[0]*val[1]*val[3];
  flux[4] = val[1]*(  (real)0.5*val[0]*( val[1]*val[1] + val[2]*val[2] + val[3]*val[3] )
                     + val[4]/Gamma_m1 + val[4]  );

} // FUNCTION : Set_Flux
#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )



#endif // #if ( MODEL == HYDRO  &&  ( RSOLVER == EXACT || CHECK_INTE == EXACT ) && ( SCHEME == MHM/MHM_RP/CTU ) )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_EXACT__
