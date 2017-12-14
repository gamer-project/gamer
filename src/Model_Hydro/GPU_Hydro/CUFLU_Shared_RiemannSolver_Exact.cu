#ifndef __CUFLU_RIEMANNSOLVER_EXACT_CU__
#define __CUFLU_RIEMANNSOLVER_EXACT_CU__



#include "CUFLU.h"
#include "CUFLU_Shared_FluUtility.cu"

static __device__ FluVar CUFLU_RiemannSolver_Exact( const int XYZ, FluVar5 &eival_out, FluVar5 &L_star_out,
                                                    FluVar5 &R_star_out, const FluVar L_In, const FluVar R_In,
                                                    const real Gamma );
static __device__ real Solve_f( const real _rho, const real p, const real p_star, const real Gamma,
                                const real Gamma_m1, const real Gamma_p1, const real _Gamma,
                                const real _Gamma_m1, const real _Gamma_p1 );
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
static __device__ void Set_Flux( FluVar &flux, const FluVar val, const real _Gamma_m1 );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_RiemmanSolver_Exact
// Description :  Exact Riemann solver
//
// Note        :  1. The input data should be primitive variables
//                2. This function is shared by WAF, MHM, MHM_RP, and CTU schemes
//                3. The "__noinline__" qualifier is added in Fermi GPUs for higher performance and
//                   faster compilation (not work in WAF scheme!?)
//
// Parameter   :  XYZ         : Target spatial direction : (0/1/2) --> (x/y/z)
//                eival_out   : Output array to store the speed of waves
//                L_star_out  : Output array to store the primitive variables in the left star region
//                R_star_out  : Output array to store the primitive variables in the right star region
//                L_In        : Input **primitive** variables in the left region
//                              --> But note that the input passive scalars should be mass density instead of mass fraction
//                R_In        : Input **primitive** variables in the right region
//                              --> But note that the input passive scalars should be mass density instead of mass fraction
//                Gamma       : Ratio of specific heats
//------------------------------------------------------------------------------------------------------
#if ( __CUDA_ARCH__ >= 200  &&  FLU_SCHEME != WAF )
__noinline__
#endif
__device__ FluVar CUFLU_RiemannSolver_Exact( const int XYZ, FluVar5 *eival_out, FluVar5 *L_star_out,
                                             FluVar5 *R_star_out, const FluVar L_In,
                                             const FluVar R_In, const real Gamma )
{

// reorder the input variables for different spatial directions
   FluVar L = CUFLU_Rotate3D( L_In, XYZ, true );
   FluVar R = CUFLU_Rotate3D( R_In, XYZ, true );

   const real  Gamma_p1 = Gamma + (real)1.0;
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma    = (real)1.0/Gamma;
   const real _Gamma_p1 = (real)1.0/Gamma_p1;
   const real _Gamma_m1 = (real)1.0/Gamma_m1;
   const real c         = Gamma_m1*_Gamma_p1;
   const real _RhoL     = (real)1.0/L.Rho;
   const real _RhoR     = (real)1.0/R.Rho;

   FluVar eival, L_star, R_star;    // don't use FluVar5 since these variables will be sent into Set_Flux()
   real   Temp;


// solution of the tangential velocity
   L_star.Py = L.Py;
   L_star.Pz = L.Pz;
   R_star.Py = R.Py;
   R_star.Pz = R.Pz;


// solution of pressure
   {
      const real du = R.Px - L.Px;

      real f;
      real f_L;
      real f_R;
      real bound[2];
      real compare[2];

      bound[0] = FMIN( L.Egy, R.Egy );
      bound[1] = FMAX( L.Egy, R.Egy );

      for (int i=0; i<2; i++)
      {
         f_L = Solve_f( _RhoL, L.Egy, bound[i], Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );
         f_R = Solve_f( _RhoR, R.Egy, bound[i], Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );

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
                 f_L = Solve_f( _RhoL, L.Egy, bound[i], Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );
                 f_R = Solve_f( _RhoR, R.Egy, bound[i], Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );

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
         L_star.Egy = (real)0.5*( bound[0] + bound[1] );

         if (  ( L_star.Egy == bound[0] )  ||  ( L_star.Egy == bound[1] )  )     break;
         else
         {
            f_L = Solve_f( _RhoL, L.Egy, L_star.Egy, Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );
            f_R = Solve_f( _RhoR, R.Egy, L_star.Egy, Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );

            f = f_L + f_R + du;

            if ( f > (real)0.0 )    bound[1] = L_star.Egy;
            else                    bound[0] = L_star.Egy;
         }
      }
      while ( FABS(f) >= MAX_ERROR );

   }

   R_star.Egy = L_star.Egy;


// solution normal velocity
   {
      real f_L = Solve_f( _RhoL, L.Egy, L_star.Egy, Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );
      real f_R = Solve_f( _RhoR, R.Egy, R_star.Egy, Gamma, Gamma_m1, Gamma_p1, _Gamma, _Gamma_m1, _Gamma_p1 );
      L_star.Px = (real)0.5*( L.Px + R.Px ) + (real)0.5*( f_R - f_L );
   }

   R_star.Px = L_star.Px;


// left complete solution
   if ( L_star.Egy > L.Egy ) // left shock
   {
      real r = L_star.Egy / L.Egy;
      L_star.Rho = L.Rho*(  ( r + c ) / ( c*r + (real)1.0 )  ); // solution of density
   }

   else                      // left rarefaction
   {
      L_star.Rho = L.Rho*POW( L_star.Egy/L.Egy, _Gamma ); //solution of density

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(L.Egy) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L.Egy,      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(L.Rho) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L.Rho,      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(L_star.Egy) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L_star.Egy, __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(L_star.Rho) )
         printf( "ERROR : negative density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 L_star.Rho, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real a_L    = SQRT( Gamma*L.Egy     /L.Rho );      // sound speed of left region
      real a_star = SQRT( Gamma*L_star.Egy/L_star.Rho ); // sound speed of left star region

      if ( L.Px < a_L  &&  L_star.Px > a_star ) // sonic rarefaction
      {
         L_star.Rho = L.Rho*POW (  (real)2.0*_Gamma_p1 + c*L.Px/a_L, (real)2.0*_Gamma_m1  );
         L_star.Px  = (real)2.0*( a_L + (real)0.5*Gamma_m1*L.Px ) * _Gamma_p1;
         L_star.Egy = L.Egy*POW (  (real)2.0*_Gamma_p1 + c*L.Px/a_L, (real)2.0*_Gamma_m1*Gamma  );
      }
   }


// right complete solution
   if ( R_star.Egy > R.Egy ) // right shock
   {
      real r = R_star.Egy / R.Egy;
      R_star.Rho = R.Rho*(  ( r + c ) / ( c*r + (real)1.0 )  ); // solution of density
   }

   else                     // right rarefaction
   {
      R_star.Rho = R.Rho * POW( R_star.Egy/R.Egy, _Gamma ); // solution of density

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(R.Egy) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R.Egy,      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R.Rho) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R.Rho,      __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R_star.Egy) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R_star.Egy, __FILE__, __LINE__, __FUNCTION__ );
      if ( CUFLU_CheckNegative(R_star.Rho) )
         printf( "ERROR : negative density(%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 R_star.Rho, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real a_R    = SQRT( Gamma*R.Egy     /R.Rho );      // sound speed of right region
      real a_star = SQRT( Gamma*R_star.Egy/R_star.Rho ); // sound speed of right star region

      if (  R.Px > -a_R  &&  R_star.Px < -a_star  )      // sonic rarefaction
      {
         R_star.Rho = R.Rho*POW(  (real)2.0*_Gamma_p1 - c*R.Px/a_R, (real)2.0*_Gamma_m1  );
         R_star.Px  = (real)2.0*( -a_R + (real)0.5*Gamma_m1*R.Px ) * _Gamma_p1;
         R_star.Egy = R.Egy*POW(  (real)2.0*_Gamma_p1 - c*R.Px/a_R, (real)2.0*_Gamma_m1*Gamma  );
      }
   }

// solve speed of waves
   eival.Px = L_star.Px;
   eival.Py = L_star.Px;
   eival.Pz = L_star.Px;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(L.Egy) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L.Egy,      __FILE__, __LINE__, __FUNCTION__ );
   if ( CUFLU_CheckNegative(L.Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L.Rho,      __FILE__, __LINE__, __FUNCTION__ );
   if ( CUFLU_CheckNegative(R.Egy) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R.Egy,      __FILE__, __LINE__, __FUNCTION__ );
   if ( CUFLU_CheckNegative(R.Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R.Rho,      __FILE__, __LINE__, __FUNCTION__ );
#  endif

   if ( L.Egy < L_star.Egy )
   {
      Temp = (real)0.5*_Gamma*( Gamma_p1*L_star.Egy/L.Egy + Gamma_m1 );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(Temp) )
         printf( "ERROR : negative value (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Temp, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      eival.Rho = L.Px - SQRT( Gamma*L.Egy/L.Rho )*SQRT( Temp );        // left shock
   }
   else
      eival.Rho = L.Px - SQRT( Gamma*L.Egy/L.Rho );                     // left rarefaction

   if ( R.Egy < R_star.Egy )
   {
      Temp = (real)0.5*_Gamma*( Gamma_p1*R_star.Egy/R.Egy + Gamma_m1 );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(Temp) )
         printf( "ERROR : negative value (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Temp, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      eival.Egy = R.Px + SQRT( Gamma*R.Egy/R.Rho )*SQRT( Temp );        // right shock
   }
   else
      eival.Egy = R.Px + SQRT( Gamma*R.Egy/R.Rho );                     // right rarefaction


// evaluate the average flux along the t axis
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

   FluVar Flux_Out;

   if (  FABS( eival.Px ) < MAX_ERROR  ) // contact wave is zero
   {
      Flux_Out.Rho = (real)0.0;
      Flux_Out.Px  = L_star.Egy;
      Flux_Out.Py  = (real)0.0;
      Flux_Out.Pz  = (real)0.0;
      Flux_Out.Egy = (real)0.0;
   }

   else
   {
      if ( eival.Px > (real)0.0 )
      {
         if ( eival.Rho > (real)0.0 )  Set_Flux( Flux_Out, L,      _Gamma_m1 );
         else                          Set_Flux( Flux_Out, L_star, _Gamma_m1 );
      }

      else
      {
         if ( eival.Egy < (real)0.0 )  Set_Flux( Flux_Out, R,      _Gamma_m1 );
         else                          Set_Flux( Flux_Out, R_star, _Gamma_m1 );
      }
   }


// evaluate the fluxes for passive scalars
// --> note that L_In and R_In are mass density instead of mass fraction for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out.Rho >= (real)0.0 )
   {
      const real vx = Flux_Out.Rho / L_In.Rho;

      for (int v=0; v<NCOMP_PASSIVE; v++)    Flux_Out.Passive[v] = L_In.Passive[v]*vx;
   }

   else
   {
      const real vx = Flux_Out.Rho / R_In.Rho;

      for (int v=0; v<NCOMP_PASSIVE; v++)    Flux_Out.Passive[v] = R_In.Passive[v]*vx;
   }
#  endif


// restore the correct order
   Flux_Out = CUFLU_Rotate3D( Flux_Out, XYZ, false );

   return Flux_Out;

#  elif ( FLU_SCHEME == WAF )

   eival_out->Rho  = eival.Rho;
   eival_out->Px   = eival.Px;
   eival_out->Py   = eival.Py;
   eival_out->Pz   = eival.Pz;
   eival_out->Egy  = eival.Egy;

   L_star_out->Rho = L_star.Rho;
   L_star_out->Px  = L_star.Px;
   L_star_out->Py  = L_star.Py;
   L_star_out->Pz  = L_star.Pz;
   L_star_out->Egy = L_star.Egy;

   R_star_out->Rho = R_star.Rho;
   R_star_out->Px  = R_star.Px;
   R_star_out->Py  = R_star.Py;
   R_star_out->Pz  = R_star.Pz;
   R_star_out->Egy = R_star.Egy;

// return an arbitray FluVar type variable since it's useless
   return L_In;

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU ) ... else ...

} // FUNCTION : CUFLU_RiemannSolver_Exact



//-------------------------------------------------------------------------------------------------------
// Function    :  Solve_f
// Description :  Solve the parameter f in Godunov's method
//
// paremater   :  rho         : Density
//                p           : Pressure
//                p_star      : Pressure in star region
//                Gamma       : Ratio of specific heats
//                Gamma_m1    : Gamma - 1
//                Gamma_p1    : Gamma + 1
//                _Gamma      : 1/Gamma
//                _Gamma_m1   : 1/(Gamma-1)
//                _Gamma_p1   : 1/(Gamma+1)
//-------------------------------------------------------------------------------------------------------
__device__ real Solve_f( const real _rho, const real p, const real p_star, const real Gamma, const real Gamma_m1,
                         const real Gamma_p1, const real _Gamma, const real _Gamma_m1, const real _Gamma_p1 )
{

   real f, Temp;

   if ( p_star > p)
   {
      real A = (real)2.0*_rho*_Gamma_p1;
      real B = p*Gamma_m1*_Gamma_p1;
      Temp   = A/(p_star+B);

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(Temp) )
         printf( "ERROR : negative value (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Temp, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      f = (p_star-p)*SQRT( Temp );
   }

   else
   {
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(p) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 p, __FILE__, __LINE__, __FUNCTION__ );

      if ( CUFLU_CheckNegative((real)1.0/_rho) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 (real)1.0/_rho, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      real a = SQRT( Gamma*p*_rho );
      real c = p_star/p;
      f = (real)2.0*a*(  POW( c, (real)0.5*Gamma_m1*_Gamma ) - (real)1.0  ) * _Gamma_m1;
   }

   return f;

} // FUNCTION : Solve_f



#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Flux
// Description :  Set the flux function evaluated at the given stat
//
// Parameter   :  flux        : Output flux
//                val         : Input primitive variables
//                _Gamma_m1   : 1/(Gamma-1)
//-------------------------------------------------------------------------------------------------------
__device__ void Set_Flux( FluVar &flux, const FluVar val, const real _Gamma_m1 )
{

// set flux
  flux.Rho = val.Rho*val.Px;
  flux.Px  = val.Rho*val.Px*val.Px + val.Egy;
  flux.Py  = val.Rho*val.Px*val.Py;
  flux.Pz  = val.Rho*val.Px*val.Pz;
  flux.Egy = val.Px*(  (real)0.5*val.Rho*( val.Px*val.Px + val.Py*val.Py + val.Pz*val.Pz )
                      + val.Egy*_Gamma_m1 + val.Egy  );

} // FUNCTION : Set_Flux
#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_EXACT_CU__
