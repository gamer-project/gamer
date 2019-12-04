#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__

#include "CUFLU.h"
#include <stdio.h>

struct Fun_params
{
  real M_Dsqr;
  real E_D;
};

GPU_DEVICE 
real VectorDotProduct( real V1, real V2, real V3 );
GPU_DEVICE 
static void Fun_DFun (real Q, void *ptr, real * f, real * df, real Gamma);
GPU_DEVICE 
static void NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma);

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_Rotate3D
// Description :  Rotate the input fluid variables properly to simplify the 3D calculation
//
// Note        :  1. x : (0,1,2,3,4) <--> (0,1,2,3,4)
//                   y : (0,1,2,3,4) <--> (0,2,3,1,4)
//                   z : (0,1,2,3,4) <--> (0,3,1,2,4)
//                2. Work if InOut includes/excludes passive scalars since they are not modified at all
//
// Parameter   :  InOut    : Array storing both the input and output data
//                XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Forward  : (true/false) <--> (forward/backward)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_Rotate3D (real InOut[], const int XYZ, const bool Forward)
{
  if (XYZ == 0)
    return;

  real Temp[3];
  for (int v = 0; v < 3; v++)
    Temp[v] = InOut[v + 1];

  if (Forward)
    {
      switch (XYZ)
	{
	case 1:
	  InOut[1] = Temp[1]; InOut[2] = Temp[2]; InOut[3] = Temp[0]; break;
	case 2:
	  InOut[1] = Temp[2]; InOut[2] = Temp[0]; InOut[3] = Temp[1]; break;
	}
    }

  else				// backward
    {
      switch (XYZ)
	{
	case 1:
	  InOut[1] = Temp[2]; InOut[2] = Temp[0]; InOut[3] = Temp[1]; break;
	case 2:
	  InOut[1] = Temp[1]; InOut[2] = Temp[2]; InOut[3] = Temp[0]; break;
	}
    }
}				// FUNCTION : SRHydro_Rotate3D


//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_GetTemperature
// Description :  Evaluate the fluid temperature
//
// Note        :  1. Currently work with the adiabatic EOS and general EOS
//                2. For simplicity, currently this function only returns kB*T/mc**2, which does NOT include normalization
//                   --> Also note that currently it only checks minimum temperature but not minimum density or pressure
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density
//                Gamma        : the adiabatic index
//
// Return      :  kB*T/mc**2
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_GetTemperature (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                             const real Gamma, const real MinTemp )
{
  real In[5] = {Dens, MomX, MomY, MomZ, Engy};

  real guess, root;
  real Msqr = VectorDotProduct(In[1], In[2], In[3]);
  real M = SQRT (Msqr); // magnitude of momentum
  real Dsqr = SQR(In[0]);
  real abc = (real)1.0 / Dsqr;
  real E_D = In[4] / In[0];
  real E_Dsqr = abc * SQR(In[4]);
  real M_Dsqr = abc * Msqr;

# if ( EOS == APPROXIMATED_GENERAL )
# if   ( CONSERVED_ENERGY == 1 )
/* initial guess  */
  real Constant = E_Dsqr - M_Dsqr;
  if ( Constant > (real)1.0 )
  {
   if ( Dsqr - (real)0.0625 * Msqr >= (real)0 )
   {
     if ( Constant > (real)2.5 ) 
          guess = SQRT( FMA( (real)1.0, E_Dsqr, - FMA ( (real)0.93750, M_Dsqr, (real)2.0 ) ) ) / (real)3.0;
     else guess = (Constant - (real)1.0) / (real)3.0;
   }
   else // 1 - (M/D)**2 < 0
   {
     if ( Constant >  (real)1.5 + SQRT( FMA ((real)0.0625, M_Dsqr, (real) -0.75 ))) 
          guess = SQRT( FMA( (real)1.0, E_Dsqr, - FMA ( (real)0.93750, M_Dsqr, (real)2.0 ) ) ) / (real)3.0;
     else guess = (Constant - (real)1.0) / (real)3.0;
   }
  } else return MinTemp;
# elif ( CONSERVED_ENERGY == 2 )
/* initial guess  */
  real Constant = E_Dsqr - M_Dsqr + (real)2.0 * E_D;
  if ( Constant > (real)0.0 )
  {
      if ( Dsqr - (real)0.0625 * Msqr >= (real)0 )
      {
        if ( Constant > (real)1.5 ) 
             guess = SQRT( E_Dsqr + (real)2.0 * E_D - (real)0.93750 * M_Dsqr - (real)1.0 ) / (real)3.0;
        else guess = Constant / (real)3.0;
      }
      else // 1 - (M/D)**2 < 0
      {
        if ( Constant >  (real)0.5 + SQRT( (real)0.0625 * M_Dsqr - (real)0.75 )) 
             guess = SQRT( E_Dsqr + (real)2.0 * E_D - (real)0.93750 * M_Dsqr - (real)1.0 ) / (real)3.0;
        else guess = Constant / (real)3.0;
      }
  } else return MinTemp;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
# elif ( EOS == CONSTANT_GAMMA )
  real Gamma_m1 = Gamma - (real) 1.0;
/* initial guess */
# if   ( CONSERVED_ENERGY == 1 )
  real Constant = E_Dsqr - M_Dsqr;

  if ( Constant > (real)1.0 )
  {
     if ( Constant > (real)1.0 + (real)2* (Gamma_m1 / Gamma ) * (M / In[0]) )
     {
       real A = (real)1.0 / SQR(Gamma_m1);
       real B = (real)2.0 /Gamma_m1;
       real C = (((real)2*Gamma-(real)1.0)/(Gamma*Gamma)) * M_Dsqr - E_Dsqr;
       real delta = SQRT( B * B - (real)4 * A * C );
       guess = -(real)2.0 *  C / ( B + delta);
     }
     else guess = (real)0.5*Gamma_m1 * ( Constant - (real)1.0);
  } else return MinTemp;
# elif ( CONSERVED_ENERGY == 2 )
  real Constant = E_Dsqr - M_Dsqr + (real)2 * E_D;

  if ( Constant > (real)0.0 )
  {
    if ( Constant > (real)2* (Gamma_m1 / Gamma ) * (M / In[0]) )
    {
      real A = (real)1.0 / SQR(Gamma_m1);
      real B = (real)2.0 /Gamma_m1;
      real C = (((real)2*Gamma-(real)1.0)/(Gamma*Gamma)) * M_Dsqr - E_Dsqr - (real)2 * E_D;
      real delta = SQRT( B * B - (real)4 * A * C );
      guess = -(real)2.0 *  C / ( B + delta);
    }
    else guess = (real)0.5*Gamma_m1 * Constant;
  } else return MinTemp;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
# else
# error: unsupported EoS!
# endif

     struct Fun_params params = { M_Dsqr, E_D };

# ifdef CHECK_FAILED_CELL_IN_FLUID
  if ( guess != guess || guess >= HUGE_NUMBER || guess <= TINY_NUMBER )
  { 
    printf ("guess root = %14.7e\n", guess);
    printf ("D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n", In[0], In[1], In[2], In[3], In[4]);
  }
# endif 

  NewtonRaphsonSolver(&params ,&root, guess, (real) TINY_NUMBER, (real) EPSILON, Gamma);
  //root = Constant / (real)3.0;
  return root;
}				// FUNCTION : SRHydro_GetTemperature

//-------------------------------------------------------------------------------------------------------
// Function    :  SpecificEnthalpy
// Description :  Evaluate the specific enthalpy
//
// Note        :  1. Currently work with the adiabatic EOS and general EOS
//                2. For simplicity, currently this function only returns h/c**2, which does NOT include normalization
//
// Parameter   :  Con[]       : conservative variables
//
// Return      :  h/c**2
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SpecificEnthalpy( const real Con[], real Temp, real Gamma )
{
   real h;

#  if ( EOS == APPROXIMATED_GENERAL )
   h = FMA( (real)2.5, Temp, SQRT( FMA( (real)2.25, Temp*Temp, (real)1.0 ) ) );

// if there is an overfolw due to ultra-high temperature, we simply set Gamma = 4/3.
   if ( h != h && Con != NULL) 
   {    
		Gamma = (real)1.3333333;
        real Gamma_m1 = Gamma - (real)1.0;

#       if ( CONSERVED_ENERGY == 1  )
        real E_Dsqr = SQR( Con[ENGY] / Con[DENS] );
#       else
        real E_Dsqr = SQR( ( Con[ENGY] + Con[DENS] ) / Con[DENS] );
#       endif

        real M_Dsqr = VectorDotProduct(Con[MOMX], Con[MOMY], Con[MOMZ]) / SQR(Con[DENS]);

        real A = (real)1.0 / SQR(Gamma_m1);
        real B = (real)2.0 /Gamma_m1;
        real C = (((real)2.0*Gamma-(real)1.0)/(Gamma*Gamma)) * M_Dsqr - E_Dsqr;

        real delta = SQRT( B * B - (real)4.0 * A * C );

        Temp = -(real)2.0 *  C / ( B + delta);

        h = (real)4.0 * Temp;
   }
   else if ( h!= h && Con == NULL )
   {
        printf("Conservative variables should be provided! function: %s: %d\n", __FUNCTION__, __LINE__);
   }

#  elif ( EOS == CONSTANT_GAMMA ) 
   h = (real)1.0 + Temp * Gamma / (Gamma-(real)1.0);
#  else
#  error: unsupported EoS!
#  endif

   return h;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_Con2Pri
// Description :  Convert the conserved variables to the primitive variables
//
// Note        :  1. This function always check if the pressure to be returned is greater than the
//                   given minimum threshold
//
// Parameter   :  In                 : Array storing the input conserved variables
//                Out                : Array to store the output primitive variables
//                Gamma              : Gamma
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_Con2Pri (const real In[], real Out[], const real Gamma, const real MinTemp)
{
   real Temp = SRHydro_GetTemperature (In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY], Gamma, MinTemp );


   real h = SpecificEnthalpy(In, Temp, Gamma);

   real factor = In[0]*h;
   Out[1] = In[1]/factor;
   Out[2] = In[2]/factor;
   Out[3] = In[3]/factor;

   real Lorentz = SQRT((real)1.0 + VectorDotProduct(Out[1], Out[2], Out[3]));

   Out[0] = In[0]/Lorentz;
   Out[4] = Out[0] * Temp; // P = rho*kB*T*c**2

   return Lorentz;
}// FUNCTION : SRHydro_Con2Pri


#ifdef __CUDACC__
GPU_DEVICE
void SRHydro_Pri2Con (const real In[], real Out[], const real Gamma)
#else
template <class T> 
void SRHydro_Pri2Con (const T In[], T Out[], const T Gamma)
#endif
{

  real h = SpecificEnthalpy( NULL, In[4]/In[0], Gamma );
  real nh   = In[0] * h;

  real Factor0 = (real)1.0 + VectorDotProduct(In[1], In[2], In[3]);
  real Factor1 = SQRT(Factor0); // Lorentz factor
  
  Out[0] = In[0] * Factor1; // mass density in inertial frame
  Out[1] = Out[0]*h*In[1]; // MomX
  Out[2] = Out[0]*h*In[2]; // MomX
  Out[3] = Out[0]*h*In[3]; // MomX
# if   ( CONSERVED_ENERGY == 1 )
  Out[4] = FMA( nh, Factor0, - In[4] ); // total_energy
# elif ( CONSERVED_ENERGY == 2 )
  Out[4] = nh * Factor0 - In[4] - Out[0]; // ( total_energy ) - ( rest_mass_energy )
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
}				// FUNCTION : SRHydro_Pri2Con

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_4Velto3Vel
// Description :  Convert 4-velocity to 3-velocity
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
GPU_DEVICE
void SRHydro_4Velto3Vel ( const real In[], real Out[])
#else
template <class T> 
void SRHydro_4Velto3Vel ( const T In[], T Out[])
#endif
{
  real Factor = (real)1.0 / SQRT ((real)1.0 + VectorDotProduct(In[1], In[2], In[3]));

  Out[0] = In[0];
  Out[1] = In[1] * Factor;
  Out[2] = In[2] * Factor;
  Out[3] = In[3] * Factor;
  Out[4] = In[4];
}				// FUNCTION : SRHydro_4Velto3Vel

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_3Velto4Vel
// Description :  Convert 3-velocity to 4-velocity
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
GPU_DEVICE
void SRHydro_3Velto4Vel (const real In[], real Out[])
#else
template <class T> 
void SRHydro_3Velto4Vel (const T In[], T Out[])
#endif

{
  real Factor = (real)1.0 / SQRT ((real)1.0 - VectorDotProduct(In[1], In[2], In[3]));

  Out[0] = In[0];
  Out[1] = In[1] * Factor;
  Out[2] = In[2] * Factor;
  Out[3] = In[3] * Factor;
  Out[4] = In[4];
}				// FUNCTION : SRHydro_4Velto3Vel

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_Con2Flux
// Description :  Evaluate the hydrodynamic fluxes by the input conserved variables
//
// Parameter   :  XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux     : Array to store the output fluxes
//                Input    : Array storing the input conserved variables
//                Gamma    : adiabatic index
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_Con2Flux (const int XYZ, real Flux[], const real Con[], const real Pri[], const real Gamma, const real MinTemp )
{
//*** we don't need to include passive scalars since they don't have to be rotated ***
  real ConVar[NCOMP_FLUID], PriVar[NCOMP_FLUID];

  for (int v=0;v<NCOMP_FLUID;v++) ConVar[v] = Con[v];
  for (int v=0;v<NCOMP_FLUID;v++) PriVar[v] = Pri[v];

  SRHydro_Rotate3D (ConVar, XYZ, true);
  SRHydro_Rotate3D (PriVar, XYZ, true);

  real Lorentz = SQRT((real)1.0 + VectorDotProduct(PriVar[1], PriVar[2], PriVar[3]));

  real Vx = PriVar[1] / Lorentz;

  Flux[0] = ConVar[0] * Vx;
  Flux[1] = FMA( ConVar[1], Vx, PriVar[4] );
  Flux[2] = ConVar[2] * Vx;
  Flux[3] = ConVar[3] * Vx;
# if ( CONSERVED_ENERGY == 1 )
  Flux[4] = ConVar[1];
# elif ( CONSERVED_ENERGY == 2 )
  Flux[4] = ConVar[1] - Flux[0];
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif

  SRHydro_Rotate3D (Flux, XYZ, false);
}				// FUNCTION : SRHydro_Con2Flux

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_CheckMinDens
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_CheckMinDens (const real InDens, const real MinDens)
{
  return FMAX (InDens, MinDens);
}// FUNCTION : SRHydro_CheckMinDens


//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_CheckMinTemp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_CheckMinTemp (const real InTemp, const real MinTemp)
{
  return FMAX (InTemp, MinTemp);
}// FUNCTION : SRHydro_CheckMinDens


//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_CheckMinTempInEngy
// Description :  Ensure that the Temp in the input total energy is greater than the given threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) temperature caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
//                   --> Currently it simply sets a minimum allowed value for temerature
//                       --> Please set MIN_TEMP in the runtime parameter file "Input__Parameter"
//                3. One must input conserved variables instead of primitive variables
//
// Parameter   :  Cons[]      : D, Mx, My, Mz, E
//                MinTemp     : Minimum allowed temperature
//                Gamma       : adiabatic index
//
// Return      :  Total energy with pressure greater than the given threshold
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_CheckMinTempInEngy (const real Cons[], const real MinTemp, const real Gamma )
{

// conservative variables should NOT be put into SpecificEnthalpy() now,
// --> In case h is overflow, Gamma will be automatically set as 4/3.
  real h_min = SpecificEnthalpy( NULL, MinTemp, Gamma );

  real D  = Cons[0];
  real Mx = Cons[1];
  real My = Cons[2];
  real Mz = Cons[3];

  real Msqr = SQR(Mx) + SQR(My) + SQR(Mz);
  real Dh = D*h_min;
  real factor = SQRT(Dh*Dh + Msqr);

# if ( CONSERVED_ENERGY == 1 )
  real E_min = factor - D*Dh*MinTemp / factor;
# elif ( CONSERVED_ENERGY == 2 )
  real E_min = factor - D*Dh*MinTemp / factor - D;
# endif

  if ( Cons[4] >= E_min) return Cons[4];
  else                   return E_min;
}



//-------------------------------------------------------------------------------------------------------
// Function    : SRHydro_CheckUnphysical
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool SRHydro_CheckUnphysical( const real Con[], const real Pri[], const real Gamma, const real MinTemp, const char s[], const int line, bool show )
{
   real discriminant;
   real Msqr;
   real ConsVar[NCOMP_FLUID];
   real Pri4Vel[NCOMP_FLUID];
   real Pri3Vel[NCOMP_FLUID];


   for (int i = 0;i < NCOMP_FLUID; i++) { ConsVar[i]=NAN; Pri4Vel[i]=NAN; Pri3Vel[i]=NAN; }

//--------------------------------------------------------------//
//------------ only check conserved variables-------------------//
//--------------------------------------------------------------//

    if ( Con != NULL && Pri == NULL){
      for(int i=0; i< NCOMP_FLUID; i++) ConsVar[i]=Con[i];


// check NaN
      if (  ConsVar[DENS] != ConsVar[DENS]
         || ConsVar[MOMX] != ConsVar[MOMX]
         || ConsVar[MOMY] != ConsVar[MOMY]
         || ConsVar[MOMZ] != ConsVar[MOMZ]
         || ConsVar[ENGY] != ConsVar[ENGY]  )                                             goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= ConsVar[DENS] || ConsVar[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= ConsVar[MOMX] || ConsVar[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= ConsVar[MOMY] || ConsVar[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= ConsVar[MOMZ] || ConsVar[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= ConsVar[ENGY] || ConsVar[ENGY]  >= (real)HUGE_NUMBER )             goto FAIL;


// check energy
      Msqr = VectorDotProduct( ConsVar[MOMX], ConsVar[MOMY], ConsVar[MOMZ] );
#     if ( CONSERVED_ENERGY == 1 )
      discriminant = ( SQR( ConsVar[ENGY] ) -  Msqr ) / SQR ( ConsVar[DENS] );
      if ( discriminant <= (real)1.0 )                                                        goto FAIL;
#     elif ( CONSERVED_ENERGY == 2 )
      discriminant = SQR(ConsVar[ENGY]/ConsVar[DENS]) + (real)2*(ConsVar[ENGY]/ConsVar[DENS]) - Msqr/SQR(ConsVar[DENS]);
      if ( discriminant <= TINY_NUMBER )                                                   goto FAIL;
#     else
#     error: CONSERVED_ENERGY must be 1 or 2!
#     endif

//      SRHydro_Con2Pri(ConsVar, Pri4Vel, Gamma, MinTemp);
//
//// check NaN
//      if (  Pri4Vel[DENS] != Pri4Vel[DENS]
//         || Pri4Vel[MOMX] != Pri4Vel[MOMX]
//         || Pri4Vel[MOMY] != Pri4Vel[MOMY]
//         || Pri4Vel[MOMZ] != Pri4Vel[MOMZ]
//         || Pri4Vel[ENGY] != Pri4Vel[ENGY]  )                                              goto FAIL;
//
//// check +inf and -inf
//      if (  (real)  TINY_NUMBER >= Pri4Vel[DENS] || Pri4Vel[DENS]  >= (real)HUGE_NUMBER
//         || (real) -HUGE_NUMBER >= Pri4Vel[MOMX] || Pri4Vel[MOMX]  >= (real)HUGE_NUMBER
//         || (real) -HUGE_NUMBER >= Pri4Vel[MOMY] || Pri4Vel[MOMY]  >= (real)HUGE_NUMBER
//         || (real) -HUGE_NUMBER >= Pri4Vel[MOMZ] || Pri4Vel[MOMZ]  >= (real)HUGE_NUMBER
//         || (real)  TINY_NUMBER >= Pri4Vel[ENGY] || Pri4Vel[ENGY]  >= (real)HUGE_NUMBER )              goto FAIL;
//
//// check whether 3-velocity is greater or equal to speed of light
//      SRHydro_4Velto3Vel(Pri4Vel,Pri3Vel);
//
//      if (VectorDotProduct( Pri3Vel[1], Pri3Vel[2], Pri3Vel[3] ) >= (real) 1.0)                      goto FAIL;

// pass all checks 
      return false;
   }

//--------------------------------------------------------------//
//------------ only check primitive variables-------------------//
//--------------------------------------------------------------//

   else if ( Con == NULL && Pri != NULL){

      for(int i=0; i< NCOMP_FLUID; i++) Pri4Vel[i]=Pri[i];


// check NaN
      if (  Pri4Vel[DENS] != Pri4Vel[DENS]
         || Pri4Vel[MOMX] != Pri4Vel[MOMX]
         || Pri4Vel[MOMY] != Pri4Vel[MOMY]
         || Pri4Vel[MOMZ] != Pri4Vel[MOMZ]
         || Pri4Vel[ENGY] != Pri4Vel[ENGY]  )                                       goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= Pri4Vel[DENS] || Pri4Vel[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri4Vel[MOMX] || Pri4Vel[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri4Vel[MOMY] || Pri4Vel[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri4Vel[MOMZ] || Pri4Vel[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Pri4Vel[ENGY] || Pri4Vel[ENGY]  >= (real)HUGE_NUMBER )       goto FAIL;


      SRHydro_4Velto3Vel(Pri4Vel,Pri3Vel);
      SRHydro_Pri2Con(Pri4Vel, ConsVar, (real) Gamma);

// check whether 3-velocity is greater or equal to speed of light
      if (VectorDotProduct( Pri3Vel[1], Pri3Vel[2], Pri3Vel[3] ) >= (real) 1.0)                      goto FAIL;
   
// check NaN
      if (  ConsVar[DENS] != ConsVar[DENS]
         || ConsVar[MOMX] != ConsVar[MOMX]
         || ConsVar[MOMY] != ConsVar[MOMY]
         || ConsVar[MOMZ] != ConsVar[MOMZ]
         || ConsVar[ENGY] != ConsVar[ENGY]  )                                       goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= ConsVar[DENS] || ConsVar[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= ConsVar[MOMX] || ConsVar[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= ConsVar[MOMY] || ConsVar[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= ConsVar[MOMZ] || ConsVar[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= ConsVar[ENGY] || ConsVar[ENGY]  >= (real)HUGE_NUMBER )        goto FAIL;


// check energy
      Msqr = VectorDotProduct( ConsVar[MOMX], ConsVar[MOMY], ConsVar[MOMZ] );
#     if ( CONSERVED_ENERGY == 1 )
      discriminant = ( SQR( ConsVar[ENGY] ) - Msqr ) / SQR ( ConsVar[DENS] );
      if ( discriminant <= (real)1.0 )                                                   goto FAIL;
#     elif ( CONSERVED_ENERGY == 2 )
      discriminant = SQR(ConsVar[ENGY]/ConsVar[DENS]) + (real)2*(ConsVar[ENGY]/ConsVar[DENS]) - Msqr/SQR(ConsVar[DENS]);
      if ( discriminant <= TINY_NUMBER )                                                   goto FAIL;
#     else
#     error: CONSERVED_ENERGY must be 1 or 2!
#     endif      

// pass all checks 
      return false;
   }

// print all variables if goto FAIL
      FAIL:
      {
        if ( show ) 
         {
           printf( "\n\nfunction: %s: %d\n", s, line);
 
           if ( Con != NULL && Pri == NULL)
           {
              printf( "D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n",
                                   ConsVar[DENS], ConsVar[MOMX], ConsVar[MOMY], ConsVar[MOMZ], ConsVar[ENGY]);
#             if ( CONSERVED_ENERGY == 1 )
              printf( "(E^2-|M|^2)/D^2=%14.7e\n", discriminant );
#             elif ( CONSERVED_ENERGY == 2 )
              printf( "E^2+2*E*D-|M|^2=%14.7e\n", discriminant );
#             endif
           }
           else
           {
              printf( "n=%14.7e, Ux=%14.7e, Uy=%14.7e, Uz=%14.7e, P=%14.7e\n", 
                                   Pri4Vel[0], Pri4Vel[1], Pri4Vel[2], Pri4Vel[3], Pri4Vel[4]);
              printf( "Vx=%14.7e, Vy=%14.7e, Vz=%14.7e, |V|=%14.7e\n",
                                   Pri3Vel[1], Pri3Vel[2], Pri3Vel[3], SQRT( VectorDotProduct( Pri3Vel[1], Pri3Vel[2], Pri3Vel[3] )));
           }
 
         }
        return true;
      }
}

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_GetPressure
// Description :  Evaluate the fluid pressure
//
// Note        :  1. Currently only work with the adiabatic EOS
//                2. Invoked by the functions "Hydro_GetTimeStep_Fluid", "Prepare_PatchData", "InterpolateGhostZone",
//                   "Hydro_Aux_Check_Negative" ...
//                3. One must input conserved variables instead of primitive variables
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density
//                Gamma      : Adiabatic index
//
// Return      :  Pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_GetPressure (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy, 
                          const real Gamma, const real MinTemp )
{
  real In[NCOMP_FLUID];
  real Out[NCOMP_FLUID];

  In[0] = Dens;
  In[1] = MomX;
  In[2] = MomY;
  In[3] = MomZ;
  In[4] = Engy;

  SRHydro_Con2Pri (In, Out, Gamma, MinTemp);

  return Out[4];
}				// FUNCTION : SRHydro_GetPressure


GPU_DEVICE
real SoundSpeedSquare( real Temp, real Gamma )
{
  real Cs_sq;

# if ( EOS == APPROXIMATED_GENERAL )
  real factor = SQRT( (real)2.25*Temp*Temp + (real)1.0 );

  real A = (real) 4.5*SQR(Temp) + (real) 5.0 * Temp * factor;
  real B = (real)18.0*SQR(Temp) + (real)12.0 * Temp * factor + (real)3.0;

  Cs_sq = A/B;

# elif ( EOS == CONSTANT_GAMMA )
  Cs_sq = Gamma * Temp; /* Mignone Eq 4 */

# endif
  
  return Cs_sq;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  Evaluate internal energy density (including rest mass energy)
//                true : lab frame
//                false: fluid frame
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_InternalEngy( real Con[], real Pri[], real Lorentz, real Gamma, bool frame)
{
  real h, Eint;

  h = SpecificEnthalpy( Con, Pri[4]/Pri[0], Gamma );

  frame ? ( Eint = Con[0] * h - Lorentz * Pri[4] )
         :( Eint = Pri[0] * h -           Pri[4] );

  return Eint;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  Evaluate thermal energy density
//                true : lab frame
//                false: fluid frame
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_ThermalEngy( real Con[], real Pri[], real Lorentz, real Gamma, bool frame )
{
  real h, E_thermal;

  h = SpecificEnthalpy( Con, Pri[4]/Pri[0], Gamma );

  
  frame ? ( E_thermal = Con[0] * (h-(real)1.0) - Lorentz * Pri[4] )
         :( E_thermal = Pri[0] * (h-(real)1.0) -           Pri[4] );

  return E_thermal;
}

//-------------------------------------------------------------------------------------------------------
// Function    :
// Description : Evaluate kinetic energy density in lab frame 
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_KineticEngy( real Con[], real Pri[], real Lorentz, real Gamma )
{
  real h;

  h = SpecificEnthalpy( Con, Pri[4]/Pri[0], Gamma );
  
  return ( Con[DENS] * h + Pri[4] ) * ( Lorentz - (real)1.0 );
}

#ifdef GRAVITY
GPU_DEVICE
real SRHydro_PoissonSource( real Con[], real Gamma, real MinTemp )
{
  real Source, Pri[NCOMP_FLUID];
 
  SRHydro_Con2Pri( Con, Pri, Gamma, MinTemp );

  Source = Pri[0];
  
  return Source;
}
#endif


GPU_DEVICE
static void 
NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma)
{
 int iter = 0;
 int max_iter = 20;

 real f, df;
 real delta;
 real tolerance;
 *root = guess;

 do
   {
     iter++;
     Fun_DFun(*root, ptr, &f, &df, Gamma);

#    ifdef CHECK_FAILED_CELL_IN_FLUID
     if ( df == (real)0.0 )                                                  printf("derivative is zero\n");
     if (  f != f  ||(real) -HUGE_NUMBER >= f  || f  >= (real)HUGE_NUMBER )  printf("function value is not finite\n");
     if ( df != df ||(real) -HUGE_NUMBER >= df || df >= (real)HUGE_NUMBER )  printf("derivative value is not finite\n");
#    endif     

      delta = f/df;
      *root = *root - delta;

      tolerance =  FMA( epsrel, FABS(*root), epsabs );

   }while ( fabs(delta) >= tolerance && iter < max_iter );

}


//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
//-------------------------------------------------------------------------------------------------------

GPU_DEVICE
static void
Fun_DFun (real Temp, void *ptr, real * f, real * df, real Gamma)
{
  struct Fun_params *params = (struct Fun_params *) ptr;

  real M_Dsqr = (params->M_Dsqr);
  real E_D    = (params->E_D);
  real Tsqr = Temp * Temp;

# if ( EOS == APPROXIMATED_GENERAL )
  real abc = SQRT(FMA( (real)9.0, Tsqr, (real)4.0 ));
  real h = FMA( (real)2.5, Temp, SQRT(FMA( (real)2.25, Tsqr, (real)1.0 ))); // approximate enthalpy
  real dh = FMA( (real)9.0, Temp / SQRT(FMA( (real)36.0, Tsqr, (real)16.0 )), (real)2.5 );
  real hsqr = SQR(h);
# if   (CONSERVED_ENERGY == 1)
  real Constant = FMA( E_D, E_D, - M_Dsqr );
  *f = (real)3.5 * Tsqr + (real)1.5 * Temp * abc + hsqr * Tsqr / (hsqr + M_Dsqr) + (real)1.0 - Constant;
# elif (CONSERVED_ENERGY == 2)
  real Constant = SQR(E_D) + (real)2*(E_D) - M_Dsqr;
  *f = (real)3.5 * Tsqr + (real)1.5 * Temp * abc + hsqr * Tsqr / (hsqr + M_Dsqr) - Constant;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
  *df = FMA( (real)7.0, Temp, FMA( (real)1.5, abc, (real)13.5 * Tsqr / abc )) + (real)2*h*Temp*((h*hsqr + M_Dsqr*h + Temp*dh*M_Dsqr) / SQR( hsqr + M_Dsqr) );

# elif ( EOS == CONSTANT_GAMMA )
  real zeta = (real)1.0 / ( Gamma - (real)1.0 );
  real alpha = Gamma * zeta;
  real h = (real)1 + alpha * Temp;
  real hsqr = SQR(h);
  real beta = ((real)2 - Gamma) * zeta * alpha;
  real theta = (real)2 * zeta;

# if   (CONSERVED_ENERGY == 1)
  real Constant = SQR(E_D) - M_Dsqr;
  *f = (real)1.0 + beta * Tsqr + theta * Temp + hsqr * Tsqr / (hsqr + M_Dsqr) - Constant; 
# elif (CONSERVED_ENERGY == 2)
  real Constant = SQR(E_D) + (real)2*(E_D) - M_Dsqr;
  *f = beta * Tsqr + theta * Temp + hsqr * Tsqr / (hsqr + M_Dsqr) - Constant; 
# endif
  real dh = alpha;
  *df = (real)2 * beta * Temp + theta + (real)2*h*Temp*((h*hsqr + M_Dsqr*h + Temp*dh*M_Dsqr) / SQR( hsqr + M_Dsqr) );
# else
# error: unsupported EoS!
# endif // #if ( EOS == APPROXIMATED_GENERAL )
}



GPU_DEVICE
real VectorDotProduct( real V1, real V2, real V3 )
{
  real Product = (real)0.0;
  
  Product = FMA( V1, V1 , Product );
  Product = FMA( V2, V2 , Product );
  Product = FMA( V3, V3 , Product );
  
  return Product;
} 

#ifndef __CUDACC__
template void SRHydro_Pri2Con(const double*, double*, const double);
template void SRHydro_Pri2Con(const float* , float* , const float);
template void SRHydro_4Velto3Vel ( const double* , double* );
template void SRHydro_4Velto3Vel ( const float* , float* );
template void SRHydro_3Velto4Vel ( const double* , double* );
template void SRHydro_3Velto4Vel ( const float* , float* );
#endif
 


#endif
