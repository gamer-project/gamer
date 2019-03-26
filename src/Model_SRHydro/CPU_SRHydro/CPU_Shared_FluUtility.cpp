#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__

#include "CUFLU.h"
#include <stdio.h>

struct Fun_params
{
  real M_Dsqr;
  real E_D;
};


#ifdef __CUDACC__
GPU_DEVICE
real SRHydro_GetTemperature (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                             const real Gamma, const real MinTemp );
GPU_DEVICE 
static real Fun (real Q, void *ptr);   // function to be solved
GPU_DEVICE 
static real DFun (real Q, void *ptr);  // the first derivative of above function
GPU_DEVICE 
static void Fun_DFun (real Q, void *ptr, real * f, real * df, real Gamma);
GPU_DEVICE 
void SRHydro_4Velto3Vel (const real In[], real Out[]);
GPU_DEVICE 
static void NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma);
GPU_DEVICE 
real VectorDotProduct( const real *V1, const real *V2, int Ini_i, int Final_i );
#else 
#include "../../../include/SRHydroPrototypes.h"
static real Fun (real Q, void *ptr);   // function to be solved
static real DFun (real Q, void *ptr);  // the first derivative of above function
static void Fun_DFun (real Q, void *ptr, real * f, real * df, real Gamma);
static void NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma);
#endif


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
void SRHydro_Con2Pri (const real In[], real Out[], const real Gamma, const real MinTemp)
{
      real Temp = SRHydro_GetTemperature (In[0], In[1], In[2], In[3], In[4], Gamma, MinTemp );
#if ( EOS == RELATIVISTIC_IDEAL_GAS )
      real h = FMA( 2.5, Temp, SQRT( FMA( 2.25, Temp*Temp, 1.0 ) ) );
#elif ( EOS == IDEAL_GAS ) 
      real h = 1 + Temp * Gamma / (Gamma-1.0);
#else
#error: unsupported EoS!
#endif
      real factor = In[0]*h;
      Out[1] = In[1]/factor;
      Out[2] = In[2]/factor;
      Out[3] = In[3]/factor;

      real factor1 = SQRT(1 + VectorDotProduct(Out, Out, 1, 3));

      Out[0] = In[0]/factor1;
      Out[4] = Out[0] * Temp; // P = nkT
}// FUNCTION : SRHydro_Con2Pri


GPU_DEVICE
void SRHydro_Pri2Con (const real In[], real Out[], const real Gamma)
{
# if ( EOS == RELATIVISTIC_IDEAL_GAS )
  real nh = FMA( 2.5, In[4], SQRT( FMA( 2.25, SQR(In[4]), SQR(In[0]) ) )); // approximate enthalpy * proper number density
# elif ( EOS == IDEAL_GAS )
  real Gamma_m1 = (real) Gamma - 1.0;
  real nh = In[0] + ( Gamma / Gamma_m1) * In[4]; // enthalpy * proper number density
# else
# error: unsupported EoS!
# endif

  real Factor0 = 1.0 + VectorDotProduct(In, In, 1, 3);
  real Factor1 = SQRT(Factor0); // Lorentz factor
  real Factor2 = nh * Factor1;
  
  Out[0] = In[0] * Factor1; // number density in inertial frame
  Out[1] = Factor2 * In[1]; // MomX
  Out[2] = Factor2 * In[2]; // MomX
  Out[3] = Factor2 * In[3]; // MomX
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
GPU_DEVICE
void SRHydro_4Velto3Vel ( const real In[], real Out[])
{
  real Factor = 1 / SQRT (1 + VectorDotProduct(In, In, 1, 3));

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
GPU_DEVICE
void SRHydro_3Velto4Vel (const real In[], real Out[])
{
  real Factor = 1 / SQRT (1 - VectorDotProduct(In, In, 1, 3));

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
void SRHydro_Con2Flux (const int XYZ, real Flux[], const real Input[], const real Gamma, const real MinTemp )
{
  real ConVar[NCOMP_FLUID];	// don't need to include passive scalars since they don't have to be rotated1
  real PriVar4[NCOMP_FLUID];	// D, Ux, Uy, Uz, P
  real PriVar3[NCOMP_FLUID];	// D, Vx, Vy, Vz, P
  real Pres, Vx;

  for (int v = 0; v < NCOMP_FLUID; v++)   ConVar[v] = Input[v];

  SRHydro_Rotate3D (ConVar, XYZ, true);

  SRHydro_Con2Pri (ConVar, PriVar4, Gamma, MinTemp);

  SRHydro_4Velto3Vel (PriVar4, PriVar3);

  Vx = PriVar3[1];
  Pres = PriVar3[4];

  Flux[0] = ConVar[0] * Vx;
  Flux[1] = FMA( ConVar[1], Vx, Pres );
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
# if ( EOS == IDEAL_GAS ) 
  real h_min = 1.0 + Gamma * MinTemp / (Gamma - 1.0);
# elif ( EOS == RELATIVISTIC_IDEAL_GAS )
  real h_min = 2.5*MinTemp + SQRT(2.25*MinTemp*MinTemp + 1);
# endif

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
      if (    TINY_NUMBER >= ConsVar[DENS] || ConsVar[DENS]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= ConsVar[MOMX] || ConsVar[MOMX]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= ConsVar[MOMY] || ConsVar[MOMY]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= ConsVar[MOMZ] || ConsVar[MOMZ]  >= HUGE_NUMBER
         ||   TINY_NUMBER >= ConsVar[ENGY] || ConsVar[ENGY]  >= HUGE_NUMBER )             goto FAIL;


// check energy
      Msqr = VectorDotProduct( ConsVar, ConsVar, MOMX, MOMZ );
#     if ( CONSERVED_ENERGY == 1 )
      discriminant = FMA(ConsVar[ENGY], ConsVar[ENGY], - Msqr - SQR(ConsVar[DENS]));
      if ( discriminant <= TINY_NUMBER )                                                   goto FAIL;
#     elif ( CONSERVED_ENERGY == 2 )
      discriminant = SQR(ConsVar[ENGY]) + 2*ConsVar[ENGY]*ConsVar[DENS] - Msqr;
      if ( discriminant <= TINY_NUMBER )                                                   goto FAIL;
#     else
#     error: CONSERVED_ENERGY must be 1 or 2!
#     endif

      SRHydro_Con2Pri(ConsVar, Pri4Vel, Gamma, MinTemp);

// check NaN
      if (  Pri4Vel[DENS] != Pri4Vel[DENS]
         || Pri4Vel[MOMX] != Pri4Vel[MOMX]
         || Pri4Vel[MOMY] != Pri4Vel[MOMY]
         || Pri4Vel[MOMZ] != Pri4Vel[MOMZ]
         || Pri4Vel[ENGY] != Pri4Vel[ENGY]  )                                              goto FAIL;

// check +inf and -inf
      if (    TINY_NUMBER >= Pri4Vel[DENS] || Pri4Vel[DENS]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= Pri4Vel[MOMX] || Pri4Vel[MOMX]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= Pri4Vel[MOMY] || Pri4Vel[MOMY]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= Pri4Vel[MOMZ] || Pri4Vel[MOMZ]  >= HUGE_NUMBER
         ||   TINY_NUMBER >= Pri4Vel[ENGY] || Pri4Vel[ENGY]  >= HUGE_NUMBER )              goto FAIL;

// check whether 3-velocity is greater or equal to speed of light
      SRHydro_4Velto3Vel(Pri4Vel,Pri3Vel);

      if (VectorDotProduct( Pri3Vel, Pri3Vel, 1, 3 ) >= 1.0)                      goto FAIL;

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
      if (    TINY_NUMBER >= Pri4Vel[DENS] || Pri4Vel[DENS]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= Pri4Vel[MOMX] || Pri4Vel[MOMX]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= Pri4Vel[MOMY] || Pri4Vel[MOMY]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= Pri4Vel[MOMZ] || Pri4Vel[MOMZ]  >= HUGE_NUMBER
         ||   TINY_NUMBER >= Pri4Vel[ENGY] || Pri4Vel[ENGY]  >= HUGE_NUMBER )       goto FAIL;


      SRHydro_4Velto3Vel(Pri4Vel,Pri3Vel);
      SRHydro_Pri2Con(Pri4Vel, ConsVar, (real) Gamma);

// check whether 3-velocity is greater or equal to speed of light
      if (VectorDotProduct( Pri3Vel, Pri3Vel, 1, 3 ) >= 1.0)                      goto FAIL;
   
// check NaN
      if (  ConsVar[DENS] != ConsVar[DENS]
         || ConsVar[MOMX] != ConsVar[MOMX]
         || ConsVar[MOMY] != ConsVar[MOMY]
         || ConsVar[MOMZ] != ConsVar[MOMZ]
         || ConsVar[ENGY] != ConsVar[ENGY]  )                                       goto FAIL;

// check +inf and -inf
      if (    TINY_NUMBER >= ConsVar[DENS] || ConsVar[DENS]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= ConsVar[MOMX] || ConsVar[MOMX]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= ConsVar[MOMY] || ConsVar[MOMY]  >= HUGE_NUMBER
         ||  -HUGE_NUMBER >= ConsVar[MOMZ] || ConsVar[MOMZ]  >= HUGE_NUMBER
         ||   TINY_NUMBER >= ConsVar[ENGY] || ConsVar[ENGY]  >= HUGE_NUMBER )        goto FAIL;


// check energy
      Msqr = VectorDotProduct( ConsVar, ConsVar, MOMX, MOMZ );
#     if ( CONSERVED_ENERGY == 1 )
      discriminant = FMA(ConsVar[ENGY], ConsVar[ENGY], - Msqr - SQR(ConsVar[DENS]));
      if ( discriminant <= TINY_NUMBER )                                              goto FAIL;
#     elif ( CONSERVED_ENERGY == 2 )
      discriminant = SQR(ConsVar[ENGY]) + 2*ConsVar[ENGY]*ConsVar[DENS] - Msqr;
      if ( discriminant <= TINY_NUMBER )                                              goto FAIL;
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
 
            printf( "D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n",
                                 ConsVar[DENS], ConsVar[MOMX], ConsVar[MOMY], ConsVar[MOMZ], ConsVar[ENGY]);
#           if ( CONSERVED_ENERGY == 1 )
            printf( "E^2-|M|^2-D^2=%14.7e\n", discriminant );
#           elif ( CONSERVED_ENERGY == 2 )
            printf( "E^2+2*E*D-|M|^2=%14.7e\n", discriminant );
#           endif
            printf( "n=%14.7e, Ux=%14.7e, Uy=%14.7e, Uz=%14.7e, P=%14.7e\n", 
                                 Pri4Vel[0], Pri4Vel[1], Pri4Vel[2], Pri4Vel[3], Pri4Vel[4]);
            printf( "Vx=%14.7e, Vy=%14.7e, Vz=%14.7e, |V|=%14.7e\n",
                                 Pri3Vel[1], Pri3Vel[2], Pri3Vel[3], SQRT( VectorDotProduct( Pri3Vel, Pri3Vel, 1, 3 )));
 
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



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_GetTemperature
// Description :  Evaluate the fluid temperature
//
// Note        :  1. Currently only work with the adiabatic EOS
//                2. For simplicity, currently this function only returns **pressure/density**, which does NOT include normalization
//                   --> For OPT__FLAG_LOHNER_TEMP only
//                   --> Also note that currently it only checks minimum pressure but not minimum density
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density
//                Gamma        : the adiabatic index
//
// Return      :  Temperature
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_GetTemperature (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                             const real Gamma, const real MinTemp )
{
  real In[5] = {Dens, MomX, MomY, MomZ, Engy};

  real guess, root;
  real Msqr = VectorDotProduct(In, In, 1, 3);
  real M = SQRT (Msqr); // magnitude of momentum
  real Dsqr = SQR(In[0]);
  real abc = 1.0 / Dsqr;
  real E_D = In[4] / In[0];
  real E_Dsqr = abc * SQR(In[4]);
  real M_Dsqr = abc * Msqr;

# if ( EOS == RELATIVISTIC_IDEAL_GAS )
# if   ( CONSERVED_ENERGY == 1 )
/* initial guess  */
  real Constant = E_Dsqr - M_Dsqr;
  if ( Constant > 1.0 )
   {
	if ( Dsqr - 0.0625 * Msqr >= 0 )
	  {
	    if ( Constant > 2.5 ) 
	      {
		 guess = SQRT(FMA( 0.1111111, E_Dsqr, - FMA ( 0.1041667, M_Dsqr, 0.2222222 ) ));
	      }
	    else guess = (Constant - 1.0) * 0.3333333;
	  }
	else // 1 - (M/D)**2 < 0
	  {
	    if ( Constant >  1.5 + SQRT( FMA (0.0625, M_Dsqr, - 0.75 ))) 
	      {
		 guess = SQRT(FMA( 0.1111111, E_Dsqr, - FMA ( 0.1041667, M_Dsqr, 0.2222222 ) ));
	      }
	    else guess = (Constant - 1.0) * 0.3333333;
	  }
    } else return MinTemp;
# elif ( CONSERVED_ENERGY == 2 )
/* initial guess  */
   real Constant = E_Dsqr - M_Dsqr + 2 * E_D;
   if ( Constant > 0.0 )
    {
	   if ( Dsqr - 0.0625 * Msqr >= 0 )
	     {
	       if ( Constant > 1.5 ) 
		 {
		    guess = SQRT( 0.1111111 * E_Dsqr + 0.2222222 * E_D - 0.1041667 * M_Dsqr - 0.1111111 );
		 }
	       else guess = Constant * 0.3333333;
	     }
	   else // 1 - (M/D)**2 < 0
	     {
	       if ( Constant >  0.5 + SQRT( 0.0625 * M_Dsqr - 0.75 )) 
		 {
		    guess = SQRT( 0.1111111 * E_Dsqr + 0.2222222 * E_D - 0.1041667 * M_Dsqr - 0.1111111 );
		 }
	       else guess = Constant * 0.3333333;
	     }
    } else return MinTemp;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
# elif ( EOS == IDEAL_GAS )
  real Gamma_m1 = Gamma - (real) 1.0;
/* initial guess */
# if   ( CONSERVED_ENERGY == 1 )
   real Constant = E_Dsqr - M_Dsqr;
   if ( Constant > 1.0 )
    {
	      if ( Constant > 1.0 + 2* (Gamma_m1 / Gamma ) * (M / In[0]) )
		{
		  real A = 1.0 / SQR(Gamma_m1);
		  real B = 2.0 /Gamma_m1;
		  real C = ((2*Gamma-1.0)/(Gamma*Gamma)) * M_Dsqr - E_Dsqr;
		  real delta = SQRT( B * B - 4 * A * C );
		  guess = -2.0 *  C / ( B + delta);
		}
	     else guess = 0.5*Gamma_m1 * ( Constant - 1.0);
    } else return MinTemp;
# elif ( CONSERVED_ENERGY == 2 )
    real Constant = E_Dsqr - M_Dsqr + 2 * E_D;
    if ( Constant > 0.0 )
     {
	  if ( Constant > 2* (Gamma_m1 / Gamma ) * (M / In[0]) )
	    {
	      real A = 1.0 / SQR(Gamma_m1);
	      real B = 2.0 /Gamma_m1;
	      real C = ((2*Gamma-1.0)/(Gamma*Gamma)) * M_Dsqr - E_Dsqr - 2 * E_D;
	      real delta = SQRT( B * B - 4 * A * C );
	      guess = -2.0 *  C / ( B + delta);
	    }
	  else guess = 0.5*Gamma_m1 * Constant;
     } else return MinTemp;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
#else
#error: unsupported EoS!
#endif

   struct Fun_params params = { M_Dsqr, E_D };

#ifdef CHECK_NEGATIVE_IN_FLUID
   if ( guess != guess || guess >= HUGE_NUMBER || guess <= TINY_NUMBER )
    { 
      printf ("guess root = %14.7e\n", guess);
      printf ("D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n", In[0], In[1], In[2], In[3], In[4]);
    }
#endif 

   NewtonRaphsonSolver(&params ,&root, guess, TINY_NUMBER, EPSILON, Gamma);

   return root;
}				// FUNCTION : SRHydro_GetTemperature




static void 
NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma)
{
 int iter = 0;
 int max_iter = 20;

 real f, df;
 real root_old;
 real tolerance;
 *root = guess;
 do
   {
     iter++;
     Fun_DFun(*root, ptr, &f, &df, Gamma);

#    ifdef CHECK_NEGATIVE_IN_FLUID
     if ( df == 0.0 )                                            printf("derivative is zero\n");
     if (  f != f  || -HUGE_NUMBER >= f  || f  >= HUGE_NUMBER )  printf("function value is not finite\n");
     if ( df != df || -HUGE_NUMBER >= df || df >= HUGE_NUMBER )  printf("derivative value is not finite\n");
#    endif     

      root_old = *root;
      *root = FMA( -f, 1.0/df, *root );
      //printf("df = %20.17e\n", df);
      //printf("f  = %20.17e\n", f);
      tolerance =  FMA( epsrel, FABS(*root), epsabs );
   }while ( fabs(root_old - *root) >= tolerance && iter < max_iter );

}


//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
//-------------------------------------------------------------------------------------------------------

static void
Fun_DFun (real Temp, void *ptr, real * f, real * df, real Gamma)
{
  struct Fun_params *params = (struct Fun_params *) ptr;

  real M_Dsqr = (params->M_Dsqr);
  real E_D    = (params->E_D);
  real Tsqr = Temp * Temp;

# if ( EOS == RELATIVISTIC_IDEAL_GAS )
  real abc = SQRT(FMA( 9, Tsqr, 4 ));
  real h = FMA( 2.5, Temp, SQRT(FMA( 2.25, Tsqr, 1.0 ))); // approximate enthalpy
  real dh = FMA( 9.0, Temp / SQRT(FMA( 36, Tsqr, 16 )), 2.5 );
  real hsqr = SQR(h);
# if   (CONSERVED_ENERGY == 1)
  real Constant = FMA( E_D, E_D, - M_Dsqr );
  *f = FMA( 3.5, Tsqr, FMA( 1.5, FMA( Temp, abc, hsqr * Tsqr / (hsqr + M_Dsqr) ),  1.0 - Constant ));
# elif (CONSERVED_ENERGY == 2)
  real Constant = SQR(E_D) + 2*(E_D) - M_Dsqr;
  *f = 3.5 * Tsqr + 1.5 * Temp * abc + hsqr * Tsqr / (hsqr + M_Dsqr) - Constant;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
  *df = FMA( 7, Temp, FMA( 1.5, abc, 13.5 * Tsqr / abc )) + 2*h*Temp*((h*hsqr + M_Dsqr*h + Temp*dh*M_Dsqr) / SQR( hsqr + M_Dsqr) );

# elif ( EOS == IDEAL_GAS )
  real zeta = 1.0 / ( Gamma - 1.0 );
  real alpha = Gamma * zeta;
  real h = 1 + alpha * Temp;
  real hsqr = SQR(h);
  real beta = (2 - Gamma) * zeta * alpha;
  real theta = 2 * zeta;

# if   (CONSERVED_ENERGY == 1)
  real Constant = SQR(E_D) - M_Dsqr;
  *f = 1.0 + beta * Tsqr + theta * Temp + hsqr * Tsqr / (hsqr + M_Dsqr) - Constant; 
# elif (CONSERVED_ENERGY == 2)
  real Constant = SQR(E_D) + 2*(E_D) - M_Dsqr;
  *f = beta * Tsqr + theta * Temp + hsqr * Tsqr / (hsqr + M_Dsqr) - Constant; 
# endif
  real dh = alpha;
  *df = 2 * beta * Temp + theta + 2*h*Temp*((h*hsqr + M_Dsqr*h + Temp*dh*M_Dsqr) / SQR( hsqr + M_Dsqr) );
# else
# error: unsupported EoS!
# endif // #if ( EOS == RELATIVISTIC_IDEAL_GAS )
}

void QuadraticSolver (real A, real B, real C, real delta, real *x_plus, real *x_minus, const int line)
{
  real tolerance = EPSILON;

  if ( FABS(A) > tolerance )
  {
       if ( delta >= 0.0 )
       {
             real factor = FMA( -0.5, B, SIGN(B) * -0.5 * SQRT(delta) );
     
           if ( B >= 0.0 )
           {
     	     *x_plus   = C/factor;
     	     *x_minus  = factor/A;      return;
           }
           else
           {
     	     *x_plus   = factor/A;
     	     *x_minus  = C/factor;      return;
           }
       }
       else if ( -tolerance < delta )
       {
         *x_plus = -0.5*B/A;
         *x_minus = *x_plus;            return;
       }
       else                             goto NO_REAL_SOLUTIONS_CASE1;
  }
  else
  {
      if ( FABS(B) >= tolerance )
      {
        *x_minus = -C/B;                return;
      }
      else                              goto NO_REAL_SOLUTIONS_CASE2;
  }


     NO_REAL_SOLUTIONS_CASE1:
     {
#       ifdef CHECK_NEGATIVE_IN_FLUID
        printf( "case1 line: %d No real solution in Quadratic Solver!\n", line);
        printf( "A=%14.7e, B=%14.7e, C=%14.7e\n", A, B, C);
        printf( "B*B-4*A*C=%14.7e\n", B*B-4*A*C);
#       endif
        return;
     }

     NO_REAL_SOLUTIONS_CASE2:
     {
#       ifdef CHECK_NEGATIVE_IN_FLUID
        printf( "case2 line: %d No real solution in Quadratic Solver!\n", line);
        printf( "A=%14.7e, B=%14.7e, C=%14.7e\n", A, B, C);
        printf( "B*B-4*A*C=%14.7e\n", B*B-4*A*C);
#       endif
        return;
     }
}


real VectorDotProduct( const real *V1, const real *V2, int Ini_i, int Final_i )
{

real Product = 0.0;

for (int i=Ini_i; i<=Final_i; i++)
 {
    Product = FMA( V1[i], V2[i] , Product );
 }

return Product;
} 

#endif
