#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__

#include "CUFLU.h"
#include <stdio.h>

#if ( CONSERVED_ENERGY == 1 )
# error: CONSERVED_ENERGY == 1 does not support!!
#endif

GPU_DEVICE 
real VectorDotProduct( real V1, real V2, real V3 );
GPU_DEVICE 
void SRHydro_HTilde_Function (real HTilde, real M_Dsqr, real Constant, real *Fun, real *DiffFun, real Gamma);
GPU_DEVICE 
static
void NewtonRaphsonSolver(void (*FunPtr)(real, real, real, real*, real*, real), real M_D, real E_D, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma);

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

GPU_DEVICE
void SRHydro_HTilde2Temperature (const real HTilde, real *Temp, real *DiffTemp, const real Gamma )
{

#  if ( EOS == APPROXIMATED_GENERAL )
   real Factor0 = (real)2.0*SQR(HTilde) + (real)4.0*HTilde;
   real Factor1 = SQRT( (real)9.0*SQR(HTilde) + (real)18.0*HTilde + (real)25.0 );
   real Factor2 = (real)5.0*HTilde + (real)5.0 + Factor1;


   if ( Temp != NULL )
   {
     *Temp = Factor0 / Factor2;
   }

   if ( DiffTemp != NULL )
   {
     *DiffTemp = ( (real)4.0*HTilde + (real)4.0 ) / Factor2
			   -  Factor0 * ( ( (real)9.0*HTilde + (real)9.0 ) /  Factor1 + (real)5.0 ) / SQR(Factor2);
   }
#  elif ( EOS == CONSTANT_GAMMA )
   if ( Temp != NULL )
   {
     *Temp = ( Gamma - (real)1.0 ) * HTilde / Gamma;
   }

   if ( DiffTemp != NULL )
   {
     *DiffTemp = ( Gamma - (real)1.0 ) / Gamma;
   }
#  endif

}// FUNCTION : SRHydro_GetTemperature

GPU_DEVICE
real SRHydro_GetHTilde( const real Con[], real Gamma )
{
  real HTilde, guess, Discrimination, Constant;

  real Msqr = VectorDotProduct(Con[1], Con[2], Con[3]);
  real Dsqr = SQR(Con[0]);
  real abc = (real)1.0 / Dsqr;
  real E_D = Con[4] / Con[0];
  real M_Dsqr = abc * Msqr;
  real M_D = SQRT( M_Dsqr );


  // (x+y)(x-y) is more accurate than x**2-y**2
  Constant = (E_D + M_D) * (E_D - M_D) + (real)2.0 * E_D;

  
# if ( EOS == APPROXIMATED_GENERAL )
  real A = (real)437.0 * M_Dsqr + (real)117.0;
  real B = (real)1.0 + M_Dsqr;
  real C = (real)43.0*M_Dsqr + (real)63.0;

  Discrimination  = (real)3240000.0 * SQR( B );
  Discrimination /= SQR( A );


  if ( Constant >= Discrimination )
  {
     guess = (real)1.3333333 * SQRT( Constant );
  }
  else
  {
     guess  = (real)11.18 * SQRT( B*C*Constant + (real)45.0*B*B );
	 guess -= (real)75.0 * B;
	 guess /= C;
  }


# elif ( EOS == CONSTANT_GAMMA )
  real A = Gamma / ( Gamma - (real)1.0 );
  real B = -A * ( A - (real)1.0 ) * ( M_Dsqr + (real)1.0 );
  real Asqr = A*A;
  real C = Asqr * M_Dsqr + Asqr - (real)2.0*A*M_Dsqr - (real)2.0*A + (real)1.0;
  

  Discrimination = B / C;

  if ( Constant >= Discrimination )
  {
     guess  = Constant * Gamma;
  }
  else
  {
     real Term1 = (real)1.0 - (real)2.0/A + (real)1.0 / ( Asqr * ( M_Dsqr + (real)1.0 ) );
     real Term2 = (real)2.0 / Gamma;
     real Term3 = -Constant;

     real Del = SQRT( Term2*Term2 - (real)4.0*Term1*Term3 );

     real guess1 = (real)2.0 * Term3 / ( -Term2 + Del );
     real guess2 = (real)2.0 * Term3 / ( -Term2 - Del );

     guess1 > (real)0.0 ?  guess = guess1 : guess = guess2 ;
  }


# endif
 
  void (*FunPtr)( real HTilde, real M_Dsqr, real Constant, real *Fun, real *DiffFun, real Gamma ) = &SRHydro_HTilde_Function;

  NewtonRaphsonSolver(FunPtr, M_Dsqr, -E_D, &HTilde, guess, (real) TINY_NUMBER, (real) EPSILON, Gamma);
  
  return HTilde;
}



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
real SRHydro_GetTemperature( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                             const real Gamma, const real MinTemp )
{
   real Temp, Cons[NCOMP_FLUID], HTilde;
   
   Cons[DENS] = Dens;
   Cons[MOMX] = MomX;
   Cons[MOMY] = MomY;
   Cons[MOMZ] = MomZ;
   Cons[ENGY] = Engy;

   HTilde = SRHydro_GetHTilde( Cons, Gamma );

   SRHydro_HTilde2Temperature( HTilde, &Temp, NULL, Gamma ); 

   return Temp;
}



GPU_DEVICE
real SRHydro_Temperature2HTilde (const real Temperature, const real Gamma )
{
    real HTilde;

#   if ( EOS == APPROXIMATED_GENERAL )
	HTilde = (real)2.5*Temperature + 2.25*SQR(Temperature)/(SQRT( (real)2.25*SQR(Temperature) + (real)1.0 ) + (real)1.0);
#   elif ( EOS == CONSTANT_GAMMA )
    HTilde = Gamma * Temperature / (Gamma - (real)1.0);
#   endif

	return HTilde;
}

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
   real h = 0.0;
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
   real LorentzFactor, HTilde, Factor0, Usqr;

   HTilde = SRHydro_GetHTilde( In, Gamma );

   Factor0 = In[DENS] * HTilde + In[DENS];

   Out[1] = In[MOMX] / Factor0;
   Out[2] = In[MOMY] / Factor0;
   Out[3] = In[MOMZ] / Factor0;

   Usqr = VectorDotProduct( Out[1], Out[2], Out[3] );

   LorentzFactor = SQRT( (real)1.0 + Usqr );

#  ifdef USE_3_VELOCITY
// 4-velocities -> 3-velocities
   Out[1] /= LorentzFactor;
   Out[2] /= LorentzFactor;
   Out[3] /= LorentzFactor;

   real Vsqr = VectorDotProduct( Out[1], Out[2], Out[3] );

   LorentzFactor = (real)1.0 / SQRT( (real)1.0 - Vsqr );

   Out[1] = In[MOMX] / Factor0 / LorentzFactor;
   Out[2] = In[MOMY] / Factor0 / LorentzFactor;
   Out[3] = In[MOMZ] / Factor0 / LorentzFactor;
#  endif

   Out[0] = In[DENS] / LorentzFactor;

   real Temperature;
   SRHydro_HTilde2Temperature ( HTilde, &Temperature, NULL, Gamma );

   Out[4] = Out[0] * Temperature;


   return LorentzFactor;
}// FUNCTION : SRHydro_Con2Pri

#ifdef __CUDACC__
GPU_DEVICE 
void SRHydro_Pri2Con (const real In[], real Out[], const real Gamma)
#else
template <class T> 
void SRHydro_Pri2Con (const T In[], T Out[], const T Gamma)
#endif
{
  real LorentzFactor, Temperature, HTilde, Factor0, Usqr;

  Usqr = VectorDotProduct( In[1], In[2], In[3] );

# ifdef USE_3_VELOCITY
  LorentzFactor = (real)1.0 / SQRT( (real)1.0 - Usqr );
# else
  LorentzFactor = SQRT( (real)1.0 + Usqr );
# endif

  Temperature = In[4]/In[0];

  HTilde = SRHydro_Temperature2HTilde( Temperature, Gamma );

  Out[DENS] = In[0] * LorentzFactor;

  Factor0 = Out[DENS]*HTilde + Out[DENS];

# ifdef USE_3_VELOCITY
  Out[MOMX] = Factor0*In[1]*LorentzFactor;
  Out[MOMY] = Factor0*In[2]*LorentzFactor;
  Out[MOMZ] = Factor0*In[3]*LorentzFactor;
# else
  Out[MOMX] = Factor0*In[1];
  Out[MOMY] = Factor0*In[2];
  Out[MOMZ] = Factor0*In[3];
# endif

  real M_Dsqr = VectorDotProduct( Out[MOMX], Out[MOMY], Out[MOMZ] ) / SQR( Out[DENS] );
  real Fun;

  SRHydro_HTilde_Function( HTilde, M_Dsqr, (real)0.0, &Fun, NULL, Gamma );

  Out[ENGY] = Out[DENS] * Fun;


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
  real ConVar[NCOMP_FLUID], PriVar[NCOMP_FLUID], Lorentz, Vx;

  for (int v=0;v<NCOMP_FLUID;v++) ConVar[v] = Con[v];
  for (int v=0;v<NCOMP_FLUID;v++) PriVar[v] = Pri[v];

  SRHydro_Rotate3D (ConVar, XYZ, true);
  SRHydro_Rotate3D (PriVar, XYZ, true);

# ifdef USE_3_VELOCITY
  Vx = PriVar[1];
# else
  Lorentz = SQRT((real)1.0 + VectorDotProduct(PriVar[1], PriVar[2], PriVar[3]));
  Vx = PriVar[1] / Lorentz;
# endif

  Flux[0] = ConVar[0] * Vx;
  Flux[1] = FMA( ConVar[1], Vx, PriVar[4] );
  Flux[2] = ConVar[2] * Vx;
  Flux[3] = ConVar[3] * Vx;
  Flux[4] = ( ConVar[4] + PriVar[4] )*Vx;

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

  real E_min = factor - D*Dh*MinTemp / factor - D;

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
   real Msqr, M, E_D, M_D;


//--------------------------------------------------------------//
//------------ only check conserved variables-------------------//
//--------------------------------------------------------------//
   if ( Con != NULL && Pri == NULL)
   {
// check NaN
      if (  Con[DENS] != Con[DENS]
         || Con[MOMX] != Con[MOMX]
         || Con[MOMY] != Con[MOMY]
         || Con[MOMZ] != Con[MOMZ]
         || Con[ENGY] != Con[ENGY]  )                                                   goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= Con[DENS] || Con[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Con[MOMX] || Con[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Con[MOMY] || Con[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Con[MOMZ] || Con[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Con[ENGY] || Con[ENGY]  >= (real)HUGE_NUMBER )       goto FAIL;


// check minimum energy
      Msqr = VectorDotProduct( Con[MOMX], Con[MOMY], Con[MOMZ] );
      M = SQRT( Msqr );
	  E_D = Con[ENGY] / Con[DENS];
	  M_D = M / Con[DENS];
      discriminant = ( E_D + M_D ) * ( E_D - M_D ) + (real)2.0 * E_D;

      if ( discriminant <= TINY_NUMBER )                                                goto FAIL;


// pass all checks 
      return false;
   }

//--------------------------------------------------------------//
//------------ only check primitive variables-------------------//
//--------------------------------------------------------------//

   else if ( Con == NULL && Pri != NULL)
   {
// check NaN
      if (  Pri[DENS] != Pri[DENS]
         || Pri[MOMX] != Pri[MOMX]
         || Pri[MOMY] != Pri[MOMY]
         || Pri[MOMZ] != Pri[MOMZ]
         || Pri[ENGY] != Pri[ENGY]  )                                                   goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= Pri[DENS] || Pri[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri[MOMX] || Pri[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri[MOMY] || Pri[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri[MOMZ] || Pri[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Pri[ENGY] || Pri[ENGY]  >= (real)HUGE_NUMBER )       goto FAIL;


#     ifdef USE_3_VELOCITY
// check whether magnitude of 3-velocity is greater or equal to speed of light
      if (VectorDotProduct( Pri[1], Pri[2], Pri[3] ) >= (real) 1.0)                     goto FAIL;
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
                                   Con[DENS], Con[MOMX], Con[MOMY], Con[MOMZ], Con[ENGY]);
              printf( "E^2+2*E*D-|M|^2=%14.7e\n", discriminant );
           }
           else
           {
#             ifdef USE_3_VELOCITY
              printf( "Vx=%14.7e, Vy=%14.7e, Vz=%14.7e, |V|=%14.7e\n",
                                   Pri[1], Pri[2], Pri[3], SQRT( VectorDotProduct( Pri[1], Pri[2], Pri[3] )));
#             else
              printf( "n=%14.7e, Ux=%14.7e, Uy=%14.7e, Uz=%14.7e, P=%14.7e\n", 
                                   Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);
#             endif
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
  real HTilde, E_thermal;

  HTilde = SRHydro_Temperature2HTilde ( Pri[4]/Pri[0], Gamma );

  
  frame ? ( E_thermal = Con[0] * HTilde - Lorentz * Pri[4] )
         :( E_thermal = Pri[0] * HTilde -           Pri[4] );

  return E_thermal;
}

//-------------------------------------------------------------------------------------------------------
// Function    :
// Description : Evaluate kinetic energy density in lab frame 
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real SRHydro_KineticEngy( real Con[], real Pri[], real Lorentz, real Gamma )
{
  real h, Usqr;

  h = SpecificEnthalpy( Con, Pri[4]/Pri[0], Gamma );

# ifdef USE_3_VELOCITY
  real Vsqr = SQR(Pri[1]) + SQR(Pri[2]) + SQR(Pri[3]);

  real LorentzFactor = (real)1.0 / SQRT( (real)1.0 - Vsqr );

  Usqr = Vsqr * SQR(LorentzFactor);
# else
  Usqr = VectorDotProduct( Pri[1], Pri[2], Pri[3]  );
# endif
  
  return ( Con[DENS] * h + Pri[4] ) * Usqr / ( Lorentz + (real)1.0 );
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
NewtonRaphsonSolver(void (*FunPtr)(real, real, real, real*, real*, real), real M_D, real E_D, real *root, const real guess, const real epsabs, const real epsrel, const real Gamma)
{
 int iter = 0;

 #ifdef FLOAT8
 int max_iter = 20;
 #else
 int max_iter = 10;
 #endif

 real Fun, DiffFun;
 real delta;
 real tolerance;
 *root = guess;

 do
   {
     iter++;
     FunPtr(*root, M_D, E_D, &Fun, &DiffFun, Gamma);

#    ifdef CHECK_FAILED_CELL_IN_FLUID
     if ( DiffFun == (real)0.0 )                                                  printf("derivative is zero\n");
     if (  Fun != Fun  ||(real) -HUGE_NUMBER >= Fun  || Fun  >= (real)HUGE_NUMBER )  printf("function value is not finite\n");
     if ( DiffFun != DiffFun ||(real) -HUGE_NUMBER >= DiffFun || DiffFun >= (real)HUGE_NUMBER )  printf("derivative value is not finite\n");
#    endif     

      delta = Fun/DiffFun;
      *root = *root - delta;

      tolerance =  FMA( epsrel, FABS(*root), epsabs );

   }while ( fabs(delta) >= tolerance && iter < max_iter );

}


//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
//-------------------------------------------------------------------------------------------------------

GPU_DEVICE
void SRHydro_HTilde_Function (real HTilde, real M_Dsqr, real Constant, real *Fun, real *DiffFun, real Gamma)
{
  real Temp, DiffTemp;


# if ( EOS == APPROXIMATED_GENERAL )
  SRHydro_HTilde2Temperature ( HTilde, &Temp, &DiffTemp, Gamma );

  real H =  HTilde + (real)1.0;
  real Hsqra  = SQR(H);

  real AA = M_Dsqr/Hsqra;

  real BB = SQRT( (real)1.0 + AA );

  real CC = (real)1.0 + BB;

  if ( Fun != NULL )
  
  *Fun = HTilde * BB + AA / CC - Temp / BB + Constant;
  
  if ( DiffFun != NULL )
  {
    real Hcube  = Hsqra * H;
    real H5     = Hcube * Hsqra;
    
    *DiffFun = SQR(M_Dsqr) / H5 / BB / SQR(CC) 
             - M_Dsqr * HTilde / Hcube / BB
             - (real)2.0 * M_Dsqr / Hcube / CC 
             + BB 
             - DiffTemp / BB 
             - Temp * M_Dsqr / ( Hcube * CUBE(BB) );
  }


# elif ( EOS == CONSTANT_GAMMA )
  SRHydro_HTilde2Temperature ( HTilde, &Temp, &DiffTemp, Gamma  );


  real H =  HTilde + (real)1.0;
  real Factor0 = SQR( H ) + M_Dsqr;

  if ( Fun != NULL )

  *Fun = SQR( HTilde ) + (real)2.0*HTilde - (real)2.0*Temp - (real)2.0*Temp*HTilde
		  + SQR( Temp * H ) / Factor0 + Constant;

  if ( DiffFun != NULL )

  *DiffFun = (real)2.0*H - (real)2.0*Temp - (real)2.0*H*DiffTemp +
		  ( (real)2.0*Temp*DiffTemp*H*H - (real)2.0*Temp*Temp*H ) / SQR( Factor0 );

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
