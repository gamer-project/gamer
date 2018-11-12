#include "GAMER.h"
#include "CUFLU.h"

struct Fun_params
{
  real D;
  real M1;
  real M2;
  real M3;
  real E;
};

// some functions in this file need to be defined even when using GPU
static real Fun (real Q, void *ptr);   // function to be solved
static real DFun (real Q, void *ptr);  // the first derivative of above function
static void Fun_DFun (real Q, void *ptr, real * f, real * df);
void CPU_4Velto3Vel (const real In[], real Out[]);
static real CPU_Con2Temperature (const real In[], const real Gamma); 
bool CPU_CheckUnphysical( const real Con[], const real Pri[]);
static void NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel);


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Rotate3D
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
void
CPU_Rotate3D (real InOut[], const int XYZ, const bool Forward)
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
	  InOut[1] = Temp[1];
	  InOut[2] = Temp[2];
	  InOut[3] = Temp[0];
	  break;
	case 2:
	  InOut[1] = Temp[2];
	  InOut[2] = Temp[0];
	  InOut[3] = Temp[1];
	  break;
	}
    }

  else				// backward
    {
      switch (XYZ)
	{
	case 1:
	  InOut[1] = Temp[2];
	  InOut[2] = Temp[0];
	  InOut[3] = Temp[1];
	  break;
	case 2:
	  InOut[1] = Temp[1];
	  InOut[2] = Temp[2];
	  InOut[3] = Temp[0];
	  break;
	}
    }
}				// FUNCTION : CPU_Rotate3D

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Pri
// Description :  Convert the conserved variables to the primitive variables
//
// Note        :  1. This function always check if the pressure to be returned is greater than the
//                   given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//                   --> But note that here we do NOT ensure "sum(mass fraction) == 1.0"
//                       --> It is done by calling CPU_NormalizePassive() in CPU_Shared_FullStepUpdate()
//
// Parameter   :  In                 : Array storing the input conserved variables
//                Out                : Array to store the output primitive variables
//                Gamma              : Gamma
//                MinPres            : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void
CPU_Con2Pri (const real In[], real Out[], const real Gamma)
{
      real Temp = CPU_GetTemperature (In[0], In[1], In[2], In[3], In[4], NAN, NAN, NAN);
#if ( EOS == RELATIVISTIC_IDEAL_GAS )
      real h = 2.5*Temp + SQRT(2.25*Temp*Temp + 1.0);
#elif ( EOS == IDEAL_GAS ) 
      real h = 1 + Temp * GAMMA / (GAMMA-1.0);
#else
#error: unsupported EoS!
#endif
      real factor = In[0]*h;
      Out[1] = In[1]/factor;
      Out[2] = In[2]/factor;
      Out[3] = In[3]/factor;

      real factor1 = SQRT(1 + SQR (Out[1]) + SQR (Out[2]) + SQR (Out[3]));

      Out[0] = In[0]/factor1;
      Out[4] = Out[0] * Temp; // P = nkT
}// FUNCTION : CPU_Con2Pri


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Temperature
// Description :  Convert the conserved variables to Q
//-------------------------------------------------------------------------------------------------------

static real 
CPU_Con2Temperature (const real In[], const real Gamma)
{
  real Msqr = SQR (In[1]) + SQR (In[2]) + SQR (In[3]);
  real M = SQRT (Msqr); // magnitude of momentum

  real guess, root;
# if ( EOS == RELATIVISTIC_IDEAL_GAS )
# if   ( CONSERVED_ENERGY == 1 )
/* initial guess  */
  real Constant = SQR(In[4]/In[0]) - Msqr/SQR(In[0]);
  if ( Constant > 1.0 )
   {
	if ( 2.0 - Msqr/(16*In[0]*In[0]) >= 0 )
	  {
	    if ( Constant > 2.5 ) 
	      {
		 guess = SQRT( SQR(In[4]/In[0]) - 0.9375*Msqr/SQR(In[0]) - 2.0 ) / 3.0;
	      }
	    else guess = (Constant - 1.0)/ 3.0;
	  }
	else // 1 - (M/D)**2 < 0
	  {
	    if ( Constant >  1.5 + SQRT( Msqr/(16*SQR(In[0])) - 3.0/4.0 )) 
	      {
		 guess = SQRT( SQR(In[4]/In[0]) - 0.9375*Msqr/SQR(In[0]) - 2.0 ) / 3.0;
	      }
	    else guess = (Constant - 1.0)/ 3.0;
	  }
    } else return NAN;
# elif ( CONSERVED_ENERGY == 2 )
/* initial guess  */
   real Constant = SQR(In[4]/In[0]) + 2*In[4]/In[0] - Msqr/SQR(In[0]);
   if ( Constant > 0.0 )
    {
	   if ( 1.0 - Msqr/(16*In[0]*In[0]) >= 0 )
	     {
	       if ( Constant > 1.5 ) 
		 {
		    guess = SQRT( SQR(In[4]/In[0]) + 2*In[4]/In[0] - 0.9375*Msqr/SQR(In[0]) - 1.0 ) / 3.0;
		 }
	       else guess = Constant / 3.0;
	     }
	   else // 1 - (M/D)**2 < 0
	     {
	       if ( Constant >  0.5 + SQRT( Msqr/(16*SQR(In[0])) - 3.0/4.0 )) 
		 {
		    guess = SQRT( SQR(In[4]/In[0]) + 2*In[4]/In[0] - 0.9375*Msqr/SQR(In[0]) - 1.0 ) / 3.0;
		 }
	       else guess = Constant / 3.0;
	     }
    } else return NAN;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
# elif ( EOS == IDEAL_GAS )
  real Gamma_m1 = Gamma - (real) 1.0;
/* initial guess */
# if   ( CONSERVED_ENERGY == 1 )
   real Constant = SQR(In[4]/In[0]) - Msqr/SQR(In[0]);
   if ( Constant > 1.0 )
    {
	      if ( Constant > 1.0 + 2* (Gamma_m1 / Gamma ) * (M / In[0]) )
		{
		  real A = 1.0/SQR(Gamma_m1);
		  real B = 2.0/Gamma_m1;
		  real C = 1 + ((2*Gamma-1.0)/(Gamma*Gamma))*Msqr/(In[0]*In[0]) - SQR(In[4]/In[0]);
		  real delta = SQRT(B*B-4*A*C);
		  guess = -2.0 *  C / ( B + delta);
		}
	     else guess = 0.5*Gamma_m1 * ( Constant - 1.0);
    } else return NAN;
# elif ( CONSERVED_ENERGY == 2 )
    real Constant = SQR(In[4]/In[0]) + 2*In[4]/In[0] - Msqr/SQR(In[0]);
    if ( Constant > 0.0 )
     {
	  if ( Constant > 2* (Gamma_m1 / Gamma ) * (M / In[0]) )
	    {
	      real A = 1.0/SQR(Gamma_m1);
	      real B = 2.0/Gamma_m1;
	      real C = 1 + ((2*Gamma-1.0)/(Gamma*Gamma))*Msqr/(In[0]*In[0]) - SQR(In[4]/In[0]) - 2*(In[4]/In[0]);
	      real delta = SQRT(B*B-4*A*C);
	      guess = -2.0 *  C / ( B + delta);
	    }
	  else guess = 0.5*Gamma_m1 * Constant;
     } else return NAN;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
#else
#error: unsupported EoS!
#endif

   struct Fun_params params = { In[0], In[1], In[2], In[3], In[4] };
  
   NewtonRaphsonSolver(&params ,&root, guess, 0.0, 1.0e-14);

   return root;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Pri2Con
// Description :  Convert the primitive variables to the conserved variables
//
// Note        :  1. This function does NOT check if the input pressure is greater than the
//                   given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//
// Parameter   : [1] In           : Array storing the input primitive variables
//               [2] Out          : Array to store the output conserved variables
//               [3] Gamma        : Adiabatic index
//-------------------------------------------------------------------------------------------------------
void
CPU_Pri2Con (const real In[], real Out[], const real Gamma)
{
#if ( EOS == RELATIVISTIC_IDEAL_GAS )
  real nh = 2.5*In[4] + SQRT(2.25*SQR(In[4]) + SQR(In[0])); // approximate enthalpy * proper number density
#elif ( EOS == IDEAL_GAS )
  real Gamma_m1 = (real) Gamma - 1.0;
  real nh = In[0] + ( GAMMA / Gamma_m1) * In[4]; // enthalpy * proper number density
#else
#error: unsupported EoS!
#endif

  real Factor0 = 1.0 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]);
  real Factor1 = SQRT(Factor0); // Lorentz factor
  real Factor2 = nh * Factor1;
  
  Out[0] = In[0] * Factor1; // number density in inertial frame
  Out[1] = Factor2 * In[1]; // MomX
  Out[2] = Factor2 * In[2]; // MomX
  Out[3] = Factor2 * In[3]; // MomX
# if   ( CONSERVED_ENERGY == 1 )
  Out[4] = nh * Factor0 - In[4]; // total_energy
# elif ( CONSERVED_ENERGY == 2 )
  Out[4] = nh * Factor0 - In[4] - Out[0]; // ( total_energy ) - ( rest_mass_energy )
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
}				// FUNCTION : CPU_Pri2Con

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_4Velto3Vel
// Description :  Convert 4-velocity to 3-velocity
//-------------------------------------------------------------------------------------------------------

void
CPU_4Velto3Vel (const real In[], real Out[])
{
  real Factor = 1 / SQRT (1 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]));

  Out[0] = In[0];
  Out[1] = In[1] * Factor;
  Out[2] = In[2] * Factor;
  Out[3] = In[3] * Factor;
  Out[4] = In[4];
}				// FUNCTION : CPU_4Velto3Vel

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_3Velto4Vel
// Description :  Convert 3-velocity to 4-velocity
//-------------------------------------------------------------------------------------------------------

void
CPU_3Velto4Vel (const real In[], real Out[])
{
  real Factor = 1 / SQRT (1 - SQR (In[1]) - SQR (In[2]) - SQR (In[3]));

  Out[0] = In[0];
  Out[1] = In[1] * Factor;
  Out[2] = In[2] * Factor;
  Out[3] = In[3] * Factor;
  Out[4] = In[4];
}				// FUNCTION : CPU_4Velto3Vel

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Flux
// Description :  Evaluate the hydrodynamic fluxes by the input conserved variables
//
// Parameter   :  XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux     : Array to store the output fluxes
//                Input    : Array storing the input conserved variables
//                Gamma_m1 : Gamma - 1
//                MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void
CPU_Con2Flux (const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres)
{
  const bool CheckMinPres_Yes = true;
  real ConVar[NCOMP_FLUID];	// don't need to include passive scalars since they don't have to be rotated1
  real PriVar4[NCOMP_FLUID];	// D, Ux, Uy, Uz, P
  real PriVar3[NCOMP_FLUID];	// D, Vx, Vy, Vz, P
  real Pres, Vx;
  real lFactor;

  for (int v = 0; v < NCOMP_FLUID; v++)   ConVar[v] = Input[v];

  CPU_Rotate3D (ConVar, XYZ, true);

  CPU_Con2Pri (ConVar, PriVar4, GAMMA);

  CPU_4Velto3Vel (PriVar4, PriVar3);

  Vx = PriVar3[1];
  Pres = PriVar3[4];

  Flux[0] = ConVar[0] * Vx;
  Flux[1] = ConVar[1] * Vx + Pres;
  Flux[2] = ConVar[2] * Vx;
  Flux[3] = ConVar[3] * Vx;
# if ( CONSERVED_ENERGY == 1 )
  Flux[4] = ConVar[1];
# elif ( CONSERVED_ENERGY == 2 )
  Flux[4] = ConVar[1] - Flux[0];
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif

  CPU_Rotate3D (Flux, XYZ, false);
}				// FUNCTION : CPU_Con2Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CheckMinPres
// Description :  Check if the input pressure is great than the minimum allowed threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
//                   --> Currently it simply sets a minimum allowed value for pressure
//                       --> Please set MIN_PRES in the runtime parameter file "Input__Parameter"
//                2. We should also support a minimum **temperature** instead of **pressure**
//                   --> NOT supported yet
//
// Parameter   :  InPres  : Input pressure to be corrected
//                MinPres : Minimum allowed pressure
//
// Return      :  max( InPres, MinPres )
//-------------------------------------------------------------------------------------------------------
real
CPU_CheckMinPres (const real InPres, const real MinPres)
{
  return FMAX (InPres, MinPres);
}				// FUNCTION : CPU_CheckMinPres

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CheckMinDens
//-------------------------------------------------------------------------------------------------------
real
CPU_CheckMinDens (const real InDens, const real MinDens)
{
  return FMAX (InDens, MinDens);
}// FUNCTION : CPU_CheckMinDens


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CheckMinDens
//-------------------------------------------------------------------------------------------------------
real
CPU_CheckMinTemp (const real InTemp, const real MinTemp)
{
  return FMAX (InTemp, MinTemp);
}// FUNCTION : CPU_CheckMinDens


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CheckMinTempInEngy
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
//                Gamma_Ratio : Gamma / (Gamma - 1)
//                MinTemp     : Minimum allowed temperature
//
// Return      :  Total energy with pressure greater than the given threshold
//-------------------------------------------------------------------------------------------------------
real
CPU_CheckMinTempInEngy (const real Cons[])
{
  real T_min = MIN_TEMP;
# if ( EOS == IDEAL_GAS ) 
  real h_min = 1.0 + GAMMA * T_min / (GAMMA - 1.0);
# elif ( EOS == RELATIVISTIC_IDEAL_GAS )
  real h_min = 2.5*T_min + SQRT(2.25*T_min*T_min + 1);
# endif

  real D  = Cons[0];
  real Mx = Cons[1];
  real My = Cons[2];
  real Mz = Cons[3];
  real E  = Cons[4];

  real Msqr = SQR(Mx) + SQR(My) + SQR(Mz);
  real Dh = D*h_min;
  real factor = SQRT(Dh*Dh + Msqr);

# if ( CONSERVED_ENERGY == 1 )
  real E_min = factor - D*Dh*T_min / factor;
# elif ( CONSERVED_ENERGY == 2 )
  real E_min = factor - D*Dh*T_min / factor - D;
# endif

  if ( Cons[4] >= E_min) return Cons[4];
  else                   return E_min;
}



//-------------------------------------------------------------------------------------------------------
// Function    : CPU_CheckUnphysical
//-------------------------------------------------------------------------------------------------------
bool CPU_CheckUnphysical( const real Con[], const real Pri[])
{
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

// check NaN, +inf and -inf
      if (  !Aux_IsFinite(ConsVar[DENS])  
	||  !Aux_IsFinite(ConsVar[MOMX])  
	||  !Aux_IsFinite(ConsVar[MOMY])  
	||  !Aux_IsFinite(ConsVar[MOMZ])  
	||  !Aux_IsFinite(ConsVar[ENGY]) )                                   goto UNPHYSICAL;

// check positivity of number density in inertial frame
      if (ConsVar[DENS] <= 0.0)                                              goto UNPHYSICAL;
 
// check energy
      Msqr = SQR(ConsVar[MOMX]) + SQR(ConsVar[MOMY]) + SQR(ConsVar[MOMZ]);
#     if ( CONSERVED_ENERGY == 1 )
      if ( SQR(ConsVar[ENGY]) <= Msqr + SQR(ConsVar[DENS]) )                 goto UNPHYSICAL;
#     elif ( CONSERVED_ENERGY == 2 )
      if ( SQR(ConsVar[ENGY]) + 2*ConsVar[ENGY]*ConsVar[DENS] - Msqr <= 0 )  goto UNPHYSICAL;
#     else
#     error: CONSERVED_ENERGY must be 1 or 2!
#     endif

      CPU_Con2Pri(ConsVar, Pri4Vel, (real) GAMMA);

// check NaN, +inf and -inf
      if (  !Aux_IsFinite(Pri4Vel[0])  
	||  !Aux_IsFinite(Pri4Vel[1])  
	||  !Aux_IsFinite(Pri4Vel[2])  
	||  !Aux_IsFinite(Pri4Vel[3])  
	||  !Aux_IsFinite(Pri4Vel[4]) )                                       goto UNPHYSICAL;

// check positivity of number density in local rest frame
      if (Pri4Vel[0] <= (real)0.0)                                            goto UNPHYSICAL;
// check positivity of pressure
      if (Pri4Vel[4] <= (real)0.0)                                            goto UNPHYSICAL;
// check whether 3-velocity is greater or equal to speed of light
      CPU_4Velto3Vel(Pri4Vel,Pri3Vel);      
      if (SQR(Pri3Vel[1]) + SQR(Pri3Vel[2]) + SQR(Pri3Vel[3]) >= 1.0)         goto UNPHYSICAL;

// pass all checks 
      return false;
   }

//--------------------------------------------------------------//
//------------ only check primitive variables-------------------//
//--------------------------------------------------------------//

   else if ( Con == NULL && Pri != NULL){
      for(int i=0; i< NCOMP_FLUID; i++) Pri4Vel[i]=Pri[i];


// check NaN, +inf and -inf
      if (  !Aux_IsFinite(Pri4Vel[0])  
	||  !Aux_IsFinite(Pri4Vel[1])  
	||  !Aux_IsFinite(Pri4Vel[2])  
	||  !Aux_IsFinite(Pri4Vel[3])  
	||  !Aux_IsFinite(Pri4Vel[4]) )                                            goto UNPHYSICAL;

// check positivity of number density in local rest frame
      if (Pri4Vel[0] <= (real)0.0)                                                 goto UNPHYSICAL;
// check positivity of pressure
      if (Pri4Vel[4] <= (real)0.0)                                                 goto UNPHYSICAL;

      CPU_4Velto3Vel(Pri4Vel,Pri3Vel);
      CPU_Pri2Con(Pri4Vel, ConsVar, (real) GAMMA);

// check whether 3-velocity is greater or equal to speed of light
      if (SQR(Pri3Vel[1]) + SQR(Pri3Vel[2]) + SQR(Pri3Vel[3]) >= 1.0)              goto UNPHYSICAL;
   
// check NaN, +inf and -inf
      if (  !Aux_IsFinite(ConsVar[DENS])  
	||  !Aux_IsFinite(ConsVar[MOMX])  
	||  !Aux_IsFinite(ConsVar[MOMY])  
	||  !Aux_IsFinite(ConsVar[MOMZ])  
	||  !Aux_IsFinite(ConsVar[ENGY]) )                                               goto UNPHYSICAL;

// check positivity of number density in inertial frame
      if (ConsVar[DENS] <= 0.0)                                                          goto UNPHYSICAL;

// check energy
      Msqr = SQR(ConsVar[MOMX]) + SQR(ConsVar[MOMY]) + SQR(ConsVar[MOMZ]);
#     if ( CONSERVED_ENERGY == 1 )
      if ( SQR(ConsVar[ENGY]) <= Msqr + SQR(ConsVar[DENS]) )                             goto UNPHYSICAL;
#     elif ( CONSERVED_ENERGY == 2 )
      if ( SQR(ConsVar[ENGY]) + 2*ConsVar[ENGY]*ConsVar[DENS] - Msqr <= 0 )              goto UNPHYSICAL;
#     else
#     error: CONSERVED_ENERGY must be 1 or 2!
#     endif      

// pass all checks 
      return false;
// print all variables if goto UNPHYSICAL
      UNPHYSICAL:
      {
        Aux_Message(stderr, "\n\nD=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n",
                             ConsVar[DENS], ConsVar[MOMX], ConsVar[MOMY], ConsVar[MOMZ], ConsVar[ENGY]);
#       if ( CONSERVED_ENERGY == 1 )
        Aux_Message(stderr, "E^2-|M|^2-D^2=%14.7e\n", SQR(ConsVar[ENGY])-Msqr-SQR(ConsVar[DENS]));
#       elif ( CONSERVED_ENERGY == 2 )
        Aux_Message(stderr, "E^2+2*E*D-|M|^2=%14.7e\n", SQR(ConsVar[ENGY]) + 2*ConsVar[ENGY]*ConsVar[DENS] - Msqr);
#       endif
        Aux_Message(stderr, "n=%14.7e, Ux=%14.7e, Uy=%14.7e, Uz=%14.7e, P=%14.7e\n", 
                             Pri4Vel[0], Pri4Vel[1], Pri4Vel[2], Pri4Vel[3], Pri4Vel[4]);
        Aux_Message(stderr, "Vx=%14.7e, Vy=%14.7e, Vz=%14.7e, |V|=%14.7e\n",
                             Pri3Vel[1], Pri3Vel[2], Pri3Vel[3], SQRT(SQR(Pri3Vel[1])+SQR(Pri3Vel[2])+SQR(Pri3Vel[3])));
        return true;
      }
   }
   else
   {
    Aux_Error(ERROR_INFO,"One of Con or Pri must be given in CPU_CheckUnphysical!\n");
    return true;
   }
}
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CheckNegative
// Description :  Check whether the input value is <= 0.0 (also check whether it's Inf or NAN)
//
// Note        :  Can be used to check whether the values of density and pressure are unphysical
//
// Parameter   :  Input : Input value
//
// Return      :  true  --> Input <= 0.0  ||  >= __FLT_MAX__  ||  != itself (Nan)
//                false --> otherwise
//-------------------------------------------------------------------------------------------------------
bool
CPU_CheckNegative (const real Input)
{
  if (Input < (real) 0.0 || Input >= HUGE_NUMBER || Input != Input)
    return true;
  else
    return false;
}				// FUNCTION : CPU_CheckNegative

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_GetPressure
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
//                Gamma_m1     : Gamma - 1, where Gamma is the adiabatic index
//                CheckMinPres : Return CPU_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                                   --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres      : Minimum allowed pressure
//
// Return      :  Pressure
//-------------------------------------------------------------------------------------------------------
real
CPU_GetPressure (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
		 const real Gamma_m1, const bool CheckMinPres, const real MinPres)
{
  real In[NCOMP_FLUID];
  real Out[NCOMP_FLUID];
  real Pres;

  In[0] = Dens;
  In[1] = MomX;
  In[2] = MomY;
  In[3] = MomZ;
  In[4] = Engy;

  CPU_Con2Pri (In, Out, GAMMA);

  Pres = Out[4];

  return Pres;
}				// FUNCTION : CPU_GetPressure



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_GetTemperature
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
//                Gamma_m1     : Gamma - 1, where Gamma is the adiabatic index
//                CheckMinPres : Return CPU_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                                   --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres      : Minimum allowed pressure
//
// Return      :  Temperature
//-------------------------------------------------------------------------------------------------------
real
CPU_GetTemperature (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                    const real Gamma_m1, const bool CheckMinPres, const real MinPres )
{
      real In[5] = {Dens, MomX, MomY, MomZ, Engy};
      real Temperature = CPU_Con2Temperature ( In, GAMMA );

 return Temperature;
}				// FUNCTION : CPU_GetTemperature



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Temperature2Pressure
// Description :  Convert gas temperature to pressure
//
// Note        :  1. Assume the ideal gas law
//                   --> P = \rho*K*T / ( mu*m_H )
//                2. Assume both input and output to be code units
//                   --> Temperature should be converted to UNIT_E in advance
//                       --> Example: T_code_unit = T_kelvin * Const_kB / UNIT_E
//                3. Pressure floor (MinPres) is applied when CheckMinPres == true
//                4. Currently this function always adopts real precision since
//                   (1) both Temp and m_H may exhibit extreme values depending on the code units, and
//                   (2) we don't really care about the performance here since this function is usually
//                       only used for constructing the initial condition
//
// Parameter   :  Dens         : Gas mass density in code units
//                Temp         : Gas temperature in code units
//                mu           : Mean molecular weight
//                m_H          : Atomic hydrogen mass in code units
//                               --> Sometimes we use the atomic mass unit (Const_amu defined in PhysicalConstant.h)
//                                   and m_H (Const_mH defined in PhysicalConstant.h) interchangeably since the
//                                   difference is small (m_H ~ 1.007825 amu)
//                CheckMinPres : Return CPU_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Minimum allowed pressure
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
real
CPU_Temperature2Pressure (const real Dens, const real Temp, const real mu, const real m_H, const bool CheckMinPres, const real MinPres)
{
  Aux_Message (stderr,"\n\nWARNING:\nfile: %s\n", __FILE__);
  Aux_Message (stderr,"\n\nPlease modify %s properly.\n", __FUNCTION__);
  abort ();

  real Pres;

  return Pres;

}				// FUNCTION : CPU_GetTemperature



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_NormalizePassive
// Description :  Normalize the target passive scalars so that the sum of their mass density is equal to
//                the gas mass density
//
// Note        :  1. Should be invoked AFTER applying the floor values to passive scalars
//                2. Invoked by CPU_Shared_FullStepUpdate(), Prepare_PatchData(), Refine(), LB_Refine_AllocateNewPatch(),
//                   Flu_FixUp(), XXX_Init_ByFunction_AssignData(), Flu_Close()
//
// Parameter   :  GasDens : Gas mass density
//                Passive : Passive scalar array (with the size NCOMP_PASSIVE)
//                NNorm   : Number of passive scalars to be normalized
//                          --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx : Target variable indices to be normalized
//                          --> Should be set to the global variable "PassiveNorm_VarIdx"
//
// Return      :  Passive
//-------------------------------------------------------------------------------------------------------
void
CPU_NormalizePassive (const real GasDens, real Passive[], const int NNorm, const int NormIdx[])
{
  Aux_Message (stderr,"\n\nWARNING:\nPlease modify %s properly.\n", __FUNCTION__);
  Aux_Message (stderr,"file: %s\n", __FILE__);
  Aux_Message (stderr,"line: %d\n", __LINE__);
  abort();

// validate the target variable indices
#ifdef GAMER_DEBUG
  const int MinIdx = 0;
#ifdef DUAL_ENERGY
  const int MaxIdx = NCOMP_PASSIVE - 2;
#else
  const int MaxIdx = NCOMP_PASSIVE - 1;
#endif

  for (int v = 0; v < NNorm; v++)
    {
      if (NormIdx[v] < MinIdx || NormIdx[v] > MaxIdx)
	Aux_Error (ERROR_INFO, "NormIdx[%d] = %d is not within the correct range ([%d <= idx <= %d]) !!\n", v, NormIdx[v], MinIdx, MaxIdx);
    }
#endif // #ifdef GAMER_DEBUG


  real Norm, PassiveDens_Sum = (real) 0.0;

  for (int v = 0; v < NNorm; v++)
    PassiveDens_Sum += Passive[NormIdx[v]];

  Norm = GasDens / PassiveDens_Sum;

  for (int v = 0; v < NNorm; v++)
    Passive[NormIdx[v]] *= Norm;

}				// FUNCTION : CPU_NormalizePassive

static void NewtonRaphsonSolver(void *ptr, real *root, const real guess, const real epsabs, const real epsrel)
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
     Fun_DFun(*root, ptr, &f, &df);

     if (df == 0.0)                printf("derivative is zero\n");
     if ( Aux_IsFinite(f) == 0 )   printf("function value is not finite\n");
     if ( Aux_IsFinite(df) == 0 )  printf("derivative value not finite\n");
     
      root_old = *root;
      *root = *root - ( f / df );
      tolerance = epsabs + epsrel * FABS(*root);
   }while ( fabs(root_old - *root) >= tolerance && iter < max_iter );
}


//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
//-------------------------------------------------------------------------------------------------------

static void
Fun_DFun (real Q, void *ptr, real * f, real * df)
{
  struct Fun_params *params = (struct Fun_params *) ptr;

  real D  = (params->D);
  real M1 = (params->M1);
  real M2 = (params->M2);
  real M3 = (params->M3);
  real E  = (params->E);

#if ( EOS == RELATIVISTIC_IDEAL_GAS )
  real h = 2.5*Q+SQRT(2.25*SQR(Q)+1.0); // approximate enthalpy
  real dh = 2.5 + 9.0*Q / SQRT(36*Q*Q+16);
  real Msqr = SQR(M1) + SQR(M2) + SQR(M3);
  real Factor = Msqr / (D*D);
  real hsqr = SQR(h);
  real hQ = h * Q;
#if   (CONSERVED_ENERGY == 1)
  real Constant = SQR(E/D) - Factor;
  //*f = hsqr- 2*hQ + (SQR(hQ)) / (hsqr + Factor) - Constant;
  *f = 3.5 * Q * Q + 1.5 * Q * SQRT(9*Q*Q+4) + SQR(hQ) / (hsqr + Factor) + 1.0 - Constant;
#elif (CONSERVED_ENERGY == 2)
  real Constant = SQR(E/D) + 2*(E/D) - Factor;
  //*f = hsqr - 1.0 - 2*hQ + (SQR(hQ)) / (hsqr + Factor) - Constant;
  *f = 3.5 * Q * Q + 1.5 * Q * SQRT(9*Q*Q+4) + SQR(hQ) / (hsqr + Factor) - Constant;
# else
# error: CONSERVED_ENERGY must be 1 or 2!
#endif
  real abc = SQRT(9*Q*Q+4);
  //*df = 2*dh*(h-Q) - 2*h + 2*hQ*((hsqr*h + h * Factor + Q*dh*Factor) / SQR(hsqr + Factor) );
  *df = 7*Q + 1.5 * abc + 13.5*Q*Q/abc + 2*h*Q*((h*hsqr + Factor*h + Q*dh*Factor) / SQR( hsqr + Factor) );

#elif ( EOS == IDEAL_GAS )
  real alpha = GAMMA / (GAMMA - 1.0);
  real h = 1 + alpha * Q;
  real hsqr = SQR(h);
  real Msqr = SQR(M1) + SQR(M2) + SQR(M3);
  real Factor = Msqr / (D*D);
  real beta = (GAMMA * (2-GAMMA)) / ((GAMMA-1.0)*(GAMMA-1.0));
  real hQ = h * Q;

#if   (CONSERVED_ENERGY == 1)
  real Constant = SQR(E/D) - Factor;
  *f = 1.0 + beta * SQR(Q) + 2 * Q / (GAMMA - 1.0) + SQR(hQ) / (hsqr + Factor) - Constant; 
#elif (CONSERVED_ENERGY == 2)
  real Constant = SQR(E/D) + 2*(E/D) - Factor;
  *f = beta * SQR(Q) + 2 * Q / (GAMMA - 1.0) + SQR(hQ) / (hsqr + Factor) - Constant; 
  //*f = h*h-1-2*h*Q + SQR(hQ) / (hsqr + Factor) - Constant; 
#endif
  real dh = alpha;
  *df = 2 * beta * Q + 2/(GAMMA-1.0) + 2*h*Q*((h*hsqr + Factor*h + Q*dh*Factor) / SQR( hsqr + Factor) );
#else
#error: unsupported EoS!
#endif // #if ( EOS == RELATIVISTIC_IDEAL_GAS )
}
