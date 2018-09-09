#include "GAMER.h"
#include "CUFLU.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>


#if ( MODEL == SR_HYDRO )
struct FUN_Q_params
{
  real d;
  real M1;
  real M2;
  real M3;
  real E;
  real Gamma;
};

// some functions in this file need to be defined even when using GPU

static real FUN_Q (real Q, void *ptr);
static real D_FUN_Q (real Q, void *ptr);
static void FDF_FUN_Q (real Q, void *ptr, real * f, real * df);


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
  real In_temp[5] = { In[0], In[1], In[2], In[3], In[4] };
  real Gamma_m1 = Gamma - (real) 1.0;

#  ifdef CHECK_NEGATIVE_IN_FLUID
  if ( CPU_CheckNegative(In_temp[0]) 
    ||     !Aux_IsFinite(In_temp[1]) 
    ||     !Aux_IsFinite(In_temp[2]) 
    ||     !Aux_IsFinite(In_temp[3]) 
    || CPU_CheckNegative(In_temp[4]))
    {
       Aux_Message (stderr, "\n\nWARNING:\n :file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
       Aux_Message (stderr, "line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, 
			    In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
    }
#  endif

  real M = SQRT (SQR (In_temp[1]) + SQR (In_temp[2]) + SQR (In_temp[3]));

#  ifdef CHECK_NEGATIVE_IN_FLUID
  if ( In_temp[4] <= M )
    {
      Aux_Message (stderr,"\n\nWARNING: |M| > E!\n");
      Aux_Message (stderr,"file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, 
                           In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
      Aux_Message (stderr,"|M|=%e, E=%e, |M|-E=%e\n\n", M, In_temp[4], M - In_temp[4]);
    }
#  endif

  if (fabs (M) > TINY_NUMBER)
    {
      int status;

      int iter = 0;
      int max_iter = 200;

      const gsl_root_fdfsolver_type *T;

      gsl_root_fdfsolver *s;

      real Q, Q0;

/* initial guess Q  */
      if (In_temp[0] > M / Gamma_m1)
	{
	  Q = M * (In_temp[4] - M) / ((1 - 1 / Gamma) * In_temp[0]);
	}
      else
	{
	  Q = In_temp[4] * Gamma;
	}

      gsl_function_fdf F;

      struct FUN_Q_params params = { In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4], Gamma };

      F.f = &FUN_Q;
      F.df = &D_FUN_Q;
      F.fdf = &FDF_FUN_Q;
      F.params = &params;

      T = gsl_root_fdfsolver_newton;
      s = gsl_root_fdfsolver_alloc (T);
      gsl_root_fdfsolver_set (s, &F, Q);

      //printf ("status = %s\n", gsl_strerror (status)); 
      do
	{
	  iter++;
	  status = gsl_root_fdfsolver_iterate (s);
	  //printf ("status = %s\n", gsl_strerror (status));
	  Q0 = Q;
	  Q = gsl_root_fdfsolver_root (s);
	  status = gsl_root_test_delta (Q, Q0, 0, 1e-16);
	  //printf ("status = %s\n", gsl_strerror (status));
	}
      while (status == GSL_CONTINUE && iter < max_iter);
      //printf ("status = %s\n", gsl_strerror (status));

#  ifdef CHECK_NEGATIVE_IN_FLUID
     if (Q < 0)
       {
	 Aux_Message (stderr,"\n\nWARNING: Q = %e < 0!\n", Q);
	 Aux_Message (stderr,"file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
	 Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, 
                              In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
       }
#  endif

      Out[1] = In_temp[1] / Q;	/*Ux */
      Out[2] = In_temp[2] / Q;	/*Uy */
      Out[3] = In_temp[3] / Q;	/*Uz */

      real Factor = SQRT (1 + SQR (Out[1]) + SQR (Out[2]) + SQR (Out[3]));

      Out[0] = In_temp[0] / Factor;	/*primitive density */
      Out[4] = (Gamma_m1 / Gamma) * FABS (Q / Factor - Out[0]);	/*pressure */

#ifdef CHECK_NEGATIVE_IN_FLUID
      if ((Out[4] < 0) || (Out[0] < 0))
	{
	  Aux_Message (stderr,"\n\nWARNING: P < 0 or d < 0!\n");
	  Aux_Message (stderr,"file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
	  Aux_Message (stderr,"Q = %e\n", Q);
          Aux_Message (stderr,"P = %e, d = %e\n", Out[4], Out[0]);
	  Aux_Message (stderr,"D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
	}
#endif

      gsl_root_fdfsolver_free (s);
    } // if (fabs (M) > TINY_NUMBER)
#ifndef CHECK_NEGATIVE_IN_FLUID
  else
#else
  else if (In_temp[4] >= In_temp[0])
#endif
    {
      Out[1] = 0.0;
      Out[2] = 0.0;
      Out[3] = 0.0;
      Out[0] = In_temp[0];
      Out[4] = Gamma_m1 * (In_temp[4] - In_temp[0]);

#ifdef CHECK_NEGATIVE_IN_FLUID
      if (Out[4] < 0)
	{
	  Aux_Message (stderr,"\n\nWARNING: P < 0!\n");
	  Aux_Message (stderr,"file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
          Aux_Message (stderr,"P = %e\n", Out[4]);
	  Aux_Message (stderr,"D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
	} // if (Out[4] < 0)
#endif

    } // else if (In_temp[4] >= In_temp[0])

#ifdef CHECK_NEGATIVE_IN_FLUID
  else
    {
      Aux_Message (stderr,"\n\nWARNING:\nToo criticle to solve! Somthing went wrong!\n");
      Aux_Message (stderr,"line:%d, D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
    } // if (fabs (M) > TINY_NUMBER)
#endif

}				// FUNCTION : CPU_Con2Pri



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
// Parameter   :  In          : Array storing the input primitive variables
//                Out         : Array to store the output conserved variables
//               _Gamma_m1    : 1 / (Gamma - 1)
//                NormPassive : true --> convert passive scalars to mass fraction
//                NNorm       : Number of passive scalars for the option "NormPassive"
//                              --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx     : Target variable indices for the option "NormPassive"
//                              --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
void
CPU_Pri2Con (const real In[], real Out[], const real Gamma)
{
  real Gamma_m1 = (real) Gamma - 1.0;
#ifdef CHECK_NEGATIVE_IN_FLUID
  if ( CPU_CheckNegative(In[0]) 
    ||     !Aux_IsFinite(In[1]) 
    ||     !Aux_IsFinite(In[2]) 
    ||     !Aux_IsFinite(In[3]) 
    || CPU_CheckNegative(In[4]))
     {
	Aux_Message (stderr,"\n\nWARNING:\nfile: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
	Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In[0], In[1], In[2], In[3], In[4]);
     }
#endif

  real Factor0 = 1 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]);
  real Factor1 = SQRT (Factor0);
  real Factor2 = Gamma / Gamma_m1;
  real Factor3 = In[0] + Factor2 * In[4]; // enthalpy* rho
  real Factor4 = Factor3 * Factor1;

  Out[0] = In[0] * Factor1;
  Out[1] = Factor4 * In[1];
  Out[2] = Factor4 * In[2];
  Out[3] = Factor4 * In[3];
  Out[4] = Factor3 * Factor0 - In[4];

#  ifdef CHECK_NEGATIVE_IN_FLUID
  if ( CPU_CheckNegative(Out[DENS]) 
    ||     !Aux_IsFinite(Out[MOMX]) 
    ||     !Aux_IsFinite(Out[MOMY]) 
    ||     !Aux_IsFinite(Out[MOMZ]) 
    || CPU_CheckNegative(Out[ENGY]))
  {
     Aux_Message (stderr,"\n\nWARNING:\nfile: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
     Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__,
                  Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY]);
  }
  real  M = SQRT (SQR (Out[MOMX]) + SQR (Out[MOMY]) + SQR (Out[MOMZ]));

  if ( Out[ENGY] <= M ) 
    {   
      Aux_Message (stderr,"\n\nWARNING: |M| > E!\n");
      Aux_Message (stderr,"file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, 
                   Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY]);
      Aux_Message (stderr,"|M|=%e, E=%e, |M|-E=%e\n\n", M, Out[ENGY], M - Out[ENGY]);
      Aux_Message (stderr,"d=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n", In[0], In[1], In[2], In[3], In[4]);
    }   
#  endif

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

#ifdef CHECK_NEGATIVE_IN_FLUID
  real V = SQRT(SQR(Out[1])+SQR(Out[2])+SQR(Out[3]));
      if (V >= 1.0)
	{
	  Aux_Message (stderr,"\n\nWARNING: |V| >= 1!\n");
	  Aux_Message (stderr,"file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
	  Aux_Message (stderr,"Ux=%e, Uy=%e, Uz=%e\n", In[1], In[2], In[3]);
	  Aux_Message (stderr,"Vx=%e, Vy=%e, Vz=%e\n", Out[1], Out[2], Out[3]);
	  Aux_Message (stderr,"|V| = %e\n", V);
	}
#endif

}				// FUNCTION : CPU_4Velto3Vel

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_3Velto4Vel
// Description :  Convert 3-velocity to 4-velocity
//-------------------------------------------------------------------------------------------------------

void
CPU_3Velto4Vel (const real In[], real Out[])
{
#ifdef CHECK_NEGATIVE_IN_FLUID
  real V = SQRT( SQR (In[1]) + SQR (In[2]) + SQR (In[3]) );
  if ( V >= 1.0 )
    {
      Aux_Message (stderr,"\n\nWARNING:%s: %d\n", __FUNCTION__, __LINE__);
      Aux_Message (stderr,"\n\n|V| >= 1!\n");
      Aux_Message (stderr,"V = %e >= 1.0 \n", V);
    }
#endif

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
#  ifdef CHECK_NEGATIVE_IN_FLUID
  if ( CPU_CheckNegative(Input[0]) 
    ||     !Aux_IsFinite(Input[1]) 
    ||     !Aux_IsFinite(Input[2]) 
    ||     !Aux_IsFinite(Input[3]) 
    || CPU_CheckNegative(Input[4]))
  {
     Aux_Message (stderr,"\n\nWARNING:\nfile: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
     Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, Input[0], Input[1], Input[2], Input[3], Input[4]);
  }
  real M = SQRT (SQR (Input[1]) + SQR (Input[2]) + SQR (Input[3]));

  if ( Input[4] <= M )
    {
      Aux_Message (stderr,"\n\nWARNING: |M| > E!\n");
      Aux_Message (stderr,"file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, Input[0], Input[1], Input[2], Input[3], Input[4]);
      Aux_Message (stderr,"|M|=%e, E=%e, |M|-E=%e\n\n", M, Input[4], M - Input[4]);
    }
#  endif
  const bool CheckMinPres_Yes = true;
  real ConVar[NCOMP_FLUID];	// don't need to include passive scalars since they don't have to be rotated1
  real PriVar4[NCOMP_FLUID];	// d, Ux, Uy, Uz, P
  real PriVar3[NCOMP_FLUID];	// d, Vx, Vy, Vz, P
  real Pres, Vx;
  real lFactor;

  for (int v = 0; v < NCOMP_FLUID; v++)
    ConVar[v] = Input[v];

  CPU_Rotate3D (ConVar, XYZ, true);

  CPU_Con2Pri (ConVar, PriVar4, GAMMA);

#ifdef CHECK_NEGATIVE_IN_FLUID
      if ((PriVar4[4] < 0) || (PriVar4[0] < 0))
	{
	  Aux_Message (stderr,"\n\nWARNING: P < 0 or d < 0!\n");
	  Aux_Message (stderr,"file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
          Aux_Message (stderr,"P = %e, d = %e\n", PriVar4[4], PriVar4[0]);
	  Aux_Message (stderr,"D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", ConVar[0], ConVar[1], ConVar[2], ConVar[3], ConVar[4]);
	}
#endif

  CPU_4Velto3Vel (PriVar4, PriVar3);

#ifdef CHECK_NEGATIVE_IN_FLUID
  real V = SQRT(SQR(PriVar3[1])+SQR(PriVar3[2])+SQR(PriVar3[3]));
      if (V >= 1.0)
	{
	  Aux_Message (stderr,"\n\nWARNING: |V| >= 1!\n");
	  Aux_Message (stderr,"file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
	  Aux_Message (stderr,"Ux=%e, Uy=%e, Uz=%e\n", PriVar4[1], PriVar4[2], PriVar4[3]);
	  Aux_Message (stderr,"Vx=%e, Vy=%e, Vz=%e\n", PriVar3[1], PriVar3[2], PriVar3[3]);
	  Aux_Message (stderr,"|V| = %e\n", V);
	}
#endif


  Vx = PriVar3[1];
  Pres = PriVar3[4];

  Flux[0] = ConVar[0] * Vx;
  Flux[1] = ConVar[1] * Vx + Pres;
  Flux[2] = ConVar[2] * Vx;
  Flux[3] = ConVar[3] * Vx;
  Flux[4] = ConVar[1];

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
#ifdef CHECK_NEGATIVE_IN_FLUID
  if ( CPU_CheckNegative(InPres))
    {
       Aux_Message (stderr,"\n\nWARNING:\nfile: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
       Aux_Message (stderr,"line:%d\nInPres = %e\n", __LINE__, InPres);
    }
#endif
  return FMAX (InPres, MinPres);
}				// FUNCTION : CPU_CheckMinPres



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CheckMinPresInEngy
// Description :  Ensure that the pressure in the input total energy is greater than the given threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
//                   --> Currently it simply sets a minimum allowed value for pressure
//                       --> Please set MIN_PRES in the runtime parameter file "Input__Parameter"
//                3. One must input conserved variables instead of primitive variables
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                Engy     : Energy density
//                Gamma_m1 : Gamma - 1
//               _Gamma_m1 : 1/(Gamma - 1)
//                MinPres  : Minimum allowed pressure
//
// Return      :  Total energy with pressure greater than the given threshold
//-------------------------------------------------------------------------------------------------------
real
CPU_ModifyEngy (const real MomX, const real MomY, const real MomZ, const real Engy, const double epsilon)
{
#  ifdef CHECK_NEGATIVE_IN_FLUID
  if (     !Aux_IsFinite(MomX) 
    ||     !Aux_IsFinite(MomY)
    ||     !Aux_IsFinite(MomZ)
    || CPU_CheckNegative(Engy))
   {
      Aux_Message (stderr,"\n\nWARNING:\nfile: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      Aux_Message (stderr,"line:%d\nMx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, MomX, MomY, MomZ, Engy);
   }
#endif
  real Msqr = SQR (MomX) + SQR (MomY) + SQR (MomZ);
  real M = SQRT (Msqr);

  if (Engy > M)
    {
      return Engy;
    }
  else
    {
      return M * (1.0 + epsilon);
    }
}

real
CPU_CheckMinDens (const real InDens, const real MinDens)
{
  return FMAX (InDens, MinDens);
}				// FUNCTION : CPU_CheckMinDens

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
#ifndef FLOAT8			// single precession
bool
CPU_CheckNegative (const real Input)
{
  if (Input < (real) 0.0 || Input >= __FLT_MAX__ || Input != Input)
    return true;
  else
    return false;

}				// FUNCTION : CPU_CheckNegative
#else				// double precession
bool
CPU_CheckNegative (const real Input)
{
  if (Input < (real) 0.0 || Input >= __DBL_MAX__ || Input != Input)
    return true;
  else
    return false;

}				// FUNCTION : CPU_CheckNegative
#endif // #ifndef FLOAT8



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
#  ifdef CHECK_NEGATIVE_IN_FLUID
  if ( CPU_CheckNegative(In[0]) 
    ||     !Aux_IsFinite(In[1]) 
    ||     !Aux_IsFinite(In[2]) 
    ||     !Aux_IsFinite(In[3]) 
    || CPU_CheckNegative(In[4]))
  {
     Aux_Message (stderr,"\n\nWARNING:\nfile: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
     Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In[0], In[1], In[2], In[3], In[4]);
  }
  real M = SQRT (SQR (In[1]) + SQR (In[2]) + SQR (In[3]));

  if ( In[4] <= M )
    {
      Aux_Message (stderr,"\n\nWARNING: |M| > E!\n");
      Aux_Message (stderr,"file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      Aux_Message (stderr,"line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In[0], In[1], In[2], In[3], In[4]);
      Aux_Message (stderr,"|M|=%e, E=%e, |M|-E=%e\n\n", M, In[4], M - In[4]);
    }
#  endif

  CPU_Con2Pri (In, Out, GAMMA);

  Pres = Out[4];

  if ( CheckMinPres )   Pres = CPU_CheckMinPres( Pres, MinPres );

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
		    const real Gamma_m1, const bool CheckMinPres, const real MinPres)
{
  Aux_Message (stderr,"\n\nWARNING:\nfile: %s\n", __FILE__);
  Aux_Message (stderr,"\n\nPlease modify %s properly.\n", __FUNCTION__);
  abort ();
  return CPU_GetPressure (Dens, MomX, MomY, MomZ, Engy, Gamma_m1, CheckMinPres, MinPres) / Dens;

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

  Pres = Dens * Temp / (mu * m_H);


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

static real
FUN_Q (real Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  real d = (params->d);
  real M1 = (params->M1);
  real M2 = (params->M2);
  real M3 = (params->M3);
  real E = (params->E);
  real Gamma = (params->Gamma);
  real Gamma_m1 = (real) Gamma - 1.0;
/*4-Velocity*/
  real U1 = M1 / Q;
  real U2 = M2 / Q;
  real U3 = M3 / Q;

  real rho = d / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3));

  real pres = (Gamma_m1 / Gamma) * (Q / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  real f = Q * SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;

  return f;
}

static real
D_FUN_Q (real Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  real d = (params->d);
  real M1 = (params->M1);
  real M2 = (params->M2);
  real M3 = (params->M3);
  real Gamma = (params->Gamma);
  real Gamma_m1 = (real) Gamma - 1.0;

  real U1 = M1 / Q;
  real U2 = M2 / Q;
  real U3 = M3 / Q;

  real dU1 = -M1 / (Q * Q);
  real dU2 = -M2 / (Q * Q);
  real dU3 = -M3 / (Q * Q);

  real dd = -d * POW (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  real dp = (Gamma_m1 / Gamma) * (1 / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3))
				  - Q * POW (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - dd);

  real df = SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) - dp;

  return df;
}

static void
FDF_FUN_Q (real Q, void *ptr, real * f, real * df)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  real d = (params->d);
  real M1 = (params->M1);
  real M2 = (params->M2);
  real M3 = (params->M3);
  real E = (params->E);
  real Gamma = (params->Gamma);
  real Gamma_m1 = (real) Gamma - 1.0;

  real U1 = M1 / Q;
  real U2 = M2 / Q;
  real U3 = M3 / Q;

  real rho = d / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3));
  real pres = (Gamma_m1 / Gamma) * (Q / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  real dU1 = -M1 / (Q * Q);
  real dU2 = -M2 / (Q * Q);
  real dU3 = -M3 / (Q * Q);

  real drho = -d * POW (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  real dp = (Gamma_m1 / Gamma) * (1 / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3))
				  - Q * POW (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - drho);

  *f = Q * SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;
  *df = SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) / SQRT (1 + SQR (U1) + SQR (U2) + SQR (U3)) - dp;
}

#endif // #if ( MODEL == SR_HYDRO )
