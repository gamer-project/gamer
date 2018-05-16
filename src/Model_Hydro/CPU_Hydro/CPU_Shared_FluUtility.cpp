#include "GAMER.h"
#include "CUFLU.h"

// some functions in this file need to be defined even when using GPU
#if ( MODEL == HYDRO )




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
void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward )
{

   if ( XYZ == 0 )   return;


   real Temp[3];
   for (int v=0; v<3; v++)    Temp[v] = InOut[v+1];

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : InOut[1] = Temp[1];  InOut[2] = Temp[2];  InOut[3] = Temp[0];     break;
         case 2 : InOut[1] = Temp[2];  InOut[2] = Temp[0];  InOut[3] = Temp[1];     break;
      }
   }

   else // backward
   {
      switch ( XYZ )
      {
         case 1 : InOut[1] = Temp[2];  InOut[2] = Temp[0];  InOut[3] = Temp[1];     break;
         case 2 : InOut[1] = Temp[1];  InOut[2] = Temp[2];  InOut[3] = Temp[0];     break;
      }
   }

} // FUNCTION : CPU_Rotate3D



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
//                Gamma_m1           : Gamma - 1
//                MinPres            : Minimum allowed pressure
//                NormPassive        : true --> convert passive scalars to mass fraction
//                NNorm              : Number of passive scalars for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                  const bool NormPassive, const int NNorm, const int NormIdx[],
                  const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const bool CheckMinPres_Yes = true;
   const real _Rho             = (real)1.0 / In[0];

   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = CPU_GetPressure( In[0], In[1], In[2], In[3], In[4], Gamma_m1, CheckMinPres_Yes, MinPres );

// pressure floor required to resolve the Jeans length
// --> note that currently we do not modify the dual-energy variable (e.g., entropy) accordingly
   if ( JeansMinPres )
   Out[4] = CPU_CheckMinPres( Out[4], JeansMinPres_Coeff*SQR(Out[0]) );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass density of target passive scalars to mass fraction
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Out[ NCOMP_FLUID + NormIdx[v] ] *= _Rho;
#  endif

} // FUNCTION : CPU_Con2Pri



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
void CPU_Pri2Con( const real In[], real Out[], const real _Gamma_m1,
                  const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];
   Out[4] = In[4]*_Gamma_m1 + (real)0.5*In[0]*( In[1]*In[1] + In[2]*In[2] + In[3]*In[3] );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass fraction of target passive scalars back to mass density
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Out[ NCOMP_FLUID + NormIdx[v] ] *= In[0];
#  endif

} // FUNCTION : CPU_Pri2Con



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
void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres )
{

   const bool CheckMinPres_Yes = true;
   real Var[NCOMP_FLUID];  // don't need to include passive scalars since they don't have to be rotated
   real Pres, Vx;

   for (int v=0; v<NCOMP_FLUID; v++)   Var[v] = Input[v];

   CPU_Rotate3D( Var, XYZ, true );

   Pres = CPU_GetPressure( Var[0], Var[1], Var[2], Var[3], Var[4], Gamma_m1, CheckMinPres_Yes, MinPres );
   Vx   = Var[1] / Var[0];

   Flux[0] = Var[1];
   Flux[1] = Vx*Var[1] + Pres;
   Flux[2] = Vx*Var[2];
   Flux[3] = Vx*Var[3];
   Flux[4] = Vx*( Var[4] + Pres );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux[v] = Input[v]*Vx;
#  endif

   CPU_Rotate3D( Flux, XYZ, false );

} // FUNCTION : CPU_Con2Flux



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
real CPU_CheckMinPres( const real InPres, const real MinPres )
{

   return FMAX( InPres, MinPres );

} // FUNCTION : CPU_CheckMinPres



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
real CPU_CheckMinPresInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                             const real Gamma_m1, const real _Gamma_m1, const real MinPres )
{

   real InPres, OutPres, Ek, _Dens;

// we didn't use CPU_GetPressure() here to void calculating kinematic energy (Ek) twice
   _Dens   = (real)1.0 / Dens;
   Ek      = (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) * _Dens;
   InPres  = Gamma_m1*( Engy - Ek );
   OutPres = CPU_CheckMinPres( InPres, MinPres );

// do not modify energy (even the round-off errors) if the input pressure passes the check of CPU_CheckMinPres()
   if ( InPres == OutPres )   return Engy;
   else                       return Ek + _Gamma_m1*OutPres;

} // FUNCTION : CPU_CheckMinPresInEngy



#ifdef CHECK_NEGATIVE_IN_FLUID
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
bool CPU_CheckNegative( const real Input )
{

   if ( Input <= (real)0.0  ||  Input >= __FLT_MAX__  ||  Input != Input )    return true;
   else                                                                       return false;

} // FUNCTION : CPU_CheckNegative
#endif // #ifdef CHECK_NEGATIVE_IN_FLUID



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
real CPU_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                      const real Gamma_m1, const bool CheckMinPres, const real MinPres )
{

   real _Dens, Pres;

  _Dens = (real)1.0 / Dens;
   Pres = Gamma_m1*(  Engy - (real)0.5*_Dens*( SQR(MomX) + SQR(MomY) + SQR(MomZ) )  );

   if ( CheckMinPres )   Pres = CPU_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : CPU_GetPressure



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
real CPU_GetTemperature( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                         const real Gamma_m1, const bool CheckMinPres, const real MinPres )
{

   return CPU_GetPressure( Dens, MomX, MomY, MomZ, Engy, Gamma_m1, CheckMinPres, MinPres ) / Dens;

} // FUNCTION : CPU_GetTemperature



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
//                4. Currently this function always adopts double precision since
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
double CPU_Temperature2Pressure( const double Dens, const double Temp, const double mu, const double m_H,
                                 const bool CheckMinPres, const double MinPres )
{

   double Pres;

   Pres = Dens*Temp/(mu*m_H);

   if ( CheckMinPres )  Pres = CPU_CheckMinPres( (real)Pres, (real)MinPres );

   return Pres;

} // FUNCTION : CPU_GetTemperature



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
void CPU_NormalizePassive( const real GasDens, real Passive[], const int NNorm, const int NormIdx[] )
{

// validate the target variable indices
#  ifdef GAMER_DEBUG
   const int MinIdx = 0;
#  ifdef DUAL_ENERGY
   const int MaxIdx = NCOMP_PASSIVE - 2;
#  else
   const int MaxIdx = NCOMP_PASSIVE - 1;
#  endif

   for (int v=0; v<NNorm; v++)
   {
      if ( NormIdx[v] < MinIdx  ||  NormIdx[v] > MaxIdx )
         Aux_Error( ERROR_INFO, "NormIdx[%d] = %d is not within the correct range ([%d <= idx <= %d]) !!\n",
                    v, NormIdx[v], MinIdx, MaxIdx );
   }
#  endif // #ifdef GAMER_DEBUG


   real Norm, PassiveDens_Sum=(real)0.0;

   for (int v=0; v<NNorm; v++)   PassiveDens_Sum += Passive[ NormIdx[v] ];

   Norm = GasDens / PassiveDens_Sum;

   for (int v=0; v<NNorm; v++)   Passive[ NormIdx[v] ] *= Norm;

} // FUNCTION : CPU_NormalizePassive



#endif // #if ( MODEL == HYDRO )
