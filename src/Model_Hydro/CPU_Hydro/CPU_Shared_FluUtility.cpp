#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

// some functions in this file need to be defined even when using GPU
#if ( MODEL == HYDRO )

real CPU_CheckMinPres( const real InPres, const real MinPres );
real CPU_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                      const real Gamma_m1, const bool CheckMinPres, const real MinPres );




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
//                XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)
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
//
// Parameter   :  In       : Array storing the input conserved variables
//                Out      : Array to store the output primitive variables
//                Gamma_m1 : Gamma - 1
//                MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres )
{

   const bool CheckMinPres_Yes = true;
   const real _Rho             = (real)1.0 / In[0];

   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = CPU_GetPressure( In[0], In[1], In[2], In[3], In[4], Gamma_m1, CheckMinPres_Yes, MinPres );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v]*_Rho;
#  endif

} // FUNCTION : CPU_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Pri2Con
// Description :  Convert the primitive variables to the conserved variables
//
// Note        :  1. This function does NOT check if the input pressure is greater than the
//                   given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//
// Parameter   :  In       : Array storing the input primitive variables
//                Out      : Array to store the output conserved variables
//               _Gamma_m1 : 1 / (Gamma - 1)
//-------------------------------------------------------------------------------------------------------
void CPU_Pri2Con( const real In[], real Out[], const real _Gamma_m1 )
{

   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];
   Out[4] = In[4]*_Gamma_m1 + (real)0.5*In[0]*( In[1]*In[1] + In[2]*In[2] + In[3]*In[3] );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[0]*In[v];
#  endif

} // FUNCTION : CPU_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Flux
// Description :  Evaluate the hydrodynamic fluxes by the input conserved variables
//
// Parameter   :  XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)
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



#endif // #if ( MODEL == HYDRO )
