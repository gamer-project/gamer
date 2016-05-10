#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

/*
#if (  !defined GPU  &&  MODEL == HYDRO  &&  \
       ( FLU_SCHEME == WAF || FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
*/

#if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
real CPU_PositivePres( const real Pres_In, const real Dens, const real _Dens );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Rotate3D
// Description :  Rotate the input 5-element fluid variables properly to simplify the 3D calculation
//
// Note        :  x : (0,1,2,3,4) <--> (0,1,2,3,4)   
//                y : (0,1,2,3,4) <--> (0,2,3,1,4)
//                z : (0,1,2,3,4) <--> (0,3,1,2,4)
//
// Parameter   :  InOut    : Array storing both the input and output data
//                XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)
//                Forward  : (true/false) <--> (forward/backward)
//-------------------------------------------------------------------------------------------------------
void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward )
{
   
   if ( XYZ == 0 )   return;


   real Temp[5];
   for (int v=0; v<5; v++)    Temp[v] = InOut[v];

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : InOut[1] = Temp[2];  InOut[2] = Temp[3];  InOut[3] = Temp[1];     break;
         case 2 : InOut[1] = Temp[3];  InOut[2] = Temp[1];  InOut[3] = Temp[2];     break;
      }
   }

   else // backward
   {
      switch ( XYZ )
      {
         case 1 : InOut[1] = Temp[3];  InOut[2] = Temp[1];  InOut[3] = Temp[2];     break;
         case 2 : InOut[1] = Temp[2];  InOut[2] = Temp[3];  InOut[3] = Temp[1];     break;
      }
   }

} // FUNCTION : CPU_Rotate3D



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Pri 
// Description :  Convert the conserved variables to the primitive variables 
//
// Parameter   :  In       : Array storing the input conserved variables
//                Out      : Array to store the output primitive variables
//                Gamma_m1 : Gamma - 1
//-------------------------------------------------------------------------------------------------------
void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1 )
{

   const real _Rho = (real)1.0 / In[0];
   
   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = (  In[4] - (real)0.5*Out[0]*( Out[1]*Out[1] + Out[2]*Out[2] + Out[3]*Out[3] )  )*Gamma_m1;

// ensure the positive pressure
#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   Out[4] = CPU_PositivePres( Out[4], Out[0], _Rho );
#  endif

} // FUNCTION : CPU_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Pri2Con 
// Description :  Convert the primitive variables to the conserved variables 
//
// Parameter   :  In       : Array storing the input primitive variables
//                Out      : Array to store the output conserved variables
//                Gamma_m1 : 1 / (Gamma - 1)
//-------------------------------------------------------------------------------------------------------
void CPU_Pri2Con( const real In[], real Out[], const real _Gamma_m1 )
{

   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];
   Out[4] = In[4]*_Gamma_m1 + (real)0.5*In[0]*( In[1]*In[1] + In[2]*In[2] + In[3]*In[3] );

} // FUNCTION : CPU_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Flux
// Description :  Evaluate the hydrodynamic fluxes by the input conserved variables
//
// Parameter   :  XYZ   : Targeted spatial direction : (0/1/2) --> (x/y/z) 
//                Flux  : Array to store the output fluxes
//                Input : Array storing the input conserved variables
//                Gamma : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma )
{

   real Var[5];
   real Pres, _Rho, Vx;

   for (int v=0; v<5; v++)    Var[v] = Input[v];

   CPU_Rotate3D( Var, XYZ, true );

   _Rho = (real)1.0 / Var[0];
   Pres = (Gamma-(real)1.0) * (  Var[4] - (real)0.5*( Var[1]*Var[1] + Var[2]*Var[2] + Var[3]*Var[3] )*_Rho  );
   Vx   = _Rho*Var[1];

   Flux[0] = Var[1];
   Flux[1] = Vx*Var[1] + Pres;
   Flux[2] = Vx*Var[2];
   Flux[3] = Vx*Var[3];
   Flux[4] = Vx*( Var[4] + Pres );

   CPU_Rotate3D( Flux, XYZ, false );

} // FUNCTION : CPU_Con2Flux



#if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PositivePres
// Description :  Ensure the positive pressure
//
// Note        :  Please set either MIN_PRES_DENS or MIN_PRES in the header CUFLU.h 
//
// Parameter   :  Pres_In  : Input pressure to be corrected
//                Dens     : Density used by MIN_PRES_DENS 
//                _Dens    : 1.0/Dens
//
// Return      :  Minimum pressure
//-------------------------------------------------------------------------------------------------------
real CPU_PositivePres( const real Pres_In, const real Dens, const real _Dens )
{

   real Pres = Pres_In;

#  ifdef MIN_PRES_DENS
   if ( Pres*_Dens < MIN_PRES_DENS )   Pres = Dens*MIN_PRES_DENS;
#  elif defined MIN_PRES
   if ( Pres       < MIN_PRES      )   Pres = MIN_PRES;
#  else
#  error : ERROR : neither MIN_PRES_DENS nor MIN_PRES are defined !!
#  endif

   return Pres;

} // FUNCTION : CPU_PositivePres



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PositivePres_In_Engy
// Description :  Ensure the positive pressure in the input total energy
//
// Note        :  Please set either MIN_PRES_DENS or MIN_PRES in the header CUFLU.h 
//
// Parameter   :  ConVar      : Conserved variable to be corrected
//                Gamma_m1    : Gamma - 1
//                _Gamma_m1   : 1/(Gamma - 1)
//
// Return      :  Total energy with the minimum pressure
//-------------------------------------------------------------------------------------------------------
real CPU_PositivePres_In_Engy( const real ConVar[], const real Gamma_m1, const real _Gamma_m1 )
{

   real TempPres, Ek, _Dens, Engy;

   _Dens    = (real)1.0 / ConVar[0];
   Ek       = (real)0.5*( SQR(ConVar[1]) + SQR(ConVar[2]) + SQR(ConVar[3]) ) * _Dens;
   TempPres = Gamma_m1*( ConVar[4] - Ek );
   TempPres = CPU_PositivePres( TempPres, ConVar[0], _Dens );
   Engy     = Ek + _Gamma_m1*TempPres;

   return Engy;

} // FUNCTION : CPU_PositivePres_In_Engy
#endif // #if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )



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
#endif



//#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == WAF || MHM || MHM_RP || CTU) )
