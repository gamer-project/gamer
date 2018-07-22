#ifndef __CUFLU_FLUUTILITY_CU__
#define __CUFLU_FLUUTILITY_CU__



#include "Macro.h"
#include "CUFLU.h"

static __device__ FluVar CUFLU_Pri2Con( const FluVar Pri, const real _Gamma_m1,
                                        const bool NormPassive, const int NNorm, const int NormIdx[] );
static __device__ FluVar CUFLU_Con2Pri( const FluVar Con, const real Gamma_m1, const real MinPres,
                                        const bool NormPassive, const int NNorm, const int NormIdx[],
                                        const bool JeansMinPres, const real JeansMinPres_Coeff );
static __device__ FluVar CUFLU_Con2Flux( const FluVar Input, const real Gamma_m1, const int XYZ, const real MinPres );
static __device__ FluVar CUFLU_Rotate3D( const FluVar In, const int XYZ, const bool Forward );
static __device__ void CUFLU_Con2Pri_AllGrids( const real g_Fluid_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                               real g_PriVar[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                               const real Gamma, const real MinPres,
                                               const bool NormPassive, const int NNorm, const int NormIdx[],
                                               const bool JeansMinPres, const real JeansMinPres_Coeff );
static __device__ real CUFLU_CheckMinPres( const real InPres, const real MinPres );
static __device__ real CUFLU_CheckMinPresInEngy( const FluVar ConVar, const real Gamma_m1, const real _Gamma_m1, const real MinPres );
#ifdef CHECK_NEGATIVE_IN_FLUID
static __device__ bool CUFLU_CheckNegative( const real Input );
#endif
static __device__ real CUFLU_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                                          const real Gamma_m1, const bool CheckMinPres, const real MinPres );
#if ( NCOMP_PASSIVE > 0 )
static __device__ void CUFLU_NormalizePassive( const real GasDens, real Passive[], const int NNorm, const int NormIdx[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Rotate3D
// Description :  Rotate the input fluid variables properly to simplify the 3D calculation
//
// Note        :  1. x : (0,1,2,3,4) <--> (0,1,2,3,4)
//                   y : (0,1,2,3,4) <--> (0,2,3,1,4)
//                   z : (0,1,2,3,4) <--> (0,3,1,2,4)
//                2. Passive scalars are not modified at all
//
// Parameter   :  In       : Input variables to be rotated
//                XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Forward  : (true/false) <--> (forward/backward)
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Rotate3D( const FluVar In, const int XYZ, const bool Forward )
{

   if ( XYZ == 0 )   return In;


   FluVar Out;

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : Out.Rho=In.Rho;  Out.Px=In.Py;  Out.Py=In.Pz;  Out.Pz=In.Px;  Out.Egy=In.Egy;    break;
         case 2 : Out.Rho=In.Rho;  Out.Px=In.Pz;  Out.Py=In.Px;  Out.Pz=In.Py;  Out.Egy=In.Egy;    break;
      }
   }

   else
   {
      switch ( XYZ )
      {
         case 1 : Out.Rho=In.Rho;  Out.Px=In.Pz;  Out.Py=In.Px;  Out.Pz=In.Py;  Out.Egy=In.Egy;    break;
         case 2 : Out.Rho=In.Rho;  Out.Px=In.Py;  Out.Py=In.Pz;  Out.Pz=In.Px;  Out.Egy=In.Egy;    break;
      }
   }

#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)    Out.Passive[v] = In.Passive[v];
#  endif

   return Out;

} // FUNCTION : CUFLU_Rotate3D



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Con2Flux
// Description :  Conserved variables --> fluxes
//
// Parameter   :  Input    : Input conserved variables
//                Gamma_m1 : Gamma - 1
//                XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Con2Flux( const FluVar Input, const real Gamma_m1, const int XYZ, const real MinPres )
{

   FluVar Temp;
   real   Pres, _Rho, Vx;

   Temp = CUFLU_Rotate3D( Input, XYZ, true );

   _Rho = (real)1.0 / Temp.Rho;
   Pres = Gamma_m1 * (  Temp.Egy - (real)0.5*( Temp.Px*Temp.Px + Temp.Py*Temp.Py + Temp.Pz*Temp.Pz )*_Rho  );
   Pres = CUFLU_CheckMinPres( Pres, MinPres );
   Vx   = _Rho*Temp.Px;

   Temp.Rho = Temp.Px;
   Temp.Px  = Vx*Temp.Px + Pres;
   Temp.Py  = Vx*Temp.Py;
   Temp.Pz  = Vx*Temp.Pz;
   Temp.Egy = Vx*( Temp.Egy + Pres );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)    Temp.Passive[v] *= Vx;
#  endif

   Temp = CUFLU_Rotate3D( Temp, XYZ, false );

   return Temp;

} // FUNCTION : CUFLU_Con2Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Pri2Con
// Description :  Primitive variables --> conserved variables
//
// Note        :  1. This function does NOT check if the input pressure is greater than the
//                   given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//
// Parameter   :  Pri         : Input primitive variables
//                _Gamma_m1   : 1.0 / (Gamma-1.0)
//                NormPassive : true --> convert passive scalars to mass fraction
//                NNorm       : Number of passive scalars for the option "NormPassive"
//                              --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx     : Target variable indices for the option "NormPassive"
//                              --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Pri2Con( const FluVar Pri, const real _Gamma_m1,
                                 const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   FluVar Con;

   Con.Rho = Pri.Rho;
   Con.Egy = _Gamma_m1*Pri.Egy + (real)0.5*Pri.Rho*( Pri.Px*Pri.Px + Pri.Py*Pri.Py + Pri.Pz*Pri.Pz );
   Con.Px  = Pri.Rho*Pri.Px;
   Con.Py  = Pri.Rho*Pri.Py;
   Con.Pz  = Pri.Rho*Pri.Pz;

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=0; v<NCOMP_PASSIVE; v++)    Con.Passive[v] = Pri.Passive[v];

// convert the mass fraction of target passive scalars back to mass density
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Con.Passive[ NormIdx[v] ] *= Pri.Rho;
#  endif

   return Con;

} // FUNCTION : CUFLU_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Con2Pri
// Description :  Conserved variables --> primitive variables
//
// Note        :  1. This function always check if the pressure to be returned is greater than the
//                   given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//                   --> But note that here we do NOT ensure "sum(mass fraction) == 1.0"
//                       --> It is done by calling CUFLU_NormalizePassive() in CUFLU_Shared_FullStepUpdate()
//
// Parameter   :  Con                : Input conserved variables
//                Gamma_m1           : Gamma - 1.0
//                MinPres            : Minimum allowed pressure
//                NormPassive        : true --> convert passive scalars to mass fraction
//                NNorm              : Number of passive scalars for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Con2Pri( const FluVar Con, const real Gamma_m1, const real MinPres,
                                 const bool NormPassive, const int NNorm, const int NormIdx[],
                                 const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const real _Rho = (real)1.0/Con.Rho;
   FluVar Pri;

   Pri.Rho = Con.Rho;
   Pri.Px  = Con.Px*_Rho;
   Pri.Py  = Con.Py*_Rho;
   Pri.Pz  = Con.Pz*_Rho;
   Pri.Egy = Gamma_m1*(  Con.Egy - (real)0.5*( Con.Px*Con.Px + Con.Py*Con.Py + Con.Pz*Con.Pz )*_Rho  );
   Pri.Egy = CUFLU_CheckMinPres( Pri.Egy, MinPres );

// pressure floor required to resolve the Jeans length
// --> note that currently we do not modify the dual-energy variable (e.g., entropy) accordingly
   if ( JeansMinPres )
   Pri.Egy = CUFLU_CheckMinPres( Pri.Egy, JeansMinPres_Coeff*SQR(Pri.Rho) );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=0; v<NCOMP_PASSIVE; v++)    Pri.Passive[v] = Con.Passive[v];

// convert the mass density of target passive scalars to mass fraction
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Pri.Passive[ NormIdx[v] ] *= _Rho;
#  endif

   return Pri;

} // FUNCTION : CUFLU_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Con2Pri_AllGrids
// Description :  Conserved variables --> primitive variables for all grids
//
// Parameter   :  g_Fluid_In         : Global memory array storing the input fluid variables
//                g_PriVar           : Global memory array to store the output primitive variables
//                Gamma              : Ratio of specific heats
//                MinPres            : Minimum allowed pressure
//                NormPassive        : true --> convert passive scalars to mass fraction
//                NNorm              : Number of passive scalars for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_Con2Pri_AllGrids( const real g_Fluid_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                        real g_PriVar[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                        const real Gamma, const real MinPres,
                                        const bool NormPassive, const int NNorm, const int NormIdx[],
                                        const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const uint bx       = blockIdx.x;
   const uint tx       = threadIdx.x;
   const uint dID      = blockDim.x;
   const real Gamma_m1 = Gamma - (real)1.0;

   uint   ID = tx;
   FluVar Var;

// loop over all cells
   while ( ID < FLU_NXT*FLU_NXT*FLU_NXT )
   {
      Var.Rho = g_Fluid_In[bx][0][ID];
      Var.Px  = g_Fluid_In[bx][1][ID];
      Var.Py  = g_Fluid_In[bx][2][ID];
      Var.Pz  = g_Fluid_In[bx][3][ID];
      Var.Egy = g_Fluid_In[bx][4][ID];

#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)    Var.Passive[v] = g_Fluid_In[bx][ NCOMP_FLUID + v ][ID];
#     endif

      Var = CUFLU_Con2Pri( Var, Gamma_m1, MinPres, NormPassive, NNorm, NormIdx, JeansMinPres, JeansMinPres_Coeff );

      g_PriVar[bx][0][ID] = Var.Rho;
      g_PriVar[bx][1][ID] = Var.Px;
      g_PriVar[bx][2][ID] = Var.Py;
      g_PriVar[bx][3][ID] = Var.Pz;
      g_PriVar[bx][4][ID] = Var.Egy;

#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)    g_PriVar[bx][ NCOMP_FLUID + v ][ID] = Var.Passive[v];
#     endif

      ID += dID;
   } // while ( ID < FLU_NXT*FLU_NXT*FLU_NXT )

} // FUNCTION : CUFLU_Con2Pri_AllGrids



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_CheckMinPres
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
__device__ real CUFLU_CheckMinPres( const real InPres, const real MinPres )
{

   return FMAX( InPres, MinPres );

} // FUNCTION : CUFLU_CheckMinPres



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_CheckMinPresInEngy
// Description :  Ensure that the pressure in the input total energy is greater than the given threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
//                   --> Currently it simply sets a minimum allowed value for pressure
//                       --> Please set MIN_PRES in the runtime parameter file "Input__Parameter"
//                3. One must input conserved variables instead of primitive variables
//
// Parameter   :  ConVar   : Conserved variable to be corrected
//                Gamma_m1 : Gamma - 1
//               _Gamma_m1 : 1/(Gamma - 1)
//                MinPres  : Minimum allowed pressure
//
// Return      :  Total energy with pressure greater than the given threshold
//-------------------------------------------------------------------------------------------------------
__device__ real CUFLU_CheckMinPresInEngy( const FluVar ConVar, const real Gamma_m1, const real _Gamma_m1, const real MinPres )
{

   real InPres, OutPres, Ek, _Dens;

   _Dens   = (real)1.0 / ConVar.Rho;
   Ek      = (real)0.5*( SQR(ConVar.Px) + SQR(ConVar.Py) + SQR(ConVar.Pz) ) * _Dens;
   InPres  = Gamma_m1*( ConVar.Egy - Ek );
   OutPres = CUFLU_CheckMinPres( InPres, MinPres );

// do not modify energy (even the round-off errors) if the input pressure passes the check of CPU_CheckMinPres()
   if ( InPres == OutPres )   return ConVar.Egy;
   else                       return Ek + _Gamma_m1*OutPres;

} // FUNCTION : CUFLU_CheckMinPresInEngy



#ifdef CHECK_NEGATIVE_IN_FLUID
//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_CheckNegative
// Description :  Check whether the input value is <= 0.0 (also check whether it's Inf or NAN)
//
// Note        :  Can be used to check whether the values of density and pressure are unphysical
//
// Parameter   :  Input : Input value
//
// Return      :  true  --> Input <= 0.0  ||  >= __FLT_MAX__  ||  != itself (Nan)
//                false --> otherwise
//-------------------------------------------------------------------------------------------------------
__device__ bool CUFLU_CheckNegative( const real Input )
{

   if ( Input <= (real)0.0  ||  Input >= __FLT_MAX__  ||  Input != Input )    return true;
   else                                                                       return false;

} // FUNCTION : CUFLU_CheckNegative
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_GetPressure
// Description :  Evaluate the fluid pressure
//
// Note        :  1. Currently only work with the adiabatic EOS
//                2. One must input conserved variables instead of primitive variables
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density
//                Gamma_m1     : Gamma - 1, where Gamma is the adiabatic index
//                CheckMinPres : Return CUFLU_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Minimum allowed pressure
//
// Return      :  Pressure
//-------------------------------------------------------------------------------------------------------
__device__ real CUFLU_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                                   const real Gamma_m1, const bool CheckMinPres, const real MinPres )
{

   real _Dens, Pres;

  _Dens = (real)1.0 / Dens;
   Pres = Gamma_m1*(  Engy - (real)0.5*_Dens*( SQR(MomX) + SQR(MomY) + SQR(MomZ) )  );

   if ( CheckMinPres )   Pres = CUFLU_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : CUFLU_GetPressure



#if ( NCOMP_PASSIVE > 0 )
//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_NormalizePassive
// Description :  Normalize the target passive scalars so that the sum of their mass density is equal to
//                the gas mass density
//
// Note        :  1. Should be invoked AFTER applying the floor values to passive scalars
//                2. Invoked by CUFLU_Shared_FullStepUpdate()
//                3. The CPU version is defined in CPU_Shared_FluUtility->CPU_NormalizePassive()
//
// Parameter   :  GasDens  : Gas mass density
//                Passive  : Passive scalar array (with the size NCOMP_PASSIVE)
//                NNorm    : Number of passive scalars to be normalized
//                           --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx  : Target variable indices to be normalized
//                           --> Should be set to the global variable "PassiveNorm_VarIdx"
//
// Return      :  Passive
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_NormalizePassive( const real GasDens, real Passive[], const int NNorm, const int NormIdx[] )
{

// validate the target variable indices
#  ifdef GAMER_DEBUG
   const int MinIdx = 0;
   const int MaxIdx = NCOMP_PASSIVE - 1;

   for (int v=0; v<NNorm; v++)
   {
      if ( NormIdx[v] < MinIdx  ||  NormIdx[v] > MaxIdx )
         printf( "ERROR : NormIdx[%d] = %d is not within the correct range ([%d <= idx <= %d]) !!\n",
                 v, NormIdx[v], MinIdx, MaxIdx );
   }
#  endif // #ifdef GAMER_DEBUG


   real Norm, PassiveDens_Sum=(real)0.0;

   for (int v=0; v<NNorm; v++)   PassiveDens_Sum += Passive[ NormIdx[v] ];

   Norm = GasDens / PassiveDens_Sum;

   for (int v=0; v<NNorm; v++)   Passive[ NormIdx[v] ] *= Norm;

} // FUNCTION : CUFLU_NormalizePassive
#endif // #if ( NCOMP_PASSIVE > 0 )



#endif // #ifndef __CUFLU_FLUUTILITY_CU__
