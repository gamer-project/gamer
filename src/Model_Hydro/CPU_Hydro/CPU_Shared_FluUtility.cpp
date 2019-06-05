#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// internal function prototypes
// --> only necessary for GPU since they are included in Prototype.h for the CPU codes
#ifdef __CUDACC__
GPU_DEVICE
static real Hydro_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                               const real Gamma_m1, const bool CheckMinPres, const real MinPres, const real EngyB );
GPU_DEVICE
static real Hydro_CheckMinPres( const real InPres, const real MinPres );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Rotate3D
// Description :  Rotate the input fluid variables properly to simplify the 3D calculation
//
// Note        :  1. x : (x,y,z) <--> (x,y,z)
//                   y : (x,y,z) <--> (y,z,x)
//                   z : (x,y,z) <--> (z,x,y)
//                2. Work no matter InOut[] includes passive scalars or not since they are not modified at all
//                   --> For MHD, specify the array offset of magnetic field by Mag_Offset
//
// Parameter   :  InOut      : Array storing both the input and output data
//                XYZ        : Target spatial direction : (0/1/2) --> (x/y/z)
//                Forward    : (true/false) <--> (forward/backward)
//                Mag_Offset : Array offset of magnetic field (for MHD only)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset )
{

   if ( XYZ == 0 )   return;


// check
#  ifdef GAMER_DEBUG
#  ifdef MHD
   if ( Mag_Offset < NCOMP_FLUID  ||  Mag_Offset > NCOMP_TOTAL_PLUS_MAG - NCOMP_MAG )
      printf( "ERROR : invalid Mag_Offset = %d !!\n", Mag_Offset );
#  endif
#  endif


   real Temp_Flu[3];
   for (int v=0; v<3; v++)    Temp_Flu[v] = InOut[ v + 1 ];
#  ifdef MHD
   real Temp_Mag[3];
   for (int v=0; v<3; v++)    Temp_Mag[v] = InOut[ v + Mag_Offset ];
#  endif

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : InOut[              1 ] = Temp_Flu[1];
                  InOut[              2 ] = Temp_Flu[2];
                  InOut[              3 ] = Temp_Flu[0];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[1];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[0];
#                 endif
                  break;

         case 2 : InOut[              1 ] = Temp_Flu[2];
                  InOut[              2 ] = Temp_Flu[0];
                  InOut[              3 ] = Temp_Flu[1];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[0];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[1];
#                 endif
                  break;
      }
   }

   else // backward
   {
      switch ( XYZ )
      {
         case 1 : InOut[              1 ] = Temp_Flu[2];
                  InOut[              2 ] = Temp_Flu[0];
                  InOut[              3 ] = Temp_Flu[1];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[0];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[1];
#                 endif
                  break;

         case 2 : InOut[              1 ] = Temp_Flu[1];
                  InOut[              2 ] = Temp_Flu[2];
                  InOut[              3 ] = Temp_Flu[0];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[1];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[0];
#                 endif
                  break;
      }
   }

} // FUNCTION : Hydro_Rotate3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Pri
// Description :  Conserved variables --> primitive variables
//
// Note        :  1. This function always check if the pressure to be returned is greater than the
//                   given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//                   --> But note that here we do NOT ensure "sum(mass fraction) == 1.0"
//                       --> It is done by calling Hydro_NormalizePassive() in Hydro_Shared_FullStepUpdate()
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//
// Parameter   :  In                 : Input conserved variables
//                Out                : Output primitive variables
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
GPU_DEVICE
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const bool CheckMinPres_Yes = true;
   const real _Rho             = (real)1.0/In[0];
#  ifdef MHD
   const real Bx               = In[ MAG_OFFSET + 0 ];
   const real By               = In[ MAG_OFFSET + 1 ];
   const real Bz               = In[ MAG_OFFSET + 2 ];
   const real EngyB            = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  else
   const real EngyB            = NULL_REAL;
#  endif

   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = Hydro_GetPressure( In[0], In[1], In[2], In[3], In[4], Gamma_m1, CheckMinPres_Yes, MinPres, EngyB );

// pressure floor required to resolve the Jeans length
// --> note that currently we do not modify the dual-energy variable (e.g., entropy) accordingly
   if ( JeansMinPres )
   Out[4] = Hydro_CheckMinPres( Out[4], JeansMinPres_Coeff*SQR(Out[0]) );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass density of target passive scalars to mass fraction
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Out[ NCOMP_FLUID + NormIdx[v] ] *= _Rho;
#  endif

// B field
#  ifdef MHD
   for (int v=NCOMP_TOTAL; v<NCOMP_TOTAL_PLUS_MAG; v++)  Out[v] = In[v];
#  endif

} // FUNCTION : Hydro_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Pri2Con
// Description :  Primitive variables --> conserved variables
//
// Note        :  1. Does NOT check if the input pressure is greater than the given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
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
GPU_DEVICE
void Hydro_Pri2Con( const real In[], real Out[], const real _Gamma_m1,
                    const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];
   Out[4] = In[4]*_Gamma_m1 + (real)0.5*In[0]*( In[1]*In[1] + In[2]*In[2] + In[3]*In[3] );

#  ifdef MHD
   const real Bx = In[ MAG_OFFSET + 0 ];
   const real By = In[ MAG_OFFSET + 1 ];
   const real Bz = In[ MAG_OFFSET + 2 ];
   Out[4] += (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  endif

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass fraction of target passive scalars back to mass density
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Out[ NCOMP_FLUID + NormIdx[v] ] *= In[0];
#  endif

// B field
#  ifdef MHD
   for (int v=NCOMP_TOTAL; v<NCOMP_TOTAL_PLUS_MAG; v++)  Out[v] = In[v];
#  endif

} // FUNCTION : Hydro_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Flux
// Description :  Evaluate hydrodynamic/MHD fluxes from the input conserved variables
//
// Note        :  1. Flux[] and In[] may point to the same array
//                2. Flux[] and In[] should have the size of NCOMP_TOTAL_PLUS_MAG
//
// Parameter   :  XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux     : Array to store the output fluxes
//                In       : Array storing the input conserved variables
//                Gamma_m1 : Gamma - 1
//                MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real Gamma_m1, const real MinPres )
{

   const bool CheckMinPres_Yes = true;
   real InRot[ NCOMP_FLUID + NCOMP_MAG ];    // no need to include passive scalars since they don't have to be rotated

   for (int v=0; v<NCOMP_FLUID; v++)   InRot[v] = In[v];

#  ifdef MHD
   for (int v=NCOMP_FLUID; v<NCOMP_FLUID+NCOMP_MAG; v++)    InRot[v] = In[ v - NCOMP_FLUID + MAG_OFFSET ];
#  endif

   Hydro_Rotate3D( InRot, XYZ, true, NCOMP_FLUID );

#  ifdef MHD
   const real Bx    = InRot[ NCOMP_FLUID + 0 ];
   const real By    = InRot[ NCOMP_FLUID + 1 ];
   const real Bz    = InRot[ NCOMP_FLUID + 2 ];
   const real EngyB = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  else
   const real EngyB = NULL_REAL;
#  endif
   const real Pres  = Hydro_GetPressure( InRot[0], InRot[1], InRot[2], InRot[3], InRot[4],
                                         Gamma_m1, CheckMinPres_Yes, MinPres, EngyB );
   const real _Rho  = (real)1.0 / InRot[0];
   const real Vx    = _Rho*InRot[1];

   Flux[0] = InRot[1];
   Flux[1] = Vx*InRot[1] + Pres;
   Flux[2] = Vx*InRot[2];
   Flux[3] = Vx*InRot[3];
   Flux[4] = Vx*( InRot[4] + Pres );

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux[v] = In[v]*Vx;
#  endif

// B field
#  ifdef MHD
   const real Vy = _Rho*InRot[2];
   const real Vz = _Rho*InRot[3];

   Flux[              1 ] += EngyB - SQR(Bx);
   Flux[              2 ] -= Bx*By;
   Flux[              3 ] -= Bx*Bz;
   Flux[              4 ] += Vx*EngyB - Bx*( Bx*Vx + By*Vy + Bz*Vz );
   Flux[ MAG_OFFSET + 0 ]  = (real)0.0;
   Flux[ MAG_OFFSET + 1 ]  = By*Vx - Bx*Vy;
   Flux[ MAG_OFFSET + 2 ]  = Bz*Vx - Bx*Vz;
#  endif

   Hydro_Rotate3D( Flux, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_Con2Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinPres
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
GPU_DEVICE
real Hydro_CheckMinPres( const real InPres, const real MinPres )
{

   return FMAX( InPres, MinPres );

} // FUNCTION : Hydro_CheckMinPres



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinPresInEngy
// Description :  Ensure that the pressure in the input total energy is greater than the given threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
//                   --> Currently it simply sets a minimum allowed value for pressure
//                       --> Please set MIN_PRES in the runtime parameter file "Input__Parameter"
//                3. One must input conserved variables instead of primitive variables
//                4. For MHD, one must provide the magnetic energy density EngyB (i.e., 0.5*B^2)
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                Engy     : Energy density
//                Gamma_m1 : Gamma - 1
//               _Gamma_m1 : 1/(Gamma - 1)
//                MinPres  : Minimum allowed pressure
//                EngyB    : Magnetic energy density (0.5*B^2)
//                           --> For MHD only
//
// Return      :  Total energy with pressure greater than the given threshold
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinPresInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                               const real Gamma_m1, const real _Gamma_m1, const real MinPres, const real EngyB )
{

   real InPres, OutPres, Ek, _Dens;

// we didn't use Hydro_GetPressure() here to avoid calculating kinematic energy (Ek) twice
   _Dens   = (real)1.0 / Dens;
   Ek      = (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) * _Dens;
#  ifdef MHD
   InPres  = Gamma_m1*( Engy - Ek - EngyB );
#  else
   InPres  = Gamma_m1*( Engy - Ek );
#  endif
   OutPres = Hydro_CheckMinPres( InPres, MinPres );

// do not modify energy (even the round-off errors) if the input pressure passes the check of Hydro_CheckMinPres()
   if ( InPres == OutPres )
      return Engy;

   else
   {
#     ifdef MHD
      return Ek + _Gamma_m1*OutPres + EngyB;
#     else
      return Ek + _Gamma_m1*OutPres;
#     endif
   }

} // FUNCTION : Hydro_CheckMinPresInEngy



#ifdef CHECK_NEGATIVE_IN_FLUID
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckNegative
// Description :  Check whether the input value is <= 0.0 (also check whether it's Inf or NAN)
//
// Note        :  Can be used to check whether the values of density and pressure are unphysical
//
// Parameter   :  Input : Input value
//
// Return      :  true  --> Input <= 0.0  ||  >= __FLT_MAX__  ||  != itself (Nan)
//                false --> otherwise
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool Hydro_CheckNegative( const real Input )
{

   if ( Input <= (real)0.0  ||  Input >= __FLT_MAX__  ||  Input != Input )    return true;
   else                                                                       return false;

} // FUNCTION : Hydro_CheckNegative
#endif // #ifdef CHECK_NEGATIVE_IN_FLUID



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetPressure
// Description :  Evaluate the fluid pressure
//
// Note        :  1. Currently only work with the adiabatic EOS
//                2. Invoked by Hydro_GetTimeStep_Fluid(), Prepare_PatchData(), InterpolateGhostZone(),
//                   Hydro_Aux_Check_Negative() ...
//                3. One must input conserved variables instead of primitive variables
//                4. For MHD, Engy is the total energy density including the magnetic energy EngyB=0.5*B^2,
//                   and thus one must provide EngyB to calculate the gas pressure
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density (including the magnetic energy density for MHD)
//                Gamma_m1     : Gamma - 1, where Gamma is the adiabatic index
//                CheckMinPres : Return Hydro_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                                   --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres      : Minimum allowed pressure
//                EngyB        : Magnetic energy density (0.5*B^2)
//                               --> For MHD only
//
// Return      :  Gas pressure (Pres)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                        const real Gamma_m1, const bool CheckMinPres, const real MinPres, const real EngyB )
{

   real _Dens, Pres;

  _Dens  = (real)1.0 / Dens;
   Pres  = Engy - (real)0.5*_Dens*( SQR(MomX) + SQR(MomY) + SQR(MomZ) );
#  ifdef MHD
   Pres -= EngyB;
#  endif
   Pres *= Gamma_m1;

   if ( CheckMinPres )   Pres = Hydro_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : Hydro_GetPressure



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetTemperature
// Description :  Evaluate the fluid temperature
//
// Note        :  1. Currently only work with the adiabatic EOS
//                2. For simplicity, currently this function only returns **pressure/density**, which does
//                   NOT include normalization
//                   --> For OPT__FLAG_LOHNER_TEMP only
//                   --> Also note that currently it only checks minimum pressure but not minimum density
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density
//                Gamma_m1     : Gamma - 1, where Gamma is the adiabatic index
//                CheckMinPres : Return Hydro_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                                   --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres      : Minimum allowed pressure
//                EngyB        : Magnetic energy density (0.5*B^2)
//                               --> For MHD only
//
// Return      :  Temperature
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_GetTemperature( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                           const real Gamma_m1, const bool CheckMinPres, const real MinPres, const real EngyB )
{

   return Hydro_GetPressure( Dens, MomX, MomY, MomZ, Engy, Gamma_m1, CheckMinPres, MinPres, EngyB ) / Dens;

} // FUNCTION : Hydro_GetTemperature



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Temperature2Pressure
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
//                CheckMinPres : Return Hydro_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Minimum allowed pressure
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
double Hydro_Temperature2Pressure( const double Dens, const double Temp, const double mu, const double m_H,
                                   const bool CheckMinPres, const double MinPres )
{

   double Pres;

   Pres = Dens*Temp/(mu*m_H);

   if ( CheckMinPres )  Pres = Hydro_CheckMinPres( (real)Pres, (real)MinPres );

   return Pres;

} // FUNCTION : Hydro_GetTemperature



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_NormalizePassive
// Description :  Normalize the target passive scalars so that the sum of their mass density is equal to
//                the gas mass density
//
// Note        :  1. Should be invoked AFTER applying the floor values to passive scalars
//                2. Invoked by Hydro_Shared_FullStepUpdate(), Prepare_PatchData(), Refine(), LB_Refine_AllocateNewPatch(),
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
GPU_DEVICE
void Hydro_NormalizePassive( const real GasDens, real Passive[], const int NNorm, const int NormIdx[] )
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
#  endif


   real Norm, PassiveDens_Sum=(real)0.0;

   for (int v=0; v<NNorm; v++)   PassiveDens_Sum += Passive[ NormIdx[v] ];

   Norm = GasDens / PassiveDens_Sum;

   for (int v=0; v<NNorm; v++)   Passive[ NormIdx[v] ] *= Norm;

} // FUNCTION : Hydro_NormalizePassive



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBField
// Description :  Calculate the cell-centered magnetic field from the input face-centered magnetic field array
//
// Note        :  1. Use the central average operator
//                2. Return all three components of the B field
//                3. Input arrays should have the following dimension:
//                      Bx_FC[]: (Nx+1)*(Ny  )*(Nz  )
//                      By_FC[]: (Nx  )*(Ny+1)*(Nz  )
//                      Bz_FC[]: (Nx  )*(Ny  )*(Nz+1)
//
// Parameter   :  B_CC      : Cell-centered B field to be returned
//                Bx/y/z_FC : Input face-centered B field array
//                Nx/y/z    : Array dimension along different directions (see Note above)
//                i/j/k     : Target cell indices
//
// Return      :  B_CC
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_GetCellCenteredBField( real B_CC[], const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                                const int Nx, const int Ny, const int Nz, const int i, const int j, const int k )
{

   const int idx_Bx = IDX321_BX( i, j, k, Nx, Ny );
   const int idx_By = IDX321_BY( i, j, k, Nx, Ny );
   const int idx_Bz = IDX321_BZ( i, j, k, Nx, Ny );

   B_CC[0] = (real)0.5*( Bx_FC[idx_Bx] + Bx_FC[ idx_Bx + 1     ] );
   B_CC[1] = (real)0.5*( By_FC[idx_By] + By_FC[ idx_By + Nx    ] );
   B_CC[2] = (real)0.5*( Bz_FC[idx_Bz] + Bz_FC[ idx_Bz + Nx*Ny ] );

} // FUNCTION : MHD_GetCellCenteredBField



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBEnergy
// Description :  Calculate the cell-centered magnetic energy (i.e., 0.5*B^2) from the input face-centered
//                magnetic field array
//
// Note        :  1. Invoke MHD_GetCellCenteredBField()
//                2. Input arrays should have the following dimension:
//                      Bx_FC[]: (Nx+1)*(Ny  )*(Nz  )
//                      By_FC[]: (Nx  )*(Ny+1)*(Nz  )
//                      Bz_FC[]: (Nx  )*(Ny  )*(Nz+1)
//
// Parameter   :  Bx/y/z_FC : Input face-centered B field array
//                Nx/y/z    : Array dimension along different directions (see Note above)
//                i/j/k     : Target cell indices
//
// Return      :  0.5*B^2 at the center of the target cell
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real MHD_GetCellCenteredBEnergy( const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                                 const int Nx, const int Ny, const int Nz, const int i, const int j, const int k )
{

// CC = cell-centered
   real B_CC[3], BEngy;

   MHD_GetCellCenteredBField( B_CC, Bx_FC, By_FC, Bz_FC, Nx, Ny, Nz, i, j, k );

   BEngy = (real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) );

   return BEngy;

} // FUNCTION : MHD_GetCellCenteredBEnergy
#endif // #ifdef MHD



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_FLUUTILITY__
