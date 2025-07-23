#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// internal function prototypes
// --> only necessary for GPU since they are included in Prototype.h for the CPU codes
#ifdef __CUDACC__
GPU_DEVICE
static real Hydro_Con2Pres( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                            const real Passive[], const bool CheckMinPres, const real MinPres, const long PassiveFloor, const real Emag,
                            const EoS_DE2P_t EoS_DensEint2Pres, const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                            const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                            const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut );
GPU_DEVICE
static real Hydro_Con2Eint( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                            const bool CheckMinEint, const real MinEint, const long PassiveFloor, const real Emag,
                            const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                            const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                            const real *const EoS_Table[EOS_NTABLE_MAX] );
GPU_DEVICE
static real Hydro_ConEint2Etot( const real Dens, const real MomX, const real MomY, const real MomZ, const real Eint,
                                const real Emag );
GPU_DEVICE
static real Hydro_CheckMinPres( const real InPres, const real MinPres );
GPU_DEVICE
static real Hydro_CheckMinEint( const real InEint, const real MinEint );
GPU_DEVICE
static real Hydro_CheckMinTemp( const real InTemp, const real MinTemp );
GPU_DEVICE
static real Hydro_CheckMinEntr( const real InEntr, const real MinEntr );
GPU_DEVICE
static bool Hydro_IsUnphysical( const IsUnphyMode_t Mode, const real Fields[],
                                const real Emag, const EoS_DE2P_t EoS_DensEint2Pres,
                                const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                                const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                                const real *const EoS_Table[EOS_NTABLE_MAX], const long PassiveFloor,
                                const char File[], const int Line, const char Function[], const IsUnphVerb_t Verbose );
GPU_DEVICE
static bool Hydro_IsUnphysical_Single( const real Field, const char SingleFieldName[], const real Min, const real Max,
                                       const char File[], const int Line, const char Function[], const IsUnphVerb_t Verbose );
#ifdef SRHD
GPU_DEVICE
static real Hydro_Con2HTilde( const real Con[], const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                              const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                              const real *const EoS_Table[EOS_NTABLE_MAX] );
#endif
#endif // #ifdef __CUDACC__

GPU_DEVICE
void NewtonRaphsonSolver( void (*FuncPtr)( real, void*, real*, real* ), void * params, const real guess,
                          const real epsabs, const real epsrel, real *root );

#ifdef SRHD
struct Hydro_HTildeFunction_params_s{
   real MSqr_DSqr;                  // (|Momentum|/Dens)**2
   real Temp;                       // Temperature
   real Constant;                   // LHS of Eq. A3 in "Tseng et al. 2021, MNRAS, 504, 3298"
   EoS_H2TEM_t EoS_HTilde2Temp;     // EoS routine to compute the temperature
   const double *EoS_AuxArray_Flt;  // Auxiliary arrays for EoS_HTilde2Temp
   const int *EoS_AuxArray_Int;     // Auxiliary arrays for EoS_HTilde2Temp
   const real *const *EoS_Table;    // EoS tables for EoS_HTilde2Temp
};
GPU_DEVICE
static void Hydro_HTildeFunction( real HTilde, void *params, real *Func, real *DiffFunc );
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
// Note        :  1. Always apply pressure floor
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when FracPassive is on
//                   --> See the input parameters "FracPassive, NFrac, FracIdx"
//                   --> But note that here we do NOT ensure "sum(mass fraction) == 1.0"
//                       --> It is done by calling Hydro_NormalizePassive() in Hydro_Shared_FullStepUpdate()
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//
// Parameter   :  In                 : Input conserved variables
//                Out                : Output primitive variables
//                MinPres            : Minimum allowed pressure
//                PassiveFloor       : Bitwise flag to specify the passive scalars to be floored
//                FracPassive        : true --> convert passive scalars to mass fraction
//                NFrac              : Number of passive scalars for the option "FracPassive"
//                FracIdx            : Target variable indices for the option "FracPassive"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS_DensEint2Pres  : EoS routine to compute the gas pressure
//                EoS_DensPres2Eint  : EoS routine to compute the gas internal energy
//                EoS_GuessHTilde    : EoS routine to compute guessed reduced enthalpy
//                EoS_HTilde2Temp    : EoS routine to compute temperature
//                EoS_AuxArray_*     : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table          : EoS tables for EoS_DensEint2Pres()
//                EintOut            : Pointer to store the output internal energy
//                                     --> Do nothing if it is NULL
//                                     --> Internal energy floor is not applied
//                LorentzFactorPtr   : Pointer to store the Lorentz factor value
//                                     --> Do nothing if it is NULL
//                                     --> Only used for SRHD
//
// Return      :  Out[], EintOut (optional), LorentzFactorPtr (optional)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Pri( const real In[], real Out[], const real MinPres, const long PassiveFloor,
                    const bool FracPassive, const int NFrac, const int FracIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut, real* LorentzFactorPtr )
{

// conserved --> primitive
#  ifdef SRHD
   real HTilde, Factor, Temp, LorentzFactor;

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical( UNPHY_MODE_CONS, In, NULL_REAL,
                       EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table,
                       PassiveFloor, ERROR_INFO, UNPHY_VERBOSE );
#  endif

   HTilde = Hydro_Con2HTilde( In, EoS_GuessHTilde, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Factor = In[0]*( (real)1.0 + HTilde );

   Out[1] = In[1]/Factor;
   Out[2] = In[2]/Factor;
   Out[3] = In[3]/Factor;

   LorentzFactor = SQRT( (real)1.0 + SQR(Out[1]) + SQR(Out[2]) + SQR(Out[3]) );

   if ( LorentzFactorPtr != NULL )  *LorentzFactorPtr = LorentzFactor;

   const real _LorentzFactor = real(1.0) / LorentzFactor;
   Out[0] = In[0]*_LorentzFactor;

   EoS_HTilde2Temp( HTilde, &Temp, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

   Out[4] = Out[0]*Temp;
   Out[4] = Hydro_CheckMinPres( Out[4], MinPres );


#  else // #ifdef SRHD


   const bool CheckMinPres_Yes = true;
   const real _Rho             = (real)1.0/In[0];

#  ifdef MHD
   const real Bx               = In[ MAG_OFFSET + 0 ];
   const real By               = In[ MAG_OFFSET + 1 ];
   const real Bz               = In[ MAG_OFFSET + 2 ];
   const real Emag             = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  else
   const real Emag             = NULL_REAL;
#  endif

   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = Hydro_Con2Pres( In[0], In[1], In[2], In[3], In[4], In+NCOMP_FLUID, CheckMinPres_Yes, MinPres, PassiveFloor, Emag,
                            EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                            EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, EintOut );

// pressure floor required to resolve the Jeans length
// --> note that currently we do not modify the dual-energy variable (e.g., entropy) accordingly
   if ( JeansMinPres )
   {
      const real Pres0 = Out[4];
      Out[4] = Hydro_CheckMinPres( Pres0, JeansMinPres_Coeff*SQR(Out[0]) );

//    recompute internal energy to be consistent with the updated pressure
      if ( EintOut != NULL  &&  Out[4] != Pres0 )
         *EintOut = EoS_DensPres2Eint( Out[0], Out[4], In+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   }
#  endif // #ifdef SRHD ... else ...


// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

#  ifdef SRHD
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] *= _LorentzFactor;

   const real _Rho = (real)1.0/Out[0];
#  endif

// convert the mass density of target passive scalars to mass fraction
   if ( FracPassive )
      for (int v=0; v<NFrac; v++)   Out[ NCOMP_FLUID + FracIdx[v] ] *= _Rho;
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
//                   when FracPassive is on
//                   --> See the input parameters "FracPassive, NFrac, FracIdx"
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//                5. Convert pressure to internal energy using the input EoS routine by default
//                   --> But one can also specify internal energy directly through *EintIn*,
//                       by which no EoS conversion is required and the input pressure will be useless
//                       --> Mainly used by the option LR_EINT in data reconstruction
//
// Parameter   :  In                : Array storing the input primitive variables
//                Out               : Array to store the output conserved variables
//                FracPassive       : true --> input passive scalars are mass fraction instead of density
//                NFrac             : Number of passive scalars for the option "FracPassive"
//                FracIdx           : Target variable indices for the option "FracPassive"
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_Temp2HTilde   : EoS routine to compute reduced enthalpy
//                EoS_HTilde2Temp   : EoS routine to compute temperature
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensPres2Eint()
//                EoS_Table         : EoS tables for EoS_DensPres2Eint()
//                EintIn            : Pointer storing the input internal energy (see the note above)
//                                    --> Do nothing if it is NULL
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Pri2Con( const real In[], real Out[], const bool FracPassive,
                    const int NFrac, const int FracIdx[], const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_TEM2H_t EoS_Temp2HTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], const real* const EintIn )
{

#  ifdef SRHD
// primitive --> conserved
   real LorentzFactor, Temperature, HTilde, MSqr_DSqr, HTildeFunction, Factor;

   LorentzFactor = SQRT( (real)1.0 + SQR(In[1]) + SQR(In[2]) + SQR(In[3]) );
   Temperature   = In[4]/In[0];
   HTilde        = EoS_Temp2HTilde( Temperature, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Out[0]        = In[0]*LorentzFactor;
   Factor        = Out[0]*HTilde + Out[0];
   Out[1]        = In[1]*Factor;
   Out[2]        = In[2]*Factor;
   Out[3]        = In[3]*Factor;
   MSqr_DSqr     = SQR(Out[1]) + SQR(Out[2]) + SQR(Out[3]);
   MSqr_DSqr    /= SQR(Out[0]);

   struct Hydro_HTildeFunction_params_s params =
   { MSqr_DSqr, Temperature, (real)0.0, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table };

   Hydro_HTildeFunction( HTilde, &params, &HTildeFunction, NULL );

   Out[4]  = MSqr_DSqr + HTildeFunction;
   Out[4] /= (real)1.0 + SQRT( (real)1.0 + MSqr_DSqr + HTildeFunction );
   Out[4] *= Out[0];


// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v]*LorentzFactor;

// convert the mass fraction of target passive scalars back to mass density
   if ( FracPassive )
      for (int v=0; v<NFrac; v++)   Out[ NCOMP_FLUID + FracIdx[v] ] *= In[0];
#  endif


#  else // #ifdef SRHD


// passive scalars
// --> do it before invoking EoS_DensPres2Eint() since the latter requires the mass density
//     instead of mass fraction of passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass fraction of target passive scalars back to mass density
   if ( FracPassive )
      for (int v=0; v<NFrac; v++)   Out[ NCOMP_FLUID + FracIdx[v] ] *= In[0];
#  endif


// primitive --> conserved
   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];

   real Eint, Emag=NULL_REAL;
#  ifdef MHD
   const real Bx = In[ MAG_OFFSET + 0 ];
   const real By = In[ MAG_OFFSET + 1 ];
   const real Bz = In[ MAG_OFFSET + 2 ];
   Emag   = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  endif
   Eint   = ( EintIn == NULL ) ? EoS_DensPres2Eint( In[0], In[4], Out+NCOMP_FLUID, EoS_AuxArray_Flt,
                                                    EoS_AuxArray_Int, EoS_Table )
                               : *EintIn;
   Out[4] = Hydro_ConEint2Etot( Out[0], Out[1], Out[2], Out[3], Eint, Emag );


// B field
#  ifdef MHD
   for (int v=NCOMP_TOTAL; v<NCOMP_TOTAL_PLUS_MAG; v++)  Out[v] = In[v];
#  endif

#  endif // #ifdef SRHD ... else ...

} // FUNCTION : Hydro_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Flux
// Description :  Evaluate hydrodynamic/MHD fluxes from the input conserved variables
//
// Note        :  1. Flux[] and In[] may point to the same array
//                2. Flux[] and In[] should have the size of NCOMP_TOTAL_PLUS_MAG
//                3. By default, it computes pressure using the input EoS routine
//                   --> But one can also specify pressure directly through *PresIn*,
//                       by which no EoS conversion is required
//
// Parameter   :  XYZ               : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux              : Array to store the output fluxes
//                In                : Array storing the input conserved variables
//                MinPres           : Minimum allowed pressure
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table         : EoS tables for EoS_DensEint2Pres()
//                AuxArray          : Array storing
//                                    (1) the input pressure in non-SRHD (see the note above)
//                                        --> Do nothing if it is NULL
//                                    (2) primitive variables in SRHD
//
// Return      :  Flux[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres, const long PassiveFloor,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real* const AuxArray )
{

#  if ( defined GAMER_DEBUG  &&  defined SRHD )
   if ( AuxArray == NULL )
      printf( "ERROR : AuxArray == NULL at file <%s>, line <%d>, function <%s> !!\n",
              ERROR_INFO );
#  endif

   real InRot [ NCOMP_FLUID + NCOMP_MAG ];   // no need to include passive scalars since they don't have to be rotated
#  ifdef SRHD
   real PriRot[ NCOMP_FLUID + NCOMP_MAG ];
#  else
   const bool CheckMinPres_Yes = true;
#  endif

   for (int v=0; v<NCOMP_FLUID; v++)
   {
      InRot [v] = In      [v];
#     ifdef SRHD
      PriRot[v] = AuxArray[v];
#     endif
   }

#  ifdef MHD
   for (int v=NCOMP_FLUID; v<NCOMP_FLUID+NCOMP_MAG; v++)    InRot[v] = In[ v - NCOMP_FLUID + MAG_OFFSET ];
#  endif

   Hydro_Rotate3D( InRot,  XYZ, true, NCOMP_FLUID );
#  ifdef SRHD
   Hydro_Rotate3D( PriRot, XYZ, true, NCOMP_FLUID );
#  endif

#  ifdef MHD
   const real Bx   = InRot[ NCOMP_FLUID + 0 ];
   const real By   = InRot[ NCOMP_FLUID + 1 ];
   const real Bz   = InRot[ NCOMP_FLUID + 2 ];
   const real Emag = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  elif ( !defined SRHD )
   const real Emag = NULL_REAL;
#  endif

#  ifdef SRHD
   const real Pres          = PriRot[4];
   const real LorentzFactor = SQRT( (real)1.0 + SQR(PriRot[1]) + SQR(PriRot[2]) + SQR(PriRot[3]) );
   const real Vx            = PriRot[1] / LorentzFactor;

   Flux[0] = Vx*InRot[0];

#  else // #ifdef SRHD

   const real Pres = ( AuxArray == NULL ) ? Hydro_Con2Pres( InRot[0], InRot[1], InRot[2], InRot[3], InRot[4], In+NCOMP_FLUID,
                                                            CheckMinPres_Yes, MinPres, PassiveFloor, Emag, EoS_DensEint2Pres,
                                                            NULL, NULL,
                                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL )
                                          : AuxArray[0];
   const real _Rho = (real)1.0 / InRot[0];
   const real Vx   = _Rho*InRot[1];

   Flux[0] = InRot[1];
#  endif // #ifdef SRHD ... else ...

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

   Flux[              1 ] += Emag - SQR(Bx);
   Flux[              2 ] -= Bx*By;
   Flux[              3 ] -= Bx*Bz;
   Flux[              4 ] += Vx*Emag - Bx*( Bx*Vx + By*Vy + Bz*Vz );
   Flux[ MAG_OFFSET + 0 ]  = (real)0.0;
   Flux[ MAG_OFFSET + 1 ]  = By*Vx - Bx*Vy;
   Flux[ MAG_OFFSET + 2 ]  = Bz*Vx - Bx*Vz;
#  endif

   Hydro_Rotate3D( Flux, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_Con2Flux



#ifdef SRHD
//-------------------------------------------------------------------------------------------------------
// Function    : Hydro_Con2HTilde
// Description : Conserved variables --> reduced enthalpy
//
// Note        : Using Eq. 15 in "Tseng et al. 2021, MNRAS, 504, 3298"
//
// Parameter   : Con             : Conserved variables
//               EoS_GuessHTilde : EoS routine to compute guessed reduced enthalpy
//               EoS_HTilde2Temp : EoS routine to compute temperature
//               EoS_AuxArray_*  : Auxiliary arrays
//               EoS_Table       : EoS tables
//
// Return      : HTilde          : Reduced enthalpy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2HTilde( const real Con[], const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                       const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                       const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   real HTilde, GuessHTilde, MSqr_DSqr, Constant;

   MSqr_DSqr  = SQR(Con[1]) + SQR(Con[2]) + SQR(Con[3]);
   MSqr_DSqr /= SQR(Con[0]);

   GuessHTilde = EoS_GuessHTilde( Con, &Constant, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

   void (*FuncPtr)( real HTilde, void *params, real *Func, real *DiffFunc ) = &Hydro_HTildeFunction;

// set params->Temp=NAN to force re-computing temperature from the given HTilde
// set params->Constant=-A3_LHS (i.e., negative of the LHS of Eq. A3) since we want to find the root of f(HTilde) - A3_LHS
   struct Hydro_HTildeFunction_params_s params =
   { MSqr_DSqr, NAN, -Constant, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table };

   NewtonRaphsonSolver( FuncPtr, &params, GuessHTilde, TINY_NUMBER, MACHINE_EPSILON, &HTilde );

#  ifdef GAMER_DEBUG
   if ( HTilde <= (real)0.0 )
      printf( "ERROR : HTilde = %14.7e <= 0.0 (Constant %14.7e, GuessHTilde %14.7e) in %s !!\n",
              HTilde, Constant, GuessHTilde, __FUNCTION__ );
#  endif

   return HTilde;

} // FUNCTION : Hydro_Con2HTilde



//-------------------------------------------------------------------------------------------------------
// Function    : Hydro_HTildeFunction
// Description : RHS of Eq. 15 in "Tseng et al. 2021, MNRAS, 504, 3298" plus Params->Constant
//
// Note        :
//
// Parameter   : HTilde         : Reduced specific enthalpy
//               Params         : Pointer to the structure `Hydro_HTildeFunction_params_s`
//               Func           : Eq. 15 in "Tseng et al. 2021, MNRAS, 504, 3298"
//               DiffFunc       : Eq. A1 in "Tseng et al. 2021, MNRAS, 504, 3298"
//               EoS_AuxArray_* : Auxiliary arrays for EoS_DensEint2Pres()
//               EoS_Table      : EoS tables for EoS_DensEint2Pres()
//
// Return      : Func, DiffFunc
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_HTildeFunction( real HTilde, void *Params, real *Func, real *DiffFunc )
{

   struct Hydro_HTildeFunction_params_s *parameters = (struct Hydro_HTildeFunction_params_s *) Params;

   const real         MSqr_DSqr        = parameters->MSqr_DSqr;
   const real         Constant         = parameters->Constant;
   const EoS_H2TEM_t  EoS_HTilde2Temp  = parameters->EoS_HTilde2Temp;
   const double      *EoS_AuxArray_Flt = parameters->EoS_AuxArray_Flt;
   const int         *EoS_AuxArray_Int = parameters->EoS_AuxArray_Int;
   const real *const *EoS_Table        = parameters->EoS_Table;
   real Temp                           = parameters->Temp;
   real DiffTemp;

// DiffTemp required by DiffFunc is only computed when Temp == NAN
#  ifdef GAMER_DEBUG
   if ( Temp == Temp  &&  DiffFunc != NULL )
      printf( "ERROR : Temp (%13.7e) != NAN but DiffFunc != NULL at file <%s>, line <%d>, function <%s> !!\n",
              Temp, ERROR_INFO );
#  endif

// recompute temperature only when the input Temp is NAN, which is used in Hydro_Con2HTilde()
   if ( Temp != Temp )
      EoS_HTilde2Temp( HTilde, &Temp, &DiffTemp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

// Eq. 15
   const real H =  HTilde + (real)1.0;
   const real Factor0 = SQR( H ) + MSqr_DSqr;

   if ( Func != NULL )
   *Func = SQR( HTilde ) + (real)2.0*HTilde - (real)2.0*Temp - (real)2.0*Temp*HTilde
           + SQR( Temp*H )/Factor0 + Constant;

// Eq. A1
   if ( DiffFunc != NULL )
   *DiffFunc = (real)2.0*H - (real)2.0*Temp - (real)2.0*H*DiffTemp +
               (real)2.0*Temp*H*( DiffTemp*H*Factor0 + Temp*MSqr_DSqr ) / SQR( Factor0 );

} // FUNCTION : Hydro_HTildeFunction
#endif // #ifdef SRHD



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinPres
// Description :  Check if the input pressure is greater than the minimum allowed threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
//                   --> Currently it simply sets a minimum allowed value for pressure
//                       --> Set MIN_PRES in the runtime parameter file "Input__Parameter"
//                2. We should also support a minimum **temperature** instead of **pressure**
//                   --> NOT supported yet
//                3. If the input pressure is NaN, return NaN in order to trigger auto-correction such as
//                   "OPT__1ST_FLUX_CORR" and "AUTO_REDUCE_DT"
//
// Parameter   :  InPres  : Input pressure to be corrected
//                MinPres : Minimum allowed pressure
//
// Return      :  InPres != NaN --> max( InPres, MinPres )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinPres( const real InPres, const real MinPres )
{

// call FMAX() only if InPres is not NaN
   if ( InPres == InPres )    return FMAX( InPres, MinPres );
   else                       return InPres;

} // FUNCTION : Hydro_CheckMinPres



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinEint
// Description :  Similar to Hydro_CheckMinPres() except that this function checks the internal energy
//                density (Eint) instead of pressure
//
// Note        :  1. See Hydro_CheckMinPres()
//
// Parameter   :  InEint  : Input Eint to be corrected
//                MinEint : Minimum allowed Eint
//
// Return      :  InEint != NaN --> max( InEint, MinEint )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinEint( const real InEint, const real MinEint )
{

// call FMAX() only if InEint is not NaN
   if ( InEint == InEint )    return FMAX( InEint, MinEint );
   else                       return InEint;

} // FUNCTION : Hydro_CheckMinEint



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinTemp
// Description :  Similar to Hydro_CheckMinPres() except that this function checks the gas temperature
//                instead of pressure
//
// Note        :  1. See Hydro_CheckMinPres()
//
// Parameter   :  InTemp  : Input temperature to be corrected
//                MinTemp : Minimum allowed temperature
//
// Return      :  InTemp != NaN --> max( InTemp, MinTemp )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinTemp( const real InTemp, const real MinTemp )
{

// call FMAX() only if InTemp is not NaN
   if ( InTemp == InTemp )    return FMAX( InTemp, MinTemp );
   else                       return InTemp;

} // FUNCTION : Hydro_CheckMinTemp



#ifndef SRHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinEntr
// Description :  Similar to Hydro_CheckMinPres() except that this function checks the gas entropy
//                instead of pressure
//
// Note        :  1. See Hydro_CheckMinPres()
//
// Parameter   :  InEntr  : Input entropy to be corrected
//                MinEntr : Minimum allowed entropy
//
// Return      :  InEntr != NaN --> max( InEntr, MinEntr )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinEntr( const real InEntr, const real MinEntr )
{

// call FMAX() only if InEntr is not NaN
   if ( InEntr == InEntr )    return FMAX( InEntr, MinEntr );
   else                       return InEntr;

} // FUNCTION : Hydro_CheckMinEntr



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinEintInEngy
// Description :  Ensure that the internal energy density in the input total energy density is greater than
//                a given threshold
//
// Note        :  1. Invoke Hydro_CheckMinEint()
//                2. Input conserved instead of primitive variables
//                3. For MHD, one must provide the magnetic energy density Emag (i.e., 0.5*B^2)
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                InEngy       : Energy density
//                MinEint      : Internal energy density floor
//                PassiveFloor : Bitwise flag to specify the passive scalars to be floored
//                Emag         : Magnetic energy density (0.5*B^2) --> For MHD only
//
// Return      :  Total energy density with internal energy density greater than a given threshold
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinEintInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real InEngy,
                               const real MinEint, const long PassiveFloor, const real Emag )
{

   const bool CheckMinEint_No = false;
   real InEint, OutEint, OutEngy;

   InEint  = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, InEngy, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                             NULL, NULL, NULL, NULL, NULL );
   OutEint = Hydro_CheckMinEint( InEint, MinEint );

// do not modify energy (even the round-off errors) if the input data pass the check
   if ( InEint == OutEint )   OutEngy = InEngy;
   else                       OutEngy = InEngy - InEint + OutEint;

   return OutEngy;

} // FUNCTION : Hydro_CheckMinEintInEngy
#endif // #ifndef SRHD



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_IsUnphysical
// Description :  Check unphysical results
//
// Note        :  1. Support various modes:
//                   UNPHY_MODE_CONS         : Check if the input conserved variables, including passive scalars, are unphysical
//                   UNPHY_MODE_PRIM         : Check if the input primitive variables, including passive scalars, are unphysical
//                   UNPHY_MODE_PASSIVE_ONLY : Check if the input passive scalars are unphysical
//                2. For UNPHY_MODE_CONS with SRHD, we also check if Eq. 15 in "Tseng et al. 2021, MNRAS, 504, 3298"
//                   has a positive root
//                3. For UNPHY_MODE_CONS:
//                   - Mass density and total energy density must be positive
//                   - Internal energy can be zero or even slightly negative if it's within machine precision
//                   - Pressure cannot be negative (only for non-trivial EoS)
//                   For UNPHY_MODE_PRIM:
//                   - Mass density must be positive
//                   - Pressure cannot be negative
//
// Parameter   :  Mode              : UNPHY_MODE_CONS, UNPHY_MODE_PRIM, UNPHY_MODE_PASSIVE_ONLY
//                                    --> See "Note" for details
//                Fields            : Field data to be checked
//                Emag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_*             : EoS parameters
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                File              : __FILE__
//                Line              : __LINE__
//                Function          : __FUNCTION__
//                Verbose           : UNPHY_VERBOSE --> Show error messages
//                                    UNPHY_SILENCE --> Show nothing
//
// Return      :  true  --> Input field is unphysical
//                false --> Otherwise
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool Hydro_IsUnphysical( const IsUnphyMode_t Mode, const real Fields[],
                         const real Emag, const EoS_DE2P_t EoS_DensEint2Pres,
                         const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                         const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                         const real *const EoS_Table[EOS_NTABLE_MAX], const long PassiveFloor,
                         const char File[], const int Line, const char Function[], const IsUnphVerb_t Verbose )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fields == NULL )   printf( "ERROR : access a NULL pointer at file <%s>, line <%d>, function <%s> !!\n",
                                   File, Line, Function );
#  endif


   bool UnphyCell = false;

   switch ( Mode )
   {
//    === check conserved variables, including passive scalars ===
      case UNPHY_MODE_CONS:
      {
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
//          check NaN
            if ( Fields[v] != Fields[v] )
                  UnphyCell = true;

//          check momentum densities
            if ( v == MOMX  ||  v == MOMY  ||  v == MOMZ )
            {
               if ( Fields[v] < -HUGE_NUMBER  ||  Fields[v] > HUGE_NUMBER )
                  UnphyCell = true;
            }

//          check mass and energy densities (which cannot be zero)
            else if ( v < NCOMP_FLUID )
            {
               if ( Fields[v] <= (real)0.0  ||  Fields[v] > HUGE_NUMBER )
                  UnphyCell = true;
            }

//          check passive scalars (which can be zero)
            else
            {
               if ( Fields[v] < (real)0.0  &&  PassiveFloor & BIDX(v) )
                  UnphyCell = true;
               if ( Fields[v] < -HUGE_NUMBER  ||  Fields[v] > HUGE_NUMBER )
                  UnphyCell = true;
            }
         } // for (int v=0; v<NCOMP_TOTAL; v++)

//       check discriminant for SRHD
//       --> positive if and only if Eq. 15 in "Tseng et al. 2021, MNRAS, 504, 3298" has a positive root
#        ifdef SRHD
         const real Msqr         = SQR(Fields[MOMX]) + SQR(Fields[MOMY]) + SQR(Fields[MOMZ]);
         const real Dsqr         = SQR(Fields[DENS]);
         const real E_D          = Fields[ENGY] / Fields[DENS];
         const real M_D          = SQRT( Msqr / Dsqr );
         const real Temp         = SQRT( E_D*E_D + (real)2.0*E_D );
         const real Discriminant = ( Temp + M_D )*( Temp - M_D ); // replace a^2-b^2 with (a+b)*(a-b) to alleviate a catastrophic cancellation

         if ( Discriminant <= (real)0.0  ||  Discriminant > HUGE_NUMBER  ||  Discriminant != Discriminant )
            UnphyCell = true;

#        else // #ifdef SRHD

#        ifndef BAROTROPIC_EOS
//       check internal energy (which can be zero or slightly negative if it's within machine precision)
         const real CheckMinEint_No = false;
         const real Eint = Hydro_Con2Eint( Fields[DENS], Fields[MOMX], Fields[MOMY], Fields[MOMZ], Fields[ENGY],
                                           CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                                           NULL, NULL, NULL, NULL, NULL );

         if ( Eint < (real)-3.0*Fields[ENGY]*MACHINE_EPSILON  ||  Eint > HUGE_NUMBER  ||  Eint != Eint )
            UnphyCell = true;

//       check pressure for non-trivial EoS (which cannot be negative)
//       --> for trivial EoS like EOS_GAMMA, checking internal energy is sufficient and pressure can be
//           slightly negative if it's within machine precision
#        if ( EOS != EOS_GAMMA )
         const real Pres = EoS_DensEint2Pres( Fields[DENS], Eint, Fields+NCOMP_FLUID,
                                              EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

         if ( Pres < (real)0.0  ||  Pres > HUGE_NUMBER  ||  Pres != Pres )
            UnphyCell = true;
#        endif // #if ( EOS != EOS_GAMMA )
#        endif // #ifndef BAROTROPIC_EOS

#        endif // #ifdef SRHD ... else ...

#        ifdef MHD
         if ( Emag < (real)0.0  ||  Emag > HUGE_NUMBER  ||  Emag != Emag )
            UnphyCell = true;
#        endif

//       print out the unphysical values
#        if ( !defined __CUDACC__  ||  defined CHECK_UNPHYSICAL_IN_FLUID )
         if ( UnphyCell  &&  Verbose )
         {
            printf( "ERROR : unphysical conserved variables at file <%s>, line <%d>, function <%s> !!\n",
                    File, Line, Function );
            printf( "D=%14.7e Mx=%14.7e My=%14.7e Mz=%14.7e E=%14.7e",
                    Fields[DENS], Fields[MOMX], Fields[MOMY], Fields[MOMZ], Fields[ENGY] );
#           ifdef SRHD
            printf( " E^2+2*E*D-|M|^2=%14.7e", Discriminant );
#           else
#           ifndef BAROTROPIC_EOS
            printf( " Eint=%14.7e", Eint );
#           if ( EOS != EOS_GAMMA )
            printf( " Pres=%14.7e", Pres );
#           endif
#           endif // #ifndef BAROTROPIC_EOS
#           endif // #ifdef SRHD ... else ...
#           ifdef MHD
            printf( " Emag=%14.7e", Emag );
#           endif
            printf( "\n" );
#           if ( NCOMP_PASSIVE > 0 )
            printf( "Passive:" );
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  printf( " [%d]=%13.7e", v-NCOMP_FLUID, Fields[v] );
            printf( "\n" );
#           endif
         }
#        endif
      } // case UNPHY_MODE_CONS
      break;


//    === check primitive variables, including passive scalars ===
      case UNPHY_MODE_PRIM:
      {
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
//          check NaN
            if ( Fields[v] != Fields[v] )
                  UnphyCell = true;

//          check velocities
            if ( v == MOMX  ||  v == MOMY  ||  v == MOMZ )
            {
               if ( Fields[v] < -HUGE_NUMBER  ||  Fields[v] > HUGE_NUMBER )
                  UnphyCell = true;
            }

//          check mass density (which cannot be zero)
            else if ( v == DENS )
            {
               if ( Fields[v] <= (real)0.0  ||  Fields[v] > HUGE_NUMBER )
                  UnphyCell = true;
            }

//          check pressure (which cannot be negative; a pressure floor should be applied prior to calling this routine)
            else if ( v == ENGY )
            {
               if ( Fields[ENGY] < (real)0.0  ||  Fields[ENGY] > HUGE_NUMBER )
                  UnphyCell = true;
            }

//          check passive scalars (which can be zero)
            else
            {
               if ( Fields[v] < (real)0.0  &&  PassiveFloor & BIDX(v) )
                  UnphyCell = true;
               if ( Fields[v] < -HUGE_NUMBER  ||  Fields[v] > HUGE_NUMBER )
                  UnphyCell = true;
            }
         } // for (int v=0; v<NCOMP_TOTAL; v++)

#        ifdef MHD
         if ( Emag < (real)0.0  ||  Emag > HUGE_NUMBER  ||  Emag != Emag )
            UnphyCell = true;
#        endif

//       print out the unphysical values
#        if ( !defined __CUDACC__  ||  defined CHECK_UNPHYSICAL_IN_FLUID )
         if ( UnphyCell  &&  Verbose )
         {
            printf( "ERROR : unphysical primitive variables at file <%s>, line <%d>, function <%s> !!\n",
                    File, Line, Function );
            printf( "D=%14.7e Vx=%14.7e Vy=%14.7e Vz=%14.7e P=%14.7e",
                    Fields[DENS], Fields[MOMX], Fields[MOMY], Fields[MOMZ], Fields[ENGY] );
#           ifdef MHD
            printf( " Emag=%14.7e", Emag );
#           endif
            printf( "\n" );
#           if ( NCOMP_PASSIVE > 0 )
            printf( "Passive:" );
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  printf( " [%d]=%13.7e", v-NCOMP_FLUID, Fields[v] );
            printf( "\n" );
#           endif
         }
#        endif
      } // case UNPHY_MODE_PRIM
      break;


//    === only check passive scalars ===
      case UNPHY_MODE_PASSIVE_ONLY:
      {
         for (int v=0; v<NCOMP_PASSIVE; v++)
         {
//          check NaN
            if ( Fields[v] != Fields[v] )
               UnphyCell = true;

//          check negative (passive scalars can be zero)
            if ( Fields[v] < (real)0.0  &&  PassiveFloor & BIDX(v+NCOMP_FLUID) )
               UnphyCell = true;

//          check infinity
            if ( Fields[v] < -HUGE_NUMBER  ||  Fields[v] > HUGE_NUMBER )
               UnphyCell = true;
         }

//       print out the unphysical values
#        if ( !defined __CUDACC__  ||  defined CHECK_UNPHYSICAL_IN_FLUID )
         if ( UnphyCell  &&  Verbose )
         {
            printf( "ERROR : unphysical passive scalars at file <%s>, line <%d>, function <%s> !!\n",
                    File, Line, Function );
#           if ( NCOMP_PASSIVE > 0 )
            printf( "Passive:" );
            for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " [%d]=%13.7e", v, Fields[v] );
            printf( "\n" );
#           endif
         }
#        endif
      } // case UNPHY_MODE_PASSIVE_ONLY
      break;


      default:
      {
#        if ( !defined __CUDACC__  ||  defined CHECK_UNPHYSICAL_IN_FLUID )
         printf( "ERROR : unsupported mode (%d) at file <%s>, line <%d>, function <%s> !!\n",
                 Mode, File, Line, Function );
#        endif
      } // default
      break;

   } // switch ( Mode )


   return UnphyCell;

} // FUNCTION : Hydro_IsUnphysical



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_IsUnphysical_Single
// Description :  Check if the input field is NAN or lies outside the accepted range
//
// Parameter   :  Field             : Field data to be checked
//                SingleFieldName   : Name of the target field
//                Min/Max           : Accepted range
//                File              : __FILE__
//                Line              : __LINE__
//                Function          : __FUNCTION__
//                Verbose           : UNPHY_VERBOSE --> Show error messages
//                                    UNPHY_SILENCE --> Show nothing
//
// Return      :  true  --> Input field is unphysical
//                false --> Otherwise
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool Hydro_IsUnphysical_Single( const real Field, const char SingleFieldName[], const real Min, const real Max,
                                const char File[], const int Line, const char Function[], const IsUnphVerb_t Verbose )
{

   bool UnphyCell = false;

// check if the input single field is NAN or lies outside the accepted range
   if ( Field < Min  ||  Field > Max  ||  Field != Field )
      UnphyCell = true;

// print out the unphysical value
#  if ( !defined __CUDACC__  ||  defined CHECK_UNPHYSICAL_IN_FLUID )
   if ( UnphyCell  &&  Verbose )
      printf( "ERROR : invalid %s = %14.7e (min %14.7e, max %14.7e) at file <%s>, line <%d>, function <%s> !!\n",
               (SingleFieldName==NULL)?"unknown field":SingleFieldName, Field, Min, Max,
               File, Line, Function );
#  endif

   return UnphyCell;

} // FUNCTION : Hydro_IsUnphysical_Single



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Pres
// Description :  Evaluate the fluid pressure
//
// Note        :  1. Invoke the EoS routine EoS_DensEint2Pres() to support different EoS
//                2. For MHD, Engy is the total energy density including the magnetic energy Emag=0.5*B^2
//                   and thus one must provide Emag to subtract it
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Energy density (including the magnetic energy density for MHD)
//                Passive           : Passive scalars
//                CheckMinPres      : Apply pressure floor by invoking Hydro_CheckMinPres()
//                                    --> In some cases we actually want to check if pressure becomes unphysical,
//                                        for which this option should be disabled
//                                        --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres           : Pressure floor
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                Emag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_GuessHTilde   : EoS routine to compute guessed reduced enthalpy
//                EoS_HTilde2Temp   : EoS routine to compute temperature
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table         : EoS tables for EoS_DensEint2Pres()
//                EintOut           : Pointer to store the output internal energy
//                                    --> Do nothing if it is NULL
//                                    --> Internal energy floor is not applied
//
// Return      :  Gas pressure (Pres), EintOut (optional)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Pres( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinPres, const real MinPres, const long PassiveFloor, const real Emag,
                     const EoS_DE2P_t EoS_DensEint2Pres, const EoS_GUESS_t EoS_GuessHTilde,
                     const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[],
                     const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut )
{

   real Pres;

#  ifdef SRHD
   real Prim[NCOMP_TOTAL], Cons[NCOMP_TOTAL];

   Cons[0] = Dens;
   Cons[1] = MomX;
   Cons[2] = MomY;
   Cons[3] = MomZ;
   Cons[4] = Engy;

   Hydro_Con2Pri( Cons, Prim, (CheckMinPres)?MinPres:-HUGE_NUMBER, PassiveFloor, false, NULL_INT, NULL,
                  NULL_BOOL, NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, NULL );
   Pres = Prim[4];

#  else // #ifdef SRHD

   const bool CheckMinEint_No = false;
   real Eint;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );
   Pres = EoS_DensEint2Pres( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

   if ( CheckMinPres )   Pres = Hydro_CheckMinPres( Pres, MinPres );

   if ( EintOut != NULL )  *EintOut = Eint;
#  endif // #ifdef SRHD ... else ...

// check unphysical results
#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   if (  Hydro_IsUnphysical_Single( Pres, "output pressure", (real)0.0, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE )  )
   {
      printf( "Input: Dens %14.7e MomX %14.7e MomY %14.7e MomZ %14.7e Engy %14.7e",
              Dens, MomX, MomY, MomZ, Engy );
#     ifndef SRHD
      printf( " Eint %14.7e", Eint );
#     endif
#     ifdef MHD
      printf( " Emag %14.7e", Emag );
#     endif
      printf( "\n" );
#     if ( NCOMP_PASSIVE > 0 )
//    note that we may have Passive==NULL even when NCOMP_PASSIVE>0
//    (e.g., when calling Hydro_Con2Pres() in the Roe solver with GAMMA EoS)
      if ( Passive != NULL )
      {
         printf( "       Passive " );
         for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " [%d] %14.7e", v, Passive[v] );
         printf( "\n" );
      }
#     endif
   }
#  endif // #ifdef CHECK_UNPHYSICAL_IN_FLUID

   return Pres;

} // FUNCTION : Hydro_Con2Pres



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Eint
// Description :  Evaluate the gas internal energy density
//
// Note        :  1. For MHD, Engy is the total energy density including the magnetic energy Emag=0.5*B^2
//                   and thus one must provide Emag to subtract it
//                2. Internal energy density is energy per volume instead of per mass
//
// Parameter   :  Dens            : Mass density
//                MomX/Y/Z        : Momentum density
//                Engy            : Energy density (including the magnetic energy density for MHD)
//                CheckMinEint    : Apply internal energy floor by invoking Hydro_CheckMinEint()
//                                  --> In some cases we actually want to check if internal energy
//                                      becomes unphysical, for which this option should be disabled
//                MinEint         : Internal energy floor
//                PassiveFloor    : Bitwise flag to specify the passive scalars to be floored
//                Emag            : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_GuessHTilde : EoS routine to compute guessed reduced enthalpy
//                EoS_HTilde2Temp : EoS routine to compute temperature
//                EoS_AuxArray_*  : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table       : EoS tables for EoS_DensEint2Pres()
//
// Return      :  Gas internal energy density (Eint)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Eint( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const bool CheckMinEint, const real MinEint, const long PassiveFloor, const real Emag,
                     const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                     const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   real Eint;

#  ifdef SRHD
   real Prim[NCOMP_TOTAL], Cons[NCOMP_TOTAL];
   real HTilde;

   Cons[0] = Dens;
   Cons[1] = MomX;
   Cons[2] = MomY;
   Cons[3] = MomZ;
   Cons[4] = Engy;

   Hydro_Con2Pri( Cons, Prim, -HUGE_NUMBER, PassiveFloor, false, NULL_INT, NULL,
                  NULL_BOOL, NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, NULL );

   HTilde  = Hydro_Con2HTilde( Cons, EoS_GuessHTilde, EoS_HTilde2Temp, EoS_AuxArray_Flt,
                               EoS_AuxArray_Int, EoS_Table );
   Eint    = HTilde*Prim[0] - Prim[4];

#  else // #ifdef SRHD

//###NOTE: assuming Etot = Eint + Ekin + Emag
   Eint    = Engy - (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) / Dens;
#  ifdef MHD
   Eint   -= Emag;
#  endif
#  endif // #ifdef SRHD ... else ...

   if ( CheckMinEint )   Eint = Hydro_CheckMinEint( Eint, MinEint );

   return Eint;

} // FUNCTION : Hydro_Con2Eint



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_ConEint2Etot
// Description :  Evaluate total energy from the input conserved variables and internal energy
//
// Note        :  1. For MHD, total energy density includes the magnetic energy Emag=0.5*B^2
//                2. Internal energy density is energy per volume instead of per mass
//                3. Does NOT support SRHD
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                Eint     : Internal energy density
//                Emag     : Magnetic energy density (0.5*B^2) --> For MHD only
//
// Return      :  Total energy density (including the magnetic energy density for MHD)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_ConEint2Etot( const real Dens, const real MomX, const real MomY, const real MomZ, const real Eint,
                         const real Emag )
{

#  if ( defined SRHD  &&  defined GAMER_DEBUG )
#  ifdef __CUDACC__
   printf( "ERROR :  SRHD does not support Hydro_ConEint2Etot at file <%s>, line <%d>, function <%s> !!\n",
           ERROR_INFO );
#  else
   Aux_Error( ERROR_INFO, "SRHD does not support Hydro_ConEint2Etot !!\n" );
#  endif
#  endif

//###NOTE: assuming Etot = Eint + Ekin + Emag
   real Etot;

   Etot  = (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) / Dens;
   Etot += Eint;
#  ifdef MHD
   Etot += Emag;
#  endif

   return Etot;

} // FUNCTION : Hydro_ConEint2Etot



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Temp
// Description :  Evaluate the fluid temperature
//
// Note        :  1. Invoke the EoS routine EoS_DensEint2Temp() to support different EoS
//                2. Temperature is in kelvin
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Energy density
//                Passive           : Passive scalars
//                CheckMinTemp      : Apply temperature floor by calling Hydro_CheckMinTemp()
//                                    --> In some cases we actually want to check if temperature becomes unphysical,
//                                        for which we don't want to enable this option
//                MinTemp           : Temperature floor
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                Emag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_DensEint2Temp : EoS routine to compute the gas temperature
//                EoS_GuessHTilde   : EoS routine to compute guessed reduced enthalpy
//                EoS_HTilde2Temp   : EoS routine to compute temperature
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Temp()
//                EoS_Table         : EoS tables for EoS_DensEint2Temp()
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Temp( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinTemp, const real MinTemp, const long PassiveFloor, const real Emag,
                     const EoS_DE2T_t EoS_DensEint2Temp, const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                     const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef SRHD
   if ( EoS_GuessHTilde == NULL  ||  EoS_HTilde2Temp == NULL )
   {
#     ifdef __CUDACC__
      printf( "ERROR :  EoS_GuessHTilde == NULL || EoS_HTilde2Temp == NULL at file <%s>, line <%d>, function <%s> !!\n",
               ERROR_INFO );
#     else
      Aux_Error( ERROR_INFO, "EoS_GuessHTilde == NULL || EoS_HTilde2Temp == NULL !!\n" );
#     endif
   }

#  else

   if ( EoS_DensEint2Temp == NULL )
   {
#     ifdef __CUDACC__
      printf( "ERROR : EoS_DensEint2Temp == NULL at file <%s>, line <%d>, function <%s> !!\n",
               ERROR_INFO );
#     else
      Aux_Error( ERROR_INFO, "EoS_DensEint2Temp == NULL !!\n" );
#     endif
   }
#  endif // #ifdef SRHD ... else ...
#  endif // #ifdef GAMER_DEBUG


   real Temp;

#  ifdef SRHD
   real Prim[NCOMP_TOTAL], Cons[NCOMP_TOTAL];

   Cons[0] = Dens;
   Cons[1] = MomX;
   Cons[2] = MomY;
   Cons[3] = MomZ;
   Cons[4] = Engy;

   Hydro_Con2Pri( Cons, Prim, -HUGE_NUMBER, PassiveFloor, false, NULL_INT, NULL,
                  NULL_BOOL, NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, NULL );
   Temp  = Prim[4]/Prim[0];
   Temp *= EoS_AuxArray_Flt[0];

#  else // #ifdef SRHD

   const bool CheckMinEint_No = false;
   real Eint;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );
   Temp = EoS_DensEint2Temp( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  endif // #ifdef SRHD ... else ...

   if ( CheckMinTemp )   Temp = Hydro_CheckMinTemp( Temp, MinTemp );

   return Temp;

} // FUNCTION : Hydro_Con2Temp



#ifndef SRHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Entr
// Description :  Evaluate the fluid entropy
//
// Note        :  1. Invoke the EoS routine EoS_DensEint2Entr() to support different EoS
//                2. We regard the entropy used in the EoS routines and that used in the dual-energy formalism
//                   as two completely separate fields
//                   --> The former is referred to as Entr/ENTR and manipulated by the EoS API, while the latter is
//                       usually referred to as Enpy/Dual/DUAL and manipulated by the routines in CPU_Shared_DualEnergy.cpp
//                   --> This routine, Hydro_Con2Entr(), belongs to the former
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Energy density
//                Passive           : Passive scalars
//                CheckMinEntr      : Apply entropy floor by calling Hydro_CheckMinEntr()
//                                    --> In some cases we actually want to check if entropy becomes unphysical,
//                                        for which we don't want to enable this option
//                MinEntr           : Entropy floor
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                Emag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_DensEint2Entr : EoS routine to compute the gas entropy
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Entr()
//                EoS_Table         : EoS tables for EoS_DensEint2Entr()
//
// Return      :  Gas entropy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Entr( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinEntr, const real MinEntr, const long PassiveFloor, const real Emag,
                     const EoS_DE2S_t EoS_DensEint2Entr, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef SRHD
#  ifdef __CUDACC__
   printf( "ERROR : SRHD does not support entropy evaluation at file <%s>, line <%d>, function <%s> !!\n",
           ERROR_INFO );
#  else
   Aux_Error( ERROR_INFO, "SRHD does not support entropy evaluation !!\n" );
#  endif
#  endif // #ifdef SRHD

   if ( EoS_DensEint2Entr == NULL )
   {
#     ifdef __CUDACC__
      printf( "ERROR : EoS_DensEint2Entr == NULL at file <%s>, line <%d>, function <%s> !!\n",
              __FILE__, __LINE__, __FUNCTION__ );
#     else
      Aux_Error( ERROR_INFO, "EoS_DensEint2Entr == NULL !!\n" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG


   const bool CheckMinEint_No = false;
   real Eint, Entr;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );
   Entr = EoS_DensEint2Entr( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

   if ( CheckMinEntr )   Entr = Hydro_CheckMinEntr( Entr, MinEntr );

   return Entr;

} // FUNCTION : Hydro_Con2Entr
#endif // #ifndef SRHD



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



//-------------------------------------------------------------------------------------------------------
// Function    : NewtonRaphsonSolver
// Description : The one-dimensional root-finder using the Newton's method
//
// Note        : 1. Solving arbitrary one-dimensional function with N parameters (a1,..,aN)
//                  --> i.e. f(a1,..,aN; x) + constant = 0, where constant = Params->Constant
//               2. Iteration stops when either |x1-x0| < EpsAbs + EpsRel*x0 or number of iterations > threshold
//                  --> x1/x0          : the estimated solution in current/previous iteration
//                  --> EpsAbs, EpsRel : See below
//
// Parameter   : FuncPtr : Pointer to the target function with the signature as follows:
//                         --> (*FuncPtr)( real Unknown, void *Params, real *Func, real *DiffFunc )
//                             --> Unknown  : Independent variables of the target function
//                             --> Params   : Pointer to a user-defined structure that groups the N parameters and the `constant`
//                             --> Func     : Evaluation of FuncPtr at `Unknown`
//                             --> DiffFunc : Evaluation of the derivative of FuncPtr w.r.t `x` at `Unknown`
//               Params  : See above
//               Guess   : Initial guess of x
//               EpsAbs  : Absolute error between the current and previous solution
//               EpsRel  : Relative error between the current and previous solution
//               Root    : Pointer to the root of the target function
//
// Return      : Root
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void NewtonRaphsonSolver( void (*FuncPtr)( real Unknown, void *Params, real *Func, real *DiffFunc ),
                          void *Params, const real Guess, const real EpsAbs, const real EpsRel, real *Root )
{

   int Iter = 0;

#  ifdef FLOAT8
   int MaxIter = 20;
#  else
   int MaxIter = 10;
#  endif

   real Func, DiffFunc, Delta, Tolerance;

   *Root = Guess;

   do
   {
      Iter++;

      FuncPtr( *Root, Params, &Func, &DiffFunc );

#     ifdef GAMER_DEBUG
      if ( DiffFunc == (real)0.0 )
         printf( "ERROR : derivative is zero at file <%s>, line <%d>, function <%s> !!\n", ERROR_INFO );

      if ( Func != Func  ||  Func < -HUGE_NUMBER  ||  Func > HUGE_NUMBER )
         printf( "ERROR : function value is not finite at file <%s>, line <%d>, function <%s> !!\n", ERROR_INFO );

      if ( DiffFunc != DiffFunc  ||  DiffFunc < -HUGE_NUMBER  ||  DiffFunc > HUGE_NUMBER )
         printf( "ERROR : derivative value is not finite at file <%s>, line <%d>, function <%s> !!\n", ERROR_INFO );
#     endif

      Delta     = Func/DiffFunc;
      Tolerance =  EpsRel*FABS(*Root) + EpsAbs;
      *Root     = *Root - Delta;

   } while ( FABS(Delta) >= Tolerance  &&  Iter < MaxIter );

} // FUNCTION : NewtonRaphsonSolver



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
