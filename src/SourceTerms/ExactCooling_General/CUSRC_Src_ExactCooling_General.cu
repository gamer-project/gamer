#include "CUFLU.h"
#include "Global.h"

#ifdef EXACT_COOLING_GENERAL

// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUDA_CheckError.h"
#include "CUDA_ConstMemory.h"
#include "CUFLU_Shared_FluUtility.cu"

extern double *d_SrcECG_TEF_Temp;
extern double *d_SrcECG_TEF_lambda;
extern double *d_SrcECG_TEF_alpha;
extern double *d_SrcECG_TEF_Y;

#endif // #ifdef __CUDACC__


// local function prototypes
#ifndef __CUDACC__

void Src_SetAuxArray_ExactCooling_General( double [], int [] );
void Src_SetConstMemory_ExactCooling_General( const double AuxArray_Flt[], const int AuxArray_Int[],
                                              double *&DevPtr_Flt, int *&DevPtr_Int );
void Src_PassData2GPU_ExactCooling_General();
void Src_SetCPUFunc_ExactCooling_General( SrcFunc_t & );
#ifdef GPU
void Src_SetGPUFunc_ExactCooling_General( SrcFunc_t & );
void CUAPI_MemFree_ExactCooling_General();
#endif
void Src_WorkBeforeMajorFunc_ExactCooling_General( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                                   double AuxArray_Flt[], int AuxArray_Int[] );
void Src_End_ExactCooling_General();
void Cooling_efficiency_table( const double Z, double TEF_Temp[], double TEF_lambda[] );
#endif // #ifdef __CUDACC__

GPU_DEVICE static
double TEF_PiecewisePowerLaw( const double Temp, const int k, const double TEF_Temp[], const double TEF_lambda[], const double TEF_alpha[], const double TEF_Y[],
                              const double AuxArray_Flt[], const int AuxArray_Int[] );
GPU_DEVICE static
double TEF_inverse_PiecewisePowerLaw( const double Y, const double TEF_Temp[], const double TEF_lambda[], const double TEF_alpha[], const double TEF_Y[],
                                      const double AuxArray_Flt[], const int AuxArray_Int[] );

/********************************************************
1. General exact-cooling source term
   --> Enabled by the compilation option "EXACT_COOLING_GENERAL" and the runtime option "SRC_EXACTCOOLING_GENERAL"

2. This file is shared by both CPU and GPU
   
   CUSRC_Src_ExactCooling_General.cu -> CPU_Src_ExactCooling_General.cpp

3. Four steps are required to implement a source term

   I.   Set auxiliary arrays
   II.  Implement the source-term function
   III. [Optional] Add the work to be done every time
        before calling the major source-term function
   IV.  Set initialization functions

4. The source-term function must be thread-safe and
   not use any global variable
********************************************************/



// =======================
// I. Set auxiliary arrays
// =======================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetAuxArray_ExactCooling_General
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_ExactCooling_General()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_ECG defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_ExactCooling_General( double AuxArray_Flt[], int AuxArray_Int[] )
{
   // ==========================
   // Gas settings
   // ==========================
   const double Z     = SrcTerms.ECG_Z;                                                 // metallicity (Z = 1.0 means solar metallicity)
   const double X_H   = 0.7;                                                            // hydrogen mass fraction
   const double X_Z   = 0.018 * Z;                                                      // metal mass fraction (X_Z = 0.018 means solar-abundace)

   // mean molecular weights (for full ionization regime)
   const double mu_e  = 2.0 / ( 1.0 + X_H );                                            // mean electron molecular weight
   const double mu_H  = 1.0 / X_H;                                                      // mean proton molecular weight
   const double mu    = 1.0 / ( 2.0 * X_H + 0.75 * ( 1.0 - X_H - X_Z ) + 0.5 * X_Z );   // mean total molecularweight
                        
   // ============================================
   // Piecewise power law cooling functions
   //=============================================
   const int    TEF_N    = SrcTerms.ECG_TEF_N;                                          // number of temperature points
   const double TEF_TN   = 1.0e8;                                                       // maximum temperature (K)
   const double TEF_Tmin = 1.0e4;                                                       // minimum temperature (K)

   // store them in the aux arrays
   AuxArray_Flt[0]  = Z;
   AuxArray_Flt[1]  = X_H;                                                                                                  
   AuxArray_Flt[2]  = X_Z;       

   AuxArray_Flt[3]  = mu_e;                                                                    
   AuxArray_Flt[4]  = mu_H;                                                                                     
   AuxArray_Flt[5]  = mu;   

   AuxArray_Flt[6]  = Const_mp / UNIT_M;           // proton mass (g)
   AuxArray_Flt[7]  = Const_kB / UNIT_E;           // Boltzmann constant (erg/K)
   AuxArray_Flt[8]  = GAMMA;                       // ratio of specific heats    
   
   AuxArray_Flt[9]  = TEF_TN;
   AuxArray_Flt[10] = TEF_Tmin;

   AuxArray_Int[0]  = TEF_N;

} // FUNCTION : Src_SetAuxArray_ExactCooling_General
#endif // #ifndef __CUDACC__



// ======================================
// II. Implement the source-term function
// ======================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_ExactCooling_General
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_ExactCooling_General() for the values stored in AuxArray_Flt/Int[]
//                3. Follow Townsend (2009, ApJS, 181, 391) to implement the exact integration of the piecewise power law cooling function
//                4. Shared by both CPU and GPU
//
// Parameter   :  fluid             : Fluid array storing both the input and updated values
//                                    --> Including both active and passive variables
//                B                 : Cell-centered magnetic field
//                SrcTerms          : Structure storing all source-term variables
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                x/y/z             : Target physical coordinates
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                EoS               : EoS object
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_ExactCooling_General( real fluid[], const real B[],
                                      const SrcTerms_t *SrcTerms, const real dt, const real dh,
                                      const double x, const double y, const double z,
                                      const double TimeNew, const double TimeOld,
                                      const real MinDens, const real MinPres, const real MinEint, const long PassiveFloor,
                                      const EoS_t *EoS, const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif
   
   // ==========================
   // (1) read parameters
   // ==========================                         
   // mean molecular weights (for full ionization regime)
   const double mu_e  = AuxArray_Flt[3];                             // mean electron molecular weight
   const double mu_H  = AuxArray_Flt[4];                             // mean proton molecular weight
   const double mu    = AuxArray_Flt[5];                             // mean total molecular weight

   // Constants
   const double mp    = AuxArray_Flt[6];                             // proton mass (g)
   const double kB    = AuxArray_Flt[7];                             // Boltzmann constant (erg/K)
   const double gamma = AuxArray_Flt[8];                             // ratio of specific heats

   // ===========================================================
   
   // Piecewise power law cooling function
   const int    TEF_N    = AuxArray_Int[0];                          // number of temperature points
   const double TEF_TN   = AuxArray_Flt[9];                          // maximum temperature (K)
   const double TEF_Tmin = AuxArray_Flt[10];                         // minimum temperature (K)

#  ifdef __CUDACC__
   const double *TEF_Temp   = SrcTerms->ECG_TEF_Temp_DevPtr;         // temperature points (K)
   const double *TEF_lambda = SrcTerms->ECG_TEF_lambda_DevPtr;       // cooling efficiency (erg cm^3 / s)
   const double *TEF_alpha  = SrcTerms->ECG_TEF_alpha_DevPtr;        // power law slope
   const double *TEF_Y      = SrcTerms->ECG_TEF_Y_DevPtr;            // exact cooling integration constant
#  else
   const double *TEF_Temp   = h_SrcECG_TEF_Temp;
   const double *TEF_lambda = h_SrcECG_TEF_lambda;
   const double *TEF_alpha  = h_SrcECG_TEF_alpha;
   const double *TEF_Y      = h_SrcECG_TEF_Y;
#  endif // #ifdef __CUDACC__


   // ==========================
   // (2) calculate parameters
   // ==========================
   // mass density
   const double rho  = fluid[DENS];           // g cm^-3                  

   // number densities
   const double n_H  = rho / (mu_H * mp);     // cm^-3       
   const double n_e  = rho / (mu_e * mp);     // cm^-3                   
   const double n    = rho / (mu   * mp);     // cm^-3                  

// please change the length of Yks to your N value
   double Tnow, Tnew, Tk, Eint, Eintf, Emag, Pres, lambdaTnow, Y_T, Y, tcool;
   int    k;

   // ===========================================================
   // (3) calculate the current temperature and internal energy
   // ===========================================================
   // current internal energy
   const bool CheckMinEint_No = false;
   Emag  = 0.0;
#  ifdef __CUDACC__
   Eint  = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                           CheckMinEint_No, NULL_REAL, PassiveFloor, Emag, 
                           EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr, 
                           EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                           EoS->Table );
#  else
   Eint  = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                           CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int,
                           h_EoS_Table );
#  endif

   // current temperature
#  ifdef __CUDACC__
   Tnow = EoS->DensEint2Temp_FuncPtr( fluid[DENS], Eint, NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table  );
#  else
   Tnow = EoS_DensEint2Temp_CPUPtr  ( fluid[DENS], Eint, NULL, EoS_AuxArray_Flt,        EoS_AuxArray_Int,        h_EoS_Table );
#  endif

   // =====================================================================
   // (4) calculate the new temperature and internal energy after cooling
   // =====================================================================
   // new temperature
   if ( Tnow >= TEF_Tmin && Tnow <= TEF_TN )
   {
      // k : the index of the temperature interval where Tnow falls into
      for (int i = TEF_N-2; i >= 0; i--)
      {
         if ( Tnow >= TEF_Temp[i] )
         {
            k = i+1;
            break;
         }
      }
      
      // power law cooling function: lambdak * ( Tnow / Tk ) ^ ( alphak )
      Tk            = TEF_Temp[k-1];                                                // K
      lambdaTnow    = TEF_lambda[k-1] * pow( Tnow / Tk, TEF_alpha[k-1] );           // erg cm^3 s^-1

      // compute the current cooling time and store it to the field TCOOLG
      tcool         = n * kB * Tnow / ( (gamma - 1.0) * n_e * n_H * lambdaTnow );
      fluid[TCOOLG] = tcool;

      // do NOT update the internal energy if dt = 0
      if ( dt == 0.0 )   return;

      // Y(T)
      Y_T      = TEF_PiecewisePowerLaw( Tnow, k, TEF_Temp, TEF_lambda, TEF_alpha, TEF_Y, AuxArray_Flt, AuxArray_Int );

      Y        = Y_T + Tnow / TEF_TN * TEF_lambda[TEF_N-1] / lambdaTnow * dt / tcool;
      
      Tnew     = TEF_inverse_PiecewisePowerLaw( Y, TEF_Temp, TEF_lambda, TEF_alpha, TEF_Y, AuxArray_Flt, AuxArray_Int );
      
      if (Tnew < TEF_Tmin)
         Tnew = TEF_Tmin;

   }
   else
   {
#     ifdef GAMER_DEBUG
      printf( "ERROR : Tnow = %24.16e is out of range (min: %24.16e, max: %24.16e) at TimeNew = %24.16e !!\n",
              Tnow, TEF_Tmin, TEF_TN, TimeNew );
#     endif
      
      fluid[ENGY]   = NAN;
      fluid[TCOOLG] = NAN;
      return;
   }

   // new internal energy
#  ifdef __CUDACC__
   Pres  = EoS->DensTemp2Pres_FuncPtr( fluid[DENS], Tnew,  NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table  );
   Eintf = EoS->DensPres2Eint_FuncPtr( fluid[DENS], Pres , NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table  );
#  else
   Pres  = EoS_DensTemp2Pres_CPUPtr  ( fluid[DENS], Tnew,  NULL, EoS_AuxArray_Flt,        EoS_AuxArray_Int,        h_EoS_Table );
   Eintf = EoS_DensPres2Eint_CPUPtr  ( fluid[DENS], Pres , NULL, EoS_AuxArray_Flt,        EoS_AuxArray_Int,        h_EoS_Table );
#  endif

   // convert back
   fluid[ENGY] = Hydro_ConEint2Etot  ( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eintf, Emag );


} // FUNCTION : Src_ExactCooling_General


//-------------------------------------------------------------------------------------------------------
// Function    : TEF_PiecewisePowerLaw
// Description : Temporal evolution function (TEF) for Piecewise Power Law Cooling Function
//
// Note        :
//
// Parameter   : Temp             : Temperature in Kelvin
//               k                : Temperature interval index
//               TEF_Temp         : Temperature points of the piecewise power law
//               TEF_lambda       : The cooling efficiency in piecewise power law
//               TEF_alpha        : The power of the cooling efficiency in piecewise power law
//               TEF_Y            : The constant TEF in the specific temparature region
//               AuxArray_Flt/Int : Auxiliary arrays (see the Note above)
//
// Return      : TEF(Temp)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE static
double TEF_PiecewisePowerLaw( const double Temp, const int k, const double TEF_Temp[], const double TEF_lambda[], const double TEF_alpha[], const double TEF_Y[],
                              const double AuxArray_Flt[], const int AuxArray_Int[] )
{
   const int    TEF_N  = AuxArray_Int[0];
   const double TEF_TN = AuxArray_Flt[9];
   const double TEF_LN = TEF_lambda[TEF_N-1];

   const double Tk = TEF_Temp[k-1];
   const double Lk = TEF_lambda[k-1];
   const double ak = TEF_alpha[k-1];

   double TEF;

   if ( ak == 1.0 )
   {
      TEF = TEF_Y[k-1] + TEF_LN / Lk * Tk / TEF_TN * log( Tk / Temp );
   }
   else
   {
      TEF = TEF_Y[k-1] + 1.0 / (1.0 - ak) * TEF_LN / Lk * Tk / TEF_TN * (1.0 - pow( Tk / Temp, (ak - 1.0)));
   }

   return TEF;

} // FUNCTION : TEF_PiecewisePowerLaw


//-------------------------------------------------------------------------------------------------------
// Function    : TEF_inverse_PiecewisePowerLaw
// Description : Inverse Temporal evolution function (TEF^-1) for Piecewise Power Law Cooling Function
//
// Note        :
//
// Parameter   : Y                : Input parameter for the inverse TEF
//               TEF_Temp         : Temperature points of the piecewise power law
//               TEF_lambda       : The cooling efficiency in piecewise power law
//               TEF_alpha        : The power of the cooling efficiency in piecewise power law
//               TEF_Y            : The constant TEF in the specific temparature region
//               AuxArray_Flt/Int : Auxiliary arrays (see the Note above)
//               
//
// Return      : inverseTEF(Y)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE static
double TEF_inverse_PiecewisePowerLaw( const double Y, const double TEF_Temp[], const double TEF_lambda[], const double TEF_alpha[], const double TEF_Y[],
                                      const double AuxArray_Flt[], const int AuxArray_Int[] )
{
   
   const int    TEF_N    = AuxArray_Int[0];
   const double TEF_TN   = AuxArray_Flt[9];
   const double TEF_Tmin = AuxArray_Flt[10];
   
   // cooling floor
   if ( Y >= TEF_Y[0] )
       return TEF_Tmin;
   
   int k = -1;
   for ( int i = 1; i < TEF_N; i++ )
   {
      if ( TEF_Y[i] <= Y && Y < TEF_Y[i-1] )
      {
         k = i;
         break;
      }
   }
   if ( k < 1 )
   {
   #  ifdef GAMER_DEBUG   
      printf( "ERROR : Invalid Y = %e in %s !!\n", Y, __FUNCTION__ );
   #  endif
      return TEF_Tmin;
   }

   const double Tk = TEF_Temp[k-1];
   const double Lk = TEF_lambda[k-1];
   const double ak = TEF_alpha[k-1];
   const double Yk = TEF_Y[k-1];

   const double TEF_LN = TEF_lambda[TEF_N-1];

   // inverseY(Y)
   double inverseY;
   if ( ak == 1.0 )
   {
      inverseY = Tk * exp(-Lk / TEF_LN * TEF_TN / Tk * (Y - Yk));
   }
   else
   {
      inverseY = Tk * pow( (1.0 - (1.0 - ak) * Lk / TEF_LN * TEF_TN / Tk * (Y - Yk)), 1.0 / (1.0 - ak));
   }
   
   return inverseY;
}  // FUNCTION : TEF_inverse_PiecewisePowerLaw



// ==================================================
// III. [Optional] Add the work to be done every time
//      before calling the major source-term function
// ==================================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc_ExactCooling_General
// Description :  Specify work to be done every time before calling the major source-term function
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc()
//                   --> By linking to "Src_WorkBeforeMajorFunc_User_Ptr" in Src_Init_ExactCooling_General()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  lv               : Target refinement level
//                TimeNew          : Target physical time to reach
//                TimeOld          : Physical time before update
//                                   --> The major source-term function will update the system from TimeOld to TimeNew
//                dt               : Time interval to advance solution
//                                   --> Physical coordinates : TimeNew - TimeOld == dt
//                                       Comoving coordinates : TimeNew - TimeOld == delta(scale factor) != dt
//                AuxArray_Flt/Int : Auxiliary arrays
//                                   --> Can be used and/or modified here
//                                   --> Must call Src_SetConstMemory_ExactCooling_General() after modification
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_WorkBeforeMajorFunc_ExactCooling_General( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                                   double AuxArray_Flt[], int AuxArray_Int[] )
{

// nothing to do here

} // FUNCTION : Src_WorkBeforeMajorFunc_ExactCooling_General
#endif // #ifndef __CUDACC__


#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_PassData2GPU_ExactCooling_General
// Description :  Transfer data to GPU
//
// Note        :  1. Invoked by Src_Init_ExactCooling_General()
//                2. Use synchronous transfer
//
// Parameter   :  None
// Return      :  None
// -------------------------------------------------------------------------------------------------------
void Src_PassData2GPU_ExactCooling_General()
{

   const long ECG_TEF_MemSize = sizeof(double) * 201;
   
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcECG_TEF_Temp,   ECG_TEF_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcECG_TEF_lambda, ECG_TEF_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcECG_TEF_alpha,  ECG_TEF_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_SrcECG_TEF_Y,      ECG_TEF_MemSize )  );

// store the device pointers in SrcTerms when using GPU
   SrcTerms.ECG_TEF_Temp_DevPtr   = d_SrcECG_TEF_Temp;
   SrcTerms.ECG_TEF_lambda_DevPtr = d_SrcECG_TEF_lambda;
   SrcTerms.ECG_TEF_alpha_DevPtr  = d_SrcECG_TEF_alpha;
   SrcTerms.ECG_TEF_Y_DevPtr      = d_SrcECG_TEF_Y;

// use synchronous transfer
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcECG_TEF_Temp,   h_SrcECG_TEF_Temp,   ECG_TEF_MemSize, cudaMemcpyHostToDevice )  ); 
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcECG_TEF_lambda, h_SrcECG_TEF_lambda, ECG_TEF_MemSize, cudaMemcpyHostToDevice )  );
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcECG_TEF_alpha,  h_SrcECG_TEF_alpha,  ECG_TEF_MemSize, cudaMemcpyHostToDevice )  );
   CUDA_CHECK_ERROR(  cudaMemcpy( d_SrcECG_TEF_Y,      h_SrcECG_TEF_Y,      ECG_TEF_MemSize, cudaMemcpyHostToDevice )  );

} // FUNCTION : Src_PassData2GPU_ExactCooling_General
#endif // #ifdef __CUDACC__


// ================================
// IV. Set initialization functions
// ================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_ExactCooling_General;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetCPU/GPUFunc_ExactCooling_General
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_ExactCooling_General()
//                2. Call-by-reference
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetGPUFunc_ExactCooling_General( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
} // FUNCTION : Src_SetGPUFunc_ExactCooling_General

#else

void Src_SetCPUFunc_ExactCooling_General( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
} // FUNCTION : Src_SetCPUFunc_ExactCooling_General

#endif // #ifdef __CUDACC__ ... else ...



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetConstMemory_ExactCooling_General
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Src_Init_ExactCooling_General() and, if necessary, Src_WorkBeforeMajorFunc_ExactCooling_General()
//                3. SRC_NAUX_EXACTCOOLING_GENERAL is defined in Macro.h
//
// Parameter   :  AuxArray_Flt/Int : Auxiliary arrays to be copied to the constant memory
//                DevPtr_Flt/Int   : Pointers to store the addresses of constant memory arrays
//
// Return      :  c_Src_ExactCooling_General_AuxArray_Flt[], c_Src_ExactCooling_General_AuxArray_Int[], DevPtr_Flt, DevPtr_Int
//---------------------------------------------------------------------------------------------------
void Src_SetConstMemory_ExactCooling_General( const double AuxArray_Flt[], const int AuxArray_Int[],
                                              double *&DevPtr_Flt, int *&DevPtr_Int )
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_ECG_AuxArray_Flt, AuxArray_Flt, SRC_NAUX_ECG*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_ECG_AuxArray_Int, AuxArray_Int, SRC_NAUX_ECG*sizeof(int   ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Flt, c_Src_ECG_AuxArray_Flt) );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Int, c_Src_ECG_AuxArray_Int) );

} // FUNCTION : Src_SetConstMemory_ExactCooling_General
#endif // #ifdef __CUDACC__



#ifndef __CUDACC__

//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_ExactCooling_General
// Description :  Initialize a user-specified source term
//
// Note        :  1. Set auxiliary arrays by invoking Src_SetAuxArray_*()
//                   --> Copy to the GPU constant memory and store the associated addresses
//                2. Set the source-term function by invoking Src_SetCPU/GPUFunc_*()
//                3. Set the function pointers "Src_WorkBeforeMajorFunc_ExactCooling_General_Ptr" and "Src_End_ExactCooling_General_Ptr"
//                4. Invoked by Src_Init()
//                   --> Enable it by linking to the function pointer "Src_Init_ExactCooling_General_Ptr"
//                5. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_Init_ExactCooling_General()
{

// set the auxiliary arrays
   Src_SetAuxArray_ExactCooling_General( Src_ECG_AuxArray_Flt, Src_ECG_AuxArray_Int );

// copy the auxiliary arrays to the GPU constant memory and store the associated addresses
#  ifdef GPU
   Src_SetConstMemory_ExactCooling_General( Src_ECG_AuxArray_Flt, Src_ECG_AuxArray_Int,
                                            SrcTerms.ECG_AuxArrayDevPtr_Flt, SrcTerms.ECG_AuxArrayDevPtr_Int );
#  else
   SrcTerms.ECG_AuxArrayDevPtr_Flt = Src_ECG_AuxArray_Flt;
   SrcTerms.ECG_AuxArrayDevPtr_Int = Src_ECG_AuxArray_Int;
#  endif

// set the major source-term function
   Src_SetCPUFunc_ExactCooling_General( SrcTerms.ECG_CPUPtr );

#  ifdef GPU
   Src_SetGPUFunc_ExactCooling_General( SrcTerms.ECG_GPUPtr );
   SrcTerms.ECG_FuncPtr = SrcTerms.ECG_GPUPtr;
#  else
   SrcTerms.ECG_FuncPtr = SrcTerms.ECG_CPUPtr;
#  endif

   if ( OPT__INIT == INIT_BY_RESTART )
      for (int i=0; i<NLEVEL; i++)   SrcTerms.ECG_TCoolGInit[i] = true;
   else
      for (int i=0; i<NLEVEL; i++)   SrcTerms.ECG_TCoolGInit[i] = false;

   if ( SrcTerms.ECG_TEF_N != 201 )
      Aux_Error( ERROR_INFO, "SRC_ECG_TEF_N must be 201 for this table !!\n" );

// allocate h_SrcECG_* arrays
   h_SrcECG_TEF_Temp   = new double [SrcTerms.ECG_TEF_N];
   h_SrcECG_TEF_lambda = new double [SrcTerms.ECG_TEF_N];
   h_SrcECG_TEF_alpha  = new double [SrcTerms.ECG_TEF_N];
   h_SrcECG_TEF_Y      = new double [SrcTerms.ECG_TEF_N];
   
   SrcTerms.ECG_TEF_Temp_DevPtr   = h_SrcECG_TEF_Temp;
   SrcTerms.ECG_TEF_lambda_DevPtr = h_SrcECG_TEF_lambda;
   SrcTerms.ECG_TEF_alpha_DevPtr  = h_SrcECG_TEF_alpha;
   SrcTerms.ECG_TEF_Y_DevPtr      = h_SrcECG_TEF_Y;

// initialize the cooling function
   const int    TEF_N    = Src_ECG_AuxArray_Int[0];   // number of temperature points
   
   const double Z        = Src_ECG_AuxArray_Flt[0];   // metallicity
   const double X_H      = Src_ECG_AuxArray_Flt[1];   // hydrogen mass fraction
   const double X_Z      = Src_ECG_AuxArray_Flt[2];   // metal mass fraction
   const double mu_e     = Src_ECG_AuxArray_Flt[3];   // mean electron molecular weight
   const double mu_H     = Src_ECG_AuxArray_Flt[4];   // mean proton molecular weight
   const double mu       = Src_ECG_AuxArray_Flt[5];   // mean total molecular weight

   const double mp       = Src_ECG_AuxArray_Flt[6];   // proton mass (g)
   const double kB       = Src_ECG_AuxArray_Flt[7];   // Boltzmann constant (erg/K)
   const double gamma    = Src_ECG_AuxArray_Flt[8];   // ratio of specific heats

   const double TEF_TN   = Src_ECG_AuxArray_Flt[9];   // maximum temperature (K)
   const double TEF_Tmin = Src_ECG_AuxArray_Flt[10];  // minimum temperature (K)

// build cooling efficiency table at metallicity Z
   Cooling_efficiency_table( Z, h_SrcECG_TEF_Temp, h_SrcECG_TEF_lambda );

// initialize the power of the cooling efficiency
// (1) k = 201
   h_SrcECG_TEF_alpha[TEF_N-1] = 0.0;   // is never required 
   
   double Tk, Tkp1, lambdak, lambdakp1;
// (2) k = 1, 2, ..., 200
   for (int k=TEF_N-1; k>0; k--)
   {
      Tk        = h_SrcECG_TEF_Temp[k-1];
      Tkp1      = h_SrcECG_TEF_Temp[k];
      lambdak   = h_SrcECG_TEF_lambda[k-1];
      lambdakp1 = h_SrcECG_TEF_lambda[k];

      h_SrcECG_TEF_alpha[k-1] = ( log10(lambdakp1) - log10(lambdak) ) / ( log10(Tkp1) - log10(Tk) );
   }

// initialize the constant TEF of integration
// (1) k = 201
   h_SrcECG_TEF_Y[TEF_N-1] = 0.0;   // YN = 0.0

   double Tk_2, Tkp1_2, alphak, lambdak_2, lambdaN;
// (2) k = 1, 2, ..., 200
   for (int k=TEF_N-1; k>0; k--)
   {
      Tk_2      = h_SrcECG_TEF_Temp[k-1];
      Tkp1_2    = h_SrcECG_TEF_Temp[k];
      alphak    = h_SrcECG_TEF_alpha[k-1];
      lambdak_2 = h_SrcECG_TEF_lambda[k-1];
      lambdaN   = h_SrcECG_TEF_lambda[TEF_N-1];

      if ( alphak == 1.0 )
         h_SrcECG_TEF_Y[k-1] = h_SrcECG_TEF_Y[k] - lambdaN / lambdak_2 * Tk_2 / TEF_TN * log( Tk_2 / Tkp1_2 );
      else
         h_SrcECG_TEF_Y[k-1] = h_SrcECG_TEF_Y[k] - 1.0 / ( 1.0 - alphak ) * lambdaN / lambdak_2 * Tk_2 / TEF_TN * ( 1.0 - pow( Tk_2 / Tkp1_2, alphak - 1.0 ) ); 
   }   

#  ifdef GPU
   Src_PassData2GPU_ExactCooling_General();
#  endif

} // FUNCTION : Src_Init_ExactCooling_General



// -------------------------------------------------------------------------------------------------------
// Function    : Cooling_efficiency_table
// Description : Build the cooling efficiency table at metallicity Z
//
// Note        :  1. Invoked by Src_Init_ExactCooling_General()
//                2. The cooling efficiency table is interpolated from the Gnat & Sternberg (2007) table for CIE
//
// Parameter   :  Z : Metallicity, 1.0 mean solar metallicity
//
// Return      :  TEF_Temp  : Temperature points
//                TEF_lambda: Cooling efficiency
//-------------------------------------------------------------------------------------------------------

void Cooling_efficiency_table( const double Z, double TEF_Temp[], double TEF_lambda[] )
{

// read the cooling efficiency table   
   double      *CoolingEfficiencyTable     = NULL;
   const char  *CoolingEfficiencyTableFile = "cooling_Gnat_Sternberg_table.txt";
   const bool   RowMajor_No                = false;
   const bool   AllocMem_Yes               = true;
   const int    Ncol                       = 6;
   const int    Col[Ncol]                  = {0, 1, 2, 3, 4, 5};
   const double ZTable[5]                  = {1.0e-3, 1.0e-2, 1.0e-1, 1.0, 2.0};
   const int    Nrow                       = Aux_LoadTable( CoolingEfficiencyTable, CoolingEfficiencyTableFile, Ncol, Col, RowMajor_No, AllocMem_Yes );

   if ( Nrow != 201 )
      Aux_Error( ERROR_INFO, "unexpected cooling table row count (%d; expected 201) !!\n", Nrow );

   if ( Z < ZTable[0] || Z > ZTable[4] )
      Aux_Error( ERROR_INFO, "metallicity Z = %13.7e is outside the table range !!\n", Z );

// cooling efficiency at metallicity Z
   for (int i=0; i<Nrow; i++)
   {
      TEF_Temp[i] = CoolingEfficiencyTable[i];   // column 0: temperature

      double LambdaByZ[5];
      for (int m=0; m<5; m++)
         LambdaByZ[m] = CoolingEfficiencyTable[(m+1)*Nrow + i];

      if ( Z == ZTable[4] )
         TEF_lambda[i] = LambdaByZ[4];
      else
      {
         TEF_lambda[i] = Mis_InterpolateFromTable( 5, ZTable, LambdaByZ, Z );

         if ( TEF_lambda[i] == NULL_REAL )
            Aux_Error( ERROR_INFO, "cooling efficiency interpolation failed at Z=%13.7e !!\n", Z );
      }
   }

   delete [] CoolingEfficiencyTable; 
}




//-----------------------------------------------------------------------------------------
// Function    :  Src_End_ExactCooling_General
// Description :  Free the resources used by a user-specified source term
//
// Note        :  1. Invoked by Src_End()
//                   --> Enable it by linking to the function pointer "Src_End_User_Ptr"
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_ExactCooling_General()
{


} // FUNCTION : Src_End_ExactCooling_General

#endif // #ifndef __CUDACC__

#endif // #ifdef EXACT_COOLING_GENERAL
