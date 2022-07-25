#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT



// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__

extern double Bondi_MassBH;
extern double Bondi_Soften_R;
extern bool   Bondi_Soliton;
extern double Bondi_Soliton_m22;
extern double Bondi_Soliton_rc;
extern int    Bondi_Soliton_type;
extern double Bondi_Soliton_t;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_Bondi
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_Bondi()
//
// Note        :  1. Invoked by Init_ExtAcc_Bondi()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//                Time     : Target physical time
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_Bondi( double AuxArray[], const double Time )
{

   AuxArray[0] = amr->BoxCenter[0];
   AuxArray[1] = amr->BoxCenter[1];
   AuxArray[2] = amr->BoxCenter[2];
   AuxArray[3] = NEWTON_G*Bondi_MassBH;   // gravitational_constant*point_source_mass
   AuxArray[4] = Bondi_Soften_R;          // soften_length (<=0.0 --> disable)
   AuxArray[5] = Bondi_Soliton;
   if( Bondi_Soliton )
   {
      AuxArray[6] = Bondi_Soliton_m22;
      AuxArray[7] = Bondi_Soliton_rc;
      AuxArray[8] = NEWTON_G*Const_Msun/UNIT_M;
      AuxArray[9] = UNIT_L/Const_kpc;
      if( Bondi_Soliton_type > 0 )
      {
          switch(Bondi_Soliton_type)
          {
              case 1: // arctan function
                  AuxArray[8] *= 2/(real)3.14159265*ATAN(Time/Bondi_Soliton_t);
                  break;
              case 2: // linear function
                  if( Time < Bondi_Soliton_t )
                      AuxArray[8] *= Time/Bondi_Soliton_t;
                  break;
              case 3: // smooth step function
                  if( Time < Bondi_Soliton_t )
                      AuxArray[8] *= 3*SQR(Time/Bondi_Soliton_t)-2*CUBE(Time/Bondi_Soliton_t);
                  break;
              case 4: // sigmoid
                  AuxArray[8] *= 1/(1+exp(-Time+Bondi_Soliton_t/2));
                  break;
                
          }
      } 
   }

} // FUNCTION : SetExtAccAuxArray_Bondi
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_Bondi
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_Bondi()
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_Bondi( real Acc[], const double x, const double y, const double z, const double Time,
                          const double UserArray[] )
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
         real GM       = (real)UserArray[3];
   const real eps      = (real)UserArray[4];
   const bool SOL      = (real)UserArray[5];
   const real dx       = (real)(x - Cen[0]);
   const real dy       = (real)(y - Cen[1]);
   const real dz       = (real)(z - Cen[2]);
   const real r        = SQRT( dx*dx + dy*dy + dz*dz );

   if( SOL )
   {
      const real m22      = (real)UserArray[6];
      const real rc       = (real)UserArray[7];  // In code unit or kpc
      const real Coeff    = (real)UserArray[8];
      const real UNIT_L   = (real)UserArray[9];

#ifdef Plummer
      double M = GM*CUBE(r)/pow(SQR(r)+SQR(rc),1.5);
#else
      real a = SQRT(POW(2.0,1.0/8.0)-1)*(r/rc);
      real M = (real)4.17e9/(SQR(m22/1e-1)*(rc*UNIT_L*1e3)*POW(SQR(a)+1, 7.0))*((real)3465*POW(a,13.0)+(real)23100*POW(a,11.0)+(real)65373*POW(a,9.0)+(real)101376*POW(a,7.0)+(real)92323*POW(a,5.0)+(real)48580*POW(a,3.0)-(real)3465*a+(real)3465*POW(SQR(a)+1, 7.0)*ATAN(a));
      M *= Coeff;
#endif
      GM += M;
   }

// Plummer
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif

   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;

} // FUNCTION : ExtAcc_Bondi



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_Bondi;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_Bondi
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_Bondi()
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_Bondi( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_Bondi( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_Bondi( double [], const double );
void SetCPUExtAcc_Bondi( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_Bondi( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_Bondi
// Description :  Initialize external acceleration
//
// Note        :  1. Set an auxiliary array by invoking SetExtAccAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external acceleration major routines by invoking SetCPU/GPUExtAcc_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtAcc_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtAcc_Bondi()
{

   SetExtAccAuxArray_Bondi( ExtAcc_AuxArray, Time[0] );

   SetCPUExtAcc_Bondi( CPUExtAcc_Ptr );
#  ifdef GPU
   SetGPUExtAcc_Bondi( GPUExtAcc_Ptr );
#  endif

} // FUNCTION : Init_ExtAcc_Bondi

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
