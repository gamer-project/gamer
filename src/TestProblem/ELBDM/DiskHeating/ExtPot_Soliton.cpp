#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY

extern double m_22;
extern double CoreRadius;
extern double Cen[3];


// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_Soliton
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_Soliton()
//
// Note        :  1. Invoked by Init_ExtPot_Soliton()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//                Time             : Target physical time
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_Soliton( double AuxArray_Flt[], int AuxArray_Int[], const double Time )
{

   const double A = 0.0019/m_22/m_22/POW(CoreRadius,4)*1.0e10*Const_Msun/Const_kpc; // soliton central density in g/cm^3
   const double B = 9.1*0.01;
   const double G = Const_NewtonG/UNIT_V/UNIT_V;   // convert output potential to code unit

   AuxArray_Flt[0] = Cen[0];  // x coordinate of the external potential center
   AuxArray_Flt[1] = Cen[1];  // y ...
   AuxArray_Flt[2] = Cen[2];  // z ...
   AuxArray_Flt[3] = CoreRadius*Const_kpc;
   AuxArray_Flt[4] = A;
   AuxArray_Flt[5] = B;
   AuxArray_Flt[6] = G;
   AuxArray_Flt[7] = UNIT_L;

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "EXT_POT_AUX_ARRAY:\n" );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  CenX       = %13.7e\n", Cen[0]     );
      Aux_Message( stdout, "  CenY       = %13.7e\n", Cen[1]     );
      Aux_Message( stdout, "  Cenz       = %13.7e\n", Cen[2]     );
      Aux_Message( stdout, "  CoreRadius = %13.7e\n", CoreRadius );
      Aux_Message( stdout, "  UNIT_L     = %13.7e\n", UNIT_L     );
      Aux_Message( stdout, "  UNIT_V     = %13.7e\n", UNIT_V     );
      Aux_Message( stdout, "  m_22       = %13.7e\n", m_22       );
      Aux_Message( stdout, "=============================================================================\n" );
   }

} // FUNCTION : SetExtPotAuxArray_Soliton
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_Soliton
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_SOliton(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = gravitational_constant*point_source_mass
//                3. Currently it does not support the soften length
//                4. GenePtr has the size of EXT_POT_NGENE_MAX defined in Macro.h (default = 6)
//
// Parameter   :  x/y/z             : Target spatial coordinates
//                Time              : Target physical time
//                UserArray_Flt/Int : User-provided floating-point/integer auxiliary arrays
//                Usage             : Different usages of external potential when computing total potential on level Lv
//                                    --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                        EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                        EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                                    --> This parameter is useless in most cases
//                PotTable          : 3D potential table used by EXT_POT_TABLE
//                GenePtr           : Array of pointers for general potential tables
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_Soliton( const double x, const double y, const double z, const double Time,
                            const double UserArray_Flt[], const int UserArray_Int[],
                            const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{

   const double Center[3] = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
   const double r_sol     = UserArray_Flt[3];
   const double dx        = (x - Center[0]);
   const double dy        = (y - Center[1]);
   const double dz        = (z - Center[2]);
   const double A         = UserArray_Flt[4];
   const double B         = UserArray_Flt[5];
   const double G         = UserArray_Flt[6];
   const double unit_l    = UserArray_Flt[7];
   const double r         = SQRT( dx*dx + dy*dy + dz*dz )*unit_l;
   const double Y         = r/r_sol;

// soliton potential, consistent with Eq. [6] in Chiang et al., PRD 103, 103019 (2021)
// --> this form may suffer from large numerical errors at r->0 as both the numerator and denominator approach 0
   double soliton_potential = -M_PI*pow(Y,2.)*A*G/53760./pow(B,1.5)
                               *( pow(B,0.5)/pow((1.+B*pow(Y,2.)),6.)
                                 *( 11895. + 36685.*B*pow(Y,2.)
                                  + 55638.*pow(B,2.)*pow(Y,4.) + 45738.*pow(B,3.)*pow(Y,6.)
                                  + 19635.*pow(B,4.)*pow(Y,8.) +  3465.*pow(B,5.)*pow(Y,10.) )
                                + 3465.*atan(pow(B,0.5)*Y)/Y );

   return (real)soliton_potential;

} // FUNCTION : ExtPot_Soliton



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_Soliton;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_Soliton
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_Soliton()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_Soliton( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_Soliton( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_Soliton( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_Soliton( double [], int [], const double );
void SetCPUExtPot_Soliton( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_Soliton( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_Soliton
// Description :  Initialize external potential
//
// Note        :  1. Set auxiliary arrays by invoking SetExtPotAuxArray_*()
//                   --> They will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_Soliton()
{

   SetExtPotAuxArray_Soliton( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int, Time[0] );
   SetCPUExtPot_Soliton( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_Soliton( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_Soliton

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
