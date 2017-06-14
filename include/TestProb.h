#ifndef __TESTPROB_H__
#define __TESTPROB_H__



// common headers
#include "ReadPara.h"


// common function prototypes
static void Validate();
static void SetParameter();
static void SetGridIC( real *fluid, const double x, const double y, const double z, const double Time );


// function pointers of various user-specified routines
extern void (*Init_Function_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_User_Ptr)();
extern bool (*Flag_User_Ptr)( const int i, const int j, const int k, const int lv, const int PID, const double Threshold );
extern void (*Mis_GetTimeStep_User_Ptr)( double &dt, double &dTime, const double dt_dTime );


// handy constants for loading the runtime parameters
const char   Useless_str   = '\0';
const bool   Useless_bool  = false;
const int    NoMax_int     = (int   )+__INT_MAX__;
const long   NoMax_long    = (long  )+__INT_MAX__;
const uint   NoMax_uint    = (uint  )+__INT_MAX__;
const ulong  NoMax_ulong   = (ulong )+__INT_MAX__;
const float  NoMax_float   = (float )+__FLT_MAX__;
const double NoMax_double  = (double)+__FLT_MAX__;
const int    NoMin_int     = (int   )-__INT_MAX__;
const long   NoMin_long    = (long  )-__INT_MAX__;
const uint   NoMin_uint    = (uint  )0;
const ulong  NoMin_ulong   = (ulong )0;
const float  NoMin_float   = (float )-__FLT_MAX__;
const double NoMin_double  = (double)-__FLT_MAX__;
const float  NoZero_float  = (float )__FLT_MIN__;
const double NoZero_double = (double)__FLT_MIN__;


// helper macro for printing warning messages when resetting parameters
#  define FORMAT_INT       %- 21d
#  define FORMAT_LONG      %- 21ld
#  define FORMAT_UINT      %- 21u
#  define FORMAT_ULONG     %- 21lu
#  define FORMAT_BOOL      %- 21d
#  define FORMAT_REAL      %- 21.14e
#  define PRINT_WARNING( name, var, format )                                                             \
   {                                                                                                     \
      if ( MPI_Rank == 0 )                                                                               \
         Aux_Message( stderr, "WARNING : parameter [%-25s] is reset to [" EXPAND_AND_QUOTE(format) "] "  \
                              "for the adopted test problem\n", name, var );                             \
   }



#endif // #ifndef __TESTPROB_H__
