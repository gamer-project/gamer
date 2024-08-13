#ifndef __MACRO__
#define __MACRO__

#include <stdint.h>

#ifndef NUM_THREADS
#  define NUM_THREADS 1
#endif

// single/double precision
#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif

//#ifdef FLOAT8
//# error: This code does not support double precision
//#endif

#  define ELECTRON_MASS         ( 9.10938356e-28         ) // g
#  define ELETRON_CHARGE        ( 4.80320425e-10         ) // statcoulombs (cm**1.5 g**0.5 s**-1)
#  define SPEED_OF_LIGHT        ( 29979245800.0          ) // cm/s
#  define MASS_PROTON_GRAM      ( 1.6726219e-24          ) // g
#  define PROTON_MASS_ENERGY    ( 9.3827208816e8         ) // eV
#  define BOLTZMANN_CONST_ERG   ( 1.38064852e-16         ) // erg/K
#  define BOLTZMANN_CONST_EV    ( 8.617333262145e-5      ) // eV/K
#  define REDUCED_PLANCK_EV     ( 6.582119569e-16        ) // eV s
#  define PLANCK_EV             ( 4.135667696e-15        ) // eV s
#  define ELECTRON_MASS_ENERGY  ( 510998.95              ) // eV
#  define THOMSON_CROSS_SECTION ( 6.652e-25              ) // cm**2
#  define ERG2EV                ( 6.2415e11              ) // 1 erg = 6.2415e11 eV
#  define EV2ERG                ( 1.6021789633902107e-12 ) // 1 / ERG2EV
#  define MILLIBARN_2_CM2       ( 1e-27                  ) // 1 mb = 1e-27 cm**-2

#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define   TAN( a )         tan( a )
#  define   POW( a , b )     pow( a , b )
#  define   LOG( a )         log( a )
#  define   EXP( a )         expf( a )
#  define   MPI_MYREAL       MPI_DOUBLE
#  define   EPSILON          DBL_EPSILON
#else
#  define  FABS( a )        fabsf( a )
#  define  SQRT( a )        sqrtf( a )
#  define   SIN( a )         sinf( a )
#  define   COS( a )         cosf( a )
#  define   TAN( a )         tanf( a )
#  define   POW( a , b )     powf( a , b )
#  define   LOG( a )         logf( a )
#  define   EXP( a )         exp( a )
#  define   MPI_MYREAL       MPI_FLOAT
#  define   EPSILON          FLT_EPSILON
#endif // #ifdef FLOAT8

#define PERSPECTIVE_PROJECTION
//#define SLICE


// distance between the sun and GC
#define R_SUN                    (-8.0) // kpc

// beta model
#define PEAK_DENS    ( 0.11  ) // cm**-3 // 0.46-0.35
#define CORE_RADIUS  ( 0.08  ) // kpc    // 0.35-0.27
#define BETA         ( 0.71  )

// Galacic halo
#define HALO_RADIUS  ( 250.0 ) // kpc
#define HALO_TEMP    ( 1.0e6 ) // Kelvien
#define MU_ELECTRON  ( 1.67  ) // molecular weight per electron in the Galactic halo

// Handy macro
#define STRING_LENGTH              50
#define MAX(a, b)             (  ( (a) > (b) ) ? (a) : (b)  )
#define SQR( x )              ( ( x ) * ( x ) )
#define CUBE( x )             ( ( x ) * ( x ) * ( x ) )
#define SIGN( a )             (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )

typedef struct Keys
{
   int UpperBound_ll;
   int UpperBound_bb;
#  ifdef PERSPECTIVE_PROJECTION
   real b_max;
   real b_min;
   real l_max;
   real l_min;
   real dt;
#  endif
   char AMRDataName[STRING_LENGTH];
   char FRBDataName[STRING_LENGTH];
   uint64_t EpochTimeStampInFRBData;
   int numAzimuthalAngle;
#  if ( defined HADRONIC_GAMMARAY || defined LEPTONIC_GAMMARAY )
   real scatteredEnergy;
#  elif ( defined SYNCHROTRON )
   real observedFreq;
#  endif
#  if ( defined HADRONIC_GAMMARAY || defined LEPTONIC_GAMMARAY || defined SYNCHROTRON )
   real gamma_max;
   real gamma_min;
   real spectral_index;
#  endif
#  ifdef SLICE
   char CuttingPlane[STRING_LENGTH];
#  endif
} Keys_t;


#define STR( x )                #x
#define SHOW_MACRO( x )         STR( x )

// print function
#define MASTER_PRINT( fmt, ... )     \
 {                                   \
    if ( MPI_Rank == 0 )             \
    {                                \
       printf( fmt, ##__VA_ARGS__ ); \
       fflush( stdout );             \
    }                                \
 }

#endif
