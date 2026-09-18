#ifndef __TURBULENCE_H__
#define __TURBULENCE_H__



#include "Macro.h"
#include <stdint.h>

//-------------------------------------------------------------------------------------------------------
// Structure   :  Turbulence_t
// Description :  Data structure of turbulence
//
// Data Member :  NMode            : Number of non-zero k modes
//                RNGState         : Current random state
//                Tdecay           : Turbulence correlation time
//                dt               : dt to update turbulence pattern
//                OUvar            : Ornstein-Uhlenbeck process variance
//                TimeLast         : Time for last turbulence field
//                TimeNext         : Time for next turbulence field
//                IdxLast          : Index for last turbulence field
//                IdxNext          : Index for next turbulence field
//                Amplitude        : Amplitude of each k mode
//                Kmode            : Non-zero k modes
//                Sin              : Array to store pre-computed sin modes
//                Cos              : Array to store pre-computed cos modes
//                OUphase          : Random phases updated by Ornstein-Uhlenbeck process, two entries for TimeLast and TimeNext
//
// Method      :  Turbulence_t     : Constructor
//               ~Turbulence_t     : Destructor
//-------------------------------------------------------------------------------------------------------
struct Turbulence_t
{

// data members
// ===================================================================================
   int     NMode;
   double  Tdecay;
   double  dt;
   double  OUvar;
   double  TimeLast;
   double  TimeNext;
   int     IdxLast;
   int     IdxNext;

   double *Amplitude;
   double *Kmode[3];
   double *Sin  [3];
   double *Cos  [3];
   double *OUphase[2];
   uint64_t RNGState;

// used for pcg32 rng process
   const uint64_t multiplier = 6364136223846793005ULL;
   const uint64_t increment  = 1442695040888963407ULL;

   //===================================================================================
   // Constructor :  Turbulence_t
   // Description :  Constructor of the structure "Turbulence_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Turbulence_t()
   {
      NMode    =  0;
      Tdecay   =  1.0;
      dt       =  1.0;
      OUvar    =  1.0;
      TimeLast = -1.0;
      TimeNext = -1.0;
      IdxLast  =  0;
      IdxNext  =  1;
      RNGState =  (uint64_t) 123;

      Amplitude  = NULL;

      for (int i=0; i<3; i++)
      {
         Kmode[i] = NULL;
         Sin  [i] = NULL;
         Cos  [i] = NULL;
      }

      for (int i=0; i<2; i++)
         OUphase[i] = NULL;

   } // METHOD : Turbulence_t



   //===================================================================================
   // Destructor  :  ~Turbulence_t
   // Description :  Destructor of the structure "Turbulence_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Turbulence_t()
   {
      if ( Amplitude != NULL ) delete [] Amplitude;

      for (int i=0; i<3; i++)
      {
         if ( Kmode[i] != NULL ) delete [] Kmode[i];
         if ( Sin  [i] != NULL ) delete [] Sin  [i];
         if ( Cos  [i] != NULL ) delete [] Cos  [i];
      }

      for (int i=0; i<2; i++)
         if ( OUphase[i] != NULL ) delete [] OUphase[i];

   } // METHOD : ~Turbulence_t

   void SetRNGState(int seed)
   {
      RNGState = (uint64_t)seed + increment;
   }

   uint32_t rotr32(uint32_t x, unsigned r)
   {
      return x >> r | x << (-r & 31);
   }

// uniformly distributed random number between (0, 1)
   double uniform( uint64_t& state )
   {
//    PCG-XSH-RR rng sequence
      uint64_t x = state;
      unsigned count = (unsigned)(x >> 59);

      state = x * multiplier + increment;
      x ^= x >> 18;
      uint32_t result = rotr32( (uint32_t)(x >> 27), count );

      return ( (double)result + 0.5 ) / 4294967296.0;
   }

// get random pair using Box–Muller transformation
   void GetRNG( double& a, double& b )
   {
      double r1 = uniform( RNGState );
      double r2 = uniform( RNGState );

      double mag = OUvar*sqrt(2.0*log(1.0/r1));
      double phi = 2*M_PI*r2;

      a = mag*cos( phi );
      b = mag*sin( phi );
   }

}; // struct Turbulence_t



#endif // #ifndef __TURBULENCE_H__
