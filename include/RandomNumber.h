#ifndef __RANDOM_NUMBER_H__
#define __RANDOM_NUMBER_H__



#if ( RANDOM_NUMBER == RNG_CPP11 )
#include <random>

#elif ( RANDOM_NUMBER == RNG_GNU_EXT )
#include <stdlib.h>

#else
#error : ERROR : unsupported RANDOM_NUMBER !!
#endif

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  RandomNumber_t
// Description :  Data structure of the random number generator (RNG)
//
// Note        :  1. Adopted RNG implementation is controlled by the Makefile option "RANDOM_NUMBER"
//                   --> Currently supported options:
//                          RNG_CPP11  : c++11 library <random>
//                          RNG_GNU_EXT: GNU extension drand48_r
//                2. All implementations must be thread-safe
//
// Data Member :  RNG          : Random number generator
//                Distribution : Random number distribution used by RNG_CP11
//                N_RNG        : Total number of RNG
//
// Method      :  RandomNumber_t : Constructor
//               ~RandomNumber_t : Destructor
//                GetValue       : Return random number
//                SetSeed        : Set random seed
//-------------------------------------------------------------------------------------------------------
struct RandomNumber_t
{

// data members
// ===================================================================================
#  if   ( RANDOM_NUMBER == RNG_CPP11 )
   mt19937 *RNG;
   uniform_real_distribution <double> Distribution;

#  elif ( RANDOM_NUMBER == RNG_GNU_EXT )
   struct drand48_data *RNG;
#  endif

   int N_RNG;


   //===================================================================================
   // Constructor :  RandomNumber_t
   // Description :  Constructor of the structure "RandomNumber_t"
   //
   // Note        :  1. Allocate RNG
   //                2. Call SetSeed() to set the random seed for each RNG
   //
   // Parameter   :  N : Number of RNG to be initialized
   //                    --> Usually set equal to the number of OpenMP threads
   //===================================================================================
#  if   ( RANDOM_NUMBER == RNG_CPP11 )
   RandomNumber_t( const int N ) : Distribution( 0.0, 1.0 )
#  else
   RandomNumber_t( const int N )
#  endif
   {

//    check
      if ( N <= 0 )  Aux_Error( ERROR_INFO, "N (%d) <= 0 !!\n", N );

//    record the number of RNG
      N_RNG = N;

//    allocate RNG
#     if   ( RANDOM_NUMBER == RNG_CPP11 )
      RNG = new mt19937      [N];
#     elif ( RANDOM_NUMBER == RNG_GNU_EXT )
      RNG = new drand48_data [N];
#     endif

   } // METHOD : RandomNumber_t



   //===================================================================================
   // Destructor  :  ~RandomNumber_t
   // Description :  Destructor of the structure "RandomNumber_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~RandomNumber_t()
   {

      delete [] RNG;

   } // METHOD : ~RandomNumber_t



   //===================================================================================
   // Constructor :  GetValue
   // Description :  Return a uniformly distributed random number in the specified range
   //
   // Note        :  1. Only return a single random number
   //                   --> Must specify the ID of the target RNG
   //
   // Parameter   :  ID  : Target RNG (0 <= ID < N_RNG)
   //                Min : Lower limit of the random number
   //                Max : Upper limit of the random number
   //
   // Return      :  Random number
   //===================================================================================
   double GetValue( const int ID, const double Min, const double Max )
   {

//    check
#     ifdef GAMER_DEBUG
      if ( ID < 0  ||  ID >= N_RNG )
         Aux_Error( ERROR_INFO, "incorrect RNG ID = %d (total number of RNG = %d) !!\n", ID, N_RNG );
#     endif

//    get a uniformly distributed random number in the range [0.0, 1.0)
      double Random;
#     if   ( RANDOM_NUMBER == RNG_CPP11 )
      Random = Distribution( RNG[ID] );
#     elif ( RANDOM_NUMBER == RNG_GNU_EXT )
      drand48_r( RNG+ID, &Random );
#     endif

//    convert the range to [Min, Max) and return
      return Random*(Max-Min) + Min;

   } // METHOD : GetValue



   //===================================================================================
   // Constructor :  SetSeed
   // Description :  Set random seed for the target RNG
   //
   // Note        :  1. Only set random seed for a single RNG
   //                   --> Must specify the ID of the target RNG
   //
   // Parameter   :  ID   : Target RNG (0 <= ID < N_RNG)
   //                Seed : Random seed
   //===================================================================================
   void SetSeed( const int ID, const long Seed )
   {

//    check
#     ifdef GAMER_DEBUG
      if ( ID < 0  ||  ID >= N_RNG )
         Aux_Error( ERROR_INFO, "incorrect RNG ID = %d (total number of RNG = %d) !!\n", ID, N_RNG );
#     endif

//    set random seed
#     if   ( RANDOM_NUMBER == RNG_CPP11 )
      RNG[ID].seed( Seed );
#     elif ( RANDOM_NUMBER == RNG_GNU_EXT )
      srand48_r( Seed, RNG + ID );
#     endif

   } // METHOD : SetSeed


}; // struct RandomNumber_t



#endif // #ifndef __RANDOM_NUMBER_H__
