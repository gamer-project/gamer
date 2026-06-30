#include "GAMER.h"

#ifdef STAR_FORMATION


static RandomNumber_t *RNG = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar
// Description :  Interface for invoking various star-formation routines
//
// Note        :  1. Invoked by EvolveLevel()
//                2. Allocate NT thread-safe random number generators the first time this function is called,
//                   where NT is the number of OpenMP threads
//                   --> Must be freed by calling SF_FreeRNG() manually
//                3. Assign UIDs to new particles at the end of this routine
//
// Parameter   :  lv      : Target refinement level
//                TimeNew : Current physical time (after advancing solution by dt)
//                TimeOld : Physical time before update
//                          --> This function updates physical time from TimeOld to TimeNew
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void SF_CreateStar( const int lv, const real TimeNew, const real TimeOld, const real dt )
{

// only form stars on levels above the given minimum level
   if ( lv < SF_CREATE_STAR_MIN_LEVEL )   return;


// initialiez the random number generators the first time this function is called
   static bool FirstTime = true;

   if ( FirstTime )
   {
//    get the number of OpenMP threads
      int NT;
#     ifdef OPENMP
#     pragma omp parallel
#     pragma omp master
      {  NT = omp_get_num_threads();  }
#     else
      {  NT = 1;                      }
#     endif

//    allocate RNG
      RNG = new RandomNumber_t( NT );

//    make sure that different threads in different MPI ranks have different random seeds
//    --> even when SF_CREATE_STAR_DET_RANDOM is enabled, we still set random seeds here just in case that not all
//        star formation methods support SF_CREATE_STAR_DET_RANDOM
      for (int t=0; t<NT; t++)
      {
         const long RSeed = SF_CREATE_STAR_RSEED + MPI_Rank*1000 + t;
         RNG->SetSeed( t, RSeed );
      }

      FirstTime = false;
   }


// determine if metallicity is included
   const bool UseMetal = ( Idx_Metal != Idx_Undefined );


#  ifdef COMOVING
// convert dt from the comoving time interval to the physical time interval (2nd order)
// --> see Equation (15) of Schive, Tsai, & Chiueh (2010)
// --> it could be improved further by mimicking Miscellaneous/Mis_dTime2dt.cpp for the given TimeOld and TimeNew.
   const real dt_sf = (real)2.0*dt*SQR(TimeOld*TimeNew)/( SQR(TimeOld) + SQR(TimeNew) );
#  else
   const real dt_sf = dt;
#  endif


// invoke the target star-formation method
   switch ( SF_CREATE_STAR_SCHEME )
   {
#     if ( MODEL == HYDRO )
      case SF_CREATE_STAR_SCHEME_AGORA:
      case SF_CREATE_STAR_SCHEME_DWARFGALAXY:
         SF_CreateStar_GeneralGalaxy( lv, TimeNew, dt_sf, RNG, UseMetal );
         break;
#     endif

      case SF_CREATE_STAR_SCHEME_NONE:    return;
      break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SF_CREATE_STAR_SCHEME", SF_CREATE_STAR_SCHEME );
   } // switch ( SF_CREATE_STAR_SCHEME )

// assign UIDs for new particles
   Par_SetParUID();

} // FUNCTION : SF_CreateStar



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_FreeRNG
// Description :  Free the random number generator allocated by SF_CreateStar()
//
// Note        :  1. Invoked by End_MemFree()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SF_FreeRNG()
{

   delete RNG;
   RNG = NULL;

} // FUNCTION : SF_FreeRNG



#endif // #ifdef STAR_FORMATION
