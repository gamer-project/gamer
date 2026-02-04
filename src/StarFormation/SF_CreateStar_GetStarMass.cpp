#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )

static real SF_CreateStar_GetStarMass_StochasticLoaclSchmidtLaw( const real GasDens, const real dv, const real dt, RandomNumber_t *RNG,
                                                                 const real Efficiency, const real MinStarMass, const int TID );
static real SF_CreateStar_GetStarMass_MaxStarM( const real GasDens, const real dv, const real MaxStarMFrac );



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_GetStarMass
// Description :  Get the mass of star according to the adopted models
//
// Note        :  1. Return StarMass = -1.0 if the star is decided not to form (e.g., in the stochastic models )
//
// Parameter   :  GasDens        : Gas density
//                dv             : Cell volume
//                dt             : Time interval to advance solution
//                                 --> Currently this function does not distinguish dt and the physical time interval (dTime)
//                                 --> Does NOT support COMOVING yet
//                RNG            : Random number generator
//                TID            : OpenMP thread ID
//
// Return      :  StarMass
//-------------------------------------------------------------------------------------------------------
real SF_CreateStar_GetStarMass( const real GasDens, const real dv, const real dt, RandomNumber_t *RNG, const int TID )
{

   real StarMass = -1.0;

   switch ( SF_CREATE_STAR_SCHEME )
   {
      case SF_CREATE_STAR_SCHEME_AGORA:
      case SF_CREATE_STAR_SCHEME_DWARFGALAXY:
         StarMass = SF_CreateStar_GetStarMass_StochasticLoaclSchmidtLaw( GasDens, dv, dt, RNG, SF_CREATE_STAR_MASS_EFF, SF_CREATE_STAR_MIN_STAR_MASS, TID );
         break;

      case SF_CREATE_STAR_SCHEME_NONE:
         break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SF_CREATE_STAR_SCHEME", SF_CREATE_STAR_SCHEME );
   } // switch ( SF_CREATE_STAR_SCHEME )


   return StarMass;

} // FUNCTION : SF_CreateStar_GetStarMass



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_GetStarMass_StochasticLoaclSchmidtLaw
// Description :  Decide the new star particles stochastically using the presription suggested by the AGORA project
//
// Note        :  1. Ref: (1) Nathan Goldbaum, et al., 2015, ApJ, 814, 131 (arXiv: 1510.08458), sec. 2.4
//                        (2) Ji-hoon Kim, et al., 2016, ApJ, 833, 202 (arXiv: 1610.03066), sec. 3.2
//
// Parameter   :  GasDens        : Gas density
//                dv             : Cell volume
//                dt             : Time interval to advance solution
//                                 --> Currently this function does not distinguish dt and the physical time interval (dTime)
//                                 --> Does NOT support COMOVING yet
//                RNG            : Random number generator
//                Efficiency     : Gas-to-star mass efficiency                                    (--> "SF_CREATE_STAR_MASS_EFF"      )
//                MinStarMass    : Minimum star particle mass for the stochastical star formation (--> "SF_CREATE_STAR_MIN_STAR_MASS" )
//                TID            : OpenMP thread ID
//
// Return      :  StarMass
//-------------------------------------------------------------------------------------------------------
real SF_CreateStar_GetStarMass_StochasticLoaclSchmidtLaw( const real GasDens, const real dv, const real dt, RandomNumber_t *RNG,
                                                          const real Efficiency, const real MinStarMass, const int TID )
{

   real StarMass = -1.0;

   const real   Coeff_FreeFall = SQRT( (32.0*NEWTON_G)/(3.0*M_PI) );
   const real  _MinStarMass    = (real)1.0 / MinStarMass;
   const real   Eff_times_dt   = Efficiency*dt;


// 1. estimate the gas free-fall time
// --> consider only the gas density under the assumption that the dark matter doesn't collapse
   real _Time_FreeFall = Coeff_FreeFall * SQRT( GasDens );


// 2. estimate the gas mass fraction to convert to stars
   real StarMFrac = Eff_times_dt*_Time_FreeFall;
   real GasMass   = GasDens*dv;
   StarMass       = GasMass*StarMFrac;


// 3. stochastic star formation
// --> if the star particle mass (StarMass) is below the minimum mass (MinStarMass), we create a
//     new star particle with a mass of MinStarMass and a probability of StarMass/MinStarMass
// --> Eq. [5] in Goldbaum et al. (2015)
   if ( StarMass < MinStarMass )
   {
      const double Min = 0.0;
      const double Max = 1.0;

      double Random = RNG->GetValue( TID, Min, Max );

      if ( (real)Random < StarMass*_MinStarMass )  StarMass = MinStarMass;
      else                                         StarMass = -1.0;
   }


   return StarMass;

} // FUNCTION : SF_CreateStar_GetStarMass_StochasticLoaclSchmidtLaw



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_GetStarMass_MaxStarM
// Description :  Maximum gas mass allowed to be converted to star
//
// Note        :  1.
//
// Parameter   :  GasDens        : Gas density
//                dv             : Cell volume
//                MaxStarMFrac   : Minimum star particle mass for the stochastical star formation (--> "SF_CREATE_STAR_MAX_STAR_MFRAC" )
//
// Return      :  StarMass
//-------------------------------------------------------------------------------------------------------
real SF_CreateStar_GetStarMass_MaxStarM( const real GasDens, const real dv, const real MaxStarMFrac )
{

   real StarMass = -1.0;

   StarMass  = GasDens * dv * MaxStarMFrac;

   return StarMass;

} // FUNCTION : SF_CreateStar_GetStarMass_MaxStarM



#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )
