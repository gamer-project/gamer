#include "GAMER.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <stdio.h>

using std::cout; using std::endl;
using std::vector; using std::copy;
using std::string;

#ifdef FEEDBACK



// function pointers to be set by FB_Init_User_Template()
extern void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                            real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[],
                            const int TID, RandomNumber_t *RNG );
extern void (*FB_End_User_Ptr)();


// internal funcitons
real UV_IonizingLuminosity( const double m );

int FB_GetNumberOfSNe( const real MinimumMass, RandomNumber_t *RNG, const int TID /*, const RandomSeed*/ );

real FB_WindEjectMass( const real age );

real FB_WindEjectEnergy( const real age );

real FB_SNeEjectMass( const real age );

void FB_UVFeedback( float (*Fluid)[PS2][PS2][PS2], const real particleIonizingLuminosity, 
                    const int pararticleWillExplode, const int idx[], const double dh );

void FB_WindFeedback( const real particleAge, real &windEjectMass, real &windEjectEnergy, 
                      const real particleSoonestExplosion, real &coldWind, const real minimumStarMass,
                      const double dt );

void FB_SNMetalFeedback( real (*Fluid)[PS2][PS2][PS2], const real windEjectDensity, 
                         const real SNeEjectDensity, const real particleAge, const int explosionFlag, 
                         const int imetal, const real metalf, const int idx[],  
                         const int imetalSNII, const real metalfSNIIn);

void FB_SNeMomentumFeedback( real (*Fluid)[PS2][PS2][PS2], const int particleExplosionFlag, 
                             const real windEjectDensity, 
                             const real SNeEjectDensity, const real SNeEjectEnergy,
                             const int imethod, const int idual, 
                             const int distrad, const int diststep, const int distcells,
                             const double EdgeL[], const int idx[],  const double dh );

//-------------------------------------------------------------------------------------------------------
// Function    :  FB_User_Template
// Description :  Template of a user-defined feedback
//
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAtt[], respectively
//                2. Must use ParSortID[] to access ParAtt[]
//                   --> ParAtt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
//                3. Particles may be outside the target region
//                4. To ensure the consistency of random numbers, one must call the random number generator for
//                   ALL particles, including those too far away to affect the target region
//                5. No need to worry about the periodic boundary condition here
//                   --> Particle positions have been remapped in FB_AdvanceDt()
//                6. CoarseFine[] records the coarse-fine boundaries along the 26 sibling directions, defined as
//                            24  13  25
//                            15  05  17     z+1 plane
//                            22  12  23
//
//                            08  03  09
//                            00  XX  01     z   plane
//                            06  02  07
//                   y
//                   ^        20  11  21
//                   |        14  04  16     z-1 plane
//                   --->x    18  10  19
//                7. Invoked by FB_AdvanceDt()
//                8. Must NOT change particle positions
//                9. Linked to FB_User_Ptr in FB_Init_User_Template()
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                NPar       : Number of particles
//                ParSortID  : Sorted particle IDs
//                ParAtt     : Particle attribute arrays
//                Fluid      : Array to store the input/output fluid data
//                             --> Array size is fixed to PS2^3
//                EdgeL      : Left edge of Fluid[]
//                             --> Right edge is given by EdgeL[]+PS2*dh
//                dh         : Cell size of Fluid[]
//                CoarseFine : Coarse-fine boundaries along the 26 sibling directions
//                TID        : Thread ID
//                RNG        : Random number generator
//                             --> Random number can be obtained by "RNG->GetValue( TID, Min, Max )",
//                                 where Min/Max specify the range of random numbers
//
// Return      :  Fluid, ParAtt
//-------------------------------------------------------------------------------------------------------


void FB_User_Template( const int lv, const double TimeNew, const double TimeOld, const double dt,
                       const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                       real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[],
                       const int TID, RandomNumber_t *RNG )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid == NULL )    Aux_Error( ERROR_INFO, "Fluid == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAtt == NULL )      Aux_Error( ERROR_INFO, "ParAtt == NULL for NPar = %d !!\n", NPar );
   }
   
   printf("feedback calculation is called.\n");
   fflush( stdout );

   // if particle in grid check (SNHandeler.C)
   
#  endif // #ifdef GAMER_DEBUG

// Needed for determine SNe number in a particle
   int *explosionFlag      = new int[NPar];
   int *explosionFlagIa    = new int[NPar];
   int *willExplode        = new int[NPar];
   real *soonestExplosion  = new double[NPar];
   real *ionizeLuminosity  = new double[NPar];
   int rand_int		   = 0;

// Polynomial coefficients for the fit to the delay time distribution
   const double p_delay[6] =
   {
        3505021.4516666,
       16621326.48066255,
        4382816.59085307,
       46194173.14420852,
      -52836941.28241706,
       22967062.02780452,
   }; // const double p_delay[6]

// Polynomial coefficients for the fit to the stellar mass distribution as a function of delay time.
   const double p_mass[10] = 
   {
        4.42035634891,
      -13.2890466089,
       26.1103296098,
      -30.1876007562,
       21.8976126631,
      -10.2544493943,
        3.09621304958,
       -0.581870413299,
        0.0618795946119,
       -0.00284352366041
   }; // const double p_mass[10] 


   int n = 0;

   for ( n = 0; n < NPar; n++ )
   {
      // constants 
      const real _dh         = 1.0 / dh;
      const double Msun = 1.9885e33;
      const double MinParMass = 3.2e3 * Msun;

      // Particle variables
      const int par_idx        = ParSortID[n];
      const double par_pos[3]  = { ParAtt[PAR_POSX][par_idx], 
                                   ParAtt[PAR_POSY][par_idx], 
                                   ParAtt[PAR_POSZ][par_idx] };
      const double par_vel[3]  = { ParAtt[PAR_VELX][par_idx], 
                                   ParAtt[PAR_VELY][par_idx],
                                   ParAtt[PAR_VELZ][par_idx] };
      const double par_creTime = ParAtt[Idx_ParCreTime][par_idx];          // creation time of particle
      const double par_age     = ( TimeOld-par_creTime )*UNIT_T/3.16e14;   // age of a particle
      printf("The time now is %lf.\n ", par_age);
      double par_mass          = ParAtt[PAR_MASS][par_idx];
      
      
      // Fluid variables
      int idx[3];    // define the coordinate of a particle_consisting cell
      for ( int d = 0; d < 3; d++ ) idx[d] = (int)FLOOR( ( par_pos[d] - EdgeL[d] )*_dh );
      real flu_dens      = Fluid[DENS][idx[0]][idx[1]][idx[2]];      // density of this cell
      real flu_energy    = Fluid[ENGY][idx[0]][idx[1]][idx[2]];      // energy of this cell
      
      
      // Feedback variables
      real fb_energy     = 0.0;                    // total energy feedback of particle
      
      // 1: UV feedback
      real ionized;
      
      // 2: Wind feedback
      real wind_mass     = 0.0;                    // total mass feedback by wind
      real wind_energy   = 0.0;                    // total energy feedback by wind
      real coldWind      = 0.0;                    // wind of old particle
      
      // 3: SNe feedback
      real SNe_mass      = 0.0;
      real SNe_dens      = 0.0;
      real wind_dens     = 0.0;
      real sn_energy     = 0.0;                    // energy emitted by supernova
      setbuf(stdout, NULL);
      printf("The supernova mass feedback = %f\n", SNe_mass);
      printf("The supernova energy feedback = %f\n", sn_energy);
 
      // 4: Metal feedback
      real ZSN           = 0.0;
      
      // 5: Momentum feedback
      real p_sn          = 0.0;                    // momentum to deposit per SNe in code units
      
      
      
      // ============================================================
      // Check the particle and the patch 
      // ============================================================
      // 1: The particle is too old, so there is no feedback
      if ( par_age >= 3.7e7*Const_yr/UNIT_T )       continue;
      
      // TODO: There is always particle outside the patch. This check might not be good for GAMER
      // 2. check bounds - if star particle is outside of this grid
      // then something is very wrong and we should crash
      // TODO : Not defined nx, ny, nz
      //if (i < 1 || i > nx || j < 1 || j > ny || k < 1 || k > nz) 
      //{
      //   Aux_Error( ERROR_INFO, "star particle out of grid !!\n" );
      //}
      
      // 3: The particle creation time 
      if ( par_creTime > 0.0 )        continue;
      
      // 4: The particle mass
      if ( par_mass > 0.0 )         continue;
      
      // TODO: down below
      // 5: The level of the patch
      
      
      
      // ============================================================
      // Calculate the nummber of the feedback events in particle
      // ============================================================
      soonestExplosion[n] = -1.0;
      
      // Seed the RNG for this particle.
      // This allows us to reproduceably calculate the SN properties of each particle each timestep.
      // TODO : No funciton called init in GAMER.
      //no need to initialize here
      //init( par_idx );
      
      explosionFlag[n] = -1;
      willExplode[n]   =  0;
      
      // figure out how many explosions this particle undergoes.
      // TODO: No variable called MinStarMass 
      const int number_of_sn = FB_GetNumberOfSNe( MinParMass, RNG, TID );
      
      // If there are no SN in this particle, go on to the next one.
      if ( number_of_sn == 0 )
      {
         explosionFlag[n] = 0;
         continue; // Go to next particle
      } // if ( number_of_sn == 0 )



      // ============================================================
      // Calculate physics feedback parameter
      // ============================================================
      // If there are any explosions, we need to know the ionizing luminosity of each progenitor
      // and whether or not that progenitor will explode soon or in this timestep.
      
      for ( int kk = 0; kk < number_of_sn; ++kk )     // Loop over expected SN.
      {
         // Draw the delay time for this SN event.
         // TODO : Use RNG->GetValue( TID, Min, Max ) in GAMER.
         // rand_int = random();
         const real x = RNG->GetValue( TID, 0, 32767.0 ) / 32768.0;
         real delayTime =
               p_delay[0] +
               p_delay[1]*x +
               p_delay[2]*x*x +
               p_delay[3]*x*x*x +
               p_delay[4]*x*x*x*x +
               p_delay[5]*x*x*x*x*x; // in years.
      
         // Convert the delay time to a progenitor mass.
         const real td7 = delayTime*1e7;
         double progenitor_mass =
               p_mass[0] +
               p_mass[1]*td7 +
               p_mass[2]*td7*td7 +
               p_mass[3]*td7*td7*td7 +
               p_mass[4]*td7*td7*td7*td7 +
               p_mass[5]*td7*td7*td7*td7*td7 +
               p_mass[6]*td7*td7*td7*td7*td7*td7 +
               p_mass[7]*td7*td7*td7*td7*td7*td7*td7 +
               p_mass[8]*td7*td7*td7*td7*td7*td7*td7*td7 +
               p_mass[9]*td7*td7*td7*td7*td7*td7*td7*td7*td7;
         progenitor_mass = POW( 10, progenitor_mass ); // in solar masses
      
         delayTime *= Const_yr/UNIT_T; // convert years -> seconds -> code units.
         const real relativeTime = TimeOld - ParAtt[Idx_ParCreTime][par_idx]; // in simulation units.
      
         if ( delayTime > relativeTime && delayTime < relativeTime+dt )
         {
            if ( explosionFlag[n] == -1 )
            {
               explosionFlag[n] = 1;
            } // if ( explosionFlag[n] == -1 )
            else if ( explosionFlag[n] == 0 )
            {
               Aux_Error( ERROR_INFO, "I thought there were no explosions but it turns out there's one right now");
            } // if ( explosionFlag[n] == -1 ) ... else if ...
            else if ( explosionFlag[n] > 0 )
            {
               explosionFlag[n] += 1;
            } // if ( explosionFlag[n] == -1 ) ... else if ... else if ...
            else 
            {
               // other behaviors possible here
               // e.g. ignore this new SN, and/or print a warning.
               Aux_Error( ERROR_INFO, "Weird value of explosionFlag");
            } // if ( explosionFlag[n] == -1 ) ... else if ... else if ... else ...
          
      	} // if ( delayTime > relativeTime && delayTime < relativeTime+dtFixed )

         if ( relativeTime < delayTime  )
         {
            // SN hasn't gone off yet, so get ionizing luminosity and add it to the total for this particle
            ionizeLuminosity[n] += UV_IonizingLuminosity( progenitor_mass );  // in 10^49 photons s^-1
            willExplode[n] = 1;
               
            if ( soonestExplosion[n] < 0 || delayTime < soonestExplosion[n] )
            {
      	      soonestExplosion[n] = delayTime;
      	   } // if ( soonestExplosion[n]<0 || delayTime<soonestExplosion[n])
          
      	} // if( relativeTime < delayTime  )
      
         if ( explosionFlag[n] == -1 )  explosionFlag[n] = 0;

      } // for ( int kk=0; kk < number_of_sn; ++kk )




      // ============================================================
      // Feedback part
      // ============================================================
      // 0: Presetting
      real (*fluid_PG)[PS2][PS2][PS2] = new real [NCOMP_TOTAL][PS2][PS2][PS2];
      memcpy( fluid_PG, Fluid, sizeof(Fluid) );
      
      // 1: UV feedback
//      FB_UVFeedback( fluid_PG, ionizeLuminosity[n], willExplode[n], idx, dh);
//      memcpy( Fluid, fluid_PG, sizeof(*fluid_PG) );
      
      // 2: Wind feedback
//      FB_WindFeedback( par_age, wind_mass, wind_energy, soonestExplosion[n], coldWind, MinParMass, dt );
      memcpy( Fluid, fluid_PG, sizeof(*fluid_PG) );
      
      // 3: SN energy feedback
      // evolution of SN ejection by time
//      if ( explosionFlag[n] > 0 )       SNe_mass = FB_SNeEjectMass( par_age );

//      SNe_mass  = POW(10.0, SNe_mass) * explosionFlag[n];

      //record SNe_mass
      printf("The supernova = %f\n", SNe_mass);

      SNe_dens  = SNe_mass * Msun / CUBE( dh*UNIT_L ) / UNIT_D;
      wind_dens = wind_mass * Msun / CUBE( dh*UNIT_L ) / UNIT_D;

      if ( (SNe_dens + wind_dens) > par_mass )
      {
         printf("WARNING: losing too much mass\n");
         printf("Setting particle mass to zero\n");
         par_mass = 0;
      }
      else // if ((SNe_dens + wind_dens) > par_mass)
      {
         par_mass = par_mass - SNe_dens - wind_dens;
      } // if ((SNe_dens + wind_dens) > par_mass) ... else ...

      // update mass change of particles
      ParAtt[PAR_MASS][par_idx] = par_mass;

      if (explosionFlag[n] > 0)
      {
         sn_energy = explosionFlag[n] * 
                     ( 1.0e51 / ( UNIT_D * CUBE( UNIT_L ) ) / UNIT_V / UNIT_V );
      }
      else // if (explosionFlag[n] > 0)
      {
         sn_energy = 0;
      } // if (explosionFlag[n] > 0) ... else ...


      // calculate energy (wind will be included later)
      fb_energy = sn_energy + wind_energy * wind_dens * CUBE( dh );
      fb_energy = fb_energy / ( (flu_dens + wind_dens + SNe_dens) * CUBE( dh ) );

      // add energy into fluid
      Fluid[ENGY][idx[0]][idx[1]][idx[2]] += fb_energy;



      // 4: SN Metal feedback (in metal fraction)
      // no metal feedback now
      // TODO
      // FB_SNMetalFeedback( Fluid, wind_energy, sn_energy, par_age, 
      //                     imetal, metalf, imetalSNII, metalfSNIIn );
      
      
      // Mass and momentum feedback from winds and SN ejecta (NOT from SN blastwaves)
      Fluid[MOMX][idx[0]][idx[1]][idx[2]] += (wind_dens + SNe_dens) * par_vel[0];
      Fluid[MOMY][idx[0]][idx[1]][idx[2]] += (wind_dens + SNe_dens) * par_vel[1];
      Fluid[MOMZ][idx[0]][idx[1]][idx[2]] += (wind_dens + SNe_dens) * par_vel[2];
      Fluid[DENS][idx[0]][idx[1]][idx[2]] +=  wind_dens + SNe_dens;
      
      // TODO: dual energy fixed
      #ifdef DUAL_ENERGY
      Fluid[ENPY][idx[0]][idx[1]][idx[2]] += fb_energy;
      #endif // #ifdef DUAL_ENERGY

      // 5: SN Momentum feedback
      // no momentum feedback now
      // TODO
      // FB_SNeMomentumFeedback( Fluid, explosionFlag[n], wind_dens, SNe_dens, sn_energy, imethod,  idual, 
      //                         distrad, diststep, distcells, dh )
      
   
   
   
   } // for (n = 0; n < NPar; n++ )

//   delete[] ionizeLuminosity;
//   delete[] explosionFlag;
//   delete[] explosionFlagIa;
//   delete[] willExplode;
//   delete[] soonestExplosion;

} // FUNCTION : FB_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_User_Template
// Description :  Free the resources used by the user-specified feedback
//
// Note        :  1. Invoked by FB_End()
//                2. Linked to FB_End_User_Ptr in FB_Init_User_Template()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_User_Template()
{


} // FUNCTION : FB_End_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_User_Template
// Description :  Initialize the user-specified feedback
//
// Note        :  1. Invoked by FB_Init()
//                   --> Enable it by linking to the function pointer "FB_Init_User_Ptr"
//                2. Set FB_User_Ptr and FB_End_User_Ptr
//
// Parameter   :  None
//
// Return      :  FB_User_Ptr and FB_End_User_Ptr
//-------------------------------------------------------------------------------------------------------
void FB_Init_User_Template()
{

   FB_User_Ptr     = FB_User_Template;
   FB_End_User_Ptr = FB_End_User_Template;

} // FUNCTION : FB_Init_User_Template

/*
//-------------------------------------------------------------------------------------------------------
// Function    :  UV_IonizingLuminosity
// Description :  Returns a broken power-law approximation to the ionizing luminosity emitted by 
//                a main sequence star of mass m, normalized to 10^49 photon/s
//
// Note        :
//
// Parameter   :  m : Mass of a main sequence star.
//
// Return      :  ..........
//
// Reference   :  Parravano et al. 2003 ApJ 584 797 (doi:10.1086/345807)
//-------------------------------------------------------------------------------------------------------
real UV_IonizingLuminosity( const real m )
{
   if ( m < 5 )         return 0;
   else if ( m < 7  )   return( 2.23e-15 * POW(m, 11.5) );
   else if ( m < 12 )   return( 3.69e-13 * POW(m, 8.87) );
   else if ( m < 20 )   return( 4.8e-12  * POW(m, 7.85) );
   else if ( m < 30 )   return( 3.12e-8  * POW(m, 4.91) );
   else if ( m < 40 )   return( 2.8e-5   * POW(m, 2.91) );
   else if ( m < 60 )   return( 3.49e-4  * POW(m, 2.23) );
   else                 return( 2.39e-3  * POW(m, 1.76) );
} // FUNCTION :  UV_IonizingLuminosity


*/

//-------------------------------------------------------------------------------------------------------
// Function    : FB_GetNumberOfSNe
// Description : 
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : http://goo.gl/sgLPcj
//-------------------------------------------------------------------------------------------------------
/*
int FB_GetNumberOfSNe( const real MinimumMass, RandomNumber_t *RNG, const int TID , const RandomSeed )
{
   // This is the number of Supernovae expected from a 10^6 Msun cluster, 
   // as predicted from Starburst99.
   const float expected_sn_s99 = 10616.955572;
   // lambda is the expected number of SN for a stellar population with mass MinimumMass
   const float lambda = expected_sn_s99 * 1.0e-6 * MinimumMass;
   const real L = EXP( -lambda );
   int k = 0;
   real p = 1.0;
   real u;
   
   // The Knuth algorithm for drawing from a Poisson Distribution.
   while( p > L )
   {
      ++k;
      // TODO : Use RNG->GetValue( TID, Min, Max ) in GAMER.
      int rand_int = RNG->GetValue( TID, 0.0, 32767.0 );
      u = rand_int / 32768.0;
      p *= u;
   } // while( p > L )
   
   // Now k-1 ~ Pois(lambda)
   k -= 1;

   return k;
} // FUNCTION : void FB_GetNumberOfSNe



//-------------------------------------------------------------------------------------------------------
// Function    : FB_WindEjectMass 
// Description :  
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
real FB_WindEjectMass( const real age )
{
   real fb_mass = 0.0;
   if ( age <= 0.29 )
   {
      fb_mass -= 2.72379457924             ;
      fb_mass += 9.03549102928 * age       ;
      fb_mass -= 349.073894935 * POW(age,2);
      fb_mass += 6133.48804337 * POW(age,3);
      fb_mass -= 45526.5891824 * POW(age,4);
      fb_mass += 160159.422053 * POW(age,5);
      fb_mass -= 254942.557778 * POW(age,6);
      fb_mass += 133227.581992 * POW(age,7);
   }
   else if ( age < 0.29 && age <= 0.447)
   {
      fb_mass += 1024395.67006             ;
      fb_mass -= 22826983.8411 * age       ;
      fb_mass += 221695566.585 * POW(age,2);
      fb_mass -= 1225731636.28 * POW(age,3);
      fb_mass += 4219889881.74 * POW(age,4);
      fb_mass -= 9263931021.04 * POW(age,5);
      fb_mass += 12664879027.3 * POW(age,6);
      fb_mass -= 9858823353.65 * POW(age,7);
      fb_mass += 3345877714.47 * POW(age,8);
   }
   else if ( age >= 0.447 && age <= 3.24)
   {
      fb_mass -= 69.9540656568               ;
      fb_mass += 540.489990705  * age        ;
      fb_mass -= 1785.45996729  * POW(age,2) ;
      fb_mass += 3267.21333796  * POW(age,3) ;
      fb_mass -= 3721.01017711  * POW(age,4) ;
      fb_mass += 2776.66619261  * POW(age,5) ;
      fb_mass -= 1381.69750895  * POW(age,6) ;
      fb_mass += 454.445793525  * POW(age,7) ;
      fb_mass -= 94.8530920783  * POW(age,8) ;
      fb_mass += 11.3775633348  * POW(age,9) ;
      fb_mass -= 0.597112042033 * POW(age,10);
   }
   else if ( age >= 3.24 && age <= 7.28)
   {
      fb_mass += 14.8304028947;
      fb_mass -= 16.0726787694   * age;
      fb_mass += 4.55554172673   * POW(age,2);
      fb_mass -= 0.522521195039  * POW(age,3);
      fb_mass += 0.0214246750892 * POW(age,4);
   }
   else if ( age >= 7.28 && age <= 7.8)
   {
      fb_mass += 3533.50268935             ;
      fb_mass -= 1430.45985176 * age       ;
      fb_mass += 192.934935096 * POW(age,2);
      fb_mass -= 8.67530777079 * POW(age,3);
   }
   else if (age >= 7.8)
   {
      fb_mass -= 3.5441266853         ;
      fb_mass -= 0.0101229134154 * age;
   } // if ( age  < 0.29 )
   return fb_mass;

} // FUNCTUION : FB_WindEjectMass



//-------------------------------------------------------------------------------------------------------
// Function    : FB_WindEjectEnergy 
// Description :  
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
real FB_WindEjectEnergy( const real age )
{
   real fb_energy = 0.0;
   if ( age < 0.29 )
   {
      fb_energy += 16.7458830136             ;
      fb_energy += 2.87170884625 * age       ;
      fb_energy -= 260.188160495 * POW(age,2);
      fb_energy += 7588.41970548 * POW(age,3);
      fb_energy -= 109128.119673 * POW(age,4);
      fb_energy += 834565.297424 * POW(age,5);
      fb_energy -= 3488638.06781 * POW(age,6);
      fb_energy += 7534432.3913  * POW(age,7);
      fb_energy -= 6577304.157   * POW(age,8);
   }
   else if ( age >= 0.29 && age <= 0.518 )
   {
      fb_energy -= 8756108.90226                 ;
      fb_energy += 249183578.797    * age        ;
      fb_energy -= 3210919258.32    * POW(age,2) ;
      fb_energy += 24729596291.2    * POW(age,3) ;
      fb_energy -= 126485776877.0   * POW(age,4) ;
      fb_energy += 451123199767.0   * POW(age,5) ;
      fb_energy -= 1.14486556168e12 * POW(age,6) ;
      fb_energy += 2.067395301e12   * POW(age,7) ;
      fb_energy -= 2.60335249368e12 * POW(age,8) ;
      fb_energy += 2.17720380795e12 * POW(age,9) ;
      fb_energy -= 1.08835588273e12 * POW(age,10);
      fb_energy += 246366864343.0   * POW(age,11);
   }
   else if ( age >= 0.518 && age <= 2.0 )
   {
      fb_energy += 300.659606389             ;
      fb_energy -= 2175.28137376 * age       ;
      fb_energy += 7038.17965731 * POW(age,2);
      fb_energy -= 12640.7809456 * POW(age,3);
      fb_energy += 13818.3936865 * POW(age,4);
      fb_energy -= 9434.54106014 * POW(age,5);
      fb_energy += 3935.34399667 * POW(age,6);
      fb_energy -= 918.140140181 * POW(age,7);
      fb_energy += 91.8268783049 * POW(age,8);
   }
   else if ( age >= 2.0 && age <= 3.23 )
   {
      fb_energy += 1955.41904193             ;
      fb_energy -= 4288.5933438  * age       ;
      fb_energy += 3935.44109106 * POW(age,2);
      fb_energy -= 1921.4747372  * POW(age,3);
      fb_energy += 526.42694795  * POW(age,4);
      fb_energy -= 76.729393462  * POW(age,5);
      fb_energy += 4.64818353202 * POW(age,6);
   }
   else if ( age >= 3.23 )
   {
      fb_energy = fb_energy + 12.6540102838;
   } // if ( age < 0.29 )
   return fb_energy;

} // FUNCTION : FB_WindEjectEnergy



//-------------------------------------------------------------------------------------------------------
// Function    : FB_SNeEjectMass 
// Description :  
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
real FB_SNeEjectMass( const real age )
{
   real fb_mass = 0.0;
   if ( age < 0.513 )
   {
      fb_mass += 3.40965833751              ;
      fb_mass -= 16.0346449798 * age        ;
      fb_mass += 31.5091825735 * POW(age, 2);
      fb_mass -= 21.3218283568 * POW(age, 3);
   }
   else if (0.513 <= age && age <= 0.918)
   {
      fb_mass -= 314538.854117               ;
      fb_mass += 4453582.08399 * age         ;
      fb_mass -= 28218211.3741 * POW(age, 2) ;
      fb_mass += 105370186.068 * POW(age, 3) ;
      fb_mass -= 256824281.305 * POW(age, 4) ;
      fb_mass += 426986197.681 * POW(age, 5) ;
      fb_mass -= 490461521.485 * POW(age, 6) ;
      fb_mass += 384394390.035 * POW(age, 7) ;
      fb_mass -= 196752045.251 * POW(age, 8) ;
      fb_mass += 59399337.5861 * POW(age, 9) ;
      fb_mass -= 8033095.66643 * POW(age, 10);
   }
   else if (0.918 <= age && age <= 3.23)	
   {
      fb_mass += 1.74261906723                 ;
      fb_mass -= 0.92589554122    * age        ;
      fb_mass += 0.551250718292   * POW(age, 2);
      fb_mass -= 0.220085806978   * POW(age, 3);
      fb_mass += 0.0510662546479  * POW(age, 4);
      fb_mass -= 0.00504400687495 * POW(age, 5);
   }
   else if (age >= 3.23)			
   {
      fb_mass += 2.67991943387                ;
      fb_mass -= 0.461075452846  * age        ;
      fb_mass -= 0.0326899620754 * POW(age, 2);
   } // if ( age < 0.513 )
   return fb_mass;

} // FUNCTION : FB_SNeEjectMass



//-------------------------------------------------------------------------------------------------------
// Function    : FB_UVFeedback 
// Description :  Add internal energy to account for photoionization heating.
//
// Note        : Energy is added in proportion to the ratio of the stromgren radius to the cell volume 
//               if the stromgren sphere is smaller than a cell.
//
// Parameter   : s49_tot : the total ionizing luminosity for this particle
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
void FB_UVFeedback( float (*Fluid)[PS2][PS2][PS2], const real particleIonizingLuminosity, 
                    const int pararticleWillExplode, const int idx[], const double dh )
{
   #ifdef DUAL_ENERGY
   const real cell_volume = POW(dh*UNIT_L, 3);
   const real mu          = 0.6;
   const real mass_h      = 1.67262171e-24;
   const real alpha       = 2.60e-13;               // Case B recombination, assuming T = 10^4 K.
   
   const real flu_dens    = Fluid[DENS][ idx[0] ][ idx[1] ][ idx[2] ];      // old density of this cell
   const real flu_energy  = Fluid[ENPY][ idx[0] ][ idx[1] ][ idx[2] ];      // old energy of this cell
   
   const real num = flu_dens * UNIT_D / mu / mass_h;
   const real stromgren_radius = POW( ( (3.0*particleIonizingLuminosity*1e49) / (4.0*M_PI*alpha*SQR(num)) )
                                    , (1.0/3.0) );
   const real stromgren_volume = (4.0/3.0) * M_PI * CUBE(stromgren_radius);
      
   real ionized = Const_kB*1e4 / mu / mass_h / SQR(UNIT_L) * SQR(UNIT_T);
      
   if ( stromgren_volume < cell_volume )  ionized *= stromgren_volume/cell_volume;

   if ( pararticleWillExplode == 1 && flu_energy < ionized )
   {
      // TODO: dual energy update
      const real theDiff = ionized - flu_energy; // thermal energy difference caused by UV ionization
      Fluid[ENGY][idx[0]][idx[1]][idx[2]] += theDiff; // Update the energy
      Fluid[ENPY][idx[0]][idx[1]][idx[2]] += theDiff;
   }
      #endif // #ifdef DUAL_ENERGY
} // FUNCTION : FB_UVFeedback



//-------------------------------------------------------------------------------------------------------
// Function    : FB_WindFeedback
// Description :  
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
void FB_WindFeedback( const real particleAge, real &windEjectMass, real &windEjectEnergy, 
                      const real particleSoonestExplosion, real &coldWind, const real minimumStarMass,
                      const double dt )
{
   double wind_mass   = FB_WindEjectMass( particleAge );
   double wind_energy = FB_WindEjectEnergy( particleAge );
      
   // If the most massive star in this particle has a delay time greater than 22 Myr
   // (i.e. it is low mass) turn on 10 km/s wind outflows.
   if ( particleSoonestExplosion < 0 )                    coldWind = 1.0;
   if ( particleSoonestExplosion*UNIT_T/3.16e14 > 2.20 )  coldWind = 1.0;
   if ( coldWind == 1 )                  windEjectEnergy = 12.6540102838;

   // Convert windEjectEnergy from cm^2/s^2 to code units
   windEjectEnergy = POW( 10.0, ( windEjectEnergy-2.0*log10(UNIT_V) ) );

   // Find mass ejection rate in solar masses per year
   windEjectMass = POW( 10.0, windEjectMass );

   // This ejection rate is correct for the 10^6 solar mass 'cluster'
   // used to compute these rates with starburst 99.
   // Reduce to account for the size of the star particle
   windEjectMass *= minimumStarMass / 1e6;

   // Convert mass ejection rate from solar masses per year to solar masses per second
   windEjectMass /= Const_yr;

   // Convert windEjectMass into a mass in solar masses rather than 
   // a mass ejection rate by multiplying through by dt in seconds
   windEjectMass *= dt * UNIT_T;

} // FUNCTION : FB_WindFeedback
*/


//-------------------------------------------------------------------------------------------------------
// Function    : FB_SNMetalFeedback
// Description :  
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
void FB_SNMetalFeedback( real (*Fluid)[PS2][PS2][PS2], const real windEjectDensity, 
                         const real SNeEjectDensity, const real particleAge, const int explosionFlag, 
                         const int imetal, const real metalf, const int idx[], 
                         const int imetalSNII, const real metalfSNIIn)
{
   if ( imetal != 1 ) return;    // metal feedback option not open
   
   const real fluidDensity = Fluid[ENGY][ idx[0] ][ idx[1] ][ idx[2] ];
   real fluidMetal; // need to be initialized
   real fluidMetalSNII; // need to be initialized
   
   real ZSN = 0.0;
   if ( explosionFlag > 0 )
   {
      if ( particleAge <= 0.365 )
      {
         ZSN -= 2.12788803787                        ;
         ZSN -= 3.08562080661 * particleAge          ;
         ZSN += 83.9782331967 * POW( particleAge, 2 );
         ZSN -= 158.867576409 * POW( particleAge, 3 );
      }
      else if ( 0.365 <= particleAge && particleAge <= 0.442)
      {
         ZSN += 0.584246828515                       ;
         ZSN -= 10.2313474777 * particleAge          ;
         ZSN += 42.3869874214 * POW( particleAge, 2 );
         ZSN -= 46.7483110349 * POW( particleAge, 3 );
      }
      else if ( 0.442 <= particleAge && particleAge <= 0.517 )
      {
         ZSN -= 6913398.58591                        ;
         ZSN += 99253357.5895 * particleAge          ;
         ZSN -= 610064398.162 * POW( particleAge, 2 );
         ZSN += 2081033583.17 * POW( particleAge, 3 );
         ZSN -= 4254712179.17 * POW( particleAge, 4 );
         ZSN += 5213617301.35 * POW( particleAge, 5 );
         ZSN -= 3545287887.67 * POW( particleAge, 6 );
         ZSN += 1032027574.24 * POW( particleAge, 7 );
      }
      else if ( 0.517 <= particleAge && particleAge <= 0.718)
      {
         ZSN -= 0.856562917033                        ;
         ZSN += 2.72919708256  * particleAge          ;
         ZSN -= 1.22818002617  * POW( particleAge, 2 );
         ZSN -= 0.536735943166 * POW( particleAge, 3 );
      }
      else if ( 0.718 <= particleAge && particleAge <= 1.7)
      {
         ZSN += 6.69150908774                        ;
         ZSN -= 51.5159032945 * particleAge          ;
         ZSN += 172.980768687 * POW( particleAge, 2 );
         ZSN -= 307.64633271  * POW( particleAge, 3 );
         ZSN += 311.333072359 * POW( particleAge, 4 );
         ZSN -= 180.28789728  * POW( particleAge, 5 );
         ZSN += 55.6987751625 * POW( particleAge, 6 );
         ZSN -= 7.12566268516 * POW( particleAge, 7 );
      }
      else if ( 1.7 <= particleAge && particleAge <= 3.2 )
      {
         ZSN += 0.173226556725                          ;
         ZSN -= 0.142758066129   * particleAge          ;
         ZSN += 0.0568361245745  * POW( particleAge, 2 );
         ZSN -= 0.00740131298973 * POW( particleAge, 3 );
      }
      else
      {
         ZSN -= 0.459335515363                         ;
         ZSN += 0.51236004354   * particleAge          ;
         ZSN -= 0.194211146544  * POW( particleAge, 2 );
         ZSN += 0.0264297153731 * POW( particleAge, 3 );
      } // if ( particleAge <= 0.365 )
   } // if ( explosionFlag[n] > 0 )

   fluidMetal *= fluidDensity;
   fluidMetal += metalf * windEjectDensity;
   fluidMetal += ZSN * SNeEjectDensity;
   fluidMetal /= ( windEjectDensity + SNeEjectDensity + fluidDensity );
   // Fluid[METAL][ idx[0] ][ idx[1] ][ idx[2] ]= fluidMetal;
      
   if (imetalSNII == 1)
   {
      fluidMetalSNII *= fluidDensity;
      fluidMetalSNII += ( metalfSNIIn * windEjectDensity );
      fluidMetalSNII += ZSN * SNeEjectDensity;
      fluidMetalSNII /= ( windEjectDensity + SNeEjectDensity + fluidDensity );
      // Fluid[METALSNII][ idx[0] ][ idx[1] ][ idx[2] ] = fluidMetalSNII;
   } // if (imetalSNII == 1)

} //  FUNCTION : FB_SNMetalFeedback



//-------------------------------------------------------------------------------------------------------
// Function    : FB_SNeMomentumFeedback
// Description :  
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
void FB_SNeMomentumFeedback( real (*Fluid)[PS2][PS2][PS2], const int particleExplosionFlag, 
                             const real windEjectDensity, 
                             const real SNeEjectDensity, const real SNeEjectEnergy,
                             const int imethod, const int idual, 
                             const int distrad, const int diststep, const int distcells,
                             const double EdgeL[], const int idx[], const double dh)
{
   // Check the condition
   if ( SNeEjectEnergy == 0 ) return; // no supernova explosion event
   if ( distrad == 0 )   return;      // not open momentum feedback option
   
   
   // Find momentum to deposit per SNe in code units.
   // TODO: Change to GAMER format
   const double Msun = 1.9885e33;
   const real SN_1CellMom = particleExplosionFlag * 3e5 * Msun / (UNIT_D*POW(UNIT_L, 3)) * (Const_km/Const_s) / UNIT_V / 26;

   // Approximate that the explosion happens in the center of the feedback zone
   // TODO: This should be wrong, I don't know how to get the position in GAMER
   // TODO: Define xstart, ystart, zstart
   const real explode_x = EdgeL[0] + ( (real) idx[0] - 0.5 ) * dh;
   const real explode_y = EdgeL[1] + ( (real) idx[1] - 0.5 ) * dh;
   const real explode_z = EdgeL[2] + ( (real) idx[2] - 0.5 ) * dh;
   
   const real centerDens = Fluid[DENS][ idx[0] ][ idx[1] ][ idx[2] ];

   // do momentum feedback over a distributed region distribute equally among zones participating in feedback
   // Limit velocity after feedback to 1000 km/s to avoid accelerating material in cells with very little mass
   real kesum = 0.0;
   real cellMom = 0.0;
   const real maxVel = 1e3*( Const_km/Const_s )/UNIT_V;  // TODO: Change to GAMER format
   
   const real cellVol = CUBE( dh );
   
   const int size = 3;
   const int size2 = SQR( size );
   int idx_temp = 0;
   for (int idx_temp = 0; CUBE( size ); idx_temp++)
   {     
      // index for the feedback zone
      // minus one is for center shift
      const int di = idx_temp % size2        - 1; // kc
      const int dj = idx_temp % size2 / size - 1;
      const int dk = idx_temp / size         - 1;
      
      // index for the fluid array
      const int flu_i = idx[0] + di; // k
      const int flu_j = idx[1] + dj; 
      const int flu_k = idx[2] + dk; 
   
      const int cellstep = abs(di) + abs(dj) + abs(dk);
      if ( cellstep <= diststep )
      {
         if ( cellstep == 0 )  continue; // Don't do momentum feedback in the cell the SN goes off in

         // center position of the cell
         // TODO: This should be wrong, I don't know how to get the position in GAMER
         const real cell_x = EdgeL[0] + ( (real) di - 0.5 ) * dh;
         const real cell_y = EdgeL[1] + ( (real) dj - 0.5 ) * dh;
         const real cell_z = EdgeL[2] + ( (real) dk - 0.5 ) * dh;

         const real radius = SQRT( SQR( cell_x - explode_x ) + 
                                   SQR( cell_y - explode_y ) + 
                                   SQR( cell_z - explode_z )  );

         const real rx = (cell_x - explode_x) / radius;
         const real ry = (cell_y - explode_y) / radius;
         const real rz = (cell_z - explode_z) / radius;

         const real cellMass = Fluid[DENS][ flu_i ][ flu_j ][ flu_k ] * cellVol;

         // Limit momentum feedback such that each cell can do no more than 
         // its fair share of the thermal energy from the SNe
         if ( SQR(SN_1CellMom) / (2.0*cellMass) > SNeEjectEnergy / (distcells-1) )
         {
            cellMom = 2.0 * cellMass * SNeEjectEnergy / (distcells-1);
            cellMom = SQRT(cellMom);
         }
         else
         {
            cellMom = SN_1CellMom;
         }

         // Limit maximum velocity boost to 1000 km/s
         if ( (SN_1CellMom/cellMass) > maxVel )
         {
            cellMom = cellMass * maxVel;
         } // if ( (SN_1CellMom/cellMass) > maxVel )

         // Add momentum feedback
         Fluid[MOMX][ flu_i ][ flu_j ][ flu_k ] += rx * cellMom / cellMass;
         Fluid[MOMY][ flu_i ][ flu_j ][ flu_k ] += ry * cellMom / cellMass;
         Fluid[MOMZ][ flu_i ][ flu_j ][ flu_k ] += rz * cellMom / cellMass;

         // TODO : check the dual energy and switch to GAMER variable
         #ifdef DUAL_ENERGY
         Fluid[ENGY][ flu_i ][ flu_j ][ flu_k ] = 
	 Fluid[ENPY][ flu_i ][ flu_j ][ flu_k ] + 
	 0.5*( POW(Fluid[MOMX][ flu_i ][ flu_j ][ flu_k ], 2) + 
	 POW(Fluid[MOMY][ flu_i ][ flu_j ][ flu_k ], 2) + 
	 POW(Fluid[MOMZ][ flu_i ][ flu_j ][ flu_k ], 2) ) /
	 Fluid[DENS][ flu_i ][ flu_j ][ flu_k ];
         #endif // #ifdef DUAL_ENERGY

         // Record kinetic energy deposited in this cell
         kesum += 0.5 / cellMass * SQR( cellMom );
      } // if ( cellstep <= diststep )

   } // for (int idx_temp = 0; CUBE( size ); idx_temp++)

   
   if ( kesum > SNeEjectEnergy )
   {  
      // TODO : change the output to GAMER function
      printf( "ssn1: Kinetic energy from momentum feedback exceeds thermal energy\n" );
      kesum = SNeEjectEnergy;
   } // if ( kesum > SNeEjectEnergy )

//   if ( kesum / (centerDens * cellVol) >= Fluid[ENPY][ idx[0] ][ idx[1] ][ idx[2] ] )
//   {
      // TODO : chenge the output to GAMER function
//      Aux_Error( ERROR_INFO, "ssn2: Kinetic energy from momentum feedback exceeds thermal energy\n" );
//   } // if ( kesum / (Fluid[DENS][ idx[0] ][ idx[1] ][ idx[2] ] * cellVol) >= ge(i,j,k) )

   // Assume kinetic energy from momentum feedback is extracted from thermal energy in the host zone
   // TODO : dual energy update
   Fluid[ENGY][ idx[0] ][ idx[1] ][ idx[2] ] -= kesum / (centerDens * cellVol);
   
} // FUNCTION : FB_SNeMomentumFeedback



#endif // #ifdef FEEDBACK
