#include "GAMER.h"
#include <fstream>

#ifdef FEEDBACK




// function pointers to be set by FB_Init_User_Template()
extern void (*FB_SNe_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                            real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[],
                            const int TID, RandomNumber_t *RNG );
extern void (*FB_End_SNe_Ptr)();


// internal funcitons
float s49Lookup(float m);

float UV_IonizingLuminosity( const double m );


int FB_GetNumberOfSNe( const real MinimumMass, RandomNumber_t *RNG, const int TID /*, const RandomSeed*/ );

real FB_WindEjectMass( const double age );

real FB_WindEjectEnergy( const double age );

real FB_SNeEjectMass( const double age );

void FB_UVFeedback( double (*Fluid)[PS2][PS2][PS2], const double particleIonizingLuminosity,
                    const int pararticleWillExplode, const int idx[], const double dh );

void FB_WindFeedback( const double particleAge, real &windEjectMass, real &windEjectEnergy,
                      const real particleSoonestExplosion, real &coldWind, const real minimumStarMass,
                      const double dt );


//-------------------------------------------------------------------------------------------------------
// Function    :  FB_SNe
// Description :  Supernova explosion feedback
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
void FB_SNe( const int lv, const double TimeNew, const double TimeOld, const double dt,
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
#  endif // #ifdef GAMER_DEBUG
  const int par_idx = *ParSortID;
  float m = ParAtt[PAR_MASS][par_idx];
  s49Lookup(m);

// Needed for determine SNe number in a particle
   int *explosionFlag      = new int[NPar];
   int *explosionFlagIa    = new int[NPar];
   int *willExplode        = new int[NPar];
   real *soonestExplosion  = new double[NPar];
   real *ionizeLuminosity  = new double[NPar];
   int rand_int            = 0;

// recording
   FILE *fout;
   FILE *Fout;
   FILE *FOut;
   FILE *FOUt;
   FILE *FOUT;

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
      const double Msun = Const_Msun;
      const double MinParMass = 3.2e3;     //in unit of M_sun

      // Particle variables
      const int par_idx        = ParSortID[n];
      const double par_pos[3]  = { ParAtt[PAR_POSX][par_idx],
                                   ParAtt[PAR_POSY][par_idx],
                                   ParAtt[PAR_POSZ][par_idx] };
      const double par_vel[3]  = { ParAtt[PAR_VELX][par_idx],
                                   ParAtt[PAR_VELY][par_idx],
                                   ParAtt[PAR_VELZ][par_idx] };
      const double par_creTime = ParAtt[Idx_ParCreTime][par_idx];          // creation time of particle
      const double par_age     = ( TimeOld-par_creTime );   // age of a particle
      double par_mass          = ParAtt[PAR_MASS][par_idx];


      // Fluid variables
      int idx[3];    // define the coordinate of a particle_consisting cell
      for ( int d = 0; d < 3; d++ ) idx[d] = (int)FLOOR( ( par_pos[d] - EdgeL[d] )*_dh );
      // check particle is in grid or not
      if ( idx[0] < 0 || idx[0] >= PS2 || idx[1] < 0 || idx[1] >= PS2 || idx[2] < 0 || idx[2] >= PS2 ) continue;

      real flu_dens      = Fluid[DENS][idx[0]][idx[1]][idx[2]];      // density of this cell
      real flu_energy    = Fluid[ENGY][idx[0]][idx[1]][idx[2]];      // energy of this cell

      // Feedback variables
      real fb_energy     = 0.0;                    // total energy feedback of particle
      real fb_energy_Msun = 0.0;

      // 1: UV feedback
      real ionized;

      // 2: Wind feedback
      real wind_mass     = 0.0;                    // total mass feedback by wind
      real wind_energy   = 0.0;                    // total energy feedback by wind
      real wind_dens     = 0.0;
      real coldWind      = 0.0;                    // wind of old particle

      // 3: SNe feedback
      real SNe_mass      = 0.0;
      real SNe_dens      = 0.0;
      real sn_energy     = 0.0;                    // energy emitted by supernova
      real SNe_dens_Msun = 0.0;

      // 4: Metal feedback
      real ZSN           = 0.0;

      // 5: Momentum feedback
      real p_sn          = 0.0;                    // momentum to deposit per SNe in code units


      // ============================================================
      // Check the particle and the patch 
      // ============================================================
      // 1: The particle is too old, so there is no feedback

//      printf( "par_age = %e\n", par_age );
//      printf( "par_creTime = %e\n", par_creTime );
//      printf( "par_mass = %e\n", par_mass );

      if ( par_age >= 3.7e7 * Const_yr / UNIT_T )       continue;

      // 2. check bounds - if star particle is outside of this grid
      // then something is very wrong and we should crash


      // 3: The particle creation time 
      if ( par_creTime > 0.0 )        continue;

      // 4: The particle mass
      if ( par_mass <= 0.0 )         continue;




      // ============================================================
      // Calculate the nummber of the feedback events in particle
      // ============================================================
      soonestExplosion[n] = -1.0;
      explosionFlag[n] = -1;
      willExplode[n]   =  0;

      // figure out how many explosions this particle undergoes.
      const int number_of_sn = FB_GetNumberOfSNe( MinParMass, RNG, TID );
//      printf("number of supernova = %d\n", number_of_sn);

      // If there are no SN in this particle, go on to the next one.
      if ( number_of_sn == 0 )
      {
         explosionFlag[n] = 0;
         continue; // Go to next particle
      } // if ( number_of_sn == 0 )
//      printf("explosion flag = %d\n", explosionFlag[n]);



      // ============================================================
      // Calculate physics feedback parameter
      // ============================================================
      // If there are any explosions, we need to know the ionizing luminosity of each progenitor
      // and whether or not that progenitor will explode soon or in this timestep.

      for ( int kk = 0; kk < number_of_sn; ++kk )     // Loop over expected SN.
      {
         // Draw the delay time for this SN event.
         const real x = RNG->GetValue( TID, 0, 32767.0 ) / 32768.0;
         real delayTime =
               p_delay[0] +
               p_delay[1]*x +
               p_delay[2]*x*x +
               p_delay[3]*x*x*x +
               p_delay[4]*x*x*x*x +
               p_delay[5]*x*x*x*x*x; // in years.
//	 printf("第%d個爆炸在%f時\n", kk, delayTime);

         // Convert the delay time to a progenitor mass.
         const real td7 = delayTime*1e7;
//	 printf("delay time = %f\n", td7);
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
//	 printf("progenitor mass scale = %e\n", progenitor_mass);
	 progenitor_mass = POW( 10, progenitor_mass ); // in solar masses

         delayTime *= Const_yr/UNIT_T; // convert years -> seconds -> code units.
         const real relativeTime = TimeOld - ParAtt[Idx_ParCreTime][par_idx]; // in simulation units.

//	 printf("explosion flag before explosion = %d\n", explosionFlag[n]);

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
//	    printf("progenitor mass = %e\n", progenitor_mass);
//	    printf("luminosity by unexploded SNe = %e\n", UV_IonizingLuminosity( progenitor_mass ));
//	    printf("total luminosity = %e\n", ionizeLuminosity[n]);

            willExplode[n] = 1;

            if ( soonestExplosion[n] < 0 || delayTime < soonestExplosion[n] )
            {
              soonestExplosion[n] = delayTime;
           } // if ( soonestExplosion[n]<0 || delayTime<soonestExplosion[n])

        } // if( relativeTime < delayTime  )

//         if ( explosionFlag[n] == -1 )  explosionFlag[n] = 0;

      } // for ( int kk=0; kk < number_of_sn; ++kk )
      if ( explosionFlag[n] == -1 )  explosionFlag[n] = 0;
//      printf("fluid density = %f\n", flu_dens);




      // ============================================================
      // Feedback part
      // ============================================================
      // 0: Presetting
//    double (*fluid_PG)[PS2][PS2][PS2] = new double [NCOMP_TOTAL][PS2][PS2][PS2];
//    memcpy( fluid_PG, Fluid, sizeof(Fluid) );

      // 1: UV feedback
//      FB_UVFeedback( fluid_PG, ionizeLuminosity[n], willExplode[n], idx, dh);

      // 2: Wind feedback
      double feedback_t = par_age * UNIT_T / Const_Myr / 10.0;
//      FB_WindFeedback( feedback_t, wind_mass, wind_energy, soonestExplosion[n], coldWind, MinParMass, dt );


      // 3: SN energy feedback
      // evolution of SN ejection by time
      if ( explosionFlag[n] > 0 )       SNe_mass = FB_SNeEjectMass( feedback_t );

      SNe_mass  = POW(10.0, SNe_mass) * explosionFlag[n];

      SNe_dens  = SNe_mass * Msun / CUBE( dh*UNIT_L ) / UNIT_D * Msun / UNIT_M;
      wind_dens = wind_mass * Msun / CUBE( dh*UNIT_L ) / UNIT_D;
//      printf("SNE mass density = %f\n", SNe_dens);
//      printf("cell width = %e\n", dh);
//      concert SNe_dens into solar mass to compare with enzo
      SNe_dens_Msun = SNe_dens * UNIT_M / Msun;

      FOUt=fopen("time.txt","a");
      if(FOUt==NULL) {
	printf("Fail To Open File time.txt!!");
                     }
      fprintf(FOUt, "%f\n", TimeOld);
      fclose(FOUt);

      FOUT=fopen("timestep.txt","a");
      if(FOUT==NULL) {
        printf("Fail To Open File timestep.txt!!");
                     }
      fprintf(FOUT, "%f\n", dt);
      fclose(FOUT);

      if(SNe_dens_Msun != 0){
      	fout=fopen("mass.txt","a");
      	if(fout==NULL) {
        	printf("Fail To Open File mass.txt!!");
                       }
      	fprintf(fout, "%f\n", SNe_dens_Msun);
      	fclose(fout);

      	FOut=fopen("timestamp.txt","a");
      	if(FOut==NULL) {
        	printf("Fail To Open File timestamp.txt!!");
                       }
      	fprintf(FOut, "%f\n", TimeOld);
      	fclose(FOut);
			   }//if(SNe_dens_Msun != 0)

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
//      printf("supernova energy = %e\n", sn_energy);
//      printf("wind energy = %e\n", wind_energy);
      fb_energy = sn_energy + wind_energy * wind_dens * CUBE( dh );
//    convert energy into solar mass result to compare with enzo
      fb_energy_Msun = fb_energy / ( (flu_dens + wind_dens + SNe_dens_Msun) * CUBE( dh ) );
      fb_energy = fb_energy / ( (flu_dens + wind_dens + SNe_dens) * CUBE( dh ) );

      if(fb_energy_Msun != 0)
	{
      	Fout=fopen("energy.txt","a");
      	if(Fout==NULL) {
        	printf("Fail To Open File energy.txt!!");
                       }
      	fprintf(Fout, "%f\n", fb_energy_Msun);
      	fclose(Fout);
	}//if(fb_energy_Msun != 0)

      // add energy into fluid
//      printf("feedback energy = %e\n", fb_energy);
      if((idx[0] < 0) || (idx[1] < 0) || (idx[2] < 0))	{
      	printf("x left edge = %e\n", EdgeL[0]);
      	printf("y left edge = %e\n", EdgeL[1]);
      	printf("z left edge = %e\n", EdgeL[2]);
      	printf("x coordinate = %2.15f\n", par_pos[0]);
      	printf("y coordinate = %2.15f\n", par_pos[1]);
      	printf("z coordinate = %2.15f\n", par_pos[2]);
      	printf("x index = %d\n", idx[0]);
      	printf("y index = %d\n", idx[1]);
      	printf("z index = %d\n", idx[2]);
      	printf("x velocity = %e\n", par_vel[0]);
      	printf("y velocity = %e\n", par_vel[1]);
      	printf("z velocity = %e\n", par_vel[2]);	}
//      printf("fluid energy = %e\n", Fluid[ENGY][idx[0]][idx[1]][idx[2]]);

      Fluid[ENGY][idx[0]][idx[1]][idx[2]] += fb_energy;

//    memcpy( Fluid, fluid_PG, sizeof(*fluid_PG) );



      // 4: SN Metal feedback (in metal fraction)
      // no metal feedback now
      // TODO
      // FB_SNMetalFeedback( Fluid, wind_energy, sn_energy, par_age, 
      //                     imetal, metalf, imetalSNII, metalfSNIIn );
      //
      //
      // Mass and momentum feedback from winds and SN ejecta (NOT from SN blastwaves)
      Fluid[MOMX][idx[0]][idx[1]][idx[2]] += (wind_dens + SNe_dens) * par_vel[0];
      Fluid[MOMY][idx[0]][idx[1]][idx[2]] += (wind_dens + SNe_dens) * par_vel[1];
      Fluid[MOMZ][idx[0]][idx[1]][idx[2]] += (wind_dens + SNe_dens) * par_vel[2];
      Fluid[DENS][idx[0]][idx[1]][idx[2]] +=  wind_dens + SNe_dens;

      #ifdef DUAL_ENERGY
      Fluid[ENPY][idx[0]][idx[1]][idx[2]] += fb_energy;
      #endif // #ifdef DUAL_ENERGY

      // 5: SN Momentum feedback
      // no momentum feedback now
      // TODO
      // FB_SNeMomentumFeedback( Fluid, explosionFlag[n], wind_dens, SNe_dens, sn_energy, imethod,  idual, 
      //                         distrad, diststep, distcells, dh )



   } // for (n = 0; n < NPar; n++ )

   delete[] ionizeLuminosity;
   delete[] explosionFlag;
   delete[] explosionFlagIa;
   delete[] willExplode;
   delete[] soonestExplosion;

}// FUNCTION : FB_SNe



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_SNe
// Description :  Free the resources used by the SNe feedback
//
// Note        :  1. Invoked by FB_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_SNe()
{

} // FUNCTION : FB_End_SNe



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_SNe
// Description :  Initialize the SNe feedback
//
// Note        :  1. Invoked by FB_Init()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init_SNe()
{

//   FB_SNe_Ptr     = FB_SNe;
//   FB_End_SNe_Ptr = FB_End_SNe;

} // FUNCTION : FB_Init_SNe



//-------------------------------------------------------------------------------------------------------
//// Function    :  s49Lookup
//// Description :  Returns a broken power-law approximation to the ionizing luminosity emitted 
////		    by a main sequence star of mass m, normalized to 10^49 photon/s
////		    See Parravano et al. 2003 ApJ 584 797 (doi:10.1086/345807)
////
//// Note        :  None
////
//// Parameter   :  Mass of the particle
////
//// Return      :  Luminosity (10^49 photon/s)
////-------------------------------------------------------------------------------------------------------
float s49Lookup(float m)
{
       if (m<5)
          {
          return 0;
          }
          else if (m<7)
          {
          return (2.23e-15*pow(m, 11.5));
          }
          else if (m<12)
          {
          return(3.69e-13*pow(m, 8.87));
          }
          else if (m<20)
          {
          return(4.8e-12*pow(m, 7.85));
          }
          else if (m<30)
          {
          return(3.12e-8*pow(m, 4.91));
          }
          else if (m<40)
          {
          return(2.8e-5*pow(m, 2.91));
          }
          else if (m<60)
          {
          return(3.49e-4*pow(m, 2.23));
          }
          else
          {
          return(2.39e-3*pow(m, 1.76));
          }
} // FUNCTION : s49Lookup


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
float UV_IonizingLuminosity( const double m )
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



//-------------------------------------------------------------------------------------------------------
// Function    : FB_GetNumberOfSNe
// Description : This function gives the number of Supernovae expected from a 10^6 Msun cluster,
// 		 as predicted from Starburst99.
//
// Note        :
//
// Parameter   : lambda : the expected number of SN for a stellar population with mass MinimumMass 
//		 k      : the estimated number of SN for this particular particle
//
// Return      : k
//
// Reference   : http://goo.gl/sgLPcj
//-------------------------------------------------------------------------------------------------------
int FB_GetNumberOfSNe( const real MinimumMass, RandomNumber_t *RNG, const int TID /*, const RandomSeed*/ )
{
   const float expected_sn_s99 = 10616.955572;
   const float lambda = expected_sn_s99 * 1.0e-6 * MinimumMass;
   const real L = EXP( -lambda );
   int k = 0;
   real p = 1.0;
   real u;
// printf("expectation factor = %e\n", L);

   // The Knuth algorithm for drawing from a Poisson Distribution.
   while( p > L )
   {
      ++k;
      int rand_int = RNG->GetValue( TID, 0.0, 32767.0 );
      u = rand_int / 32768.0;
      p *= u;
//    printf("random dies = %f\n", u);
//    printf("quota = %f\n", p);
//    printf("number of SNe goes to %d\n", k);
   } // while( p > L )

   // Now k-1 ~ Pois(lambda)
   k -= 1;

   return k;
} // FUNCTION : void FB_GetNumberOfSNe



//-------------------------------------------------------------------------------------------------------
// Function    : FB_WindEjectMass 
// Description : This function estimates the mass of gas ejected by stellar wind
// 		 as a function of the star particle's age.
//
// Note        :
//
// Parameter   : age     : The time this particle has been through from its creation.
// 		 fb_mass : The mass ejected by wind.
//
// Return      : fb_mass
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
real FB_WindEjectMass( const double age )
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
// Description : This function estimates the energy of gas ejected by stellar wind
//               as a function of the star particle's age. 
//
// Note        :
//
// Parameter   : age       : The time this particle has been through from its creation.
//               fb_energy : The energy ejected by wind.
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
real FB_WindEjectEnergy( const double age )
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
real FB_SNeEjectMass( const double age )
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
void FB_UVFeedback( double (*Fluid)[PS2][PS2][PS2], const real particleIonizingLuminosity,
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
      const real theDiff = ionized - flu_energy; // thermal energy difference caused by UV ionization
      Fluid[ENGY][idx[0]][idx[1]][idx[2]] += theDiff; // Update the energy
      Fluid[ENPY][idx[0]][idx[1]][idx[2]] += theDiff;
   }
      #endif // #ifdef DUAL_ENERGY
} // FUNCTION : FB_UVFeedback



//-------------------------------------------------------------------------------------------------------
// Function    : FB_WindFeedback
// Description : Convert mass and energy feedback to feedback rate of time. 
//
// Note        :
//
// Parameter   : 
//
// Return      :  
//
// Reference   : 
//-------------------------------------------------------------------------------------------------------
void FB_WindFeedback( const double particleAge, real &windEjectMass, real &windEjectEnergy,
                      const real particleSoonestExplosion, real &coldWind, const real minimumStarMass,
                      const double dt )
{
   double wind_Mass   = FB_WindEjectMass( particleAge );
   double wind_Energy = FB_WindEjectEnergy( particleAge );

   // If the most massive star in this particle has a delay time greater than 22 Myr
   // (i.e. it is low mass) turn on 10 km/s wind outflows.
   if ( particleSoonestExplosion < 0 )                    coldWind = 1.0;
   if ( particleSoonestExplosion*UNIT_T/3.16e14 > 2.20 )  coldWind = 1.0;
   if ( coldWind == 1 )                  windEjectEnergy = 12.6540102838;

   // Convert windEjectEnergy from cm^2/s^2 to code units
   windEjectEnergy = POW( 10.0, ( wind_Energy-2.0*log10(UNIT_V) ) );

   // Find mass ejection rate in solar masses per year
   windEjectMass = POW( 10.0, wind_Mass );

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




#endif // #ifdef FEEDBACK
