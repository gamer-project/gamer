#include "Copyright.h"
#include "GAMER.h"

#if (  defined PARTICLE  &&  defined STAR_FORMATION  &&  ( MODEL==HYDRO || MODEL==MHD )  )


#ifdef GRAVITY
#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_AGORA
// Description :  Create new star particles stochastically using the presription suggested by the AGORA project
//
// Note        :  1. Ref: (1) Nathan Goldbaum, et al., 2015, ApJ, 814, 131 (arXiv: 1510.08458), sec. 2.4
//                        (2) Ji-hoon Kim, et al., 2016, ApJ, 833, 202 (arXiv: 1610.03066), sec. 3.2
//                2. One must turn on STORE_POT_GHOST When adopting STORE_PAR_ACC
//                   --> It is because, currently, this function always uses the pot_ext[] array of each patch
//                       to calculate the gravitationally acceleration of the new star particles
//                3. One must invoke Buf_GetBufferData( ..., _TOTAL, ... ) after calling this function
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Current physical time (after advancing solution by dt)
//                dt           : Time interval to advance solution
//                               --> Currently this function does not distinguish dt and the physical time interval (dTime)
//                               --> Does NOT support COMOVING yet
//                drand_buf    : Buffer for the reentrant and thread-safe random number generator drand48_r()
//                GasDensThres : Minimum gas density for creating star particles                (--> "SF_CREATE_STAR_MIN_GAS_DENS"  )
//                Efficiency   : Gas-to-star mass efficiency                                    (--> "SF_CREATE_STAR_MASS_EFF"      )
//                MinStarMass  : Minimum star particle mass for the stochastical star formation (--> "SF_CREATE_STAR_MIN_STAR_MASS" )
//                MaxStarMFrac : Maximum gas mass fraction allowed to convert to stars          (--> "SF_CREATE_STAR_MAX_STAR_MFRAC")
//                UseMetal     : Store the metal mass fraction in star particles
//
// Return      :  1. Particle repository will be updated
//                2. fluid[] array of gas will be updated
//-------------------------------------------------------------------------------------------------------
void SF_CreateStar_AGORA( const int lv, const real TimeNew, const real dt, struct drand48_data *drand_buf,
                          const real GasDensThres, const real Efficiency, const real MinStarMass, const real MaxStarMFrac,
                          const bool UseMetal )
{

// check
#  if ( defined STORE_PAR_ACC  &&  !defined STORE_POT_GHOST )
#     error : STAR_FORMATION + STORE_PAR_ACC must work with STORE_POT_GHOST !!
#  endif

#  ifndef GRAVITY
#     error : must turn on GRAVITY for SF_CreateStar_AGORA() !!
#  endif

#  ifdef COMOVING
#     error : SF_CreateStar_AGORA() does not support COMOVING yet !!
#  endif

   if ( UseMetal )
   {
      if ( PAR_NPASSIVE != 2 )
         Aux_Error( ERROR_INFO, "currently the metal field is hard coded for UseMetal and thus PAR_NPASSIVE (%d) must be 2 !!\n",
                    PAR_NPASSIVE );
   }

   else
   {
      if ( PAR_NPASSIVE != 1 )
         Aux_Error( ERROR_INFO, "currently the passive particle attributes must be hard coded when adding the new particles !!\n" );
   }


   const double dh             = amr->dh[lv];
   const real   dv             = CUBE( dh );
   const int    FluSg          = amr->FluSg[lv];
   const int    PotSg          = amr->PotSg[lv];
   const real   Coeff_FreeFall = SQRT( (32.0*NEWTON_G)/(3.0*M_PI) );
   const real  _MinStarMass    = (real)1.0 / MinStarMass;
   const real   Eff_times_dt   = Efficiency*dt;

   double x0, y0, z0, x, y, z, Random;
   real   GasDens, _GasDens, GasMass, _Time_FreeFall, StarMFrac, GasMFracLeft;
   real   NewParVar[PAR_NVAR], NewParPassive[PAR_NPASSIVE];
   real   (*fluid)[PS1][PS1][PS1] = NULL;
   real   (*pot)[PS1][PS1]        = NULL;
   long   ParID;


// loop over all real patches
//###: OPENMP
//#  pragma omp parallel for
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;


      fluid = amr->patch[FluSg][lv][PID]->fluid;
      pot   = amr->patch[PotSg][lv][PID]->pot;
      x0    = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0    = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0    = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {

//       1. check the star formation criteria
//       ===========================================================================================================
         GasDens = fluid[DENS][k][j][i];
         GasMass = GasDens*dv;

//       1-1. create star particles only if the gas density exceeds the given threshold
         if ( GasDens < GasDensThres )    continue;


//       1-2. estimate the gas free-fall time
//       --> consider only the gas density under the assumption that the dark matter doesn't collapse
         _Time_FreeFall = Coeff_FreeFall * SQRT( GasDens );


//       1-3. estimate the gas mass fraction to convert to stars
         StarMFrac = Eff_times_dt*_Time_FreeFall;
         StarMass  = GasMass*StarMFrac;


//       1-4. stochastic star formation
//       --> if the star particle mass (StarMass) is below the minimum mass (MinStarMass), we create a
//           new star particle with a mass of MinStarMass and a probability of StarMass/MinStarMass
//       --> Eq. [5] in Goldbaum et al. (2015)
         if ( StarMass < MinStarMass )
         {
//###: OPENMP
//          drand48_r( drand_buf[thread_id], &Random )
            drand48_r( drand_buf[0], &Random );

            if ( (real)Random < StarMass*_MinStarMass )  StarMFrac = MinStarMass / GasMass;
            else                                         continue;
         }


//       1-5. check the maximum gas mass fraction allowed to convert to stars
         StarMFrac = MIN( StarMFrac, MaxStarMFrac );
         StarMass  = GasMass*StarMFrac;



//       2. create the star particles
//       ===========================================================================================================
//###: OPENMP
//       2-1. calculate the new particle attributes
//       2-1-1. intrinsic attributes
         _GasDens = (real)1.0 / GasDens;
         x        = x0 + i*dh;
         y        = y0 + j*dh;
         z        = z0 + k*dh;

         NewParVar[PAR_MASS] = StarMass;
         NewParVar[PAR_POSX] = x;
         NewParVar[PAR_POSY] = y;
         NewParVar[PAR_POSZ] = z;
         NewParVar[PAR_VELX] = fluid[MOMX][k][j][i]*_GasDens;
         NewParVar[PAR_VELY] = fluid[MOMY][k][j][i]*_GasDens;
         NewParVar[PAR_VELZ] = fluid[MOMZ][k][j][i]*_GasDens;
         NewParVar[PAR_TIME] = TimeNew;

#        ifdef STORE_PAR_ACC
         NewParVar[PAR_ACCX] = NULL_REAL;
         NewParVar[PAR_ACCY] = NULL_REAL;
         NewParVar[PAR_ACCZ] = NULL_REAL;
#        endif

//       2-1-2. passive attributes
//###: HARD-CODED FIELDS
//       note that we store the metal mass **fraction** instead of density in particles
//       currently the metal field is hard coded ... ugh!
         if ( UseMetal )
         NewParPassive[PAR_METAL_FRAC   ] = fluid[METAL] * _GasDens;

         NewParPassive[PAR_CREATION_TIME] = TimeNew;


//       2-2. add particles to the particle repository
         ParID = amr->Par->AddOneParticle( NewParVar, NewParPassive );


//       2-3. add particles to the patch
#        ifdef DEBUG_PARTICLE
//       do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//       may change after calling amr->Par->AddOneParticle()
         const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         char Comment[100];
         sprintf( Comment, "%s", __FUNCTION__ );

         amr->patch[0][lv][PID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[lv],
                                              ParPos, amr->Par->NPar_AcPlusInac, Comment );
#        else
         amr->patch[0][lv][PID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[lv] );
#        endif



//       3. remove the gas that has been converted to stars
//       ===========================================================================================================
         GasMFracLeft = (real)1.0 - StarMFrac;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v][k][j][i] *= GasMFracLeft;
      } // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : SF_CreateStar_AGORA



#endif // #if (  defined PARTICLE  &&  defined STAR_FORMATION  &&  ( MODEL==HYDRO || MODEL==MHD )  )
