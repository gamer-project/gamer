#include "GAMER.h"
#ifdef SUPPORT_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#else
#error : ERROR : please turn on SUPPORT_GSL for the cluster merger test problem !!
#endif

#ifdef PARTICLE


double MassProf_DM( const double r );
static void RanVec_FixRadius( const double r, double RanVec[] );

extern int     ClusterMerger_RanSeed;
extern double  ClusterMerger_Rcut;
extern int     ClusterMerger_NBin_MassProf;
extern int     ClusterMerger_NBin_SigmaProf;

extern double *ClusterMerger_DM_MassProf;
extern double *ClusterMerger_DM_SigmaProf;
extern double *ClusterMerger_DM_MassProf_R;
extern double *ClusterMerger_DM_SigmaProf_R;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize the particle position and velocity
//
// Note        :  1. Invoked by "Init_GAMER"
//                2. Assuming all particles are active initially
//                3. To give reproducible results, currently we let the master rank to contruct the initial
//                   the particle initial condition and then broadcast to all ranks
//                   --> It can be time-consuming, and there must be a better way (e.g., use SPRNG?) ...
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass/Pos(X/Y/Z)/Vel(X/Y/Z)/Time
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };

// currently only the master rank will construct the initial condition
// --> for bit-wise reproducibility between serial and parallel runs
// --> consider using SPRNG if the performance becomes unacceptably slow
   if ( MPI_Rank == 0 )
   {
      const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
      double TotM, ParM, RanM, RanR, RanVec[3], Sigma, RcutInBox;

      Mass_AllRank   = new real [amr->Par->NPar_Active_AllRank];
      for (int d=0; d<3; d++) {
      Pos_AllRank[d] = new real [amr->Par->NPar_Active_AllRank];
      Vel_AllRank[d] = new real [amr->Par->NPar_Active_AllRank]; }


//    set the random seed
      srand( ClusterMerger_RanSeed );


//    initialize GSL random number generator
      gsl_rng *GSL_RNG = gsl_rng_alloc( gsl_rng_mt19937 );
      gsl_rng_set( GSL_RNG, ClusterMerger_RanSeed );


//    reset the dark matter cut-off radius to be within the simulation box since we do not allow particles to lie outside the
//    simulation box at the beginning
//    --> but we may want the cut-off radius for gas to be larger than the box (i.e., no cut-off for gas)
      RcutInBox = ClusterMerger_Rcut;
      for (int d=0; d<3; d++)    RcutInBox = fmin( RcutInBox, 0.5*amr->BoxSize[d] );

      if ( RcutInBox != ClusterMerger_Rcut )
         Aux_Message( stderr, "WARNING : dark matter cut-off radius is reset to %13.7e !!\n", RcutInBox );


//    set particle mass (note that it's set by total dark matter mass within **Rcut** instead of **Rvir**)
      TotM = MassProf_DM( RcutInBox );
      ParM = TotM / amr->Par->NPar_Active_AllRank;

      Aux_Message( stdout, "NOTE : total dark matter mass within the cut-off radius = %13.7e Msun\n",      TotM*UNIT_M/Const_Msun );
      Aux_Message( stdout, "       --> dark matter particle mass is set to          = %13.7e Msun\n",      ParM*UNIT_M/Const_Msun );
      Aux_Message( stdout, "                                                        = %13.7e code unit\n", ParM );


//    set particle attributes
      for (long p=0; p<amr->Par->NPar_Active_AllRank; p++)
      {
//       (1) mass
         Mass_AllRank[p] = ParM;


//       (2) position (sample from the cumulative mass profile)
         RanM = ( (double)rand()/RAND_MAX )*TotM;
         RanR = Mis_InterpolateFromTable( ClusterMerger_NBin_MassProf, ClusterMerger_DM_MassProf,
                                          ClusterMerger_DM_MassProf_R, RanM );
         if ( RanR == NULL_REAL )
         {
            if ( RanM < ClusterMerger_DM_MassProf[0] )
            {
               Aux_Message( stderr, "WARNING : random mass sample (%20.14e) < minimum mass in the interpolation table (%20.14e) !!\n",
                            RanM, ClusterMerger_DM_MassProf[0] );
               Aux_Message( stderr, "          --> Particle radial distance is default to %20.14e\n",
                            ClusterMerger_DM_MassProf_R[0] );
               Aux_Message( stderr, "          --> You might want to set ClusterMerger_IntRmin smaller\n" );

               RanR = ClusterMerger_DM_MassProf_R[0];
            }

            else
               Aux_Error( ERROR_INFO, "random mass sample (%20.14e) >= maximum mass in the interpolation table (%20.14e) !!\n",
                          RanM, ClusterMerger_DM_MassProf[ ClusterMerger_NBin_MassProf - 1 ] );
         }

//       randomly set the position vector with a given radius
         RanVec_FixRadius( RanR, RanVec );
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + BoxCenter[d];


//       (3) velocity
//       interpolate velocity dispersion at the particle position
         Sigma = Mis_InterpolateFromTable( ClusterMerger_NBin_SigmaProf, ClusterMerger_DM_SigmaProf_R,
                                           ClusterMerger_DM_SigmaProf, RanR );

//       RanR must already lie within the interpolation table since it is returned from the table of mass profile
         if ( Sigma == NULL_REAL )
            Aux_Error( ERROR_INFO, "cannot determine velocity dispersion (target radius = %20.14e) !!\n", RanR );

//       assume velocity distribution is isotropic and Gaussian with a standard deviation of Sigma
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = gsl_ran_gaussian( GSL_RNG, Sigma );

      } // for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)

//    free GSL resource
      gsl_rng_free( GSL_RNG );
   } // if ( MPI_Rank == 0 )


// synchronize all particles to the physical time at the base level
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   amr->Par->Time[p] = Time[0];


// get the number of particles in each rank and set the corresponding offsets
   if ( amr->Par->NPar_Active_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                 amr->Par->NPar_Active_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank], NPar_MyRank=(int)amr->Par->NPar_AcPlusInac;

   MPI_Gather( &NPar_MyRank, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particle attributes from the master rank to all ranks
   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

#  ifdef FLOAT8
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, amr->Par->NPar_AcPlusInac, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], amr->Par->NPar_AcPlusInac, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], amr->Par->NPar_AcPlusInac, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

#  else
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, amr->Par->NPar_AcPlusInac, MPI_FLOAT,  0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], amr->Par->NPar_AcPlusInac, MPI_FLOAT,  0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], amr->Par->NPar_AcPlusInac, MPI_FLOAT,  0, MPI_COMM_WORLD );
   }
#  endif


   if ( MPI_Rank == 0 )
   {
      delete [] Mass_AllRank;

      for (int d=0; d<3; d++)
      {
         delete [] Pos_AllRank[d];
         delete [] Vel_AllRank[d];
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r        : Input radius
//                RanVec   : Array to store the random 3D vector
//
// Return      :  RanVec
//-------------------------------------------------------------------------------------------------------
void RanVec_FixRadius( const double r, double RanVec[] )
{

   double Norm, RanR2;

   do
   {
      RanR2 = 0.0;

      for (int d=0; d<3; d++)
      {
         RanVec[d]  = ( (double)rand()/RAND_MAX )*2.0 - 1.0;
         RanR2     += SQR( RanVec[d] );
      }
   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



#endif // #ifdef PARTICLE
