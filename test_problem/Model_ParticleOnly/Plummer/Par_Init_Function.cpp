#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE

extern int  Plummer_RSeed;
extern real Plummer_MaxR;
extern real Plummer_Rho0;
extern real Plummer_R0;
extern int  Plummer_NBinR;
extern bool Plummer_Collision;
extern real Plummer_Collision_D;
extern real Plummer_Center[3];
extern real Plummer_BulkVel[3];

double MassProf_Plummer( const double r );
static void RanVec_FixRadius( const double r, double RanVec[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Function
// Description :  Initialize the particle position and velocity 
//
// Note        :  Invoked by "Init_GAMER"
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass, amr->Par->PosX/Y/Z, amr->Par->VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_Function()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// synchronize all particles to the physical time at the base level
   for (long p=0; p<amr->Par->NPar; p++)  amr->Par->Time[p] = Time[0];


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

   const double TotM_Inf    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
   const double Vmax_Fac    = sqrt( 2.0*NEWTON_G*TotM_Inf );
   const double Coll_Offset = 0.5*Plummer_Collision_D/sqrt(3.0);

   double *Table_MassProf_r = NULL;
   double *Table_MassProf_M = NULL;
   double  TotM, ParM, dr, RanM, RanR, rL, rR, ML, MR, EstM, ErrM, ErrM_Max, RanVec[3];
   double  Vmax, RanV, RanProb, Prob;
   int     IdxL, IdxR;


// ============================================================================================================
// set the random seed
   srand( Plummer_RSeed + MPI_Rank*1000000 );


// determine the total enclosed mass within the maximum radius
   TotM = MassProf_Plummer( Plummer_MaxR );   
   ParM = TotM / amr->Par->NPar;

   if ( Plummer_Collision )   ParM *= 2.0;


// construct the mass profile table
   Table_MassProf_r = new double [Plummer_NBinR];
   Table_MassProf_M = new double [Plummer_NBinR];

   dr = Plummer_MaxR / (Plummer_NBinR-1);

   for (int b=0; b<Plummer_NBinR; b++)
   {
      Table_MassProf_r[b] = dr*b;
      Table_MassProf_M[b] = MassProf_Plummer( Table_MassProf_r[b] );
   }


// set particle attributes
   for (int p=0; p<amr->Par->NPar; p++)
   {
//    mass
      Mass[p] = ParM;


//    position (sample from the cumulative mass profile and perform linear interpolation)
      RanM = ( (double)rand()/RAND_MAX )*TotM;
      IdxL = Mis_BinarySearch_Real( Table_MassProf_M, 0, Plummer_NBinR-1, RanM );
      IdxR = IdxL + 1;
      rL   = Table_MassProf_r[IdxL];
      rR   = Table_MassProf_r[IdxR];
      ML   = Table_MassProf_M[IdxL];
      MR   = Table_MassProf_M[IdxR];

//    linear interpolation
      RanR = rL + (rR-rL)/(MR-ML)*(RanM-ML);

//    record the maximum error
      EstM     = MassProf_Plummer( RanR );
      ErrM     = fabs( (EstM-RanM)/RanM );
      ErrM_Max = fmax( ErrM, ErrM_Max );

//    randomly set the position vector with a given radius
      RanVec_FixRadius( RanR, RanVec );
      for (int d=0; d<3; d++)    Pos[d][p] = RanVec[d] + Plummer_Center[d];

//    set position offset for the Plummer collision test
      if ( Plummer_Collision )
      for (int d=0; d<3; d++)    Pos[d][p] += Coll_Offset*( (p<amr->Par->NPar/2)?-1.0:+1.0 );

//    check periodicity
      if ( OPT__BC_POT == BC_POT_PERIODIC )
      for (int d=0; d<3; d++)    Pos[d][p] = FMOD( Pos[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );


//    velocity
//    determine the maximum velocity (the escaping velocity)
      Vmax = Vmax_Fac*pow( SQR(Plummer_R0) + SQR(RanR), -0.25 );

//    randomly determine the velocity amplitude (ref: Aarseth, S. et al. 1974, A&A, 37, 183)
      do
      {
         RanV    = ( (double)rand()/RAND_MAX );          // 0.0 ... 1.0
         RanProb = ( (double)rand()/RAND_MAX )*0.1;      // 0.0 ... 0.1
         Prob    = SQR(RanV)*pow( 1.0-SQR(RanV), 3.5 );  // < 0.1
      }
      while ( RanProb > Prob );

//    randomly set the velocity vector with the given amplitude (RanV*Vmax)
      RanVec_FixRadius( RanV*Vmax, RanVec );
      for (int d=0; d<3; d++)    Vel[d][p] = RanVec[d] + Plummer_BulkVel[d];

   } // for (int p=0; p<amr->Par->NPar; p++)

   Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
   Aux_Message( stdout, "   Total enclosed mass to inifinity = %13.7e\n",  TotM_Inf );
   Aux_Message( stdout, "   Enclosed ratio                   = %6.2f%%\n", 100.0*TotM/TotM_Inf );
   Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
   Aux_Message( stdout, "   Maximum mass error               = %13.7e\n",  ErrM_Max );


// free memory
   delete [] Table_MassProf_r;
   delete [] Table_MassProf_M;
// ============================================================================================================
 

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_Function



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Plummer
// Description :  Mass profile of the Plummer model
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_Plummer( const double r )
{

   const double x = r / Plummer_R0;

   return 4.0/3.0*M_PI*Plummer_Rho0*CUBE(r)*pow( 1.0+x*x, -1.5 );

} // FUNCTION : MassProf_Plummer



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
