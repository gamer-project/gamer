#include "GAMER.h"

#ifdef PARTICLE

extern int  NParAcc_Mode;
extern int  NParAcc_RSeed;
extern real NParAcc_MaxR;
extern real NParAcc_Rho0;
extern real NParAcc_R0;
extern int  NParAcc_NBinR;

double MassProf_Plummer( const double r );
double MassProf_Burkert( const double r );
double MassProf_NFW    ( const double r );
static void RanVec_FixRadius( const double r, double RanVec[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize the particle position and velocity
//
// Note        :  Invoked by "Init_GAMER"
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass, amr->Par->PosX/Y/Z, amr->Par->VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   for (long p=0; p<amr->Par->NPar_Active_AllRank; p++)  amr->Par->Time[p] = Time[0];


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

   const double Cen[3] = { 0.5*amr->BoxSize[0],
                           0.5*amr->BoxSize[1],
                           0.5*amr->BoxSize[2] };

   double (*MassProf)( const double r ) = NULL;
   double *Table_MassProf_r = NULL;
   double *Table_MassProf_M = NULL;
   double  TotM, ParM, dr, RanM, RanR, rL, rR, ML, MR, EstM, ErrM, ErrM_Max=1.0, RanVec[3];
   int     IdxL, IdxR;


// ============================================================================================================
// set the random seed
   srand( NParAcc_RSeed );


// set the mass profile of the target model
   switch ( NParAcc_Mode )
   {
      case 1:  MassProf = MassProf_Plummer;   break;
      case 2:  MassProf = MassProf_Burkert;   break;
      case 3:  MassProf = MassProf_NFW;       break;

      default: Aux_Error( ERROR_INFO, "unsupported mode (%d) !!\n", NParAcc_Mode );
   }


// determine the total enclosed mass within the maximum radius
   TotM = MassProf( NParAcc_MaxR );
   ParM = TotM / amr->Par->NPar_Active_AllRank;


// construct the mass profile table
   Table_MassProf_r = new double [NParAcc_NBinR];
   Table_MassProf_M = new double [NParAcc_NBinR];

   dr = NParAcc_MaxR / (NParAcc_NBinR-1);

   for (int b=0; b<NParAcc_NBinR; b++)
   {
      Table_MassProf_r[b] = dr*b;
      Table_MassProf_M[b] = MassProf( Table_MassProf_r[b] );
   }


// set particle attributes
   for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)
   {
//    mass
      Mass[p] = ParM;


//    position (sample from the cumulative mass profile and perform linear interpolation)
      RanM = ( (double)rand()/RAND_MAX )*TotM;
      IdxL = Mis_BinarySearch_Real( Table_MassProf_M, 0, NParAcc_NBinR-1, RanM );
      IdxR = IdxL + 1;
      rL   = Table_MassProf_r[IdxL];
      rR   = Table_MassProf_r[IdxR];
      ML   = Table_MassProf_M[IdxL];
      MR   = Table_MassProf_M[IdxR];

//    linear interpolation
      RanR = rL + (rR-rL)/(MR-ML)*(RanM-ML);

//    record the maximum error
      EstM     = MassProf( RanR );
      ErrM     = fabs( (EstM-RanM)/RanM );
      ErrM_Max = fmax( ErrM, ErrM_Max );

//    randomly set the position vector with a given radius
      RanVec_FixRadius( RanR, RanVec );
      for (int d=0; d<3; d++)    Pos[d][p] = RanVec[d] + Cen[d];


//    velocity must be initialized as zero since we use the updated velocity array to store the acceleration
      for (int d=0; d<3; d++)    Vel[d][p] = 0.0;
   } // for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)

   Aux_Message( stdout, "   Total enclosed mass = %13.7e\n", TotM );
   Aux_Message( stdout, "   Particle mass       = %13.7e\n", ParM );
   Aux_Message( stdout, "   Maximum mass error  = %13.7e\n", ErrM_Max );


// free memory
   delete [] Table_MassProf_r;
   delete [] Table_MassProf_M;
// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



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

   const double x = r / NParAcc_R0;

   return 4.0/3.0*M_PI*NParAcc_Rho0*CUBE(r)*pow( 1.0+x*x, -1.5 );

} // FUNCTION : MassProf_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Burkert
// Description :  Mass profile of the Burkert model
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_Burkert( const double r )
{

   const double x = r / NParAcc_R0;

   return M_PI*CUBE(NParAcc_R0)*NParAcc_Rho0*(  -2.0*atan(x) + log( SQR(1.0+x)*(1.0+x*x) )  );

} // FUNCTION : MassProf_Burkert



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_NFW
// Description :  Mass profile of the NFW model
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_NFW( const double r )
{

   const double x = r / NParAcc_R0;

   return 4.0*M_PI*CUBE(NParAcc_R0)*NParAcc_Rho0*( -x/(1.0+x) + log(1.0+x) );

} // FUNCTION : MassProf_NFW



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
