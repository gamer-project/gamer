/***Total deep copy***/
#include "GAMER.h"
#include "Burkert_calculator.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#define DEBUG
#ifdef PARTICLE

extern int    Burkert_RSeed;
extern double Burkert_Rho0;
extern double Burkert_R0;
extern double Burkert_MaxR;
extern bool   Burkert_Collision;
extern double Burkert_Collision_D;
extern double Burkert_Center[3];
extern double Burkert_BulkVel[3];
extern double Burkert_GasMFrac;
extern int    Burkert_MassProfNBin;

static RandomNumber_t *RNG = NULL;


static double MassProf_Burkert( const double r );
static void   RanVec_FixRadius( const double r, double RanVec[] );

double f(double x,void *params);
double function01(double x,void *params);
double function02(double x,void *params);

/***gsl function***/
double Burkert_potential(double a3){
   double s=(4*M_PI*NEWTON_G*Burkert_Rho0*Burkert_R0*Burkert_R0)*(atan(a3)/(2*a3)-log(a3*a3+1)/(4*a3)-log(a3+1)/(2*a3)+atan(a3)/2+log(a3*a3+1)/4-log(a3+1)/2-M_PI/4);
   return s;
}/*
double f (double x, void * parameters) {
   double a3=x;
   double e= *(double *)parameters;
   double fs= -potential(x)-e;
   return fs;
}
double function01 (double x, void * parameters) {
      double rho_dx=2*x*Burkert_Rho0/((x*x+1)*(x*x+1)*(x+1))+Burkert_Rho0/((x*x+1)*(x+1)*(x+1));
      double psi_dx=4*M_PI*NEWTON_G*Burkert_R0*Burkert_Rho0*(-0.5*atan(x)+0.25*log(x*x+1)+0.5*log(x+1))/(x*x);
      double fs = rho_dx/psi_dx;
      return fs;
      }
double de_rho_over_de_psi(double x) {
      gsl_function F;
      double a3=x;
      double result, abserr;
      F.function = &function01;
      gsl_deriv_central (&F, x, 1e-8, &result, &abserr);
      return result;
}
double function02(double x, void * parameters) {
      
      double e= *(double *)parameters;
            
      double fs=de_rho_over_de_psi(x)*pow(e+potential(x),-0.5);
      return fs;
      
} 
double function02(double x, double e) {    
      double fs=de_rho_over_de_psi(x)*pow(e+potential(x),-0.5);
      Burkert_calculator a;
      double *b;
      double s=a.test2(3,b);
      return fs;
     
} */ 
       
/***gsl function***/

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Burkert
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
Burkert_calculator a;
void Par_Init_ByFunction_Burkert( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] )
{
   static bool flag=0;
   if(flag==0){
      a.init();
      flag=1;
      cout<<"done"<<endl;
   }
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      //const double TotM_Inf    = 4.0/3.0*M_PI*CUBE(Burkert_R0)*Burkert_Rho0;
      //const double Vmax_Fac    = sqrt( 2.0*NEWTON_G*TotM_Inf );
      const double Coll_Offset = 0.5*Burkert_Collision_D/sqrt(3.0);

      double *Table_MassProf_r = NULL;
      double *Table_MassProf_M = NULL;
      double  TotM, ParM, dr, RanM, RanR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3];
      double  Vmax, RanV, RanProb, Prob;

      Mass_AllRank = new real [NPar_AllRank];
      for (int d=0; d<3; d++)
      {
         Pos_AllRank[d] = new real [NPar_AllRank];
         Vel_AllRank[d] = new real [NPar_AllRank];
      }


//    initialize the random number generator
      RNG = new RandomNumber_t( 1 );
      RNG->SetSeed( 0, Burkert_RSeed );


//    determine the total enclosed mass within the maximum radius
      TotM = MassProf_Burkert( Burkert_MaxR );
      ParM = TotM / NPar_AllRank;

      if ( Burkert_Collision )   ParM *= 2.0;

//    rescale particle mass to account for the gas contribution
      ParM *= 1.0 - Burkert_GasMFrac;


//    construct the mass profile table
      Table_MassProf_r = new double [Burkert_MassProfNBin];
      Table_MassProf_M = new double [Burkert_MassProfNBin];

      dr = Burkert_MaxR / (Burkert_MassProfNBin-1);

      for (int b=0; b<Burkert_MassProfNBin; b++)
      {
         Table_MassProf_r[b] = dr*b;
         Table_MassProf_M[b] = MassProf_Burkert( Table_MassProf_r[b] );
      }

      double max=0.0;
//    set particle attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
//       mass
         Mass_AllRank[p] = ParM;


//       position
//       --> sample from the cumulative mass profile with linear interpolation
         RanM = RNG->GetValue( 0, 0.0, 1.0 )*TotM;
         RanR = Mis_InterpolateFromTable( Burkert_MassProfNBin, Table_MassProf_M, Table_MassProf_r, RanM );

//       record the maximum error
         EstM     = MassProf_Burkert( RanR );
         ErrM     = fabs( (EstM-RanM)/RanM );
         ErrM_Max = fmax( ErrM, ErrM_Max );

//       randomly set the position vector with a given radius
         RanVec_FixRadius( RanR, RanVec );
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + Burkert_Center[d];

//       set position offset for the Burkert collision test
         if ( Burkert_Collision )
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] += Coll_Offset*( (p<NPar_AllRank/2)?-1.0:+1.0 );

//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
         }


//       velocity

         double a3=RanR/Burkert_R0;
         RanV = a.set_vel(a3);      

//       randomly set the velocity vector with the given amplitude (RanV*Vmax)
         RanVec_FixRadius( RanV, RanVec );
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = RanVec[d] + Burkert_BulkVel[d];

      } // for (long p=0; p<NPar_AllRank; p++)
      

      Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
      //Aux_Message( stdout, "   Total enclosed mass to inifinity = %13.7e\n",  TotM_Inf );
      //Aux_Message( stdout, "   Enclosed mass ratio              = %6.2f%%\n", 100.0*TotM/TotM_Inf );
      Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
      Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  ErrM_Max );


//    free memory
      delete [] Table_MassProf_r;
      delete [] Table_MassProf_M;
   } // if ( MPI_Rank == 0 )


// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)   ParTime[p] = Time[0];


// get the number of particles in each rank and set the corresponding offsets
   if ( NPar_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                 NPar_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank];
   int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                             // --> (ii) is necessary for OpenMPI version < 1.7

   MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particle attributes from the master rank to all ranks
   real *Mass   =   ParMass;
   real *Pos[3] = { ParPosX, ParPosY, ParPosZ };
   real *Vel[3] = { ParVelX, ParVelY, ParVelZ };

#  ifdef FLOAT8
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

#  else
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
   }
#  endif


   if ( MPI_Rank == 0 )
   {
      delete RNG;
      delete [] Mass_AllRank;

      for (int d=0; d<3; d++)
      {
         delete [] Pos_AllRank[d];
         delete [] Vel_AllRank[d];
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Burkert



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Burkert
// Description :  Mass profile of the Burkert model
//
// Note        :  Calculate the enclosed mass within the given radius
//
// Parameter   :  r : Input radius
//
// Return      :  Enclosed mass
//-------------------------------------------------------------------------------------------------------

double MassProf_Burkert( const double r )
{
   
   const double x = r / Burkert_R0;
   //const double m0=1.6*Burkert_Rho0*pow(Burkert_R0,3);

   //return 4.0*m0*(log(1+x)+0.5*log(1+x*x)-atan(x));
   
   return M_PI *Burkert_Rho0 *pow(Burkert_R0,3) *(log(1+x*x) + 2 * log(1+x) -2 *atan(x));


} // FUNCTION : MassProf_Burkert



//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r      : Input radius
//                RanVec : Array to store the random 3D vector
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
         RanVec[d]  = RNG->GetValue( 0, -1.0, +1.0 );
         RanR2     += SQR( RanVec[d] );
      }
   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



#endif // #ifdef PARTICLE
