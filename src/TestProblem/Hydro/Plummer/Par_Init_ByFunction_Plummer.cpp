#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

extern int    Plummer_RSeed;
extern double Plummer_Rho0;
extern double Plummer_R0;
extern double Plummer_MaxR;
extern bool   Plummer_Collision;
extern double Plummer_Collision_D;
extern double Plummer_Center[3];
extern double Plummer_BulkVel[3];
extern double Plummer_GasMFrac;
extern double Plummer_ParMFrac;
extern double Plummer_ExtAccMFrac;
extern double Plummer_ExtPotMFrac;
extern int    Plummer_MassProfNBin;

static RandomNumber_t *RNG = NULL;


static double MassProf_Plummer( const double r );
static void   RanVec_FixRadius( const double r, double RanVec[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Plummer
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
// Parameter   :  NPar_ThisRank   : Number of particles to be set by this MPI rank
//                NPar_AllRank    : Total Number of particles in all MPI ranks
//                ParMass         : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z     : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z     : Particle velocity array with the size of NPar_ThisRank
//                ParTime         : Particle time     array with the size of NPar_ThisRank
//                ParType         : Particle type     array with the size of NPar_ThisRank
//                AllAttributeFlt : Pointer array for all particle floating-point attributes
//                                  --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                      to access the data
//                AllAttributeInt : Pointer array for all particle integer attributes
//                                  --> Dimension = [PAR_NATT_INT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttributeFlt, AllAttributeInt
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Plummer( const long NPar_ThisRank, const long NPar_AllRank,
                                  real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                  real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                  long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                  long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real_par *ParFltData_AllRank[PAR_NATT_FLT_TOTAL];
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   ParFltData_AllRank[v] = NULL;
   long_par *ParIntData_AllRank[PAR_NATT_INT_TOTAL];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   ParIntData_AllRank[v] = NULL;

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      const double TotM_Inf    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
      const double Vmax_Fac    = sqrt( 2.0*NEWTON_G*TotM_Inf );
      const double Coll_Offset = 0.5*Plummer_Collision_D/sqrt(3.0);

      double *Table_MassProf_r = NULL;
      double *Table_MassProf_M = NULL;
      double  TotM, ParM, dr, RanM, RanR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3];
      double  Vmax, RanV, RanProb, Prob;

      ParFltData_AllRank[PAR_MASS] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSZ] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELZ] = new real_par [NPar_AllRank];


//    initialize the random number generator
      RNG = new RandomNumber_t( 1 );
      RNG->SetSeed( 0, Plummer_RSeed );


//    determine the total enclosed mass within the maximum radius
      TotM = MassProf_Plummer( Plummer_MaxR );
      ParM = TotM / NPar_AllRank;

      if ( Plummer_Collision )   ParM *= 2.0;

//    rescale particle mass to account for the contribution from gas and external gravity
      ParM *= Plummer_ParMFrac;


//    construct the mass profile table
      Table_MassProf_r = new double [Plummer_MassProfNBin];
      Table_MassProf_M = new double [Plummer_MassProfNBin];

      dr = Plummer_MaxR / (Plummer_MassProfNBin-1);

      for (int b=0; b<Plummer_MassProfNBin; b++)
      {
         Table_MassProf_r[b] = dr*b;
         Table_MassProf_M[b] = MassProf_Plummer( Table_MassProf_r[b] );
      }


//    set particle attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
//       mass
         ParFltData_AllRank[PAR_MASS][p] = ParM;


//       position
//       --> sample from the cumulative mass profile with linear interpolation
         RanM = RNG->GetValue( 0, 0.0, 1.0 )*TotM;
         RanR = Mis_InterpolateFromTable( Plummer_MassProfNBin, Table_MassProf_M, Table_MassProf_r, RanM );

//       record the maximum error
         EstM     = MassProf_Plummer( RanR );
         ErrM     = fabs( (EstM-RanM)/RanM );
         ErrM_Max = fmax( ErrM, ErrM_Max );

//       randomly set the position vector with a given radius
         RanVec_FixRadius( RanR, RanVec );
         for (int d=0; d<3; d++)    ParFltData_AllRank[PAR_POSX+d][p] = real_par(RanVec[d] + Plummer_Center[d]);

//       set position offset for the Plummer collision test
         if ( Plummer_Collision )
         for (int d=0; d<3; d++)    ParFltData_AllRank[PAR_POSX+d][p] += real_par(Coll_Offset*( (p<NPar_AllRank/2)?-1.0:+1.0 ));

//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               ParFltData_AllRank[PAR_POSX+d][p] = FMOD( ParFltData_AllRank[PAR_POSX+d][p]+(real_par)amr->BoxSize[d], (real_par)amr->BoxSize[d] );
         }


//       velocity
//       determine the maximum velocity (i.e., the escaping velocity)
         Vmax = Vmax_Fac*pow( SQR(Plummer_R0) + SQR(RanR), -0.25 );

//       randomly determine the velocity amplitude (ref: Aarseth, S. et al. 1974, A&A, 37, 183: Eq. [A4,A5])
         do
         {
            RanV    = RNG->GetValue( 0, 0.0, 1.0 );         // (0.0, 1.0)
            RanProb = RNG->GetValue( 0, 0.0, 0.1 );         // (0.0, 0.1)
            Prob    = SQR(RanV)*pow( 1.0-SQR(RanV), 3.5 );  // < 0.1
         }
         while ( RanProb > Prob );

//       randomly set the velocity vector with the given amplitude (RanV*Vmax)
         RanVec_FixRadius( RanV*Vmax, RanVec );
         for (int d=0; d<3; d++)    ParFltData_AllRank[PAR_VELX+d][p] = real_par(RanVec[d] + Plummer_BulkVel[d]);

      } // for (long p=0; p<NPar_AllRank; p++)

      Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
      Aux_Message( stdout, "   Total enclosed mass to inifinity = %13.7e\n",  TotM_Inf );
      Aux_Message( stdout, "   Enclosed mass ratio              = %6.2f%%\n", 100.0*TotM/TotM_Inf );
      Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
      Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  ErrM_Max );


//    free memory
      delete [] Table_MassProf_r;
      delete [] Table_MassProf_M;
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, _NONE,
                            ParFltData_AllRank, ParIntData_AllRank, AllAttributeFlt, AllAttributeInt );

// synchronize all particles to the physical time on the base level,
// set generic particle type
   for (long p=0; p<NPar_ThisRank; p++) {
      ParTime[p] = (real_par)Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }

// free resource
   if ( MPI_Rank == 0 )
   {
      delete RNG;
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] ParFltData_AllRank[v];
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] ParIntData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Plummer
// Description :  Mass profile of the Plummer model
//
// Note        :  Calculate the enclosed mass within the given radius
//
// Parameter   :  r : Input radius
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



#endif // #ifdef MASSIVE_PARTICLES
