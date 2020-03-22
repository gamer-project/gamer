/***Total deep copy***/
#include "GAMER.h"
#include "NFW_calculator.h"
#include "Plummer_calculator.h"
#include"Burkert_calculator.h"
#include"Jaffe_calculator.h"
#include "Hernquist_calculator.h"
#include "Einasto_calculator.h"
#include"UNKNOWN_calculator.h"

#define DEBUG
#ifdef PARTICLE

extern int*    Models_RSeed;
extern double* Models_Rho0;
extern double* Models_R0;
extern double* Models_MaxR;
extern bool*   Models_Collision;
extern double* Models_Collision_D;
extern double** Models_Center;
extern double** Models_BulkVel;
extern double* Models_GasMFrac;
extern int*    Models_MassProfNBin;

extern int Models_num;
extern string* Models_Paras;
extern string* Models_Type;
extern string* Models_Profile;
extern double* Models_Alpha;
extern int*    Models_r_col;
extern int*    Models_rho_col;
extern bool* Models_truncation;

static RandomNumber_t *RNG = NULL;


static double MassProf_Models( const double r ,string model_name,int k);
static void   RanVec_FixRadius( const double r, double RanVec[] );


//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Models
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
NFW_calculator cal_NFW;//main
Plummer_calculator cal_Plummer;
Burkert_calculator cal_Burkert;
Jaffe_calculator cal_Jaffe;
Hernquist_calculator cal_Hernquist;
Einasto_calculator cal_Einasto;
UNKNOWN_calculator cal_UNKNOWN;
void Par_Init_ByFunction_Models( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] )
{
// Input Model names
   static bool good=1;
   string *model_name;
   static bool *flag;
   if(good){
      model_name=new string[Models_num];
      flag =new bool[Models_num];
      for(int k=0;k<Models_num;k++)flag[k]=0;
      good=0;
   }
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };
   
   if ( MPI_Rank == 0 )
      {  
         Mass_AllRank = new real [NPar_AllRank];
         for (int d=0; d<3; d++)
         {
            Pos_AllRank[d] = new real [NPar_AllRank];
            Vel_AllRank[d] = new real [NPar_AllRank];
         }
      }

   
   for(int k=0;k<Models_num;k++){
   // Initialize calculators
      model_name[k]=Models_Type[k];
      if(flag[k]==0){
         const char* profile=Models_Profile[k].c_str();
         if(model_name[k]=="Plummer")cal_Plummer.init(NEWTON_G,Models_Rho0[k],Models_R0[k]);
         else if(model_name[k]=="NFW")cal_NFW.init(NEWTON_G,Models_Rho0[k],Models_R0[k],Models_MassProfNBin[k],Models_MaxR[k],Models_RSeed[k],Models_truncation[k],0.7);
         else if(model_name[k]=="Burkert")cal_Burkert.init(NEWTON_G,Models_Rho0[k],Models_R0[k],Models_MassProfNBin[k],Models_MaxR[k],Models_RSeed[k],Models_truncation[k],0.7);
         else if(model_name[k]=="Jaffe")cal_Jaffe.init(NEWTON_G,Models_Rho0[k],Models_R0[k],Models_MassProfNBin[k],Models_MaxR[k],Models_RSeed[k],Models_truncation[k],0.5);
         else if(model_name[k]=="Hernquist")cal_Hernquist.init(NEWTON_G,Models_Rho0[k],Models_R0[k],Models_MassProfNBin[k],Models_MaxR[k],Models_RSeed[k],Models_truncation[k],0.7);
         else if(model_name[k]=="Einasto")cal_Einasto.init(Models_Alpha[k],NEWTON_G,Models_Rho0[k],Models_R0[k],Models_MaxR[k]);
         else if(model_name[k]=="UNKNOWN")cal_UNKNOWN.init(NEWTON_G,Models_r_col[k],Models_rho_col[k],profile);
         flag[k]=1;
      }
      

   // only the master rank will construct the initial condition
      if ( MPI_Rank == 0 )
      {  
         const double Coll_Offset = 0.5*Models_Collision_D[k]/sqrt(3.0);

         double *Table_MassProf_r = NULL;
         double *Table_MassProf_M = NULL;
         double  TotM, ParM, dr, RanM, RanR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3];
         double  Vmax, RanV, RanProb, Prob;

   //    initialize the random number generator
         RNG = new RandomNumber_t( 1 );
         RNG->SetSeed( 0, Models_RSeed[k]);


   //    determine the total enclosed mass within the maximum radius
         TotM = MassProf_Models( Models_MaxR[k] ,model_name[k],k);
         ParM = TotM / NPar_AllRank;

         if ( Models_Collision[k] )   ParM *= 2.0;

   //    rescale particle mass to account for the gas contribution
         ParM *= 1.0 - Models_GasMFrac[k];


   //    construct the mass profile table
         Table_MassProf_r = new double [Models_MassProfNBin[k]];
         Table_MassProf_M = new double [Models_MassProfNBin[k]];

         dr = Models_MaxR[k] / (Models_MassProfNBin[k]-1);

         for (int b=0; b<Models_MassProfNBin[k]; b++)
         {
            Table_MassProf_r[b] = dr*b;
            Table_MassProf_M[b] = MassProf_Models( Table_MassProf_r[b] ,model_name[k],k);
         }

         double max=0.0;
   //    set particle attributes
         for (long p=int(NPar_AllRank*k/Models_num); p<int(NPar_AllRank*(k+1)/Models_num); p++)
         {
   //       mass
            Mass_AllRank[p] = ParM;


   //       position
   //       --> sample from the cumulative mass profile with linear interpolation
            RanM = RNG->GetValue( 0, 0.0, 1.0 )*TotM;
            RanR = Mis_InterpolateFromTable( Models_MassProfNBin[k], Table_MassProf_M, Table_MassProf_r, RanM );

   //       record the maximum error
            EstM     = MassProf_Models( RanR ,model_name[k],k);
            ErrM     = fabs( (EstM-RanM)/RanM );
            ErrM_Max = fmax( ErrM, ErrM_Max );

   //       randomly set the position vector with a given radius
            RanVec_FixRadius( RanR, RanVec );
            for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + Models_Center[k][d];

   //       set position offset for the Models collision test
            if ( Models_Collision[k] )
            for (int d=0; d<3; d++)    Pos_AllRank[d][p] += Coll_Offset*( (p<NPar_AllRank/2)?-1.0:+1.0 );

   //       check periodicity
            for (int d=0; d<3; d++)
            {
               if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
                  Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
            }


   //       velocity//main

            double a3=RanR/Models_R0[k];
            if(model_name[k]=="NFW")RanV = cal_NFW.set_vel(a3); 
            else if(model_name[k]=="Plummer")RanV = cal_Plummer.set_vel(a3);   
            else if(model_name[k]=="Burkert")RanV = cal_Burkert.set_vel(a3);
            else if(model_name[k]=="Jaffe")RanV = cal_Jaffe.set_vel(a3);    
            else if(model_name[k]=="Hernquist")RanV = cal_Hernquist.set_vel(a3);      
            else if(model_name[k]=="Einasto")RanV = cal_Einasto.set_vel(a3);     
            else if(model_name[k]=="UNKNOWN")RanV = cal_UNKNOWN.set_vel(RanR);            

   //       randomly set the velocity vector with the given amplitude (RanV*Vmax)
            RanVec_FixRadius( RanV, RanVec );
            for (int d=0; d<3; d++)    Vel_AllRank[d][p] = RanVec[d] + Models_BulkVel[k][d];

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

   }
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
} // FUNCTION : Par_Init_ByFunction_Models



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Models
// Description :  Mass profile of the Models model
//
// Note        :  Calculate the enclosed mass within the given radius
//
// Parameter   :  r : Input radius
//
// Return      :  Enclosed mass
//-------------------------------------------------------------------------------------------------------

double MassProf_Models( const double r ,string model_name,int k)
{
   
   const double x = r / Models_R0[k];
   if     (model_name=="NFW"      )       return cal_NFW.set_mass(x);
   else if(model_name=="Plummer"  )       return 4.0/3.0*M_PI*Models_Rho0[k]*CUBE(r)*pow( 1.0+x*x, -1.5 );
   else if(model_name=="Burkert"  )       return cal_Burkert.set_mass(x);
   else if(model_name=="Jaffe"    )       return cal_Jaffe.set_mass(x);
   else if(model_name=="Hernquist")       return Models_Rho0[k] *pow(Models_R0[k],3) *pow(x/(1+x),2);
   else if(model_name=="Einasto"  )       return cal_Einasto.set_mass(x);
   else if(model_name=="UNKNOWN"  )       return cal_UNKNOWN.set_mass(r);
   
} // FUNCTION : MassProf_Models



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
