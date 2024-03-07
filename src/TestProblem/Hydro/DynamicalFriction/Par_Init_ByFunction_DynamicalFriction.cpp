#include "GAMER.h"
#include "Par_EquilibriumIC.h"
#include <tuple>
#include <vector>

using namespace std;



#ifdef MASSIVE_PARTICLES

extern double GC_xx;
extern double GC_yy;
extern double GC_zz;

static RandomNumber_t *RNG = NULL;

static char Table_Name[MAX_STRING];
static double GC_MASS;
static double GC_POSX;
static double GC_POSY;
static double GC_POSZ;
static double GC_VELX;
static double GC_VELY;
static double GC_VELZ;

double GC_R;
bool Pure_Table;

static char HaloType[MAX_STRING];
double Halo_Rho0;
double Halo_Rs;
double Halo_Rt;




vector<double> GC_pos, GC_vel, Halo_vel;

tuple<vector<double>, vector<double>, vector<double>> Calculate_IC( const char* HaloType, const double GC_MASS, const double GC_R,
								    const double Halo_Rho0, const double Halo_Rs, const double Halo_Rt,
								    const char* Table_Name, const bool Pure_Table);

//#include "global_var.h"

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_DynamicalFriction
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
//                ParType       : Particle type     array with the size of NPar_ThisRan
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------

void Par_Init_ByFunction_DynamicalFriction( const long NPar_ThisRank, const long NPar_AllRank,
                                   real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                   real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                   real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// define the particle attribute arrays
   real *ParData_AllRank[PAR_NATT_TOTAL];
   real *new_ParData_AllRank[PAR_NATT_TOTAL];
   for (int v=0; v<PAR_NATT_TOTAL; v++){
         ParData_AllRank[v] = NULL;
         new_ParData_AllRank[v] = NULL;}


// define the particle IC constructor
   Par_EquilibriumIC Filename_Loader;


// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 ) {
      

// load run-time parameters
       const char* FileName = "Input__TestProb";
       ReadPara_t *ReadPara  = new ReadPara_t;
    // ********************************************************************************************************************************
    // ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
    // ********************************************************************************************************************************

    // For GC
       ReadPara->Add( "GC_MASS",                 &GC_MASS,               NoDef_double,  NoMin_double,     NoMax_double      );
       ReadPara->Add( "GC_initial_R",            &GC_R,                  NoDef_double,  NoMin_double,     NoMax_double      );

       ReadPara->Add( "GC_POSX",                 &GC_POSX,               NoDef_double,  NoMin_double,     NoMax_double      );  
       ReadPara->Add( "GC_POSY",                 &GC_POSY,               NoDef_double,  NoMin_double,     NoMax_double      );
       ReadPara->Add( "GC_POSZ",                 &GC_POSZ,               NoDef_double,  NoMin_double,     NoMax_double      );

       ReadPara->Add( "GC_VELX",                 &GC_VELX,               NoDef_double,  NoMin_double,     NoMax_double      );
       ReadPara->Add( "GC_VELY",                 &GC_VELY,               NoDef_double,  NoMin_double,     NoMax_double      ); 
       ReadPara->Add( "GC_VELZ",                 &GC_VELZ,               NoDef_double,  NoMin_double,     NoMax_double      );
    // Halo
       ReadPara->Add( "HALO_TYPE",               HaloType,            "None",        Useless_str,   Useless_str          );
       ReadPara->Add( "HALO_RHO_0",              &Halo_Rho0,             NoDef_double,  NoMin_double,     NoMax_double      );  
       ReadPara->Add( "HALO_Rs",                 &Halo_Rs,               NoDef_double,  NoMin_double,     NoMax_double      );  
       ReadPara->Add( "HALO_Rt",                 &Halo_Rt,               NoDef_double,  NoMin_double,     NoMax_double      );  
 
    // PURE_TABLE
       ReadPara->Add( "PURE_TABLE",              &Pure_Table,            false,          Useless_bool,     Useless_bool      );  
       ReadPara->Add( "Density_Table_Name",      Table_Name,             "Profile_Table.txt",        Useless_str,   Useless_str          );
      
       ReadPara->Read( FileName );

    // call the function to calculate the velocity, generate density table,        
      
       tie(GC_pos, GC_vel, Halo_vel) = Calculate_IC(HaloType, GC_MASS, GC_R,Halo_Rho0,Halo_Rs,Halo_Rt,Table_Name,Pure_Table);

       Aux_Message(stdout, " GC_pos %f, %f, %f \n", GC_pos[0],GC_pos[1],GC_pos[2]);
       Aux_Message(stdout, " GC_vel %f, %f, %f \n", GC_vel[0],GC_vel[1],GC_vel[2]);
       Aux_Message(stdout, " Halo_vel %f, %f, %f \n", Halo_vel[0],Halo_vel[1],Halo_vel[2]);
       Aux_Message(stdout, "%f, %s \n", GC_MASS,HaloType);












//    allocate memory for particle attribute arrays
      ParData_AllRank[PAR_MASS] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSZ] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELZ] = new real [NPar_AllRank];
      ParData_AllRank[PAR_TYPE] = new real [NPar_AllRank];
//    input filenames as parameters into Filename_Loader
      Filename_Loader.Read_Filenames( "Input__TestProb" );
      long Par_Idx0 = 0;


//       initialize Par_EquilibriumIC for each cloud // There is one cloud only, we're just using the ParticleEquilibriem Testprob and modify it
         Par_EquilibriumIC Cloud_Constructor;
         Cloud_Constructor.Load_Physical_Params( Filename_Loader.filenames, 0, NPar_AllRank-1 );
         Cloud_Constructor.Init();

//      check whether the particle number of each cloud is reasonable
        if ( (Par_Idx0 + Cloud_Constructor.params.Cloud_Par_Num) > NPar_AllRank-1 ) {
           Aux_Error( ERROR_INFO, "particle number doesn't match (%ld + %ld = %ld > %ld) !!\n",
                       Par_Idx0, Cloud_Constructor.params.Cloud_Par_Num, Par_Idx0+Cloud_Constructor.params.Cloud_Par_Num, NPar_AllRank-1 );
        }
   
//       set an equilibrium initial condition for each cloud
         Cloud_Constructor.Par_SetEquilibriumIC( ParData_AllRank[PAR_MASS], ParData_AllRank+PAR_POSX, ParData_AllRank+PAR_VELX, Par_Idx0 );

   
//	set the particle type of the cloud

//       update the particle index offset for the next cloud
         Par_Idx0 += Cloud_Constructor.params.Cloud_Par_Num; 
//      set the particle type of the cloud
	for (int k=0;k<Par_Idx0;k++)ParData_AllRank[PAR_TYPE][k]=PTYPE_GENERIC_MASSIVE;    


//       ParData_AllRank[PAR_MASS][Par_Idx0] = GC_MASS;
//       ParData_AllRank[PAR_POSX][Par_Idx0] = GC_POSX;
//       ParData_AllRank[PAR_POSY][Par_Idx0] = GC_POSY;
//       ParData_AllRank[PAR_POSZ][Par_Idx0] = GC_POSZ;
//       ParData_AllRank[PAR_VELX][Par_Idx0] = GC_VELX;
//       ParData_AllRank[PAR_VELY][Par_Idx0] = GC_VELY;
//       ParData_AllRank[PAR_VELZ][Par_Idx0] = GC_VELZ;
       ParData_AllRank[PAR_MASS][Par_Idx0] = GC_MASS;
       ParData_AllRank[PAR_POSX][Par_Idx0] = GC_pos[0];
       ParData_AllRank[PAR_POSY][Par_Idx0] = GC_pos[1];
       ParData_AllRank[PAR_POSZ][Par_Idx0] = GC_pos[2];
       ParData_AllRank[PAR_VELX][Par_Idx0] = GC_vel[0];
       ParData_AllRank[PAR_VELY][Par_Idx0] = GC_vel[1];
       ParData_AllRank[PAR_VELZ][Par_Idx0] = GC_vel[2];
       ParData_AllRank[PAR_TYPE][Par_Idx0] = PTYPE_GC;
       // Broadcast the GC information from the root rank
       //GC_xx = GC_POSX;
       //GC_yy = GC_POSY;
       //GC_zz = GC_POSZ;


      

       

           

  } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL|_PAR_TYPE, ParData_AllRank, AllAttribute );


// synchronize all particles to the physical time on the base level
// and assign particle type
   for (long p=0; p<NPar_ThisRank; p++) {
      ParTime[p] = Time[0];
   }


// free resource
   if ( MPI_Rank == 0 )
   {
      delete RNG;
      for (int v=0; v<PAR_NATT_TOTAL; v++){
         delete [] ParData_AllRank[v];
      }
   }

// Broadcast the GC position from the root rank

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&GC_xx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&GC_yy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&GC_zz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
   

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_DynamicalFriction

#endif // #ifdef MASSIVE_PARTICLES

