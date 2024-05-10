#include "GAMER.h"
#include "Par_EquilibriumIC.h"

#ifdef MASSIVE_PARTICLES

static RandomNumber_t *RNG = NULL;

extern int      ParEqmIC_NumCloud;
extern char   (*ParEqmIC_Cloud_ParaFilenames)[MAX_STRING];

extern double (*ParEqmIC_Cloud_Center)[3];
extern double (*ParEqmIC_Cloud_BulkVel)[3];
extern char   (*ParEqmIC_Cloud_Type)[MAX_STRING];
extern double  *ParEqmIC_Cloud_Rho0;
extern double  *ParEqmIC_Cloud_R0;
extern double  *ParEqmIC_Cloud_EinastoPowerFactor;
extern char   (*ParEqmIC_Cloud_DensityTable)[MAX_STRING];
extern long    *ParEqmIC_Cloud_ParNum;
extern double  *ParEqmIC_Cloud_MaxR;
extern int     *ParEqmIC_Cloud_DensProfNBin;
extern int     *ParEqmIC_Cloud_EnergyNBin;
extern int     *ParEqmIC_Cloud_RSeed;
extern int     *ParEqmIC_Cloud_AddExtPotAnaly;
extern int     *ParEqmIC_Cloud_AddExtPotTable;
extern char   (*ParEqmIC_Cloud_ExtPotTable)[MAX_STRING];




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_ParEqmIC
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

void Par_Init_ByFunction_ParEqmIC( const long NPar_ThisRank, const long NPar_AllRank,
                                   real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                   real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                   real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// define the particle attribute arrays
   real *ParData_AllRank[PAR_NATT_TOTAL];
   for (int v=0; v<PAR_NATT_TOTAL; v++)   ParData_AllRank[v] = NULL;


// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {

//    allocate memory for particle attribute arrays
      ParData_AllRank[PAR_MASS] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSZ] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELZ] = new real [NPar_AllRank];

      long Par_Idx0 = 0;

      for (int i=0; i<ParEqmIC_NumCloud; i++)
      {

         // Check whether the particle number of each cloud is reasonable
         if ( (Par_Idx0 + ParEqmIC_Cloud_ParNum[i]) > NPar_AllRank )
         {
            Aux_Error( ERROR_INFO, "particle number doesn't match (%ld + %ld = %ld > %ld) !!\n",
                       Par_Idx0, ParEqmIC_Cloud_ParNum[i], Par_Idx0+ParEqmIC_Cloud_ParNum[i], NPar_AllRank );
         }

         // Initialize Par_EquilibriumIC for each cloud
         Par_EquilibriumIC Cloud_Constructor( ParEqmIC_Cloud_Type[i] );

         // Set the parameters for each particle cloud
         Cloud_Constructor.setCenterAndBulkVel(      ParEqmIC_Cloud_Center            [i][0],
                                                     ParEqmIC_Cloud_Center            [i][1],
                                                     ParEqmIC_Cloud_Center            [i][2],
                                                     ParEqmIC_Cloud_BulkVel           [i][0],
                                                     ParEqmIC_Cloud_BulkVel           [i][1],
                                                     ParEqmIC_Cloud_BulkVel           [i][2] );
         Cloud_Constructor.setModelParameters(       ParEqmIC_Cloud_Rho0              [i],
                                                     ParEqmIC_Cloud_R0                [i] );
         Cloud_Constructor.setEinastoPowerFactor(    ParEqmIC_Cloud_EinastoPowerFactor[i] );
         Cloud_Constructor.setDensProfTableFilename( ParEqmIC_Cloud_DensityTable      [i] );
         Cloud_Constructor.setParticleParameters(    ParEqmIC_Cloud_ParNum            [i],
                                                     ParEqmIC_Cloud_MaxR              [i],
                                                     ParEqmIC_Cloud_DensProfNBin      [i],
                                                     ParEqmIC_Cloud_EnergyNBin        [i],
                                                     ParEqmIC_Cloud_RSeed             [i] );
         Cloud_Constructor.setExternalPotential(     ParEqmIC_Cloud_AddExtPotAnaly    [i],
                                                     ParEqmIC_Cloud_AddExtPotTable    [i],
                                                     ParEqmIC_Cloud_ExtPotTable       [i] );

         // initialize the particle cloud
         Cloud_Constructor.initialize();

//       set an equilibrium initial condition for each cloud
         Cloud_Constructor.constructParticles( ParData_AllRank[PAR_MASS], ParData_AllRank+PAR_POSX, ParData_AllRank+PAR_VELX, Par_Idx0 );

         Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  Cloud_Constructor.getTotCloudMass() );
         Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  Cloud_Constructor.getParticleMass() );
         Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  Cloud_Constructor.getMaxMassError() );

//       update the particle index offset for the next cloud
         Par_Idx0 += ParEqmIC_Cloud_ParNum[i];

      } // for (int i=0; i<ParEqmIC_CloudNum; i++)

      // check whether the total particle number is reasonable
      if ( Par_Idx0 != NPar_AllRank )
      {
         Aux_Error( ERROR_INFO, "total particle number doesn't match (total = %ld > NPar_AllRank = %ld) !!\n", Par_Idx0, NPar_AllRank );
      }

   } // if ( MPI_Rank == 0 )


// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, ParData_AllRank, AllAttribute );


// synchronize all particles to the physical time on the base level
// and assign particle type
   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParTime[p] = Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }


// free resource
   if ( MPI_Rank == 0 )
   {
      delete RNG;
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_ParEqmIC

#endif // #ifdef MASSIVE_PARTICLES
