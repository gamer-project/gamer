#include "GAMER.h"

#ifdef PARTICLE

extern int    ParTest_NPar[3];
extern double ParTest_Point_Mass;
extern double ParTest_Par_Sep;
extern bool   ParTest_Use_Tracers;
extern bool   ParTest_Use_Massive;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_ParticleTest
// Description :  Initialize all particle attributes for the particle test problem
//                --> Modified from "Par_Init_ByFile.cpp"
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
//                4. File format: plain C binary in the format [Number of particles][Particle attributes]
//                   --> [Particle 0][Attribute 0], [Particle 0][Attribute 1], ...
//                   --> Note that it's different from the internal data format in the particle repository,
//                       which is [Particle attributes][Number of particles]
//                   --> Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                       (and exactly in this order)
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_ParticleTest( const long NPar_ThisRank, const long NPar_AllRank,
                                       real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                       real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                       real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   long NPar_All = 0;
   if ( ParTest_Use_Massive ) NPar_All += 2;
   if ( ParTest_Use_Tracers ) NPar_All += ParTest_NPar[0]*ParTest_NPar[1]*ParTest_NPar[2];

   if ( NPar_All != NPar_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found [%ld] != expect [%ld] !!\n",
                 NPar_All, NPar_AllRank );

// define the particle attribute arrays
   real *ParData_AllRank[PAR_NATT_TOTAL];
   for (int v=0; v<PAR_NATT_TOTAL; v++)   ParData_AllRank[v] = NULL;

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 ) {

//    allocate memory for particle attribute arrays
      ParData_AllRank[PAR_MASS] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSZ] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELZ] = new real [NPar_AllRank];
      ParData_AllRank[PAR_TYPE] = new real [NPar_AllRank];

      long p = 0;

      if ( ParTest_Use_Massive ) {

         const double v = 0.5*SQRT(Const_NewtonG*ParTest_Point_Mass/(0.5*ParTest_Par_Sep));

         for (int ii=0; ii<2; ii++) {

            const double dir = 2.0*ii-1.0;

            ParData_AllRank[PAR_MASS][p] = real( ParTest_Point_Mass );

            ParData_AllRank[PAR_POSX][p] = real( 0.5*amr->BoxSize[0] +
                 0.5*ParTest_Par_Sep*dir );
            ParData_AllRank[PAR_POSY][p] = real( 0.5*amr->BoxSize[1] );
            ParData_AllRank[PAR_POSZ][p] = real( 0.5*amr->BoxSize[2] );

            ParData_AllRank[PAR_VELX][p] = 0.0;
            ParData_AllRank[PAR_VELY][p] = real( v*dir );
            ParData_AllRank[PAR_VELZ][p] = 0.0;

//          set the particle type to be generic massive
            ParData_AllRank[PAR_TYPE][p] = PTYPE_GENERIC_MASSIVE;

            p++;
         }
      } // if ( ParTest_Use_Massive )

      if ( ParTest_Use_Tracers ) {

         const double delta_p[3] = { 0.5*amr->BoxSize[0]/(ParTest_NPar[0]+1),
                                     0.5*amr->BoxSize[1]/(ParTest_NPar[1]+1),
                                     0.5*amr->BoxSize[2]/(ParTest_NPar[2]+1) };

         for (long kk=0; kk<ParTest_NPar[2]; kk++)
         for (long jj=0; jj<ParTest_NPar[1]; jj++)
         for (long ii=0; ii<ParTest_NPar[0]; ii++)
         {

//          tracer particles have no mass
            ParData_AllRank[PAR_MASS][p] = 0.0;

//          assign positions
            ParData_AllRank[PAR_POSX][p] = real(
                 (ii+1)*delta_p[0]+0.25*amr->BoxSize[0] );
            ParData_AllRank[PAR_POSY][p] = real(
                 (jj+1)*delta_p[1]+0.25*amr->BoxSize[1] );
            ParData_AllRank[PAR_POSZ][p] = real(
                 (kk+1)*delta_p[2]+0.25*amr->BoxSize[2] );

//          set velocities to zero (these will be updated from the grid later)
            ParData_AllRank[PAR_VELX][p] = 0.0;
            ParData_AllRank[PAR_VELY][p] = 0.0;
            ParData_AllRank[PAR_VELZ][p] = 0.0;

//          set the particle type to be tracer
            ParData_AllRank[PAR_TYPE][p] = PTYPE_TRACER;

            p++;
         }
      } // if ( ParTest_Use_Tracers )
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank,
                            _PAR_MASS|_PAR_POS|_PAR_VEL|_PAR_TYPE,
                            ParData_AllRank, AllAttribute );

// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)
      ParTime[p] = Time[0];

// free resource
   if ( MPI_Rank == 0 )
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParData_AllRank[v];


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_ParticleTest



#endif // #ifdef PARTICLE
