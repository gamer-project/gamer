#include "GAMER.h"

#ifdef PARTICLE

extern int ParFlag_MinLv;
extern int ParFlag_MaxLv;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_ParFlag
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
void Par_Init_ByFunction_ParFlag( const long NPar_ThisRank, const long NPar_AllRank,
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
      ParFltData_AllRank[PAR_MASS] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSZ] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELZ] = new real_par [NPar_AllRank];
      ParIntData_AllRank[PAR_FLAG] = new long_par [NPar_AllRank];

//    set general attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
         ParFltData_AllRank[PAR_MASS][p] = 1.0;
         if ( p < 3 ) {
         ParFltData_AllRank[PAR_POSX][p] = amr->BoxEdgeL [0] + 1.5*amr->dh[0];
         ParFltData_AllRank[PAR_POSY][p] = amr->BoxEdgeL [1] + (p*2.0+1.0)/6.0*amr->BoxSize[1];
         }
         ParFltData_AllRank[PAR_POSZ][p] = amr->BoxCenter[2] + 0.5*amr->dh[0];
         ParFltData_AllRank[PAR_VELX][p] = 1.0;
         ParFltData_AllRank[PAR_VELY][p] = 0.0;
         ParFltData_AllRank[PAR_VELZ][p] = 0.0;
      }

//    set particle flags
//    particle 0: refine to the minimum level
      ParIntData_AllRank[PAR_FLAG][0] = ( OPT__FLAG_PAR_TARGET == FLAG_PAR_CAN ) ? -ParFlag_MinLv
                                                                                 : +ParFlag_MinLv;

//    particle 1: refine to the maximum level
      ParIntData_AllRank[PAR_FLAG][1] = ( OPT__FLAG_PAR_TARGET == FLAG_PAR_CAN ) ? -ParFlag_MaxLv
                                                                                 : +ParFlag_MaxLv;

//    particle 2: no refinement
      ParIntData_AllRank[PAR_FLAG][2] = PFLAG_NO;

//    particles 3/4: place around particle 1 with PAR_FLAG < 0 to allow refinement
      if ( OPT__FLAG_PAR_TARGET == FLAG_PAR_BOTH ) {
      ParFltData_AllRank[PAR_POSX][3] = ParFltData_AllRank[PAR_POSX][1] - amr->dh[0];
      ParFltData_AllRank[PAR_POSX][4] = ParFltData_AllRank[PAR_POSX][1] + amr->dh[0];
      ParFltData_AllRank[PAR_POSY][3] = ParFltData_AllRank[PAR_POSY][1];
      ParFltData_AllRank[PAR_POSY][4] = ParFltData_AllRank[PAR_POSY][1];

      ParIntData_AllRank[PAR_FLAG][3] = -ParFlag_MaxLv;
      ParIntData_AllRank[PAR_FLAG][4] = -ParFlag_MaxLv;
      }
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, _PAR_FLAG,
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
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] ParFltData_AllRank[v];
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] ParIntData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_ParFlag



#endif // #ifdef PARTICLE
