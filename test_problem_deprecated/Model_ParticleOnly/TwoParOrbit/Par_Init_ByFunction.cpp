#include "GAMER.h"

#ifdef PARTICLE

extern real TwoParOrbit_M;
extern real TwoParOrbit_v;
extern real TwoParOrbit_R;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize the particle position and velocity 
//
// Note        :  Invoked by "Init_GAMER"
//
// Parameter   :  None
//
// Return      :  amr->Par->PosX/Y/Z, amr->Par->VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// synchronize all particles to the physical time at the base level
   for (long p=0; p<amr->Par->NPar_Active_AllRank; p++)  amr->Par->Time[p] = Time[0];


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

   const real Cen[3] = { 0.5*amr->BoxSize[0],
                         0.5*amr->BoxSize[1],
                         0.5*amr->BoxSize[2] };

   for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)
   {
//    set particle mass to extremely small values when using external potential so that self-gravity can be ignored
      Mass  [p] = ( OPT__EXTERNAL_POT ) ? TwoParOrbit_M*1.0e-20 : TwoParOrbit_M;

      Pos[0][p] = Cen[0] + (1.0-p*2.0)*TwoParOrbit_R;
      Pos[1][p] = Cen[1];
      Pos[2][p] = Cen[2];

      Vel[0][p] = 0.0;
      Vel[1][p] =          (1.0-p*2.0)*TwoParOrbit_v;
      Vel[2][p] = 0.0;
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
