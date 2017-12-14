#include "GAMER.h"

#ifdef PARTICLE

extern int  TwoParAcc_Mode;
extern real TwoParAcc_MinR;
extern real TwoParAcc_MaxR;
extern int  TwoParAcc_LogSample;




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

   for (long p=0; p<amr->Par->NPar_Active_AllRank; p++)  amr->Par->Time[p] = Time[0];


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

   const real Cen[3] = { 0.5*amr->BoxSize[0],
                         0.5*amr->BoxSize[1],
                         0.5*amr->BoxSize[2] };

   real PosMin, PosMax, r, theta, phi, sign;


// ============================================================================================================
// set particle mass = 1.0 except for mode 3 (for which the mass of particle 1 is set to zero for the self-force test)
   switch ( TwoParAcc_Mode ) 
   {
      case 1:  Mass[0] = 1.0;   Mass[1] = 1.0;   break;
      case 2:  Mass[0] = 1.0;   Mass[1] = 1.0;   break;
      case 3:  Mass[0] = 1.0;   Mass[1] = 0.0;   break;

      default: Aux_Error( ERROR_INFO, "unsupported mode (%d) !!\n", TwoParAcc_Mode );
   }


// velocity must be initialized as zero since we use the updated velocity array to store the acceleration
   for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)
   for (int d=0; d<3; d++)
      Vel[d][p] = 0.0;


// set the particle position
   switch ( TwoParAcc_Mode ) 
   {
      case 1:
//       put the first particle randomly in the eight base-level patches around the box center
         for (int d=0; d<3; d++)
         {
            PosMin = Cen[d] - PS1*amr->dh[0];
            PosMax = Cen[d] + PS1*amr->dh[0];
            
            Pos[d][0] = ( (real)rand()/RAND_MAX )*( PosMax - PosMin ) + PosMin;
         }

//       put the second particle randomly with a distance MinR <= R <= MaxR from the first particle
         if ( TwoParAcc_LogSample )
         r     = exp(  ( (real)rand()/RAND_MAX )*( log(TwoParAcc_MaxR) - log(TwoParAcc_MinR) ) + log(TwoParAcc_MinR)  );
         else
         r     = ( (real)rand()/RAND_MAX )*( TwoParAcc_MaxR - TwoParAcc_MinR ) + TwoParAcc_MinR;

         theta = ( (real)rand()/RAND_MAX )*( M_PI           - 0.0            ) + 0.0;
         phi   = ( (real)rand()/RAND_MAX )*( 2.0*M_PI       - 0.0            ) + 0.0;

         Pos[0][1] = Pos[0][0] + r*sin(theta)*cos(phi);
         Pos[1][1] = Pos[1][0] + r*sin(theta)*sin(phi);
         Pos[2][1] = Pos[2][0] + r*cos(theta);
      break;


      case 2: case 3:
//       put particle 1 at the center of a MAX_LEVEL patch near the box center
         for (int d=0; d<3; d++)    Pos[d][1] = Cen[d] + 0.5*PS1*amr->dh[MAX_LEVEL];

//       displace particle 0 along the x direction
         if ( TwoParAcc_LogSample )
         r         = exp(  ( (real)rand()/RAND_MAX )*( log(TwoParAcc_MaxR) - log(TwoParAcc_MinR) ) + log(TwoParAcc_MinR)  );
         else
         r         = ( (real)rand()/RAND_MAX )*( TwoParAcc_MaxR - TwoParAcc_MinR ) + TwoParAcc_MinR;

         sign      = ( rand()%2 )*2 - 1;
         Pos[0][0] = Pos[0][1] + sign*r;

         for (int d=1; d<3; d++)    Pos[d][0] = Pos[d][1];
      break;

      default: Aux_Error( ERROR_INFO, "unsupported mode (%d) !!\n", TwoParAcc_Mode );
   } // switch ( TwoParAcc_Mode )
// ============================================================================================================

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
