#include "GAMER.h"

#if ( !defined GPU  &&  MODEL == HYDRO )




//-----------------------------------------------------------------------------------------
// Function    :  CPU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) required for the hydro solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
//
// Parameter   :  dt_Array  : Array to store the minimum dt in each target patch
//                Flu_Array : Array storing the prepared fluid data of each target patch
//                NPG       : Number of target patch groups
//                dh        : Grid size
//                Safety    : dt safety factor
//                Gamma     : Ratio of specific heats
//                MinPres   : Minimum allowed pressure
//
// Return      :  dt_Array
//-----------------------------------------------------------------------------------------
void CPU_dtSolver_HydroCFL( real dt_Array[], const real Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                            const int NPG, const real dh, const real Safety, const real Gamma, const real MinPres )
{

   const bool CheckMinPres_Yes = true;
   const int  NPatch           = 8*NPG;
   const real Gamma_m1         = Gamma - (real)1.0;

   real fluid[NCOMP_FLUID], _Rho, Vx, Vy, Vz, Pres, Cs, MaxV, MaxCFL;

// loop over all patches
#  pragma omp parallel for private( fluid, _Rho, Vx, Vy, Vz, Pres, Cs, MaxV, MaxCFL ) schedule( runtime )
   for (int p=0; p<NPatch; p++)
   {
      MaxCFL = (real)0.0;

      for (int t=0; t<CUBE(PS1); t++)
      {
         for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = Flu_Array[p][v][t];

        _Rho  = (real)1.0 / fluid[DENS];
         Vx   = FABS( fluid[MOMX] )*_Rho;
         Vy   = FABS( fluid[MOMY] )*_Rho;
         Vz   = FABS( fluid[MOMZ] )*_Rho;
         Pres = CPU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                 Gamma_m1, CheckMinPres_Yes, MinPres );
         Cs   = SQRT( Gamma*Pres*_Rho );

#        if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU  ||  FLU_SCHEME == WAF )
         MaxV   = FMAX( Vx, Vy );
         MaxV   = FMAX( Vz, MaxV );
         MaxCFL = FMAX( MaxV+Cs, MaxCFL );

#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
         MaxV   = Vx + Vy + Vz;
         MaxCFL = FMAX( MaxV+(real)3.0*Cs, MaxCFL );
#        endif
      } // for (int t=0; t<CUBE(PS1); t++)

      dt_Array[p] = Safety*dh/MaxCFL;

   } // for (int p=0; p<NPatch; p++)

} // FUNCTION : CPU_dtSolver_HydroCFL



#endif // #if ( !defined GPU  &&  MODEL == HYDRO )
