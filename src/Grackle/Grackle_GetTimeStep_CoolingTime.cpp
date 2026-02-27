#include "GAMER.h"

#ifdef SUPPORT_GRACKLE

static real GetMinCoolingTime( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_GetTimeStep_CoolingTime
// Description :  Estimate the evolution time-step from Grackle cooling time
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Time-step is set to restrict the cooling process
//                   --> dt = DT__GRACKLE_COOLING * Min(CoolingTime)
//                3. This is usually for the purpose of outputting data for analysis;
//                   Inside the Grackle solver, subcycling time-steps will be adopted for the correctness of evolution
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Grackle_GetTimeStep_CoolingTime( const int lv )
{

// get the minimum cooling time
   const real   MinCoolingTime = GetMinCoolingTime( lv );

// get the time-step
   const double dt = DT__GRACKLE_COOLING * MinCoolingTime;

   return dt;

} // FUNCTION : Grackle_GetTimeStep_CoolingTime



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMinCoolingTime
// Description :  Get the minimum Grackle cooling time at the target level among all MPI ranks
//
// Note        :  1. Invoked by Grackle_GetTimeStep_CoolingTime()
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MinCoolingTime_AllRank
//-------------------------------------------------------------------------------------------------------
real GetMinCoolingTime( const int lv )
{

   real   MinCoolingTime = INFINITY;
   bool   AnyCell        = false;

#  pragma omp parallel
   {

      real *Grackle_TCool = new real [ CUBE(PS2) ];   // array storing ONE patch group of grackle cooling time

//    get the minimum Grckle cooling time in this rank
#     pragma omp for reduction( min:MinCoolingTime ) reduction( ||:AnyCell ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {

//       calculate the Grackle cooling time
         Grackle_Calculate( Grackle_TCool, _GRACKLE_TCOOL, lv, 1, &PID0 );

         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int PID = PID0 + LocalID;

//          if OPT__FIXUP_RESTRICT is enabled, skip all non-leaf patches
//          becausethey are later overwritten by the refined patches
//          note that this leads to the timestep being "inf" when a level is completely refined
            if ( OPT__FIXUP_RESTRICT  &&  amr->patch[0][lv][PID]->son != -1 )    continue;

            AnyCell = true;

//          calculate the minimum
            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
               // remember to take the absolute value for the cooling time, which could be negative or positive
               MinCoolingTime = MIN( MinCoolingTime,
                                     FABS( Grackle_TCool[ (LocalID*CUBE(PS1) + k*SQR(PS1) + j*PS1 + i) ] ) );
            } // k,j,i
         } // for (int LocalID=0; LocalID<8; LocalID++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

      delete [] Grackle_TCool;

   } // OpenMP parallel region


// get the minimum Grackle cooling time in all ranks
   real MinCoolingTime_AllRank;
   bool AnyCell_AllRank;

   MPI_Allreduce( &MinCoolingTime, &MinCoolingTime_AllRank, 1, MPI_GAMER_REAL, MPI_MIN,    MPI_COMM_WORLD );
   MPI_Reduce(    &AnyCell,        &AnyCell_AllRank,        1, MPI_CXX_BOOL,   MPI_LOR, 0, MPI_COMM_WORLD );


// check
   if ( MPI_Rank == 0  &&  MinCoolingTime_AllRank == INFINITY  &&  AnyCell_AllRank )
      Aux_Error( ERROR_INFO, "MinCoolingTime == INFINITY at lv %d !!\n", lv );


   return MinCoolingTime_AllRank;

} // FUNCTION : GetMinCoolingTime



#endif // #ifdef SUPPORT_GRACKLE
