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
//                   inside the Grackle solver, subcycling time-steps are adopted to ensure correct evolution
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
   const double dt_phy = DT__GRACKLE_COOLING * MinCoolingTime;

#  ifdef COMOVING
// convert the time-step size to comoving coordinates
   const double dt = dt_phy / SQR(Time[lv]);
#  else
   const double dt = dt_phy;
#  endif

   return dt;

} // FUNCTION : Grackle_GetTimeStep_CoolingTime



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMinCoolingTime
// Description :  Get the minimum Grackle cooling time at the target level among all MPI ranks
//
// Note        :  1. Invoked by Grackle_GetTimeStep_CoolingTime()
//                2. OpenMP parallelization should be performed within Grackle, not outside of Grackle_Calculate()
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MinCoolingTime_AllRank
//-------------------------------------------------------------------------------------------------------
real GetMinCoolingTime( const int lv )
{

   real  MinCoolingTime = INFINITY;
   bool  AnyCell        = false;

   const int NTotal  = amr->NPatchComma[lv][1]/8;          // total number of patch groups to calculate the cooling time
   const int NPG_Max = CHE_GPU_NPGROUP;                    // maximum number of patch groups to calculate the cooling time at a time

   int  *PID0_List = new int [NTotal];                     // list recording the patch indices with LocalID==0 to calculate the cooling time
   for (int t=0; t<NTotal; t++)  PID0_List[t] = 8*t;

   real *Grackle_TCool = new real [ CUBE(PS2)*NPG_Max ];   // array storing NPG_Max patch groups of grackle cooling time

// get the minimum Grackle cooling time in this rank
   for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)
   {
      const int NPG = ( NPG_Max < NTotal-Disp ) ? NPG_Max : NTotal-Disp; // number of patch groups to calculate the cooling time at a time

//    calculate the Grackle cooling time
      Grackle_Calculate( Grackle_TCool, _GRACKLE_TCOOL, lv, NPG, PID0_List+Disp );

//    loop over NPG patch groups to find the minimum with OpenMP parallelization
#     pragma omp parallel for reduction( min:MinCoolingTime ) reduction( ||:AnyCell ) schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         const int PID0 = (PID0_List+Disp)[TID];

//       loop over 8 patches in the patch group
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int PID = PID0  + LocalID;
            const int N   = TID*8 + LocalID;

//          if OPT__FIXUP_RESTRICT is enabled, skip all non-leaf patches
//          because they are later overwritten by the refined patches
//          --> note that this leads to the timestep being "inf" when a level is completely refined
            if ( OPT__FIXUP_RESTRICT  &&  amr->patch[0][lv][PID]->son != -1 )    continue;

            AnyCell = true;

//          loop over cells in the patch
            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
//             remember to take the absolute value for the cooling time, which could be negative or positive
               const real AbsTCool = FABS( Grackle_TCool[ (N*CUBE(PS1) + k*SQR(PS1) + j*PS1 + i) ] );

//             take the minimum
//             --> it is fine if Grackle returns an abnormally NaN cooling time; FMIN will not return it
               MinCoolingTime = FMIN( MinCoolingTime, AbsTCool );

            } // k,j,i
         } // for (int LocalID=0; LocalID<8; LocalID++)
      } // for (int TID=0; TID<NPG; TID++)
   } // for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)

   delete [] PID0_List;
   delete [] Grackle_TCool;


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
