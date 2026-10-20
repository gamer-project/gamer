#include "GAMER.h"

#ifdef EXACT_COOLING_GENERAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_ExactCooling_General
// Description :  Estimate the evolution time-step constrained by the general exact-cooling source term
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Invoked by Mis_GetTimeStep()
//                3. Enabled by setting the runtime option SRC_ECG_DTCOEF >= 0.0
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep_ExactCooling_General( const int lv, const double dTime_dt )
{

   if ( !SrcTerms.ExactCooling_General )   return HUGE_NUMBER;


   const double dh = amr->dh[lv];

// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double  dt_ECG     = HUGE_NUMBER;
   double *OMP_dt_ECG = new double [NT];


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize the array
      OMP_dt_ECG[TID] = HUGE_NUMBER;

#     pragma omp for schedule( static )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         for (int k=0; k<PS1; k++)  {
         for (int j=0; j<PS1; j++)  {
         for (int i=0; i<PS1; i++)  {

            real tcool_Code = NULL_REAL;

            if ( SrcTerms.ECG_TCoolGInit[lv] ) {
//             use the stored cooling time
#              ifdef TCOOLG
               tcool_Code = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[TCOOLG][k][j][i];
#              endif
            }

            else
            {
#              ifdef MHD
               real B[NCOMP_MAG];
               MHD_GetCellCenteredBFieldInPatch( B, lv, PID, i, j, k, amr->MagSg[lv] );
#              else
               real *B = NULL;
#              endif

//             call Src_ExactCooling() to compute the cooling time if not initialized yet
               const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh;
               const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
               const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;

//             get the input arrays
               real fluid[FLU_NIN_S];
               for (int v=0; v<FLU_NIN_S; v++)  fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

//             call the CPU version of the exact-cooling function
//             --> must set dt=0.0 since we only want to compute the cooling time fluid[TCOOL] without evolving the gas
//             --> must NOT pass the EoS object since it's only used by the GPU code
               SrcTerms.ECG_CPUPtr( fluid, B, &SrcTerms, 0.0, NULL_REAL, x, y, z, NULL_REAL, NULL_REAL,
                                    MIN_DENS, MIN_PRES, MIN_EINT, PassiveFloorMask, NULL,
                                    Src_ECG_AuxArray_Flt, Src_ECG_AuxArray_Int );
#              ifdef TCOOLG
               tcool_Code = fluid[TCOOLG];
#              endif
            } // if ( SrcTerms.ECG_TCoolGInit[lv] ) ... else ...

//          compare the cooling time and store the minimum value
            OMP_dt_ECG[TID] = fmin( OMP_dt_ECG[TID], (double)tcool_Code );

         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // OpenMP parallel region

// find the minimum over all OpenMP threads
   for (int TID=0; TID<NT; TID++)   dt_ECG = fmin( dt_ECG, OMP_dt_ECG[TID] );

// free per-thread arrays
   delete [] OMP_dt_ECG;

// find the minimum over all MPI processes
#  ifndef SERIAL
   MPI_Allreduce( MPI_IN_PLACE, &dt_ECG, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
#  endif


   return SrcTerms.ECG_dtCoef * dt_ECG;

} // FUNCTION : Mis_GetTimeStep_ExactCooling_General



#endif // #ifdef EXACT_COOLING_GENERAL
