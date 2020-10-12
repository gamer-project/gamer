#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  CPU_ExtPotSolver_BaseLevel
// Description :  Add external potential on the base level
//
// Note        :  1. External potential is specified by the input function Func()
//                2. Set PotIsInit to false if the base-level potential has not been initialized
//                   --> Useful when self-gravity is disabled
//                3. Invoked by Gra_AdvanceDt()
//
// Parameter   :  Func             : Function pointer to the external potential routine
//                AuxArray_Flt/Int : Auxiliary floating-point/integer arrays for adding external potential
//                Table            : 3D potential table for EXT_POT_TABLE
//                Time             : Target physical time
//                PotIsInit        : Whether patch->pot[] has been initialized
//                                   --> true : **add** to the original data
//                                       false: **overwrite** the original data
//                SaveSg           : Sandglass to store the updated potential
//
// Return      :  amr->patch->pot[]
//-----------------------------------------------------------------------------------------
void CPU_ExtPotSolver_BaseLevel( const ExtPot_t Func, const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real Table[], const double Time, const bool PotIsInit, const int SaveSg )
{

// check
#  ifdef GAMER_DEBUG
   if ( Func == NULL )
      Aux_Error( ERROR_INFO, "Func == NULL !!\n" );

   if ( AuxArray_Flt == NULL )
      Aux_Error( ERROR_INFO, "AuxArray_Flt == NULL !!\n" );

   if ( AuxArray_Int == NULL )
      Aux_Error( ERROR_INFO, "AuxArray_Int == NULL !!\n" );

   if ( SaveSg != 0  &&  SaveSg != 1 )
      Aux_Error( ERROR_INFO, "incorrect SaveSg (%d) !!\n", SaveSg );

   if ( Time < 0.0 )
      Aux_Error( ERROR_INFO, "Time (%14.7e) < 0.0 !!\n", Time );
#  endif


   const int    lv   = 0;
   const double dh   = amr->dh[lv];
   const double dh_2 = 0.5*dh;

#  pragma omp parallel
   {
//    thread-private variables
      double x0, y0, z0, x, y, z;
      real   ExtPot;

#     pragma omp for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         x0 = amr->patch[0][lv][PID]->EdgeL[0] + dh_2;
         y0 = amr->patch[0][lv][PID]->EdgeL[1] + dh_2;
         z0 = amr->patch[0][lv][PID]->EdgeL[2] + dh_2;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

            ExtPot = Func( x, y, z, Time, AuxArray_Flt, AuxArray_Int, EXT_POT_USAGE_ADD, Table );

            if ( PotIsInit )  amr->patch[SaveSg][lv][PID]->pot[k][j][i] += ExtPot;  // add
            else              amr->patch[SaveSg][lv][PID]->pot[k][j][i]  = ExtPot;  // overwrite
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // end of OpenMP parallel region

} // FUNCTION : CPU_ExtPotSolver_BaseLevel



#endif // #ifdef GRAVITY
