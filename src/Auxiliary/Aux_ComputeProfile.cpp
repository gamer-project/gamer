#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ComputeProfile
// Description :  Compute the average radial profile of a target field
//
// Note        :  1. Maximum radius adopted when actually computing the profile may be larger than the
//                   input "r_max"
//                   --> Because "r_max" in general does not coincide with the right edge of the maximum bin
//
// Parameter   :
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_ComputeProfile( Profile_t *Prof, const double Center[], const double r_max_input, const double dr_min,
                         const bool LogBin, const double LogBinRatio, const bool RemoveEmpty )
{

// check
#  ifdef GAMER_DEBUG
   if ( r_max_input <= 0.0 )
      Aux_Error( ERROR_INFO, "r_max_input (%14.7e) <= 0.0 !!\n", r_max_input );

   if ( dr_min <= 0.0 )
      Aux_Error( ERROR_INFO, "dr_min (%14.7e) <= 0.0 !!\n", dr_min );

   if ( LobBin  &&  LogBinRatio <= 1.0 )
      Aux_Error( ERROR_INFO, "LogBinRatio (%14.7e) <= 1.0 !!\n", LogBinRatio );
#  endif


// get the total number of radial bins and the corresponding maximum radius
   if ( LogBin )
   {
      Prof->NBin      = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 2;
      Prof->MaxRadius = dr_min*pow( LogBinRatio, Prof->NBin-1 );
   }

   else // linear bin
   {
      Prof->NBin      = (int)ceil( r_max_input / dr_min );
      Prof->MaxRadius = dr_min*Prof->NBin;
   }


// record profile parameters
   for (int d=0; d<3; d++)    Prof->Center[d] = Center[d];

   Prof->LogBin = LogBin;

   if ( LogBin )  Prof->LogBinRatio = LogBinRatio;


// allocate and initialize memory
   Prof->AllocateMemory();


// record radial coordinates
   if ( LogBin )
      for (int b=0; b<Prof->NBin; b++)    Prof->Radius[b] = dr_min*pow( LogBinRatio, b-0.5 );
   else
      for (int b=0; b<Prof->NBin; b++)    Prof->Radius[b] = (b+0.5)*dr_min;


// collect profile dat in this rank
   const double r_max2 = SQR( Prof->MaxRadius );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

         for (int k=0; k<PS1; k++)  {  const double dz = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  const double dy = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  const double dx = x0 + i*dh;

            const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

            if ( r2 < r_max2 )
            {
               const double r   = sqrt( r2 );
               const int    bin = ( LogBin ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(LogBinRatio) ) + 1  )
                                             : int( r/dr_min );
//             prevent from round-off errors
               if ( bin >= Prof->NBin )   continue;

//             check
#              ifdef GAMER_DEBUG
               if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#              endif

//###WORK      replace it with a user-specified function
               Prof->Data  [bin] += amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i]*dv;
               Prof->Weight[bin] += dv;
               Prof->NCell [bin] ++;
            } // if ( r2 < r_max2 )
         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


//###WORK
// in-place MPI reduce for the profile and weight


// compute profile by the root rank
   if ( MPI_Rank == 0 )
   {
      for (int b=0; b<Prof->NBin; b++)
         Prof->Data[b] /= Prof->Weight[b];

//    reset or remove the empty bins
      for (int b=0; b<Prof->NBin; b++)
      {
         if ( Prof->NCell[b] != 0L )   continue;

         if ( RemoveEmpty )
         {
            for (int b_up=b+1; b_up<Prof->NBin; b_up++)
            {
               const int b_up_m1 = b_up - 1;

               Prof->Radius[b_up_m1] = Prof->Radius[b_up];
               Prof->Data  [b_up_m1] = Prof->Data  [b_up];
               Prof->Weight[b_up_m1] = Prof->Weight[b_up];
               Prof->NCell [b_up_m1] = Prof->NCell [b_up];
            }

//          reset the total number of bins
            Prof->NBin --;

//          reset the maximum radius if we are removing the last bin
            if ( b == Prof->NBin )
            {
               if ( LogBin )
                  Prof->MaxRadius = dr_min*pow( LogBinRatio, Prof->NBin-1 );
               else
                  Prof->MaxRadius = dr_min*Prof->NBin;
            }

//          reduce counter since all bins above b have been shifted downed
            b --;
         } // if ( RemoveEmpty )

         else
         {
            Prof->Data[b] = 0.0;
         } // if ( RemoveEmpty ) ... else ...
      } // for (int b=0; b<Prof->NBin; b++)
   } // if ( MPI_Rank == 0 )


//###WORK
// broadcast Data and Weight

} // FUNCTION : Aux_ComputeProfile
