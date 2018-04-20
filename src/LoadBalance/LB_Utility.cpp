#include "GAMER.h"


void  LB_Hilbert_i2c( ulong index, ulong coord[], const uint NBits );
ulong LB_Hilbert_c2i( ulong const coord[], const uint NBits );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Corner2Index
// Description :  Transform the 3D patch corner coordinatest to the 1D load-balance indices
//
// Note        :  1. "amr->ResPower2[lv]" must be prepared in advance
//                2. Periodic B.C. is applied
//                3. By inputting "amr->ResPower2[lv]" for LB_Hilbert_c2i, we ensure the following property
//                   "FaLBIdx*8 == SonLBIdx - SonLBIdx%8"
//                   --> Hilbert curve in different levels have the similar indexing order in space
//                   --> Son patches are more likely to be put in the same MPI rank as their father patches
//                4. This function is defined even when LOAD_BALANCE is off since we want to store LB_Idx for
//                   all patches in any case
//                5. Experiments show that "LB_Hilbert_c2i( Coord, NBits1 )" and ""LB_Hilbert_c2i( Coord, NBits2 )"
//                   return the same value if NBits1%3 = NBits2%3
//                   --> LB_Hilbert_c2i( Coord, NBits1 ) = LB_Hilbert_c2i( Coord, NBits1+3 ) = LB_Hilbert_c2i( Coord, NBits1+6 ) ...
//
// Parameter   :  Check :  Check whether the input corner lies in the simulation box
//                         --> effective only in the DEBUG mode
//
// Return      :  1D load-balance indices
//-------------------------------------------------------------------------------------------------------
long LB_Corner2Index( const int lv, const int Corner[], const Check_t Check )
{

#  ifdef GAMER_DEBUG
   if ( Check == CHECK_ON )
   for (int d=0; d<3; d++)
   {
      if ( Corner[d] < 0  ||  Corner[d] >= amr->BoxScale[d] )
         Aux_Error( ERROR_INFO, "Corner = (%8d,%8d,%8d) lies outside the simulation box !!\n",
                    Corner[0], Corner[1], Corner[2] );
   }
#  endif


   const int PatchScale  = amr->scale[lv]*PATCH_SIZE;
   const int PatchPower2 = (int)log2( PS1 );

#  ifdef GAMER_DEBUG
   if ( 1<<PatchPower2 != PS1 )  Aux_Error( ERROR_INFO, "2^%d != %d !!\n", PatchPower2, PS1 );
#  endif

   ulong Coord[3];
   int   Cr_Periodic[3];

   for (int d=0; d<3; d++)
   {
//###PERIODIC B.C. (ok if the non-periodic B.C. is adopted since other parts of the code have been carefully revised)
      Cr_Periodic[d] = ( Corner[d] + amr->BoxScale[d] ) % amr->BoxScale[d];
      Coord      [d] = Cr_Periodic[d] / PatchScale;
   }

   return LB_Hilbert_c2i( Coord, amr->ResPower2[lv]-PatchPower2 );

} // FUNCTION : LB_Corner2Index



#ifdef LOAD_BALANCE

//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Index2Corner
// Description :  Transform the 1D load-balance indices to the 3D patch corner coordinates
//
// Note        :  1. "amr->ResPower2[lv]" must be prepared in advance
//                2. By inputting "amr->ResPower2[lv]" for LB_Hilbert_i2c, we ensure the following property
//                   "FaLBIdx*8 == SonLBIdx - SonLBIdx%8"
//                   --> Hilbert curve in different levels have the similar indexing order in space
//                   --> Son patches are more likely to be put in the same MPI rank as their father patches
//
// Parameter   :  Check :  Check whether the output corner lies in the simulation box
//                         --> effective only in the DEBUG mode
//-------------------------------------------------------------------------------------------------------
void LB_Index2Corner( const int lv, const long LB_Idx, int Corner[], const Check_t Check )
{

   const int PatchScale = amr->scale[lv]*PATCH_SIZE;
   const int PatchPower2 = (int)log2( PS1 );
   ulong Coord[3];

#  ifdef GAMER_DEBUG
   if ( 1<<PatchPower2 != PS1 )  Aux_Error( ERROR_INFO, "2^%d != %d !!\n", PatchPower2, PS1 );
#  endif

   LB_Hilbert_i2c( LB_Idx, Coord, amr->ResPower2[lv]-PatchPower2 );

   for (int d=0; d<3; d++)    Corner[d] = Coord[d]*PatchScale;

#  ifdef GAMER_DEBUG
   if ( Check == CHECK_ON )
   for (int d=0; d<3; d++)
   {
      if ( Corner[d] < 0  ||  Corner[d] >= NX0_TOT[d]*amr->scale[0] )
         Aux_Error( ERROR_INFO, "Corner = (%8d,%8d,%8d) lies outside the simulation box !!\n",
                    Corner[0], Corner[1], Corner[2] );
   }
#  endif

} // FUNCTION : LB_Index2Corner



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Index2Rank
// Description :  Return the MPI rank which the input LB_Idx belongs to
//
// Note        :  "LB_CutPoint[lv]" must be prepared in advance
//
// Return      :  success  : target MPI rank
//                fail     : -1
//
// Parameter   :  lv       : Refinement level of the input LB_Idx
//                LB_Idx   : Space-filling-curve index for load balance
//                Check    : Check whether or not the input LB_Idx lies in the range of LB_CutPoint
//                           (this check may fail during the grid refinement --> if that's the case,
//                            turn this check off during the grid refinement)
//
// Return      :  MPI rank
//-------------------------------------------------------------------------------------------------------
int LB_Index2Rank( const int lv, const long LB_Idx, const Check_t Check )
{

//###OPTIMIZATION: better search algorithm
   for (int r=0; r<MPI_NRank; r++)
      if ( LB_Idx >= amr->LB->CutPoint[lv][r]  &&  LB_Idx < amr->LB->CutPoint[lv][r+1] )  return r;

   if ( Check == CHECK_ON )
      Aux_Error( ERROR_INFO, "no target rank was found for lv %d, LB_Idx %ld !!\n",
                 lv, LB_Idx );

   return -1;

} // FUNCTION : LB_Index2Rank

#endif // #ifdef LOAD_BALANCE
