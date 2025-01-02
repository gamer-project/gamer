#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined FEEDBACK )

extern bool   Plummer_FB_Exp;
extern double Plummer_FB_ExpEMin;
extern double Plummer_FB_ExpEMax;
extern double Plummer_FB_ExpMMin;
extern double Plummer_FB_ExpMMax;
extern bool   Plummer_FB_Acc;
extern double Plummer_FB_AccMMin;
extern double Plummer_FB_AccMMax;
extern double Plummer_FB_Like;




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Plummer
// Description :  Example of explosion and mass accretion feedbacks
//                --> Just an example and is by no means to be physically correct
//
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAttFlt/Int[], respectively
//                   --> This function is responsible for updating gas and particles within
//                       ** FB_GHOST_SIZE <= cell indices i,j,k < FB_GHOST_SIZE+PS2 **
//                   --> Updating gas and particles outside this range is fine but will have no effect at all
//                2. Must use ParSortID[] to access ParAttFlt[] and ParAttInt[]
//                   --> ParAttFlt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
//                   --> ParAttInt[PAR_TYPE/etc][ ParSortID[...] ]
//                3. Particles may be outside the target region
//                4. To ensure the consistency of random numbers, one must call the random number generator for
//                   ALL particles, including those too far away to affect the target region
//                5. No need to worry about the periodic boundary condition here
//                   --> Particle positions have been remapped in FB_AdvanceDt()
//                6. CoarseFine[] records the coarse-fine boundaries along the 26 sibling directions, defined as
//                            24  13  25
//                            15  05  17     z+1 plane
//                            22  12  23
//
//                            08  03  09
//                            00  XX  01     z   plane
//                            06  02  07
//                   y
//                   ^        20  11  21
//                   |        14  04  16     z-1 plane
//                   --->x    18  10  19
//                7. Invoked by FB_AdvanceDt()
//                8. Must NOT change particle positions
//                9. Since Fluid[] stores both the input and output data, the order of particles may affect the
//                   final output results
//                   --> For example, particle 2 may use the data updated by particle 1 as the input data
//                   --> Actually, even if we separate Fluid[] to input and output arrays, the final output results
//                       may still depend on the order of particles for non-local feedback since different particles
//                       may update the same cell
//                10. In general, it is recommended to have the maximum feedback radius no larger than half of the patch size
//                    (i.e., PATCH_SIZE/2=4 cells for PATCH_SIZE=8)
//                    --> Increase PATCH_SIZE if necessary
//                11. Linked to FB_User_Ptr in FB_Init_Plummer()
//                12. Must have FB_GHOST_SIZE>=1 for the mass accretion feedback
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                NPar       : Number of particles
//                ParSortID  : Sorted particle IDs
//                ParAttFlt  : Particle floating-point attribute arrays
//                ParAttInt  : Particle integer        attribute arrays
//                Fluid      : Array to store the input/output fluid data
//                             --> Array size is fixed to (FB_NXT)^3=(PS2+2*FB_GHOST_SIZE)^3
//                EdgeL      : Left edge of Fluid[]
//                             --> Right edge is given by EdgeL[]+FB_NXT*dh
//                dh         : Cell size of Fluid[]
//                CoarseFine : Coarse-fine boundaries along the 26 sibling directions
//                TID        : Thread ID
//                RNG        : Random number generator
//                             --> Random number can be obtained by "RNG->GetValue( TID, Min, Max )",
//                                 where Min/Max specify the range of random numbers
//
// Return      :  Fluid, ParAttFlt, ParAttInt
//-------------------------------------------------------------------------------------------------------
int FB_Plummer( const int lv, const double TimeNew, const double TimeOld, const double dt,
                const int NPar, const long *ParSortID, real_par *ParAttFlt[PAR_NATT_FLT_TOTAL], long_par *ParAttInt[PAR_NATT_INT_TOTAL],
                real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[],
                const int TID, RandomNumber_t *RNG )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid == NULL )    Aux_Error( ERROR_INFO, "Fluid == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAttFlt == NULL )   Aux_Error( ERROR_INFO, "ParAttFlt == NULL for NPar = %d !!\n", NPar );
      if ( ParAttInt == NULL )   Aux_Error( ERROR_INFO, "ParAttInt == NULL for NPar = %d !!\n", NPar );
   }
#  endif // #ifdef GAMER_DEBUG


   const double _dh  = 1.0 / dh;
   const double dv   = CUBE( dh );
   const double _dv  = 1.0 / dv;
   const int    MaxR = 1;

   if ( Plummer_FB_Acc  &&  FB_GHOST_SIZE < MaxR )
      Aux_Error( ERROR_INFO, "FB_GHOST_SIZE (%d) < MaxR (%d) for Plummer_FB_Acc !!\n", FB_GHOST_SIZE, MaxR );

   bool CheckCF = false;
#  if ( FB_GHOST_SIZE > 0 )
   for (int s=0; s<26; s++) {
      if ( CoarseFine[s] ) {
         CheckCF = true;
         break;
      }
   }
#  endif

   for (int t=0; t<NPar; t++)
   {
      const long   p      = ParSortID[t];
      const double xyz[3] = { (double)ParAttFlt[PAR_POSX][p], (double)ParAttFlt[PAR_POSY][p], (double)ParAttFlt[PAR_POSZ][p] };

      int idx[3];
      for (int d=0; d<3; d++)    idx[d] = (int)floor( ( xyz[d] - EdgeL[d] )*_dh );


//    skip particles too close to coarse-fine boundaries
      bool Skip = false;

      if ( CheckCF )
      {
//       only need to check the 27 extreme cases
         int ijk[3];
         for (int dk=-1; dk<=1; dk++)  {  if ( Skip )  break;  ijk[2] = idx[2] + dk*MaxR;
         for (int dj=-1; dj<=1; dj++)  {  if ( Skip )  break;  ijk[1] = idx[1] + dj*MaxR;
         for (int di=-1; di<=1; di++)  {  if ( Skip )  break;  ijk[0] = idx[0] + di*MaxR;

            const int CellPatchRelPos = FB_Aux_CellPatchRelPos( ijk );
            if ( CellPatchRelPos != -1  &&  CoarseFine[CellPatchRelPos] )  Skip = true;   // cell is in a coarse patch

         }}}
      }


      if ( RNG->GetValue(TID,0.0,1.0) < Plummer_FB_Like )
      {
         const real ExpMassFac = 1.0 - RNG->GetValue( TID, Plummer_FB_ExpMMin, Plummer_FB_ExpMMax );
         const real ExpEngyFac = 1.0 + RNG->GetValue( TID, Plummer_FB_ExpEMin, Plummer_FB_ExpEMax );
         const real AccMassFac =       RNG->GetValue( TID, Plummer_FB_AccMMin, Plummer_FB_AccMMax );

//       to ensure the consistency of random numbers, it's better to skip particles **after** generating all random numbers
//       --> since the same particle may be skipped when updating one patch group but not skipped when updating another patch group
//           (e.g., if we skip particles in nearby patch groups too far away from a target patch group)
         if ( Skip )    continue;

         if ( OPT__VERBOSE ) {
            Aux_Message( stdout, "\n" );
            Aux_Message( stdout, "---------------------------------------------------\n" );
            Aux_Message( stdout, "Rank %d, Time %21.14e, ParID % 5d: ExpMassFac %13.7e, ExpEngyFac %13.7e, AccMassFac %13.7e\n",
                         MPI_Rank, TimeNew, p, ExpMassFac, ExpEngyFac, AccMassFac );
            Aux_Message( stdout, "---------------------------------------------------\n" );
         }

//       explosion
         if ( Plummer_FB_Exp )
         {
//          inject explosion energy to the nearby 27 cells including itself
//          --> skip cells in the ghost zones
            for (int dk=-1; dk<=1; dk++)  {  const int k = idx[2] + dk;  if ( k < FB_GHOST_SIZE || k >= FB_GHOST_SIZE+PS2 )   continue;
            for (int dj=-1; dj<=1; dj++)  {  const int j = idx[1] + dj;  if ( j < FB_GHOST_SIZE || j >= FB_GHOST_SIZE+PS2 )   continue;
            for (int di=-1; di<=1; di++)  {  const int i = idx[0] + di;  if ( i < FB_GHOST_SIZE || i >= FB_GHOST_SIZE+PS2 )   continue;

               Fluid[ENGY][k][j][i] *= ExpEngyFac;

            }}}

//          reduce particle mass
//          --> skip particles in the ghost zones
            if ( idx[0] >= FB_GHOST_SIZE  &&  idx[0] < FB_GHOST_SIZE+PS2  &&
                 idx[1] >= FB_GHOST_SIZE  &&  idx[1] < FB_GHOST_SIZE+PS2  &&
                 idx[2] >= FB_GHOST_SIZE  &&  idx[2] < FB_GHOST_SIZE+PS2   )
               ParAttFlt[PAR_MASS][p] *= (real_par)ExpMassFac;
         } // if ( Plummer_FB_Exp )


//       mass accretion
         if ( Plummer_FB_Acc )
         {
            real   dM;
            double dM_sum = 0.0;

//          accrete mass from the nearby 27 cells including itself
//          --> include cells in the ghost zones in order to
//              (i) ensure the following particles see the updated gas density
//                  (since particles in the ghost zones may affect particles inside the patch group)
//              (ii) get the correct total accreted mass
            for (int dk=-1; dk<=1; dk++)  {  const int k = idx[2] + dk;  if ( k < 0 || k >= FB_NXT )  continue;
            for (int dj=-1; dj<=1; dj++)  {  const int j = idx[1] + dj;  if ( j < 0 || j >= FB_NXT )  continue;
            for (int di=-1; di<=1; di++)  {  const int i = idx[0] + di;  if ( i < 0 || i >= FB_NXT )  continue;

               dM                    = AccMassFac*Fluid[DENS][k][j][i]*dv;
               dM_sum               += dM;
               Fluid[DENS][k][j][i] -= dM*_dv;

            }}}

//          increase particle mass
//          --> skip particles in the ghost zones
            if ( idx[0] >= FB_GHOST_SIZE  &&  idx[0] < FB_GHOST_SIZE+PS2  &&
                 idx[1] >= FB_GHOST_SIZE  &&  idx[1] < FB_GHOST_SIZE+PS2  &&
                 idx[2] >= FB_GHOST_SIZE  &&  idx[2] < FB_GHOST_SIZE+PS2   )
               ParAttFlt[PAR_MASS][p] += (real_par)dM_sum;
         } // if ( Plummer_FB_Acc )

      } // if ( RNG->GetValue(TID,0.0,1.0) < Plummer_FB_Like )
   } // for (int t=0; t<NPar; t++)


   return GAMER_SUCCESS;

} // FUNCTION : FB_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_Plummer
// Description :  Free the resources used by the user-specified feedback
//
// Note        :  1. Invoked by FB_End()
//                2. Linked to FB_End_User_Ptr in FB_Init_Plummer()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_Plummer()
{


} // FUNCTION : FB_End_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_Plummer
// Description :  Initialize the user-specified feedback
//
// Note        :  1. Invoked by FB_Init()
//                   --> Enable it by linking to the function pointer "FB_Init_User_Ptr"
//                2. Set FB_User_Ptr and FB_End_User_Ptr
//
// Parameter   :  None
//
// Return      :  FB_User_Ptr and FB_End_User_Ptr
//-------------------------------------------------------------------------------------------------------
void FB_Init_Plummer()
{

   FB_User_Ptr     = FB_Plummer;
   FB_End_User_Ptr = FB_End_Plummer;

} // FUNCTION : FB_Init_Plummer



#endif // #if ( MODEL == HYDRO  &&  defined FEEDBACK )
