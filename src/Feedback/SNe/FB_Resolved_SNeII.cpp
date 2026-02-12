#include "GAMER.h"

#ifdef FEEDBACK


static const int      maxfbDiameter  =   FB_GHOST_SIZE + 1;  // maximum diameter to apply feedback, constraint by ghost zone size
static       real  ***fbDepositWeighting[FB_GHOST_SIZE + 1]; // array of weighting for each feedback diameter

static const int      nVarRecSNeII   = 23;                   // number of variables to be record for each SNII
static const int      maxNumRecSNeII = 100;                  // maximum  number of recorded SNeII in each OpenMP thread
static       int     *numRecSNeII    = NULL;                 // number of recorded SNeII in each OpenMP thread
static       double **recordSNeII    = NULL;                 // array of stored variables to record the SNeII info




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Resolved_SNeII
// Description :  Resolved type II supernovae explosion feedback
//
// Note        :  1. Input fluid, output fluid and particle data are stored in Fluid_In[], Fluid_Out[] and ParAttFlt/Int[], respectively
//                   --> This function is responsible for updating gas and particles within
//                       ** FB_GHOST_SIZE <= cell indices i,j,k < FB_GHOST_SIZE+PS2 **
//                   --> Updating gas and particles outside this range is fine but will have no effect at all
//                2. Must use ParSortID[] to access ParAttFlt/Int[]
//                   --> ParAttFlt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
//                   --> ParAttInt[PAR_TYPE/etc][ ParSortID[...] ]
//                3. Particles may be outside the target region
//                4. To ensure the consistency of random numbers, one must call the random number generator for
//                   ALL particles, including those too far away to affect the target region
//                5. No need to worry about the periodic boundary condition here
//                   --> Particle positions have been remapped in FB_AdvanceDt()
//                6. Non-local feedback should not cross the boundary to either a coarser or a finer level
//                   to provide complete feedback and avoid inconsistency between levels
//                   --> Use CoarseFine[] and NonLeaf[] to check whether the cells are in the coarse patches and non-leaf patches, respectively.
//                   (a) CoarseFine[] records the coarse-fine boundaries along the 26 sibling directions, defined as
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
//                   (b) NonLeaf[] records the non-leaf patches among the 64 local patches, defined as
//                          62  46  47  63
//                          51  30  31  55
//                          50  28  29  54   z=3 plane
//                          60  44  45  61
//
//                          37  22  23  39
//                          11  05  07  15
//                          10  03  06  14   z=2 plane
//                          33  18  19  35
//
//                          36  20  21  38
//                          09  02  04  13
//                          08  00  01  12   z=1 plane
//                          32  16  17  34
//
//                   y      58  42  43  59
//                   ^      49  25  27  53
//                   |      48  24  26  52   z=0 plane
//                   --->x  56  40  41  57
//                7. Invoked by FB_AdvanceDt()
//                8. Must NOT change particle positions
//                9. Although the input and output data are separated as Fluid_In[] and Fluid_Out[],
//                   the final output results may still depend on the order of particles for non-local feedback
//                   since different particles may update the same cell
//                10. In general, it is recommended to have the maximum feedback radius no larger than half of the patch size
//                    (i.e., PATCH_SIZE/2=4 cells for PATCH_SIZE=8)
//                    --> Increase PATCH_SIZE if necessary
//                11. This feedback method assumes that the particle mass resolution is sufficiently high,
//                    so there will be at most one SNII per particle
//                    Whether a star particle has a SNII progenitor and when will it explode
//                    are already determined when the star formed
//                    --> Please add the attribute ParSNIITime in the star formation scheme accordingly
//                        (e.g. in StarFormation/SF_CreateStar_AGORA.cpp)
//                    --> After the explosion, its ParSNIITime will be set as negative so it will not explode again
//                12. This feedback method assumes the grid resolution is sufficiently high
//                    to resolve the Sedov phase blast wave caused by the supernova explosion
//                    --> Thermal energy (along with the mass and metal from SN explosion) is
//                        only injected into one cell where the particle is located
//                        and thre is no kinetic (outward momentum) feedback
//                13. Ref: Sec. 2.6 of Chia-Yu Hu, et al., 2023, ApJ, 950, 132 (https://doi.org/10.3847/1538-4357/accf9e)
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
//                Fluid_In   : Array to store the input fluid data
//                             --> Array size is fixed to (FB_NXT)^3=(PS2+2*FB_GHOST_SIZE)^3
//                Fluid_Out  : Array to store the output fluid data
//                             --> Array size is fixed to (FB_NXT)^3=(PS2+2*FB_GHOST_SIZE)^3
//                EdgeL      : Left edge of Fluid_In[] and Fluid_Out[]
//                             --> Right edge is given by EdgeL[]+FB_NXT*dh
//                dh         : Cell size of Fluid_In[] and Fluid_Out[]
//                CoarseFine : Coarse-fine boundaries along the 26 sibling directions
//                NonLeaf    : Non-leaf patches among the 64 patches of the Fluid array
//                TID        : Thread ID
//                RNG        : Random number generator
//                             --> Random number can be obtained by "RNG->GetValue( TID, Min, Max )",
//                                 where Min/Max specify the range of random numbers
//
// Return      :  Fluid_Out, ParAttFlt, ParAttInt
//-------------------------------------------------------------------------------------------------------
int FB_Resolved_SNeII( const int lv, const double TimeNew, const double TimeOld, const double dt,
                       const int NPar, const long *ParSortID, real_par *ParAttFlt[PAR_NATT_FLT_TOTAL], long_par *ParAttInt[PAR_NATT_INT_TOTAL],
                       const real (*Fluid_In)[FB_NXT][FB_NXT][FB_NXT],
                       real (*Fluid_Out)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[], bool NonLeaf[],
                       const int TID, RandomNumber_t *RNG )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid_In  == NULL )      Aux_Error( ERROR_INFO, "Fluid_In  == NULL !!\n" );
   if ( Fluid_Out == NULL )      Aux_Error( ERROR_INFO, "Fluid_Out == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAttFlt == NULL )   Aux_Error( ERROR_INFO, "ParAttFlt == NULL for NPar = %d !!\n", NPar );
      if ( ParAttInt == NULL )   Aux_Error( ERROR_INFO, "ParAttInt == NULL for NPar = %d !!\n", NPar );
   }
#  endif // #ifdef GAMER_DEBUG

// cell volume
   const real dv = CUBE( dh );

// determine if metallicity is included
   const bool UseMetal = ( Idx_Metal != Idx_Undefined  &&  Idx_ParMetalFrac != Idx_Undefined );

// whether to check the coarse-fine boundary
   bool CheckCF = false;
   if ( maxfbDiameter > 1 )
   for (int s=0; s<26; s++)
   {
      if ( CoarseFine[s] )
      {
         CheckCF = true;
         break;
      }
   } // for (int s=0; s<26; s++)

// whether to check the non-leaf patch
   bool CheckNL = false;
   if ( maxfbDiameter > 1 )
   for (int nbp=0; nbp<64; nbp++)
   {
      if ( NonLeaf[nbp] )
      {
         CheckNL = true;
         break;
      }
   } // for (int nbp=0; nbp<64; nbp++)

   for (int p=0; p<NPar; p++)
   {

//    1. prepare the input
      const long     par_idx       =   ParSortID[p];                                        // index of particle

      const real_par par_mass      =   ParAttFlt[        PAR_MASS][par_idx];                // mass of particle
      const real_par par_pos[3]    = { ParAttFlt[        PAR_POSX][par_idx],
                                       ParAttFlt[        PAR_POSY][par_idx],
                                       ParAttFlt[        PAR_POSZ][par_idx] };              // position of particle
      const real_par par_vel[3]    = { ParAttFlt[        PAR_VELX][par_idx],
                                       ParAttFlt[        PAR_VELY][par_idx],
                                       ParAttFlt[        PAR_VELZ][par_idx] };              // velocity of particle
      const long_par par_type      =   ParAttInt[        PAR_TYPE][par_idx];                // type of particle
#     ifdef STAR_FORMATION
      const real_par par_creTime   =   ParAttFlt[  Idx_ParCreTime][par_idx];                // creation time of particle
#     else
      const real_par par_creTime   =   0.0;
#     endif
      const real_par par_SNIITime  =   ParAttFlt[ Idx_ParSNIITime][par_idx];                // SNII explosion time of particle
      const real_par par_metalMass = ( UseMetal ) ?
                                       ParAttFlt[Idx_ParMetalFrac][par_idx]*par_mass : 0.0; // metal fraction of particle

      const int  cell_idx[3]       = { (int)floor( ( par_pos[0] - EdgeL[0] )/dh ),          // indices of cell where the particle is in
                                       (int)floor( ( par_pos[1] - EdgeL[1] )/dh ),          // Note it could be negative for particles in nearby patches
                                       (int)floor( ( par_pos[2] - EdgeL[2] )/dh ) };
      const bool isOnLeftCell[3]   = { ( par_pos[0] < EdgeL[0] + (cell_idx[0]+0.5)*dh ),    // whether the particle is on the left side of the cell
                                       ( par_pos[1] < EdgeL[1] + (cell_idx[1]+0.5)*dh ),    // true  = 1 -> on the left side
                                       ( par_pos[2] < EdgeL[2] + (cell_idx[2]+0.5)*dh ) };  // false = 0 -> on the right side

#     ifdef GAMER_DEBUG
//    sanity check for the particle
//    particle should not be in the coarse patch
      if ( CheckCF )
      {
         const int CellPatchRelPos = FB_Aux_CellPatchRelPos( cell_idx );
         if ( CellPatchRelPos != -1  &&  CoarseFine[CellPatchRelPos] )
            Aux_Error( ERROR_INFO, "The particle in the feedback routine should not be in the coarse patch !!\n" );
      }

//    particle should not be in the non-leaf patch
      if ( CheckNL )
      {
         const int PatchLocalID = FB_Aux_PatchLocalID( cell_idx );
         if ( NonLeaf[PatchLocalID] )
            Aux_Error( ERROR_INFO, "The particle in the feedback routine should not be in the non-leaf patch !!\n" );
      }
#     endif // #ifdef GAMER_DEBUG


//    2. check the feedback condition

//    2.1 to give feedback, the particle has to be a star particle
      if ( par_type != PTYPE_STAR )                               continue;

//    2.2 to give feedback, the particle has to be created in the simulation
      if ( par_creTime < 0.0 )                                    continue;

//    2.3 to give feedback, the particle mass has to be positive
      if ( par_mass <= 0.0 )                                      continue;

//    2.4 to give feedback, the particle explosion time has to be in this step;
      if ( par_SNIITime < TimeOld  ||  par_SNIITime > TimeNew )   continue; // both par_SNIITime == TimeOld and par_SNIITime ==  TimeNew are allowed

//    2.5 to give feedback, the particle has to be in the range to update fluid
      if ( cell_idx[0] + (maxfbDiameter  -isOnLeftCell[0])/2 <        FB_GHOST_SIZE  ||
           cell_idx[0] - (maxfbDiameter-1+isOnLeftCell[0])/2 >= PS2 + FB_GHOST_SIZE  ||
           cell_idx[1] + (maxfbDiameter  -isOnLeftCell[1])/2 <        FB_GHOST_SIZE  ||
           cell_idx[1] - (maxfbDiameter-1+isOnLeftCell[1])/2 >= PS2 + FB_GHOST_SIZE  ||
           cell_idx[2] + (maxfbDiameter  -isOnLeftCell[2])/2 <        FB_GHOST_SIZE  ||
           cell_idx[2] - (maxfbDiameter-1+isOnLeftCell[2])/2 >= PS2 + FB_GHOST_SIZE )  continue;


//    3. calculate the feedback

//    3.1 energy feedback
      real SNII_DepositedEnergy = FB_RESOLVED_SNEII_EJECT_ENGY;

//    3.2 mass feedback
      real SNII_DepositedMass   = MIN( par_mass, FB_RESOLVED_SNEII_EJECT_MASS );                   // particle mass cannot be less than zero

//    3.3 metal feedback
      real SNII_DepositedMetal  = MIN( par_metalMass,
                                       MIN( FB_RESOLVED_SNEII_EJECT_METAL, SNII_DepositedMass ) ); // ejected metal mass cannot be larger than ejected total mass

//    3.4 region to apply feedback
      int    fbDiameter  = 0;        // diameter of the region to apply feedback
      double fbFluidMass = 0.0;      // enclosed mass in the region to apply feedback

//    increase the feedback region to include enough gas mass
      for (int diameter=1; diameter<=maxfbDiameter; diameter++)
      {
#        ifdef GAMER_DEBUG
         if ( fbDepositWeighting[diameter-1] == NULL )
            Aux_Error( ERROR_INFO, "fbDepositWeighting[%d] == NULL !!\n", diameter-1 );
#        endif // #ifdef GAMER_DEBUG

//       calculate offsets based on cell position
         const int diameterMinus1 = diameter-1;
         const int offsets[3][2]  = { { (diameterMinus1+isOnLeftCell[0])/2, (diameter-isOnLeftCell[0])/2 },
                                      { (diameterMinus1+isOnLeftCell[1])/2, (diameter-isOnLeftCell[1])/2 },
                                      { (diameterMinus1+isOnLeftCell[2])/2, (diameter-isOnLeftCell[2])/2 } };

//       whether any cell in the region are in the forbidden patch (coarse patch and non-leaf patch)
         bool   touchForbiddenPat = false;

//       sum the total enclosed mass of the cells in the range
         double enclosedMass      = 0.0;

//       loop through cells in the region
         for (int dk=-offsets[2][0]; dk<=offsets[2][1] && !touchForbiddenPat; dk++) { const int k = cell_idx[2] + dk;   const int k_w = offsets[2][0] + dk;
         for (int dj=-offsets[1][0]; dj<=offsets[1][1] && !touchForbiddenPat; dj++) { const int j = cell_idx[1] + dj;   const int j_w = offsets[1][0] + dj;
         for (int di=-offsets[0][0]; di<=offsets[0][1] && !touchForbiddenPat; di++) { const int i = cell_idx[0] + di;   const int i_w = offsets[0][0] + di;

//          skip cells with non-positive weighting
            if ( fbDepositWeighting[diameterMinus1][k_w][j_w][i_w] <= 0.0 )   continue;

//          add cell mass to total
            enclosedMass += Fluid_In[DENS][k][j][i]*dv;

//          check if the cell included is in the coarse patch
            if ( CheckCF )
            {
               const int ijk[3] = { i, j, k };
               const int CellPatchRelPos = FB_Aux_CellPatchRelPos( ijk );
               if ( CellPatchRelPos != -1  &&  CoarseFine[CellPatchRelPos] )   touchForbiddenPat = true;

            } // if ( CheckCF )

//          check if the cell included is in the non-leaf patch
            if ( CheckNL )
            {
               const int ijk[3] = { i, j, k };
               const int PatchLocalID = FB_Aux_PatchLocalID( ijk );
               if ( NonLeaf[PatchLocalID] )   touchForbiddenPat = true;

            } // if ( CheckNL )

         }}} // i,j,k

//       stop expanding the region and go back to the previous one when it touches the coarse patch or non-leaf patch
         if ( touchForbiddenPat )   break;

//       store the current diameter and mass
         fbDiameter  = diameter;
         fbFluidMass = enclosedMass;

//       stop expanding the region when the enclosed mass is larger than the given thredshold
         if ( enclosedMass >= FB_RESOLVED_SNEII_MIN_M_GAS )   break;

      } // for (int diameter=1; diameter<=maxfbDiameter; diameter++)

//    to give feedback, the particle has to be in the range to update fluid
      if ( cell_idx[0] + (fbDiameter  -isOnLeftCell[0])/2 <        FB_GHOST_SIZE  ||
           cell_idx[0] - (fbDiameter-1+isOnLeftCell[0])/2 >= PS2 + FB_GHOST_SIZE  ||
           cell_idx[1] + (fbDiameter  -isOnLeftCell[1])/2 <        FB_GHOST_SIZE  ||
           cell_idx[1] - (fbDiameter-1+isOnLeftCell[1])/2 >= PS2 + FB_GHOST_SIZE  ||
           cell_idx[2] + (fbDiameter  -isOnLeftCell[2])/2 <        FB_GHOST_SIZE  ||
           cell_idx[2] - (fbDiameter-1+isOnLeftCell[2])/2 >= PS2 + FB_GHOST_SIZE )  continue;


//    4. record the information of SNeII event before applying feedback
      if ( FB_RESOLVED_SNEII_RECORD  &&  FB_Aux_CellPatchRelPos( cell_idx ) == -1 )
      {
#        ifdef GAMER_DEBUG
         if ( recordSNeII      == NULL )   Aux_Error( ERROR_INFO, "recordSNeII == NULL !!\n" );
         if ( numRecSNeII      == NULL )   Aux_Error( ERROR_INFO, "numRecSNeII == NULL !!\n" );
         if ( recordSNeII[TID] == NULL )   Aux_Error( ERROR_INFO, "recordSNeII[TID] == NULL !!\n" );
#        endif // #ifdef GAMER_DEBUG

         if ( numRecSNeII[TID] >= maxNumRecSNeII )   Aux_Error( ERROR_INFO, "Number of SNeII to record is larger than maximum (%d)!!\n", maxNumRecSNeII );

         int nVar = 0;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = TimeOld;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = TimeNew;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_SNIITime;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = SNII_DepositedEnergy;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = SNII_DepositedMass;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = SNII_DepositedMetal;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = fbDiameter*dh;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = fbFluidMass;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_idx;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_mass;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_pos[0];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_pos[1];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_pos[2];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_vel[0];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_vel[1];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_vel[2];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = par_metalMass;
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = Fluid_In[DENS][cell_idx[2]][cell_idx[1]][cell_idx[0]];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = Fluid_In[MOMX][cell_idx[2]][cell_idx[1]][cell_idx[0]];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = Fluid_In[MOMY][cell_idx[2]][cell_idx[1]][cell_idx[0]];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = Fluid_In[MOMZ][cell_idx[2]][cell_idx[1]][cell_idx[0]];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = Fluid_In[ENGY][cell_idx[2]][cell_idx[1]][cell_idx[0]];
         recordSNeII[TID][ nVarRecSNeII*numRecSNeII[TID] + (nVar++) ] = ( UseMetal ) ? Fluid_In[Idx_Metal][cell_idx[2]][cell_idx[1]][cell_idx[0]] : 0.0;
         numRecSNeII[TID] += 1;

         if ( nVar != nVarRecSNeII )   Aux_Error( ERROR_INFO, "Number of variables to record (%d) is not expected (%d)!!\n", nVar, nVarRecSNeII );

      } // if ( FB_RESOLVED_SNEII_RECORD  &&  FB_Aux_CellPatchRelPos( cell_idx ) == -1 )


//    5. update

//    5.1 update the particle
      ParAttFlt[Idx_ParSNIITime ][par_idx] *= -1;                // negative ParSNIITime will be recoreded and indicates the SN has exploded
      ParAttFlt[PAR_MASS        ][par_idx] -= (real_par)SNII_DepositedMass;
      if ( UseMetal )
      ParAttFlt[Idx_ParMetalFrac][par_idx]  = (par_metalMass - (real_par)SNII_DepositedMetal)/ParAttFlt[PAR_MASS][par_idx];

//    5.2 update the fluid
      const int fbDiameterMinus1 = fbDiameter-1;
      const int fbOffsets[3][2]  = { { (fbDiameterMinus1+isOnLeftCell[0])/2, (fbDiameter-isOnLeftCell[0])/2 },
                                     { (fbDiameterMinus1+isOnLeftCell[1])/2, (fbDiameter-isOnLeftCell[1])/2 },
                                     { (fbDiameterMinus1+isOnLeftCell[2])/2, (fbDiameter-isOnLeftCell[2])/2 } };

      for (int dk=-fbOffsets[2][0]; dk<=fbOffsets[2][1]; dk++) { const int k = cell_idx[2] + dk;   const int k_w = fbOffsets[2][0] + dk;
      for (int dj=-fbOffsets[1][0]; dj<=fbOffsets[1][1]; dj++) { const int j = cell_idx[1] + dj;   const int j_w = fbOffsets[1][0] + dj;
      for (int di=-fbOffsets[0][0]; di<=fbOffsets[0][1]; di++) { const int i = cell_idx[0] + di;   const int i_w = fbOffsets[0][0] + di;
#        ifdef DUAL_ENERGY
#        if   ( DUAL_ENERGY == DE_ENPY )
         const real flu_Dens = Fluid_Out[DENS][k][j][i];
         const real flu_Dual = Fluid_Out[DUAL][k][j][i];
         const real flu_Pres = Hydro_DensDual2Pres( flu_Dens, flu_Dual,
                                                    EoS_AuxArray_Flt[1], false, NULL_REAL );
         const real flu_Eint = EoS_DensPres2Eint_CPUPtr( flu_Dens, flu_Pres,
                                                         NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#        endif
#        endif // #ifdef DUAL_ENERGY ... else ...

         Fluid_Out[DENS     ][k][j][i] += SNII_DepositedMass                    * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;
         Fluid_Out[MOMX     ][k][j][i] += SNII_DepositedMass * (real)par_vel[0] * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;
         Fluid_Out[MOMY     ][k][j][i] += SNII_DepositedMass * (real)par_vel[1] * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;
         Fluid_Out[MOMZ     ][k][j][i] += SNII_DepositedMass * (real)par_vel[2] * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;
         Fluid_Out[ENGY     ][k][j][i] += SNII_DepositedEnergy                  * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;
         if ( UseMetal )
         Fluid_Out[Idx_Metal][k][j][i] += SNII_DepositedMetal                   * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;

#        ifdef DUAL_ENERGY
#        if   ( DUAL_ENERGY == DE_ENPY )
         const real Eint                = flu_Eint + SNII_DepositedEnergy       * fbDepositWeighting[fbDiameterMinus1][k_w][j_w][i_w] / dv;
         const real Pres                = EoS_DensEint2Pres_CPUPtr( Fluid_Out[DENS][k][j][i], Eint, NULL,
                                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         Fluid_Out[DUAL     ][k][j][i]  = Hydro_DensPres2Dual( Fluid_Out[DENS][k][j][i], Pres, EoS_AuxArray_Flt[1] );
#        endif
#        endif // #ifdef DUAL_ENERGY

      }}} // i,j,k

   } // for (int p=0; p<NPar; p++)


   return GAMER_SUCCESS;

} // FUNCTION : FB_Resolved_SNeII



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_Resolved_SNeII
// Description :  Free the resources used by the resolved Type II SNe feedback
//
// Note        :  1. Invoked by FB_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_Resolved_SNeII()
{

// free memory
   for (int diameter=1; diameter<=maxfbDiameter; diameter++)
      Aux_DeallocateArray3D( fbDepositWeighting[diameter-1] );


   if ( FB_RESOLVED_SNEII_RECORD )
   {
      for (int tid=0; tid<OMP_NTHREAD; tid++)   delete [] recordSNeII[tid];

      delete [] numRecSNeII;
      delete [] recordSNeII;
   }

} // FUNCTION : FB_End_Resolved_SNeII



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_Resolved_SNeII
// Description :  Initialize the resolved Type II SNe feedback
//
// Note        :  1. Invoked by FB_Init()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init_Resolved_SNeII()
{

// check
   if ( Idx_ParSNIITime == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParSNIITime is undefined !!\n" );


// set the corresponding weighting distribution for each case of feedback region diameter
   for (int diameter=1; diameter<=maxfbDiameter; diameter++)
   {
//    allocate a memory of cube for each diameter
      Aux_AllocateArray3D( fbDepositWeighting[diameter-1], diameter, diameter, diameter );

      double cubeCent       = 0.5*diameter;   // center coordinate of the cube of width = diameter
      double fbRadius       = 0.5*diameter;   // target radius to apply feedback within
      real   totalWeighting = 0.0;            // for normalization of the feedback

//    set the weighting for each cell in the cube
      for (int k=0; k<diameter; k++)
      for (int j=0; j<diameter; j++)
      for (int i=0; i<diameter; i++)
      {
//       set a uniform distribution in the inscribed sphere
         fbDepositWeighting[diameter-1][k][j][i] = ( SQR( i+0.5-cubeCent ) + SQR( j+0.5-cubeCent ) + SQR( k+0.5-cubeCent ) <= SQR( fbRadius ) )
                                                   ? 1.0 : 0.0;
//       sum the total weighting
         totalWeighting += fbDepositWeighting[diameter-1][k][j][i];

      } // i,j,k

//    normalize the weighting
      for (int k=0; k<diameter; k++)
      for (int j=0; j<diameter; j++)
      for (int i=0; i<diameter; i++)
      {
         fbDepositWeighting[diameter-1][k][j][i] /= totalWeighting;

      } // i,j,k

   } // for (int diameter=1; diameter<=maxfbDiameter; diameter++)


// prepare the buffer array to collect the information to record
   if ( FB_RESOLVED_SNEII_RECORD )
   {
      numRecSNeII = new int     [OMP_NTHREAD];
      recordSNeII = new double* [OMP_NTHREAD];

      for (int tid=0; tid<OMP_NTHREAD; tid++)
      {
         numRecSNeII[tid] = 0;
         recordSNeII[tid] = new double [maxNumRecSNeII*nVarRecSNeII];

      } // for (int tid=0; tid<OMP_NTHREAD; tid++)
   }

} // FUNCTION : FB_Init_Resolved_SNeII



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_FB_Resolved_SNeII
// Description :  Record the information of FB_Resolved_SNeII
//
// Note        :  1. Invoked by FB_AdvanceDt()
//                2. The information is gathered parallelly by different MPI ranks and OpenMP threads
//                3. After write to file, the arrays will be reset and rewritten for the next update
//
// Parameter   :  lv         : Target refinement level
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Record_FB_Resolved_SNeII( const int lv )
{

// nothing to do if FB_RESOLVED_SNEII_RECORD is disabled
   if ( !FB_RESOLVED_SNEII_RECORD )   return;

// set the filename
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__FB_Resolved_SNeII", OUTPUT_DIR );

// write to the file rank by rank
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
//       open file
         FILE *File = fopen( FileName, "a" );

//       output header
         static bool FirstTime = true;
         if ( TargetMPIRank == 0  &&  FirstTime )
         {
            fprintf( File, "#%5s%6s%6s%16s%16s",
                     "Rank", "TID", "lv", "TimeOld", "TimeNew" );
            fprintf( File, "%16s%16s%16s%16s%16s%16s",
                     "SNII_Time", "SNII_Energy", "SNII_Mass", "SNII_Metal", "FB_Diameter", "FB_Flu_Mass" );
            fprintf( File, "%16s%16s%16s%16s%16s",
                     "Par_ID", "Par_Mass", "Par_PosX", "Par_PosY", "Par_PosZ" );
            fprintf( File, "%16s%16s%16s%16s",
                     "Par_VelX", "Par_VelY", "Par_VelZ", "Par_MetalMass" );
            fprintf( File, "%16s%16s%16s%16s%16s%16s\n",
                     "Flu_Dens", "Flu_MomX", "Flu_MomY", "Flu_MomZ", "Flu_Engy", "Flu_Metal" );

            FirstTime = false;

         } // if ( TargetMPIRank == 0  &&  FirstTime )

//       output content
         for (int tid=0; tid<OMP_NTHREAD; tid++)
         {
            for (int iSNII=0; iSNII<numRecSNeII[tid]; iSNII++)
            {
               fprintf( File, "%6d%6d%6d", MPI_Rank, tid, lv );

               for (int iVar=0; iVar<nVarRecSNeII; iVar++)
                  fprintf( File, "%16.8e", recordSNeII[tid][ nVarRecSNeII*iSNII + iVar ] );

               fprintf( File, "\n" );

            } // for (int iSNII=0; iSNII<numRecSNeII[tid]; iSNII++)
         } // for (int tid=0; tid<OMP_NTHREAD; tid++)

//       close file
         fclose( File );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)

// reset the number of SNeII to zero
   for (int tid=0; tid<OMP_NTHREAD; tid++)   numRecSNeII[tid] = 0;

} // FUNCTION : Record_FB_Resolved_SNeII



#endif // #ifdef FEEDBACK
