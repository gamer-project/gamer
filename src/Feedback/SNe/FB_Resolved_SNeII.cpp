#include "GAMER.h"

#ifdef FEEDBACK




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Resolved_SNeII
// Description :  Resolved type II supernovae explosion feedback
//
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAttFlt/Int[], respectively
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
int FB_Resolved_SNeII( const int lv, const double TimeNew, const double TimeOld, const double dt,
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

// cell volume
   const real dv = CUBE( dh );

// determine if metallicity is included
   const bool UseMetal = ( Idx_Metal != Idx_Undefined  &&  Idx_ParMetalFrac != Idx_Undefined );

   for (int p=0; p<NPar; p++)
   {

//    1. Prepare the input
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

      const int ii = (int)FLOOR( ( par_pos[0] - EdgeL[0] )/dh );                            // indices of cell where the particle is in
      const int jj = (int)FLOOR( ( par_pos[1] - EdgeL[1] )/dh );                            // Note it could be negative for particles in nearby patches
      const int kk = (int)FLOOR( ( par_pos[2] - EdgeL[2] )/dh );


//    2. Check the condition

//    2.1 The particle is in the range to update fluid
      if ( ii < FB_GHOST_SIZE  ||  ii >= (PS2 + FB_GHOST_SIZE)  ||
           jj < FB_GHOST_SIZE  ||  jj >= (PS2 + FB_GHOST_SIZE)  ||
           kk < FB_GHOST_SIZE  ||  kk >= (PS2 + FB_GHOST_SIZE) )  continue;

//    2.2 The particle is a star particle
      if ( par_type != PTYPE_STAR )                               continue;

//    2.3 The particle is created in the simulation
      if ( par_creTime < 0.0 )                                    continue;

//    2.4 The particle mass is positive
      if ( par_mass <= 0.0 )                                      continue;

//    2.5 The particle explosion time is in this step;
      if ( par_SNIITime < TimeOld  ||  par_SNIITime > TimeNew )   continue; // Both par_SNIITime == TimeOld and par_SNIITime ==  TimeNew are allowed


//    3. Calculate the feedback

//    3.1 Energy feedback
      real SNII_DepositedEnergy = FB_RESOLVED_SNEII_EJECT_ENGY;

//    3.2 Mass feedback
      real SNII_DepositedMass   = MIN( par_mass, FB_RESOLVED_SNEII_EJECT_MASS );                   // Particle mass cannot be less than zero

//    3.3 Metal feedback
      real SNII_DepositedMetal  = MIN( par_metalMass,
                                       MIN( FB_RESOLVED_SNEII_EJECT_METAL, SNII_DepositedMass ) ); // Ejected metal mass cannot be larger than Ejected total mass


//    4. Update

//    4.1 Update the particle
      ParAttFlt[Idx_ParSNIITime ][par_idx] *= -1;                // Negative ParSNIITime will be recoreded and indicates the SN has exploded
      ParAttFlt[PAR_MASS        ][par_idx] -= SNII_DepositedMass;
      if ( UseMetal )
      ParAttFlt[Idx_ParMetalFrac][par_idx]  = (par_metalMass - SNII_DepositedMetal)/ParAttFlt[PAR_MASS][par_idx];

//    4.2 Update the fluid
#     ifdef DUAL_ENERGY
#     if   ( DUAL_ENERGY == DE_ENPY )
      const real flu_Dens = Fluid[DENS][kk][jj][ii];
      const real flu_Dual = Fluid[DUAL][kk][jj][ii];
      const real flu_Pres = Hydro_DensDual2Pres( flu_Dens, flu_Dual,
                                                 EoS_AuxArray_Flt[1], false, NULL_REAL );
      const real flu_Eint = EoS_DensPres2Eint_CPUPtr( flu_Dens, flu_Pres,
                                                      NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#     endif
#     endif // #ifdef DUAL_ENERGY ... else ...

      Fluid[DENS     ][kk][jj][ii] += SNII_DepositedMass              / dv;
      Fluid[MOMX     ][kk][jj][ii] += SNII_DepositedMass * par_vel[0] / dv;
      Fluid[MOMY     ][kk][jj][ii] += SNII_DepositedMass * par_vel[1] / dv;
      Fluid[MOMZ     ][kk][jj][ii] += SNII_DepositedMass * par_vel[2] / dv;
      Fluid[ENGY     ][kk][jj][ii] += SNII_DepositedEnergy            / dv;
      if ( UseMetal )
      Fluid[Idx_Metal][kk][jj][ii] += SNII_DepositedMetal             / dv;

#     ifdef DUAL_ENERGY
#     if   ( DUAL_ENERGY == DE_ENPY )
      const real Eint               = flu_Eint + SNII_DepositedEnergy / dv;
      const real Pres               = EoS_DensEint2Pres_CPUPtr( Fluid[DENS][kk][jj][ii], Eint, NULL,
                                                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      Fluid[DUAL][kk][jj][ii]       = Hydro_DensPres2Dual( Fluid[DENS][kk][jj][ii], Pres,
                                                           EoS_AuxArray_Flt[1] );
#     endif
#     endif // #ifdef DUAL_ENERGY

   }


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
   if ( Idx_ParSNIITime == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParSNIITime is undefined !!\n" );

} // FUNCTION : FB_Init_Resolved_SNeII



#endif // #ifdef FEEDBACK
