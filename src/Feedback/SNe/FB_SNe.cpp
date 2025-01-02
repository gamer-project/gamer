#include "GAMER.h"

#ifdef FEEDBACK




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_SNe
// Description :  Supernova explosion feedback
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
int FB_SNe( const int lv, const double TimeNew, const double TimeOld, const double dt,
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


   return GAMER_SUCCESS;

} // FUNCTION : FB_SNe



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_SNe
// Description :  Free the resources used by the SNe feedback
//
// Note        :  1. Invoked by FB_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_SNe()
{

} // FUNCTION : FB_End_SNe



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_SNe
// Description :  Initialize the SNe feedback
//
// Note        :  1. Invoked by FB_Init()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init_SNe()
{

} // FUNCTION : FB_Init_SNe



#endif // #ifdef FEEDBACK
