#include "GAMER.h"

#ifdef FEEDBACK



// function pointers to be set by FB_Init_User_Template()
extern void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                            real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[] );
extern void (*FB_End_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_User_Template
// Description :  Template of a user-defined feedback
//
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAtt[], respectively
//                2. Must use ParSortID[] to access ParAtt[]
//                   --> ParAtt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
//                3. To ensure the consistency of random numbers, one must call the random number generator for
//                   ALL particles, including those too far away to affect the target patch group
//                4. No need to worry about the periodic boundary condition here
//                   --> Particle positions have been remapped in FB_AdvanceDt()
//                5. CoarseFine[] records the coarse-fine boundaries along the 26 sibling directions, defined as
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
//                6. Invoked by FB_AdvanceDt()
//                7. Linked to FB_User_Ptr in FB_Init_User_Template()
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                NPar       : Number of particles
//                ParSortID  : Sorted particle IDs
//                ParAtt     : Particle attribute arrays
//                Fluid      : Array to store the input/output fluid data
//                             --> Array size is fixed to PS2^3
//                EdgeL      : Left edge of Fluid[]
//                             --> Right edge is given by EdgeL[]+PS2*dh
//                dh         : Cell size of Fluid[]
//                CoarseFine : Coarse-fine boundaries along the 26 sibling directions
//
// Return      :  Fluid, ParAtt
//-------------------------------------------------------------------------------------------------------
void FB_User_Template( const int lv, const double TimeNew, const double TimeOld, const double dt,
                       const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                       real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid == NULL )    Aux_Error( ERROR_INFO, "Fluid == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAtt == NULL )      Aux_Error( ERROR_INFO, "ParAtt == NULL for NPar = %d !!\n", NPar );
   }
#  endif // #ifdef GAMER_DEBUG




} // FUNCTION : FB_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_User_Template
// Description :  Free the resources used by the user-specified feedback
//
// Note        :  1. Invoked by FB_End()
//                2. Linked to FB_End_User_Ptr in FB_Init_User_Template()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_User_Template()
{


} // FUNCTION : FB_End_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_User_Template
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
void FB_Init_User_Template()
{

   FB_User_Ptr     = FB_User_Template;
   FB_End_User_Ptr = FB_End_User_Template;

} // FUNCTION : FB_Init_User_Template



#endif // #ifdef FEEDBACK
