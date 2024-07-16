#include "GAMER.h"

// declare as static so that other functions cannot invoke them directly and must use the function pointers
static void Init_User_Template();
static void Init_User_AfterPoisson_Template();

// these function pointers must be set by a test problem initializer
void (*Init_User_Ptr)()              = NULL;
void (*Init_User_AfterPoisson_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_Template
// Description :  Template of user-defined initialization
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Gravitational potential on grids and particle acceleration have not been computed
//                   at this stage
//                   --> Use Init_User_AfterPoisson_Ptr() instead if this information is required
//                   --> It's OK to modify mass density on grids and particle mass/position here
//                       --> But grid distribution won't change unless you manually call the corresponding
//                           grid refinement routines here
//                       --> Modifying particle position requires special attention in order to ensure that
//                           all particles still reside in leaf patches. Do this only if you know what
//                           you are doing.
//                3. To add new particles, remember to call Par_AddParticleAfterInit()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_User_Template()
{

} // FUNCTION : Init_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_AfterPoisson_Template
// Description :  Template of user-defined initialization after invoking the Poisson solver
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_AfterPoisson_Ptr",
//                   which must be set by a test problem initializer
//                2. Unlike Init_User_Ptr(), this routine is invoked after the Poisson solver
//                   --> Do not modify mass density on grids and particle mass/position unless
//                        you manually invoke the Poisson solver again
//                   --> One typical usage of this routine is to reset particle velocity according to
//                       the gravity field
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_User_AfterPoisson_Template()
{

} // FUNCTION : Init_User_AfterPoisson_Template
