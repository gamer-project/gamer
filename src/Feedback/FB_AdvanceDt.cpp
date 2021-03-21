#include "GAMER.h"

#ifdef FEEDBACK




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_AdvanceDt
// Description :  Feedback from particles to grids
//
// Note        :  1. Invoked by EvolveLevel()
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                SaveSg_Flu : Sandglass to store the updated fluid data
//                SaveSg_Mag : Sandglass to store the updated B field
//
// Return      :  Update both grids and particles
//-------------------------------------------------------------------------------------------------------
void FB_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const int SaveSg_Flu, const int SaveSg_Mag )
{


} // FUNCTION : FB_AdvanceDt



#endif // #ifdef FEEDBACK
