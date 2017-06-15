#include "Copyright.h"
#include "GAMER.h"


// ******************************************************************
// add the new test problem function prototypes here
// ******************************************************************
void Init_TestProb_BlastWave();
void Init_TestProb_AcousticWave();




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize the target test problem
//
// Note        :  1. Use TESTPROB_ID to choose the target test problem
//                2. All test problem IDs are defined in "include/Typedef.h"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

// ******************************************************************
// add the new test problem IDs here
// ******************************************************************
   switch ( TESTPROB_ID )
   {
      case TESTPROB_NONE :                                           break;
      case TESTPROB_BLAST_WAVE :    Init_TestProb_BlastWave();       break;
      case TESTPROB_ACOUSTIC_WAVE : Init_TestProb_AcousticWave();    break;

      default: Aux_Error( ERROR_INFO, "unsupported TESTPROB_ID (%d) !!\n", TESTPROB_ID );
   } // switch( TESTPROB_ID )

} // FUNCTION : Init_TestProb
