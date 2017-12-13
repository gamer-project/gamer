#include "GAMER.h"

#if ( MODEL == ELBDM )

extern real  *Soliton_DensProfileR;
extern real  *Soliton_DensProfile;
extern real  *Soliton_Scale;
extern real (*Soliton_Center)[3];




//-------------------------------------------------------------------------------------------------------
// Function    :  End_TestProb
// Description :  Release memory for the ELBDM soliton test
//
// Note        :  None
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void End_TestProb()
{  

   if ( Soliton_DensProfileR != NULL )
   {
      delete [] Soliton_DensProfileR;
      Soliton_DensProfileR = NULL;
   }

   if ( Soliton_DensProfile != NULL )
   {
      delete [] Soliton_DensProfile;
      Soliton_DensProfile = NULL;
   }

   if ( Soliton_Scale != NULL )
   {
      delete [] Soliton_Scale;
      Soliton_Scale = NULL;
   }

   if ( Soliton_Center != NULL )
   {
      delete [] Soliton_Center;
      Soliton_Center = NULL;
   }

} // FUNCTION : End_TestProb



#endif // #if ( MODEL == ELBDM )
