#include "GAMER.h"

#if ( MODEL == HYDRO )

extern double  *Jet_Radius;
extern double  *Jet_HalfHeight;
extern double  *Jet_SrcVel;
extern double  *Jet_SrcDens;
extern double  *Jet_SrcTemp;
extern double (*Jet_Vec)[3];
extern double (*Jet_CenOffset)[3];
extern double  *Jet_SrcEint;
extern double (*Jet_Cen)[3];
extern double  *Jet_WaveK;




//-------------------------------------------------------------------------------------------------------
// Function    :  End_TestProb
// Description :  Release memory
//
// Note        :  None
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_TestProb()
{

   delete [] Jet_Radius;
   delete [] Jet_HalfHeight;
   delete [] Jet_SrcVel;
   delete [] Jet_SrcDens;
   delete [] Jet_SrcTemp;
   delete [] Jet_Vec;
   delete [] Jet_CenOffset;
   delete [] Jet_SrcEint;
   delete [] Jet_Cen;
   delete [] Jet_WaveK;

} // FUNCTION : End_TestProb



#endif // #if ( MODEL == HYDRO )
