#include "GAMER.h"

#ifdef PARTICLE


extern double *ClusterMerger_Gas_PresProf;
extern double *ClusterMerger_DM_MassProf;
extern double *ClusterMerger_DM_SigmaProf;

extern double *ClusterMerger_Gas_PresProf_R;
extern double *ClusterMerger_DM_MassProf_R;
extern double *ClusterMerger_DM_SigmaProf_R;




//-------------------------------------------------------------------------------------------------------
// Function    :  End_TestProb
// Description :  Release memory for the N particles force test
//
// Note        :  None
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_TestProb()
{

   delete [] ClusterMerger_Gas_PresProf;
   delete [] ClusterMerger_DM_MassProf;
   delete [] ClusterMerger_DM_SigmaProf;

   delete [] ClusterMerger_Gas_PresProf_R;
   delete [] ClusterMerger_DM_MassProf_R;
   delete [] ClusterMerger_DM_SigmaProf_R;

} // FUNCTION : End_TestProb



#endif // #ifdef PARTICLE
