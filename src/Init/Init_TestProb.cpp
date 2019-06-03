#include "GAMER.h"


// ******************************************************************
// add the new test problem function prototypes here
// ******************************************************************
void Init_TestProb_Hydro_BlastWave();
void Init_TestProb_Hydro_AcousticWave();
void Init_TestProb_Hydro_Bondi();
void Init_TestProb_Hydro_ClusterMerger_vs_Flash();
void Init_TestProb_Hydro_AGORA_IsolatedGalaxy();
void Init_TestProb_Hydro_Caustic();
void Init_TestProb_Hydro_SphericalCollapse();
void Init_TestProb_Hydro_KelvinHelmholtzInstability();
void Init_TestProb_Hydro_Riemann();
void Init_TestProb_Hydro_CollidingJets();
void Init_TestProb_Hydro_Plummer();
void Init_TestProb_Hydro_Gravity();

void Init_TestProb_SRHydro_AcousticWave();
void Init_TestProb_SRHydro_BlastWave();
void Init_TestProb_SRHydro_Riemann();
void Init_TestProb_SRHydro_DoubleMachReflection();
void Init_TestProb_SRHydro_Jets();
void Init_TestProb_SRHydro_PrecessedJet();
void Init_TestProb_SRHydro_DiffPrecessedJet();
void Init_TestProb_SRHydro_WScaling_BlastWave();

void Init_TestProb_ELBDM_ExtPot();




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
      case TESTPROB_NONE :                                                                                  break;
#     if ( MODEL == HYDRO )
      case TESTPROB_HYDRO_BLAST_WAVE :                   Init_TestProb_Hydro_BlastWave();                   break;
      case TESTPROB_HYDRO_ACOUSTIC_WAVE :                Init_TestProb_Hydro_AcousticWave();                break;
//      case TESTPROB_HYDRO_BONDI :                        Init_TestProb_Hydro_Bondi();                       break;
      case TESTPROB_HYDRO_CLUSTER_MERGER_VS_FLASH :      Init_TestProb_Hydro_ClusterMerger_vs_Flash();      break;
      case TESTPROB_HYDRO_AGORA_ISOLATED_GALAXY :        Init_TestProb_Hydro_AGORA_IsolatedGalaxy();        break;
      case TESTPROB_HYDRO_CAUSTIC :                      Init_TestProb_Hydro_Caustic();                     break;
      case TESTPROB_HYDRO_SPHERICAL_COLLAPSE :           Init_TestProb_Hydro_SphericalCollapse();           break;
      case TESTPROB_HYDRO_KELVIN_HELMHOLTZ_INSTABILITY : Init_TestProb_Hydro_KelvinHelmholtzInstability();  break;
      case TESTPROB_HYDRO_RIEMANN :                      Init_TestProb_Hydro_Riemann();                     break;
      case TESTPROB_HYDRO_COLLIDING_JETS :               Init_TestProb_Hydro_CollidingJets();               break;
      case TESTPROB_HYDRO_PLUMMER :                      Init_TestProb_Hydro_Plummer();                     break;
      case TESTPROB_HYDRO_GRAVITY :                      Init_TestProb_Hydro_Gravity();                     break;
#     elif ( MODEL == SR_HYDRO )
      case TESTPROB_SRHYDRO_ACOUSTIC_WAVE :              Init_TestProb_SRHydro_AcousticWave();                break;
      case TESTPROB_SRHYDRO_BLAST_WAVE :                 Init_TestProb_SRHydro_BlastWave();                 break;
      case TESTPROB_SRHYDRO_WEAK_SCALING_BLAST_WAVE :    Init_TestProb_SRHydro_WScaling_BlastWave();        break;
      case TESTPROB_SRHYDRO_RIEMANN :                    Init_TestProb_SRHydro_Riemann();                   break;
      case TESTPROB_SRHYDRO_DOUBLE_MACH_REFLECTION :     Init_TestProb_SRHydro_DoubleMachReflection();      break;
      case TESTPROB_SRHYDRO_JETS:                        Init_TestProb_SRHydro_Jets();                      break;
      case TESTPROB_SRHYDRO_PRECESSED_JET:               Init_TestProb_SRHydro_PrecessedJet();              break;
      case TESTPROB_SRHYDRO_DIFFERENTIAL_PRECESSED_JET : Init_TestProb_SRHydro_DiffPrecessedJet();          break;
#     elif ( MODEL == ELBDM )
      case TESTPROB_ELBDM_EXTPOT :                       Init_TestProb_ELBDM_ExtPot();                      break;
#     endif
      default: Aux_Error( ERROR_INFO, "unsupported TESTPROB_ID (%d) !!\n", TESTPROB_ID );
   } // switch( TESTPROB_ID )

} // FUNCTION : Init_TestProb
