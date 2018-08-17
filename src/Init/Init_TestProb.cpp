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

void Init_TestProb_ELBDM_ExtPot();
void Init_TestProb_ELBDM_JeansInstabilityComoving();
void Init_TestProb_ELBDM_JeansInstabilityPhysical();
void Init_TestProb_ELBDM_Soliton();
void Init_TestProb_ELBDM_SelfSimilarHalo();
void Init_TestProb_ELBDM_VortexPairRotating();




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

      case TESTPROB_HYDRO_BLAST_WAVE :                   Init_TestProb_Hydro_BlastWave();                   break;
      case TESTPROB_HYDRO_ACOUSTIC_WAVE :                Init_TestProb_Hydro_AcousticWave();                break;
//    case TESTPROB_HYDRO_BONDI :                        Init_TestProb_Hydro_Bondi();                       break;
      case TESTPROB_HYDRO_CLUSTER_MERGER_VS_FLASH :      Init_TestProb_Hydro_ClusterMerger_vs_Flash();      break;
      case TESTPROB_HYDRO_AGORA_ISOLATED_GALAXY :        Init_TestProb_Hydro_AGORA_IsolatedGalaxy();        break;
      case TESTPROB_HYDRO_CAUSTIC :                      Init_TestProb_Hydro_Caustic();                     break;
      case TESTPROB_HYDRO_SPHERICAL_COLLAPSE :           Init_TestProb_Hydro_SphericalCollapse();           break;
      case TESTPROB_HYDRO_KELVIN_HELMHOLTZ_INSTABILITY : Init_TestProb_Hydro_KelvinHelmholtzInstability();  break;
      case TESTPROB_HYDRO_RIEMANN :                      Init_TestProb_Hydro_Riemann();                     break;
      case TESTPROB_HYDRO_COLLIDING_JETS :               Init_TestProb_Hydro_CollidingJets();               break;
      case TESTPROB_HYDRO_PLUMMER :                      Init_TestProb_Hydro_Plummer();                     break;
      case TESTPROB_HYDRO_GRAVITY :                      Init_TestProb_Hydro_Gravity();                     break;

      case TESTPROB_ELBDM_EXTPOT :                       Init_TestProb_ELBDM_ExtPot();                      break;
      case TESTPROB_ELBDM_JEANS_INSTABILITY_COMOVING :   Init_TestProb_ELBDM_JeansInstabilityComoving();    break;
//    case TESTPROB_ELBDM_JEANS_INSTABILITY_PHYSICAL :   Init_TestProb_ELBDM_JeansInstabilityPhysical();    break;
      case TESTPROB_ELBDM_SOLITON :                      Init_TestProb_ELBDM_Soliton();                     break;
      case TESTPROB_ELBDM_SELF_SIMILAR_HALO :            Init_TestProb_ELBDM_SelfSimilarHalo();             break;
      case TESTPROB_ELBDM_VORTEX_PAIR_ROTATING :         Init_TestProb_ELBDM_VortexPairRotating();          break;

      default: Aux_Error( ERROR_INFO, "unsupported TESTPROB_ID (%d) !!\n", TESTPROB_ID );
   } // switch( TESTPROB_ID )

} // FUNCTION : Init_TestProb
