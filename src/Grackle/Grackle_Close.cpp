#include "Copyright.h"
#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Close
// Description :  Copy the specific internal energy updated by the CPU/GPU Grackle solver back to the
//                patch pointers
//
// Note        :  1. Use SaveSg to determine where to store the data
//                   --> Currently it's set to the same Sg as the fluid data when calling
//                       Grackle_AdvanceDt() in EvolveLevel()
//
// Parameter   :  lv          : Target refinement level
//                SaveSg      : Sandglass to store the updated data
//                h_Che_Array : Host array storing the updated data
//                NPG         : Number of patch groups to store the updated data
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Grackle_Close( const int lv, const int SaveSg, const real h_Che_Array[][CHE_NPREP][ CUBE(PS1) ],
                    const int NPG, const int *PID0_List )
{

   const int Idx_Dens   = 0;
   const int Idx_sEint  = 1;
   const int Idx_Ek     = 2;
   const real  Gamma_m1 = GAMMA - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;

   int  N, PID, PID0;
   real Dens, Pres;


#  pragma omp parallel for private( N, PID, PID0, Dens, Pres ) schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

         for (int t=0; t<CUBE(PS1); t++)
         {
//          apply the minimum pressure check
            Dens = h_Che_Array[N][Idx_Dens ][t];
            Pres = h_Che_Array[N][Idx_sEint][t]*Dens*Gamma_m1;
            Pres = CPU_CheckMinPres( Pres, MIN_PRES );

//          update the total energy density
            *( amr->patch[SaveSg][lv][PID]->fluid[ENGY][0][0] + t ) = Pres*_Gamma_m1 + h_Che_Array[N][Idx_Ek][t];

//          update the dual-energy variable to be consistent with the updated pressure
#           ifdef DUAL_ENERGY
#           if   ( DUAL_ENERGY == DE_ENPY )
            *( amr->patch[SaveSg][lv][PID]->fluid[ENPY][0][0] + t ) = CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );

#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif
#           endif // #ifdef DUAL_ENERGY
         } // for (int t=0; t<CUBE(PS1); t++)
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Grackle_Close



#endif // #ifdef SUPPORT_GRACKLE
