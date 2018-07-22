#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Prepare
// Description :  Fill up the input host array "h_Che_Array" for the CPU/GPU Grackle solver
//
// Note        :  1. Prepare CHE_NPREP variables
//                   --> CHE_NPREP = 4 currently
//                   --> [mass density, specific internal energy, kinematic energy density, metal density]
//                2. This function always prepares the latest FluSg data
//
// Parameter   :  lv          : Target refinement level
//                h_Che_Array : Host array to store the prepared data
//                NPG         : Number of patch groups prepared at a time
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Grackle_Prepare( const int lv, real h_Che_Array[], const int NPG, const int *PID0_List )
{

   const int  hCheIdx_Dens    = 0;
   const int  hCheIdx_sEint   = 1;
   const int  hCheIdx_Ek      = 2;
   const int  hCheIdx_Metal   = 3;
   const real dh              = (real)amr->dh[lv];
   const int  Size1pg         = CUBE(PS2);
   const int  Size1v          = NPG*Size1pg;
#  ifdef DUAL_ENERGY
   const real  Gamma_m1       = GAMMA - (real)1.0;
   const real _Gamma_m1       = (real)1.0 / Gamma_m1;
   const bool CheckMinPres_No = false;
#  endif

   int  idx_pg, PID, PID0;    // idx_pg: array indices within a patch group
   real Dens, Px, Py, Pz, Etot, _Dens, Ek, sEint;

   real *Ptr_Dens0  = h_Che_Array + hCheIdx_Dens *Size1v;
   real *Ptr_sEint0 = h_Che_Array + hCheIdx_sEint*Size1v;
   real *Ptr_Ek0    = h_Che_Array + hCheIdx_Ek   *Size1v;
   real *Ptr_Metal0 = h_Che_Array + hCheIdx_Metal*Size1v;

   real *Ptr_Dens=NULL, *Ptr_sEint=NULL, *Ptr_Ek=NULL, *Ptr_Metal=NULL;


// check
#  ifdef GAMER_DEBUG
   if ( GRACKLE_METAL  &&  Idx_Metal == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_Metal is undefined for \"GRACKLE_METAL\" !!\n" );
#  endif


#  pragma omp parallel for private( idx_pg, PID, PID0, Dens, Px, Py, Pz, Etot, _Dens, Ek, sEint, \
                                    Ptr_Dens, Ptr_sEint, Ptr_Ek, Ptr_Metal ) schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0      = PID0_List[TID];
      idx_pg    = 0;
      Ptr_Dens  = Ptr_Dens0  + TID*Size1pg;
      Ptr_sEint = Ptr_sEint0 + TID*Size1pg;
      Ptr_Ek    = Ptr_Ek0    + TID*Size1pg;
      Ptr_Metal = Ptr_Metal0 + TID*Size1pg;

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;

         for (int idx_p=0; idx_p<CUBE(PS1); idx_p++)
         {
            Dens  = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][0][0] + idx_p );
            Px    = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][0][0] + idx_p );
            Py    = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][0][0] + idx_p );
            Pz    = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][0][0] + idx_p );
            Etot  = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][0][0] + idx_p );
            _Dens = (real)1.0 / Dens;
            Ek    = (real)0.5*( SQR(Px) + SQR(Py) + SQR(Pz) )*_Dens;

//          use the dual-energy variable to calculate the internal energy if applicable
#           ifdef DUAL_ENERGY

#           if   ( DUAL_ENERGY == DE_ENPY )
            sEint = CPU_DensEntropy2Pres( Dens, *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENPY][0][0] + idx_p ), Gamma_m1,
                                          CheckMinPres_No, NULL_REAL )*_Dens*_Gamma_m1;
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

#           else
            sEint = ( Etot - Ek )*_Dens;
#           endif // #ifdef DUAL_ENERGY ... else

            Ptr_Dens [idx_pg] = Dens;
            Ptr_sEint[idx_pg] = sEint;
            Ptr_Ek   [idx_pg] = Ek;

//          prepare metallicity if metal cooling is enabled
            if ( GRACKLE_METAL )
            Ptr_Metal[idx_pg] = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[Idx_Metal][0][0] + idx_p );

            idx_pg ++;
         } // for (int idx_p=0; idx_p<CUBE(PS1); idx_p++)

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

   if ( GRACKLE_MODE == GRACKLE_MODE_ORI )
   {
      Che_FieldData->density         = Ptr_Dens0;
      Che_FieldData->internal_energy = Ptr_sEint0;
      Che_FieldData->grid_dx         = dh;

      if ( GRACKLE_METAL )
      Che_FieldData->metal_density   = Ptr_Metal0;
   }

} // FUNCTION : Grackle_Prepare



#endif // #ifdef SUPPORT_GRACKLE
