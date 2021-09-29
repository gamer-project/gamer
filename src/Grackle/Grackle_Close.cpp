#include "GAMER.h"

#ifdef SUPPORT_GRACKLE


// global variables for accessing h_Che_Array[]
// --> declared in Init_MemAllocate_Grackle.cpp
extern int Che_NField;
extern int CheIdx_Dens;
extern int CheIdx_sEint;
extern int CheIdx_Ent;
extern int CheIdx_e;
extern int CheIdx_HI;
extern int CheIdx_HII;
extern int CheIdx_HeI;
extern int CheIdx_HeII;
extern int CheIdx_HeIII;
extern int CheIdx_HM;
extern int CheIdx_H2I;
extern int CheIdx_H2II;
extern int CheIdx_DI;
extern int CheIdx_DII;
extern int CheIdx_HDI;
extern int CheIdx_Metal;




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Close
// Description :  Copy the specific internal energy updated by the Grackle solver back to the
//                patch pointers
//
// Note        :  1. Use SaveSg to determine where to store the data
//                   --> Currently it's set to the same Sg as the fluid data when calling
//                       Grackle_AdvanceDt() in EvolveLevel()
//                2. Che_NField and the corresponding array indices in h_Che_Array[] (e.g., CheIdx_Dens)
//                   are declared and set by Init_MemAllocate_Grackle()
//
// Parameter   :  lv          : Target refinement level
//                SaveSg      : Sandglass to store the updated data
//                h_Che_Array : Host array storing the updated data
//                NPG         : Number of patch groups to store the updated data
//                PID0_List   : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Grackle_Close( const int lv, const int SaveSg, const real h_Che_Array[], const int NPG, const int *PID0_List )
{

   const int  Size1pg      = CUBE(PS2);
   const int  Size1v       = NPG*Size1pg;
   const real MassRatio_ep = Const_me / Const_mp;

   const real *Ptr_Dens0  = h_Che_Array + CheIdx_Dens *Size1v;
   const real *Ptr_sEint0 = h_Che_Array + CheIdx_sEint*Size1v;
   const real *Ptr_Ent0   = h_Che_Array + CheIdx_Ent  *Size1v;
   const real *Ptr_e0     = h_Che_Array + CheIdx_e    *Size1v;
   const real *Ptr_HI0    = h_Che_Array + CheIdx_HI   *Size1v;
   const real *Ptr_HII0   = h_Che_Array + CheIdx_HII  *Size1v;
   const real *Ptr_HeI0   = h_Che_Array + CheIdx_HeI  *Size1v;
   const real *Ptr_HeII0  = h_Che_Array + CheIdx_HeII *Size1v;
   const real *Ptr_HeIII0 = h_Che_Array + CheIdx_HeIII*Size1v;
   const real *Ptr_HM0    = h_Che_Array + CheIdx_HM   *Size1v;
   const real *Ptr_H2I0   = h_Che_Array + CheIdx_H2I  *Size1v;
   const real *Ptr_H2II0  = h_Che_Array + CheIdx_H2II *Size1v;
   const real *Ptr_DI0    = h_Che_Array + CheIdx_DI   *Size1v;
   const real *Ptr_DII0   = h_Che_Array + CheIdx_DII  *Size1v;
   const real *Ptr_HDI0   = h_Che_Array + CheIdx_HDI  *Size1v;


#  pragma omp parallel
   {

// thread-private variables
   int  idx_p, idx_pg, PID, PID0, offset;    // idx_p/idx_pg: array indices within a patch/patch group
   real Dens, Eint;
#  ifdef DUAL_ENERGY
   real Pres;
#  endif
   real (*fluid)[PS1][PS1][PS1]=NULL;

   const real *Ptr_Dens=NULL, *Ptr_sEint=NULL, *Ptr_Ent=NULL, *Ptr_e=NULL, *Ptr_HI=NULL, *Ptr_HII=NULL;
   const real *Ptr_HeI=NULL, *Ptr_HeII=NULL, *Ptr_HeIII=NULL, *Ptr_HM=NULL, *Ptr_H2I=NULL, *Ptr_H2II=NULL;
   const real *Ptr_DI=NULL, *Ptr_DII=NULL, *Ptr_HDI=NULL;

#  pragma omp for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0      = PID0_List[TID];
      idx_pg    = 0;
      offset    = TID*Size1pg;

      Ptr_Dens  = Ptr_Dens0  + offset;
      Ptr_sEint = Ptr_sEint0 + offset;
      Ptr_Ent   = Ptr_Ent0   + offset;
      Ptr_e     = Ptr_e0     + offset;
      Ptr_HI    = Ptr_HI0    + offset;
      Ptr_HII   = Ptr_HII0   + offset;
      Ptr_HeI   = Ptr_HeI0   + offset;
      Ptr_HeII  = Ptr_HeII0  + offset;
      Ptr_HeIII = Ptr_HeIII0 + offset;
      Ptr_HM    = Ptr_HM0    + offset;
      Ptr_H2I   = Ptr_H2I0   + offset;
      Ptr_H2II  = Ptr_H2II0  + offset;
      Ptr_DI    = Ptr_DI0    + offset;
      Ptr_DII   = Ptr_DII0   + offset;
      Ptr_HDI   = Ptr_HDI0   + offset;

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID   = PID0 + LocalID;
         idx_p = 0;
         fluid = amr->patch[SaveSg][lv][PID]->fluid;

         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
//          apply internal energy floor
            Dens = Ptr_Dens [idx_pg];
            Eint = Ptr_sEint[idx_pg]*Dens;
            Eint = Hydro_CheckMinEint( Eint, MIN_EINT );

//          update the total energy density
            *( fluid[ENGY     ][0][0] + idx_p ) = Eint + Ptr_Ent[idx_pg];

//          update the dual-energy variable to be consistent with the updated pressure
#           ifdef DUAL_ENERGY
#           if   ( DUAL_ENERGY == DE_ENPY )
//          DE_ENPY only works with EOS_GAMMA, which does not involve passive scalars
            Pres = EoS_DensEint2Pres_CPUPtr( Dens, Eint, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
            *( fluid[ENPY     ][0][0] + idx_p ) = Hydro_DensPres2Entropy( Dens, Pres, EoS_AuxArray_Flt[1] );

#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif
#           endif // #ifdef DUAL_ENERGY

//          update all chemical species
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
            *( fluid[Idx_e    ][0][0] + idx_p ) = Ptr_e    [idx_pg] * MassRatio_ep;
            *( fluid[Idx_HI   ][0][0] + idx_p ) = Ptr_HI   [idx_pg];
            *( fluid[Idx_HII  ][0][0] + idx_p ) = Ptr_HII  [idx_pg];
            *( fluid[Idx_HeI  ][0][0] + idx_p ) = Ptr_HeI  [idx_pg];
            *( fluid[Idx_HeII ][0][0] + idx_p ) = Ptr_HeII [idx_pg];
            *( fluid[Idx_HeIII][0][0] + idx_p ) = Ptr_HeIII[idx_pg];
            }

//          9-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
            *( fluid[Idx_HM   ][0][0] + idx_p ) = Ptr_HM   [idx_pg];
            *( fluid[Idx_H2I  ][0][0] + idx_p ) = Ptr_H2I  [idx_pg];
            *( fluid[Idx_H2II ][0][0] + idx_p ) = Ptr_H2II [idx_pg];
            }

//          12-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
            *( fluid[Idx_DI   ][0][0] + idx_p ) = Ptr_DI   [idx_pg];
            *( fluid[Idx_DII  ][0][0] + idx_p ) = Ptr_DII  [idx_pg];
            *( fluid[Idx_HDI  ][0][0] + idx_p ) = Ptr_HDI  [idx_pg];
            }

            idx_p  ++;
            idx_pg ++;
         } // i,j,k
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

   } // end of OpenMP parallel region

} // FUNCTION : Grackle_Close



#endif // #ifdef SUPPORT_GRACKLE
