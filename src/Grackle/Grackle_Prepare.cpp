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
extern int CheIdx_vHeatingRate;
extern int CheIdx_sHeatingRate;


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static real_che Grackle_vHeatingRate_User_Template( const double x, const double y, const double z, const double Time, const double n_H );
static real_che Grackle_sHeatingRate_User_Template( const double x, const double y, const double z, const double Time );


// these function pointers must be set by a test problem initializer
real_che (*Grackle_vHeatingRate_User_Ptr)( const double x, const double y, const double z, const double Time, const double n_H ) = NULL;
real_che (*Grackle_sHeatingRate_User_Ptr)( const double x, const double y, const double z, const double Time )                   = NULL;



//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Prepare
// Description :  Fill up the input host array h_Che_Array[] for the Grackle solver
//
// Note        :  1. Prepare Che_NField variables
//                   --> Che_NField and the corresponding array indices in h_Che_Array[] (e.g., CheIdx_Dens)
//                       are declared and set by Init_MemAllocate_Grackle()
//                2. This function always prepares the latest FluSg data
//
// Parameter   :  lv          : Target refinement level
//                h_Che_Array : Host array to store the prepared data
//                NPG         : Number of patch groups prepared at a time
//                PID0_List   : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Grackle_Prepare( const int lv, real_che h_Che_Array[], const int NPG, const int *PID0_List )
{

// check
#  ifdef GAMER_DEBUG
   if ( CheIdx_Dens == Idx_Undefined )
      Aux_Error( ERROR_INFO, "CheIdx_Dens is undefined !!\n" );
   if ( CheIdx_sEint == Idx_Undefined )
      Aux_Error( ERROR_INFO, "CheIdx_sEint is undefined !!\n" );
   if ( CheIdx_Ent == Idx_Undefined )
      Aux_Error( ERROR_INFO, "CheIdx_Ent is undefined !!\n" );

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
      if (  Idx_e == Idx_Undefined  ||  CheIdx_e == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_e is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HI == Idx_Undefined  ||  CheIdx_HI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HI is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HII == Idx_Undefined  ||  CheIdx_HII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HII is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HeI == Idx_Undefined  ||  CheIdx_HeI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HeI is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HeII == Idx_Undefined  ||  CheIdx_HeII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HeII is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HeIII == Idx_Undefined  ||  CheIdx_HeIII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HeIII is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
      if (  Idx_HM == Idx_Undefined  ||  CheIdx_HM == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HM is undefined for \"GRACKLE_PRI_CHE_NSPE9\" !!\n" );
      if (  Idx_H2I == Idx_Undefined  ||  CheIdx_H2I == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_H2I is undefined for \"GRACKLE_PRI_CHE_NSPE9\" !!\n" );
      if (  Idx_H2II == Idx_Undefined  ||  CheIdx_H2II == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_H2II is undefined for \"GRACKLE_PRI_CHE_NSPE9\" !!\n" );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
      if (  Idx_DI == Idx_Undefined  ||  CheIdx_DI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_DI is undefined for \"GRACKLE_PRI_CHE_NSPE12\" !!\n" );
      if (  Idx_DII == Idx_Undefined  ||  CheIdx_DII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_DII is undefined for \"GRACKLE_PRI_CHE_NSPE12\" !!\n" );
      if (  Idx_HDI == Idx_Undefined  ||  CheIdx_HDI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HDI is undefined for \"GRACKLE_PRI_CHE_NSPE12\" !!\n" );
   }

   if ( GRACKLE_METAL ) {
      if (  Idx_Metal == Idx_Undefined  ||  CheIdx_Metal == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_Metal is undefined for \"GRACKLE_METAL\" !!\n" );
   }

   if ( GRACKLE_USE_V_HEATING_RATE ) {
      if ( Grackle_vHeatingRate_User_Ptr == NULL )
         Aux_Error( ERROR_INFO, "Grackle_vHeatingRate_User_Ptr == NULL !!\n" );
      if ( CheIdx_vHeatingRate == Idx_Undefined )
         Aux_Error( ERROR_INFO, "CheIdx_vHeatingRate is undefined for \"GRACKLE_USE_V_HEATING_RATE\" !!\n" );
   }

   if ( GRACKLE_USE_S_HEATING_RATE ) {
      if ( Grackle_sHeatingRate_User_Ptr == NULL )
         Aux_Error( ERROR_INFO, "Grackle_sHeatingRate_User_Ptr == NULL !!\n" );
      if ( CheIdx_sHeatingRate == Idx_Undefined )
         Aux_Error( ERROR_INFO, "CheIdx_sHeatingRate is undefined for \"GRACKLE_USE_S_HEATING_RATE\" !!\n" );
   }
#  endif // #ifdef GAMER_DEBUG


   const int  Size1pg          = CUBE(PS2);
   const int  Size1v           = NPG*Size1pg;
   const real MassRatio_pe    = Const_mp / Const_me;
#  ifdef DUAL_ENERGY
   const bool CheckMinPres_No  = false;
#  else
   const bool CheckMinEint_Yes = true;
#  endif

   const double dh = amr->dh[lv];

// Floor value for the mass density of each species applied by Grackle
// --> See ceiling_species_g() in grackle/src/clib/solve_rate_cool_g.F
//     and the definition of tiny in grackle/src/include/grackle_fortran_types.def
// --> Note this value is a quantity with **code density units UNIT_D**
   const real_che Grackle_tiny = 1.0e-20;

// Minimum of total density for the Grackle solver
// --> When the total density is too low (close to Grackle_tiny),
//     the species floor value in Grackle will be a significant mass fraction
// --> The factor here is decided such that the species with low mass fraction
//     remains negligible after being applied the Grackle_tiny floor
// --> Note that to maintain accuracy, this should have a physically small density
   const real_che Che_MinDens  = 1.0e+04 * Grackle_tiny;

   real_che *Ptr_Dens0         = h_Che_Array + CheIdx_Dens        *Size1v;
   real_che *Ptr_sEint0        = h_Che_Array + CheIdx_sEint       *Size1v;
   real_che *Ptr_Ent0          = h_Che_Array + CheIdx_Ent         *Size1v;
   real_che *Ptr_e0            = h_Che_Array + CheIdx_e           *Size1v;
   real_che *Ptr_HI0           = h_Che_Array + CheIdx_HI          *Size1v;
   real_che *Ptr_HII0          = h_Che_Array + CheIdx_HII         *Size1v;
   real_che *Ptr_HeI0          = h_Che_Array + CheIdx_HeI         *Size1v;
   real_che *Ptr_HeII0         = h_Che_Array + CheIdx_HeII        *Size1v;
   real_che *Ptr_HeIII0        = h_Che_Array + CheIdx_HeIII       *Size1v;
   real_che *Ptr_HM0           = h_Che_Array + CheIdx_HM          *Size1v;
   real_che *Ptr_H2I0          = h_Che_Array + CheIdx_H2I         *Size1v;
   real_che *Ptr_H2II0         = h_Che_Array + CheIdx_H2II        *Size1v;
   real_che *Ptr_DI0           = h_Che_Array + CheIdx_DI          *Size1v;
   real_che *Ptr_DII0          = h_Che_Array + CheIdx_DII         *Size1v;
   real_che *Ptr_HDI0          = h_Che_Array + CheIdx_HDI         *Size1v;
   real_che *Ptr_Metal0        = h_Che_Array + CheIdx_Metal       *Size1v;
   real_che *Ptr_vHeatingRate0 = h_Che_Array + CheIdx_vHeatingRate*Size1v;
   real_che *Ptr_sHeatingRate0 = h_Che_Array + CheIdx_sHeatingRate*Size1v;


#  pragma omp parallel
   {

// thread-private variables
   int  idx_p, idx_pg, PID, PID0, offset;    // idx_p/idx_pg: array indices within a patch/patch group
   real Dens, Etot, Eint, Ent;
#  ifdef DUAL_ENERGY
   real Pres;
#  else
   real Px, Py, Pz, Emag=NULL_REAL;
#  endif // #ifdef DUAL_ENERGY ... else ...
   real (*fluid)[PS1][PS1][PS1]=NULL;

   real_che *Ptr_Dens=NULL, *Ptr_sEint=NULL, *Ptr_Ent=NULL, *Ptr_e=NULL, *Ptr_HI=NULL, *Ptr_HII=NULL;
   real_che *Ptr_HeI=NULL, *Ptr_HeII=NULL, *Ptr_HeIII=NULL, *Ptr_HM=NULL, *Ptr_H2I=NULL, *Ptr_H2II=NULL;
   real_che *Ptr_DI=NULL, *Ptr_DII=NULL, *Ptr_HDI=NULL, *Ptr_Metal=NULL;
   real_che *Ptr_vHeatingRate=NULL, *Ptr_sHeatingRate=NULL;
   real_che  Ratio_FloorDens;

#  pragma omp for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0   = PID0_List[TID];
      idx_pg = 0;
      offset = TID*Size1pg;

      Ptr_Dens         = Ptr_Dens0         + offset;
      Ptr_sEint        = Ptr_sEint0        + offset;
      Ptr_Ent          = Ptr_Ent0          + offset;
      Ptr_e            = Ptr_e0            + offset;
      Ptr_HI           = Ptr_HI0           + offset;
      Ptr_HII          = Ptr_HII0          + offset;
      Ptr_HeI          = Ptr_HeI0          + offset;
      Ptr_HeII         = Ptr_HeII0         + offset;
      Ptr_HeIII        = Ptr_HeIII0        + offset;
      Ptr_HM           = Ptr_HM0           + offset;
      Ptr_H2I          = Ptr_H2I0          + offset;
      Ptr_H2II         = Ptr_H2II0         + offset;
      Ptr_DI           = Ptr_DI0           + offset;
      Ptr_DII          = Ptr_DII0          + offset;
      Ptr_HDI          = Ptr_HDI0          + offset;
      Ptr_Metal        = Ptr_Metal0        + offset;
      Ptr_vHeatingRate = Ptr_vHeatingRate0 + offset;
      Ptr_sHeatingRate = Ptr_sHeatingRate0 + offset;

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID   = PID0 + LocalID;
         idx_p = 0;
         fluid = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;

         const double x0  = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0  = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0  = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            Dens  = *( fluid[DENS][0][0] + idx_p );
            Etot  = *( fluid[ENGY][0][0] + idx_p );

//          use the dual-energy variable to calculate the internal energy if applicable
#           ifdef DUAL_ENERGY

#           if   ( DUAL_ENERGY == DE_ENPY )
            Pres  = Hydro_DensDual2Pres( Dens, *(fluid[DUAL][0][0]+idx_p), EoS_AuxArray_Flt[1], CheckMinPres_No, NULL_REAL );
//          EOS_GAMMA does not involve passive scalars
            Eint  = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

#           else // #ifdef DUAL_ENERGY

            Px    = *( fluid[MOMX][0][0] + idx_p );
            Py    = *( fluid[MOMY][0][0] + idx_p );
            Pz    = *( fluid[MOMZ][0][0] + idx_p );
#           ifdef MHD
            Emag  = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#           endif
            Eint  = Hydro_Con2Eint( Dens, Px, Py, Pz, Etot, CheckMinEint_Yes, MIN_EINT, PassiveFloorMask, Emag,
                                    NULL, NULL, NULL, NULL, NULL );
#           endif // #ifdef DUAL_ENERGY ... else

//          Grackle doesn't know cosmic rays so we must exclude the cosmic-ray energy from the input gas internal energy
#           ifdef COSMIC_RAY
            Eint -= *( fluid[CRAY][0][0] + idx_p );
#           endif

//          set the ratio to the density floor
            if ( Dens < Che_MinDens )
            {
//             when the density is too low, we will pre-floor the total density as Che_MinDens
//             while keeping the fraction of each species fixed
//             such that the mass fraction of the species being floored in Grackle is
//             less than "Grackle_tiny/Che_MinDens"
               Ratio_FloorDens = Che_MinDens/Dens;

               Aux_Message( stderr, "\nWARNING : density = %16.8e is too low and will be scaled to %16.8e as the input to the Grackle solver !!\n",
                                    Dens, Che_MinDens );
            }
            else
            {
               Ratio_FloorDens = 1.0;
            }

//          mandatory fields
            Ptr_Dens [idx_pg] = Dens * Ratio_FloorDens;
            Ptr_sEint[idx_pg] = Eint / Dens;
            Ptr_Ent  [idx_pg] = (Etot - Eint) * Ratio_FloorDens; // non-thermal energy density

//          6-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
            Ptr_e    [idx_pg] = *( fluid[Idx_e    ][0][0] + idx_p ) * Ratio_FloorDens * MassRatio_pe;
            Ptr_HI   [idx_pg] = *( fluid[Idx_HI   ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_HII  [idx_pg] = *( fluid[Idx_HII  ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_HeI  [idx_pg] = *( fluid[Idx_HeI  ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_HeII [idx_pg] = *( fluid[Idx_HeII ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_HeIII[idx_pg] = *( fluid[Idx_HeIII][0][0] + idx_p ) * Ratio_FloorDens;
            }

//          9-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
            Ptr_HM   [idx_pg] = *( fluid[Idx_HM   ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_H2I  [idx_pg] = *( fluid[Idx_H2I  ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_H2II [idx_pg] = *( fluid[Idx_H2II ][0][0] + idx_p ) * Ratio_FloorDens;
            }

//          12-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
            Ptr_DI   [idx_pg] = *( fluid[Idx_DI   ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_DII  [idx_pg] = *( fluid[Idx_DII  ][0][0] + idx_p ) * Ratio_FloorDens;
            Ptr_HDI  [idx_pg] = *( fluid[Idx_HDI  ][0][0] + idx_p ) * Ratio_FloorDens;
            }

//          metallicity for metal cooling
            if ( GRACKLE_METAL )
            Ptr_Metal[idx_pg] = *( fluid[Idx_Metal][0][0] + idx_p ) * Ratio_FloorDens;

//          user-provided array of volumetric heating rates
            if ( GRACKLE_USE_V_HEATING_RATE ) {
            const double n_H = Ptr_Dens[idx_pg] * UNIT_D * GRACKLE_HYDROGEN_MFRAC / Const_mH; // hydrogen number density in units of cm^-3
            Ptr_vHeatingRate[idx_pg] = Grackle_vHeatingRate_User_Ptr( x0+i*dh, y0+j*dh, z0+k*dh, Time[lv], n_H );
            }

//          user-provided array of specific heating rates
            if ( GRACKLE_USE_S_HEATING_RATE )
            Ptr_sHeatingRate[idx_pg] = Grackle_sHeatingRate_User_Ptr( x0+i*dh, y0+j*dh, z0+k*dh, Time[lv] );

            idx_p  ++;
            idx_pg ++;
         } // i,j,k

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

   } // end of OpenMP parallel region


// set cell size and link pointers for different fields
   Che_FieldData->grid_dx         = amr->dh[lv];

   Che_FieldData->density         = Ptr_Dens0;
   Che_FieldData->internal_energy = Ptr_sEint0;

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   Che_FieldData->e_density       = Ptr_e0;
   Che_FieldData->HI_density      = Ptr_HI0;
   Che_FieldData->HII_density     = Ptr_HII0;
   Che_FieldData->HeI_density     = Ptr_HeI0;
   Che_FieldData->HeII_density    = Ptr_HeII0;
   Che_FieldData->HeIII_density   = Ptr_HeIII0;
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   Che_FieldData->HM_density      = Ptr_HM0;
   Che_FieldData->H2I_density     = Ptr_H2I0;
   Che_FieldData->H2II_density    = Ptr_H2II0;
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   Che_FieldData->DI_density      = Ptr_DI0;
   Che_FieldData->DII_density     = Ptr_DII0;
   Che_FieldData->HDI_density     = Ptr_HDI0;
   }

   if ( GRACKLE_METAL )
   Che_FieldData->metal_density   = Ptr_Metal0;

   if ( GRACKLE_USE_V_HEATING_RATE )
   Che_FieldData->volumetric_heating_rate = Ptr_vHeatingRate0;

   if ( GRACKLE_USE_S_HEATING_RATE )
   Che_FieldData->specific_heating_rate   = Ptr_sHeatingRate0;

} // FUNCTION : Grackle_Prepare



//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_vHeatingRate_User_Template
// Description :  Function template to set Grackle's volumetric heating rate
//
// Note        :  1. Invoked by Grackle_Prepare() using the function pointer
//                   "Grackle_vHeatingRate_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Returned rate should be in the unit of erg s^-1 cm^-3
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                n_H      : Hydrogen number density, should be in the unit of cm^-3
//
// Return      :  volumetric_heating_rate
//-------------------------------------------------------------------------------------------------------
real_che Grackle_vHeatingRate_User_Template( const double x, const double y, const double z, const double Time, const double n_H )
{

   const double   Center[3]                 = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double   radius_to_center          = sqrt( SQR(x-Center[0]) + SQR(y-Center[1]) );
   const double   ScaleLength               = 0.125*amr->BoxSize[0];
   const double   R_core                    = 0.125*amr->BoxSize[0];
   const double   volumetric_heating_rate_0 = n_H * GRACKLE_PE_HEATING_RATE; // GRACKLE_PE_HEATING_RATE is in units of erg cm^-3 s^-1 n_H^-1
   const real_che volumetric_heating_rate   = volumetric_heating_rate_0 * exp( (R_core - MAX(radius_to_center, R_core)) / ScaleLength );

   return volumetric_heating_rate;

} // FUNCTION : Grackle_vHeatingRate_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_sHeatingRate_User_Template
// Description :  Function template to set Grackle's specific heating rate
//
// Note        :  1. Invoked by Grackle_Prepare() using the function pointer
//                   "Grackle_sHeatingRate_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Returned rate should be in the units of erg s^-1 g^-1
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//
// Return      :  specific_heating_rate
//-------------------------------------------------------------------------------------------------------
real_che Grackle_sHeatingRate_User_Template( const double x, const double y, const double z, const double Time )
{

   const double   Center[3]                 = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double   radius_to_center          = sqrt( SQR(x-Center[0]) + SQR(y-Center[1]) );
   const double   ScaleLength               = 0.125*amr->BoxSize[0];
   const double   R_core                    = 0.125*amr->BoxSize[0];
   const double   specific_heating_rate_0   = GRACKLE_PE_HEATING_RATE * GRACKLE_HYDROGEN_MFRAC / Const_mH; // GRACKLE_PE_HEATING_RATE is in units of erg cm^-3 s^-1 n_H^-1
   const real_che specific_heating_rate     = specific_heating_rate_0 * exp( (R_core - MAX(radius_to_center, R_core)) / ScaleLength );

   return specific_heating_rate;

} // FUNCTION : Grackle_sHeatingRate_User_Template



#endif // #ifdef SUPPORT_GRACKLE
