#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Init_FieldData
// Description :  Initialize the "Che_FieldData" grackle_field_data object of Grackle
//
// Note        :  1. Invoked by Grackle_Init()
//                   --> "Che_FieldData" is freed by End_MemFree()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Grackle_Init_FieldData()
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


// allocate memory
   Che_FieldData = new grackle_field_data;


// initialization
   const int NDim = 3;

// fields not evolving with time
   Che_FieldData->grid_rank               = NDim;
   Che_FieldData->grid_dimension          = new int [NDim];
   Che_FieldData->grid_start              = new int [NDim];
   Che_FieldData->grid_end                = new int [NDim];

// grid_dimension, grid_start, and grid_end are set by CPU_GrackleSolver() since the number
// of patch groups advanced at a time is not a constant
   /*
   for (int d=0; d<NDim; d++)
   {
      Che_FieldData->grid_dimension[d]    = PS2;
      Che_FieldData->grid_start    [d]    = 0;
      Che_FieldData->grid_end      [d]    = PS2 - 1;
   }
   */

// fields set by Grackle_Prepare() during each time-step
   Che_FieldData->density                 = NULL;
   Che_FieldData->internal_energy         = NULL;
   Che_FieldData->grid_dx                 = NULL_REAL;
   Che_FieldData->e_density               = NULL;
   Che_FieldData->HI_density              = NULL;
   Che_FieldData->HII_density             = NULL;
   Che_FieldData->HeI_density             = NULL;
   Che_FieldData->HeII_density            = NULL;
   Che_FieldData->HeIII_density           = NULL;
   Che_FieldData->HM_density              = NULL;
   Che_FieldData->H2I_density             = NULL;
   Che_FieldData->H2II_density            = NULL;
   Che_FieldData->DI_density              = NULL;
   Che_FieldData->DII_density             = NULL;
   Che_FieldData->HDI_density             = NULL;
   Che_FieldData->metal_density           = NULL;

// fields not supported yet
   Che_FieldData->x_velocity              = NULL;
   Che_FieldData->y_velocity              = NULL;
   Che_FieldData->z_velocity              = NULL;
   Che_FieldData->volumetric_heating_rate = NULL;
   Che_FieldData->specific_heating_rate   = NULL;
   Che_FieldData->RT_HI_ionization_rate   = NULL;
   Che_FieldData->RT_HeI_ionization_rate  = NULL;
   Che_FieldData->RT_HeII_ionization_rate = NULL;
   Che_FieldData->RT_H2_dissociation_rate = NULL;
   Che_FieldData->RT_heating_rate         = NULL;

} // FUNCTION : Grackle_Init_FieldData



#endif // #ifdef SUPPORT_GRACKLE
