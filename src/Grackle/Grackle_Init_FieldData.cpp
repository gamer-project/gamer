#include "Copyright.h"
#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Init_FieldData
// Description :  Initialize the "Che_FieldData" grackle_field_data objects of Grackle
//
// Note        :  1. Must be called AFTER Init_Load_Parameter()
//                2. Invoked by Init_MemAllocate()
//                   --> "Che_FieldData" is freed by End_MemFree()
//                3. Useless if GRACKLE_MODE == GRACKLE_MODE_GAMER
//
// Parameter   :  Che_NPG : Number of patch groups to be evaluated at a time
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Grackle_Init_FieldData( const int Che_NPG )
{

// nothing to do if we are not using the original Grackle
   if ( GRACKLE_MODE != GRACKLE_MODE_ORI )   return;


// allocate memory
   const int Che_NP = 8*Che_NPG;

   Che_FieldData = new grackle_field_data [Che_NP];


// initialization
   const int NDim = 3;

   for (int t=0; t<Che_NP; t++)
   {
//    fields not evolving with time
      Che_FieldData[t].grid_rank               = NDim;
      Che_FieldData[t].grid_dimension          = new int [NDim];
      Che_FieldData[t].grid_start              = new int [NDim];
      Che_FieldData[t].grid_end                = new int [NDim];

      for (int d=0; d<NDim; d++)
      {
         Che_FieldData[t].grid_dimension[d]    = PS1;
         Che_FieldData[t].grid_start    [d]    = 0;
         Che_FieldData[t].grid_end      [d]    = PS1 - 1;
      }

//    fields set by Grackle_Prepare() during each time-step
      Che_FieldData[t].density                 = NULL;
      Che_FieldData[t].internal_energy         = NULL;
      Che_FieldData[t].grid_dx                 = NULL_REAL;

//    fields not supported yet
      Che_FieldData[t].x_velocity              = NULL;
      Che_FieldData[t].y_velocity              = NULL;
      Che_FieldData[t].z_velocity              = NULL;
      Che_FieldData[t].metal_density           = NULL;

//    fields not supported yet
      Che_FieldData[t].HI_density              = NULL;
      Che_FieldData[t].HII_density             = NULL;
      Che_FieldData[t].HM_density              = NULL;
      Che_FieldData[t].HeI_density             = NULL;
      Che_FieldData[t].HeII_density            = NULL;
      Che_FieldData[t].HeIII_density           = NULL;
      Che_FieldData[t].H2I_density             = NULL;
      Che_FieldData[t].H2II_density            = NULL;
      Che_FieldData[t].DI_density              = NULL;
      Che_FieldData[t].DII_density             = NULL;
      Che_FieldData[t].HDI_density             = NULL;
      Che_FieldData[t].e_density               = NULL;
      Che_FieldData[t].volumetric_heating_rate = NULL;
      Che_FieldData[t].specific_heating_rate   = NULL;
      Che_FieldData[t].RT_HI_ionization_rate   = NULL;
      Che_FieldData[t].RT_HeI_ionization_rate  = NULL;
      Che_FieldData[t].RT_HeII_ionization_rate = NULL;
      Che_FieldData[t].RT_H2_dissociation_rate = NULL;
      Che_FieldData[t].RT_heating_rate         = NULL;
   } // for (int t=0; t<Che_NP; t++)

} // FUNCTION : Grackle_Init_FieldData



#endif // #ifdef SUPPORT_GRACKLE
