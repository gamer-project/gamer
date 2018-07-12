#include "GAMER.h"

static int NDefinedField = 0;    // total number of defined fields --> for debug only




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Field
// Description :  Set all fields (including both predefined and user-defined fields)
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Total number of fields is determined by NCOMP_TOTAL = NCOMP_FLUID + NCOMP_PASSIVE
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_Field()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. add main predefined fields
//    --> must not change the following order of their declaration
#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   Idx_Dens = AddField( "Dens" );
   Idx_MomX = AddField( "MomX" );
   Idx_MomY = AddField( "MomY" );
   Idx_MomZ = AddField( "MomZ" );
   Idx_Engy = AddField( "Engy" );

#  elif ( MODEL == ELBDM )

#  else
#    error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 2. add other predefined fields
// entropy
// grackle-related field


// 3. add user-defined fields


// 4. validate if all fields have been set properly
   if ( NDefinedField != NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "Total number of defined fields (%d) != expectation (%d) !!\n"
                 "        --> Modify NCOMP_PASSIVE_USER in the Makefile properly\n",
                 NDefinedField, NCOMP_TOTAL );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Field



//-------------------------------------------------------------------------------------------------------
// Function    :  AddField
// Description :  Add a new field to the field list
//
// Note        :  1. This function will (i) set the label and (ii) return the index of the newly added field
//                   --> Field label will be used as the name of the output field
//                   --> Field index can be used to access the field data (e.g., amr->patch->fluid[])
//                2. One must invoke AddField() exactly NCOMP_TOTAL times to set the labels of all fields
//                3. Invoked by Init_Field() and various test problem initializers
//
// Parameter   :  InputLabel : Label (i.e., name) of the new field
//
// Return      :  Set FieldLabel[] and return the index of the newly added field
//-------------------------------------------------------------------------------------------------------
int AddField( const char *InputLabel )
{

// check
   if ( NDefinedField+1 > NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "Total number of defined fields (%d) exceeds expectation (%d) !!\n"
                 "        --> Modify NCOMP_PASSIVE_USER in the Makefile properly\n",
                 NDefinedField+1, NCOMP_TOTAL );

// set field label
   FieldLabel[NDefinedField] = InputLabel;

// return field index
   return NDefinedField ++;

} // FUNCTION : AddField
