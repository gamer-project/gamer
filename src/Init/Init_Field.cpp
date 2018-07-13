#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
void Init_Field_User();

// this function pointer may be overwritten by various test problem initializers
void (*Init_Field_User_Ptr)() = Init_Field_User;


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
//    --> must not change the following order of declaration since they must be consistent
//        with the symbolic constants defined in Macro.h (e.g., DENS)
#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   Idx_Dens    = AddField( "Dens" );
   Idx_MomX    = AddField( "MomX" );
   Idx_MomY    = AddField( "MomY" );
   Idx_MomZ    = AddField( "MomZ" );
   Idx_Engy    = AddField( "Engy" );

#  elif ( MODEL == ELBDM )

#  else
#    error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 2. add other predefined fields
//    --> must declare the dual-energy variable first in order to be consistent with the symbolic
//        constant ENPY (or EINT) defined in Macro.h
#  if   ( DUAL_ENERGY == DE_ENPY )
   Idx_Enpy    = AddField( "Entropy" );
#  elif ( DUAL_ENERGY == DE_EINT )
   Idx_Eint    = AddField( "Eint" );
#  endif

#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_METAL )
   Idx_Metal   = AddField( "Metal" );
#  endif


// 3. add user-defined fields
   if ( Init_Field_User_Ptr != NULL )  Init_Field_User_Ptr();


// 4. validate if all fields have been set properly
   if ( NDefinedField != NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "Total number of defined fields (%d) != expectation (%d) !!\n"
                 "        --> Modify NCOMP_PASSIVE_USER in the Makefile or invoke AddField() properly\n",
                 NDefinedField, NCOMP_TOTAL );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Field



//-------------------------------------------------------------------------------------------------------
// Function    :  AddField
// Description :  Add a new field to the field list
//
// Note        :  1. This function will (i) set the label and (ii) return the index of the newly added field
//                   --> Field label will be used as the output name of the field
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
      Aux_Error( ERROR_INFO, "Total number of defined fields (%d) exceeds expectation (%d) after adding the field \"%s\" !!\n"
                 "        --> Modify NCOMP_PASSIVE_USER in the Makefile properly\n",
                 NDefinedField+1, NCOMP_TOTAL, InputLabel );

// set field label
   FieldLabel[NDefinedField] = InputLabel;

// return field index
   return NDefinedField ++;

} // FUNCTION : AddField




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Field_User
// Description :  Add user-defined field
//
// Note        :  1. Invoked by Init_Field() using the function pointer "Init_Field_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_Field_User()
{

// example
// Idx_NewField = AddField( "NewFieldLabel" );

} // FUNCTION : Init_Field_User
