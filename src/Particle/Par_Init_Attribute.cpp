#include "GAMER.h"

#ifdef PARTICLE


// declare as static so that other functions cannot invoke it directly and must use the function pointer
void Par_Init_Attribute_User_Template();

// this function pointer must be set by a test problem initializer
void (*Par_Init_Attribute_User_Ptr)() = NULL;


static int NDefinedAttFlt;    // total number of defined floating-point attributes
static int NDefinedAttInt;    // total number of defined integer        attributes



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Attribute
// Description :  Set all particle attributes (including both built-in and user-defined attributes)
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Total number of attributes is determined by PAR_NATT_FLT/INT_TOTAL = PAR_NATT_FLT/INT_BUILTIN + PAR_NATT_FLT/INT_USER
//                3. To initialize user-defined particle attributes, the function pointer "Par_Init_Attribute_User_Ptr"
//                   must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_Attribute()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 0. initialize counter as zero
   NDefinedAttFlt = 0;
   NDefinedAttInt = 0;


// 1. add built-in intrinsic attributes
//    --> must not change the following order of declaration since they must be consistent
//        with the symbolic constants defined in Macro.h (e.g., PAR_MASS)
   Idx_ParMass = AddParticleAttributeFlt( "ParMass" );
   Idx_ParPosX = AddParticleAttributeFlt( "ParPosX" );
   Idx_ParPosY = AddParticleAttributeFlt( "ParPosY" );
   Idx_ParPosZ = AddParticleAttributeFlt( "ParPosZ" );
   Idx_ParVelX = AddParticleAttributeFlt( "ParVelX" );
   Idx_ParVelY = AddParticleAttributeFlt( "ParVelY" );
   Idx_ParVelZ = AddParticleAttributeFlt( "ParVelZ" );

   Idx_ParType = AddParticleAttributeInt( "ParType" );


// 2. add other built-in attributes
#  ifdef STAR_FORMATION
   Idx_ParCreTime = AddParticleAttributeFlt( "ParCreTime" );
#  endif


// 3. add user-defined attributes
   if ( Par_Init_Attribute_User_Ptr != NULL )   Par_Init_Attribute_User_Ptr();


// 4. must put built-in attributes not to be stored on disk at the END of the attribute list
//    --> make it easier to discard them when storing data on disk (see Output_DumpData_Total(_HDF5).cpp)
//    --> must also be consistent with the symbolic constant (e.g., PAR_TIME and PAR_ACC*) defined in Macro.h
//    --> total number of attributes not to be stored on disk is set by PAR_NATT_FLT/INT_UNSTORED
//        --> currently including time and acceleration*3
#  ifdef STORE_PAR_ACC
   Idx_ParAccX = AddParticleAttributeFlt( "ParAccX" );
   Idx_ParAccY = AddParticleAttributeFlt( "ParAccY" );
   Idx_ParAccZ = AddParticleAttributeFlt( "ParAccZ" );
#  endif
   Idx_ParTime = AddParticleAttributeFlt( "ParTime" );



// 5. validate if all attributes have been set properly
   if ( NDefinedAttFlt != PAR_NATT_FLT_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined floating-point attributes (%d) != expectation (%d) !!\n"
                 "        --> Modify PAR_NATT_FLT_USER in the Makefile or invoke AddParticleAttributeFlt() properly\n",
                 NDefinedAttFlt, PAR_NATT_FLT_TOTAL );

   if ( NDefinedAttInt != PAR_NATT_INT_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined integer attributes (%d) != expectation (%d) !!\n"
                 "        --> Modify PAR_NATT_INT_USER in the Makefile or invoke AddParticleAttributeInt() properly\n",
                 NDefinedAttInt, PAR_NATT_INT_TOTAL );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_Attribute



//-------------------------------------------------------------------------------------------------------
// Function    :  AddParticleAttributeFlt
// Description :  Add a new particle floating-point attribute to the attribute list
//
// Note        :  1. This function will
//                   (1) set the attribute label, which will be used as the output name of the attribute
//                   (2) return the index of the new attribute, which can be used to access the attribute
//                       data (e.g., amr->Par->AttributeFlt[])
//                2. One must invoke AddParticleAttributeFlt() exactly PAR_NATT_FLT_TOTAL times to set the labels
//                   of all attributes
//                3. Invoked by Par_Init_Attribute() and various test problem initializers
//
// Parameter   :  InputLabel : Label (i.e., name) of the new attribute
//
// Return      :  (1) ParAttFltLabel[]
//                (2) Index of the newly added attribute
//-------------------------------------------------------------------------------------------------------
FieldIdx_t AddParticleAttributeFlt( const char *InputLabel )
{

   const FieldIdx_t AttIdx = NDefinedAttFlt ++;


// check
   if ( InputLabel == NULL )
      Aux_Error( ERROR_INFO, "InputLabel == NULL !!\n" );

   if ( NDefinedAttFlt > PAR_NATT_FLT_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined particle floating-point attributes (%d) exceeds expectation (%d) after adding \"%s\" !!\n"
                 "        --> Modify PAR_NATT_FLT_USER in the Makefile properly\n",
                 NDefinedAttFlt, PAR_NATT_FLT_TOTAL, InputLabel );

   for (int v=0; v<NDefinedAttFlt-1; v++)
      if (  strcmp( ParAttFltLabel[v], InputLabel ) == 0  )
         Aux_Error( ERROR_INFO, "duplicate particle attribute label \"%s\" !!\n", InputLabel );


// set attribute label
   strcpy( ParAttFltLabel[AttIdx], InputLabel );


// return attribute index
   return AttIdx;

} // FUNCTION : AddParticleAttributeFlt



//-------------------------------------------------------------------------------------------------------
// Function    :  AddParticleAttributeInt
// Description :  Add a new particle integer attribute to the attribute list
//
// Note        :  1. This function will
//                   (1) set the attribute label, which will be used as the output name of the attribute
//                   (2) return the index of the new attribute, which can be used to access the attribute
//                       data (e.g., amr->Par->AttributeInt[])
//                2. One must invoke AddParticleAttributeInt() exactly PAR_NATT_INT_TOTAL times to set the labels
//                   of all attributes
//                3. Invoked by Par_Init_Attribute() and various test problem initializers
//
// Parameter   :  InputLabel : Label (i.e., name) of the new attribute
//
// Return      :  (1) ParAttIntLabel[]
//                (2) Index of the newly added attribute
//-------------------------------------------------------------------------------------------------------
FieldIdx_t AddParticleAttributeInt( const char *InputLabel )
{

   const FieldIdx_t AttIdx = NDefinedAttInt ++;


// check
   if ( InputLabel == NULL )
      Aux_Error( ERROR_INFO, "InputLabel == NULL !!\n" );

   if ( NDefinedAttInt > PAR_NATT_INT_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined particle integer attributes (%d) exceeds expectation (%d) after adding \"%s\" !!\n"
                 "        --> Modify PAR_NATT_INT_USER in the Makefile properly\n",
                 NDefinedAttInt, PAR_NATT_INT_TOTAL, InputLabel );

   for (int v=0; v<NDefinedAttInt-1; v++)
      if (  strcmp( ParAttIntLabel[v], InputLabel ) == 0  )
         Aux_Error( ERROR_INFO, "duplicate particle attribute label \"%s\" !!\n", InputLabel );


// set attribute label
   strcpy( ParAttIntLabel[AttIdx], InputLabel );


// return attribute index
   return AttIdx;

} // FUNCTION : AddParticleAttributeInt



//-------------------------------------------------------------------------------------------------------
// Function    :  GetParticleAttributeFltIndex
// Description :  Return the index of the target particle floating-point attribute
//
// Note        :  1. Usage: AttIdx = GetParticleAttributeFltIndex( ParAttFltLabel );
//                2. Return Idx_Undefined if the target attribute cannot be found
//
// Parameter   :  InputLabel : Target attribute label
//                Check      : Whether or not to terminate the program if the target attribute cannot be found
//                             --> Accepted options: CHECK_ON / CHECK_OFF
//
// Return      :  Sucess: index of the target attribute
//                Failed: Idx_Undefined
//-------------------------------------------------------------------------------------------------------
FieldIdx_t GetParticleAttributeFltIndex( const char *InputLabel, const Check_t Check )
{

   FieldIdx_t Idx_Out = Idx_Undefined;

   for (int v=0; v<NDefinedAttFlt; v++)
   {
      if (  strcmp( ParAttFltLabel[v], InputLabel ) == 0  )
      {
         Idx_Out = v;
         break;
      }
   }

   if ( Check == CHECK_ON  &&  Idx_Out == Idx_Undefined )
      Aux_Error( ERROR_INFO, "cannot find the target particle floating-point attribute \"%s\" !!\n", InputLabel );

   return Idx_Out;

} // FUNCTION : GetParticleAttributeFltIndex



//-------------------------------------------------------------------------------------------------------
// Function    :  GetParticleAttributeIntIndex
// Description :  Return the index of the target particle integer attribute
//
// Note        :  1. Usage: AttIdx = GetParticleAttributeIntIndex( ParAttIntLabel );
//                2. Return Idx_Undefined if the target attribute cannot be found
//
// Parameter   :  InputLabel : Target attribute label
//                Check      : Whether or not to terminate the program if the target attribute cannot be found
//                             --> Accepted options: CHECK_ON / CHECK_OFF
//
// Return      :  Sucess: index of the target attribute
//                Failed: Idx_Undefined
//-------------------------------------------------------------------------------------------------------
FieldIdx_t GetParticleAttributeIntIndex( const char *InputLabel, const Check_t Check )
{

   FieldIdx_t Idx_Out = Idx_Undefined;

   for (int v=0; v<NDefinedAttInt; v++)
   {
      if (  strcmp( ParAttIntLabel[v], InputLabel ) == 0  )
      {
         Idx_Out = v;
         break;
      }
   }

   if ( Check == CHECK_ON  &&  Idx_Out == Idx_Undefined )
      Aux_Error( ERROR_INFO, "cannot find the target particle integer attribute \"%s\" !!\n", InputLabel );

   return Idx_Out;

} // FUNCTION : GetParticleAttributeIntIndex



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Attribute_User_Template
// Description :  Template of adding user-defined particle attributes
//
// Note        :  1. Invoked by Par_Init_Attribute() using the function pointer "Par_Init_Attribute_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_Attribute_User_Template()
{

// example
// Idx_NewAttFlt = AddParticleAttributeFlt( "NewAttFltLabel" );
// Idx_NewAttInt = AddParticleAttributeInt( "NewAttIntLabel" );

} // FUNCTION : Par_Init_Attribute_User_Template



#endif // #ifdef PARTICLE
