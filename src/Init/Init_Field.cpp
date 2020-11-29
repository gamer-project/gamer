#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_Field_User_Template();

// this function pointer must be set by a test problem initializer
void (*Init_Field_User_Ptr)() = NULL;


static int NDefinedField;  // total number of defined fields --> for debug only




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Field
// Description :  Set all fields (including both predefined and user-defined fields)
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Total number of fields is determined by NCOMP_TOTAL = NCOMP_FLUID + NCOMP_PASSIVE
//                3. To initialize user-defined fields, the function pointer "Init_Field_User_Ptr" must be
//                   set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_Field()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 0. initialize counters as zeros
   NDefinedField    = 0;
   PassiveNorm_NVar = 0;


// 1. add main predefined fields
//    --> must not change the following order of declaration since they must be consistent
//        with the symbolic constants defined in Macro.h (e.g., DENS)
#  if   ( MODEL == HYDRO )
   Idx_Dens    = AddField( "Dens",     NORMALIZE_NO );
   Idx_MomX    = AddField( "MomX",     NORMALIZE_NO );
   Idx_MomY    = AddField( "MomY",     NORMALIZE_NO );
   Idx_MomZ    = AddField( "MomZ",     NORMALIZE_NO );
   Idx_Engy    = AddField( "Engy",     NORMALIZE_NO );

   if ( Idx_Dens != DENS )    Aux_Error( ERROR_INFO, "inconsistent Idx_Dens (%d != %d) !!\n", Idx_Dens, DENS );
   if ( Idx_MomX != MOMX )    Aux_Error( ERROR_INFO, "inconsistent Idx_MomX (%d != %d) !!\n", Idx_MomX, MOMX );
   if ( Idx_MomY != MOMY )    Aux_Error( ERROR_INFO, "inconsistent Idx_MomY (%d != %d) !!\n", Idx_MomY, MOMY );
   if ( Idx_MomZ != MOMZ )    Aux_Error( ERROR_INFO, "inconsistent Idx_MomZ (%d != %d) !!\n", Idx_MomZ, MOMZ );
   if ( Idx_Engy != ENGY )    Aux_Error( ERROR_INFO, "inconsistent Idx_Engy (%d != %d) !!\n", Idx_Engy, ENGY );

#  ifdef MHD
   MagLabel[MAGX] = "MagX";
   MagLabel[MAGY] = "MagY";
   MagLabel[MAGZ] = "MagZ";
#  endif

#  elif ( MODEL == ELBDM )

#  else
#    error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 2. add other predefined fields
#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   Idx_e       = AddField( "Electron", NORMALIZE_YES );
   Idx_HI      = AddField( "HI",       NORMALIZE_YES );
   Idx_HII     = AddField( "HII",      NORMALIZE_YES );
   Idx_HeI     = AddField( "HeI",      NORMALIZE_YES );
   Idx_HeII    = AddField( "HeII",     NORMALIZE_YES );
   Idx_HeIII   = AddField( "HeIII",    NORMALIZE_YES );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   Idx_HM      = AddField( "HM",       NORMALIZE_YES );
   Idx_H2I     = AddField( "H2I",      NORMALIZE_YES );
   Idx_H2II    = AddField( "H2II",     NORMALIZE_YES );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   Idx_DI      = AddField( "DI",       NORMALIZE_YES );
   Idx_DII     = AddField( "DII",      NORMALIZE_YES );
   Idx_HDI     = AddField( "HDI",      NORMALIZE_YES );
   }

// normalize the metallicity field only when adopting the non-equilibrium chemistry
// --> may need a machanism to allow users to overwrite this default setup
   if ( GRACKLE_METAL )
   Idx_Metal   = AddField( "Metal",    (GRACKLE_PRIMORDIAL==GRACKLE_PRI_CHE_CLOUDY)?NORMALIZE_NO:NORMALIZE_YES );
#  endif // #ifdef SUPPORT_GRACKLE


// 3. add user-defined fields
   if ( Init_Field_User_Ptr != NULL )  Init_Field_User_Ptr();


// 4. must put all built-in scalars at the END of the field list and with the same order as their
//    corresponding symbolic constants (e.g., ENPY/EINT/CRAY) defined in Macro.h
//    --> as we still rely on these constants (e.g., DENS, ENPY) in the fluid solvers
#  ifdef COSMIC_RAY
   Idx_CRay    = AddField( "CRay",     NORMALIZE_NO );
   if ( Idx_CRay != CRAY )    Aux_Error( ERROR_INFO, "inconsistent Idx_CRay (%d != %d) !!\n", Idx_CRay, CRAY );
#  endif

#  if   ( DUAL_ENERGY == DE_ENPY )
   Idx_Enpy    = AddField( "Entropy",  NORMALIZE_NO );
   if ( Idx_Enpy != ENPY )    Aux_Error( ERROR_INFO, "inconsistent Idx_Enpy (%d != %d) !!\n", Idx_Enpy, ENPY );
#  elif ( DUAL_ENERGY == DE_EINT )
   Idx_Eint    = AddField( "Eint",     NORMALIZE_NO );
   if ( Idx_Eint != EINT )    Aux_Error( ERROR_INFO, "inconsistent Idx_Eint (%d != %d) !!\n", Idx_Eint, EINT );
#  endif



// 5. validate if all fields have been set properly
   if ( NDefinedField != NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined fields (%d) != expectation (%d) !!\n"
                 "        --> Modify NCOMP_PASSIVE_USER in the Makefile or invoke AddField() properly\n",
                 NDefinedField, NCOMP_TOTAL );


// 6. turn off OPT__NORMALIZE_PASSIVE if there are no passive scalars to be normalized
   if ( OPT__NORMALIZE_PASSIVE  &&  PassiveNorm_NVar == 0 )
   {
      OPT__NORMALIZE_PASSIVE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : disabling \"OPT__NORMALIZE_PASSIVE\" since there are no passive scalars to be normalized !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Field



//-------------------------------------------------------------------------------------------------------
// Function    :  AddField
// Description :  Add a new field to the field list
//
// Note        :  1. This function will
//                   (1) set the field label, which will be used as the output name of the field
//                   (2) return the index of the new field, which can be used to access the field data
//                       (e.g., amr->patch->fluid[])
//                   (3) add the new field to the normalization list PassiveNorm_VarIdx[] if requested
//                       --> The sum of passive scalars in this list will be normalized to the gas density
//                       --> sum(passive_scalar_density) == gas_density
//                2. One must invoke AddField() exactly NCOMP_TOTAL times to set the labels of all fields
//                3. Invoked by Init_Field() and various test problem initializers
//
// Parameter   :  InputLabel : Label (i.e., name) of the new field
//                Norm       : whether or not to normalize the new field
//
// Return      :  (1) FieldLabel[]
//                (2) Index of the newly added field
//                (3) PassiveNorm_NVar & PassiveNorm_VarIdx[]
//-------------------------------------------------------------------------------------------------------
FieldIdx_t AddField( char *InputLabel, const NormPassive_t Norm )
{

   const FieldIdx_t FieldIdx = NDefinedField ++;


// check
   if ( InputLabel == NULL )
      Aux_Error( ERROR_INFO, "InputLabel == NULL !!\n" );

   if ( NDefinedField > NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined fields (%d) exceeds expectation (%d) after adding the field \"%s\" !!\n"
                 "        --> Modify NCOMP_PASSIVE_USER in the Makefile properly\n",
                 NDefinedField, NCOMP_TOTAL, InputLabel );

   for (int v=0; v<NDefinedField-1; v++)
      if (  strcmp( FieldLabel[v], InputLabel ) == 0  )
         Aux_Error( ERROR_INFO, "duplicate field label \"%s\" !!\n", InputLabel );


// set field label
   strcpy( FieldLabel[FieldIdx], InputLabel );


// set the normalization list
// --> note that PassiveNorm_VarIdx[] starts from 0 instead of NCOMP_FLUID
// --> currently we set PassiveNorm_VarIdx[] no mater OPT__NORMALIZE_PASSIVE is on or off
//     --> this allows Aux_Check_Conservation() to work on passive scalars even when OPT__NORMALIZE_PASSIVE is off
   if ( Norm )
   {
      const int NormIdx = PassiveNorm_NVar ++;

      if ( PassiveNorm_NVar > NCOMP_PASSIVE )
         Aux_Error( ERROR_INFO, "total number of normalized passive scalars (%d) exceeds expectation (%d) after adding the field \"%s\" !!\n",
                    PassiveNorm_NVar, NCOMP_PASSIVE, InputLabel );

      if ( FieldIdx < NCOMP_FLUID )
      Aux_Error( ERROR_INFO, "Field index to be normalized (%d) < NCOMP_FLUID (%d) !!\n"
                 "        --> One should not add any main field to the normalization list\n",
                 FieldIdx, NCOMP_FLUID );

      PassiveNorm_VarIdx[NormIdx] = FieldIdx - NCOMP_FLUID;
   }


// return field index
   return FieldIdx;

} // FUNCTION : AddField



//-------------------------------------------------------------------------------------------------------
// Function    :  GetFieldIndex
// Description :  Return the index of the target field
//
// Note        :  1. Usage: FieldIdx = GetFieldIndex( FieldLabel );
//                2. Return Idx_Undefined if the target field cannot be found
//
// Parameter   :  InputLabel : Target field label
//                Check      : Whether or not to terminate the program if the target field cannot be found
//                             --> Accepted options: CHECK_ON / CHECK_OFF
//
// Return      :  Sucess: index of the target field
//                Failed: Idx_Undefined
//-------------------------------------------------------------------------------------------------------
FieldIdx_t GetFieldIndex( char *InputLabel, const Check_t Check )
{

   FieldIdx_t Idx_Out = Idx_Undefined;

   for (int v=0; v<NDefinedField; v++)
   {
      if (  strcmp( FieldLabel[v], InputLabel ) == 0  )
      {
         Idx_Out = v;
         break;
      }
   }

   if ( Check == CHECK_ON  &&  Idx_Out == Idx_Undefined )
      Aux_Error( ERROR_INFO, "cannot find the target field \"%s\" !!\n", InputLabel );

   return Idx_Out;

} // FUNCTION : GetFieldIndex



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Field_User_Template
// Description :  Template of adding user-defined fields
//
// Note        :  1. Invoked by Init_Field() using the function pointer "Init_Field_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_Field_User_Template()
{

// example
// Idx_NewField = AddField( "NewFieldLabel", NORMALIZE_YES );

} // FUNCTION : Init_Field_User_Template
