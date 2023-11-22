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


// 0. initialize variables
   NDefinedField       = 0;
   PassiveNorm_NVar    = 0;
   PassiveIntFrac_NVar = 0;

   for (int v=0; v<NCOMP_PASSIVE; v++)
   {
      PassiveNorm_VarIdx   [v] = -1;
      PassiveIntFrac_VarIdx[v] = -1;
   }


// 1. add main predefined fields
//    --> must not change the following order of declaration since they must be consistent
//        with the symbolic constants defined in Macro.h (e.g., DENS)
#  if   ( MODEL == HYDRO )
   Idx_Dens    = AddField( "Dens",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_MomX    = AddField( "MomX",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_MomY    = AddField( "MomY",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_MomZ    = AddField( "MomZ",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_Engy    = AddField( "Engy",     NORMALIZE_NO, INTERP_FRAC_NO );

   if ( Idx_Dens != DENS )    Aux_Error( ERROR_INFO, "inconsistent Idx_Dens (%d != %d) !!\n", Idx_Dens, DENS );
   if ( Idx_MomX != MOMX )    Aux_Error( ERROR_INFO, "inconsistent Idx_MomX (%d != %d) !!\n", Idx_MomX, MOMX );
   if ( Idx_MomY != MOMY )    Aux_Error( ERROR_INFO, "inconsistent Idx_MomY (%d != %d) !!\n", Idx_MomY, MOMY );
   if ( Idx_MomZ != MOMZ )    Aux_Error( ERROR_INFO, "inconsistent Idx_MomZ (%d != %d) !!\n", Idx_MomZ, MOMZ );
   if ( Idx_Engy != ENGY )    Aux_Error( ERROR_INFO, "inconsistent Idx_Engy (%d != %d) !!\n", Idx_Engy, ENGY );

#  ifdef MHD
   strcpy( MagLabel[MAGX], "MagX" );
   strcpy( MagLabel[MAGY], "MagY" );
   strcpy( MagLabel[MAGZ], "MagZ" );
#  endif

#  elif ( MODEL == ELBDM )


#  if ( ELBDM_SCHEME == ELBDM_WAVE )
   Idx_Dens    = AddField( "Dens",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_Real    = AddField( "Real",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_Imag    = AddField( "Imag",     NORMALIZE_NO, INTERP_FRAC_NO );

   if ( Idx_Dens != DENS )    Aux_Error( ERROR_INFO, "inconsistent Idx_Dens (%d != %d) !!\n", Idx_Dens, DENS );
   if ( Idx_Real != REAL )    Aux_Error( ERROR_INFO, "inconsistent Idx_Real (%d != %d) !!\n", Idx_Real, REAL );
   if ( Idx_Imag != IMAG )    Aux_Error( ERROR_INFO, "inconsistent Idx_Imag (%d != %d) !!\n", Idx_Imag, IMAG );

#  elif ( ELBDM_SCHEME == ELBDM_HYBRID )
   Idx_Dens    = AddField( "Dens",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_Phas    = AddField( "Phase",    NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_Stub    = AddField( "Stub",     NORMALIZE_NO, INTERP_FRAC_NO );
   Idx_Real    = Idx_Phas;
   Idx_Imag    = Idx_Stub;

   if ( Idx_Dens != DENS )    Aux_Error( ERROR_INFO, "inconsistent Idx_Dens  (%d != %d) !!\n", Idx_Dens, DENS );
   if ( Idx_Real != REAL )    Aux_Error( ERROR_INFO, "inconsistent Idx_Real  (%d != %d) !!\n", Idx_Real, REAL );
   if ( Idx_Imag != IMAG )    Aux_Error( ERROR_INFO, "inconsistent Idx_Imag  (%d != %d) !!\n", Idx_Imag, IMAG );
   if ( Idx_Phas != PHAS )    Aux_Error( ERROR_INFO, "inconsistent Idx_Phase (%d != %d) !!\n", Idx_Phas, PHAS );
   if ( Idx_Stub != STUB )    Aux_Error( ERROR_INFO, "inconsistent Idx_Stub  (%d != %d) !!\n", Idx_Stub, STUB );
#  else
#    error : ERROR : unsupported ELBDM_SCHEME !!
#  endif

#  else
#    error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 2. add other predefined fields
#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   Idx_e       = AddField( "Electron", NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_HI      = AddField( "HI",       NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_HII     = AddField( "HII",      NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_HeI     = AddField( "HeI",      NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_HeII    = AddField( "HeII",     NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_HeIII   = AddField( "HeIII",    NORMALIZE_YES, INTERP_FRAC_YES );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   Idx_HM      = AddField( "HM",       NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_H2I     = AddField( "H2I",      NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_H2II    = AddField( "H2II",     NORMALIZE_YES, INTERP_FRAC_YES );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   Idx_DI      = AddField( "DI",       NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_DII     = AddField( "DII",      NORMALIZE_YES, INTERP_FRAC_YES );
   Idx_HDI     = AddField( "HDI",      NORMALIZE_YES, INTERP_FRAC_YES );
   }

// normalize the metallicity field only when adopting the non-equilibrium chemistry
// --> may need a machanism to allow users to overwrite this default setup
   if ( GRACKLE_METAL )
   Idx_Metal   = AddField( "Metal",    (GRACKLE_PRIMORDIAL==GRACKLE_PRI_CHE_CLOUDY)?NORMALIZE_NO:NORMALIZE_YES,
                                                      INTERP_FRAC_YES );
#  endif // #ifdef SUPPORT_GRACKLE


// 3. add user-defined fields
   if ( Init_Field_User_Ptr != NULL )  Init_Field_User_Ptr();


// 4. must put all built-in scalars at the END of the field list and with the same order as their
//    corresponding symbolic constants (e.g., DUAL/CRAY) defined in Macro.h
//    --> as we still rely on these constants (e.g., DENS, DUAL) in the fluid solvers
#  ifdef COSMIC_RAY
   Idx_CRay = AddField( "CRay", NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Idx_CRay != CRAY )    Aux_Error( ERROR_INFO, "inconsistent Idx_CRay (%d != %d) !!\n", Idx_CRay, CRAY );
#  endif

#  ifdef DUAL_ENERGY
   Idx_Dual = AddField( "Dual", NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Idx_Dual != DUAL )    Aux_Error( ERROR_INFO, "inconsistent Idx_Dual (%d != %d) !!\n", Idx_Dual, DUAL );
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
         Aux_Message( stderr, "WARNING : disabling \"OPT__NORMALIZE_PASSIVE\" since there are no passive scalars\n"
                              "          to be normalized !!\n" );
   }


// 7. turn off OPT__INT_FRAC_PASSIVE_* if there are no passive scalars to be interpolated in fractional form
   if ( OPT__INT_FRAC_PASSIVE_LR  &&  PassiveIntFrac_NVar == 0 )
   {
      OPT__INT_FRAC_PASSIVE_LR = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : disabling \"OPT__INT_FRAC_PASSIVE_*\" since there are no passive scalars\n"
                              "          to be interpolated in fractional form !!\n" );
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
//                   (3) add the new field to PassiveNorm_VarIdx[] if Norm==true
//                       --> The sum of passive scalars in this list will be normalized to the gas density
//                       --> sum(passive_scalar_density) == gas_density
//                   (4) add the new field to PassiveIntFrac_VarIdx[] if IntFrac==true
//                       --> These passive scalars will be converted to fracion form during interpolation
//                2. One must invoke AddField() exactly NCOMP_TOTAL times to set the labels of all fields
//                3. Invoked by Init_Field() and various test problem initializers
//
// Parameter   :  InputLabel : Label (i.e., name) of the new field
//                Norm       : whether or not to normalize the new field
//                IntFrac    : whether or not to convert the new field to fraction form during interpolation
//
// Return      :  (1) FieldLabel[]
//                (2) Index of the newly added field
//                (3) PassiveNorm_NVar & PassiveNorm_VarIdx[]
//                (4) PassiveIntFrac_NVar & PassiveIntFrac_VarIdx[]
//-------------------------------------------------------------------------------------------------------
FieldIdx_t AddField( const char *InputLabel, const NormPassive_t Norm, const IntFracPassive_t IntFrac )
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
#  if ( NCOMP_PASSIVE > 0 )
   if ( Norm )
   {
      const int NormIdx = PassiveNorm_NVar ++;

      if ( PassiveNorm_NVar > NCOMP_PASSIVE )
         Aux_Error( ERROR_INFO, "total number of normalized passive scalars (%d) exceeds expectation (%d) after adding the field \"%s\" !!\n",
                    PassiveNorm_NVar, NCOMP_PASSIVE, InputLabel );

      if ( FieldIdx < NCOMP_FLUID )
         Aux_Error( ERROR_INFO, "Field index to be normalized (%d) < NCOMP_FLUID (%d) !!\n"
                    "        --> Set NORMALIZE_NO when adding a built-in fluid field by AddField()\n",
                    FieldIdx, NCOMP_FLUID );

      PassiveNorm_VarIdx[NormIdx] = FieldIdx - NCOMP_FLUID;
   }
#  endif


// set the fractional-form interpolation list
// --> note that PassiveIntFrac_VarIdx[] starts from 0 instead of NCOMP_FLUID
#  if ( NCOMP_PASSIVE > 0 )
   if ( IntFrac )
   {
      const int IntFracIdx = PassiveIntFrac_NVar ++;

      if ( PassiveIntFrac_NVar > NCOMP_PASSIVE )
         Aux_Error( ERROR_INFO, "total number of fractional passive scalars (%d) exceeds expectation (%d) after adding the field \"%s\" !!\n",
                    PassiveIntFrac_NVar, NCOMP_PASSIVE, InputLabel );

      if ( FieldIdx < NCOMP_FLUID )
         Aux_Error( ERROR_INFO, "Field index for interpolation in fractional form (%d) < NCOMP_FLUID (%d) !!\n"
                    "        --> Set INTERP_FRAC_NO when adding a built-in fluid field by AddField()\n",
                    FieldIdx, NCOMP_FLUID );

      PassiveIntFrac_VarIdx[IntFracIdx] = FieldIdx - NCOMP_FLUID;
   }
#  endif


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
FieldIdx_t GetFieldIndex( const char *InputLabel, const Check_t Check )
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
// Idx_NewField = AddField( "NewFieldLabel", NORMALIZE_YES, INTERP_FRAC_YES );

} // FUNCTION : Init_Field_User_Template
