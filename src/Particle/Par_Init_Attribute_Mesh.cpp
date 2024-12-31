#include "GAMER.h"

#ifdef PARTICLE

static long Par_Get_MeshIndex( const char *InputLabel );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Get_MeshIndex
// Description :  Retrieve the bitwise index corresponding to the input field label
//
// Note        :  1. Invoked by Par_Init_Attribute_Mesh()
//                2. Similar to GetFieldIndex(), but supports derived fields
//                3. Return Idx_Undefined if the target field cannot be found
//
// Parameter   :  InputLabel : Target field label
//
// Return      :  Success: bitwise index of the target field
//                Failed : Idx_Undefined
//-------------------------------------------------------------------------------------------------------
long Par_Get_MeshIndex( const char *InputLabel )
{

// predefined active and passive fields
   for (int v=0; v<NCOMP_TOTAL; v++)
      if (  strcmp( FieldLabel[v], InputLabel ) == 0  )
         return BIDX(v);

// derived fields
#  if ( MODEL == HYDRO )
   if (  strcmp( "VelX", InputLabel ) == 0  )   return _VELX;
   if (  strcmp( "VelY", InputLabel ) == 0  )   return _VELY;
   if (  strcmp( "VelZ", InputLabel ) == 0  )   return _VELZ;
   if (  strcmp( "Pres", InputLabel ) == 0  )   return _PRES;
   if (  strcmp( "Temp", InputLabel ) == 0  )   return _TEMP;
   if (  strcmp( "Entr", InputLabel ) == 0  )   return _ENTR;
   if (  strcmp( "Eint", InputLabel ) == 0  )   return _EINT;
#  ifdef MHD
   if (  strcmp( "MagX", InputLabel ) == 0  )   return _MAGX_CC;
   if (  strcmp( "MagY", InputLabel ) == 0  )   return _MAGY_CC;
   if (  strcmp( "MagZ", InputLabel ) == 0  )   return _MAGZ_CC;
   if (  strcmp( "MagE", InputLabel ) == 0  )   return _MAGE_CC;
#  endif
#  ifdef GRAVITY
   if (  strcmp( "Pote", InputLabel ) == 0  )   return _POTE;
#  endif
//###REVISE: support derived fields defined in Flu_DerivedField_BuiltIn() and Flu_DerivedField_User()
#  endif // #if ( MODEL == HYDRO )


   return Idx_Undefined;

} // FUNCTION : Par_Get_MeshIndex



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Attribute_Mesh
// Description :  Parse the fields specified in the file Input__Par_Mesh and initialize
//                associated variables
//
// Note        :  1. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_Attribute_Mesh()
{

   char FileName[] = "Input__Par_Mesh";

// retrieve the number of mesh quantities specified in Input__Par_Mesh
   const int Mesh_NAttr = Aux_CountRow( FileName );

   amr->Par->Mesh_Attr_Num = Mesh_NAttr;

   if ( !Mesh_NAttr )
   {
      Aux_Message( stderr, "WARNING : no fields have been specified in \"%s\" !!\n", FileName );

      return;
   }


// allocate the memory
   amr->Par->Mesh_Attr       = (real_par**) malloc( Mesh_NAttr*sizeof(real_par*) );
   amr->Par->Mesh_Attr_Label = (char    **) malloc( Mesh_NAttr*sizeof(char    *) );
   amr->Par->Mesh_Attr_Idx   = (long    * ) malloc( Mesh_NAttr*sizeof(long     ) );

   for (int v=0; v<Mesh_NAttr; v++)
   {
      amr->Par->Mesh_Attr      [v] = NULL;
      amr->Par->Mesh_Attr_Label[v] = (char*) malloc( MAX_STRING*sizeof(char) );
   }


// load data in Input__Par_Mesh and convert the field name to field index
   char FirstItem[MAX_STRING];
   int  NRow=0, NItem;

   char *Line = new char [MAX_STRING];
   FILE *File = fopen( FileName, "r" );

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      NItem = sscanf( Line, "%s", FirstItem );

//    skip empty lines and lines starting with #
      if ( NItem <= 0  ||  FirstItem[0] == '#' )   continue;

//    initialize the bitwise indices and labels of mesh quantities to be mapped from
      const long FieldIdx = Par_Get_MeshIndex( FirstItem );

      if ( FieldIdx == Idx_Undefined )
         Aux_Error( ERROR_INFO, "unknown input field label (%s) in %s !!\n", FirstItem, FileName );

      amr->Par->Mesh_Attr_Idx[NRow] = FieldIdx;

      sprintf( amr->Par->Mesh_Attr_Label[NRow], "Mesh%s", FirstItem );

      NRow++;
   }


   fclose( File );
   delete [] Line;

} // FUNCTION : Par_Init_Attribute_Mesh



#endif // #ifdef PARTICLE
