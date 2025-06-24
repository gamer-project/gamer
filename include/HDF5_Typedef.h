#ifndef __HDF5_TYPEDEF_H__
#define __HDF5_TYPEDEF_H__


#include "hdf5.h"
#include "Macro.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif

#ifdef PARTICLE
# ifdef FLOAT8_PAR
#  define H5T_GAMER_REAL_PAR H5T_NATIVE_DOUBLE
# else
#  define H5T_GAMER_REAL_PAR H5T_NATIVE_FLOAT
# endif

# ifdef INT8_PAR
#  define H5T_GAMER_LONG_PAR H5T_NATIVE_LONG
# else
#  define H5T_GAMER_LONG_PAR H5T_NATIVE_INT
# endif
#endif // #ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    :  SyncHDF5File
// Description :  Force NFS to synchronize (dump data to the disk before return)
//
// Note        :  1. Synchronization is accomplished by first opening the target file with the "appending" mode
//                   and then close it immediately
//                2. It's an inline function
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
inline void SyncHDF5File( const char *FileName )
{

   FILE *File      = fopen( FileName, "ab" );
   int CloseStatus = fclose( File );

   if ( CloseStatus != 0 )
      Aux_Message( stderr, "WARNING : failed to close the file \"%s\" for synchronization !!\n", FileName );

} // FUNCTION : SyncHDF5File



#define NFIELD_MAX   1000     // maximum number of parameters
#define TYPE_INT        1     // various data types
#define TYPE_LONG       2
#define TYPE_UINT       3
#define TYPE_ULONG      4
#define TYPE_BOOL       5
#define TYPE_FLOAT      6
#define TYPE_DOUBLE     7
#define TYPE_STRING     8



//-------------------------------------------------------------------------------------------------------
// Structure   :  HDF5_Output_t
// Description :  Data structure for outputting fields as compound dataset in HDF5 snapshot
//
// Note        :  1. Boolean is stored as integer in the HDF5 snapshot
//                2. For the non-array field, `ArrDim` must be `0` and `ArrLen` must be `NULL`
//                3. The string size MUST be `MAX_STRING`
//                4. The structure can NOT be empty when writing the HDF5 snapshot
//
// Data Member :  NField           : Total number of fields
//                TotalSize        : Size of the structure
//                Key              : Parameter names
//                Ptr              : Pointer to the field variables
//                ArrSize          : The total number of elements in the array
//                ArrDim           : The array dimension
//                ArrLen           : The array length in each dimension
//                Type             : Parameter data types
//                TypeSize         : Size of field data type
//                RS_Compare       : Whether or not to do the comparision during restart
//                RS_FatalNonexist : Whether or not the nonexistence of the target field is fatal
//                RS_FatalCompare  : Whether or not the comparison result is fatal
//                H5_TypeID_VarStr : String type in HDF5
//
// Method      :  HDF5_Output_t         : Constructor
//               ~HDF5_Output_t         : Destructor
//                Add                   : Add a new field
//                CheckAndSetType       : Check the input data type and record it
//                GetH5Type             : Get the HDF5 type
//                Write2CompoundDataset : Write the data to HDF5 compund dataset
//                Copy2Dataset          : Copy the data in the structure to the dataset
//                GetCompound           : Get the HDF5 compound dataset of this structure
//                Compare2File          : Compare the data with the given restart file
//                GetIndexByname        : Get the field array index in the structure
//                CheckNewKey           : Varify the new key name
//-------------------------------------------------------------------------------------------------------
struct HDF5_Output_t
{

// data members
// ===================================================================================
   int       NField;
   size_t    TotalSize;
   char    (*Key)[MAX_STRING];
   void    **Ptr;
   size_t   *ArrSize;
   int      *ArrDim;
   hsize_t **ArrLen;
   int      *Type;
   size_t   *TypeSize;
   bool     *RS_Compare;
   bool     *RS_FatalNonexist;
   bool     *RS_FatalCompare;
   hid_t     H5_TypeID_VarStr;

   //===================================================================================
   // Constructor :  HDF5_Output_t
   // Description :  Constructor of the structure "HDF5_Output_t"
   //
   // Note        :  Initialize variables and allocate memory
   //===================================================================================
   HDF5_Output_t()
   {

      NField           = 0;
      TotalSize        = (size_t)0;
      Key              = new char     [NFIELD_MAX][MAX_STRING];
      Ptr              = new void*    [NFIELD_MAX];
      ArrSize          = new size_t   [NFIELD_MAX];
      ArrDim           = new int      [NFIELD_MAX];
      ArrLen           = new hsize_t* [NFIELD_MAX];
      Type             = new int      [NFIELD_MAX];
      TypeSize         = new size_t   [NFIELD_MAX];
      RS_Compare       = new bool     [NFIELD_MAX];
      RS_FatalNonexist = new bool     [NFIELD_MAX];
      RS_FatalCompare  = new bool     [NFIELD_MAX];

      H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
      H5Tset_size( H5_TypeID_VarStr, MAX_STRING ); // H5T_VARIABLE will cause segmentation fault

   } // METHOD : HDF5_Output_t



   //===================================================================================
   // Constructor :  ~HDF5_Output_t
   // Description :  Destructor of the structure "HDF5_Output_t"
   //
   // Note        :  Deallocate memory
   //===================================================================================
   ~HDF5_Output_t()
   {

      delete [] Key;
      for (int i=0; i<NField; i++)   delete [] (char*)Ptr[i];
      delete [] Ptr;
      delete [] ArrSize;
      delete [] ArrDim;
      for (int i=0; i<NField; i++)   delete [] ArrLen[i];
      delete [] ArrLen;
      delete [] Type;
      delete [] TypeSize;
      delete [] RS_Compare;
      delete [] RS_FatalNonexist;
      delete [] RS_FatalCompare;

      H5Tclose( H5_TypeID_VarStr );

   } // METHOD : ~HDF5_Output_t



   //===================================================================================
   // Constructor :  Add (pointer input)
   // Description :  Add a new field in the structure
   //
   // Note        :  1. This function stores the name, address, data type, array size,
   //                   and dimension of the field
   //                2. Data type (e.g., integer, float, ...) is determined by the input pointer
   //                3. The string size MUST be `MAX_STRING`
   //                4. The memory space of the array MUST be continuous
   //                5. For the non-array field, `Dim` must be `0` and `Len` must be `NULL`
   //
   // Parameter   :  NewKey         : Key name in HDF5 snapshot
   //                NewPtr         : Pointer to the field
   //                Dim            : The array dimension
   //                Len            : The array length in each dimension
  //                 Compare        : Whether or not to do the comparision during restart
   //                Fatal_Nonexist : Whether or not the nonexistence of the target field is fatal
   //                                 --> true  : terminate the program     if the target key cannot be found
   //                                     false : display a warning message if the target key cannot be found
   //                Fatal_Compr    : Whether or not the comparison result is fatal
   //                                 --> true  : terminate the program     if the values are different
   //                                     false : display a warning message if the values are different
   //===================================================================================
   template <typename T>
   void Add( const char NewKey[], T* NewPtr, const int Dim, int *Len, const bool Compare,
             const bool Fatal_Nonexist, const bool Fatal_Compr )
   {

//    check
      if ( NField >= NFIELD_MAX )  Aux_Error( ERROR_INFO, "exceed the maximum number of fields (%d) !!\n", NFIELD_MAX );

      if ( Dim >= __INT_MAX__  ||  Dim < 0 )  Aux_Error( ERROR_INFO, "invalid array dimension (%d) of \"%s\"!!\n", Dim, NewKey );

      if ( Len != NULL  &&  Dim == 0 )  Aux_Error( ERROR_INFO, "key \"%s\" Len must be NULL for single field !!\n", NewKey );
      if ( Len == NULL  &&  Dim != 0 )  Aux_Error( ERROR_INFO, "key \"%s\" Len is NULL for non-single field !!\n", NewKey );

      for (int d=0; d<Dim; d++)
         if ( Len[d] >= __INT_MAX__  ||  Len[d] <= 0 )  Aux_Error( ERROR_INFO, "invalid array length (%d) of \"%s\"!!\n", Len[d], NewKey );

      CheckNewKey( NewKey );

      if (  !CheckAndSetType( typeid(T), Type[NField], TypeSize[NField] )  )
         Aux_Error( ERROR_INFO, "unsupported data type %s for \"%s\"(float*, double*, int*, long*, uint*, ulong*, bool*, and char* only) !!\n",
                    typeid(T).name(), NewKey );

//    field name
      strncpy( Key[NField], NewKey, MAX_STRING );

      RS_Compare      [NField] = Compare;
      RS_FatalNonexist[NField] = Fatal_Nonexist;
      RS_FatalCompare [NField] = Fatal_Compr;
      ArrSize         [NField] = (size_t)1;
      ArrDim          [NField] = Dim;

      if ( Dim == 0 )   ArrLen[NField] = NULL; // must be NULL to avoid segmentation fault when freeing
      else              ArrLen[NField] = new hsize_t [ Dim ];

      for (int d=0; d<Dim; d++)
      {
         ArrLen [NField][d] = (hsize_t)Len[d];
         ArrSize[NField]   *= (size_t )Len[d];
      }

//    copy field from the pointer
      Ptr[NField] = new char [ (int)ArrSize[NField] * (int)TypeSize[NField] ];
      memcpy( Ptr[NField], NewPtr, ArrSize[NField] * TypeSize[NField] );

      TotalSize += ArrSize[NField] * TypeSize[NField];
      NField++;

   } // METHOD : Add



   //===================================================================================
   // Constructor :  Add (value input)
   // Description :  Add a new field in the structure
   //
   // Note        :  1. This function stores the name, value, data type, array size,
   //                   and dimension of the field
   //                2. Data type (e.g., integer, float, ...) is determined by the input pointer
   //                3. `char` data type is not supported
   //                4. `Dim` MUST be `0` and `Len` MUST be `NULL`
   //
   // Parameter   :  NewKey         : Key name in HDF5 snapshot
   //                KeyVal         : Value of the field
   //                Dim            : The array dimension
   //                Len            : The array length in each dimension
   //                Compare        : Whether or not to do the comparision during restart
   //                Fatal_Nonexist : Whether or not the nonexistence of the target field is fatal
   //                                 --> true  : terminate the program     if the target key cannot be found
   //                                     false : display a warning message if the target key cannot be found
   //                Fatal_Compr    : Whether or not the comparison result is fatal
   //                                 --> true  : terminate the program     if the values are different
   //                                     false : display a warning message if the values are different
   //===================================================================================
   template <typename T>
   void Add( const char NewKey[], const T KeyVal, const int Dim, int *Len, const bool Compare,
             const bool Fatal_Nonexist, const bool Fatal_Compr )
   {

//    check
      if ( NField >= NFIELD_MAX )  Aux_Error( ERROR_INFO, "exceed the maximum number of fields (%d) !!\n", NFIELD_MAX );

      if ( Dim != 0  ||  Len != NULL )  Aux_Error( ERROR_INFO, "Len must be NULL and Dim must be zero for the single constant field \"%s\"!!\n", NewKey );

      CheckNewKey( NewKey );

      if (  !CheckAndSetType( typeid(T), Type[NField], TypeSize[NField] )  )
         Aux_Error( ERROR_INFO, "unsupported data type %s for \"%s\"(float, double, int, long, uint, ulong, bool, and char only) !!\n",
                    typeid(T).name(), NewKey );

      if ( typeid(T) == typeid(char) )
         Aux_Error( ERROR_INFO, "passing a constant char value is not supported !!\n" );

//    field name
      strncpy( Key[NField], NewKey, MAX_STRING );

      RS_Compare      [NField] = Compare;
      RS_FatalNonexist[NField] = Fatal_Nonexist;
      RS_FatalCompare [NField] = Fatal_Compr;
      ArrSize         [NField] = (size_t)1;
      ArrDim          [NField] = Dim;
      ArrLen          [NField] = NULL;

//    copy field value
      Ptr[NField] = new T [ (int)ArrSize[NField] ];
      memcpy( Ptr[NField], &KeyVal, ArrSize[NField] * TypeSize[NField] );

      TotalSize += ArrSize[NField] * TypeSize[NField];
      NField++;

   } // METHOD : Add



   //===================================================================================
   // Constructor :  CheckAndSetType
   // Description :  Check if the input data type is valid and set the corresponding type fields
   //
   // Note        :  1. Support int, long, uint, ulong, bool, float, double, and string data types
   //                2. Call-by-reference
   //                3. Boolean is stored as integer, so the type size is `sizeof(int)`
   //                4. String size is always assigned as `MAX_STRING`
   //
   // Parameter   :  type_info : `typeid()` of the input data type
   //                type      : Data member `Type`     of the field
   //                type_size : Data member `TypeSize` of the field
   //
   // Return      :  "true " if data type is     supported
   //                "false" if data type is NOT supported
   //===================================================================================
   bool CheckAndSetType( const std::type_info &type_info, int &type, size_t &type_size )
   {

      if      ( type_info == typeid(int   ) ) { type = TYPE_INT;    type_size = sizeof(int   );     }
      else if ( type_info == typeid(long  ) ) { type = TYPE_LONG;   type_size = sizeof(long  );     }
      else if ( type_info == typeid(uint  ) ) { type = TYPE_UINT;   type_size = sizeof(uint  );     }
      else if ( type_info == typeid(ulong ) ) { type = TYPE_ULONG;  type_size = sizeof(ulong );     }
      else if ( type_info == typeid(bool  ) ) { type = TYPE_BOOL;   type_size = sizeof(int   );     }
      else if ( type_info == typeid(float ) ) { type = TYPE_FLOAT;  type_size = sizeof(float );     }
      else if ( type_info == typeid(double) ) { type = TYPE_DOUBLE; type_size = sizeof(double);     }
      else if ( type_info == typeid(char  ) ) { type = TYPE_STRING; type_size = (size_t)MAX_STRING; }
      else   return false;

      return true;

   } // METHOD : CheckAndSetType



   //===================================================================================
   // Constructor :  GetH5Type
   // Description :  Get the HDF5 type ID
   //
   // Note        :  1. Call-by-reference
   //                2. Boolean type is stored as integer type
   //
   // Parameter   :  type    : The type macro of the field
   //                H5_Type : HDF5 type ID for storing the field
   //===================================================================================
   void GetH5Type( const int type, hid_t &H5_Type )
   {

      switch ( type )
      {
         case TYPE_INT:    H5_Type = H5T_NATIVE_INT;    break;
         case TYPE_LONG:   H5_Type = H5T_NATIVE_LONG;   break;
         case TYPE_UINT:   H5_Type = H5T_NATIVE_UINT;   break;
         case TYPE_ULONG:  H5_Type = H5T_NATIVE_ULONG;  break;
         case TYPE_BOOL:   H5_Type = H5T_NATIVE_INT;    break;
         case TYPE_FLOAT:  H5_Type = H5T_NATIVE_FLOAT;  break;
         case TYPE_DOUBLE: H5_Type = H5T_NATIVE_DOUBLE; break;
         case TYPE_STRING: H5_Type = H5_TypeID_VarStr;  break;
         default: Aux_Error( ERROR_INFO, "unrecognized type: %d !!\n", type ); break;
      } // switch ( type )

   } // METHOD : GetH5Type



   //===================================================================================
   // Constructor :  Write2CompoundDataset
   // Description :  Write all fields in HDF5_Output to an HDF5 compound dataset
   //
   // Note        :  1. Should be called only on master rank
   //
   // Parameter   :  H5_SetID    : HDF5 dataset ID
   //                H5_TypeID   : HDF5 compound type ID associated with HDF5_Output
   //
   // Return      :  H5_Status_write : Status of the write operation
   //===================================================================================
   herr_t Write2CompoundDataset( const hid_t H5_SetID, const hid_t H5_TypeID )
   {

      if ( TotalSize == 0 )   Aux_Error( ERROR_INFO, "TotalSize can not be zero !!\n" );
      if ( MPI_Rank != 0 )   Aux_Message( stderr, "WARNING : H5Dwrite might not work correctly on multiple ranks !!\n");

      herr_t H5_Status;
      char *dataset = new char [ TotalSize ];

      Copy2Dataset( dataset );

      H5_Status = H5Dwrite( H5_SetID, H5_TypeID, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset );

      delete [] dataset;

      return H5_Status;

   } // METHOD : Write2CompoundDataset



   //===================================================================================
   // Constructor :  Copy2Dataset
   // Description :  Copy all the data in the structure to the dataset pointer
   //
   // Note        :  1. The size of the dataset MUST be consistent with `TotalSize`
   //
   // Parameter   :  dataset : the dataset pointer to be stored
   //
   // Return      :  dataset
   //===================================================================================
   void Copy2Dataset( char *dataset )
   {

      size_t offset = (size_t)0;

      for (int i=0; i<NField; i++)
      {
         if ( Type[i] == TYPE_BOOL ) // bool is stored as int
         {
            int *temp = new int [ (int)ArrSize[i] ];

            for (int j=0; j<(int)ArrSize[i]; j++)
               temp[j] = ( *((bool *)(Ptr[i])+j) ) ? 1 : 0;

            memcpy( dataset+offset, temp, ArrSize[i] * TypeSize[i] );

            delete [] temp;
         }
         else
            memcpy( dataset+offset, Ptr[i], ArrSize[i] * TypeSize[i] );

         offset += ArrSize[i] * TypeSize[i];
      } // for (int i=0; i<NField; i++)

   } // METHOD : Copy2Dataset



   //===================================================================================
   // Constructor :  GetCompound
   // Description :  Create the compound dataset for the structure
   //
   // Note        :  1. The returned H5_TypeID must be closed manually
   //                2. Call-by-reference
   //
   // Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
   //===================================================================================
   void GetCompound( hid_t &H5_TypeID )
   {

      if ( TotalSize == (size_t)0 )   Aux_Error( ERROR_INFO, "HDF5_Output_t structure must not be empty !!\n" );

      hid_t  H5_Type, H5_ArrType;

      H5_TypeID = H5Tcreate( H5T_COMPOUND, TotalSize );

      size_t offset = (size_t)0;

      for (int i=0; i<NField; i++)
      {
         GetH5Type( Type[i], H5_Type );

         if ( ArrDim[i] != 0 )
         {
            H5_ArrType = H5Tarray_create( H5_Type, ArrDim[i], ArrLen[i] );

            H5Tinsert( H5_TypeID, Key[i], offset, H5_ArrType );

            H5Tclose( H5_ArrType );
         }
         else
            H5Tinsert( H5_TypeID, Key[i], offset, H5_Type );

         offset += ArrSize[i] * TypeSize[i];
      } // for (int i=0; i<NField; i++)

   } // METHOD : GetCompound



   //===================================================================================
   // Constructor :  Compare2File
   // Description :  Compare the field value in the structure with the field value in
   //                the HDF5 snapshot
   //
   // Parameter   :  KeyName          : Key name of the target field
   //                H5_SetID_Target  : HDF5 dataset  ID of the target compound variable
   //                H5_TypeID_Target : HDF5 datatype ID of the target compound variable
   //                Fatal_Nonexist   : Whether or not the nonexistence of the target field is fatal
   //                                   --> true  : terminate the program     if the target key cannot be found
   //                                       false : display a warning message if the target key cannot be found
   //                Fatal_Compr      : Whether or not the comparison result is fatal
   //                                   --> true  : terminate the program     if the values are different
   //                                       false : display a warning message if the values are different
   //
   // Return      : Success/fail <-> 0/<0
   //===================================================================================
   herr_t Compare2File( const char *KeyName, const hid_t H5_SetID_Target, const hid_t H5_TypeID_Target,
                        const bool Fatal_Nonexist, const bool Fatal_Compr )
   {

//    get field index in structure
      int FieldIdx = GetFieldIdxByKey( KeyName );
      if ( FieldIdx < 0 )
      {
         if ( Fatal_Nonexist )
            Aux_Error( ERROR_INFO, "key \"%s\" does not exist in HDF5_Output_t structure !!\n", KeyName );
         else if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : key \"%s\" does not exist in the HDF5_Output_t structure !!\n", KeyName );

         return -4;
      }

      bool   Check_Pass = true;
      int    H5_FieldIdx;
      size_t H5_FieldSize, MinSize;
      hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
      hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
      herr_t H5_Status;
      char  *FieldPtr;

//    get field index in file
      H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, KeyName );

      if ( H5_FieldIdx < 0 )
      {
         if ( Fatal_Nonexist )
            Aux_Error( ERROR_INFO, "target key \"%s\" does not exist in the restart file !!\n", KeyName );

         else if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : target key \"%s\" does not exist in the restart file !!\n", KeyName );

         return -1;
      } // if ( H5_FieldIdx < 0 )

//    load
      H5_TypeID_Field = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
      H5_FieldSize    = H5Tget_size( H5_TypeID_Field );
      FieldPtr        = new char [ H5_FieldSize ]; // assuming the data type is the same
      MinSize         = MIN( H5_FieldSize, ArrSize[FieldIdx] * TypeSize[FieldIdx] );

      H5_TypeID_Load  = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
      H5_Status       = H5Tinsert( H5_TypeID_Load, KeyName, 0, H5_TypeID_Field );

      H5_Status       = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)FieldPtr );
      if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", KeyName );

      H5_Status       = H5Tclose( H5_TypeID_Field );
      H5_Status       = H5Tclose( H5_TypeID_Load  );

//    check size
      if ( H5_FieldSize != ArrSize[FieldIdx] * TypeSize[FieldIdx]  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" comparison might be incomplete due to the size mismatch (RS: %zu Bytes != RT: %zu Bytes)!!\n",
                      KeyName, H5_FieldSize, ArrSize[FieldIdx] * TypeSize[FieldIdx] );

//    comparison
      char KeyNameWithIdx[MAX_STRING];
      for (int offset=0; offset<(int)ArrSize[FieldIdx]; offset++)
      {
         if ( offset*TypeSize[FieldIdx] >= MinSize )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "offset reaches the minimum size (%zu Bytes), "
                                      "not comparing the rest of the field \"%s\" !!\n",
                          MinSize, KeyName );
               return -5;
            }
            else
            {
               if ( MPI_Rank == 0 )
               Aux_Message( stdout, "WARNING : offset reaches the minimum size (%zu Bytes), "
                                    "not comparing the rest of the field \"%s\" !!\n",
                            MinSize, KeyName );
               break;
            }
         } // if ( offset*TypeSize[FieldIdx] >= MinSize )

         sprintf( KeyNameWithIdx, "%s", KeyName );
         int skip = (int)ArrSize[FieldIdx];
         for (int d=0; d<ArrDim[FieldIdx]; d++)
         {
            skip /= ArrLen[FieldIdx][d];
            int idx = (offset / skip) % ArrLen[FieldIdx][d];
            snprintf( KeyNameWithIdx+strlen(KeyNameWithIdx), sizeof(KeyNameWithIdx)-strlen(KeyNameWithIdx), "[%d]", idx );
         }

#        define COMPARE_NON_STRING_TYPE( field_name, format, fatal_compare, RS_Val, RT_Val )          \
         {                                                                                            \
            if ( RS_Val == RT_Val )   continue;                                                       \
                                                                                                      \
            if ( fatal_compare )                                                                      \
            {                                                                                         \
               Aux_Error( ERROR_INFO, "\"%30s\" : RESTART file (" EXPAND_AND_QUOTE(format) ")"        \
                                      " != runtime (" EXPAND_AND_QUOTE(format) ") !!\n",              \
                          field_name, RS_Val, RT_Val );                                               \
               return -2;                                                                             \
            }                                                                                         \
            else                                                                                      \
            {                                                                                         \
               Check_Pass = false;                                                                    \
               if ( MPI_Rank == 0 )                                                                   \
               Aux_Message( stderr, "WARNING : \"%30s\" : RESTART file (" EXPAND_AND_QUOTE(format) ")"\
                                    " != runtime (" EXPAND_AND_QUOTE(format) ") !!\n",                \
                            field_name, RS_Val, RT_Val );                                             \
            }                                                                                         \
         }

         switch ( Type[FieldIdx] )
         {
            case TYPE_INT:    COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_INT,  Fatal_Compr,       *((int*   )FieldPtr+offset), *((int*   )Ptr[FieldIdx]+offset) ); break;
            case TYPE_LONG:   COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_LONG, Fatal_Compr,       *((long*  )FieldPtr+offset), *((long*  )Ptr[FieldIdx]+offset) ); break;
            case TYPE_UINT:   COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_INT,  Fatal_Compr,       *((uint*  )FieldPtr+offset), *((uint*  )Ptr[FieldIdx]+offset) ); break;
            case TYPE_ULONG:  COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_LONG, Fatal_Compr,       *((ulong* )FieldPtr+offset), *((ulong* )Ptr[FieldIdx]+offset) ); break;
            case TYPE_BOOL:   COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_BOOL, Fatal_Compr, (bool)*((int*   )FieldPtr+offset), *((bool*  )Ptr[FieldIdx]+offset) ); break;
            case TYPE_FLOAT:  COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_REAL, Fatal_Compr,       *((float* )FieldPtr+offset), *((float* )Ptr[FieldIdx]+offset) ); break;
            case TYPE_DOUBLE: COMPARE_NON_STRING_TYPE( KeyNameWithIdx, FORMAT_REAL, Fatal_Compr,       *((double*)FieldPtr+offset), *((double*)Ptr[FieldIdx]+offset) ); break;
            case TYPE_STRING:
               if (  strcmp( (char*)FieldPtr+offset*MAX_STRING, (char*)Ptr[FieldIdx]+offset*MAX_STRING ) == 0 )   continue;
               if ( Fatal_Compr )
               {
                  Aux_Error( ERROR_INFO, "\"%30s\" : RESTART file (%s) != runtime (%s) !!\n",
                             KeyNameWithIdx, (char*)FieldPtr+offset*MAX_STRING, (char*)Ptr[FieldIdx]+offset*MAX_STRING );
                  return -2;
               }
               else
               {
                  Check_Pass = false;
                  if ( MPI_Rank == 0 )
                  Aux_Message( stderr, "WARNING : \"%30s\" : RESTART file (%s) != runtime (%s) !!\n",
                               KeyNameWithIdx, (char*)FieldPtr+offset*MAX_STRING, (char*)Ptr[FieldIdx]+offset*MAX_STRING );
               }
               break;
            default:
               Aux_Error( ERROR_INFO, "Unrecognized type: %d !!\n", Type[FieldIdx] ); break;
         } // switch ( Type[FieldIdx] )

#        undef COMPARE_NON_STRING_TYPE
      } // for (int offset=0; offset<(int)ArrSize[FieldIdx]; offset++)

      delete [] FieldPtr;

      return (Check_Pass) ? 0 : -3;

   } // METHOD : Compare2File



   //===================================================================================
   // Constructor :  Compare2File_All
   // Description :  Compare all field values in the structure with the HDF5 snapshot
   //
   // Parameter   :  H5_SetID_Target  : HDF5 dataset  ID of the target compound variable
   //                H5_TypeID_Target : HDF5 datatype ID of the target compound variable
   //
   // Return      : Success/fail <-> 0/<0
   //===================================================================================
   void Compare2File_All( const hid_t H5_SetID_Target, const hid_t H5_TypeID_Target )
   {

      for (int i=0; i<NField; i++)
      {
         if ( !RS_Compare[i] )   continue;

         Compare2File( Key[i], H5_SetID_Target, H5_TypeID_Target, RS_FatalNonexist[i],
                       RS_FatalCompare[i] );
      } // for (int i=0; i<NField; i++)

   } // METHOD : Compare2File_All



   //===================================================================================
   // Constructor :  GetFieldIdxByKey
   // Description :  Get the index of the key in the data member arrays
   //
   // Note        :  1. Assuming there is no dupilcated key name
   //
   // Parameter   :  KeyName : Key name
   //
   // Return      :  idx : The index of the key in the data member arrays
   //===================================================================================
   int GetFieldIdxByKey( const char *KeyName )
   {

      int idx = -1;

      for (int i=0; i<NField; i++)
      {
         if (  strcmp( Key[i], KeyName ) != 0  )   continue;

         idx = i;
         break;
      } // for (int i=0; i<NField; i++)

      return idx;

   } // METHOD : GetFieldIdxByKey



   //===================================================================================
   // Constructor :  CheckNewKey
   // Description :  Check the new input key
   //
   // Note        :  1. Invoked by the method Add()
   //                2. Check if the key is empty
   //                3. Check if the new key exists in the current structure
   //
   // Parameter   :  NewKey : New key to be added
   //===================================================================================
   void CheckNewKey( const char *NewKey )
   {

      if ( NewKey == NULL )   Aux_Error( ERROR_INFO, "the key string pointer is NULL !!\n" );
      if ( NewKey[0] == '\0' )   Aux_Error( ERROR_INFO, "the key of the field can not be empty !!\n" );

      for (int i=0; i<NField; i++)
      {
         if (  strcmp( Key[i], NewKey ) != 0  )   continue;

         Aux_Error( ERROR_INFO, "Key \"%s\" exists in HDF5_Output_t structure. Please change the key name !!\n", NewKey );
      } // for (int i=0; i<NField; i++)

   } // METHOD : CheckNewKey

}; // struct HDF5_Output_t



#undef NFIELD_MAX
#undef TYPE_INT
#undef TYPE_LONG
#undef TYPE_UINT
#undef TYPE_ULONG
#undef TYPE_FLOAT
#undef TYPE_DOUBLE
#undef TYPE_BOOL
#undef TYPE_STRING



#endif // #ifndef __HDF5_TYPEDEF_H__
