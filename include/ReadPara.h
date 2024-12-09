#ifndef __READPARA_H__
#define __READPARA_H__


#include <typeinfo>
#include "Macro.h"
#include "Global.h"

bool Aux_CheckFileExist( const char *FileName );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );

#define NPARA_MAX    1000     // maximum number of parameters
#define COMMENT_SYM   '#'     // comment symbol
#define TYPE_INT        1     // various data types
#define TYPE_LONG       2
#define TYPE_UINT       3
#define TYPE_ULONG      4
#define TYPE_BOOL       5
#define TYPE_FLOAT      6
#define TYPE_DOUBLE     7
#define TYPE_STRING     8

#define GET_VOID( type, var )    ( *( (type*)(var) ) )   // helper function to convert the (void*) pointer


// handy constants for setting the allowed range of the input parameters
const char   Useless_str[1] = "";
const bool   Useless_bool   = false;
const int    NoMax_int      = (int   )+__INT_MAX__;
const long   NoMax_long     = (long  )+__INT_MAX__;
const uint   NoMax_uint     = (uint  )+__INT_MAX__;
const ulong  NoMax_ulong    = (ulong )+__INT_MAX__;
const float  NoMax_float    = (float )+__FLT_MAX__;
const double NoMax_double   = (double)+__DBL_MAX__;
const int    NoMin_int      = (int   )-__INT_MAX__;
const long   NoMin_long     = (long  )-__INT_MAX__;
const uint   NoMin_uint     = (uint  )0;
const ulong  NoMin_ulong    = (ulong )0;
const float  NoMin_float    = (float )-__FLT_MAX__;
const double NoMin_double   = (double)-__DBL_MAX__;
const float  Eps_float      = (float )+__FLT_MIN__;
const double Eps_double     = (double)+__DBL_MIN__;
const float  NoDef_float    = NoMax_float;
const double NoDef_double   = NoMax_double;
const char   NoDef_str[1]   = "";




//-------------------------------------------------------------------------------------------------------
// Structure   :  ReadPara_t
// Description :  Data structure for loading runtime parameters
//
// Data Member :  NPara  : Total number of parameters
//                Key    : Parameter names
//                Ptr    : Pointer to the parameter variables
//                Type   : Parameter data types
//                Loaded : Record whether each parameter has been properly set or not (true --> set already)
//                Def    : Default values
//                Min    : Minimum allowed values
//                Max    : Maximum allowed values
//
// Method      :  ReadPara_t : Constructor
//               ~ReadPara_t : Destructor
//                Add        : Add a new parameter
//                Read       : Read parameters from the given parameter file
//                SetDefault : Set the default values
//                Validate   : Check if all parameters are within their allowed range
//-------------------------------------------------------------------------------------------------------
struct ReadPara_t
{

// data members
// ===================================================================================
   int    NPara;
   char (*Key)[MAX_STRING];
   void **Ptr;
   int   *Type;
   bool  *Loaded;
   void **Def;
   void **Min;
   void **Max;


   //===================================================================================
   // Constructor :  ReadPara_t
   // Description :  Constructor of the structure "ReadPara_t"
   //
   // Note        :  Initialize variables and allocate memory
   //===================================================================================
   ReadPara_t()
   {

      NPara  = 0;
      Key    = new char  [NPARA_MAX][MAX_STRING];
      Ptr    = new void* [NPARA_MAX];
      Type   = new int   [NPARA_MAX];
      Loaded = new bool  [NPARA_MAX];
      Def    = new void* [NPARA_MAX];
      Min    = new void* [NPARA_MAX];
      Max    = new void* [NPARA_MAX];

      for (int t=0; t<NPARA_MAX; t++)  Loaded[t] = false;

   } // METHOD : ReadPara_t



   //===================================================================================
   // Constructor :  ~ReadPara_t
   // Description :  Destructor of the structure "ReadPara_t"
   //
   // Note        :  Deallocate memory
   //===================================================================================
   ~ReadPara_t()
   {

      for (int t=0; t<NPara; t++)
      {
         free( Def[t] );
         free( Min[t] );
         free( Max[t] );
      }

      delete [] Key;
      delete [] Ptr;
      delete [] Type;
      delete [] Loaded;
      delete [] Def;
      delete [] Min;
      delete [] Max;

   } // METHOD : ~ReadPara_t



   //===================================================================================
   // Constructor :  Add
   // Description :  Add a new parameter to be loaded later
   //
   // Note        :  1. This function stores the name, address, and data type of the new parameter
   //                2. Data type (e.g., integer, float, ...) is determined by the input pointer
   //                3. NewPtr, NewDef, NewMin, and NewMax must have the same data type
   //                4. String parameters are handled by a separate overloaded function
   //===================================================================================
   template <typename T>
   void Add( const char NewKey[], T* NewPtr, T NewDef, T NewMin, T NewMax )
   {

      if ( NPara >= NPARA_MAX )  Aux_Error( ERROR_INFO, "exceed the maximum number of parameters (%d) !!\n", NPARA_MAX );


//    parameter name
      strncpy( Key[NPara], NewKey, MAX_STRING );

//    parameter address
      Ptr[NPara] = NewPtr;

//    parameter data type
      if      ( typeid(T) == typeid(int   ) )   Type[NPara] = TYPE_INT;
      else if ( typeid(T) == typeid(long  ) )   Type[NPara] = TYPE_LONG;
      else if ( typeid(T) == typeid(uint  ) )   Type[NPara] = TYPE_UINT;
      else if ( typeid(T) == typeid(ulong ) )   Type[NPara] = TYPE_ULONG;
      else if ( typeid(T) == typeid(bool  ) )   Type[NPara] = TYPE_BOOL;
      else if ( typeid(T) == typeid(float ) )   Type[NPara] = TYPE_FLOAT;
      else if ( typeid(T) == typeid(double) )   Type[NPara] = TYPE_DOUBLE;
      else
         Aux_Error( ERROR_INFO, "unsupported data type for \"%s\" (float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n",
                    NewKey );

//    set the default value and the allowed range of each variable
//    --> use malloc() and free() since they are more appropriate for the void* pointers
      Def[NPara] = (T*)malloc( sizeof(T) );
      Min[NPara] = (T*)malloc( sizeof(T) );
      Max[NPara] = (T*)malloc( sizeof(T) );

      GET_VOID( T, Def[NPara] ) = NewDef;
      GET_VOID( T, Min[NPara] ) = NewMin;
      GET_VOID( T, Max[NPara] ) = NewMax;

      NPara ++;

   } // METHOD : Add



   //===================================================================================
   // Constructor :  Add (string)
   // Description :  Add a new string parameter to be loaded later
   //
   // Note        :  1. Overloaded function for strings
   //===================================================================================
   void Add( const char NewKey[], char* NewPtr, const char* NewDef, const char* NewMin, const char* NewMax )
   {

      if ( NPara >= NPARA_MAX )  Aux_Error( ERROR_INFO, "exceed the maximum number of parameters (%d) !!\n", NPARA_MAX );


//    parameter name
      strncpy( Key[NPara], NewKey, MAX_STRING );

//    parameter address
      Ptr[NPara] = NewPtr;

//    parameter data type
      Type[NPara] = TYPE_STRING;

//    set the default value
//    --> use malloc() and free() since they are more appropriate for the void* pointers
//    --> string parameter doesn't need range
      Def[NPara] = (char*)malloc( MAX_STRING*sizeof(char) );
      Min[NPara] = NULL;
      Max[NPara] = NULL;

      strncpy( (char*)Def[NPara], NewDef, MAX_STRING );

      NPara ++;

   } // METHOD : Add (string)



   //===================================================================================
   // Constructor :  Read
   // Description :  Read all parameters added by Add()
   //
   // Note        :  1. Format:   KEY   VALUE
   //                2. Use # to comment out lines
   //===================================================================================
   void Read( const char *FileName )
   {

      if ( !Aux_CheckFileExist(FileName) )
         Aux_Error( ERROR_INFO, "runtime parameter file \"%s\" does not exist !!\n", FileName );


      char LoadKey[MAX_STRING], LoadValue[MAX_STRING];
      int  MatchIdx, LineNum=0, NLoad;

      char *Line = new char [MAX_STRING];
      FILE *File = fopen( FileName, "r" );


//    loop over all lines in the target file
      while ( fgets( Line, MAX_STRING, File ) != NULL )
      {
         LineNum ++;

//       load the key and value at the target line
         NLoad = sscanf( Line, "%s%s", LoadKey, LoadValue );

//       skip lines with incorrect format (e.g., empty lines)
         if ( NLoad < 2 )
         {
            if ( NLoad == 1  &&  LoadKey[0] != COMMENT_SYM )
               Aux_Error( ERROR_INFO, "cannot find the value to assign to the key \"%s\" at line %d !!\n",
                          LoadKey, LineNum );
            continue;
         }

//       skip lines starting with the comment symbol
         if ( LoadKey[0] == COMMENT_SYM )    continue;

//       locate the target key
         MatchIdx = -1;
         for (int k=0; k<NPara; k++)
         {
            if (  strcmp( Key[k], LoadKey ) == 0  )
            {
               if ( Loaded[k]  &&  MPI_Rank == 0 )
                  Aux_Message( stderr, "WARNING : duplicate      parameter [%-30s] at line %4d !!\n", LoadKey, LineNum );

               MatchIdx  = k;
               Loaded[k] = true;
               break;
            }
         }

         if ( MatchIdx >= 0 )
         {
            switch ( Type[MatchIdx] )
            {
               case TYPE_INT    :   GET_VOID( int,    Ptr[MatchIdx] ) = (int   )atol( LoadValue );    break;
               case TYPE_LONG   :   GET_VOID( long,   Ptr[MatchIdx] ) = (long  )atol( LoadValue );    break;
               case TYPE_UINT   :   GET_VOID( uint,   Ptr[MatchIdx] ) = (uint  )atol( LoadValue );    break;
               case TYPE_ULONG  :   GET_VOID( ulong,  Ptr[MatchIdx] ) = (ulong )atol( LoadValue );    break;
               case TYPE_BOOL   :   GET_VOID( bool,   Ptr[MatchIdx] ) = (bool  )atol( LoadValue );    break;
               case TYPE_FLOAT  :   GET_VOID( float,  Ptr[MatchIdx] ) = (float )atof( LoadValue );    break;
               case TYPE_DOUBLE :   GET_VOID( double, Ptr[MatchIdx] ) = (double)atof( LoadValue );    break;
               case TYPE_STRING :   strncpy( (char*)Ptr[MatchIdx], LoadValue, MAX_STRING );           break;

               default: Aux_Error( ERROR_INFO, "unsupported data type (char*, float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n" );
            }
         }

         else if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : unrecognizable parameter [%-30s] at line %4d !!\n", LoadKey, LineNum );
      } // while ( fgets( Line, MAX_STRING, File ) != NULL )


      fclose( File );
      delete [] Line;


//    set the parameters missing in the runtime parameter file to their default values
      SetDefault();


//    validate the parameters
      Validate();

   } // METHOD : Read


   //===================================================================================
   // Constructor :  SetDefault
   // Description :  Set parameters missing in the runtime parameter file to their default values
   //
   // Note        :  1. Parameters already set by the runtime parameter file will have Loaded[] == true
   //                2. We do NOT set default values for strings. An error will be raised instead.
   //===================================================================================
   void SetDefault()
   {

      for (int t=0; t<NPara; t++)
      {
         if ( !Loaded[t] )
         {
            long   def_int = NULL_INT;
            double def_flt = NULL_REAL;

            switch ( Type[t] )
            {
               case TYPE_INT    :   def_int = GET_VOID( int,    Ptr[t] ) = GET_VOID( int,    Def[t] );   break;
               case TYPE_LONG   :   def_int = GET_VOID( long,   Ptr[t] ) = GET_VOID( long,   Def[t] );   break;
               case TYPE_UINT   :   def_int = GET_VOID( uint,   Ptr[t] ) = GET_VOID( uint,   Def[t] );   break;
               case TYPE_ULONG  :   def_int = GET_VOID( ulong,  Ptr[t] ) = GET_VOID( ulong,  Def[t] );   break;
               case TYPE_BOOL   :   def_int = GET_VOID( bool,   Ptr[t] ) = GET_VOID( bool,   Def[t] );   break;
               case TYPE_FLOAT  :   def_flt = GET_VOID( float,  Ptr[t] ) = GET_VOID( float,  Def[t] );   break;
               case TYPE_DOUBLE :   def_flt = GET_VOID( double, Ptr[t] ) = GET_VOID( double, Def[t] );   break;
               case TYPE_STRING :   strncpy( (char*)Ptr[t], (char*)Def[t], MAX_STRING );                 break;

               default: Aux_Error( ERROR_INFO, "no default value for this data type (char*, float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n" );
            }

            if ( MPI_Rank == 0 )
            {
               if      ( def_int != NULL_INT )
                  Aux_Message( stdout, "NOTE : parameter [%-30s] is set to the default value [%- 21ld]\n",
                               Key[t], def_int );
               else if ( def_flt != NULL_REAL )
                  Aux_Message( stdout, "NOTE : parameter [%-30s] is set to the default value [%- 21.14e]\n",
                               Key[t], def_flt );
               else
                  Aux_Message( stdout, "NOTE : parameter [%-30s] is set to the default value [%- 21s]\n",
                               Key[t], (char*)Ptr[t] );
            }
         }
      } // for (int t=0; t<NPara; t++)

   } // METHOD : SetDefault


   //===================================================================================
   // Constructor :  Validate
   // Description :  Validate if all parameters are within their correct range
   //
   // Note        :  1. Correct range: Min <= Value <= Max
   //                2. No check for bool and string variables
   //===================================================================================
   void Validate()
   {

#     define CHECK_RANGE( type, format )                                                                                   \
      {                                                                                                                    \
         if (  GET_VOID(type,Ptr[t]) < GET_VOID(type,Min[t])  ||  GET_VOID(type,Ptr[t]) > GET_VOID(type,Max[t])  )         \
            Aux_Error( ERROR_INFO, "parameter \"%s\" = " EXPAND_AND_QUOTE(format) " is not within the correct range ["     \
                       EXPAND_AND_QUOTE(format) " ... " EXPAND_AND_QUOTE(format) "] !!\n",                                 \
                       Key[t], GET_VOID(type,Ptr[t]), GET_VOID(type,Min[t]), GET_VOID(type,Max[t]) );                      \
      }


      for (int t=0; t<NPara; t++)
      {
         if ( Type[t] != TYPE_BOOL  &&  Type[t] != TYPE_STRING )
         {
            switch ( Type[t] )
            {
               case TYPE_INT    :   CHECK_RANGE( int,    %d     );   break;
               case TYPE_LONG   :   CHECK_RANGE( long,   %ld    );   break;
               case TYPE_UINT   :   CHECK_RANGE( uint,   %u     );   break;
               case TYPE_ULONG  :   CHECK_RANGE( ulong,  %lu    );   break;
               case TYPE_FLOAT  :   CHECK_RANGE( float,  %14.7e );   break;
               case TYPE_DOUBLE :   CHECK_RANGE( double, %14.7e );   break;

               default: Aux_Error( ERROR_INFO, "no check for this data type (float*, double*, int*, long*, unit*, ulong* only) !!\n" );
            }
         }
      } // for (int t=0; t<NPara; t++)

#     undef CHECK_RANGE

   } // METHOD : Validate


}; // struct ReadPara_t



// remove symbolic constants and macros only used in this structure
#undef NPARA_MAX
#undef COMMENT_SYM
#undef TYPE_INT
#undef TYPE_LONG
#undef TYPE_UINT
#undef TYPE_ULONG
#undef TYPE_FLOAT
#undef TYPE_DOUBLE
#undef TYPE_BOOL
#undef TYPE_STRING
#undef GET_VOID



#endif // #ifndef __READPARA_H__
