#ifndef __READPARA_H__
#define __READPARA_H__


#include <typeinfo>
#include "Global.h"

bool Aux_CheckFileExist( const char *FileName );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );

#define NPARA_MAX    1000     // maximum number of parameters
#define TYPE_INT        1     // different data types
#define TYPE_LONG       2
#define TYPE_UINT       3
#define TYPE_ULONG      4
#define TYPE_FLOAT      5
#define TYPE_DOUBLE     6
#define TYPE_BOOL       7
#define TYPE_STRING     8




//-------------------------------------------------------------------------------------------------------
// Structure   :  ReadPara_t
// Description :  Data structure ...
//
// Data Member :
//
// Method      :  ReadPara_t : Constructor
//               ~ReadPara_t : Destructor
//                Add        : Add a new parameter
//                Read       : Read parameters from the parameter file
//-------------------------------------------------------------------------------------------------------
struct ReadPara_t
{

// data members
// ===================================================================================
   int    NPara;
   char (*Key)[MAX_STRING];
   void **Ptr;
   int   *Type;


   //===================================================================================
   // Constructor :  ReadPara_t
   // Description :  Constructor of the structure "ReadPara_t"
   //
   // Note        :  Initialize variables and allocate memory
   //===================================================================================
   ReadPara_t()
   {

      NPara = 0;
      Key   = new char  [NPARA_MAX][MAX_STRING];
      Ptr   = new void* [NPARA_MAX];
      Type  = new int   [NPARA_MAX];

   } // METHOD : ReadPara_t



   //===================================================================================
   // Constructor :  ~ReadPara_t
   // Description :  Destructor of the structure "ReadPara_t"
   //
   // Note        :  Deallocate memory
   //===================================================================================
   ~ReadPara_t()
   {

      delete [] Key;
      delete [] Ptr;
      delete [] Type;

   } // METHOD : ~ReadPara_t



   //===================================================================================
   // Constructor :  Add
   // Description :  Add a new parameter to be loaded later
   //
   // Note        :  1. This function stores the name, address, and data type of the new parameter
   //                2. Data type (e.g., integer, float, ...) is determined by the input pointer
   //===================================================================================
   template <typename T>
   void Add( const char NewKey[], T* NewPtr )
   {

      if ( NPara >= NPARA_MAX )  Aux_Error( ERROR_INFO, "exceed the maximum number of parameters (%d) !!\n", NPARA_MAX );


//    parameter name
      strcpy( Key[NPara], NewKey );

//    parameter address
      Ptr[NPara] = NewPtr;

//    parameter data type
      if      ( typeid(T) == typeid(int   ) )   Type[NPara] = TYPE_INT;
      else if ( typeid(T) == typeid(long  ) )   Type[NPara] = TYPE_LONG;
      else if ( typeid(T) == typeid(uint  ) )   Type[NPara] = TYPE_UINT;
      else if ( typeid(T) == typeid(ulong ) )   Type[NPara] = TYPE_ULONG;
      else if ( typeid(T) == typeid(float ) )   Type[NPara] = TYPE_FLOAT;
      else if ( typeid(T) == typeid(double) )   Type[NPara] = TYPE_DOUBLE;
      else if ( typeid(T) == typeid(bool  ) )   Type[NPara] = TYPE_BOOL;
      else if ( typeid(T) == typeid(char  ) )   Type[NPara] = TYPE_STRING;
      else
         Aux_Error( ERROR_INFO, "unsupported data type for \"%s\" (char*, float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n",
                    NewKey );

      NPara ++;

   } // METHOD : Add


   //===================================================================================
   // Constructor :  Read
   // Description :  Read all parameters added by Add()
   //
   // Note        :  1. Format:   KEY   VALUE
   //                2. Use # to comment out lines
   //===================================================================================
   void Read( const char *FileName )
   {

      if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );


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

//       skip lines with incorrect format (e.g., empyt lines)
         if ( NLoad < 2  &&  MPI_Rank == 0 )
         {
            if ( NLoad == 1 )
               Aux_Error( ERROR_INFO, "cannot find the value to assign to the key \"%s\" at line %d !!\n",
                          LoadKey, LineNum );
            continue;
         }

//       skip lines starting with #
         if ( LoadKey[0] == '#' )   continue;

//       locate the target key
         MatchIdx = -1;
         for (int k=0; k<NPara; k++)
         {
            if (  strcmp( Key[k], LoadKey ) == 0  )
            {
               MatchIdx  = k;
               Key[k][0] = '%';  // reset the key to any strange format to detect duplicate parameters
               break;
            }
         }

         if ( MatchIdx >= 0 )
         {
            switch ( Type[MatchIdx] )
            {
               case TYPE_INT    :   *( (int*   ) Ptr[MatchIdx] ) = (int   )atol( LoadValue );   break;
               case TYPE_LONG   :   *( (long*  ) Ptr[MatchIdx] ) = (long  )atol( LoadValue );   break;
               case TYPE_UINT   :   *( (uint*  ) Ptr[MatchIdx] ) = (uint  )atol( LoadValue );   break;
               case TYPE_ULONG  :   *( (ulong* ) Ptr[MatchIdx] ) = (ulong )atol( LoadValue );   break;
               case TYPE_FLOAT  :   *( (float* ) Ptr[MatchIdx] ) = (float )atof( LoadValue );   break;
               case TYPE_DOUBLE :   *( (double*) Ptr[MatchIdx] ) = (double)atof( LoadValue );   break;
               case TYPE_BOOL   :   *( (bool*  ) Ptr[MatchIdx] ) = (bool  )atol( LoadValue );   break;
               case TYPE_STRING :   strcpy( (char*)Ptr[MatchIdx], LoadValue );                  break;

               default: Aux_Error( ERROR_INFO, "unsupported data type (char*, float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n" );
            }
         }

         else
         {
            Aux_Error( ERROR_INFO, "unrecognizable parameter \"%s\" at line %d (either non-existing or duplicate) !!\n",
                       LoadKey, LineNum );
         }
      } // while ( fgets( Line, MAX_STRING, File ) != NULL )


      fclose( File );
      delete [] Line;


//    set default


//    validate the loaded parameters


   } // METHOD : Read

}; // struct ReadPara_t



// remove symbolic constants since they are only used in this structure
#undef NPARA_MAX
#undef TYPE_INT
#undef TYPE_LONG
#undef TYPE_UINT
#undef TYPE_ULONG
#undef TYPE_FLOAT
#undef TYPE_DOUBLE
#undef TYPE_STRING



#endif // #ifndef __READPARA_H__
