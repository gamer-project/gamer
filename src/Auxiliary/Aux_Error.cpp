#include "GAMER.h"
#include <cstdarg>




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Error
// Description :  Output the error messages and force the program to be terminated
//
// Note        :  Use the variable argument lists provided in "cstdarg" 
// 
// Parameter   :  File     : Name of the file where error occurs 
//                Line     : Line number where error occurs
//                Func     : Name of the function where error occurs
//                Format   : Output format
//                ...      : Arguments in vfprintf
//-------------------------------------------------------------------------------------------------------
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... )
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );


// output error messages
   va_list Arg;
   va_start( Arg, Format );

   Aux_Message ( stderr, "********************************************************************************\n" );
   Aux_Message ( stderr, "ERROR : " );
   vfprintf    ( stderr, Format, Arg );
   Aux_Message ( stderr, "        Rank <%d>, file <%s>, line <%d>, function <%s>\n", 
                 MPI_Rank, File, Line, Func );
   Aux_Message ( stderr, "********************************************************************************\n" );

   va_end( Arg );


// terminate the program   
   MPI_Exit();

} // FUNCTION : Aux_Error
