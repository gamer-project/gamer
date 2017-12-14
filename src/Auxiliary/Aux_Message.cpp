#include "GAMER.h"
#include <cstdarg>




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Message
// Description :  Output the message and flush it 
//
// Note        :  Use the variable argument lists provided in "cstdarg" 
// 
// Parameter   :  Type     : stdout/stderr/file stream 
//                Format   : Output format
//                ...      : Arguments in vfprintf
//                           --> It is equivalent to call "fprintf( Type, Format, ... );   fflush( Type );"
//-------------------------------------------------------------------------------------------------------
void Aux_Message( FILE *Type, const char *Format, ... )
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );

   va_list Arg;
   va_start( Arg, Format );

   vfprintf( Type, Format, Arg );
   fflush( Type ); fflush( Type ); fflush( Type );

   va_end( Arg );

} // FUNCTION : Aux_Message

