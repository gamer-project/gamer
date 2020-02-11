#include "GAMER.h"
#include <sys/stat.h>




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CheckFileExist
// Description :  Check whether or not the target file exists
//
// Note        :  Use the "stat" function to query the existence of the target file
//
// Parameter   :  FileName : Name of the target file
//
// Return      :  true/false <-> file exists/not exists
//-------------------------------------------------------------------------------------------------------
bool Aux_CheckFileExist( const char *FileName )
{

   struct stat Buf;
   return ( stat(FileName,&Buf) == 0 );

} // FUNCTION : Aux_CheckFileExist
