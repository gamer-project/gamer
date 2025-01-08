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



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CheckFolderExist
// Description :  Check whether or not the target folder exists
//
// Note        :  Use the "stat" function to query the existence of the target folder
//
// Parameter   :  FileName : Name of the target folder
//
// Return      :  true/false <-> file exists/not exists
//-------------------------------------------------------------------------------------------------------
bool Aux_CheckFolderExist( const char *FolderName )
{

   struct stat Buf;

   if( stat( FolderName, &Buf ) != 0 )   return false; // not exist
   else if( Buf.st_mode & S_IFDIR )      return true;  // is directory
   else                                  return false; // not a diretory

} // FUNCTION : Aux_CheckFolderExist
