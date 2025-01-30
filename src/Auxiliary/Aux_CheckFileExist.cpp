#include "GAMER.h"
#include <sys/stat.h>
#include <unistd.h>
#include <grp.h>




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
// Parameter   :  FolderName : Name of the target folder
//
// Return      :  true/false <-> folder exists/not exists
//-------------------------------------------------------------------------------------------------------
bool Aux_CheckFolderExist( const char *FolderName )
{

   struct stat Buf;

   if ( stat(FolderName, &Buf) != 0 )    return false; // not exist
   else if ( !(Buf.st_mode & S_IFDIR) )  return false; // not a directoy

   return true;

} // FUNCTION : Aux_CheckFolderExist



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CheckPermission
// Description :  Check whether or not the target file has requested permissions
//
// Note        :  Use the "stat" function to query the permission of the target file
//
// Parameter   :  FileName : Name of the target file
//             :  perms    : Permissions code of the target file (sum of 4->read, 2->write, 1->execute)
//                           --> Examples: 1. read and write   permissions => 6 (4+2)
//                                         2. read and execute permissions => 5 (4+1)
//
// Return      :  true/false <-> you does/does not have the file permissions
//-------------------------------------------------------------------------------------------------------
bool Aux_CheckPermission( const char *FileName, const int perms )
{

   if ( perms < 0  ||  perms > 7 )   Aux_Error( ERROR_INFO, "Incorrect file permission code %d (0~7) !!\n", perms );

   const uid_t curUserId  = getuid(); // get user  id
   const gid_t curGroupId = getgid(); // get group id

   struct stat Buf;

   if ( stat(FileName, &Buf) != 0 )    Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   int perm_r, perm_w, perm_x;
   if      ( Buf.st_uid == curUserId  ) { perm_r = S_IRUSR; perm_w = S_IWUSR; perm_x = S_IXUSR; } // user
   else if ( Buf.st_gid == curGroupId ) { perm_r = S_IRGRP; perm_w = S_IWGRP; perm_x = S_IXGRP; } // group
   else                                 { perm_r = S_IROTH; perm_w = S_IWOTH; perm_x = S_IXOTH; } // other

   if ( (perms & 4)  &&  !(Buf.st_mode & perm_r) )   return false;
   if ( (perms & 2)  &&  !(Buf.st_mode & perm_w) )   return false;
   if ( (perms & 1)  &&  !(Buf.st_mode & perm_x) )   return false;

   return true;

} // FUNCTION : Aux_CheckPermission
