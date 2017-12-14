#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemoryPool
// Description :  Initialize the memory pool for the patch data
//
// Note        :  1. Load the number of preallocated patches at each level from the table "Input__MemoryPool"
//                   --> Set the numbers at higher levels to zero if they are not specified in the table
//                   --> Currently the table must have one header line
//                2. Preallocate patches with both fluid and pot data allocated
//                3. Controlled by the option "OPT__MEMORY_POOL"
//                   --> Must turn on "OPT__REUSE_MEMORY" as well
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_MemoryPool()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


   if ( ! OPT__REUSE_MEMORY )    Aux_Error( ERROR_INFO, "Please turn on OPT__REUSE_MEMORY for OPT__MEMORY_POOL !!\n" );


   const char FileName[] = "Input__MemoryPool";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   char  *input_line = NULL;
   size_t len = 0;
   int    Trash, n;

   int *NPatchInPool = new int [MAX_LEVEL+1];

   for (int lv=0; lv<=MAX_LEVEL; lv++)    NPatchInPool[lv] = 0;

// skip the header
   getline( &input_line, &len, File );

// begin to read
   for (int lv=0; lv<=MAX_LEVEL; lv++)
   {
      n = getline( &input_line, &len, File );

//    check
      if ( n <= 1 )
      {
         if ( MPI_Rank == 0 )
         {
            Aux_Message( stderr, "WARNING : cannot load the number of patches to be preallocated at lv %d ...\n", lv );
            Aux_Message( stderr, "          --> no patches will be preallocated at levels >= %d\n", lv );
         }

         break;
      }

      sscanf( input_line, "%d%d", &Trash, NPatchInPool+lv );
   }

   fclose( File );

   if ( input_line != NULL )     free( input_line );


// allocate the memory pool
   const bool WithFluData_Yes = true;
   const bool WithPotData_Yes = true;
   const bool ReuseMemory_Yes = true;

   for (int lv=0; lv<=MAX_LEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Preallocating %8d patches at level %2d ... ", NPatchInPool[lv], lv );

//    preallocate patches (with corner arbitrarily set)
      for (int PID=0; PID<NPatchInPool[lv]; PID++)    amr->pnew( lv, 0, 0, 0, NULL_INT, WithFluData_Yes, WithPotData_Yes );

//    deactivate all patches (which does not deallocate memory)
      for (int PID=0; PID<NPatchInPool[lv]; PID++)    amr->pdelete( lv, PID, ReuseMemory_Yes );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }

   delete [] NPatchInPool;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_MemoryPool
