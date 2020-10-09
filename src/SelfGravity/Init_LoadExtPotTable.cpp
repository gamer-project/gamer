#include "GAMER.h"
#include "ReadPara.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_LoadExtPotTable
// Description :  Load the external potential table
//
// Note        :  1. Invoked by Init_ExtAccPot()
//                2. Enabled by the runtime option "OPT__EXT_POT=EXT_POT_TABLE"
//                3. Floating-point type is set by FLOAT8 for now
//                   --> EXT_POT_TABLE_FLOAT8 is NOT supported yet
//                4. The loaded table will be sent to GPU by invoking CUAPI_SendExtPotTable2GPU()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_LoadExtPotTable()
{

// only for EXT_POT_TABLE
   if ( OPT__EXT_POT != EXT_POT_TABLE )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ...\n", __FUNCTION__ );


// checks
// parameters
   for (int d=0; d<3; d++) {
      if ( EXT_POT_TABLE_NCELL[d] < 2 )
         Aux_Error( ERROR_INFO, "EXT_POT_TABLE_NCELL[%d] = %d < 2 !!\n", d, EXT_POT_TABLE_NCELL[d] );
   }

   if ( EXT_POT_TABLE_DH <= 0.0 )
      Aux_Error( ERROR_INFO, "EXT_POT_TABLE_DH = %14.7e <= 0.0 !!\n", EXT_POT_TABLE_DH );

#  ifdef FLOAT8
   if ( ! EXT_POT_TABLE_FLOAT8 )
      Aux_Error( ERROR_INFO, "must enable EXT_POT_TABLE_FLOAT8 since FLOAT8 is on !!\n" );
#  else
   if (   EXT_POT_TABLE_FLOAT8 )
      Aux_Error( ERROR_INFO, "must disable EXT_POT_TABLE_FLOAT8 since FLOAT8 is off !!\n" );
#  endif

// table must cover the simulation domain
   for (int d=0; d<3; d++) {
      if ( EXT_POT_TABLE_EDGEL[d] == NoDef_double )
         Aux_Error( ERROR_INFO, "EXT_POT_TABLE_EDGEL[%d] has not been set !!\n", d );

      if ( EXT_POT_TABLE_EDGEL[d] > amr->BoxEdgeL[d] )
         Aux_Error( ERROR_INFO, "incorrect left boundary of the external potential table (table=%13.7e > simulation=%13.7e) !!\n",
                    EXT_POT_TABLE_EDGEL[d], amr->BoxEdgeL[d] );

      if ( EXT_POT_TABLE_EDGEL[d]+EXT_POT_TABLE_NCELL[d]*EXT_POT_TABLE_DH < amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "incorrect right boundary of the external potential table (table=%13.7e < simulation=%13.7e) !!\n",
                    EXT_POT_TABLE_EDGEL[d]+EXT_POT_TABLE_NCELL[d]*EXT_POT_TABLE_DH, amr->BoxEdgeR[d] );
   }

// file existence
   if ( ! Aux_CheckFileExist(EXT_POT_TABLE_NAME) )
      Aux_Error( ERROR_INFO, "external potential table \"%s\" does not exist !!\n", EXT_POT_TABLE_NAME );

// file size
   FILE *FileTemp = fopen( EXT_POT_TABLE_NAME, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long NCell3D    = (long)EXT_POT_TABLE_NCELL[0]*EXT_POT_TABLE_NCELL[1]*EXT_POT_TABLE_NCELL[2];
   const long ExpectSize = NCell3D*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the external potential table <%s> (%ld) != expect (%ld) !!\n",
                 EXT_POT_TABLE_NAME, FileSize, ExpectSize );

   fclose( FileTemp );

// memory allocation
   if ( h_ExtPotTable == NULL )
      Aux_Error( ERROR_INFO, "h_ExtPotTable[] has not been allocated !!\n" );

   MPI_Barrier( MPI_COMM_WORLD );


// load table to CPU
   FILE *File = fopen( EXT_POT_TABLE_NAME, "rb" );
   fread( h_ExtPotTable, sizeof(real), NCell3D, File );
   fclose( File );


// transfer table to GPU
#  ifdef GPU
   CUAPI_SendExtPotTable2GPU( h_ExtPotTable );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_LoadExtPotTable



#endif // #ifdef GRAVITY
