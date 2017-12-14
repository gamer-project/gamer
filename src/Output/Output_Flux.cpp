#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_Flux
// Description :  Output the flux of a single patch
//
// Parameter   :  lv       : Target refinement level
//                PID      : Target patch ID 
//                Sib      : Target sibling direction of the flux : ( 0,1,2,3,4,5 ) <--> ( -x,+x,-y,+y,-z,+z )
//                comment  : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_Flux( const int lv, const int PID, const int Sib, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )     
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( PID < 0  ||  PID >= MAX_PATCH )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (MAX_PATCH = %d) !!\n", "PID", PID, MAX_PATCH );

   if ( !amr->WithFlux )
      Aux_Message( stderr, "WARNING : invoking %s is useless since no flux is required !!\n", __FUNCTION__ );

   if ( amr->patch[0][lv][PID] == NULL )
   {
      Aux_Message( stderr, "WARNING : level = %d, PID = %d does NOT exist !!\n", lv, PID );
      return;
   }


   patch_t *Relation  = amr->patch[0][lv][PID];

   char FileName[100];
   sprintf( FileName, "Flux_r%d_lv%d_p%d%c%c", MPI_Rank, lv, PID, 45-2*(Sib%2), 120+Sib/2 );
   if ( comment != NULL )       
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = amr->patch[0][lv][PID]->flux[Sib];
   char label[2] = { '?', '?' };

   switch ( Sib )
   {
      case 0: case 1:   label[0] = 'y';  label[1] = 'z';    break;
      case 2: case 3:   label[0] = 'z';  label[1] = 'x';    break;
      case 4: case 5:   label[0] = 'x';  label[1] = 'y';    break;
      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Sib", Sib );
   }


// output header
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Rank  %d    Level  %d    Patch  %d    Local ID  %d    Time  %13.7e    Counter  %ld\n", 
            MPI_Rank, lv, PID, PID%8, Time[lv], AdvanceCounter[lv] );
   fprintf( File, "Father  %d     Son  %d     Corner  (%10d,%10d,%10d)\n\n", Relation->father, Relation->son, 
            Relation->corner[0], Relation->corner[1], Relation->corner[2] );
   fprintf( File, "Flux at %c%c surface\n\n", 45-2*(Sib%2), 120+Sib/2  );

// output flux 
   if ( FluxPtr != NULL )
   {
#     if   ( MODEL == HYDRO )
      fprintf( File, "( %c, %c )%16s%16s%16s%16s%16s\n", label[0], label[1], "Density Flux", "Px Flux", 
                                                         "Py Flux", "Pz Flux", "Energy Flux" ); 
#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      fprintf( File, "( %c, %c )%16s\n",                 label[0], label[1], "Density Flux" );

#     else
#     warning : WARNING : DO YOU WANT TO ADD the FILE HEADER HERE FOR THE NEW MODEL ??
#     endif // MODEL

      for (int m=0; m<PATCH_SIZE; m++)
      for (int n=0; n<PATCH_SIZE; n++)
      {
         fprintf( File, "( %d, %d )", n, m );

         for (int v=0; v<NFLUX_TOTAL; v++)   fprintf( File, "  %14.7e", FluxPtr[v][m][n] );

         fprintf( File, "\n" );
      }
   } // if ( FluxPtr != NULL )

   else
   {
      fprintf( File, "\nNULL\n" );
      Aux_Message( stderr, "WARNING : lv %d, PID %d, Sib %d, the requested flux does NOT exist !!\n",
                   lv, PID, Sib );
   }

   fclose( File );

} // FUNCTION : Output_Flux
