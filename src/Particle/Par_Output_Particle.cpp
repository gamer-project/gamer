#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Output_Particle
// Description :  Output the particle position and velocity
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Output_Particle( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check
   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// output particles
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "#Time %20.14e   Step %13ld   Active / Inactive Particles %13ld / %13ld\n\n",
            Time[0], Step, amr->Par->NPar_Active, amr->Par->NPar-amr->Par->NPar_Active );
   fprintf( File, "#%9s  %21s  %21s  %21s  %21s  %21s  %21s  %21s  %21s", "ID", "Mass", "X", "Y", "Z", "Vx", "Vy", "Vz", "Time" );
#  ifdef STORE_PAR_ACC
   fprintf( File, "  %21s  %21s  %21s", "AccX", "AccY", "AccZ" );
#  endif
   fprintf( File, "\n" );

   for (long p=0; p<amr->Par->NPar; p++)
   {
//    skip inactive particles
      if ( amr->Par->Mass[p] < 0.0 )   continue;

      fprintf( File, "%10ld", p );

      for (int v=0; v<NPAR_VAR;     v++)     fprintf(  File, "  %21.14e", amr->Par->ParVar [v][p] );
      for (int v=0; v<NPAR_PASSIVE; v++)     fprintf(  File, "  %21.14e", amr->Par->Passive[v][p] );

      fprintf( File, "\n" );
   }

   fclose( File );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Par_Output_Particle



#endif // #ifdef PARTICLE
