#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_MemFree
// Description :  Check the total free memory and terminate the program if it is below a given threshold
//
// Note        :  1. The minimum free memory is specified by the variable "OPT__CK_MEMFREE"
//                2. This check will be performed every "global step"
//                   --> included in the function "Aux_Check"
//                3. The total free memory is estimated as the sum of the free, buffer, and cached memories
//
// Parameter   :  MinMemFree_Total  : Minimum total free memory (in GB)
//                comment           : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_MemFree( const double MinMemFree_Total, const char *comment )
{

   const int  StrSize               = 128;
   const char FileName_Mem[StrSize] = "/proc/meminfo";

   char   Useless[2][StrSize], *line=NULL;
   bool   GetMemTotal=false, GetMemFree=false, GetBuffers=false, GetCached=false;
   char   MemTotal_c[StrSize], MemFree_c[StrSize], Buffers_c[StrSize], Cached_c[StrSize];
   double MemTotal_f, MemFree_f, Buffers_f, Cached_f;
   double MemFree_Total;
   size_t len=0;


// 1. read the memory information
   if ( !Aux_CheckFileExist(FileName_Mem) )
   {
      Aux_Message( stderr, "WARNING : memory information file \"%s\" does not exist (Rank %d) !!\n",
                   FileName_Mem, MPI_Rank );
      return;
   }

   FILE *MemFile = fopen( FileName_Mem, "r" );

   while (  !GetMemTotal  ||  !GetMemFree  ||  !GetBuffers  ||  !GetCached  )
   {
      if ( getline( &line, &len, MemFile ) == -1 )
      {
         Aux_Message( stderr, "WARNING : some memory information is not found at Rank %d ", MPI_Rank );
         Aux_Message( stderr, "(MemTotal: %s, MemFree: %s, Buffers: %s, Cached: %s)\n",
                      GetMemTotal? "OK":"NO", GetMemFree? "OK":"NO", GetBuffers? "OK":"NO", GetCached? "OK":"NO" );
         break;
      }

      if      ( strncmp( line, "MemTotal:", 9 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], MemTotal_c, Useless[1] );
         GetMemTotal = true;
      }

      else if ( strncmp( line, "MemFree:", 8 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], MemFree_c, Useless[1] );
         GetMemFree = true;
      }

      else if ( strncmp( line, "Buffers:", 8 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], Buffers_c, Useless[1] );
         GetBuffers = true;
      }

      else if ( strncmp( line, "Cached:", 7 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], Cached_c, Useless[1] );
         GetCached = true;
      }
   } // while (  !GetMemTotal  ||  !GetMemFree  ||  !GetBuffers  ||  !GetCached  )

   fclose( MemFile );

   if ( line != NULL )  free( line );


// 2. estimate the total free memory
   MemTotal_f = (double)atof( MemTotal_c )*1.0e-6;  // char (kB) --> double (GB)
   MemFree_f  = (double)atof( MemFree_c  )*1.0e-6;
   Buffers_f  = (double)atof( Buffers_c  )*1.0e-6;
   Cached_f   = (double)atof( Cached_c   )*1.0e-6;

   MemFree_Total = MemFree_f + Buffers_f + Cached_f;


// 3. terminate the program if the total free memory is below the threshold
   int Terminate_global, Terminate_local;

   if ( MemFree_Total < MinMemFree_Total )
   {
      Terminate_local = true;

      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n", comment, __FUNCTION__, Time[0], Step );
      Aux_Message( stderr, "   Total free memory (%5.2f GB) is below the threshold (%5.2f GB) at Rank %d !!\n",
                   (double)MemFree_Total, (double)MinMemFree_Total, MPI_Rank );
      Aux_Message( stderr, "   MemTotal %5.2f GB, MemFree %5.2f GB, Buffers %5.2f GB, Cached %5.2f GB\n",
                   (double)MemTotal_f, (double)MemFree_f, (double)Buffers_f, (double)Cached_f );
      Aux_Message( stderr, "\n   The program is going to be terminated automatically ...\n\n" );
   }
   else
      Terminate_local = false;

// the program will be terminated as long as ONE process has detected the low free memory
   MPI_Allreduce( &Terminate_local, &Terminate_global, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD );

   if ( Terminate_global )
   {
//    output data
      Output_DumpData( 2 );

      End_GAMER();
   }

} // FUNCTION : Aux_Check_MemFree


