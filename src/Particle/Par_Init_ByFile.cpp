#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFile
// Description :  Initialize particle attributes from a file
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box even for the periodic BC
//                3. Particles lying outside the active region will be removed later by Par_Aux_InitCheck()
//                   if non-periodic B.C. is adopted
//                4. Particles loaded here are only temporarily stored in this rank
//                   --> They will be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> So there is no constraint on which particles should be set by this function
//                5. Currently the target file name is fixed to "PAR_IC"
//                6. The data format of the PAR_IC file is controlled by the runtime parameter "PAR_IC_FORMAT"
//                   --> PAR_IC_FORMAT_ATT_ID: [particle attribute][particle id] in a row-major order
//                       PAR_IC_FORMAT_ID_ATT: [particle id][particle attribute] in a row-major order
//                7  Load the following attributes with the specified order:
//
//                      mass, position x/y/z, velocity x/y/z, type,
//                      [, creation time (when enabling STAR_FORMATION)]
//                      [, user-specified attributes (when PAR_NATT_USER>0)]
//
//                   --> The mass of all particles can be set to PAR_IC_MASS instead (by having PAR_IC_MASS>=0.0),
//                       in which case PAR_IC should exclude partice mass
//                   --> The type of all particles can be set to PAR_IC_TYPE instead (by having PAR_IC_TYPE>=0),
//                       in which case PAR_IC should exclude partice type
//                   --> No need to provide particle acceleration and time
//                8. For LOAD_BALANCE, the number of particles in each rank must be set in advance
//                   --> Currently it's set by Init_Parallelization()
//
// Parameter   :  None
//
// Return      :  amr->Par->Attribute[]
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( amr->Par->ParICFormat != PAR_IC_FORMAT_ID_ATT  &&  amr->Par->ParICFormat != PAR_IC_FORMAT_ATT_ID )
      Aux_Error( ERROR_INFO, "unknown data format in PAR_IC (%d) !!\n", amr->Par->ParICFormat );


   const char FileName[]    = "PAR_IC";
   const long NParAllRank   = amr->Par->NPar_Active_AllRank;
         long NParThisRank  = amr->Par->NPar_AcPlusInac;       // cannot be "const" due to MPI_Allgather()
   const bool SingleParMass = amr->Par->ParICMass >= 0.0;
   const bool SingleParType = amr->Par->ParICType >= 0;

// determine the number of attributes to be loaded
   int NParAtt = PAR_NATT_TOTAL - 1;   // exclude time
#  ifdef STORE_PAR_ACC
   NParAtt -= 3;                       // exclude acceleration
#  endif
   if ( SingleParMass )    NParAtt --; // exclude mass
   if ( SingleParType )    NParAtt --; // exclude type

   const int NParAttPerLoad = ( amr->Par->ParICFormat == PAR_IC_FORMAT_ID_ATT ) ? NParAtt : 1;


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist for PAR_INIT == PAR_INIT_BY_FILE !!\n", FileName );

   FILE *FileTemp = fopen( FileName, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long ExpectSize = long(NParAtt)*NParAllRank*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                 FileName, FileSize, ExpectSize );

   fclose( FileTemp );

   MPI_Barrier( MPI_COMM_WORLD );


// set the file offset for this rank
   long NPar_EachRank[MPI_NRank], NPar_Check=0, FileOffset=0;

   MPI_Allgather( &NParThisRank, 1, MPI_LONG, NPar_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

// check if the total number of particles is correct
   for (int r=0; r<MPI_NRank; r++)  NPar_Check += NPar_EachRank[r];
   if ( NPar_Check != NParAllRank )
      Aux_Error( ERROR_INFO, "total number of particles found (%ld) != expect (%ld) !!\n", NPar_Check, NParAllRank );

   for (int r=0; r<MPI_Rank; r++)   FileOffset = FileOffset + long(NParAttPerLoad)*NPar_EachRank[r]*sizeof(real);


// load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading data ... " );

   real *ParData_ThisRank = new real [ NParThisRank*NParAtt ];

// note that fread() may fail for large files if sizeof(size_t) == 4 instead of 8
   FILE *File = fopen( FileName, "rb" );

   for (int v=0; v<NParAtt; v+=NParAttPerLoad)
   {
      fseek( File, FileOffset+v*NParAllRank*sizeof(real), SEEK_SET );
      fread( ParData_ThisRank+v*NParThisRank, sizeof(real), long(NParAttPerLoad)*NParThisRank, File );
   }

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing data into particle repository ... " );

   real *ParData1 = new real [NParAtt];

   for (long p=0; p<NParThisRank; p++)
   {
//    collect data for the target particle
//    [id][att]
      if ( amr->Par->ParICFormat == PAR_IC_FORMAT_ID_ATT )
         memcpy( ParData1, ParData_ThisRank+p*NParAtt, NParAtt*sizeof(real) );

//    [att][id]
      else
      {
         for (int v=0; v<NParAtt; v++)
            ParData1[v] = ParData_ThisRank[ v*NParThisRank + p ];
      }

//    assuming that the orders of the particle attributes stored on the disk and in Par->Attribute[] are the same
//    --> no need to skip acceleration and time since they are always put at the end of the attribute list
      for (int v_in=0, v_out=0; v_in<NParAtt; v_in++, v_out++)
      {
         if ( SingleParMass  &&  v_out == PAR_MASS )  v_out ++;
         if ( SingleParType  &&  v_out == PAR_TYPE )  v_out ++;

         amr->Par->Attribute[v_out][p] = ParData1[v_in];
      }

      if ( SingleParMass )    amr->Par->Mass[p] = amr->Par->ParICMass;
      if ( SingleParType )    amr->Par->Type[p] = amr->Par->ParICType;

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   }

   delete [] ParData_ThisRank;
   delete [] ParData1;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFile



#endif // #ifdef PARTICLE
