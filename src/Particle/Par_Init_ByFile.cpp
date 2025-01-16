#include "GAMER.h"

#ifdef PARTICLE
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Par_Init_ByFile_Default();

// this function pointer may be overwritten by various test problem initializers
// --> link to Par_Init_ByFile_Default() by default
void (*Par_Init_ByFile_User_Ptr)() = Par_Init_ByFile_Default;




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
//                      [, user-specified attributes (when PAR_NATT_FLT_USER>0 or PAR_NATT_INT_USER>0)]
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
// Return      :  amr->Par->AttributeFlt[], amr->Par->AttributeInt[]
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( amr->Par->ParICFormat != PAR_IC_FORMAT_ID_ATT  &&  amr->Par->ParICFormat != PAR_IC_FORMAT_ATT_ID )
      Aux_Error( ERROR_INFO, "unknown data format in PAR_IC (%d) !!\n", amr->Par->ParICFormat );


   Par_Init_ByFile_User_Ptr();


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCITON : Par_Init_ByFile



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFile_Default
// Description :  Initialize particle attributes from a file
//
// Note        :  1. Refer to the "Note" section in "Particle/Par_Init_ByFile.cpp -> Par_Init_ByFile()"
//
// Parameter   :  None
//
// Return      :  amr->Par->Attribute[]
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFile_Default()
{

   const char FileName[]    = "PAR_IC";
   const long NParAllRank   = amr->Par->NPar_Active_AllRank;
         long NParThisRank  = amr->Par->NPar_AcPlusInac;       // cannot be "const" due to MPI_Allgather()
   const bool SingleParMass = amr->Par->ParICMass >= 0.0;
   const bool SingleParType = amr->Par->ParICType >= 0;

// determine the number of attributes to be loaded
   int NParAttFlt = PAR_NATT_FLT_TOTAL - 1;   // exclude time
   int NParAttInt = PAR_NATT_INT_TOTAL;
#  ifdef STORE_PAR_ACC
   NParAttFlt -= 3;                       // exclude acceleration
#  endif
   if ( SingleParMass )    NParAttFlt --; // exclude mass
   if ( SingleParType )    NParAttInt --; // exclude type


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist for PAR_INIT == PAR_INIT_BY_FILE !!\n", FileName );

// determine the load_data_size for loading PAR_IC
   size_t load_data_size_flt = ( PAR_IC_FLOAT8 ) ? sizeof(double) : sizeof(float);
   size_t load_data_size_int = ( PAR_IC_INT8   ) ? sizeof(long)   : sizeof(int)  ;

   FILE *FileTemp = fopen( FileName, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long ExpectSize = NParAllRank*( long(NParAttFlt)*load_data_size_flt + long(NParAttInt)*load_data_size_int );
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


// load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading data ... " );

// allocate buffer for loading PAR_IC
   char* ParFltData_ThisRank = new char [ NParThisRank*NParAttFlt*load_data_size_flt ];
   char* ParIntData_ThisRank = new char [ NParThisRank*NParAttInt*load_data_size_int ];

// note that fread() may fail for large files if sizeof(size_t) == 4 instead of 8
   FILE *File = fopen( FileName, "rb" );

   if ( amr->Par->ParICFormat == PAR_IC_FORMAT_ID_ATT )
   {
      for (int r=0; r<MPI_Rank; r++)   FileOffset += NPar_EachRank[r]*( long(NParAttFlt)*load_data_size_flt + long(NParAttInt)*load_data_size_int );

      fseek( File, FileOffset, SEEK_SET );

      for (long p=0; p<NParThisRank; p++)
      {
         fread( ParFltData_ThisRank+p*NParAttFlt*load_data_size_flt, load_data_size_flt, long(NParAttFlt), File );
         fread( ParIntData_ThisRank+p*NParAttInt*load_data_size_int, load_data_size_int, long(NParAttInt), File );
      }
   }

   else
   {
      for (int r=0; r<MPI_Rank; r++)   FileOffset += NPar_EachRank[r]*load_data_size_flt;

      for (int v=0; v<NParAttFlt; v++)
      {
         fseek( File, FileOffset, SEEK_SET );
         fread( ParFltData_ThisRank+v*NParThisRank*load_data_size_flt, load_data_size_flt, NParThisRank, File );
         FileOffset += NParAllRank*load_data_size_flt;
      }

      for (int v=0; v<NParAttInt; v++)
      {
         fseek( File, FileOffset, SEEK_SET );
         fread( ParIntData_ThisRank+v*NParThisRank*load_data_size_int, load_data_size_int, NParThisRank, File );
         FileOffset += NParAllRank*load_data_size_int;
      }
   } // if ( amr->Par->ParICFormat == PAR_IC_FORMAT_ID_ATT ) ... else ...

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing data into particle repository ... " );

   real_par *ParFltData1 = new real_par [NParAttFlt];
   long_par *ParIntData1 = new long_par [NParAttInt];

   for (long p=0; p<NParThisRank; p++)
   {
//    collect data for the target particle
//    [id][att]
      if ( amr->Par->ParICFormat == PAR_IC_FORMAT_ID_ATT )
      {
         if ( PAR_IC_FLOAT8 )
            for (int v=0; v<NParAttFlt; v++) ParFltData1[v] = (real_par)( *((double*)(ParFltData_ThisRank+(p*NParAttFlt+v)*load_data_size_flt)) );
         else
            for (int v=0; v<NParAttFlt; v++) ParFltData1[v] = (real_par)( *((float* )(ParFltData_ThisRank+(p*NParAttFlt+v)*load_data_size_flt)) );

         if ( PAR_IC_INT8 )
            for (int v=0; v<NParAttInt; v++) ParIntData1[v] = (long_par)( *((long*  )(ParIntData_ThisRank+(p*NParAttInt+v)*load_data_size_int)) );
         else
            for (int v=0; v<NParAttInt; v++) ParIntData1[v] = (long_par)( *((int*   )(ParIntData_ThisRank+(p*NParAttInt+v)*load_data_size_int)) );
      }

//    [att][id]
      else
      {
         if ( PAR_IC_FLOAT8 )
            for (int v=0; v<NParAttFlt; v++) ParFltData1[v] = (real_par)( *((double*)(ParFltData_ThisRank+(v*NParThisRank+p)*load_data_size_flt)) );
         else
            for (int v=0; v<NParAttFlt; v++) ParFltData1[v] = (real_par)( *((float* )(ParFltData_ThisRank+(v*NParThisRank+p)*load_data_size_flt)) );

         if ( PAR_IC_INT8 )
            for (int v=0; v<NParAttInt; v++) ParIntData1[v] = (long_par)( *((long*  )(ParIntData_ThisRank+(v*NParThisRank+p)*load_data_size_int)) );
         else
            for (int v=0; v<NParAttInt; v++) ParIntData1[v] = (long_par)( *((int*   )(ParIntData_ThisRank+(v*NParThisRank+p)*load_data_size_int)) );
      }

//    assuming that the orders of the particle attributes stored on the disk and in Par->AttributeFlt/Int[] are the same
//    --> no need to skip acceleration and time since they are always put at the end of the attribute list
      for (int v_in=0, v_out=0; v_in<NParAttFlt; v_in++, v_out++)
      {
         if ( SingleParMass  &&  v_out == PAR_MASS )  v_out ++;

         amr->Par->AttributeFlt[v_out][p] = ParFltData1[v_in];
      }

      for (int v_in=0, v_out=0; v_in<NParAttInt; v_in++, v_out++)
      {
         if ( SingleParType  &&  v_out == PAR_TYPE )  v_out ++;

         amr->Par->AttributeInt[v_out][p] = ParIntData1[v_in];
      }

      if ( SingleParMass )    amr->Par->Mass[p] = amr->Par->ParICMass;
      if ( SingleParType )    amr->Par->Type[p] = amr->Par->ParICType;

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   } // for (long p=0; p<NParThisRank; p++)

   delete [] ParFltData_ThisRank;
   delete [] ParIntData_ThisRank;
   delete [] ParFltData1;
   delete [] ParIntData1;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Par_Init_ByFile_Default



#endif // #ifdef PARTICLE
