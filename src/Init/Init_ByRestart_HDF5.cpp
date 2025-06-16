#ifdef SUPPORT_HDF5

#include "GAMER.h"
#include "HDF5_Typedef.h"
#include <typeinfo>

void FillIn_HDF5_KeyInfo  ( HDF5_Output_t *HDF5_KeyInfo, int NFieldStored, const bool Load_RS,
                            const int RS_FormatVersion, const bool ReenablePar );
void FillIn_HDF5_Makefile ( HDF5_Output_t *HDF5_Makefile, const int RS_FormatVersion );
void FillIn_HDF5_SymConst ( HDF5_Output_t *HDF5_SymConst, const int RS_FormatVersion );
void FillIn_HDF5_InputPara( HDF5_Output_t *HDF5_InputPara, const int NFieldStored, char FieldLabelOut[][MAX_STRING],
                            const bool Load_RS );

template <typename T>
size_t GetFieldData( const char *FieldName, const hid_t H5_SetID_Target, const hid_t H5_TypeID_Target,
                     T *Ptr, const bool Allocate );
static void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive,
                          const int *SonList, const int (*CrList)[3],
                          const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field,
                          const hid_t *H5_SetID_FCMag, const hid_t *H5_SpaceID_FCMag, const hid_t *H5_MemID_FCMag,
                          const int *NParList, real_par **ParFltBuf, long_par **ParIntBuf, long *NewParList,
                          const hid_t *H5_SetID_ParFltData, const hid_t *H5_SetID_ParIntData,
                          const hid_t H5_SpaceID_ParData, const long *GParID_Offset, const long NParThisRank,
                          const int FormatVersion );



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByRestart_HDF5
// Description :  Reload a previous HDF5 output as the initial condition
//
// Note        :  1. This function will be invoked by "Init_ByRestart" automatically if the restart file
//                   is in the HDF5 format
//                2. Only work for format version >= 2100 (PARTICLE only works for version >= 2200)
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
void Init_ByRestart_HDF5( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( MPI_Rank == 0 )
   {
      if ( !Aux_CheckFileExist(FileName) )
         Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" does not exist !!\n", FileName );

      if ( !H5Fis_hdf5(FileName) )
         Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" is not in the HDF5 format !!\n", FileName );
   }

#  if ( defined PARTICLE  &&  !defined SERIAL  &&  !defined LOAD_BALANCE )
#     error : ERROR : PARTICLE must work with either SERIAL or LOAD_BALANCE !!
#  endif

   MPI_Barrier( MPI_COMM_WORLD );


// 1. load the simulation info
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ...\n" );

   const bool    Fatal             = true;
   const bool NonFatal             = false;
   const bool NotAllocate          = false;
   const bool Load_RS_Yes          = true;

   hid_t  H5_FileID, H5_SetID_KeyInfo, H5_TypeID_KeyInfo, H5_SetID_Cr;
#  ifdef LOAD_BALANCE
   hid_t  H5_SetID_LBIdx;
#  else
   hid_t  H5_SetID_Son;
#  endif
   herr_t H5_Status;
   int    NLvRescale, NPatchAllLv, GID_LvStart[NLEVEL];
   bool   ReenablePar = false;


// 1-1. open the HDF5 file
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 1-2. load the KeyInfo dataset and datatype
   H5_SetID_KeyInfo  = H5Dopen( H5_FileID, "Info/KeyInfo", H5P_DEFAULT );
   if ( H5_SetID_KeyInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/KeyInfo" );

   H5_TypeID_KeyInfo = H5Dget_type( H5_SetID_KeyInfo );
   if ( H5_TypeID_KeyInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the datatype of \"%s\" !!\n", "Info/KeyInfo" );


// storing restart parameters for later usage
   int FormatVersion_RS, Particle_RS, NLevel_RS, DumpID_RS;
   GetFieldData( "FormatVersion", H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &FormatVersion_RS, NotAllocate );
   GetFieldData( "Particle",      H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &Particle_RS,      NotAllocate );
   GetFieldData( "NLevel",        H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &NLevel_RS,        NotAllocate );
   GetFieldData( "DumpID",        H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &DumpID_RS,        NotAllocate );

// format version for HDF5 output
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "      The format version of the HDF5 RESTART file = %ld\n", FormatVersion_RS );

      if ( FormatVersion_RS < 2100 )
         Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2100) !!\n" );

#     ifdef PARTICLE
      if ( FormatVersion_RS < 2200 )
         Aux_Error( ERROR_INFO, "unsupported data format version for PARTICLE (only support version >= 2200) !!\n" );
#     endif

      if ( FormatVersion_RS < 2300 )
         Aux_Message( stderr, "WARNING : loading user-defined fields or particle attributes from version < 2300 will likely fail !!\n" );

#     ifdef MHD
      if ( FormatVersion_RS < 2400 )
         Aux_Error( ERROR_INFO, "unsupported data format version for MHD (only support version >= 2400) !!\n" );
#     endif

#     ifdef SRHD
      if ( FormatVersion_RS < 2473 )
         Aux_Error( ERROR_INFO, "unsupported data format version for SRHD (only support version >= 2473) !!\n" );
#     endif

   }

   MPI_Barrier( MPI_COMM_WORLD );

// support re-enabling PARTICLE from a snapshot without particles, but not vice-versa
#  ifdef PARTICLE
   if ( ! Particle_RS )   ReenablePar = true;
#  else
   if (   Particle_RS )
      Aux_Error( ERROR_INFO, "cannot disable PARTICLE when restarting from a snapshot with particles !!\n" );
#  endif

// runtime NLEVEL must be >= loaded NLEVEL
   if ( NLevel_RS > NLEVEL )
      Aux_Error( ERROR_INFO, "%s : RESTART file (%d) > runtime (%d) !!\n",
                 "NLEVEL", NLevel_RS, NLEVEL );
   else
   {
      NLvRescale = 1 << ( NLEVEL - NLevel_RS );

      if ( MPI_Rank == 0  &&  NLvRescale != 1 )
         Aux_Message( stderr, "WARNING : the NLEVEL rescale factor is set to %d\n", NLvRescale );
   }

   MPI_Barrier( MPI_COMM_WORLD );


// 1-3. compare runtime parameters with the file one-by-one (by all ranks)
   HDF5_Output_t HDF5_KeyInfo_RT;
   FillIn_HDF5_KeyInfo( &HDF5_KeyInfo_RT, NULL_INT, Load_RS_Yes, FormatVersion_RS, ReenablePar );
   HDF5_KeyInfo_RT.Compare2File_All( H5_SetID_KeyInfo, H5_TypeID_KeyInfo );


// 1-4. reset
//      --> set internal parameters
//          1. parameters must be reset
//          2. parameters reset only when OPT__RESTART_RESET is disabled

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   GetFieldData( "UseWaveScheme", H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &amr->use_wave_flag[0], NotAllocate );
#  endif
   GetFieldData( "NPatch",        H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &NPatchTotal[0],        NotAllocate );

// overwrite the current ConRef
   if ( FormatVersion_RS >= 2502 )
   {
      GetFieldData( "ConRef",     H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &ConRef[0],             NotAllocate );
      ConRefInitialized = true;
   }

   if ( ! OPT__RESTART_RESET )
   {
      GetFieldData( "Time",           H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &Time[0],           NotAllocate );
      GetFieldData( "Step",           H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &Step,              NotAllocate );
      GetFieldData( "AdvanceCounter", H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &AdvanceCounter[0], NotAllocate );
      if ( FormatVersion_RS >= 2250 )
      GetFieldData( "dTime_AllLv",    H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &dTime_AllLv[0],    NotAllocate );
#     ifdef GRAVITY
      GetFieldData( "AveDens_Init",   H5_SetID_KeyInfo, H5_TypeID_KeyInfo, &AveDensity_Init,   NotAllocate );
#     endif
   }


// 1-5. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Fclose( H5_FileID );


// 1-6. set parameters in levels that do not exist in the input file
// --> assuming dTime_AllLv[] has been initialized as 0.0 properly
   for (int lv=NLevel_RS; lv<NLEVEL; lv++)
   {
      Time              [lv] = 0.0;
      NPatchTotal       [lv] = 0;
      AdvanceCounter    [lv] = 0;
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      amr->use_wave_flag[lv] = amr->use_wave_flag[ NLevel_RS - 1 ];
#     endif
   }


// 1-7. set SgTime
   for (int lv=0; lv<NLEVEL; lv++)
   {
      amr->FluSgTime[lv][ amr->FluSg[lv] ] = Time[lv];
#     ifdef MHD
      amr->MagSgTime[lv][ amr->MagSg[lv] ] = Time[lv];
#     endif
#     ifdef GRAVITY
      amr->PotSgTime[lv][ amr->PotSg[lv] ] = Time[lv];
#     endif
   }


// 1-8. set the next dump ID
   if ( INIT_DUMPID < 0 )
      DumpID = ( OPT__RESTART_RESET ) ? 0 : DumpID_RS + 1;
   else
      DumpID = INIT_DUMPID;


// 1-10. check all other simulation information (by rank 0 only)
   if ( MPI_Rank == 0 )
   {
      hid_t FID; // file ID
      hid_t SID_Makefile, SID_SymConst, SID_InputPara; // dataset ID
      hid_t TID_Makefile, TID_SymConst, TID_InputPara; // datatype ID
      HDF5_Output_t HDF5_Makefile_RT, HDF5_SymConst_RT, HDF5_InputPara_RT; // RT = RunTime

//    1-10-1. fill in the runtime parameters
      FillIn_HDF5_Makefile ( &HDF5_Makefile_RT, FormatVersion_RS );
      FillIn_HDF5_SymConst ( &HDF5_SymConst_RT, FormatVersion_RS );
      FillIn_HDF5_InputPara( &HDF5_InputPara_RT, NCOMP_TOTAL, FieldLabel, Load_RS_Yes );   // no need to include all output fields here

//    1-10-2. open the HDF5 file
      FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );

      if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

//    1-10-3. load the dataset and datatype
      SID_Makefile  = H5Dopen( FID, "Info/Makefile",  H5P_DEFAULT );
      if ( SID_Makefile  < 0 )
         Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/Makefile"  );
      TID_Makefile  = H5Dget_type( SID_Makefile  );

      SID_SymConst  = H5Dopen( FID, "Info/SymConst",  H5P_DEFAULT );
      if ( SID_SymConst  < 0 )
         Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/SymConst"  );
      TID_SymConst  = H5Dget_type( SID_SymConst  );

      SID_InputPara = H5Dopen( FID, "Info/InputPara", H5P_DEFAULT );
      if ( SID_InputPara < 0 )
         Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/InputPara" );
      TID_InputPara = H5Dget_type( SID_InputPara );

//    1-10-4. compare all target fields one-by-one
      HDF5_Makefile_RT.Compare2File_All ( SID_Makefile , TID_Makefile  );
      HDF5_SymConst_RT.Compare2File_All ( SID_SymConst , TID_SymConst  );
      HDF5_InputPara_RT.Compare2File_All( SID_InputPara, TID_InputPara );

//    1-10-5. close all objects
      H5_Status = H5Tclose( TID_Makefile  );
      H5_Status = H5Tclose( TID_SymConst  );
      H5_Status = H5Tclose( TID_InputPara );
      H5_Status = H5Dclose( SID_Makefile  );
      H5_Status = H5Dclose( SID_SymConst  );
      H5_Status = H5Dclose( SID_InputPara );
      H5_Status = H5Fclose( FID           );
   }


// 1-11. set the GID offset at different levels
   NPatchAllLv = 0;
   for (int lv=0; lv<NLEVEL; lv++)
   {
      NPatchAllLv += NPatchTotal[lv];

      GID_LvStart[lv] = ( lv == 0 ) ? 0 : GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }

   MPI_Barrier( MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ... done\n" );



// 2. load the tree information (load-balance indices, corner, son, ... etc) of all patches (by all ranks)
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

// 2-1. corner
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading corner table ...\n" );

// allocate memory
   int (*CrList_AllLv)[3] = new int [ NPatchAllLv ][3];

// load data
   H5_SetID_Cr = H5Dopen( H5_FileID, "Tree/Corner", H5P_DEFAULT );
   if ( H5_SetID_Cr < 0 )     Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Corner" );
   H5_Status = H5Dread( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CrList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Cr );

// rescale the loaded corner (necessary when NLevel_RS != NLEVEL)
   if ( NLvRescale != 1 )
   for (int GID=0; GID<NPatchAllLv; GID++)
   for (int d=0; d<3; d++)
      CrList_AllLv[GID][d] *= NLvRescale;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading corner table ... done\n" );


// 2-2. LBIdx
#  ifdef LOAD_BALANCE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading load-balance index table ...\n" );

   long *LBIdxList_AllLv;
   long *LBIdxList_EachLv         [NLEVEL];
   int  *LBIdxList_EachLv_IdxTable[NLEVEL];
   int   LoadIdx_Start[NLEVEL], LoadIdx_Stop[NLEVEL];

// 2-2-1. allocate memory
   LBIdxList_AllLv = new long [ NPatchAllLv ];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( lv < NLevel_RS )
      {
         LBIdxList_EachLv         [lv] = LBIdxList_AllLv + GID_LvStart[lv];
         LBIdxList_EachLv_IdxTable[lv] = new int [ NPatchTotal[lv] ];
      }

      else
      {
         LBIdxList_EachLv         [lv] = NULL;
         LBIdxList_EachLv_IdxTable[lv] = NULL;
      }
   }


// 2-2-2. load LBIdx list (sorted by GID)
   H5_SetID_LBIdx = H5Dopen( H5_FileID, "Tree/LBIdx", H5P_DEFAULT );

   if ( H5_SetID_LBIdx < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/LBIdx" );

   H5_Status = H5Dread( H5_SetID_LBIdx, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, LBIdxList_AllLv );
   H5_Status = H5Dclose( H5_SetID_LBIdx );


   for (int lv=0; lv<NLevel_RS; lv++)
   {
//    2-2-3. sort the LBIdx list at each level and set the load-balance cut points
#     if ( LOAD_BALANCE != HILBERT )
      if ( NLvRescale != 1 )
      Aux_Message( stderr, "WARNING : please make sure that the patch LBIdx doesn't change when NLvRescale != 1 !!\n" );
#     endif

//    Mis_Heapsort must be called before LB_SetCutPoint in order to get LBIdxList_EachLv_IdxTable
//    (since LB_SetCutPoint will sort LBIdxList_EachLv as well)
//    --> Actually it's not necessary anymore since we now send LBIdx0_AllRank instead of LBIdxList_EachLv into LB_SetCutPoint()
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_EachLv[lv], LBIdxList_EachLv_IdxTable[lv] );

//    prepare LBIdx and load-balance weighting of each **patch group** for LB_SetCutPoint()
      const bool InputLBIdx0AndLoad_Yes = true;
      long   *LBIdx0_AllRank = NULL;
      double *Load_AllRank   = NULL;

      if ( MPI_Rank == 0 )
      {
         LBIdx0_AllRank = new long   [ NPatchTotal[lv] / 8 ];
         Load_AllRank   = new double [ NPatchTotal[lv] / 8 ];

         for (int t=0; t<NPatchTotal[lv]/8; t++)
         {
            LBIdx0_AllRank[t]  = LBIdxList_EachLv[lv][t*8];
            LBIdx0_AllRank[t] -= LBIdx0_AllRank[t] % 8;     // assuming LBIdxList_EachLv is NOT sorted yet
            Load_AllRank  [t]  = 8.0;                       // assuming all patches have the same weighting == 1.0
         }
      }

//    do NOT consider load-balance weighting of particles since at this point we don't have that information
      const double ParWeight_Zero = 0.0;
      LB_SetCutPoint( lv, NPatchTotal[lv]/8, amr->LB->CutPoint[lv], InputLBIdx0AndLoad_Yes, LBIdx0_AllRank, Load_AllRank,
                      ParWeight_Zero );

//    free memory
      if ( MPI_Rank == 0 )
      {
         delete [] LBIdx0_AllRank;
         delete [] Load_AllRank;
      }


//    2-2-4. get the target LBIdx range of each rank
      LoadIdx_Start[lv] = -1;    // -1 --> indicate that it has not been set
      LoadIdx_Stop [lv] = -1;    // must be initialized as <= LoadIdx_Start[lv]

//    skip levels with no patch (otherwise LoadIdx_Start and LoadIdx_Stop will be set incorrectly)
      if ( NPatchTotal[lv] == 0 )   continue;

      for (int t=0; t<NPatchTotal[lv]; t++)
      {
//       set LoadIdx_Start to "the first patch belonging to this rank"
         if (  LoadIdx_Start[lv] == -1  &&  LB_Index2Rank( lv, LBIdxList_EachLv[lv][t], CHECK_ON ) == MPI_Rank  )
            LoadIdx_Start[lv] = t;

//       set LoadIdx_Stop to "the last patch belonging to this rank + 1"
         if ( LoadIdx_Start[lv] != -1 )
         {
            if (  LB_Index2Rank( lv, LBIdxList_EachLv[lv][t], CHECK_ON ) > MPI_Rank  )
            {
               LoadIdx_Stop[lv] = t;
               break;
            }

//          rank owning the last patch needs to be treated separately
            else if ( t == NPatchTotal[lv] - 1 )
            {
               LoadIdx_Stop[lv] = NPatchTotal[lv];
               break;
            }
         }
      }

#     ifdef DEBUG_HDF5
      if ( LoadIdx_Start[lv]%8 != 0  &&  LoadIdx_Start[lv] != -1 )
         Aux_Error( ERROR_INFO, "LoadIdx_Start[%d] = %d --> %%8 != 0 !!\n", lv, LoadIdx_Start[lv]%8 );

      if ( LoadIdx_Stop [lv]%8 != 0  &&  LoadIdx_Stop[lv] != -1 )
         Aux_Error( ERROR_INFO, "LoadIdx_Stop [%d] = %d --> %%8 != 0 !!\n", lv, LoadIdx_Stop [lv]%8 );

      if (  ( LoadIdx_Start[lv] == -1 && LoadIdx_Stop[lv] != -1 )  ||
            ( LoadIdx_Start[lv] != -1 && LoadIdx_Stop[lv] == -1 )   )
         Aux_Error( ERROR_INFO, "LoadIdx_Start/Stop[%d] = %d/%d !!\n", lv, LoadIdx_Start[lv], LoadIdx_Stop[lv] );
#     endif
   } // for (int lv=0; lv<NLevel_RS; lv++)

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading load-balance index table ... done\n" );


#  else // #ifdef LOAD_BALANCE


// 2-3. son
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading son table ...\n" );

// allocate memory
   int *SonList_AllLv = new int [ NPatchAllLv ];

// load data
   H5_SetID_Son = H5Dopen( H5_FileID, "Tree/Son", H5P_DEFAULT );
   if ( H5_SetID_Son < 0 )    Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Son" );
   H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Son );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading son table ... done\n" );
#  endif // #ifdef LOAD_BALANCE ... else ...


// 2-4. number of particles in each patch
#  ifdef PARTICLE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading particle counts ...\n" );

   hid_t H5_SetID_NPar;

// allocate memory
   int *NParList_AllLv = new int [ NPatchAllLv ];

// load data
   if ( ReenablePar )   for (int t=0; t<NPatchAllLv; t++)   NParList_AllLv[t] = 0;
   else {
      H5_SetID_NPar = H5Dopen( H5_FileID, "Tree/NPar", H5P_DEFAULT );
      if ( H5_SetID_NPar < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/NPar" );
      H5_Status = H5Dread( H5_SetID_NPar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, NParList_AllLv );
      H5_Status = H5Dclose( H5_SetID_NPar );
   } // if ( ReenablePar ) ... else ...

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading particle counts ... done\n" );
#  endif // #ifdef PARTICLE


   H5_Status = H5Fclose( H5_FileID );


// 2-5. initialize particle variables
#  ifdef PARTICLE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing particle repository ...\n" );

// 2-5-1. get the total number of paticles in this rank
   long NParThisRank;

#  ifdef LOAD_BALANCE
   NParThisRank = 0;

   for (int lv=0; lv<NLevel_RS; lv++)
   for (int t=LoadIdx_Start[lv]; t<LoadIdx_Stop[lv]; t++)
   {
#     ifdef DEBUG_HDF5
      if ( t < 0  ||  t >= NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "incorrect load index (%d) !!\n", t );
#     endif

      NParThisRank += NParList_AllLv[ LBIdxList_EachLv_IdxTable[lv][t] + GID_LvStart[lv] ];
   }

#  ifdef DEBUG_HDF5
   long NParAllRank;
   MPI_Reduce( &NParThisRank, &NParAllRank, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   if ( MPI_Rank == 0  &&  NParAllRank != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles (%ld) != expect (%ld) !!\n",
                 NParAllRank, amr->Par->NPar_Active_AllRank );
#  endif

#  else // #ifdef LOAD_BALANCE

   NParThisRank = amr->Par->NPar_Active_AllRank;
#  endif // #ifdef LOAD_BALANCE ... else ...


// 2-5-2. initialize particle repository
   amr->Par->InitRepo( NParThisRank, MPI_NRank );

// reset the total number of particles to be zero
// --> so particle repository is pre-allocated, but it contains no active particle yet
// --> particles will be added later by calling Par->AddOneParticle()
   amr->Par->NPar_AcPlusInac = 0;
   amr->Par->NPar_Active     = 0;


// 2-5-3. calculate the starting global particle indices (i.e., GParID_Offset) for all patches
   long *GParID_Offset = new long [ NPatchAllLv ];

   GParID_Offset[0] = 0;
   for (int t=1; t<NPatchAllLv; t++)   GParID_Offset[t] = GParID_Offset[t-1] + NParList_AllLv[t-1];

#  ifdef DEBUG_HDF5
   if ( GParID_Offset[ NPatchAllLv-1 ] + NParList_AllLv[ NPatchAllLv-1 ] != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles (%ld) != expect (%ld) !!\n",
                 GParID_Offset[ NPatchAllLv-1 ] + NParList_AllLv[ NPatchAllLv-1 ], amr->Par->NPar_Active_AllRank );
#  endif


// 2-5-4. get the maximum number of particles in one patch and allocate an I/O buffer accordingly
   long MaxNParInOnePatch = 0;
   long *NewParList       = NULL;
   real_par **ParFltBuf   = NULL;
   long_par **ParIntBuf   = NULL;

   for (int t=0; t<NPatchAllLv; t++)   MaxNParInOnePatch = MAX( MaxNParInOnePatch, NParList_AllLv[t] );

   NewParList = new long [MaxNParInOnePatch];

// be careful about using ParFlt/IntBuf returned from Aux_AllocateArray2D, which is set to NULL if MaxNParInOnePatch == 0
// --> for example, accessing ParFlt/IntBuf[0...PAR_NATT_FLT/INT_STORED-1] will be illegal when MaxNParInOnePatch == 0
   Aux_AllocateArray2D( ParFltBuf, PAR_NATT_FLT_STORED, MaxNParInOnePatch );
   Aux_AllocateArray2D( ParIntBuf, PAR_NATT_INT_STORED, MaxNParInOnePatch );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing particle repository ... done\n" );
#  endif // #ifdef PARTICLE



// 3. load and allocate patches (load particles as well if PARTICLE is on)
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading patches and particles ...\n" );

   int NCompStore = NCOMP_TOTAL;

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
// do not load STUB field (the hybrid scheme always stores density and phase)
   NCompStore -= 1;
#  endif

#  ifdef LOAD_BALANCE
   const bool Recursive_No  = false;
#  else
   const bool Recursive_Yes = true;
   int  TRange_Min[3], TRange_Max[3];
#  endif

   char (*FieldName)[MAX_STRING] = new char [NCompStore][MAX_STRING];
   hsize_t H5_SetDims_Field[4], H5_MemDims_Field[4];
   hid_t   H5_SetID_Field[NCompStore], H5_MemID_Field, H5_SpaceID_Field, H5_GroupID_GridData;

#  ifdef MHD
   char (*FCMagName)[MAX_STRING] = new char [NCOMP_MAG][MAX_STRING];
   hsize_t H5_SetDims_FCMag[4], H5_MemDims_FCMag[4];
   hid_t   H5_SetID_FCMag[NCOMP_MAG], H5_MemID_FCMag[NCOMP_MAG], H5_SpaceID_FCMag[NCOMP_MAG];
#  else
   hid_t *H5_SetID_FCMag   = NULL;
   hid_t *H5_MemID_FCMag   = NULL;
   hid_t *H5_SpaceID_FCMag = NULL;
#  endif // #ifdef MHD ... else ...

#  ifdef PARTICLE
   char (*ParAttFltName)[MAX_STRING] = new char [PAR_NATT_FLT_STORED][MAX_STRING];
   char (*ParAttIntName)[MAX_STRING] = new char [PAR_NATT_INT_STORED][MAX_STRING];
   hsize_t H5_SetDims_ParData[1];
   hid_t   H5_SetID_ParFltData[PAR_NATT_FLT_STORED], H5_SetID_ParIntData[PAR_NATT_INT_STORED], H5_SpaceID_ParData, H5_GroupID_Particle;
#  else
// define useless variables when PARTICLE is off
   int       *NParList_AllLv      = NULL;
   real_par **ParFltBuf           = NULL;
   long_par **ParIntBuf           = NULL;
   long      *NewParList          = NULL;
   long      *GParID_Offset       = NULL;
   hid_t     *H5_SetID_ParFltData = NULL;
   hid_t     *H5_SetID_ParIntData = NULL;
   hid_t      H5_SpaceID_ParData  = NULL_INT;
   long       NParThisRank        = NULL_INT;
#  endif // #ifdef PARTICLE ... else ...


// 3-1. set the names of all grid fields and particle attributes
   for (int v=0; v<NCompStore; v++)       sprintf( FieldName[v], "%s", FieldLabel[v] );

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)        sprintf( FCMagName[v], "%s", MagLabel[v] );
#  endif

#  ifdef PARTICLE
// skip the last PAR_NATT_FLT/INT_UNSTORED attributes since we do not store them on disk
   for (int v=0; v<PAR_NATT_FLT_STORED; v++)  sprintf( ParAttFltName[v], "%s", ParAttFltLabel[v] );
   for (int v=0; v<PAR_NATT_INT_STORED; v++)  sprintf( ParAttIntName[v], "%s", ParAttIntLabel[v] );
#  endif


// 3-2. initialize relevant HDF5 objects
   H5_SetDims_Field[0] = NPatchAllLv;
   H5_SetDims_Field[1] = PS1;
   H5_SetDims_Field[2] = PS1;
   H5_SetDims_Field[3] = PS1;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

   H5_MemDims_Field[0] = 1;
   H5_MemDims_Field[1] = PS1;
   H5_MemDims_Field[2] = PS1;
   H5_MemDims_Field[3] = PS1;

   H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
   if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemID_Field" );

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
      H5_SetDims_FCMag[0] = NPatchAllLv;
      for (int t=1; t<4; t++)
      H5_SetDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_SpaceID_FCMag[v] = H5Screate_simple( 4, H5_SetDims_FCMag, NULL );
      if ( H5_SpaceID_FCMag[v] < 0 )
         Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n", "H5_SpaceID_FCMag", v );

      H5_MemDims_FCMag[0] = 1;
      for (int t=1; t<4; t++)
      H5_MemDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_MemID_FCMag[v] = H5Screate_simple( 4, H5_MemDims_FCMag, NULL );
      if ( H5_MemID_FCMag[v] < 0 )
         Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n", "H5_MemID_FCMag", v );
   }
#  endif // #ifdef MHD

#  ifdef PARTICLE
   if ( ! ReenablePar ) {
      H5_SetDims_ParData[0] = amr->Par->NPar_Active_AllRank;
      H5_SpaceID_ParData    = H5Screate_simple( 1, H5_SetDims_ParData, NULL );
      if ( H5_SpaceID_ParData < 0 )    Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_ParData" );
   }
#  endif


// load data with RESTART_LOAD_NRANK ranks at a time
   for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)
   {
      if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )
      {
//       3-3. open the target datasets just once
         H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
         if ( H5_FileID < 0 )
            Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
         if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );

         for (int v=0; v<NCompStore; v++)
         {
            H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );
            if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
         }

#        ifdef MHD
         for (int v=0; v<NCOMP_MAG; v++)
         {
            H5_SetID_FCMag[v] = H5Dopen( H5_GroupID_GridData, FCMagName[v], H5P_DEFAULT );
            if ( H5_SetID_FCMag[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FCMagName[v] );
         }
#        endif

#        ifdef PARTICLE
         if ( ! ReenablePar ) {
            H5_GroupID_Particle = H5Gopen( H5_FileID, "Particle", H5P_DEFAULT );
            if ( H5_GroupID_Particle < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Particle" );

            for (int v=0; v<PAR_NATT_FLT_STORED; v++)
            {
               H5_SetID_ParFltData[v] = H5Dopen( H5_GroupID_Particle, ParAttFltName[v], H5P_DEFAULT );
               if ( H5_SetID_ParFltData[v] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", ParAttFltName[v] );
            }
            for (int v=0; v<PAR_NATT_INT_STORED; v++)
            {
               H5_SetID_ParIntData[v] = H5Dopen( H5_GroupID_Particle, ParAttIntName[v], H5P_DEFAULT );
               if ( H5_SetID_ParIntData[v] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", ParAttIntName[v] );
            }
         } // if ( ! ReenablePar )
#        endif


//       3-4. begin to load data
//       3-4-1. load-balance data
#        ifdef LOAD_BALANCE
         int GID0;

         for (int lv=0; lv<NLevel_RS; lv++)
         {
            if ( MPI_Rank == TRanks )
            Aux_Message( stdout, "      Loading ranks %4d -- %4d, lv %2d ... ",
                         TRanks, MIN(TRanks+RESTART_LOAD_NRANK-1, MPI_NRank-1), lv );

//          loop over all target LBIdx
            for (int t=LoadIdx_Start[lv]; t<LoadIdx_Stop[lv]; t+=8)
            {
#              ifdef DEBUG_HDF5
               if ( t < 0  ||  t >= NPatchTotal[lv]  ||  t%8 != 0 )
                 Aux_Error( ERROR_INFO, "incorrect load index (%d) !!\n", t );
#              endif

//             make sure that we load patch from LocalID == 0
               GID0 = LBIdxList_EachLv_IdxTable[lv][t] - LBIdxList_EachLv_IdxTable[lv][t]%8 + GID_LvStart[lv];

               for (int GID=GID0; GID<GID0+8; GID++)
                  LoadOnePatch( H5_FileID, lv, GID, Recursive_No, NULL, CrList_AllLv,
                                H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field,
                                H5_SetID_FCMag, H5_SpaceID_FCMag, H5_MemID_FCMag,
                                NParList_AllLv, ParFltBuf, ParIntBuf, NewParList,
                                H5_SetID_ParFltData, H5_SetID_ParIntData, H5_SpaceID_ParData,
                                GParID_Offset, NParThisRank, FormatVersion_RS );
            }

//          check if LocalID matches corner
#           ifdef DEBUG_HDF5
            const int PatchScale = PS1*amr->scale[lv];
            int LocalID;

            for (int PID=0; PID<amr->num[lv]; PID++)
            {
               LocalID = PID%8;

               for (int d=0; d<3; d++)
               {
                  if (  TABLE_02( LocalID, 'x'+d, 0, 1 ) != (amr->patch[0][lv][PID]->corner[d]/PatchScale)%2  )
                  {
                     Output_Patch( lv, PID, 0, 0, 0, NULL );
                     Aux_Error( ERROR_INFO, "lv %d, PID %d, LocalID does not match patch corner (patch has been output) !!\n",
                                 lv, PID );
                  }
               }
            }
#           endif

            if ( MPI_Rank == TRanks )  Aux_Message( stdout, "done\n" );
         } // for (int lv=0; lv<NLevel_RS; lv++)


//       3-4.2. non load-balance data
#        else // #ifdef LOAD_BALANCE

//       set the target range of each rank
         for (int d=0; d<3; d++)
         {
            TRange_Min[d] = MPI_Rank_X[d]*NX0[d]*amr->scale[0];
            TRange_Max[d] = TRange_Min[d] + NX0[d]*amr->scale[0];
         }

//       loop over the corners of all root-level patches (rescale in advance if NLevel_RS != NLEVEL)
         const int TenPercent = MAX( NPatchTotal[0]/10, 1 );

         for (int GID=0; GID<NPatchTotal[0]; GID++)
         {
            if ( GID % TenPercent == 0 )
            Aux_Message( stdout, "      Rank %2d: %5.1lf%% completed ...\n",
                         MPI_Rank, 100.0*GID/NPatchTotal[0] );

//          load the info and data of the entire patch family recursively if the root patch is within the target range
            if (  CrList_AllLv[GID][0] >= TRange_Min[0]  &&  CrList_AllLv[GID][0] < TRange_Max[0]  &&
                  CrList_AllLv[GID][1] >= TRange_Min[1]  &&  CrList_AllLv[GID][1] < TRange_Max[1]  &&
                  CrList_AllLv[GID][2] >= TRange_Min[2]  &&  CrList_AllLv[GID][2] < TRange_Max[2]     )
               LoadOnePatch( H5_FileID, 0, GID, Recursive_Yes, SonList_AllLv, CrList_AllLv,
                             H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field,
                             H5_SetID_FCMag, H5_SpaceID_FCMag, H5_MemID_FCMag,
                             NParList_AllLv, ParFltBuf, ParIntBuf, NewParList,
                             H5_SetID_ParFltData, H5_SetID_ParIntData, H5_SpaceID_ParData,
                             GParID_Offset, NParThisRank, FormatVersion_RS );
         } // for (int GID=0; GID<NPatchTotal[0]; GID++)

#        endif // #ifdef LOAD_BALANCE ... else ...

//       free resource
         for (int v=0; v<NCompStore; v++)       H5_Status = H5Dclose( H5_SetID_Field[v] );
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG;   v++)      H5_Status = H5Dclose( H5_SetID_FCMag[v] );
#        endif
         H5_Status = H5Gclose( H5_GroupID_GridData );

#        ifdef PARTICLE
         if ( ! ReenablePar ) {
            for (int v=0; v<PAR_NATT_FLT_STORED; v++)  H5_Status = H5Dclose( H5_SetID_ParFltData[v] );
            for (int v=0; v<PAR_NATT_INT_STORED; v++)  H5_Status = H5Dclose( H5_SetID_ParIntData[v] );
            H5_Status = H5Gclose( H5_GroupID_Particle );
         }
#        endif

         H5_Status = H5Fclose( H5_FileID );
      } // if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)

// free HDF5 objects
   H5_Status = H5Sclose( H5_SpaceID_Field );
   H5_Status = H5Sclose( H5_MemID_Field );
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
      H5_Status = H5Sclose( H5_SpaceID_FCMag[v] );
      H5_Status = H5Sclose( H5_MemID_FCMag  [v] );
   }
#  endif
#  ifdef PARTICLE
   if ( ! ReenablePar )    H5_Status = H5Sclose( H5_SpaceID_ParData );
#  endif


// 3-5. record the number of real patches (and LB_IdxList_Real)
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = amr->num[lv];

#     ifdef LOAD_BALANCE
      if ( amr->LB->IdxList_Real         [lv] != NULL )   delete [] amr->LB->IdxList_Real         [lv];
      if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )   delete [] amr->LB->IdxList_Real_IdxTable[lv];

      amr->LB->IdxList_Real         [lv] = new long [ amr->NPatchComma[lv][1] ];
      amr->LB->IdxList_Real_IdxTable[lv] = new int  [ amr->NPatchComma[lv][1] ];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

      Mis_Heapsort( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );
#     endif

//    get the total number of real patches at all ranks
      Mis_GetTotalPatchNumber( lv );
   }


// 3-6. verify that all patches and particles are loaded
#  ifdef DEBUG_HDF5
   int NLoadPatch[NLevel_RS];
   MPI_Reduce( amr->num, NLoadPatch, NLevel_RS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   for (int lv=0; lv<NLevel_RS; lv++)
   {
      if ( NLoadPatch[lv] != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "Lv %d: total number of loaded patches (%d) != expectation (%d) !!\n",
                       lv, NLoadPatch[lv], NPatchTotal[lv] );
   }

#  ifdef PARTICLE
   if ( amr->Par->NPar_AcPlusInac != NParThisRank )
      Aux_Error( ERROR_INFO, "total number of particles in the repository (%ld) != expect (%ld) !!\n",
                 amr->Par->NPar_AcPlusInac, NParThisRank );
#  endif
#  endif // #ifdfe DEBUG_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading patches and particles ... done\n" );



// 4. close all HDF5 objects and free memory (especially for the variable-length string)

   delete [] FieldName;
   delete [] CrList_AllLv;
#  ifdef MHD
   delete [] FCMagName;
#  endif
#  ifdef LOAD_BALANCE
   delete [] LBIdxList_AllLv;
   for (int lv=0; lv<NLEVEL; lv++)  delete [] LBIdxList_EachLv_IdxTable[lv];
#  else
   delete [] SonList_AllLv;
#  endif
#  ifdef PARTICLE
   delete [] ParAttFltName;
   delete [] ParAttIntName;
   delete [] NParList_AllLv;
   delete [] GParID_Offset;
   delete [] NewParList;
   Aux_DeallocateArray2D( ParFltBuf );
   Aux_DeallocateArray2D( ParIntBuf );
#  endif



// 5-1. improve load balance
// ===================================================================================================================
#  ifdef LOAD_BALANCE

// no need to redistribute all patches again since we already did that when loading data from disks
// --> but note that we have not considerer load-balance weighting of particles yet

// we don't have enough information to calculate the load-balance weighting of particles when
// calling LB_Init_LoadBalance() for the first time
// --> for example, LB_EstimateWorkload_AllPatchGroup()->Par_CollectParticle2OneLevel()->Par_LB_CollectParticle2OneLevel()
//     needs amr->LB->IdxList_Real[], which will be constructed only AFTER calling LB_Init_LoadBalance()
// --> must disable particle weighting (by setting ParWeight==0.0) first

// must not reset load-balance variables (i.e., must adopt ResetLB_No) when calling LB_Init_LoadBalance() for the first time
// since we MUST NOT overwrite IdxList_Real[] and IdxList_Real_IdxList[] already set above
   const double ParWeight_Zero    = 0.0;
   const bool   Redistribute_Yes  = true;
   const bool   Redistribute_No   = false;
   const bool   SendGridData_Yes  = true;
   const bool   SendGridData_No   = false;
   const bool   ResetLB_Yes       = true;
   const bool   ResetLB_No        = false;
   const int    AllLv             = -1;

   LB_Init_LoadBalance( Redistribute_No,  SendGridData_No,  ParWeight_Zero,      ResetLB_No,  OPT__SORT_PATCH_BY_LBIDX,  AllLv );

// redistribute patches again if we want to take into account the load-balance weighting of particles
#  ifdef PARTICLE
   if ( amr->LB->Par_Weight > 0.0  &&  !ReenablePar )
   LB_Init_LoadBalance( Redistribute_Yes, SendGridData_Yes, amr->LB->Par_Weight, ResetLB_Yes, OPT__SORT_PATCH_BY_LBIDX,  AllLv );
#  endif



// 5-2. complete all levels for the case without LOAD_BALANCE
// ===================================================================================================================
#  else // #ifdef LOAD_BALANCE

   for (int lv=0; lv<NLEVEL; lv++)
   {
//    construct the relation "father <-> son"
      if ( lv > 0 )     FindFather( lv, 1 );

//    allocate the buffer patches
      Buf_AllocateBufferPatch( amr, lv );

//    set up the BaseP List
      if ( lv == 0 )    Init_RecordBasePatch();

//    set up the BounP_IDMap
      Buf_RecordBoundaryPatch( lv );

//    construct the sibling relation
      SiblingSearch( lv );

//    get the IDs of patches for sending and receiving data between neighbor ranks
      Buf_RecordExchangeDataPatchID( lv );

//    allocate flux arrays on level "lv-1"
      if ( lv > 0  &&  amr->WithFlux )    Flu_AllocateFluxArray( lv-1 );

//    allocate electric arrays on level "lv-1"
#     ifdef MHD
      if ( lv > 0  &&  amr->WithElectric )   MHD_AllocateElectricArray( lv-1 );
#     endif

//    get data for all buffer patches
      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_NO );

   } // for (int lv=0; lv<NLEVEL; lv++)

#  endif // #ifdef LOAD_BALANCE ... else ...


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ByRestart_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadOnePatch
// Description :  Allocate and load all fields (and particles if PARTICLE is on) for one patch
//
// Note        :  1. Both leaf and non-leaf patches will store data in the HDF5 output
//                   --> But only leaf patches have particles associated with them
//                2. If "Recursive == true", this function will be invoked recursively to find all children
//                   (and children's children, ...) patches
//
// Parameter   :  H5_FileID           : HDF5 file ID of the restart file
//                lv                  : Target level
//                GID                 : Target GID
//                Recursive           : Find all children (and childrens' children, ...) recuresively
//                SonList             : List of son indices
//                                      --> Set only when LOAD_BALANCE is not defined
//                CrList              : List of patch corners
//                H5_SetID_Field      : HDF5 dataset ID for cell-centered grid data
//                H5_SpaceID_Field    : HDF5 dataset dataspace ID for cell-centered grid data
//                H5_MemID_Field      : HDF5 memory dataspace ID for cell-centered grid data
//                H5_SetID_FCMag      : HDF5 dataset ID for face-centered magnetic field
//                H5_SpaceID_FCMag    : HDF5 dataset dataspace ID for face-centered magnetic field
//                H5_MemID_FCMag      : HDF5 memory dataspace ID for face-centered magnetic field
//                NParList            : List of particle counts
//                ParFlt/IntBuf       : I/O buffer for loading particle data from the disk
//                                      --> It must be preallocated with a size equal to the maximum number of
//                                          particles in one patch times the number of particles attributes
//                                          stored on disk
//                                      --> Be careful about using ParFlt/IntBuf, which is set to NULL if it has no elements
//                                          (because of the current implementation of Aux_AllocateArray2D)
//                                          --> For example, accessing ParFlt/IntBuf[0...PAR_NATT_FLT/INT_STORED-1] will be illegal when there
//                                              are no particles
//                NewParList          : Array to store the new particle indices
//                                      --> It must be preallocated with a size equal to the maximum number of
//                                          particles in one patch
//                H5_SetID_ParFltData : HDF5 dataset ID for particle floating-point data
//                H5_SetID_ParIntData : HDF5 dataset ID for particle integer        data
//                H5_SpaceID_ParData  : HDF5 dataset dataspace ID for particle data
//                GParID_Offset       : Starting global particle indices for all patches
//                NParThisRank        : Total number of particles in this rank (for check only)
//                FormatVersion       : HDF5 snapshot format version
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive,
                   const int *SonList, const int (*CrList)[3],
                   const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field,
                   const hid_t *H5_SetID_FCMag, const hid_t *H5_SpaceID_FCMag, const hid_t *H5_MemID_FCMag,
                   const int *NParList, real_par **ParFltBuf, long_par **ParIntBuf, long *NewParList,
                   const hid_t *H5_SetID_ParFltData, const hid_t *H5_SetID_ParIntData,
                   const hid_t H5_SpaceID_ParData, const long *GParID_Offset, const long NParThisRank,
                   const int FormatVersion )
{

   const bool WithData_Yes = true;

   hsize_t H5_Count_Field[4], H5_Offset_Field[4];
#  ifdef MHD
   hsize_t H5_Count_FCMag[4], H5_Offset_FCMag[4];
#  endif
   herr_t  H5_Status;
   int     SonGID0, PID;

   int NCompStore = NCOMP_TOTAL;

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
// do not load STUB field (the hybrid scheme always stores density and phase)
   NCompStore -= 1 ;
#  endif

// allocate patch
   amr->pnew( lv, CrList[GID][0], CrList[GID][1], CrList[GID][2], -1, WithData_Yes, WithData_Yes, WithData_Yes );

   PID = amr->num[lv] - 1;


// determine the subset of dataspace for grid data
   H5_Offset_Field[0] = GID;
   H5_Offset_Field[1] = 0;
   H5_Offset_Field[2] = 0;
   H5_Offset_Field[3] = 0;

   H5_Count_Field [0] = 1;
   H5_Count_Field [1] = PS1;
   H5_Count_Field [2] = PS1;
   H5_Count_Field [3] = PS1;

   H5_Status = H5Sselect_hyperslab( H5_SpaceID_Field, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
   if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the grid data !!\n" );


// load cell-centered intrinsic variables from disk
// --> excluding all derived variables such as gravitational potential and cell-centered B field
   for (int v=0; v<NCompStore; v++)
   {
      H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                           amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load a field variable (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   }


// convert phase/density to real and imaginary parts
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
      real Dens, Phas, Im, Re;

      for (int k=0; k<PS1; k++) {
      for (int j=0; j<PS1; j++) {
      for (int i=0; i<PS1; i++) {
         Dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
         Phas = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[PHAS][k][j][i];
         Re   = SQRT(Dens) * COS(Phas);
         Im   = SQRT(Dens) * SIN(Phas);

         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i] = Re;
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i] = Im;
      }}}
   }
#  endif


// load face-centered magnetic field from disk
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
//    determine the subset of dataspace for grid data
      H5_Offset_FCMag[0] = GID;
      H5_Offset_FCMag[1] = 0;
      H5_Offset_FCMag[2] = 0;
      H5_Offset_FCMag[3] = 0;

      H5_Count_FCMag [0] = 1;
      for (int t=1; t<4; t++)
      H5_Count_FCMag [t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_Status = H5Sselect_hyperslab( H5_SpaceID_FCMag[v], H5S_SELECT_SET, H5_Offset_FCMag, NULL, H5_Count_FCMag, NULL );
      if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the magnetic field %d !!\n", v );

//    load data
      H5_Status = H5Dread( H5_SetID_FCMag[v], H5T_GAMER_REAL, H5_MemID_FCMag[v], H5_SpaceID_FCMag[v], H5P_DEFAULT,
                           amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load magnetic field (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   } // for (int v=0; v<NCOMP_MAG; v++)
#  endif // #ifdef MHD


// load particle data
#  ifdef PARTICLE
   const int NParThisPatch = NParList[GID];

   hsize_t     H5_Offset_ParData[1], H5_Count_ParData[1], H5_MemDims_ParData[1];
   hid_t       H5_MemID_ParData;
   real_par    NewParAttFlt[PAR_NATT_FLT_TOTAL];
   long_par    NewParAttInt[PAR_NATT_INT_TOTAL];

   if ( NParThisPatch > 0 )
   {
//    check: only leaf patches can have particles (but note that SonList is NOT set for LOAD_BALANCE)
#     ifndef LOAD_BALANCE
      if ( SonList[GID] != -1 )
         Aux_Error( ERROR_INFO, "particles in a non-leaf patch (lv %d, PID %d, GID %d, SonPID %d, NPar %d) !!\n",
                    lv, PID, GID, SonList[GID], NParThisPatch );
#     endif

//    determine the memory space and the subset of dataspace for particle data
      H5_Offset_ParData [0] = GParID_Offset[GID];
      H5_Count_ParData  [0] = NParThisPatch;
      H5_MemDims_ParData[0] = NParThisPatch;

      H5_Status = H5Sselect_hyperslab( H5_SpaceID_ParData, H5S_SELECT_SET, H5_Offset_ParData, NULL, H5_Count_ParData, NULL );
      if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the particle data !!\n" );

      H5_MemID_ParData = H5Screate_simple( 1, H5_MemDims_ParData, NULL );
      if ( H5_MemID_ParData < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemID_ParData" );

//    load particle data from disk
      if ( FormatVersion < 2500 )
      {
         const int ParTypeIdx_old = 7;
         int skip_type = 0;
         for (int v=0; v<PAR_NATT_FLT_STORED+1; v++)
         {
//          using ParFltBuf[v] here is safe since it's NOT called when NParThisPatch == 0
            if ( v == ParTypeIdx_old )
            {
               real_par *ParType_Buf = new real_par [NParThisPatch];
               H5_Status = H5Dread( H5_SetID_ParIntData[PAR_TYPE], H5T_GAMER_REAL_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT,
                                    ParType_Buf );
               for (int p=0; p<NParThisPatch; p++)   ParIntBuf[PAR_TYPE][p] = (long_par)ParType_Buf[p];
               delete [] ParType_Buf;
               skip_type = 1;
            }
            else
            {
               H5_Status = H5Dread( H5_SetID_ParFltData[v-skip_type], H5T_GAMER_REAL_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT,
                                    ParFltBuf[v-skip_type] );
            }
            if ( H5_Status < 0 )
               Aux_Error( ERROR_INFO, "failed to load a particle floating-point attribute (lv %d, GID %d, v %d) !!\n", lv, GID, v );
         }
      } // if ( FormatVersion < 2500 )
      else
      {
         for (int v=0; v<PAR_NATT_FLT_STORED; v++)
         {
//          using ParFltBuf[v] here is safe since it's NOT called when NParThisPatch == 0
            H5_Status = H5Dread( H5_SetID_ParFltData[v], H5T_GAMER_REAL_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT,
                                 ParFltBuf[v] );
            if ( H5_Status < 0 )
               Aux_Error( ERROR_INFO, "failed to load a particle floating-point attribute (lv %d, GID %d, v %d) !!\n", lv, GID, v );
         }
         for (int v=0; v<PAR_NATT_INT_STORED; v++)
         {
//          using ParIntBuf[v] here is safe since it's NOT called when NParThisPatch == 0
            H5_Status = H5Dread( H5_SetID_ParIntData[v], H5T_GAMER_LONG_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT,
                                 ParIntBuf[v] );
            if ( H5_Status < 0 )
               Aux_Error( ERROR_INFO, "failed to load a particle integer attribute (lv %d, GID %d, v %d) !!\n", lv, GID, v );
         }
      } // if ( FormatVersion < 2500 ) ... else ...

//    store particles to the particle repository (one particle at a time)
      NewParAttFlt[PAR_TIME] = Time[0];   // all particles are assumed to be synchronized with the base level

      for (int p=0; p<NParThisPatch; p++)
      {
//       skip the last PAR_NATT_FLT/INT_UNSTORED attributes since we do not store them on disk
         for (int v=0; v<PAR_NATT_FLT_STORED; v++)  NewParAttFlt[v] = ParFltBuf[v][p];
         for (int v=0; v<PAR_NATT_INT_STORED; v++)  NewParAttInt[v] = ParIntBuf[v][p];

         NewParList[p] = amr->Par->AddOneParticle( NewParAttFlt, NewParAttInt );

//       check
         if ( NewParList[p] >= NParThisRank )
            Aux_Error( ERROR_INFO, "New particle ID (%ld) >= maximum allowed value (%ld) !!\n",
                       NewParList[p], NParThisRank );
      } // for (int p=0; p<NParThisPatch )

//    link particles to this patch
      const long_par *PType = amr->Par->Type;
#     ifdef DEBUG_PARTICLE
      const real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
      char Comment[MAX_STRING];
      sprintf( Comment, "%s, lv %d, PID %d, GID %d, NPar %d", __FUNCTION__, lv, PID, GID, NParThisPatch );
      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParList, &amr->Par->NPar_Lv[lv],
                                           PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#     else
      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParList, &amr->Par->NPar_Lv[lv],
                                           PType );
#     endif

//    free resource
      H5_Status = H5Sclose( H5_MemID_ParData );

   } // if ( NParList[GID] > 0 )
#  endif // #ifdef PARTICLE


// enter the next level
   if ( Recursive )
   {
      if ( SonList == NULL )  Aux_Error( ERROR_INFO, "SonList == NULL (lv %d, GID %d) !!\n", lv, GID );

      SonGID0 = SonList[GID];

      if ( SonGID0 != -1 )
      {
         for (int SonGID=SonGID0; SonGID<SonGID0+8; SonGID++)
            LoadOnePatch( H5_FileID, lv+1, SonGID, Recursive, SonList, CrList,
                          H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field,
                          H5_SetID_FCMag, H5_SpaceID_FCMag, H5_MemID_FCMag,
                          NParList, ParFltBuf, ParIntBuf, NewParList,
                          H5_SetID_ParFltData, H5_SetID_ParIntData, H5_SpaceID_ParData,
                          GParID_Offset, NParThisRank, FormatVersion );
      }
   } // if ( Recursive )

} // FUNCTION : LoadOnePatch



//-------------------------------------------------------------------------------------------------------
// Function    :  GetFieldData
// Description :  Get the target field data from the restart file
//
// Note        :  1. The memory size of `Ptr` must not smaller than the data size in the restart file
//                2. `Ptr` must be NULL to allocate the memory
//                3. The allocated memory must be closed manually
//
// Parameter   :  FileName         : Restart file name
//                H5_SetID_Target  : HDF5 dataset  ID of the target compound variable
//                H5_TypeID_Target : HDF5 datatype ID of the target compound variable
//                Ptr              : The target field pointer to be stored
//                Allocate         : Whether to allocate the memory for `Ptr` or not
//
// Return      :  H5_FieldSize : the size of the target field
//-------------------------------------------------------------------------------------------------------
template <typename T>
size_t GetFieldData( const char *FieldName, const hid_t H5_SetID_Target, const hid_t H5_TypeID_Target,
                     T *Ptr, const bool Allocate )
{

   if ( Ptr != NULL  &&  Allocate )
      Aux_Error( ERROR_INFO, "Pointer must be NULL to allocate memory !!\n" );

   int    H5_FieldIdx;
   size_t H5_FieldSize;
   hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
   hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
   herr_t H5_Status;

// get field index in file
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldName );

   if ( H5_FieldIdx < 0 )
      Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the restart file !!\n", FieldName );

// load
   H5_TypeID_Field = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
   H5_FieldSize    = H5Tget_size( H5_TypeID_Field );

// allocate
   if ( Allocate )   Ptr = new T [ (int)H5_FieldSize ];

   H5_TypeID_Load = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
   H5_Status      = H5Tinsert( H5_TypeID_Load, FieldName, 0, H5_TypeID_Field );

   H5_Status      = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)Ptr );
   if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", FieldName );

   H5_Status = H5Tclose( H5_TypeID_Field );
   H5_Status = H5Tclose( H5_TypeID_Load  );

   return H5_FieldSize;

} // FUNCTION : GetFieldData



#endif // #ifdef SUPPORT_HDF5
