#ifdef SUPPORT_HDF5

#include "Copyright.h"
#include "GAMER.h"
#include "HDF5_Typedef.h"
#include <typeinfo>

void FillIn_Makefile (  Makefile_t &Makefile  );
void FillIn_SymConst (  SymConst_t &SymConst  );
void FillIn_InputPara( InputPara_t &InputPara );

template <typename T>
static herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                         const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                         const T *ComprPtr, const int NCompr, const bool Fatal_Compr );
static void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive, const int *SonList,
                          const int (*CrList)[3], const hid_t *H5_SetID, const hid_t H5_SpaceID, const hid_t H5_MemID );
static void Check_Makefile ( const char *FileName );
static void Check_SymConst ( const char *FileName );
static void Check_InputPara( const char *FileName );
static void ResetParameter( const char *FileName, double *EndT, long *EndStep );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Restart_HDF5
// Description :  Reload a previous HDF5 output as the initial condition
//
// Note        :  1. "OPT__RESTART_HEADER == RESTART_HEADER_CHECK"   
//                   --> Check if the parameters loaded from the RESTART file are consistent with the 
//                       parameters loaded from the Input__Parameter file
//
//                   "OPT__RESTART_HEADER == RESTART_HEADER_SKIP"
//                   --> Skip the header information in the RESTART file
//
//                2. This function will be invokde by "Init_Restart" automatically if the restart file
//                   is determined to a HDF5 file
//                3. Only work for format version >= 2100
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
void Init_Restart_HDF5( const char *FileName )
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

   MPI_Barrier( MPI_COMM_WORLD );


// 1. load the simulation info
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ...\n" );

   const bool    Fatal  = true;
   const bool NonFatal  = false;
   const int  Model     = MODEL;
   const int  PatchSize = PATCH_SIZE;
#  ifdef FLOAT8
   const int  Float8    = 1;
#  else
   const int  Float8    = 0;
#  endif
#  ifdef GRAVITY  
   const int  Gravity   = 1;
#  else
   const int  Gravity   = 0;
#  endif
#  ifdef PARTICLE
   const int  Particle  = 1;
#  else
   const int  Particle  = 0;
#  endif

   KeyInfo_t KeyInfo;

   hid_t  H5_FileID, H5_SetID_KeyInfo, H5_TypeID_KeyInfo, H5_SetID_Cr;
#  ifdef LOAD_BALANCE
   hid_t  H5_SetID_LBIdx;
#  else
   hid_t  H5_SetID_Son;
#  endif
   herr_t H5_Status;
   int    NLvRescale, NPatchAllLv, GID_LvStart[NLEVEL];
   int   *NullPtr = NULL;


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


// 1-3. load all target fields in KeyInfo one-by-one (by all ranks)
   LoadField( "FormatVersion",  &KeyInfo.FormatVersion,  H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );

// format version for HDF5 output
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "      The format version of the HDF5 RESTART file = %ld\n", KeyInfo.FormatVersion );

      if ( KeyInfo.FormatVersion < 2100 )
         Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2100) !!\n" );
   }

   MPI_Barrier( MPI_COMM_WORLD );

   LoadField( "Model",          &KeyInfo.Model,          H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal, &Model,          1,    Fatal );
   LoadField( "Float8",         &KeyInfo.Float8,         H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal, &Float8,         1,    Fatal );
   LoadField( "Gravity",        &KeyInfo.Gravity,        H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal, &Gravity,        1,    Fatal );
   LoadField( "Particle",       &KeyInfo.Particle,       H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal, &Particle,       1,    Fatal );
   LoadField( "NLevel",         &KeyInfo.NLevel,         H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
   LoadField( "PatchSize",      &KeyInfo.PatchSize,      H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal, &PatchSize,      1,    Fatal );

// runtime NLEVEL must be >= loaded NLEVEL
   if      ( KeyInfo.NLevel > NLEVEL )
      Aux_Error( ERROR_INFO, "%s : RESTART file (%d) > runtime (%d) !!\n", 
                 "NLEVEL", KeyInfo.NLevel, NLEVEL );
   else
   {
      NLvRescale = 1 << ( NLEVEL - KeyInfo.NLevel );

      if ( MPI_Rank == 0  &&  NLvRescale != 1 )
         Aux_Message( stderr, "WARNING : the NLEVEL rescale factor is set to %d\n", NLvRescale );
   }

   MPI_Barrier( MPI_COMM_WORLD );

   LoadField( "DumpID",             &KeyInfo.DumpID,             H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
#  ifdef GRAVITY
   LoadField( "OutputPot",          &KeyInfo.OutputPot,          H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
#  endif
   LoadField( "NX0",                 KeyInfo.NX0,                H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NX0_TOT,        3,    Fatal );
   LoadField( "BoxScale",            KeyInfo.BoxScale,           H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
   LoadField( "NPatch",              KeyInfo.NPatch,             H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
   LoadField( "CellScale",           KeyInfo.CellScale,          H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );

   LoadField( "Step",               &KeyInfo.Step,               H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
   LoadField( "AdvanceCounter",      KeyInfo.AdvanceCounter,     H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );

   LoadField( "BoxSize",             KeyInfo.BoxSize,            H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  amr->BoxSize,   3,    Fatal );
   LoadField( "Time",                KeyInfo.Time,               H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
   LoadField( "CellSize",            KeyInfo.CellSize,           H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
#  ifdef GRAVITY
   LoadField( "AveDens",            &KeyInfo.AveDens,            H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
#  endif

   LoadField( "CodeVersion",        &KeyInfo.CodeVersion,        H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );
   LoadField( "DumpWallTime",       &KeyInfo.DumpWallTime,       H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );


// 1-4. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Fclose( H5_FileID );


// 1-5. set internal parameters
   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   {
      Time          [lv] = KeyInfo.Time          [lv];
      NPatchTotal   [lv] = KeyInfo.NPatch        [lv];
      AdvanceCounter[lv] = KeyInfo.AdvanceCounter[lv];
   }

   Step       = KeyInfo.Step;
#  ifdef GRAVITY
   AveDensity = KeyInfo.AveDens;
#  endif


// 1-6. set parameters in levels that do not exist in the input file
   for (int lv=KeyInfo.NLevel; lv<NLEVEL; lv++)
   {
      Time          [lv] = 0.0;
      NPatchTotal   [lv] = 0;
      AdvanceCounter[lv] = 0;
   }


// 1-7. set Flu(Pot)SgTime
   for (int lv=0; lv<NLEVEL; lv++)
   {
      amr->FluSgTime[lv][ amr->FluSg[lv] ] = Time[lv];
#     ifdef GRAVITY
      amr->PotSgTime[lv][ amr->PotSg[lv] ] = Time[lv];
#     endif
   }


// 1-8. set the next dump ID
   if ( INIT_DUMPID < 0 )  DumpID = KeyInfo.DumpID + 1;
   else                    DumpID = INIT_DUMPID;


// 1-9. reset parameters from the restart file
   ResetParameter( FileName, &END_T, &END_STEP );


// 1-10. check all other simulation information (by rank 0 only)
   if ( OPT__RESTART_HEADER  &&  MPI_Rank == 0 )
   {
      Check_Makefile ( FileName );
      Check_SymConst ( FileName );
      Check_InputPara( FileName );
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

// rescale the loaded corner (necessary when KeyInfo.NLevel != NLEVEL)
   if ( NLvRescale != 1 )
   for (int GID=0; GID<NPatchAllLv; GID++)
   for (int d=0; d<3; d++)
      CrList_AllLv[GID][d] *= NLvRescale;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading corner table ... done\n" );


// 2-2. LBIdx
#  ifdef LOAD_BALANCE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading load-balance index table ...\n" );

   const bool InputLBIdxList_Yes = true;

   long *LBIdxList_AllLv;
   long *LBIdxList_EachLv         [NLEVEL];
   int  *LBIdxList_EachLv_IdxTable[NLEVEL];
   int   LoadIdx_Start[NLEVEL], LoadIdx_Stop[NLEVEL];

// 2-2-1. allocate memory
   LBIdxList_AllLv = new long [ NPatchAllLv ];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( lv < KeyInfo.NLevel )   
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


   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   {
//    2-2-3. sort the LBIdx list at each level and set the load-balance cut points
#     if ( LOAD_BALANCE != HILBERT )
      if ( NLvRescale != 1 )
      Aux_Message( stderr, "WARNING : please make sure that the patch LBIdx doesn't change when NLvRescale != 1 !!\n" );
#     endif

//    Mis_Heapsort must be called before LB_SetCutPoint in order to get LBIdxList_EachLv_IdxTable
//    (since LB_SetCutPoint will sort LBIdxList_EachLv as well)
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_EachLv[lv], LBIdxList_EachLv_IdxTable[lv] );

      LB_SetCutPoint( lv, amr->LB->CutPoint[lv], InputLBIdxList_Yes, LBIdxList_EachLv[lv] );


//    2-2-4. get the target LBIdx range of each rank
      LoadIdx_Start[lv] = -1;    // -1 --> indicate that it has not been set
      LoadIdx_Stop [lv] = -1;    // must be initialized as <= LoadIdx_Start[lv]

//    skip levels with no patch (otherwise LoadIdx_Start and LoadIdx_Stop will be set incorrectly)
      if ( NPatchTotal[lv] == 0 )   continue;

      for (int t=0; t<NPatchTotal[lv]; t++)
      {
         if (  LoadIdx_Start[lv] == -1  &&  LB_Index2Rank( lv, LBIdxList_EachLv[lv][t], CHECK_ON ) == MPI_Rank  )
            LoadIdx_Start[lv] = t;

         if (                               LB_Index2Rank( lv, LBIdxList_EachLv[lv][t], CHECK_ON ) >  MPI_Rank  )
         {
            LoadIdx_Stop [lv] = t;
            break;
         }
      }

//    take care of the last rank with patches (ranks without any patch will have Start=Stop=-1, which is fine)
      if ( LoadIdx_Start[lv] != -1  &&  LoadIdx_Stop[lv] == -1 )  LoadIdx_Stop[lv] = NPatchTotal[lv];

#     ifdef DEBUG_HDF5
      if ( LoadIdx_Start[lv]%8 != 0  &&  LoadIdx_Start[lv] != -1 )
         Aux_Error( ERROR_INFO, "LoadIdx_Start[%d] = %d --> %%8 != 0 !!\n", lv, LoadIdx_Start[lv]%8 );

      if ( LoadIdx_Stop [lv]%8 != 0  &&  LoadIdx_Stop[lv] != -1 )
         Aux_Error( ERROR_INFO, "LoadIdx_Stop [%d] = %d --> %%8 != 0 !!\n", lv, LoadIdx_Stop [lv]%8 );

      if (  ( LoadIdx_Start[lv] == -1 && LoadIdx_Stop[lv] != -1 )  ||
            ( LoadIdx_Start[lv] != -1 && LoadIdx_Stop[lv] == -1 )   )
         Aux_Error( ERROR_INFO, "LoadIdx_Start/Stop[%d] =  %d/%d !!\n", lv, LoadIdx_Start[lv], LoadIdx_Stop[lv] );
#     endif
   } // for (int lv=0; lv<KeyInfo.NLevel; lv++)

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

   H5_Status = H5Fclose( H5_FileID );


// 3. load and allocate patches (one rank at a time)
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading patches ...\n" );

#  ifdef LOAD_BALANCE
   const bool Recursive_No  = false;
#  else
   const bool Recursive_Yes = true;
   int  TRange_Min[3], TRange_Max[3];
#  endif

   char FieldName[NCOMP][100];

   hsize_t H5_SetDims_Field[4], H5_MemDims_Field[4];
   hid_t   H5_SetID_Field[NCOMP], H5_MemID_Field, H5_SpaceID_Field, H5_GroupID_Data;


// 3-1. set the field names
#  if   ( MODEL == HYDRO )
   sprintf( FieldName[DENS], "Dens" );
   sprintf( FieldName[MOMX], "MomX" );
   sprintf( FieldName[MOMY], "MomY" );
   sprintf( FieldName[MOMZ], "MomZ" );
   sprintf( FieldName[ENGY], "Engy" );

#  elif ( MODEL == ELBDM )
   sprintf( FieldName[DENS], "Dens" );
   sprintf( FieldName[REAL], "Real" );
   sprintf( FieldName[IMAG], "Imag" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif


// 3-2. initialize relevant HDF5 objects
   H5_SetDims_Field[0] = NPatchAllLv;
   H5_SetDims_Field[1] = PATCH_SIZE;
   H5_SetDims_Field[2] = PATCH_SIZE;
   H5_SetDims_Field[3] = PATCH_SIZE;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

   H5_MemDims_Field[0] = 1;
   H5_MemDims_Field[1] = PATCH_SIZE;
   H5_MemDims_Field[2] = PATCH_SIZE;
   H5_MemDims_Field[3] = PATCH_SIZE;

   H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
   if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_Field" );


// load data in one rank at a time
   for (int TRank=0; TRank<MPI_NRank; TRank++)
   {
      if ( MPI_Rank == TRank )
      {
//       3-3. open the target datasets just once
         H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
         if ( H5_FileID < 0 )
            Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_Data = H5Gopen( H5_FileID, "Data", H5P_DEFAULT );
         if ( H5_GroupID_Data < 0 )    Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Data" );

         for (int v=0; v<NCOMP; v++)
         {
            H5_SetID_Field[v] = H5Dopen( H5_GroupID_Data, FieldName[v], H5P_DEFAULT );
            if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
         }


//       3-4. begin to load data
//       3-4-1. load-balance data
#        ifdef LOAD_BALANCE
         int GID0;

         for (int lv=0; lv<KeyInfo.NLevel; lv++)
         {
            Aux_Message( stdout, "      Loading rank %2d, lv %2d ... ", MPI_Rank, lv );

//          loop over all targeted LBIdx
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
                                H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field );
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
                     Output_Patch( lv, PID, 0, 0, NULL );
                     Aux_Error( ERROR_INFO, "lv %d, PID %d, LocalID does not match patch corner (patch has been output) !!\n",
                                 lv, PID );
                  }
               }
            }
#           endif

            Aux_Message( stdout, "done\n" );
         } // for (int lv=0; lv<KeyInfo.NLevel; lv++)


//       3-4.2. non load-balance data
#        else // #ifdef LOAD_BALANCE

//       set the target range of each rank
         for (int d=0; d<3; d++)
         {
            TRange_Min[d] = MPI_Rank_X[d]*NX0[d]*amr->scale[0];
            TRange_Max[d] = TRange_Min[d] + NX0[d]*amr->scale[0];
         }

//       loop over the corners of all root-level patches (rescale in advance if KeyInfo.NLevel != NLEVEL)
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
                             H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field );
         } // for (int GID=0; GID<NPatchTotal[0]; GID++)

#        endif // #ifdef LOAD_BALANCE ... else ...

         for (int v=0; v<NCOMP; v++)
         H5_Status = H5Dclose( H5_SetID_Field[v] );
         H5_Status = H5Gclose( H5_GroupID_Data );
         H5_Status = H5Fclose( H5_FileID );
      } // if ( MyRank == TargetRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)

// free HDF5 objects
   H5_Status = H5Sclose( H5_SpaceID_Field );
   H5_Status = H5Sclose( H5_MemID_Field );


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
   }


// 3-6. verify that all patches are loaded
#  ifdef DEBUG_HDF5
   int NLoadPatch[KeyInfo.NLevel];
   MPI_Reduce( amr->num, NLoadPatch, KeyInfo.NLevel, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   {
      if ( NLoadPatch[lv] != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "Lv %d: total number of loaded patches (%d) != expectation (%d) !!\n",
                       lv, NLoadPatch[lv], NPatchTotal[lv] );
   }
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading patches ... done\n" );



// 4. close all HDF5 objects and free memory (especially for the variable-length string)
   free ( KeyInfo.CodeVersion );
   free ( KeyInfo.DumpWallTime );

   if ( CrList_AllLv != NULL )   delete [] CrList_AllLv;
#  ifdef LOAD_BALANCE
   if ( LBIdxList_AllLv != NULL )   delete [] LBIdxList_AllLv;
   for (int lv=0; lv<NLEVEL; lv++)
      if ( LBIdxList_EachLv_IdxTable[lv] != NULL )   delete [] LBIdxList_EachLv_IdxTable[lv];
#  else
   if ( SonList_AllLv != NULL )  delete [] SonList_AllLv;
#  endif



// 5. complete the tree structure (only necessary when LOAD_BALANCE is OFF)
#  ifndef LOAD_BALANCE
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

//    allocate the flux arrays at the level "lv-1"
      if ( lv > 0  &&  amr->WithFlux )    Flu_AllocateFluxArray( lv-1 );

//    get data for all buffer patches
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _FLU, Flu_ParaBuf, USELB_NO );

   } // for (int lv=0; lv<NLEVEL; lv++)
#  endif // #ifndef LOAD_BALANCE


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ ); 

} // FUNCTION : Init_Restart_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadField
// Description :  Load a single field from the input compound dataset
//
// Note        :  1. This function works for arbitary datatype (int, float, char, 1D array ...)
//                2. Memory must be allocated for FieldPtr in advance with sufficent size (except for "char *")
//                3. For loading a string, which has (type(FieldPtr) = (char *)), the memory must be freed
//                   manually by calling "free()"
//                4. It can also compare the loaded variables (FieldPtr) with the reference value (ComprPtr)
//                   (perform comparison only if "NCompr > 0")
//                   --> Please make sure that "FieldPtr" and "ComprPtr" point to the same type since we 
//                       use the type of "ComprPtr" to typecast "FieldPtr"
//
// Parameter   :  FieldName         : Name of the target field 
//                FieldPtr          : Pointer to store the retrieved data
//                H5_SetID_Target   : HDF5 dataset  ID of the target compound variable
//                H5_TypeID_Target  : HDF5 datatype ID of the target compound variable
//                Fatal_Nonexist    : Whether or not the nonexistence of the target field is fatal
//                                    --> true  : terminate the program     if the target field cannot be found
//                                        false : display a warning message if the target field cannot be found
//                ComprPtr          : Pointer to store the reference values for comparison
//                NCompr            : Number of elements to be compared
//                Fatal_Compr       : Whether or not the comparison result is fatal
//                                    --> true  : terminate the program     if "FieldPtr[X] != ComprPtr[X]"
//                                        false : display a warning message if "FieldPtr[X] != ComprPtr[X]" 
//
// Return      :  Success/fail <-> 0/-1
//-------------------------------------------------------------------------------------------------------
template <typename T>
herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                  const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                  const T *ComprPtr, const int NCompr, const bool Fatal_Compr )
{

#  ifdef DEBUG_HDF5
   if ( NCompr > 0  &&  ComprPtr == NULL )
      Aux_Error( ERROR_INFO, "ComprPtr == NULL for NCompr = %d > 0 !!\n", NCompr );
#  endif


   bool   Check_Pass = true;
   int    H5_FieldIdx;
   size_t H5_FieldSize;
   hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
   hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
   herr_t H5_Status;


// load
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldName );

   if ( H5_FieldIdx >= 0 )
   {
      H5_TypeID_Field  = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
      H5_FieldSize     = H5Tget_size( H5_TypeID_Field );

      H5_TypeID_Load   = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
      H5_Status        = H5Tinsert( H5_TypeID_Load, FieldName, 0, H5_TypeID_Field );

      H5_Status        = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, FieldPtr );
      if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", FieldName );

      H5_Status        = H5Tclose( H5_TypeID_Field );
      H5_Status        = H5Tclose( H5_TypeID_Load  );
   } // if ( H5_FieldIdx >= 0 )

   else
   {
      if ( Fatal_Nonexist )
         Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the restart file !!\n", FieldName );

      else if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : target field \"%s\" does not exist in the restart file !!\n", FieldName );

      return -1;
   } // if ( H5_FieldIdx >= 0 ) ... else ...


// comparison
   char ArrayIdx[10];

   for (int t=0; t<NCompr; t++)
   {
      if ( NCompr == 1 )   sprintf( ArrayIdx, "%s", "" );
      else                 sprintf( ArrayIdx, "[%d]", t );

      if ( ((T*)FieldPtr)[t] != ComprPtr[t] )
      {
         if (  typeid(T) == typeid( int)  ||  typeid(T) == typeid( long)  ||
               typeid(T) == typeid(uint)  ||  typeid(T) == typeid(ulong)    )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : RESTART file (%ld) != runtime (%ld) !!\n", 
                          FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               return -2;
            }

            else
            {
               if ( MPI_Rank == 0 )
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%ld) != runtime (%ld) !!\n", 
                            FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else if (  typeid(T) == typeid(float)  ||  typeid(T) == typeid(double)  )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                          FieldName, ArrayIdx,  ((T*)FieldPtr)[t], ComprPtr[t] );
               return -2;
            }

            else
            {
               if ( MPI_Rank == 0 )
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                            FieldName, ArrayIdx, ((T*)FieldPtr)[t], ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else
         {
            if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : \"%s%s\" : unsupported data type !!\n",
                         FieldName, ArrayIdx );
            Check_Pass = false;
         }

      } // if ( ((T*)FieldPtr)[t] != ComprPtr[t] )
   } // for (int t=0; t<NCompr; t++)

   return (Check_Pass) ? 0 : -3;

} // FUNCTION : LoadField



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadOnePatch
// Description :  Allocate and load all fields for one patch
//
// Note        :  1. Both leaf and non-leaf patches will store data in the HDF5 output
//                2. If "Recursive == true", this function will be invoked recursively to find all children
//                   (and children's children, ...) patches
//
// Parameter   :  H5_FileID   : HDF5 file ID of the restart file
//                lv          : Target level
//                GID         : Target GID
//                Recursive   : Find all children (and childrens' children, ...) recuresively
//                SonList     : List of son indices
//                CrList      : List of patch corners
//                H5_SetID    : HDF5 dataset ID for NCOMP fields
//                H5_SpaceID  : HDF5 dataset dataspace ID
//                H5_MemID    : HDF5 memory dataspace ID
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive, const int *SonList,
                   const int (*CrList)[3], const hid_t *H5_SetID, const hid_t H5_SpaceID, const hid_t H5_MemID )
{

   const bool WithData_Yes = true;

   hsize_t H5_Count_Field[4], H5_Offset_Field[4];
   herr_t  H5_Status;
   int     SonGID0, PID;

// allocate patch
   amr->pnew( lv, CrList[GID][0], CrList[GID][1], CrList[GID][2], -1, WithData_Yes, WithData_Yes );

   PID = amr->num[lv] - 1;


// determine the subset of the dataspace
   H5_Offset_Field[0] = GID;
   H5_Offset_Field[1] = 0;
   H5_Offset_Field[2] = 0;
   H5_Offset_Field[3] = 0;

   H5_Count_Field [0] = 1;
   H5_Count_Field [1] = PATCH_SIZE;
   H5_Count_Field [2] = PATCH_SIZE;
   H5_Count_Field [3] = PATCH_SIZE;

   H5_Status = H5Sselect_hyperslab( H5_SpaceID, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
   if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab !!\n" );


// load the field data (potential data, if presented, are ignored and will be recalculated)
   for (int v=0; v<NCOMP; v++)
   {
      H5_Status = H5Dread( H5_SetID[v], H5T_GAMER_REAL, H5_MemID, H5_SpaceID, H5P_DEFAULT,
                           amr->patch[0][lv][PID]->fluid[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load a field variable (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   }


// enter the next level
   if ( Recursive )
   {
      if ( SonList == NULL )  Aux_Error( ERROR_INFO, "SonList == NULL (lv %d, GID %d) !!\n", lv, GID );

      SonGID0 = SonList[GID];

      if ( SonGID0 != -1 )
      {
         for (int SonGID=SonGID0; SonGID<SonGID0+8; SonGID++)
            LoadOnePatch( H5_FileID, lv+1, SonGID, Recursive, SonList, CrList,
                          H5_SetID, H5_SpaceID, H5_MemID );
      }
   }

} // FUNCTION : LoadOnePatch



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Makefile
// Description :  Load and compare the Makefile_t structure (runtime vs. restart file)
//
// Note        :  1. Unnecessary comparisons are avoided
//                2. Program may be terminated if any fatal difference is found
//
// Parameter   :  FileName : Restart file name
//-------------------------------------------------------------------------------------------------------
void Check_Makefile( const char *FileName )
{

   const bool    Fatal = true;
   const bool NonFatal = false;

   herr_t Status;


// 1. fill in the runtime parameters
   Makefile_t RT;    // RT = RunTime
   FillIn_Makefile( RT );


// 2. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 3. load the Makefile dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/Makefile", H5P_DEFAULT );  // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/Makefile" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 4. load and compare all target fields one-by-one
   Makefile_t RS;    // RS = ReStart

   LoadField( "Model",              &RS.Model,              SID, TID, NonFatal, &RT.Model,               1,    Fatal );
   LoadField( "Gravity",            &RS.Gravity,            SID, TID, NonFatal, &RT.Gravity,             1,    Fatal );
   LoadField( "IndividualDt",       &RS.IndividualDt,       SID, TID, NonFatal, &RT.IndividualDt,        1, NonFatal );
   LoadField( "Comoving",           &RS.Comoving,           SID, TID, NonFatal, &RT.Comoving,            1,    Fatal );
   LoadField( "Particle",           &RS.Particle,           SID, TID, NonFatal, &RT.Particle,            1,    Fatal );

   LoadField( "UseGPU",             &RS.UseGPU,             SID, TID, NonFatal, &RT.UseGPU,              1, NonFatal );
   LoadField( "GAMER_Optimization", &RS.GAMER_Optimization, SID, TID, NonFatal, &RT.GAMER_Optimization,  1, NonFatal );
   LoadField( "GAMER_Debug",        &RS.GAMER_Debug,        SID, TID, NonFatal, &RT.GAMER_Debug,         1, NonFatal );
   LoadField( "Timing",             &RS.Timing,             SID, TID, NonFatal, &RT.Timing,              1, NonFatal );
   LoadField( "TimingSolver",       &RS.TimingSolver,       SID, TID, NonFatal, &RT.TimingSolver,        1, NonFatal );
   LoadField( "Intel",              &RS.Intel,              SID, TID, NonFatal, &RT.Intel,               1, NonFatal );
   LoadField( "Float8",             &RS.Float8,             SID, TID, NonFatal, &RT.Float8,              1,    Fatal );
   LoadField( "Serial",             &RS.Serial,             SID, TID, NonFatal, &RT.Serial,              1, NonFatal );
   LoadField( "LoadBalance",        &RS.LoadBalance,        SID, TID, NonFatal, &RT.LoadBalance,         1, NonFatal );
   LoadField( "OverlapMPI",         &RS.OverlapMPI,         SID, TID, NonFatal, &RT.OverlapMPI,          1, NonFatal );
   LoadField( "OpenMP",             &RS.OpenMP,             SID, TID, NonFatal, &RT.OpenMP,              1, NonFatal );
   LoadField( "GPU_Arch",           &RS.GPU_Arch,           SID, TID, NonFatal, &RT.GPU_Arch,            1, NonFatal );
   LoadField( "Laohu",              &RS.Laohu,              SID, TID, NonFatal, &RT.Laohu,               1, NonFatal );
   LoadField( "SupportHDF5",        &RS.SupportHDF5,        SID, TID, NonFatal, &RT.SupportHDF5,         1, NonFatal );

   LoadField( "NLevel",             &RS.NLevel,             SID, TID, NonFatal, &RT.NLevel,              1, NonFatal );
   LoadField( "MaxPatch",           &RS.MaxPatch,           SID, TID, NonFatal, &RT.MaxPatch,            1, NonFatal );

#  ifdef GRAVITY
   LoadField( "PotScheme",          &RS.PotScheme,          SID, TID, NonFatal, &RT.PotScheme,           1, NonFatal );
   LoadField( "StorePotGhost",      &RS.StorePotGhost,      SID, TID, NonFatal, &RT.StorePotGhost,       1, NonFatal );
   LoadField( "UnsplitGravity",     &RS.UnsplitGravity,     SID, TID, NonFatal, &RT.UnsplitGravity,      1, NonFatal );
#  endif

#  if   ( MODEL == HYDRO )
   LoadField( "FluScheme",          &RS.FluScheme,          SID, TID, NonFatal, &RT.FluScheme,           1, NonFatal );
#  ifdef LR_SCHEME
   LoadField( "LRScheme",           &RS.LRScheme,           SID, TID, NonFatal, &RT.LRScheme,            1, NonFatal );
#  endif
#  ifdef RSOLVER
   LoadField( "RSolver",            &RS.RSolver,            SID, TID, NonFatal, &RT.RSolver,             1, NonFatal );
#  endif
   LoadField( "NPassive",           &RS.NPassive,           SID, TID, NonFatal, &RT.NPassive,            1,    Fatal );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   LoadField( "ConserveMass",       &RS.ConserveMass,       SID, TID, NonFatal, &RT.ConserveMass,        1, NonFatal );
   LoadField( "Laplacian4th",       &RS.Laplacian4th,       SID, TID, NonFatal, &RT.Laplacian4th,        1, NonFatal );
   LoadField( "SelfInteraction4",   &RS.SelfInteraction4,   SID, TID, NonFatal, &RT.SelfInteraction4,    1, NonFatal );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
   LoadField( "StoreParAcc",        &RS.StoreParAcc,        SID, TID, NonFatal, &RT.StoreParAcc,         1, NonFatal );
#  endif


// 5. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : Check_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_SymConst
// Description :  Load and compare the SymConst_t structure (runtime vs. restart file)
//
// Note        :  1. Unnecessary comparisons are avoided
//                2. Program may be terminated if any fatal difference is found
//
// Parameter   :  FileName : Restart file name
//-------------------------------------------------------------------------------------------------------
void Check_SymConst( const char *FileName )
{

   const bool    Fatal = true;
   const bool NonFatal = false;

   herr_t Status;


// 1. fill in the runtime parameters
   SymConst_t RT;    // RT = RunTime
   FillIn_SymConst( RT );


// 2. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 3. load the SymConst dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/SymConst", H5P_DEFAULT );  // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/SymConst" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 4. load and compare all target fields one-by-one
   SymConst_t RS;    // RS = ReStart


   LoadField( "NComp",                &RS.NComp,                SID, TID, NonFatal, &RT.NComp,                 1,    Fatal );
   LoadField( "PatchSize",            &RS.PatchSize,            SID, TID, NonFatal, &RT.PatchSize,             1,    Fatal );
   LoadField( "Flu_NIn",              &RS.Flu_NIn,              SID, TID, NonFatal, &RT.Flu_NIn,               1, NonFatal );
   LoadField( "Flu_NOut",             &RS.Flu_NOut,             SID, TID, NonFatal, &RT.Flu_NOut,              1, NonFatal );
   LoadField( "NFlux",                &RS.NFlux,                SID, TID, NonFatal, &RT.NFlux,                 1, NonFatal );
   LoadField( "Flu_GhostSize",        &RS.Flu_GhostSize,        SID, TID, NonFatal, &RT.Flu_GhostSize,         1, NonFatal );
   LoadField( "Flu_Nxt",              &RS.Flu_Nxt,              SID, TID, NonFatal, &RT.Flu_Nxt,               1, NonFatal );
   LoadField( "Debug_HDF5",           &RS.Debug_HDF5,           SID, TID, NonFatal, &RT.Debug_HDF5,            1, NonFatal );
   LoadField( "SibOffsetNonperiodic", &RS.SibOffsetNonperiodic, SID, TID, NonFatal, &RT.SibOffsetNonperiodic,  1, NonFatal );
#  ifdef LOAD_BALANCE
   LoadField( "SonOffsetLB",          &RS.SonOffsetLB,          SID, TID, NonFatal, &RT.SonOffsetLB,           1, NonFatal );
#  endif
   LoadField( "TinyValue",            &RS.TinyValue,            SID, TID, NonFatal, &RT.TinyValue,             1, NonFatal );

#  ifdef GRAVITY
   LoadField( "Gra_NIn",              &RS.Gra_NIn,              SID, TID, NonFatal, &RT.Gra_NIn,               1, NonFatal );
   LoadField( "Pot_GhostSize",        &RS.Pot_GhostSize,        SID, TID, NonFatal, &RT.Pot_GhostSize,         1, NonFatal );
   LoadField( "Gra_GhostSize",        &RS.Gra_GhostSize,        SID, TID, NonFatal, &RT.Gra_GhostSize,         1, NonFatal );
   LoadField( "Rho_GhostSize",        &RS.Rho_GhostSize,        SID, TID, NonFatal, &RT.Rho_GhostSize,         1, NonFatal );
   LoadField( "Pot_Nxt",              &RS.Pot_Nxt,              SID, TID, NonFatal, &RT.Pot_Nxt,               1, NonFatal );
   LoadField( "Gra_Nxt",              &RS.Gra_Nxt,              SID, TID, NonFatal, &RT.Gra_Nxt,               1, NonFatal );
   LoadField( "Rho_Nxt",              &RS.Rho_Nxt,              SID, TID, NonFatal, &RT.Rho_Nxt,               1, NonFatal );
#  ifdef UNSPLIT_GRAVITY
   LoadField( "USG_GhostSize",        &RS.USG_GhostSize,        SID, TID, NonFatal, &RT.USG_GhostSize,         1, NonFatal );
   LoadField( "USG_NxtF",             &RS.USG_NxtF,             SID, TID, NonFatal, &RT.USG_NxtF,              1, NonFatal );
   LoadField( "USG_NxtG",             &RS.USG_NxtG,             SID, TID, NonFatal, &RT.USG_NxtG,              1, NonFatal );
#  endif
   LoadField( "Gra_BlockSize_z",      &RS.Gra_BlockSize_z,      SID, TID, NonFatal, &RT.Gra_BlockSize_z,       1, NonFatal );
   LoadField( "ExtPotNAuxMax",        &RS.ExtPotNAuxMax,        SID, TID, NonFatal, &RT.ExtPotNAuxMax,         1, NonFatal );
   LoadField( "ExtAccNAuxMax",        &RS.ExtAccNAuxMax,        SID, TID, NonFatal, &RT.ExtAccNAuxMax,         1, NonFatal );
#  if   ( POT_SCHEME == SOR )
   LoadField( "Pot_BlockSize_z",      &RS.Pot_BlockSize_z,      SID, TID, NonFatal, &RT.Pot_BlockSize_z,       1, NonFatal );
   LoadField( "UsePSolver_10to14",    &RS.UsePSolver_10to14,    SID, TID, NonFatal, &RT.UsePSolver_10to14,     1, NonFatal );
#  elif ( POT_SCHEME == MG )
   LoadField( "Pot_BlockSize_x",      &RS.Pot_BlockSize_x,      SID, TID, NonFatal, &RT.Pot_BlockSize_x,       1, NonFatal );
#  endif
#  endif // #ifdef GRAVITY

#  ifdef PARTICLE
   LoadField( "NPar_Var",             &RS.NPar_Var,             SID, TID, NonFatal, &RT.NPar_Var,              1,    Fatal );
   LoadField( "NPar_Passive",         &RS.NPar_Passive,         SID, TID, NonFatal, &RT.NPar_Passive,          1,    Fatal );
   LoadField( "PM_GhostSize",         &RS.PM_GhostSize,         SID, TID, NonFatal, &RT.PM_GhostSize,          1, NonFatal );
   LoadField( "PM_Nxt",               &RS.PM_Nxt,               SID, TID, NonFatal, &RT.PM_Nxt,                1, NonFatal );
   LoadField( "Debug_Particle",       &RS.Debug_Particle,       SID, TID, NonFatal, &RT.Debug_Particle,        1, NonFatal );
   LoadField( "ParList_GrowthFactor", &RS.ParList_GrowthFactor, SID, TID, NonFatal, &RT.ParList_GrowthFactor,  1, NonFatal );
   LoadField( "ParList_ReduceFactor", &RS.ParList_ReduceFactor, SID, TID, NonFatal, &RT.ParList_ReduceFactor,  1, NonFatal );
#  endif

#  if   ( MODEL == HYDRO )
   LoadField( "Flu_BlockSize_x",      &RS.Flu_BlockSize_x,      SID, TID, NonFatal, &RT.Flu_BlockSize_x,       1, NonFatal );
   LoadField( "Flu_BlockSize_y",      &RS.Flu_BlockSize_y,      SID, TID, NonFatal, &RT.Flu_BlockSize_y,       1, NonFatal );
   LoadField( "CheckNegativeInFluid", &RS.CheckNegativeInFluid, SID, TID, NonFatal, &RT.CheckNegativeInFluid,  1, NonFatal );
   LoadField( "CharReconstruction",   &RS.CharReconstruction,   SID, TID, NonFatal, &RT.CharReconstruction,    1, NonFatal );
   LoadField( "CheckIntermediate",    &RS.CheckIntermediate,    SID, TID, NonFatal, &RT.CheckIntermediate,     1, NonFatal );
   LoadField( "HLL_NoRefState",       &RS.HLL_NoRefState,       SID, TID, NonFatal, &RT.HLL_NoRefState,        1, NonFatal );
   LoadField( "HLL_IncludeAllWaves",  &RS.HLL_IncludeAllWaves,  SID, TID, NonFatal, &RT.HLL_IncludeAllWaves,   1, NonFatal );
   LoadField( "WAF_Dissipate",        &RS.WAF_Dissipate,        SID, TID, NonFatal, &RT.WAF_Dissipate,         1, NonFatal );
   LoadField( "PositiveDensInFixUp",  &RS.PositiveDensInFixUp,  SID, TID, NonFatal, &RT.PositiveDensInFixUp,   1, NonFatal );
#  ifdef N_FC_VAR
   LoadField( "N_FC_Var",             &RS.N_FC_Var,             SID, TID, NonFatal, &RT.N_FC_Var,              1, NonFatal );
#  endif
#  ifdef N_SLOPE_PPM
   LoadField( "N_Slope_PPM",          &RS.N_Slope_PPM,          SID, TID, NonFatal, &RT.N_Slope_PPM,           1, NonFatal );
#  endif
#  ifdef MIN_PRES_DENS
   LoadField( "Min_Pres_Dens",        &RS.Min_Pres_Dens,        SID, TID, NonFatal, &RT.Min_Pres_Dens,         1, NonFatal );
#  endif
#  ifdef MIN_PRES
   LoadField( "Min_Pres",             &RS.Min_Pres,             SID, TID, NonFatal, &RT.Min_Pres,              1, NonFatal );
#  endif
#  ifdef MAX_ERROR
   LoadField( "MaxError",             &RS.MaxError,             SID, TID, NonFatal, &RT.MaxError,              1, NonFatal );
#  endif

#  elif ( MODEL == MHD )
   LoadField( "Flu_BlockSize_x",      &RS.Flu_BlockSize_x,      SID, TID, NonFatal, &RT.Flu_BlockSize_x,       1, NonFatal );
   LoadField( "Flu_BlockSize_y",      &RS.Flu_BlockSize_y,      SID, TID, NonFatal, &RT.Flu_BlockSize_y,       1, NonFatal );
#  warning : WAIT MHD !!!

#  elif  ( MODEL == ELBDM )
   LoadField( "Flu_BlockSize_x",      &RS.Flu_BlockSize_x,      SID, TID, NonFatal, &RT.Flu_BlockSize_x,       1, NonFatal );
   LoadField( "Flu_BlockSize_y",      &RS.Flu_BlockSize_y,      SID, TID, NonFatal, &RT.Flu_BlockSize_y,       1, NonFatal );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 5. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : Check_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_InputPara
// Description :  Load and compare the InputPara_t structure (runtime vs. restart file)
//
// Note        :  1. Unnecessary comparisons are avoided
//                2. Program may be terminated if any fatal difference is found
//
// Parameter   :  FileName : Restart file name
//-------------------------------------------------------------------------------------------------------
void Check_InputPara( const char *FileName )
{

   const bool    Fatal = true;
   const bool NonFatal = false;
   const int  N1       = MAX_LEVEL;
   const int *NullPtr  = NULL;

   herr_t Status;


// 1. fill in the runtime parameters
   InputPara_t RT;   // RT = RunTime
   FillIn_InputPara( RT );


// 2. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 3. load the InputPara dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/InputPara", H5P_DEFAULT ); // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/InputPara" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 4. load and compare all target fields one-by-one
   InputPara_t RS;   // RS = ReStart

// simulation scale
   LoadField( "BoxSize",                 &RS.BoxSize,                 SID, TID, NonFatal, &RT.BoxSize,                  1,    Fatal );
   LoadField( "NX0_Tot",                  RS.NX0_Tot,                 SID, TID, NonFatal,  RT.NX0_Tot,                  3,    Fatal );
   LoadField( "MPI_NRank",               &RS.MPI_NRank,               SID, TID, NonFatal, &RT.MPI_NRank,                1, NonFatal );
   LoadField( "MPI_NRank_X",              RS.MPI_NRank_X,             SID, TID, NonFatal,  RT.MPI_NRank_X,              3, NonFatal );
   LoadField( "OMP_NThread",             &RS.OMP_NThread,             SID, TID, NonFatal, &RT.OMP_NThread,              1, NonFatal );
   LoadField( "EndT",                    &RS.EndT,                    SID, TID, NonFatal, &RT.EndT,                     1, NonFatal );
   LoadField( "EndStep",                 &RS.EndStep,                 SID, TID, NonFatal, &RT.EndStep,                  1, NonFatal );

// boundary condition
   LoadField( "Opt__BC_Flu",              RS.Opt__BC_Flu,             SID, TID, NonFatal,  RT.Opt__BC_Flu,              6, NonFatal );
#  ifdef GRAVITY
   LoadField( "Opt__BC_Pot",             &RS.Opt__BC_Pot,             SID, TID, NonFatal, &RT.Opt__BC_Pot,              1, NonFatal );
   LoadField( "GFunc_Coeff0",            &RS.GFunc_Coeff0,            SID, TID, NonFatal, &RT.GFunc_Coeff0,             1, NonFatal );
#  endif

// particle
#  ifdef PARTICLE
   LoadField( "Par_NPar",                &RS.Par_NPar,                SID, TID, NonFatal, &RT.Par_NPar,                 1, NonFatal );
   LoadField( "Par_Init",                &RS.Par_Init,                SID, TID, NonFatal, &RT.Par_Init,                 1, NonFatal );
   LoadField( "Par_Interp",              &RS.Par_Interp,              SID, TID, NonFatal, &RT.Par_Interp,               1, NonFatal );
   LoadField( "Par_Integ",               &RS.Par_Integ,               SID, TID, NonFatal, &RT.Par_Integ,                1, NonFatal );
   LoadField( "Par_SyncDump",            &RS.Par_SyncDump,            SID, TID, NonFatal, &RT.Par_SyncDump,             1, NonFatal );
   LoadField( "Par_ImproveAcc",          &RS.Par_ImproveAcc,          SID, TID, NonFatal, &RT.Par_ImproveAcc,           1, NonFatal );
   LoadField( "Par_PredictPos",          &RS.Par_PredictPos,          SID, TID, NonFatal, &RT.Par_PredictPos,           1, NonFatal );
   LoadField( "Par_RemoveCell",          &RS.Par_RemoveCell,          SID, TID, NonFatal, &RT.Par_RemoveCell,           1, NonFatal );
#  endif

// cosmology
#  ifdef COMOVING
   LoadField( "A_Init",                  &RS.A_Init,                  SID, TID, NonFatal, &RT.A_Init,                   1, NonFatal );
   LoadField( "OmegaM0",                 &RS.OmegaM0,                 SID, TID, NonFatal, &RT.OmegaM0,                  1, NonFatal );
#  endif

// time-step determination
   LoadField( "Dt__Fluid",               &RS.Dt__Fluid,               SID, TID, NonFatal, &RT.Dt__Fluid,                1, NonFatal );
   LoadField( "Dt__FluidInit",           &RS.Dt__FluidInit,           SID, TID, NonFatal, &RT.Dt__FluidInit,            1, NonFatal );
#  ifdef GRAVITY
   LoadField( "Dt__Gravity",             &RS.Dt__Gravity,             SID, TID, NonFatal, &RT.Dt__Gravity,              1, NonFatal );
#  endif
#  if ( MODEL == ELBDM )
   LoadField( "Dt__Phase",               &RS.Dt__Phase,               SID, TID, NonFatal, &RT.Dt__Phase,                1, NonFatal );
#  endif
#  ifdef PARTICLE 
   LoadField( "Dt__ParVel",              &RS.Dt__ParVel,              SID, TID, NonFatal, &RT.Dt__ParVel,               1, NonFatal );
   LoadField( "Dt__ParVelMax",           &RS.Dt__ParVelMax,           SID, TID, NonFatal, &RT.Dt__ParVelMax,            1, NonFatal );
   LoadField( "Dt__ParAcc",              &RS.Dt__ParAcc,              SID, TID, NonFatal, &RT.Dt__ParAcc,               1, NonFatal );
#  endif
#  ifdef COMOVING
   LoadField( "Dt__MaxDeltaA",           &RS.Dt__MaxDeltaA,           SID, TID, NonFatal, &RT.Dt__MaxDeltaA,            1, NonFatal );
#  endif
   LoadField( "Opt__AdaptiveDt",         &RS.Opt__AdaptiveDt,         SID, TID, NonFatal, &RT.Opt__AdaptiveDt,          1, NonFatal );
   LoadField( "Opt__RecordDt",           &RS.Opt__RecordDt,           SID, TID, NonFatal, &RT.Opt__RecordDt,            1, NonFatal );
   LoadField( "Opt__DtUser",             &RS.Opt__DtUser,             SID, TID, NonFatal, &RT.Opt__DtUser,              1, NonFatal );
   

// domain refinement
   LoadField( "RegridCount",             &RS.RegridCount,             SID, TID, NonFatal, &RT.RegridCount,              1, NonFatal );
   LoadField( "FlagBufferSize",          &RS.FlagBufferSize,          SID, TID, NonFatal, &RT.FlagBufferSize,           1, NonFatal );
   LoadField( "MaxLevel",                &RS.MaxLevel,                SID, TID, NonFatal, &RT.MaxLevel,                 1, NonFatal );
   LoadField( "Opt__Flag_Rho",           &RS.Opt__Flag_Rho,           SID, TID, NonFatal, &RT.Opt__Flag_Rho,            1, NonFatal );
   LoadField( "Opt__Flag_RhoGradient",   &RS.Opt__Flag_RhoGradient,   SID, TID, NonFatal, &RT.Opt__Flag_RhoGradient,    1, NonFatal );
#  if ( MODEL == HYDRO ) 
   LoadField( "Opt__Flag_PresGradient",  &RS.Opt__Flag_PresGradient,  SID, TID, NonFatal, &RT.Opt__Flag_PresGradient,   1, NonFatal );
#  endif
#  if ( MODEL == ELBDM ) 
   LoadField( "Opt__Flag_EngyDensity",   &RS.Opt__Flag_EngyDensity,   SID, TID, NonFatal, &RT.Opt__Flag_EngyDensity,    1, NonFatal );
#  endif
   LoadField( "Opt__Flag_LohnerDens",    &RS.Opt__Flag_LohnerDens,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerDens,     1, NonFatal );
#  if ( MODEL == HYDRO ) 
   LoadField( "Opt__Flag_LohnerEngy",    &RS.Opt__Flag_LohnerEngy,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerEngy,     1, NonFatal );
   LoadField( "Opt__Flag_LohnerPres",    &RS.Opt__Flag_LohnerPres,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerPres,     1, NonFatal );
#  endif
   LoadField( "Opt__Flag_LohnerForm",    &RS.Opt__Flag_LohnerForm,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerForm,     1, NonFatal );
   LoadField( "Opt__Flag_User",          &RS.Opt__Flag_User,          SID, TID, NonFatal, &RT.Opt__Flag_User,           1, NonFatal );
   LoadField( "Opt__Flag_Region",        &RS.Opt__Flag_Region,        SID, TID, NonFatal, &RT.Opt__Flag_Region,         1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__Flag_NParPatch",     &RS.Opt__Flag_NParPatch,     SID, TID, NonFatal, &RT.Opt__Flag_NParPatch,      1, NonFatal );
#  endif
   LoadField( "Opt__PatchCount",         &RS.Opt__PatchCount,         SID, TID, NonFatal, &RT.Opt__PatchCount,          1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__ParLevel",           &RS.Opt__ParLevel,           SID, TID, NonFatal, &RT.Opt__ParLevel,            1, NonFatal );
#  endif

// load balance
#  ifdef LOAD_BALANCE
   LoadField( "LB_Input__WLI_Max",       &RS.LB_Input__WLI_Max,       SID, TID, NonFatal, &RT.LB_Input__WLI_Max,        1, NonFatal );
#  endif

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   LoadField( "Gamma",                   &RS.Gamma,                   SID, TID, NonFatal, &RT.Gamma,                    1, NonFatal );
   LoadField( "MinMod_Coeff",            &RS.MinMod_Coeff,            SID, TID, NonFatal, &RT.MinMod_Coeff,             1, NonFatal );
   LoadField( "EP_Coeff",                &RS.EP_Coeff,                SID, TID, NonFatal, &RT.EP_Coeff,                 1, NonFatal );
   LoadField( "Opt__LR_Limiter",         &RS.Opt__LR_Limiter,         SID, TID, NonFatal, &RT.Opt__LR_Limiter,          1, NonFatal );
   LoadField( "Opt__WAF_Limiter",        &RS.Opt__WAF_Limiter,        SID, TID, NonFatal, &RT.Opt__WAF_Limiter,         1, NonFatal );
   LoadField( "Opt__CorrUnphyScheme",    &RS.Opt__CorrUnphyScheme,    SID, TID, NonFatal, &RT.Opt__CorrUnphyScheme,     1, NonFatal );
#  endif

// ELBDM solvers
#  if ( MODEL == ELBDM )
   LoadField( "ELBDM_Mass",              &RS.ELBDM_Mass,              SID, TID, NonFatal, &RT.ELBDM_Mass,               1, NonFatal );
   LoadField( "ELBDM_PlanckConst",       &RS.ELBDM_PlanckConst,       SID, TID, NonFatal, &RT.ELBDM_PlanckConst,        1, NonFatal );
#  ifdef QUARTIC_SELF_INTERACTION
   LoadField( "ELBDM_Lambda",            &RS.ELBDM_Lambda,            SID, TID, NonFatal, &RT.ELBDM_Lambda,             1, NonFatal );
#  endif
   LoadField( "ELBDM_Taylor3_Coeff",     &RS.ELBDM_Taylor3_Coeff,     SID, TID, NonFatal, &RT.ELBDM_Taylor3_Coeff,      1, NonFatal );
   LoadField( "ELBDM_Taylor3_Auto",      &RS.ELBDM_Taylor3_Auto,      SID, TID, NonFatal, &RT.ELBDM_Taylor3_Auto,       1, NonFatal );
#  endif

// fluid solvers in both HYDRO/MHD/ELBDM
   LoadField( "Flu_GPU_NPGroup",         &RS.Flu_GPU_NPGroup,         SID, TID, NonFatal, &RT.Flu_GPU_NPGroup,          1, NonFatal );
   LoadField( "GPU_NStream",             &RS.GPU_NStream,             SID, TID, NonFatal, &RT.GPU_NStream,              1, NonFatal );
   LoadField( "Opt__FixUp_Flux",         &RS.Opt__FixUp_Flux,         SID, TID, NonFatal, &RT.Opt__FixUp_Flux,          1, NonFatal );
   LoadField( "Opt__FixUp_Restrict",     &RS.Opt__FixUp_Restrict,     SID, TID, NonFatal, &RT.Opt__FixUp_Restrict,      1, NonFatal );
   LoadField( "Opt__OverlapMPI",         &RS.Opt__OverlapMPI,         SID, TID, NonFatal, &RT.Opt__OverlapMPI,          1, NonFatal );
   LoadField( "Opt__ResetFluid",         &RS.Opt__ResetFluid,         SID, TID, NonFatal, &RT.Opt__ResetFluid,          1, NonFatal );
   LoadField( "Opt__CorrUnphy",          &RS.Opt__CorrUnphy,          SID, TID, NonFatal, &RT.Opt__CorrUnphy,           1, NonFatal );

// self-gravity
#  ifdef GRAVITY
   LoadField( "NewtonG",                 &RS.NewtonG,                 SID, TID, NonFatal, &RT.NewtonG,                  1, NonFatal );
#  if   ( POT_SCHEME == SOR )
   LoadField( "SOR_Omega",               &RS.SOR_Omega,               SID, TID, NonFatal, &RT.SOR_Omega,                1, NonFatal );
   LoadField( "SOR_MaxIter",             &RS.SOR_MaxIter,             SID, TID, NonFatal, &RT.SOR_MaxIter,              1, NonFatal );
   LoadField( "SOR_MinIter",             &RS.SOR_MinIter,             SID, TID, NonFatal, &RT.SOR_MinIter,              1, NonFatal );
#  elif ( POT_SCHEME == MG )
   LoadField( "MG_MaxIter",              &RS.MG_MaxIter,              SID, TID, NonFatal, &RT.MG_MaxIter,               1, NonFatal );
   LoadField( "MG_NPreSmooth",           &RS.MG_NPreSmooth,           SID, TID, NonFatal, &RT.MG_NPreSmooth,            1, NonFatal );
   LoadField( "MG_NPostSmooth",          &RS.MG_NPostSmooth,          SID, TID, NonFatal, &RT.MG_NPostSmooth,           1, NonFatal );
   LoadField( "MG_ToleratedError",       &RS.MG_ToleratedError,       SID, TID, NonFatal, &RT.MG_ToleratedError,        1, NonFatal );
#  endif
   LoadField( "Pot_GPU_NPGroup",         &RS.Pot_GPU_NPGroup,         SID, TID, NonFatal, &RT.Pot_GPU_NPGroup,          1, NonFatal );
   LoadField( "Opt__GraP5Gradient",      &RS.Opt__GraP5Gradient,      SID, TID, NonFatal, &RT.Opt__GraP5Gradient,       1, NonFatal );
   LoadField( "Opt__GravityType",        &RS.Opt__GravityType,        SID, TID, NonFatal, &RT.Opt__GravityType,         1, NonFatal );
   LoadField( "Opt__ExternalPot",        &RS.Opt__ExternalPot,        SID, TID, NonFatal, &RT.Opt__ExternalPot,         1, NonFatal );
#  endif

// initialization
   LoadField( "Opt__Init",               &RS.Opt__Init,               SID, TID, NonFatal, &RT.Opt__Init,                1, NonFatal );
   LoadField( "Opt__RestartHeader",      &RS.Opt__RestartHeader,      SID, TID, NonFatal, &RT.Opt__RestartHeader,       1, NonFatal );
   LoadField( "Opt__UM_Start_Level",     &RS.Opt__UM_Start_Level,     SID, TID, NonFatal, &RT.Opt__UM_Start_Level,      1, NonFatal );
   LoadField( "Opt__UM_Start_NVar",      &RS.Opt__UM_Start_NVar,      SID, TID, NonFatal, &RT.Opt__UM_Start_NVar,       1, NonFatal );
   LoadField( "Opt__UM_Start_Downgrade", &RS.Opt__UM_Start_Downgrade, SID, TID, NonFatal, &RT.Opt__UM_Start_Downgrade,  1, NonFatal );
   LoadField( "Opt__UM_Start_Refine",    &RS.Opt__UM_Start_Refine,    SID, TID, NonFatal, &RT.Opt__UM_Start_Refine,     1, NonFatal );
   LoadField( "Opt__UM_Factor_5over3",   &RS.Opt__UM_Factor_5over3,   SID, TID, NonFatal, &RT.Opt__UM_Factor_5over3,    1, NonFatal );
   LoadField( "Opt__InitRestrict",       &RS.Opt__InitRestrict,       SID, TID, NonFatal, &RT.Opt__InitRestrict,        1, NonFatal );
   LoadField( "Opt__GPUID_Select",       &RS.Opt__GPUID_Select,       SID, TID, NonFatal, &RT.Opt__GPUID_Select,        1, NonFatal );
   LoadField( "Init_Subsampling_NCell",  &RS.Init_Subsampling_NCell,  SID, TID, NonFatal, &RT.Init_Subsampling_NCell,   1, NonFatal );

// interpolation schemes
   LoadField( "Opt__Int_Time",           &RS.Opt__Int_Time,           SID, TID, NonFatal, &RT.Opt__Int_Time,            1, NonFatal );
#  if ( MODEL == ELBDM ) 
   LoadField( "Opt__Int_Phase",          &RS.Opt__Int_Phase,          SID, TID, NonFatal, &RT.Opt__Int_Phase,           1, NonFatal );
#  endif
   LoadField( "Opt__Flu_IntScheme",      &RS.Opt__Flu_IntScheme,      SID, TID, NonFatal, &RT.Opt__Flu_IntScheme,       1, NonFatal );
#  ifdef GRAVITY
   LoadField( "Opt__Pot_IntScheme",      &RS.Opt__Pot_IntScheme,      SID, TID, NonFatal, &RT.Opt__Pot_IntScheme,       1, NonFatal );
   LoadField( "Opt__Rho_IntScheme",      &RS.Opt__Rho_IntScheme,      SID, TID, NonFatal, &RT.Opt__Rho_IntScheme,       1, NonFatal );
   LoadField( "Opt__Gra_IntScheme",      &RS.Opt__Gra_IntScheme,      SID, TID, NonFatal, &RT.Opt__Gra_IntScheme,       1, NonFatal );
#  endif
   LoadField( "Opt__RefFlu_IntScheme",   &RS.Opt__RefFlu_IntScheme,   SID, TID, NonFatal, &RT.Opt__RefFlu_IntScheme,    1, NonFatal );
#  ifdef GRAVITY
   LoadField( "Opt__RefPot_IntScheme",   &RS.Opt__RefPot_IntScheme,   SID, TID, NonFatal, &RT.Opt__RefPot_IntScheme,    1, NonFatal );
#  endif
   LoadField( "IntMonoCoeff",            &RS.IntMonoCoeff,            SID, TID, NonFatal, &RT.IntMonoCoeff,             1, NonFatal );

// data dump
   LoadField( "Opt__Output_Total",       &RS.Opt__Output_Total,       SID, TID, NonFatal, &RT.Opt__Output_Total,        1, NonFatal );
   LoadField( "Opt__Output_Part",        &RS.Opt__Output_Part,        SID, TID, NonFatal, &RT.Opt__Output_Part,         1, NonFatal );
   LoadField( "Opt__Output_TestError",   &RS.Opt__Output_TestError,   SID, TID, NonFatal, &RT.Opt__Output_TestError,    1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__Output_Particle",    &RS.Opt__Output_Particle,    SID, TID, NonFatal, &RT.Opt__Output_Particle,     1, NonFatal );
#  endif
   LoadField( "Opt__Output_BasePS",      &RS.Opt__Output_BasePS,      SID, TID, NonFatal, &RT.Opt__Output_BasePS,       1, NonFatal );
   if ( OPT__OUTPUT_PART )
   LoadField( "Opt__Output_Base",        &RS.Opt__Output_Base,        SID, TID, NonFatal, &RT.Opt__Output_Base,         1, NonFatal );
#  ifdef GRAVITY
   LoadField( "Opt__Output_Pot",         &RS.Opt__Output_Pot,         SID, TID, NonFatal, &RT.Opt__Output_Pot,          1, NonFatal );
#  endif
#  ifdef PARTICLE
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_TEST_ERROR || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PARTICLE ) {
#  else
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_TEST_ERROR || OPT__OUTPUT_BASEPS ) {
#  endif
   LoadField( "Opt__Output_Mode",        &RS.Opt__Output_Mode,        SID, TID, NonFatal, &RT.Opt__Output_Mode,         1, NonFatal );
   LoadField( "Opt__Output_Step",        &RS.Opt__Output_Step,        SID, TID, NonFatal, &RT.Opt__Output_Step,         1, NonFatal );
   LoadField( "Opt__Output_Dt",          &RS.Opt__Output_Dt,          SID, TID, NonFatal, &RT.Opt__Output_Dt,           1, NonFatal );
   }
   if ( OPT__OUTPUT_PART ) {
   LoadField( "Output_PartX",            &RS.Output_PartX,            SID, TID, NonFatal, &RT.Output_PartX,             1, NonFatal );
   LoadField( "Output_PartY",            &RS.Output_PartY,            SID, TID, NonFatal, &RT.Output_PartY,             1, NonFatal );
   LoadField( "Output_PartZ",            &RS.Output_PartZ,            SID, TID, NonFatal, &RT.Output_PartZ,             1, NonFatal );
   }
   LoadField( "InitDumpID",              &RS.InitDumpID,              SID, TID, NonFatal, &RT.InitDumpID,               1, NonFatal );

// miscellaneous
   LoadField( "Opt__Verbose",            &RS.Opt__Verbose,            SID, TID, NonFatal, &RT.Opt__Verbose,             1, NonFatal );
   LoadField( "Opt__TimingBalance",      &RS.Opt__TimingBalance,      SID, TID, NonFatal, &RT.Opt__TimingBalance,       1, NonFatal );
   LoadField( "Opt__TimingMPI",          &RS.Opt__TimingMPI,          SID, TID, NonFatal, &RT.Opt__TimingMPI,           1, NonFatal );
   LoadField( "Opt__RecordMemory",       &RS.Opt__RecordMemory,       SID, TID, NonFatal, &RT.Opt__RecordMemory,        1, NonFatal );
   LoadField( "Opt__RecordPerformance",  &RS.Opt__RecordPerformance,  SID, TID, NonFatal, &RT.Opt__RecordPerformance,   1, NonFatal );
   LoadField( "Opt__ManualControl",      &RS.Opt__ManualControl,      SID, TID, NonFatal, &RT.Opt__ManualControl,       1, NonFatal );
   LoadField( "Opt__RecordUser",         &RS.Opt__RecordUser,         SID, TID, NonFatal, &RT.Opt__RecordUser,          1, NonFatal );

// simulation checks
   LoadField( "Opt__Ck_Refine",          &RS.Opt__Ck_Refine,          SID, TID, NonFatal, &RT.Opt__Ck_Refine,           1, NonFatal );
   LoadField( "Opt__Ck_ProperNesting",   &RS.Opt__Ck_ProperNesting,   SID, TID, NonFatal, &RT.Opt__Ck_ProperNesting,    1, NonFatal );
   LoadField( "Opt__Ck_Conservation",    &RS.Opt__Ck_Conservation,    SID, TID, NonFatal, &RT.Opt__Ck_Conservation,     1, NonFatal );
   LoadField( "Opt__Ck_Restrict",        &RS.Opt__Ck_Restrict,        SID, TID, NonFatal, &RT.Opt__Ck_Restrict,         1, NonFatal );
   LoadField( "Opt__Ck_Finite",          &RS.Opt__Ck_Finite,          SID, TID, NonFatal, &RT.Opt__Ck_Finite,           1, NonFatal );
   LoadField( "Opt__Ck_PatchAllocate",   &RS.Opt__Ck_PatchAllocate,   SID, TID, NonFatal, &RT.Opt__Ck_PatchAllocate,    1, NonFatal );
   LoadField( "Opt__Ck_FluxAllocate",    &RS.Opt__Ck_FluxAllocate,    SID, TID, NonFatal, &RT.Opt__Ck_FluxAllocate,     1, NonFatal );
#  if ( MODEL == HYDRO )
   LoadField( "Opt__Ck_Negative",        &RS.Opt__Ck_Negative,        SID, TID, NonFatal, &RT.Opt__Ck_Negative,         1, NonFatal );
#  endif
   LoadField( "Opt__Ck_MemFree",         &RS.Opt__Ck_MemFree,         SID, TID, NonFatal, &RT.Opt__Ck_MemFree,          1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__Ck_Particle",        &RS.Opt__Ck_Particle,        SID, TID, NonFatal, &RT.Opt__Ck_Particle,         1, NonFatal );
#  endif


// flag tables
#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   const bool Opt__FlagLohner = ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES );
#  elif ( MODEL == ELBDM )
   const bool Opt__FlagLohner = OPT__FLAG_LOHNER_DENS;
#  endif

// initialize as -1 (to work with NLvRestart < NLEVEL)
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      RS.FlagTable_Rho         [lv]    = -1.0; 
      RS.FlagTable_RhoGradient [lv]    = -1.0; 

      for (int t=0; t<3; t++)
      RS.FlagTable_Lohner      [lv][t] = -1.0; 

      RS.FlagTable_User        [lv]    = -1.0; 

#     if   ( MODEL == HYDRO )
      RS.FlagTable_PresGradient[lv]    = -1.0; 

#     elif ( MODEL == ELBDM )
      for (int t=0; t<2; t++)
      RS.FlagTable_EngyDensity [lv][t] = -1.0; 
#     endif

#     ifdef PARTICLE
      RS.FlagTable_NParPatch   [lv]    = -1; 
#     endif
   }

   if ( OPT__FLAG_RHO )
   LoadField( "FlagTable_Rho",            RS.FlagTable_Rho,           SID, TID, NonFatal,  RT.FlagTable_Rho,           N1, NonFatal );
   
   if ( OPT__FLAG_RHO_GRADIENT )
   LoadField( "FlagTable_RhoGradient",    RS.FlagTable_RhoGradient,   SID, TID, NonFatal,  RT.FlagTable_RhoGradient,   N1, NonFatal );

   if ( Opt__FlagLohner ) {
   LoadField( "FlagTable_Lohner",         RS.FlagTable_Lohner,        SID, TID, NonFatal,  NullPtr,                    -1, NonFatal );

   for (int lv=0; lv<MAX_LEVEL; lv++)
   for (int t=0; t<3; t++)
   {
      if ( RS.FlagTable_Lohner[lv][t] != RT.FlagTable_Lohner[lv][t] )
         Aux_Message( stderr, "WARNING : \"%s[%d][%d]\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                       "FlagTable_Lohner", lv, t, RS.FlagTable_Lohner[lv][t],  RT.FlagTable_Lohner[lv][t] );
   }}

   if ( OPT__FLAG_USER )
   LoadField( "FlagTable_User",           RS.FlagTable_User,          SID, TID, NonFatal,  RT.FlagTable_User,          N1, NonFatal );

#  if   ( MODEL == HYDRO )
   if ( OPT__FLAG_PRES_GRADIENT )
   LoadField( "FlagTable_PresGradient",   RS.FlagTable_PresGradient,  SID, TID, NonFatal,  RT.FlagTable_PresGradient,  N1, NonFatal );

#  elif ( MODEL == ELBDM )
   if ( OPT__FLAG_ENGY_DENSITY ) {
   LoadField( "FlagTable_EngyDensity",    RS.FlagTable_EngyDensity,   SID, TID, NonFatal,  NullPtr,                    -1, NonFatal );

   for (int lv=0; lv<MAX_LEVEL; lv++)
   for (int t=0; t<2; t++)
   {
      if ( RS.FlagTable_EngyDensity[lv][t] != RT.FlagTable_EngyDensity[lv][t] )
         Aux_Message( stderr, "WARNING : \"%s[%d][%d]\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                       "FlagTable_EngyDensity", lv, t, RS.FlagTable_EngyDensity[lv][t],  RT.FlagTable_EngyDensity[lv][t] );
   }}
#  endif

#  ifdef PARTICLE
   if ( OPT__FLAG_NPAR_PATCH )
   LoadField( "FlagTable_NParPatch",      RS.FlagTable_NParPatch,     SID, TID, NonFatal,  RT.FlagTable_NParPatch,     N1, NonFatal );
#  endif


// 5. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : Check_InputPara



//-------------------------------------------------------------------------------------------------------
// Function    :  ResetParameter
// Description :  Reset simulation parameters from the restart file
//
// Note        :  1. Replace some runtime parameters by the values loaded from the restart file
//                2. This function must be called by ALL ranks
//
// Parameter   :  FileName : Restart file name
//                EndT     : Runtime parameter "END_T" 
//                EndStep  : Runtime parameter "END_STEP"
//-------------------------------------------------------------------------------------------------------
void ResetParameter( const char *FileName, double *EndT, long *EndStep )
{

   const bool    Fatal = true;
   const bool NonFatal = false;
   const int *NullPtr  = NULL;

   herr_t Status;


// 1. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 2. load the InputPar dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/InputPara", H5P_DEFAULT ); // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/InputPara" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 3. reset the runtime parameters if required
   InputPara_t RS;   // RS = ReStart

   if ( *EndT < 0.0 )
   {
      LoadField( "EndT",      EndT,    SID, TID, Fatal, NullPtr, -1, NonFatal );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "      NOTE : parameter %s is reset to %14.7e\n", "END_T", *EndT );
   }

   if ( *EndStep < 0 )
   {
      LoadField( "EndStep",   EndStep, SID, TID, Fatal, NullPtr, -1, NonFatal );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "      NOTE : parameter %s is reset to %ld\n", "END_STEP", *EndStep );
   }

// 4. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : ResetParameter



#endif // #ifdef SUPPORT_HDF5
