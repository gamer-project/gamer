#ifdef SUPPORT_HDF5

#include "GAMER.h"
#include "HDF5_Typedef.h"
#include <ctime>

void FillIn_KeyInfo  (   KeyInfo_t &KeyInfo, const int NFieldStored );
void FillIn_Makefile (  Makefile_t &Makefile  );
void FillIn_SymConst (  SymConst_t &SymConst  );
void FillIn_InputPara( InputPara_t &InputPara, const int NFieldStored, char FieldLabelOut[][MAX_STRING] );

static void GetCompound_KeyInfo  ( hid_t &H5_TypeID );
static void GetCompound_Makefile ( hid_t &H5_TypeID );
static void GetCompound_SymConst ( hid_t &H5_TypeID );
static void GetCompound_InputPara( hid_t &H5_TypeID, const int NFieldStored );
static void GetCompound_General  ( hid_t &H5_TypeID, const HDF5_Output_t *HDF5_Output );

void (*Output_HDF5_InputTest_Ptr)( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest ) = NULL;
void (*Output_HDF5_UserPara_Ptr)( HDF5_Output_t *HDF5_UserPara ) = NULL;
static herr_t H5_write_compound( const hid_t H5_SetID, const hid_t H5_Type_ID, const HDF5_Output_t *HDF5_Output );

static void Output_HDF5_UserPara_Template( HDF5_Output_t *HDF5_UserPara );



/*======================================================================================================
Data structure:
/ -> |
     | -> Info group     -> | -> InputPara dset (compound)
     |                      | -> InputTest dset (compound)
     |                      | -> KeyInfo   dset (compound)
     |                      | -> Makefile  dset (compound)
     |                      | -> SymConst  dset (compound)
     |
     | -> User group     -> | -> UserPara dset (compound)
     |
     | -> Tree group     -> | -> Corner  dset -> Cvt2Phy attrs
     |                      | -> LBIdx   dset
     |                      | -> Father  dset
     |                      | -> Son     dset
     |                      | -> Sibling dset
     |                      | -> NPar    dset
     |
     | -> GridData group -> | -> Dens dset
     |                      | -> ...
     |                      | -> ...
     |
     | -> Particle group -> | -> ParMass dset
                            | -> ...
                            | -> ...
======================================================================================================*/



/*======================================================================================================
h5py usage (with python 2):
1. Load file: f=h5py.File("Data_000000", "r")
2. Shows the names of all groups: list(f) or f.keys()
3. Show the names of all attributes: list(f['Tree']['Corner'].attrs) or f['Tree']['Corner'].attrs.keys()
4. Show a specific attribute: f['Tree']['Corner'].attrs['Cvt2Phy']
5. Show all variables in a compound variable: f['Info']['KeyInfo'].dtype
6. Show the value of a dataset: f['Tree']['Corner'].value
7. Show density of a patch with a global ID (GID) 1234: f['Data']['Dens'][1234]
8. Show density at a specific cell [1][2][3]: f['Data']['Dens'][1234][1][2][3]
======================================================================================================*/



/*======================================================================================================
Procedure for outputting new variables:
1. Add the new variable into one of the data structures (XXX) defined in "HDF5_Typedef.h"
2. Edit "GetCompound_XXX" to insert the new variables into the compound datatype
3. Edit "FillIn_XXX" to fill in the new variables
4. Edit "Check_XXX" in "Init_ByRestart_HDF5.cpp" to load and compare the new variables
5. Update FormatVersion
======================================================================================================*/




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Total_HDF5 (FormatVersion = 2503)
// Description :  Output all simulation data in the HDF5 format, which can be used as a restart file
//                or loaded by YT
//
// Note        :  1. Please refer to the "Data structure" described on the top of this file
//                2. Patch IDs stored in the HDF5 output are always GID (global identification) instead of
//                   PID (patch identification). Unlike PID, which always starts from 0 at different ranks
//                   and different levels, GID is unique among all patches at all ranks and all levels
//                   --> Each patch has an unique GID
//                3. Both "Father, Son, and Sibling[26]" are GID instead of PID
//                4. Currently we always use HDF5 NATVIE datatypes for both memory and dataset
//                5. All arrays in the "Tree" group (e.g., Corner, LBIdx, ...) and the "GridData" group (e.g., Dens, MomX, ...)
//                   have been sorted according to GID
//                   --> Moreover, currently we store all patches at the same level together
//                   --> A higher-level patch always has a larger GID than a lower-level patch
//                5. LBIdx dataset stores the LB_Idx of all patches sorted by their GIDs
//                   --> This list will be created even when LOAD_BALANCE is not turned on
//                       so that a serial output can be loaded by a parallel run easily
//                6. All C structures (e.g., "KeyInfo, SymConst, ..." are stored as a single (scalar)
//                   compound datetype
//                   --> Better update h5py to version >= 2.3.0 to properly read it in python
//                   --> Corresponding C structures are defined in "HDF5_Typedef.h"
//                7. It seems that h5py still have problem for "modifying" the loaded data. But reading data is fine.
//                8. The "H5T_GAMER_REAL" datatype will be mapped to "H5T_NATIVE_DOUBLE / H5T_NATIVE_FLOAT" if
//                   FLOAT8 is on / off
//                9. It is found that in the parallel environment each rank must try to "synchronize" the HDF5 file
//                   before opening the existed file and add data
//                   --> To achieve that, always invoke "SyncHDF5File" before calling "H5Fopen"
//                   --> "SyncHDF5File" is defined in "HDF5_Typedef.h", which simply opens the file
//                       with the appending mode and then closes it immediately
//                10. With PARTICLE on, two additional particle information will be recorded:
//                    --> "NPar" dataset under "Tree" records the number of active particles in all patches
//                         sorted by their GIDs
//                    --> "Particle" dataset under "/" stores all particle attributes (e.g., mass, velocity, ...)
//                        --> Currently we store different attributes in separate datasets
//                        --> Particles are stored in the order of their associated GIDs as well, but the order of
//                            particles in the same patch is not specified
//
// Parameter   :  FileName : Name of the output file
//
// Revision    :  2210 : 2016/10/03 --> output HUBBLE0, OPT__UNIT, UNIT_L/M/T/V/D/E, MOLECULAR_WEIGHT
//                2216 : 2016/11/27 --> output OPT__FLAG_LOHNER_TEMP
//                2217 : 2017/01/25 --> output RESTART_LOAD_NRANK, set CodeVersion to "gamer"
//                2218 : 2017/01/28 --> output OPT__FLAG_VORTICITY and the corresponding flag table
//                2220 : 2017/02/14 --> output passive grid and particle variables
//                2221 : 2017/02/20 --> output TINY_NUMBER and HUGE_NUMBER
//                2222 : 2017/02/20 --> output OPT__NORMALIZE_PASSIVE
//                2223 : 2017/02/22 --> output NormalizePassive_NVar and NormalizePassive_VarIdx
//                2224 : 2017/02/25 --> output OPT__CK_NORMALIZE_PASSIVE
//                2225 : 2017/03/01 --> output LB_Par_Weight, rename LB_Input__WLI_Max as LB_WLI_Max
//                2226 : 2017/03/03 --> output Opt__RecordLoadBalance
//                2227 : 2017/03/21 --> output PassiveFieldName_Grid and PassiveFieldName_Par
//                2228 : 2017/03/21 --> output NCOMP_FLUID, NCOMP_PASSIVE, and PAR_NPASSIVE in KeyInfo_t
//                2229 : 2017/04/06 --> output DUAL_ENERGY and DUAL_ENERGY_SWITCH
//                2230 : 2017/05/08 --> output OPT__FLAG_PAR_MASS_CELL and FlagTable_ParMassCell
//                2231 : 2017/05/08 --> output OPT__FLAG_JEANS and FlagTable_Jeans
//                2232 : 2017/06/13 --> output TESTPROB_ID
//                2233 : 2017/06/13 --> rename Opt__Output_TestError as Opt__Output_User
//                2234 : 2017/06/25 --> output Grackle variables
//                2235 : 2017/07/11 --> replace IndividualDt and Opt__AdaptiveDt by Opt__DtLevel
//                2236 : 2017/07/16 --> output DT_FLU_BLOCK_SIZE and DT_FLU_USE_SHUFFLE
//                2237 : 2017/07/17 --> output DT__FLEXIBLE_RANGE
//                2238 : 2017/07/19 --> output DT_GRA_BLOCK_SIZE_Z and DT_FLU_USE_SHUFFLE
//                2239 : 2017/07/25 --> output JEANS_MIN_PRES, JEANS_MIN_PRES_LEVEL, JEANS_MIN_PRES_NCELL
//                2240 : 2017/07/26 --> output AUTO_REDUCE_DT*
//                2241 : 2017/08/05 --> output FLAG_BUFFER_SIZE_MAXM1_LV
//                2242 : 2017/08/09 --> output DT__SYNC_PARENT_LV and DT__SYNC_CHILDREN_LV; remove DT__FLEXIBLE_RANGE
//                2250 : 2017/08/09 --> output dTime_AllLv[]
//                2251 : 2017/08/12 --> output OPT__RESTART_RESET
//                2252 : 2017/08/17 --> output GRACKLE_PE_HEATING and GRACKLE_PE_HEATING_RATE
//                2253 : 2017/08/27 --> output STAR_FORMATION
//                2254 : 2017/08/29 --> output star formation parameters
//                2255 : 2017/09/01 --> output SF_CREATION_STAR_DET_RANDOM
//                2256 : 2017/09/10 --> output FLAG_BUFFER_SIZE_MAXM2_LV
//                2257 : 2017/09/17 --> output OPT__OPTIMIZE_AGGRESSIVE
//                2258 : 2017/09/21 --> output OPT__MINIMIZE_MPI_BARRIER
//                2259 : 2017/10/10 --> output OPT__INIT_GRID_WITH_OMP
//                2261 : 2017/12/05 --> no longer define INTEL
//                2262 : 2017/12/27 --> rename all UM variables
//                2263 : 2017/12/27 --> remove OPT__RESTART_HEADER
//                2264 : 2018/02/28 --> add RANDOM_NUMBER
//                2265 : 2018/04/02 --> add OPT__NO_FLAG_NEAR_BOUNDARY
//                2266 : 2018/05/11 --> add OPT__UM_IC_LOAD_NRANK
//                2300 : 2018/07/15 --> replace PAR_NVAR and PAR_NPASSIVE by PAR_NATT_STORED and PAR_NATT_USER;
//                                      use the new infrastructure for adding user-defined grid fields and
//                                      particle attributes
//                                      --> imcompatible with version 2266 for the data with user-defined grid fields
//                                          and particle attributes as their labels may have changed
//                2301 : 2018/07/24 --> add OPT__UM_IC_FORMAT, PAR_IC_FORMAT, and PAR_IC_MASS
//                2302 : 2018/07/24 --> replace GRACKLE_MODE by GRACKLE_ACTIVATE
//                2303 : 2018/10/04 --> set "CodeVersion" to VERSION defined in Macro.h
//                2304 : 2018/12/10 --> remove EP_Coeff that no longer exists
//                2305 : 2018/12/15 --> remove variables related to the WAF scheme
//                2306 : 2018/12/25 --> replace DT_GRA_BLOCK_SIZE_Z by DT_GRA_BLOCK_SIZE
//                2307 : 2018/12/27 --> replace GRA_BLOCK_SIZE_Z by GRA_BLOCK_SIZE
//                2308a: 2019/01/24 --> add ELBDM_REMOVE_MOTION_CM
//                2308b: 2019/03/14 --> add OPT__RECORD_NOTE and OPT__RECORD_UNPHY
//                2309 : 2019/03/27 --> add OPT__FIXUP_ELECTRIC
//                2310 : 2019/04/20 --> add OPT__CK_INTERFACE_B
//                2311 : 2019/05/22 --> add OPT__CK_DIVERGENCE_B
//                2312 : 2019/05/31 --> add OPT__GRAVITY_EXTRA_MASS
//                2400 : 2019/06/08 --> output magnetic field for MHD
//                2401 : 2019/06/30 --> output OPT__FLAG_CURRENT and FlagTable_Current
//                2402 : 2019/07/17 --> replace USG_GhostSize by USG_GhostSizeF and USG_GhostSizeG
//                2403 : 2019/09/20 --> add BIT_REP_FLUX and BIT_REP_ELECTRIC defined in CUFLU.h
//                2404 : 2019/10/16 --> add DT__MAX
//                2405 : 2019/12/29 --> output GRACKLE_THREE_BODY_RATE, GRACKLE_CIE_COOLING, GRACKLE_H2_OPA_APPROX
//                2406 : 2020/02/26 --> output EOS
//                2407 : 2020/02/27 --> output MIN_EINT
//                2408 : 2020/05/10 --> output EOS_NAUX_MAX
//                2409 : 2020/07/05 --> output HLLC_WAVESPEED
//                2410 : 2020/07/31 --> output HLLE_WAVESPEED
//                2411 : 2020/08/11 --> output LR_EINT
//                2412 : 2020/08/17 --> output FLU_NIN_T
//                2413 : 2020/08/23 --> output HLLD_WAVESPEED
//                2414 : 2020/09/06 --> output INT_OPP_SIGN_0TH_ORDER
//                2415 : 2020/09/08 --> output OPT__LAST_RESORT_FLOOR
//                2416 : 2020/09/08 --> output BAROTROPIC_EOS
//                2417 : 2020/09/09 --> output ISO_TEMP
//                2418 : 2020/09/21 --> replace (OPT__GRAVITY_TYPE, OPT__EXTERNAL_POT) by
//                                      (OPT__SELF_GRAVITY, OPT__EXT_ACC, OPT__EXT_POT)
//                2419 : 2020/09/25 --> output EXTPOT_BLOCK_SIZE
//                2420 : 2020/10/12 --> output OPT__FLAG_USER_NUM and use variable-length datatype for FlagTable_User
//                2421 : 2020/10/26 --> output COSMIC_RAY
//                2422 : 2020/10/29 --> output the parameters of external potential table
//                2423 : 2020/11/01 --> output EOS_NTABLE_MAX
//                2424 : 2020/12/22 --> output SRC_USER
//                2425 : 2020/12/24 --> output SRC_DELEPTONIZATION
//                2426 : 2020/12/24 --> output FLU_NIN_S, FLU_NOUT_S, and SRC_GPU_NPGROUP
//                2427 : 2020/12/26 --> output SRC_BLOCK_SIZE and SRC_GHOST_SIZE
//                2428 : 2020/12/27 --> output SRC_NAUX_DLEP and SRC_NAUX_USER
//                2429 : 2021/01/26 --> output SRC_DLEP_PROF_NVAR and SRC_DLEP_PROF_NBINMAX
//                2430 : 2021/02/05 --> output EXT_POT_NGENE_MAX
//                2431 : 2021/02/13 --> output DER_GHOST_SIZE, DER_NXT, DER_NOUT_MAX, SRC_NXT
//                2432 : 2021/02/13 --> output OPT__OUTPUT_* of various derived fields
//                2433 : 2021/02/14 --> output MIN_TEMP
//                2434 : 2021/03/12 --> output OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, and PassiveIntFrac_VarIdx
//                2435 : 2021/04/06 --> output OPT__UM_IC_NLEVEL
//                2436 : 2021/04/06 --> output UM_IC_RefineRegion
//                2437 : 2021/05/12 --> output OPT__CHECK_PRES_AFTER_FLU
//                2438 : 2021/06/05 --> output git information
//                2439 : 2021/06/05 --> output UniqueDataID
//                2440 : 2021/06/17 --> output NFieldStored, NMagStored, and NFieldStoredMax
//                2441 : 2021/10/20 --> output OPT__FREEZE_FLUID
//                2442 : 2022/01/22 --> output OPT__FREEZE_PAR
//                2443 : 2022/01/30 --> output MINMOD_MAX_ITER and MONO_MAX_ITER
//                2444 : 2022/03/16 --> output OPT__FLAG_LOHNER_ENTR and MIN_ENTR
//                2445 : 2022/03/25 --> output OPT__OUTPUT_ENTR
//                2446 : 2022/05/10 --> output SUPPORT_LIBYT and LIBYT_USE_PATCH_GROUP
//                2447 : 2022/05/11 --> output MASSIVE_PARTICLES, TRACER, PAR_NTYPE, GhostSizeTracer
//                2448 : 2022/05/18 --> output PAR_IC_TYPE
//                2449 : 2021/03/21 --> output FEEDBACK
//                2450 : 2021/03/21 --> output FB_LEVEL
//                2451 : 2021/04/02 --> output FB_SNE and FB_USER
//                2452 : 2021/04/03 --> output FB_RSEED
//                2453 : 2022/07/08 --> output OPT__OUTPUT_RESTART
//                2454 : 2022/07/13 --> output OPT__INT_PRIM
//                2455 : 2022/10/10 --> output OPT__SAME_INTERFACE_B
//                2456 : 2022/10/17 --> output INTERP_MASK, OPT__CK_INPUT_FLUID
//                2457 : 2022/10/20 --> output RSOLVER_RESCUE
//                2458 : 2022/10/24 --> output AUTO_REDUCE_MINMOD_FACTOR, AUTO_REDUCE_MINMOD_MIN,
//                                             AUTO_REDUCE_INT_MONO_FACTOR, AUTO_REDUCE_INT_MONO_MIN,
//                                             INT_MONO_COEFF_B
//                2459 : 2022/11/04 --> output REFINE_NLEVEL
//                2460 : 2022/12/15 --> output SUPPORT_FFTW
//                2461 : 2023/01/28 --> output OPT__RESET_FLUID_INIT
//                2462 : 2023/03/19 --> output FB_GHOST_SIZE, FB_NXT
//                2463 : 2023/03/20 --> output FB_SEP_FLUOUT
//                2464 : 2023/04/27 --> output LIBYT_INTERACTIVE
//                2465 : 2023/04/29 --> output MU_NORM
//                2466 : 2023/05/08 --> output OPT__FFTW_STARTUP
//                2467 : 2023/05/18 --> replace OPT__INIT_BFIELD_BYFILE by OPT__INIT_BFIELD_BYVECPOT
//                2468 : 2023/06/24 --> output OPT__SORT_PATCH_BY_LBIDX
//                2469 : 2023/09/09 --> output MHM_CHECK_PREDICT
//                2470 : 2023/10/16 --> output OPT__OUTPUT_TEXT_FORMAT_FLT
//                2471 : 2023/11/09 --> output cosmic-ray options
//                2472 : 2023/11/11 --> output FixUpVar_Flux and FixUpVar_Restrict
//                2473 : 2023/11/29 --> output SRHD options and fields
//                2474 : 2023/11/22 --> output OPT__UM_IC_FLOAT8 and PAR_IC_FLOAT8
//                2475 : 2024/03/28 --> output YT_JUPYTER_USE_CONNECTION_FILE, LIBYT_RELOAD, LIBYT_INTERACTIVE,
//                                      LIBYT_JUPYTER
//                2476 : 2024/03/31 --> output particle attribute with assigned precision set by FLOAT8_PAR in Makefile;
//                                      record value of FLOAT8_PAR as Makefile.Float8_Par and KeyInfo.Float8_Par
//                2477 : 2024/04/05 --> output OPT__RECORD_CENTER, COM_CEN_X, COM_CEN_Y, COM_CEN_Z,
//                                             COM_MAX_R, COM_MIN_RHO, COM_TOLERR_R, COM_MAX_ITER
//                2478 : 2024/04/09 --> output ANGMOM_ORIGIN_X, ANGMOM_ORIGIN_Y, ANGMOM_ORIGIN_Z
//                2479 : 2024/04/23 --> output OPT__RES_PHASE, ELBDM_BASE_SPECTRAL, OPT__FLAG_INTERFERENCE, ELBDM_MATCH_PHASE,
//                                             ELBDM_FIRST_WAVE_LEVEL, OPT__LB_EXCHANGE_FATHER, FlagTable_Interference
//                                      output DENS and PHAS for the hybrid scheme (discard STUB)
//                                      output use_wave_flag[lv] for the hybrid scheme
//                2480 : 2024/07/17 --> output OPT__OUTPUT_PAR_MESH and particle attributes mapped from mesh quantities
//                2481 : 2024/12/11 --> output OPT__FLAG_ANGULAR, FlagTable_Angular, FLAG_ANGULAR_CEN_X, FLAG_ANGULAR_CEN_Y, FLAG_ANGULAR_CEN_Z
//                                             OPT__FLAG_RADIAL,  FlagTable_Radial,  FLAG_RADIAL_CEN_X,  FLAG_RADIAL_CEN_Y,  FLAG_RADIAL_CEN_Z
//                2500 : 2024/07/01 --> output particle integer attributes
//                2501 : 2025/01/15 --> output OPT__OUTPUT_TEXT_LENGTH_INT
//                2502 : 2025/01/16 --> output ConRef[]
//                2503 : 2025/01/17 --> output user-defined parameters in "User/UserPara" and
//                                             Input__TestProb parameters in "Info/InputTest"
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Total_HDF5( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d)     ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// check if the target file already exists
   if ( Aux_CheckFileExist(FileName)  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );



// 0. determine all the fields to be stored
//    --> must do it before calling GetCompound_* and FillIn_*
   const int NoDump = -__INT_MAX__; // must set the dump index to an extremely negative value (not -1) to disable it
   char FieldLabelOut[NFIELD_STORED_MAX][MAX_STRING];
   int  NFieldStored = 0;

   const int FluDumpIdx0 = NFieldStored;

   int NCompFluSkip = 0;
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
#     if (  ELBDM_SCHEME == ELBDM_HYBRID  &&  !defined( GAMER_DEBUG )  )
      if ( v == STUB )  // do not store STUB field unless we are in debug mode
      {
         NCompFluSkip += 1;
         continue;
      }
#     endif

      const int FluDumpIdx = NFieldStored++;
      if ( FluDumpIdx >= NFIELD_STORED_MAX )
         Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
      sprintf( FieldLabelOut[ FluDumpIdx ], "%s", FieldLabel[v] );
   }
   const int NCompStore  = NCOMP_TOTAL - NCompFluSkip;

#  ifdef GRAVITY
   const int PotDumpIdx = ( OPT__OUTPUT_POT ) ? NFieldStored++ : NoDump;
   if ( PotDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_POT )  sprintf( FieldLabelOut[PotDumpIdx], "%s", PotLabel );
#  endif

#  ifdef MASSIVE_PARTICLES
   const int ParDensDumpIdx = ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE ) ? NFieldStored++ : NoDump;
   if ( ParDensDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if      ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_PAR_ONLY )   sprintf( FieldLabelOut[ParDensDumpIdx], "%s", "ParDens"   );
   else if ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_TOTAL    )   sprintf( FieldLabelOut[ParDensDumpIdx], "%s", "TotalDens" );
#  endif

#  ifdef MHD
   const int CCMagDumpIdx0 = ( OPT__OUTPUT_CC_MAG ) ? NFieldStored : NoDump;
   if ( CCMagDumpIdx0+NCOMP_MAG-1 >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_CC_MAG )
   {
      NFieldStored += NCOMP_MAG;
      sprintf( FieldLabelOut[ CCMagDumpIdx0 + MAGX ], "%s", "CCMagX" );
      sprintf( FieldLabelOut[ CCMagDumpIdx0 + MAGY ], "%s", "CCMagY" );
      sprintf( FieldLabelOut[ CCMagDumpIdx0 + MAGZ ], "%s", "CCMagZ" );
   }
#  endif

#  if ( MODEL == HYDRO )
   const int PresDumpIdx   = ( OPT__OUTPUT_PRES ) ? NFieldStored++ : NoDump;
   if ( PresDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_PRES   )  sprintf( FieldLabelOut[PresDumpIdx  ], "%s", "Pres"   );

   const int TempDumpIdx   = ( OPT__OUTPUT_TEMP ) ? NFieldStored++ : NoDump;
   if ( TempDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_TEMP   )  sprintf( FieldLabelOut[TempDumpIdx  ], "%s", "Temp"   );

   const int EntrDumpIdx   = ( OPT__OUTPUT_ENTR ) ? NFieldStored++ : NoDump;
   if ( EntrDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_ENTR   )  sprintf( FieldLabelOut[EntrDumpIdx  ], "%s", "Entr"   );

   const int CsDumpIdx     = ( OPT__OUTPUT_CS ) ? NFieldStored++ : NoDump;
   if ( CsDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_CS     )  sprintf( FieldLabelOut[CsDumpIdx    ], "%s", "Cs"     );

   const int DivVelDumpIdx = ( OPT__OUTPUT_DIVVEL ) ? NFieldStored++ : NoDump;
   if ( DivVelDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_DIVVEL )  sprintf( FieldLabelOut[DivVelDumpIdx], "%s", "DivVel" );

   const int MachDumpIdx   = ( OPT__OUTPUT_MACH ) ? NFieldStored++ : NoDump;
   if ( MachDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_MACH   )  sprintf( FieldLabelOut[MachDumpIdx  ], "%s", "Mach"   );
#  endif

#  ifdef MHD
   const int DivMagDumpIdx = ( OPT__OUTPUT_DIVMAG ) ? NFieldStored++ : NoDump;
   if ( DivMagDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_DIVMAG )  sprintf( FieldLabelOut[DivMagDumpIdx], "%s", "DivMag" );
#  endif

#  ifdef SRHD
   const int LorentzDumpIdx = ( OPT__OUTPUT_LORENTZ ) ? NFieldStored++ : NoDump;
   if ( LorentzDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_LORENTZ )  sprintf( FieldLabelOut[LorentzDumpIdx], "%s", "Lrtz" );

   const int VelDumpIdx0 = ( OPT__OUTPUT_3VELOCITY ) ? NFieldStored : NoDump;
   if ( VelDumpIdx0+2 >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_3VELOCITY )
   {
      NFieldStored += 3;
      sprintf( FieldLabelOut[ VelDumpIdx0     ], "%s", "VelX" );
      sprintf( FieldLabelOut[ VelDumpIdx0 + 1 ], "%s", "VelY" );
      sprintf( FieldLabelOut[ VelDumpIdx0 + 2 ], "%s", "VelZ" );
   }

   const int EnthalpyDumpIdx = ( OPT__OUTPUT_ENTHALPY ) ? NFieldStored++ : NoDump;
   if ( EnthalpyDumpIdx >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_ENTHALPY )  sprintf( FieldLabelOut[EnthalpyDumpIdx], "%s", "Enth" );
#  endif // #ifdef SRHD

   const int UserDumpIdx0 = ( OPT__OUTPUT_USER_FIELD ) ? NFieldStored : NoDump;
   if ( UserDumpIdx0+UserDerField_Num-1 >= NFIELD_STORED_MAX )
      Aux_Error( ERROR_INFO, "exceed NFIELD_STORED_MAX (%d) !!\n", NFIELD_STORED_MAX );
   if ( OPT__OUTPUT_USER_FIELD )
   {
      NFieldStored += UserDerField_Num;
      for (int v=0; v<UserDerField_Num; v++)    sprintf( FieldLabelOut[ UserDumpIdx0 + v ], "%s", UserDerField_Label[v] );
   }



// 1. gather the number of patches at different MPI ranks and set the corresponding GID offset
   LB_PatchCount pc;
   LB_AllgatherPatchCount( pc );



// 2. prepare all HDF5 variables
   hsize_t H5_SetDims_LBIdx, H5_SetDims_Cr[2], H5_SetDims_Fa, H5_SetDims_Son, H5_SetDims_Sib[2], H5_SetDims_Field[4];
   hsize_t H5_MemDims_Field[4], H5_Count_Field[4], H5_Offset_Field[4];
   hid_t   H5_MemID_Field;
   hid_t   H5_FileID, H5_GroupID_Info, H5_GroupID_Tree, H5_GroupID_GridData, H5_GroupID_User;
   hid_t   H5_SetID_LBIdx, H5_SetID_Cr, H5_SetID_Fa, H5_SetID_Son, H5_SetID_Sib, H5_SetID_Field;
   hid_t   H5_SetID_KeyInfo, H5_SetID_Makefile, H5_SetID_SymConst, H5_SetID_InputPara, H5_SetID_InputTest, H5_SetID_UserPara;
   hid_t   H5_SpaceID_Scalar, H5_SpaceID_LBIdx, H5_SpaceID_Cr, H5_SpaceID_Fa, H5_SpaceID_Son, H5_SpaceID_Sib, H5_SpaceID_Field;
   hid_t   H5_TypeID_Com_KeyInfo, H5_TypeID_Com_Makefile, H5_TypeID_Com_SymConst, H5_TypeID_Com_InputPara, H5_TypeID_Com_InputTest, H5_TypeID_Com_UserPara;
   hid_t   H5_DataCreatePropList;
   hid_t   H5_AttID_Cvt2Phy;
   herr_t  H5_Status;
#  ifdef PARTICLE
   hsize_t H5_SetDims_NPar, H5_SetDims_ParData[1], H5_MemDims_ParData[1],  H5_Count_ParData[1], H5_Offset_ParData[1];
   hid_t   H5_SetID_NPar, H5_SpaceID_NPar, H5_SpaceID_ParData, H5_GroupID_Particle, H5_SetID_ParFltData, H5_SetID_ParIntData, H5_MemID_ParData;
#  endif
#  ifdef MHD
   hsize_t H5_SetDims_FCMag[4], H5_MemDims_FCMag[4], H5_Count_FCMag[4], H5_Offset_FCMag[4];
   hid_t   H5_MemID_FCMag, H5_SetID_FCMag, H5_SpaceID_FCMag[NCOMP_MAG];
#  endif

// 2-1. do NOT write fill values to any dataset for higher I/O performance
   H5_DataCreatePropList = H5Pcreate( H5P_DATASET_CREATE );
   H5_Status             = H5Pset_fill_time( H5_DataCreatePropList, H5D_FILL_TIME_NEVER );

// 2-2. create the "scalar" dataspace
   H5_SpaceID_Scalar = H5Screate( H5S_SCALAR );



// 3. output the simulation information
   if ( MPI_Rank == 0 )
   {
//    3-1. collect all information to be recorded
      KeyInfo_t      KeyInfo;
      Makefile_t     Makefile;
      SymConst_t     SymConst;
      InputPara_t    InputPara;
      HDF5_Output_t  HDF5_InputTest;
      HDF5_Output_t  HDF5_UserPara;

      FillIn_KeyInfo  ( KeyInfo, NFieldStored );
      FillIn_Makefile ( Makefile );
      FillIn_SymConst ( SymConst );
      FillIn_InputPara( InputPara, NFieldStored, FieldLabelOut );
      if ( Output_HDF5_InputTest_Ptr != NULL )  Output_HDF5_InputTest_Ptr( LOAD_HDF5_OUTPUT, NULL, &HDF5_InputTest );
      if ( Output_HDF5_UserPara_Ptr  != NULL )  Output_HDF5_UserPara_Ptr( &HDF5_UserPara );


//    3-2. create the "compound" datatype
      GetCompound_KeyInfo  ( H5_TypeID_Com_KeyInfo  );
      GetCompound_Makefile ( H5_TypeID_Com_Makefile );
      GetCompound_SymConst ( H5_TypeID_Com_SymConst );
      GetCompound_InputPara( H5_TypeID_Com_InputPara, NFieldStored );
      if ( Output_HDF5_InputTest_Ptr != NULL )  GetCompound_General( H5_TypeID_Com_InputTest, &HDF5_InputTest );
      if ( Output_HDF5_UserPara_Ptr  != NULL )  GetCompound_General( H5_TypeID_Com_UserPara, &HDF5_UserPara );


//    3-3. create the HDF5 file (overwrite the existing file)
      H5_FileID = H5Fcreate( FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to create the HDF5 file \"%s\" !!\n", FileName );


//    3-4. write the simulation info (note: dataset doesn't support VL datatype when the fill value is not defined)
      H5_GroupID_Info = H5Gcreate( H5_FileID, "Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_Info < 0 )    Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Info" );

//    3-4-1. KeyInfo
      H5_SetID_KeyInfo   = H5Dcreate( H5_GroupID_Info, "KeyInfo", H5_TypeID_Com_KeyInfo, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_KeyInfo < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "KeyInfo" );
      H5_Status          = H5Dwrite( H5_SetID_KeyInfo, H5_TypeID_Com_KeyInfo, H5S_ALL, H5S_ALL, H5P_DEFAULT, &KeyInfo );
      H5_Status          = H5Dclose( H5_SetID_KeyInfo );

//    3-4-2. Makefile
      H5_SetID_Makefile  = H5Dcreate( H5_GroupID_Info, "Makefile", H5_TypeID_Com_Makefile, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_Makefile < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Makefile" );
      H5_Status          = H5Dwrite( H5_SetID_Makefile, H5_TypeID_Com_Makefile, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Makefile );
      H5_Status          = H5Dclose( H5_SetID_Makefile );

//    3-4-3. SymConst
      H5_SetID_SymConst  = H5Dcreate( H5_GroupID_Info, "SymConst", H5_TypeID_Com_SymConst, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_SymConst < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "SymConst" );
      H5_Status          = H5Dwrite( H5_SetID_SymConst, H5_TypeID_Com_SymConst, H5S_ALL, H5S_ALL, H5P_DEFAULT, &SymConst );
      H5_Status          = H5Dclose( H5_SetID_SymConst );

//    3-4-4. InputPara
      H5_SetID_InputPara = H5Dcreate( H5_GroupID_Info, "InputPara", H5_TypeID_Com_InputPara, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_InputPara < 0 ) Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "InputPara" );
      H5_Status          = H5Dwrite( H5_SetID_InputPara, H5_TypeID_Com_InputPara, H5S_ALL, H5S_ALL, H5P_DEFAULT, &InputPara );
      H5_Status          = H5Dclose( H5_SetID_InputPara );

//    3-4-5. InputTest
      H5_SetID_InputTest = H5Dcreate( H5_GroupID_Info, "InputTest", H5_TypeID_Com_InputTest, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_InputTest < 0 ) Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "InputTest" );
      H5_Status          = H5_write_compound( H5_SetID_InputTest, H5_TypeID_Com_InputTest, &HDF5_InputTest );
      H5_Status          = H5Dclose( H5_SetID_InputTest );

      H5_Status = H5Gclose( H5_GroupID_Info );


//    3-5. write the user info
      H5_GroupID_User = H5Gcreate( H5_FileID, "User", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_User < 0 )     Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "User" );

//    3-5-1. UserPara
      if ( Output_HDF5_UserPara_Ptr != NULL )
      {
         H5_SetID_UserPara = H5Dcreate( H5_GroupID_User, "UserPara", H5_TypeID_Com_UserPara, H5_SpaceID_Scalar,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
         if ( H5_SetID_UserPara < 0 ) Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "UserPara" );
         H5_Status         = H5_write_compound( H5_SetID_UserPara, H5_TypeID_Com_UserPara, &HDF5_UserPara );
         H5_Status         = H5Dclose( H5_SetID_UserPara );
      } // if ( Output_HDF5_UserPara_Ptr != NULL )

      H5_Status = H5Gclose( H5_GroupID_User );
      H5_Status = H5Fclose( H5_FileID );

//    3-6. free memory
      for (int lv=0; lv<NLEVEL-1; lv++)   free( InputPara.FlagTable_User[lv].p );
   } // if ( MPI_Rank == 0 )



// 4. output the AMR tree structure (father, son, sibling, LBIdx, corner, and the number of particles --> sorted by GID)
   int root = 0;

// 4-1. allocate lists
   LB_LocalPatchExchangeList  lel;
   LB_GlobalPatchExchangeList gel( pc, root );

// 4-2. collect and sort LBIdx from all ranks
   LB_AllgatherLBIdx( pc, lel, &gel );

// 4-3. store the local tree
   LB_FillLocalPatchExchangeList( pc, lel );

// 4-4. gather data from all ranks
   LB_FillGlobalPatchExchangeList( pc, lel, gel, root );

// 4-5. dump the tree info
   if ( MPI_Rank == 0 )
   {
//    reopen file
      H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

      H5_GroupID_Tree = H5Gcreate( H5_FileID, "Tree", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_Tree < 0 )    Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Tree" );

//    4-5-1. LBIdx
      H5_SetDims_LBIdx = pc.NPatchAllLv;
      H5_SpaceID_LBIdx = H5Screate_simple( 1, &H5_SetDims_LBIdx, NULL );
      H5_SetID_LBIdx   = H5Dcreate( H5_GroupID_Tree, "LBIdx", H5T_NATIVE_LONG, H5_SpaceID_LBIdx,
                                    H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_LBIdx < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "LBIdx" );

      H5_Status = H5Dwrite( H5_SetID_LBIdx, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.LBIdxList_AllLv );
      H5_Status = H5Dclose( H5_SetID_LBIdx );
      H5_Status = H5Sclose( H5_SpaceID_LBIdx );

//    4-5-2. corner
      H5_SetDims_Cr[0] = pc.NPatchAllLv;
      H5_SetDims_Cr[1] = 3;
      H5_SpaceID_Cr    = H5Screate_simple( 2, H5_SetDims_Cr, NULL );
      H5_SetID_Cr      = H5Dcreate( H5_GroupID_Tree, "Corner", H5T_NATIVE_INT, H5_SpaceID_Cr,
                                    H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Cr < 0 )    Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Corner" );

//    attach the attribute for converting corner to physical coordinates
      H5_AttID_Cvt2Phy = H5Acreate( H5_SetID_Cr, "Cvt2Phy", H5T_NATIVE_DOUBLE, H5_SpaceID_Scalar,
                                    H5P_DEFAULT, H5P_DEFAULT );

      if ( H5_AttID_Cvt2Phy < 0 )   Aux_Error( ERROR_INFO, "failed to create the attribute \"%s\" !!\n", "Cvt2Phy" );

      H5_Status = H5Awrite( H5_AttID_Cvt2Phy, H5T_NATIVE_DOUBLE, &amr->dh[TOP_LEVEL] );
      H5_Status = H5Aclose( H5_AttID_Cvt2Phy );

      H5_Status = H5Dwrite( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.CrList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Cr );
      H5_Status = H5Sclose( H5_SpaceID_Cr );

//    4-5-3. father
      H5_SetDims_Fa = pc.NPatchAllLv;
      H5_SpaceID_Fa = H5Screate_simple( 1, &H5_SetDims_Fa, NULL );
      H5_SetID_Fa   = H5Dcreate( H5_GroupID_Tree, "Father", H5T_NATIVE_INT, H5_SpaceID_Fa,
                                 H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Fa < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Father" );

      H5_Status = H5Dwrite( H5_SetID_Fa, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.FaList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Fa );
      H5_Status = H5Sclose( H5_SpaceID_Fa );

//    4-5-4. son
      H5_SetDims_Son = pc.NPatchAllLv;
      H5_SpaceID_Son = H5Screate_simple( 1, &H5_SetDims_Son, NULL );
      H5_SetID_Son   = H5Dcreate( H5_GroupID_Tree, "Son", H5T_NATIVE_INT, H5_SpaceID_Son,
                                  H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Son < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Son" );

      H5_Status = H5Dwrite( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.SonList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Son );
      H5_Status = H5Sclose( H5_SpaceID_Son );

//    4-5-5. sibling
      H5_SetDims_Sib[0] = pc.NPatchAllLv;
      H5_SetDims_Sib[1] = 26;
      H5_SpaceID_Sib    = H5Screate_simple( 2, H5_SetDims_Sib, NULL );
      H5_SetID_Sib      = H5Dcreate( H5_GroupID_Tree, "Sibling", H5T_NATIVE_INT, H5_SpaceID_Sib,
                                     H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Sib < 0 )    Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Sibling" );

      H5_Status = H5Dwrite( H5_SetID_Sib, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.SibList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Sib );
      H5_Status = H5Sclose( H5_SpaceID_Sib );

//    4-5-6. NPar
#     ifdef PARTICLE
      H5_SetDims_NPar = pc.NPatchAllLv;
      H5_SpaceID_NPar = H5Screate_simple( 1, &H5_SetDims_NPar, NULL );
      H5_SetID_NPar   = H5Dcreate( H5_GroupID_Tree, "NPar", H5T_NATIVE_INT, H5_SpaceID_NPar,
                                   H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_NPar < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "NPar" );

      H5_Status = H5Dwrite( H5_SetID_NPar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.NParList_AllLv );
      H5_Status = H5Dclose( H5_SetID_NPar );
      H5_Status = H5Sclose( H5_SpaceID_NPar );
#     endif

//    close file
      H5_Status = H5Gclose( H5_GroupID_Tree );
      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )



// 5. output the simulation grid data (density, momentum, ... etc)
   const int FieldSizeOnePatch = sizeof(real)*CUBE(PS1);
   real (*FieldData)[PS1][PS1][PS1]  = NULL;

#  ifdef MHD
   const int FCMagSizeOnePatch = sizeof(real)*PS1P1*SQR(PS1);
   real (*FCMagData)[PS1P1*SQR(PS1)] = NULL;
#  endif

// 5-1. initialize the "GridData" group and the datasets of all fields and magnetic field
   H5_SetDims_Field[0] = pc.NPatchAllLv;
   H5_SetDims_Field[1] = PS1;
   H5_SetDims_Field[2] = PS1;
   H5_SetDims_Field[3] = PS1;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
      H5_SetDims_FCMag[0] = pc.NPatchAllLv;
      for (int t=1; t<4; t++)
      H5_SetDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_SpaceID_FCMag[v] = H5Screate_simple( 4, H5_SetDims_FCMag, NULL );
      if ( H5_SpaceID_FCMag[v] < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n",
                                                  "H5_SpaceID_FCMag", v );
   }
#  endif

   if ( MPI_Rank == 0 )
   {
//    HDF5 file must be synchronized before being written by the next rank
      SyncHDF5File( FileName );

      H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

//    create the "GridData" group
      H5_GroupID_GridData = H5Gcreate( H5_FileID, "GridData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "GridData" );

//    create the datasets of all fields
      for (int v=0; v<NFieldStored; v++)
      {
         H5_SetID_Field = H5Dcreate( H5_GroupID_GridData, FieldLabelOut[v], H5T_GAMER_REAL, H5_SpaceID_Field,
                                     H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", FieldLabelOut[v] );
         H5_Status = H5Dclose( H5_SetID_Field );
      }

//    create the datasets of all magnetic field components
#     ifdef MHD
      for (int v=0; v<NCOMP_MAG; v++)
      {
         H5_SetID_FCMag = H5Dcreate( H5_GroupID_GridData, MagLabel[v], H5T_GAMER_REAL, H5_SpaceID_FCMag[v],
                                     H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_FCMag < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", MagLabel[v] );
         H5_Status = H5Dclose( H5_SetID_FCMag );
      }
#     endif

//    close the file and group
      H5_Status = H5Gclose( H5_GroupID_GridData );
      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )


// 5-2. start to dump data (serial instead of parallel)
   const bool IntPhase_No         = false;
   const bool DE_Consistency_No   = false;
   const real MinDens_No          = -1.0;
   const real MinPres_No          = -1.0;
   const real MinTemp_No          = -1.0;
   const real MinEntr_No          = -1.0;
#  ifndef MHD
   const int  OPT__MAG_INT_SCHEME = INT_NONE;
#  endif
#  ifdef PARTICLE
   const bool TimingSendPar_No    = false;
   const bool PredictParPos_No    = false;   // particles synchronization is done in "Flu_CorrAfterAllSync()"
   const bool JustCountNPar_No    = false;
#  ifdef LOAD_BALANCE
   const bool SibBufPatch         = true;
   const bool FaSibBufPatch       = true;
#  else
   const bool SibBufPatch         = NULL_BOOL;
   const bool FaSibBufPatch       = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE

// for the derived fields
   const int Der_NP = 8;

   real (*Der_FluIn)[NCOMP_TOTAL][ CUBE(DER_NXT)            ] = new real [Der_NP][NCOMP_TOTAL ][ CUBE(DER_NXT)            ];
   real (*Der_Out  )             [ CUBE(PS1)                ] = new real         [DER_NOUT_MAX][ CUBE(PS1)                ];
#  ifdef MHD
   real (*Der_MagFC)[NCOMP_MAG  ][ (DER_NXT+1)*SQR(DER_NXT) ] = new real [Der_NP][NCOMP_MAG   ][ (DER_NXT+1)*SQR(DER_NXT) ];
   real (*Der_MagCC)             [ CUBE(DER_NXT)            ] = new real         [NCOMP_MAG   ][ CUBE(DER_NXT)            ];
#  else
   real (*Der_MagFC)[NCOMP_MAG  ][ (DER_NXT+1)*SQR(DER_NXT) ] = NULL;
   real (*Der_MagCC)             [ CUBE(DER_NXT)            ] = NULL;
#  endif
// allocate the maximum required memory (i.e., with NCOMP_TOTAL fields) for the temporary array Der_FluInTmp[]
   real (*Der_FluInTmp) = new real [ Der_NP*NCOMP_TOTAL*CUBE(DER_NXT) ];


// output one level at a time so that data at the same level are consecutive on disk (even for multiple ranks)
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    set the PID0 list for general purposes
      int *PID0List = new int [ amr->NPatchComma[lv][1]/8 ];
      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

//    5-2-0. initialize the particle density array (rho_ext) and collect particles from higher levels for outputting particle density
#     ifdef MASSIVE_PARTICLES
      if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
      {
         Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ, _PAR_TYPE, PredictParPos_No,
                                       NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

         Prepare_PatchData_InitParticleDensityArray( lv, Time[lv] );
      }
#     endif

      for (int TRank=0; TRank<MPI_NRank; TRank++)
      {
         if ( MPI_Rank == TRank )
         {
//          HDF5 file must be synchronized before being written by the next rank
            SyncHDF5File( FileName );

//          reopen the file and group
            H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
            if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

            H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
            if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );


//          5-2-1. dump cell-centered data
//          5-2-1-1. determine the memory space
            H5_MemDims_Field[0] = amr->NPatchComma[lv][1];
            H5_MemDims_Field[1] = PS1;
            H5_MemDims_Field[2] = PS1;
            H5_MemDims_Field[3] = PS1;

            H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
            if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_Field" );


//          5-2-1-2. determine the subset of the dataspace
            H5_Offset_Field[0] = pc.GID_Offset[lv];
            H5_Offset_Field[1] = 0;
            H5_Offset_Field[2] = 0;
            H5_Offset_Field[3] = 0;

            H5_Count_Field [0] = amr->NPatchComma[lv][1];
            H5_Count_Field [1] = PS1;
            H5_Count_Field [2] = PS1;
            H5_Count_Field [3] = PS1;

            H5_Status = H5Sselect_hyperslab( H5_SpaceID_Field, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
            if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the grid data !!\n" );


//          output one field at one level in one rank at a time
            FieldData = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];

            for (int v=0; v<NFieldStored; v++)
            {
//             5-2-1-3. collect the target field from all patches at the current target level
//             a. gravitational potential
#              ifdef GRAVITY
               if ( v == PotDumpIdx )
               {
                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                     memcpy( FieldData[PID], amr->patch[ amr->PotSg[lv] ][lv][PID]->pot, FieldSizeOnePatch );
               }
#              else
               if ( false ) {}
#              endif

//             b. particle density on grids
#              ifdef MASSIVE_PARTICLES
               else if ( v == ParDensDumpIdx )
               {
//                we do not check minimum density here (just because it's unnecessary)
                  Prepare_PatchData( lv, Time[lv], FieldData[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List,
                                     ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_PAR_ONLY ) ? _PAR_DENS : _TOTAL_DENS, _NONE,
                                     OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                                     MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
               }
#              endif

//             c. cell-centered magnetic field
#              ifdef MHD
               else if ( v >= CCMagDumpIdx0  &&  v < CCMagDumpIdx0+NCOMP_MAG )
               {
                  const int Bv = v - CCMagDumpIdx0;
                  real CCMag_1Cell[NCOMP_MAG];

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
//                   actually we only need CCMag_1Cell[Bv] here
//                   --> but the overhead of computing the other two B components is probably acceptable
                     MHD_GetCellCenteredBFieldInPatch( CCMag_1Cell, lv, PID, i, j, k, amr->MagSg[lv] );
                     FieldData[PID][k][j][i] = CCMag_1Cell[Bv];
                  }
               }
#              endif

//             d. derived fields
#              if ( MODEL == HYDRO )
//             d-1. gas pressure
               else if ( v == PresDumpIdx )
               {
//                we do not check minimum pressure here
                  Prepare_PatchData( lv, Time[lv], FieldData[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List,
                                     _PRES, _NONE, OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00,
                                     IntPhase_No, OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,
                                     DE_Consistency_No );
               }

//             d-2. gas temperature
               else if ( v == TempDumpIdx )
               {
                  const bool CheckMinTemp_No = false;

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     real u[NCOMP_TOTAL], Temp, Emag=NULL_REAL;

                     for (int v=0; v<NCOMP_TOTAL; v++)   u[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#                    endif
                     Temp = Hydro_Con2Temp( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                                            CheckMinTemp_No, NULL_REAL, Emag, EoS_DensEint2Temp_CPUPtr,
                                            EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     FieldData[PID][k][j][i] = Temp;
                  }
               } // if ( v == TempDumpIdx )

#              ifndef SRHD
//             d-3. gas entropy
               else if ( v == EntrDumpIdx )
               {
                  const bool CheckMinEntr_No = false;

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     real u[NCOMP_TOTAL], Entr, Emag=NULL_REAL;

                     for (int v=0; v<NCOMP_TOTAL; v++)   u[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#                    endif
                     Entr = Hydro_Con2Entr( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                                            CheckMinEntr_No, NULL_REAL, Emag, EoS_DensEint2Entr_CPUPtr,
                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     FieldData[PID][k][j][i] = Entr;
                  }
               } // if ( v == EntrDumpIdx )
#              endif

//             d-4. sound speed
               else if ( v == CsDumpIdx )
               {
                  const bool CheckMinPres_No = false;

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     real u[NCOMP_TOTAL], Pres, Cs2, Emag=NULL_REAL;

                     for (int v=0; v<NCOMP_TOTAL; v++)   u[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#                    endif

#                    ifdef SRHD
                     real Prim[NCOMP_TOTAL];
                     Hydro_Con2Pri( u, Prim, (real)-HUGE_NUMBER, NULL_BOOL, NULL_INT, NULL,
                                    NULL_BOOL, NULL_REAL, EoS_DensEint2Pres_CPUPtr,
                                    EoS_DensPres2Eint_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL, NULL );

                     Cs2 = EoS_DensPres2CSqr_CPUPtr( Prim[0], Prim[4], NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#                    else
                     Pres = Hydro_Con2Pres( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                                            CheckMinPres_No, NULL_REAL, Emag, EoS_DensEint2Pres_CPUPtr,
                                            EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
                     Cs2  = EoS_DensPres2CSqr_CPUPtr( u[DENS], Pres, u+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                      h_EoS_Table );
#                    endif
                     FieldData[PID][k][j][i] = SQRT( Cs2 );
                  }
               } // if ( v == CsDumpIdx )

//             d-5. divergence(velocity)
               else if ( v == DivVelDumpIdx )
               {
                  for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
                  {
//                   prepare the input density and momentum
//                   --> no need to prepare other fields
//                   --> store in Der_FluInTmp[] first and then copy to Der_FluIn[] since the shape of the latter
//                       is fixed to "[Der_NP][NCOMP_TOTAL][CUBE(DER_NXT)]" even though some fields may be useless
                     Prepare_PatchData( lv, Time[lv], Der_FluInTmp, NULL, DER_GHOST_SIZE, 1, &PID0,
                                        _DENS|_MOMX|_MOMY|_MOMZ, _NONE, OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26,
                                        IntPhase_No, OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,
                                        DE_Consistency_No );

//                   type casting for convenience
                     real (*Der_FluInTmp3D)[4][ CUBE(DER_NXT) ] = ( real (*)[4][ CUBE(DER_NXT) ] )Der_FluInTmp;

                     for (int LocalID=0; LocalID<8; LocalID++)
                     {
//                      copy data from Der_FluInTmp[] to Der_FluIn[]
                        const int Size1v = sizeof(real)*CUBE(DER_NXT);
                        memcpy( Der_FluIn[LocalID][DENS], Der_FluInTmp3D[LocalID][0], Size1v );
                        memcpy( Der_FluIn[LocalID][MOMX], Der_FluInTmp3D[LocalID][1], Size1v );
                        memcpy( Der_FluIn[LocalID][MOMY], Der_FluInTmp3D[LocalID][2], Size1v );
                        memcpy( Der_FluIn[LocalID][MOMZ], Der_FluInTmp3D[LocalID][3], Size1v );

//                      compute and store the target derived field
                        const int PID  = PID0 + LocalID;
                        const int NDer = 1;
                        Flu_DerivedField_DivVel( FieldData[PID][0][0], Der_FluIn[LocalID][0], NULL,
                                                 NDer, DER_NXT, DER_NXT, DER_NXT, DER_GHOST_SIZE, amr->dh[lv] );
                     }
                  } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
               } // if ( v == DivVelDumpIdx )

//             d-6. Mach number
               else if ( v == MachDumpIdx )
               {
                  for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
                  {
//                   prepare the input fields
//                   --> must prepare all NCOMP_TOTAL and NCOMP_MAG fields
                     Prepare_PatchData( lv, Time[lv], Der_FluIn[0][0], Der_MagFC[0][0], DER_GHOST_SIZE, 1, &PID0,
                                        _TOTAL, _MAG, OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCH, NSIDE_26,
                                        IntPhase_No, OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,
                                        DE_Consistency_No );

                     for (int LocalID=0; LocalID<8; LocalID++)
                     {
//                      convert B field from face-centered to cell-centered
#                       ifdef MHD
                        for (int k=0; k<DER_NXT; k++)
                        for (int j=0; j<DER_NXT; j++)
                        for (int i=0; i<DER_NXT; i++)
                        {
                           const int IdxCC = IDX321( i, j, k, DER_NXT, DER_NXT );
                           real B_CC[NCOMP_MAG];

                           MHD_GetCellCenteredBField( B_CC, Der_MagFC[LocalID][MAGX], Der_MagFC[LocalID][MAGY],
                                                      Der_MagFC[LocalID][MAGZ], DER_NXT, DER_NXT, DER_NXT, i, j, k );

                           Der_MagCC[MAGX][IdxCC] = B_CC[MAGX];
                           Der_MagCC[MAGY][IdxCC] = B_CC[MAGY];
                           Der_MagCC[MAGZ][IdxCC] = B_CC[MAGZ];
                        }
#                       endif // #ifdef MHD

//                      compute and store the target derived field
                        const int PID  = PID0 + LocalID;
                        const int NDer = 1;
                        Flu_DerivedField_Mach( FieldData[PID][0][0], Der_FluIn[LocalID][0], Der_MagCC[0],
                                               NDer, DER_NXT, DER_NXT, DER_NXT, DER_GHOST_SIZE, amr->dh[lv] );
                     } // for (int LocalID=0; LocalID<8; LocalID++)
                  } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
               } // if ( v == MachDumpIdx )

//             d-7. divergence(B field)
#              ifdef MHD
               else if ( v == DivMagDumpIdx )
               {
                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     const real DivB = MHD_GetCellCenteredDivBInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
                     FieldData[PID][k][j][i] = DivB;
                  }
               }
#              endif

#              ifdef SRHD
               else if (  ( v >= VelDumpIdx0 && v < VelDumpIdx0+3 )  ||  v == LorentzDumpIdx )
               {
                  const int vv = v - VelDumpIdx0 + 1;
                  real Prim[NCOMP_TOTAL], Cons[NCOMP_TOTAL], LorentzFactor;

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     for (int fv=0; fv<NCOMP_TOTAL; fv++)  Cons[fv] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[fv][k][j][i];

                     Hydro_Con2Pri( Cons, Prim, (real)-HUGE_NUMBER, false, NULL_INT, NULL,
                                    NULL_BOOL, (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL, &LorentzFactor );

//                   d-8. Lorentz factor
                     if      ( v == LorentzDumpIdx )
                        FieldData[PID][k][j][i] = LorentzFactor;

//                   d-9. 3-velocity
                     else if ( v >= VelDumpIdx0  &&  v < VelDumpIdx0+3 )
                        FieldData[PID][k][j][i] = Prim[vv] / LorentzFactor;
                  }
               }

//             d-10. reduced enthalpy
               else if ( v == EnthalpyDumpIdx )
               {
                  real Cons[NCOMP_TOTAL], HTilde;

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     for (int fv=0; fv<NCOMP_TOTAL; fv++)  Cons[fv] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[fv][k][j][i];

                     HTilde = Hydro_Con2HTilde( Cons, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     FieldData[PID][k][j][i] = HTilde;
                  }
               }
#              endif // #ifdef SRHD
#              endif // #if ( MODEL == HYDRO )

//             d-11. user-defined derived fields
//             the following check also works for OPT__OUTPUT_USER_FIELD==false since UserDerField_Num is initialized as 0
               else if ( v >= UserDumpIdx0  &&  v < UserDumpIdx0 + UserDerField_Num )
               {
                  for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
                  {
//                   prepare the input fields
//                   --> must prepare all NCOMP_TOTAL and NCOMP_MAG fields
                     Prepare_PatchData( lv, Time[lv], Der_FluIn[0][0], Der_MagFC[0][0], DER_GHOST_SIZE, 1, &PID0,
                                        _TOTAL, _MAG, OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCH, NSIDE_26,
                                        IntPhase_No, OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,
                                        DE_Consistency_No );

                     for (int LocalID=0; LocalID<8; LocalID++)
                     {
//                      convert B field from face-centered to cell-centered
#                       ifdef MHD
                        for (int k=0; k<DER_NXT; k++)
                        for (int j=0; j<DER_NXT; j++)
                        for (int i=0; i<DER_NXT; i++)
                        {
                           const int IdxCC = IDX321( i, j, k, DER_NXT, DER_NXT );
                           real B_CC[NCOMP_MAG];

                           MHD_GetCellCenteredBField( B_CC, Der_MagFC[LocalID][MAGX], Der_MagFC[LocalID][MAGY],
                                                      Der_MagFC[LocalID][MAGZ], DER_NXT, DER_NXT, DER_NXT, i, j, k );

                           Der_MagCC[MAGX][IdxCC] = B_CC[MAGX];
                           Der_MagCC[MAGY][IdxCC] = B_CC[MAGY];
                           Der_MagCC[MAGZ][IdxCC] = B_CC[MAGZ];
                        }
#                       endif // #ifdef MHD

//                      compute and store the target derived field
//                      --> store in Der_Out[] first and then copy to FieldData[] since we only output one field at a time
//                      --> there are redundant calculations in Flu_DerivedField_User_Ptr since it always computes
//                          UserDerField_Num fields while we only need one field at a time
//###OPTIMIZATION: only compute the derived field being dumped
                        const int PID  = PID0 + LocalID;
                        const int NDer = UserDerField_Num;
                        Flu_DerivedField_User_Ptr( Der_Out[0], Der_FluIn[LocalID][0], Der_MagCC[0],
                                                   NDer, DER_NXT, DER_NXT, DER_NXT, DER_GHOST_SIZE, amr->dh[lv] );

//                      copy data from Der_Der[] to FieldData[]
                        const int DerIdx = v - UserDumpIdx0;
                        memcpy( FieldData[PID], Der_Out[DerIdx], FieldSizeOnePatch );
                     } // for (int LocalID=0; LocalID<8; LocalID++)
                  } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
               } // if ( v >= UserDumpIdx0  &&  v < UserDumpIdx0 + UserDerField_Num )

//             e. fluid variables
               else if ( v >= FluDumpIdx0  &&  v < FluDumpIdx0+NCOMP_FLUID-NCompFluSkip )
               {
//                convert real/imag to density/phase in hybrid scheme
//                bitwise reproducibility currently fails in hybrid scheme because of conversion from RE/IM to DENS/PHAS when storing fields in HDF5
//                possible solution could be to convert RE/IM <-> DENS/PHAS using high-precision routines to ensure bitwise identity for significant digits
#                 if ( ELBDM_SCHEME == ELBDM_HYBRID )
                  if (  amr->use_wave_flag[lv]  &&  ( v == REAL || v == IMAG )  ) {
                     real Re, Im;
                     for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                     {
                        for (int k=0; k<PS1; k++)
                        for (int j=0; j<PS1; j++)
                        for (int i=0; i<PS1; i++)
                        {
                           Re = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i];
                           Im = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i];

                           if        ( v == REAL ) {
                              FieldData[PID][k][j][i] = SATAN2( Im, Re );
                           } else if ( v == IMAG ) {
                              FieldData[PID][k][j][i] = (real)0.0;
                           }
                        }
                     }

                  } else
#                 endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )
                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                     memcpy( FieldData[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
               } // if ( v >= FluDumpIdx0  &&  v < FluDumpIdx0+NCOMP_FLUID-NCompFluSkip )

//             f. passive fluid variables
               else if ( v >= FluDumpIdx0+NCOMP_FLUID-NCompFluSkip  &&  v < FluDumpIdx0+NCompStore )
               {
                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                     memcpy( FieldData[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v+NCompFluSkip], FieldSizeOnePatch );
               } // if ( v >= FluDumpIdx0+NCOMP_FLUID-NCompFluSkip  &&  v < FluDumpIdx0+NCompStore )

               else
                  Aux_Error( ERROR_INFO, "incorrect index (%d) !!\n", v );


//             5-2-1-4. write data to disk
               H5_SetID_Field = H5Dopen( H5_GroupID_GridData, FieldLabelOut[v], H5P_DEFAULT );

               H5_Status = H5Dwrite( H5_SetID_Field, H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT, FieldData );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to write a field (lv %d, v %d) !!\n", lv, v );

               H5_Status = H5Dclose( H5_SetID_Field );
            } // for (int v=0; v<NFieldStored; v++)


//          5-2-1-5.free resource before dumping magnetic field to save memory
            delete [] FieldData;

            H5_Status = H5Sclose( H5_MemID_Field );

//          free memory used for outputting particle density
#           ifdef MASSIVE_PARTICLES
            if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
            {
               Prepare_PatchData_FreeParticleDensityArray( lv );

               Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );
            }
#           endif


//          5-2-2. dump magnetic field
#           ifdef MHD
//          5-2-2-0. allocate memory
//                   --> output one B component at one level in one rank at a time
            FCMagData = new real [ amr->NPatchComma[lv][1] ][ PS1P1*SQR(PS1) ];

            for (int v=0; v<NCOMP_MAG; v++)
            {
//             5-2-2-1. determine the memory space
               H5_MemDims_FCMag[0] = amr->NPatchComma[lv][1];
               for (int t=1; t<4; t++)
               H5_MemDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

               H5_MemID_FCMag = H5Screate_simple( 4, H5_MemDims_FCMag, NULL );
               if ( H5_MemID_FCMag < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_FCMag" );


//             5-2-2-2. determine the subset of the dataspace
               H5_Offset_FCMag[0] = pc.GID_Offset[lv];
               H5_Offset_FCMag[1] = 0;
               H5_Offset_FCMag[2] = 0;
               H5_Offset_FCMag[3] = 0;

               H5_Count_FCMag [0] = amr->NPatchComma[lv][1];
               for (int t=1; t<4; t++)
               H5_Count_FCMag [t] = ( 3-t == v ) ? PS1P1 : PS1;

               H5_Status = H5Sselect_hyperslab( H5_SpaceID_FCMag[v], H5S_SELECT_SET, H5_Offset_FCMag, NULL, H5_Count_FCMag, NULL );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the magnetic field !!\n" );


//             5-2-2-3. collect the target B component from all patches at the current target level
               for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  memcpy( FCMagData[PID], amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[v], FCMagSizeOnePatch );


//             5-2-2-4. write data to disk
               H5_SetID_FCMag = H5Dopen( H5_GroupID_GridData, MagLabel[v], H5P_DEFAULT );

               H5_Status = H5Dwrite( H5_SetID_FCMag, H5T_GAMER_REAL, H5_MemID_FCMag, H5_SpaceID_FCMag[v], H5P_DEFAULT, FCMagData );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to write magnetic field (lv %d, v %d) !!\n", lv, v );

               H5_Status = H5Dclose( H5_SetID_FCMag );
               H5_Status = H5Sclose( H5_MemID_FCMag );
            } // for (int v=0; v<NCOMP_MAG; v++)


//          5-2-2-5.free resource
            delete [] FCMagData;
#           endif // #ifdef MHD

            H5_Status = H5Gclose( H5_GroupID_GridData );
            H5_Status = H5Fclose( H5_FileID );
         } // if ( MPI_Rank == TRank )

         MPI_Barrier( MPI_COMM_WORLD );

      } // for (int TRank=0; TRank<MPI_NRank; TRank++)

      delete [] PID0List;
   } // for (int lv=0; lv<NLEVEL; lv++)

   H5_Status = H5Sclose( H5_SpaceID_Field );
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   H5_Status = H5Sclose( H5_SpaceID_FCMag[v] );
#  endif

   delete [] Der_FluIn;
   delete [] Der_Out;
#  ifdef MHD
   delete [] Der_MagFC;
   delete [] Der_MagCC;
#  endif
   delete [] Der_FluInTmp;



// 6. output particles
#  ifdef PARTICLE
//###ISSUE: currently we output all particles at the same level at once (although one attribute at a time),
//          which may introduce a large memory overhead
//          --> solution: we can output a fixed number of particles at a time (see Output_DumpData_Total.cpp)
   long     (*NParLv_EachRank)[NLEVEL] = new long [MPI_NRank][NLEVEL];   // number of particles at each level in each rank
   real_par (*ParFltBuf1v1Lv)          = NULL;   // buffer storing the data of one particle floating-point attribute at one level
   long_par (*ParIntBuf1v1Lv)          = NULL;   // buffer storing the data of one particle integer        attribute at one level
   int        Par_NAtt_Mesh            = amr->Par->Mesh_Attr_Num;

   long  GParID_Offset[NLEVEL];  // GParID = global particle index (==> unique for each particle)
   long  NParLv_AllRank[NLEVEL];
   long  MaxNPar1Lv, NParInBuf, ParID;

// prepare particle attributes mapped from mesh quantities
   if ( OPT__OUTPUT_PAR_MESH )   Par_Output_TracerParticle_Mesh();


// 6-1. initialize variables
// 6-1-1. allocate I/O buffer for storing particle data
   MaxNPar1Lv = 0;
   for (int lv=0; lv<NLEVEL; lv++)  MaxNPar1Lv = MAX( MaxNPar1Lv, amr->Par->NPar_Lv[lv] );

   ParFltBuf1v1Lv = new real_par [MaxNPar1Lv];
   ParIntBuf1v1Lv = new long_par [MaxNPar1Lv];

// 6-1-2. get the starting global particle index (i.e., GParID_Offset[NLEVEL]) for particles at each level in this rank
   MPI_Allgather( amr->Par->NPar_Lv, NLEVEL, MPI_LONG, NParLv_EachRank[0], NLEVEL, MPI_LONG, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      NParLv_AllRank[lv] = 0;
      for (int r=0; r<MPI_NRank; r++)     NParLv_AllRank[lv] += NParLv_EachRank[r][lv];

      GParID_Offset[lv] = 0;
      for (int FaLv=0; FaLv<lv; FaLv++)   GParID_Offset[lv] += NParLv_AllRank[FaLv];

      for (int r=0; r<MPI_Rank; r++)      GParID_Offset[lv] += NParLv_EachRank[r][lv];
   }


// 6-2. initialize the "Particle" group and the datasets of all particle attributes
   H5_SetDims_ParData[0] = amr->Par->NPar_Active_AllRank;
   H5_SpaceID_ParData    = H5Screate_simple( 1, H5_SetDims_ParData, NULL );
   if ( H5_SpaceID_ParData < 0 )    Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_ParData" );

   if ( MPI_Rank == 0 )
   {
//    HDF5 file must be synchronized before being written by the next rank
      SyncHDF5File( FileName );

      H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

//    create the "Particle" group
      H5_GroupID_Particle = H5Gcreate( H5_FileID, "Particle", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_Particle < 0 )   Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Particle" );

//    create the datasets of all particle attributes
      for (int v=0; v<PAR_NATT_FLT_STORED+Par_NAtt_Mesh; v++)
      {
         char *ParLabel = ( v < PAR_NATT_FLT_STORED ) ? ParAttFltLabel[v] : amr->Par->Mesh_Attr_Label[v - PAR_NATT_FLT_STORED];

         H5_SetID_ParFltData = H5Dcreate( H5_GroupID_Particle, ParLabel, H5T_GAMER_REAL_PAR, H5_SpaceID_ParData,
                                          H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_ParFltData < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", ParLabel );
         H5_Status = H5Dclose( H5_SetID_ParFltData );
      }

      for (int v=0; v<PAR_NATT_INT_STORED; v++)
      {
         H5_SetID_ParIntData = H5Dcreate( H5_GroupID_Particle, ParAttIntLabel[v], H5T_GAMER_LONG_PAR, H5_SpaceID_ParData,
                                          H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_ParIntData < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", ParAttIntLabel[v] );
         H5_Status = H5Dclose( H5_SetID_ParIntData );
      }

//    close the file and group
      H5_Status = H5Gclose( H5_GroupID_Particle );
      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )


// 6-3. start to dump particle data (one level, one rank, and one attribute at a time)
//      --> note that particles must be outputted in the same order as their associated patches
   for (int lv=0; lv<NLEVEL; lv++)
   for (int TRank=0; TRank<MPI_NRank; TRank++)
   {
      if ( MPI_Rank == TRank )
      {
//       HDF5 file must be synchronized before being written by the next rank
         SyncHDF5File( FileName );

//       reopen the file and group
         H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
         if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_Particle = H5Gopen( H5_FileID, "Particle", H5P_DEFAULT );
         if ( H5_GroupID_Particle < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Particle" );


//       6-3-1. determine the memory space
         H5_MemDims_ParData[0] = amr->Par->NPar_Lv[lv];
         H5_MemID_ParData      = H5Screate_simple( 1, H5_MemDims_ParData, NULL );
         if ( H5_MemID_ParData < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_ParData" );


//       6-3-2. determine the subset of the dataspace
         H5_Offset_ParData[0] = GParID_Offset[lv];
         H5_Count_ParData [0] = amr->Par->NPar_Lv[lv];

         H5_Status = H5Sselect_hyperslab( H5_SpaceID_ParData, H5S_SELECT_SET, H5_Offset_ParData, NULL, H5_Count_ParData, NULL );
         if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the particle data !!\n" );


//       output one particle floating-point attribute at one level in one rank at a time
//       --> skip the last PAR_NATT_FLT_UNSTORED floating-point attributes since we do not want to store them on disk
         for (int v=0; v<PAR_NATT_FLT_STORED+Par_NAtt_Mesh; v++)
         {
//          6-3-3. collect particle data from all patches at the current target level
            NParInBuf = 0;

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            {
               ParID = amr->patch[0][lv][PID]->ParList[p];

#              ifdef DEBUG_PARTICLE
               if ( NParInBuf >= amr->Par->NPar_Lv[lv] )
                  Aux_Error( ERROR_INFO, "lv %d, NParInBuf (%ld) >= NPar_Lv (%ld) !!\n", lv, NParInBuf, amr->Par->NPar_Lv[lv] );
#              endif

               ParFltBuf1v1Lv[ NParInBuf ++ ] = ( v < PAR_NATT_FLT_STORED )
                                              ? amr->Par->AttributeFlt[v]                      [ParID]
                                              : amr->Par->Mesh_Attr   [v - PAR_NATT_FLT_STORED][ParID];
            }


//          6-3-4. write data to disk
            char *ParLabel = ( v < PAR_NATT_FLT_STORED ) ? ParAttFltLabel[v] : amr->Par->Mesh_Attr_Label[v - PAR_NATT_FLT_STORED];

            H5_SetID_ParFltData = H5Dopen( H5_GroupID_Particle, ParLabel, H5P_DEFAULT );

            H5_Status = H5Dwrite( H5_SetID_ParFltData, H5T_GAMER_REAL_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT, ParFltBuf1v1Lv );
            if ( H5_Status < 0 )
               Aux_Error( ERROR_INFO, "failed to write a particle floating-point attribute (lv %d, v %d) !!\n", lv, v );

            H5_Status = H5Dclose( H5_SetID_ParFltData );
         } // for (int v=0; v<PAR_NATT_FLT_STORED+Par_NAtt_Mesh; v++)

//       output one particle integer attribute at one level in one rank at a time
//       --> skip the last PAR_NATT_INT_UNSTORED integer attributes since we do not want to store them on disk
         for (int v=0; v<PAR_NATT_INT_STORED; v++)
         {
//          6-3-5. collect particle data from all patches at the current target level
            NParInBuf = 0;

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            {
               ParID = amr->patch[0][lv][PID]->ParList[p];

#              ifdef DEBUG_PARTICLE
               if ( NParInBuf >= amr->Par->NPar_Lv[lv] )
                  Aux_Error( ERROR_INFO, "lv %d, NParInBuf (%ld) >= NPar_Lv (%ld) !!\n", lv, NParInBuf, amr->Par->NPar_Lv[lv] );
#              endif

               ParIntBuf1v1Lv[ NParInBuf ++ ] = amr->Par->AttributeInt[v][ParID];
            }


//          6-3-6. write data to disk
            H5_SetID_ParIntData = H5Dopen( H5_GroupID_Particle, ParAttIntLabel[v], H5P_DEFAULT );

            H5_Status = H5Dwrite( H5_SetID_ParIntData, H5T_GAMER_LONG_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT, ParIntBuf1v1Lv );
            if ( H5_Status < 0 )
               Aux_Error( ERROR_INFO, "failed to write a particle integer attribute (lv %d, v %d) !!\n", lv, v );

            H5_Status = H5Dclose( H5_SetID_ParIntData );
         } // for (int v=0; v<PAR_NATT_INT_STORED; v++)

//       free resource
         H5_Status = H5Sclose( H5_MemID_ParData );
         H5_Status = H5Gclose( H5_GroupID_Particle );
         H5_Status = H5Fclose( H5_FileID );
      } // if ( MPI_Rank == TRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TRank=0; TRank<MPI_NRank; TRank++) ... for (int lv=0; lv<NLEVEL; lv++)

   H5_Status = H5Sclose( H5_SpaceID_ParData );

   delete [] ParFltBuf1v1Lv;
   delete [] ParIntBuf1v1Lv;
   delete [] NParLv_EachRank;
#  endif // #ifdef PARTICLE



// 7. check
#  ifdef DEBUG_HDF5
   if ( MPI_Rank == 0 )
   {
      const int MirrorSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };

      H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

//    7-1. validate the father-son relation
//    7-1-1. load data
      char SetName[MAX_STRING];
      sprintf( SetName, "%s", "Tree/Father" );
      H5_SetID_Fa = H5Dopen( H5_FileID, SetName, H5P_DEFAULT );

      if ( H5_SetID_Fa < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", SetName );

      H5_Status = H5Dread( H5_SetID_Fa, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.FaList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Fa );

      sprintf( SetName, "%s", "Tree/Son" );
      H5_SetID_Son = H5Dopen( H5_FileID, SetName, H5P_DEFAULT );

      if ( H5_SetID_Son < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", SetName );

      H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.SonList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Son );

      sprintf( SetName, "%s", "Tree/Sibling" );
      H5_SetID_Sib = H5Dopen( H5_FileID, SetName, H5P_DEFAULT );

      if ( H5_SetID_Sib < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", SetName );

      H5_Status = H5Dread( H5_SetID_Sib, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gel.SibList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Sib );


      for (int lv=0; lv<NLEVEL; lv++)
      for (int GID=pc.GID_LvStart[lv]; GID<pc.GID_LvStart[lv]+NPatchTotal[lv]; GID++)
      {
//       7-1-2. root patches have no father
         if ( lv == 0 )
         if ( gel.FaList_AllLv[GID] != -1 )
            Aux_Error( ERROR_INFO, "Lv %d, GID %d, FaGID %d != -1 !!\n", lv, GID, gel.FaList_AllLv[GID] );

//       7-1-3. all patches at refinement levels have fathers
         if ( lv > 0 )
         if ( gel.FaList_AllLv[GID] < 0  ||  gel.FaList_AllLv[GID] >= pc.GID_LvStart[lv] )
            Aux_Error( ERROR_INFO, "Lv %d, GID %d, FaGID %d < 0 (or > max = %d) !!\n",
                       lv, GID, gel.FaList_AllLv[GID], pc.GID_LvStart[lv]-1 );

//       7-1-4. father->son == itself
         if ( lv > 0 )
         if ( gel.SonList_AllLv[ gel.FaList_AllLv[GID] ] + GID%8 != GID )
            Aux_Error( ERROR_INFO, "Lv %d, GID %d, FaGID %d, FaGID->Son %d ==> inconsistent !!\n",
                       lv, GID, gel.FaList_AllLv[GID], gel.SonList_AllLv[ gel.FaList_AllLv[GID] ] );

//       7-1-5. son->father == itself
         const int SonGID = gel.SonList_AllLv[GID];
         if ( SonGID != -1 )
         {
            if ( lv >= MAX_LEVEL )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d != -1 !!\n", lv, GID, SonGID );

            if ( SonGID < -1 )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d < -1 !!\n", lv, GID, SonGID );

            if ( lv < NLEVEL-1  &&  SonGID >= pc.GID_LvStart[lv+1]+NPatchTotal[lv+1] )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d > max %d !!\n", lv, GID, SonGID,
                          pc.GID_LvStart[lv+1]+NPatchTotal[lv+1]-1 );

            for (int LocalID=0; LocalID<8; LocalID++)
            if ( gel.FaList_AllLv[SonGID+LocalID] != GID )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d, SonGID->Father %d ==> inconsistent !!\n",
                          lv, GID, SonGID+LocalID, gel.FaList_AllLv[SonGID+LocalID] );
         }

//       7-1-6. sibling->sibling_mirror = itself
         for (int s=0; s<26; s++)
         {
            const int SibGID = gel.SibList_AllLv[GID][s];

            if ( SibGID >= 0 )
            {
               if ( SibGID < pc.GID_LvStart[lv]  ||  SibGID >= pc.GID_LvStart[lv]+NPatchTotal[lv] )
                  Aux_Error( ERROR_INFO, "Lv %d, GID %d, sib %d, SibGID %d lies outside the correct range (%d <= SibGID < %d) !!\n",
                             lv, GID, s, SibGID, pc.GID_LvStart[lv], pc.GID_LvStart[lv]+NPatchTotal[lv] );

               if ( gel.SibList_AllLv[SibGID][ MirrorSib[s] ] != GID )
                  Aux_Error( ERROR_INFO, "Lv %d, GID %d, sib %d, SibGID %d != SibGID->sibling %d !!\n",
                             lv, GID, s, SibGID, gel.SibList_AllLv[SibGID][ MirrorSib[s] ] );
            }
         }
      } // for (int lv=0; lv<NLEVEL; lv++)

      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )
#  endif // #ifdef DEBUG_HDF5



// 8. close all HDF5 objects and free memory
   if ( MPI_Rank == 0 )
   {
      H5_Status = H5Tclose( H5_TypeID_Com_KeyInfo   );
      H5_Status = H5Tclose( H5_TypeID_Com_Makefile  );
      H5_Status = H5Tclose( H5_TypeID_Com_SymConst  );
      H5_Status = H5Tclose( H5_TypeID_Com_InputPara );
      if ( Output_HDF5_InputTest_Ptr != NULL )   H5_Status = H5Tclose( H5_TypeID_Com_InputTest );
      if ( Output_HDF5_UserPara_Ptr  != NULL )   H5_Status = H5Tclose( H5_TypeID_Com_UserPara );
   } // if ( MPI_Rank == 0 )
   H5_Status = H5Sclose( H5_SpaceID_Scalar );
   H5_Status = H5Pclose( H5_DataCreatePropList );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d)     ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Total_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_KeyInfo
// Description :  Fill in the KeyInfo_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  KeyInfo      : Pointer to be filled in
//                NFieldStored : Number of grid fields to be stored on disk
//-------------------------------------------------------------------------------------------------------
void FillIn_KeyInfo( KeyInfo_t &KeyInfo, const int NFieldStored )
{

   const time_t CalTime = time( NULL );   // calendar time

   KeyInfo.FormatVersion        = 2503;
   KeyInfo.Model                = MODEL;
   KeyInfo.NLevel               = NLEVEL;
   KeyInfo.NCompFluid           = NCOMP_FLUID;
   KeyInfo.NCompPassive         = NCOMP_PASSIVE;
   KeyInfo.PatchSize            = PS1;
   KeyInfo.DumpID               = DumpID;
   KeyInfo.Step                 = Step;
#  ifdef GRAVITY
   KeyInfo.AveDens_Init         = AveDensity_Init;
   KeyInfo.Gravity              = 1;
#  else
   KeyInfo.Gravity              = 0;
#  endif
#  ifdef PARTICLE
   KeyInfo.Particle             = 1;
#  else
   KeyInfo.Particle             = 0;
#  endif  // #ifdef PARTICLE ... else ...
#  ifdef FLOAT8
   KeyInfo.Float8               = 1;
#  else
   KeyInfo.Float8               = 0;
#  endif
   KeyInfo.NFieldStored         = NFieldStored;
   KeyInfo.NMagStored           = NCOMP_MAG;
#  ifdef PARTICLE
   KeyInfo.Par_NPar             = amr->Par->NPar_Active_AllRank;
   KeyInfo.Par_NAttFltStored    = PAR_NATT_FLT_STORED;
   KeyInfo.Par_NAttIntStored    = PAR_NATT_INT_STORED;
#  ifdef FLOAT8_PAR
   KeyInfo.Float8_Par           = 1;
#  else
   KeyInfo.Float8_Par           = 0;
#  endif
#  ifdef INT8_PAR
   KeyInfo.Int8_Par             = 1;
#  else
   KeyInfo.Int8_Par             = 0;
#  endif
#  endif // #ifdef PARTICLE
#  if ( MODEL == HYDRO )
#  ifdef MHD
   KeyInfo.Magnetohydrodynamics = 1;
#  else
   KeyInfo.Magnetohydrodynamics = 0;
#  endif
#  ifdef SRHD
   KeyInfo.SRHydrodynamics      = 1;
#  else
   KeyInfo.SRHydrodynamics      = 0;
#  endif
#  ifdef COSMIC_RAY
   KeyInfo.CosmicRay            = 1;
#  ifdef CR_DIFFUSION
   KeyInfo.CR_Diffusion         = 1;
#  else
   KeyInfo.CR_Diffusion         = 0;
#  endif
#  else // #ifdef COSMIC_RAY
   KeyInfo.CosmicRay            = 0;
#  endif // #ifdef COSMIC_RAY .. else ...
#  endif // #if ( MODEL == HYDRO )

   for (int d=0; d<3; d++)
   {
      KeyInfo.NX0     [d] = NX0_TOT      [d];
      KeyInfo.BoxScale[d] = amr->BoxScale[d];
      KeyInfo.BoxSize [d] = amr->BoxSize [d];
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {
      KeyInfo.Time          [lv] = Time              [lv];
      KeyInfo.CellSize      [lv] = amr->dh           [lv];
      KeyInfo.CellScale     [lv] = amr->scale        [lv];
      KeyInfo.NPatch        [lv] = NPatchTotal       [lv];
      KeyInfo.AdvanceCounter[lv] = AdvanceCounter    [lv];
      KeyInfo.dTime_AllLv   [lv] = dTime_AllLv       [lv];
#     if ( MODEL == ELBDM  &&  ELBDM_SCHEME == ELBDM_HYBRID )
      KeyInfo.UseWaveScheme [lv] = amr->use_wave_flag[lv];
#     endif
   }

   KeyInfo.CodeVersion  = (char*)VERSION;
   KeyInfo.DumpWallTime = ctime( &CalTime );
   KeyInfo.DumpWallTime[ strlen(KeyInfo.DumpWallTime)-1 ] = '\0';  // remove the last character '\n'
   KeyInfo.GitBranch    = (char*)EXPAND_AND_QUOTE( GIT_BRANCH );
   KeyInfo.GitCommit    = (char*)EXPAND_AND_QUOTE( GIT_COMMIT );

//###REVISE: replace rand() by UUID
   srand( time(NULL) );
   KeyInfo.UniqueDataID = rand();

   if ( !ConRefInitialized )   Aux_Error( ERROR_INFO, "Reference values for conserved variables have not been assigned yet !!\n" );

   for (int v=0; v<1+NCONREF_MAX; v++)   KeyInfo.ConRef[v] = ConRef[v];

} // FUNCTION : FillIn_KeyInfo



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_Makefile
// Description :  Fill in the Makefile_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  Makefile : Pointer to be filled in
//-------------------------------------------------------------------------------------------------------
void FillIn_Makefile( Makefile_t &Makefile )
{

// model-independent options
   Makefile.Model                  = MODEL;

#  ifdef GRAVITY
   Makefile.Gravity                = 1;
#  else
   Makefile.Gravity                = 0;
#  endif

#  ifdef COMOVING
   Makefile.Comoving               = 1;
#  else
   Makefile.Comoving               = 0;
#  endif

#  ifdef PARTICLE
   Makefile.Particle               = 1;
#  else
   Makefile.Particle               = 0;
#  endif

#  ifdef GPU
   Makefile.UseGPU                 = 1;
#  else
   Makefile.UseGPU                 = 0;
#  endif

#  ifdef GAMER_DEBUG
   Makefile.GAMER_Debug            = 1;
#  else
   Makefile.GAMER_Debug            = 0;
#  endif

#  ifdef BITWISE_REPRODUCIBILITY
   Makefile.BitwiseReproducibility = 1;
#  else
   Makefile.BitwiseReproducibility = 0;
#  endif

#  ifdef TIMING
   Makefile.Timing                 = 1;
#  else
   Makefile.Timing                 = 0;
#  endif

#  ifdef TIMING_SOLVER
   Makefile.TimingSolver           = 1;
#  else
   Makefile.TimingSolver           = 0;
#  endif

#  ifdef FLOAT8
   Makefile.Float8                 = 1;
#  else
   Makefile.Float8                 = 0;
#  endif

#  ifdef SERIAL
   Makefile.Serial                 = 1;
#  else
   Makefile.Serial                 = 0;
#  endif

#  ifdef LOAD_BALANCE
   Makefile.LoadBalance            = LOAD_BALANCE;
#  else
   Makefile.LoadBalance            = 0;
#  endif

#  ifdef OVERLAP_MPI
   Makefile.OverlapMPI             = 1;
#  else
   Makefile.OverlapMPI             = 0;
#  endif

#  ifdef OPENMP
   Makefile.OpenMP                 = 1;
#  else
   Makefile.OpenMP                 = 0;
#  endif

#  ifdef GPU
   Makefile.GPU_Arch               = GPU_ARCH;
#  else
   Makefile.GPU_Arch               = NULL_INT;
#  endif

#  ifdef LAOHU
   Makefile.Laohu                  = 1;
#  else
   Makefile.Laohu                  = 0;
#  endif

#  ifdef SUPPORT_HDF5
   Makefile.SupportHDF5            = 1;
#  else
   Makefile.SupportHDF5            = 0;
#  endif

#  ifdef SUPPORT_GSL
   Makefile.SupportGSL             = 1;
#  else
   Makefile.SupportGSL             = 0;
#  endif

#  ifdef SUPPORT_SPECTRAL_INT
   Makefile.SupportSpectralInt     = 1;
#  else
   Makefile.SupportSpectralInt     = 0;
#  endif

#  ifdef SUPPORT_FFTW
   Makefile.SupportFFTW            = SUPPORT_FFTW;
#  else
   Makefile.SupportFFTW            = 0;
#  endif

#  ifdef SUPPORT_LIBYT
   Makefile.SupportLibYT           = 1;

#  ifdef LIBYT_USE_PATCH_GROUP
   Makefile.LibYTUsePatchGroup     = 1;
#  else
   Makefile.LibYTUsePatchGroup     = 0;
#  endif

#  ifdef LIBYT_INTERACTIVE
   Makefile.LibYTInteractive       = 1;
#  else
   Makefile.LibYTInteractive       = 0;
#  endif

#  ifdef LIBYT_RELOAD
   Makefile.LibYTReload            = 1;
#  else
   Makefile.LibYTReload            = 0;
#  endif

#  ifdef LIBYT_JUPYTER
   Makefile.LibYTJupyter           = 1;
#  else
   Makefile.LibYTJupyter           = 0;
#  endif

#  else  // #ifdef SUPPORT_LIBYT
   Makefile.SupportLibYT           = 0;
#  endif // #ifdef SUPPORT_LIBYT ... else ...

#  ifdef SUPPORT_GRACKLE
   Makefile.SupportGrackle         = 1;
#  else
   Makefile.SupportGrackle         = 0;
#  endif

   Makefile.RandomNumber           = RANDOM_NUMBER;

   Makefile.NLevel                 = NLEVEL;
   Makefile.MaxPatch               = MAX_PATCH;


// model-dependent options
#  ifdef GRAVITY
   Makefile.PotScheme              = POT_SCHEME;

#  ifdef STORE_POT_GHOST
   Makefile.StorePotGhost          = 1;
#  else
   Makefile.StorePotGhost          = 0;
#  endif

#  ifdef UNSPLIT_GRAVITY
   Makefile.UnsplitGravity         = 1;
#  else
   Makefile.UnsplitGravity         = 0;
#  endif

#  endif // #ifdef GRAVITY

#  if   ( MODEL == HYDRO )
   Makefile.FluScheme              = FLU_SCHEME;

#  ifdef LR_SCHEME
   Makefile.LRScheme               = LR_SCHEME;
#  endif

#  ifdef RSOLVER
   Makefile.RSolver                = RSOLVER;
#  endif

#  ifdef DUAL_ENERGY
   Makefile.DualEnergy             = DUAL_ENERGY;
#  else
   Makefile.DualEnergy             = 0;
#  endif

#  ifdef MHD
   Makefile.Magnetohydrodynamics   = 1;
#  else
   Makefile.Magnetohydrodynamics   = 0;
#  endif

#  ifdef SRHD
   Makefile.SRHydrodynamics        = 1;
#  else
   Makefile.SRHydrodynamics        = 0;
#  endif

#  ifdef COSMIC_RAY
   Makefile.CosmicRay              = 1;
#  ifdef CR_DIFFUSION
   Makefile.CR_Diffusion           = 1;
#  else
   Makefile.CR_Diffusion           = 0;
#  endif
#  else // #ifdef COSMIC_RAY
   Makefile.CosmicRay              = 0;
#  endif // #ifdef COSMIC_RAY .. else ...

   Makefile.EoS                    = EOS;

#  ifdef BAROTROPIC_EOS
   Makefile.BarotropicEoS          = 1;
#  else
   Makefile.BarotropicEoS          = 0;
#  endif


#  elif ( MODEL == ELBDM )

   Makefile.ELBDMScheme            = ELBDM_SCHEME;
   Makefile.WaveScheme             = WAVE_SCHEME;

#  ifdef CONSERVE_MASS
   Makefile.ConserveMass           = 1;
#  else
   Makefile.ConserveMass           = 0;
#  endif

#  ifdef LAPLACIAN_4TH
   Makefile.Laplacian4th           = 1;
#  else
   Makefile.Laplacian4th           = 0;
#  endif

#  ifdef QUARTIC_SELF_INTERACTION
   Makefile.SelfInteraction4       = 1;
#  else
   Makefile.SelfInteraction4       = 0;
#  endif


#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
#  ifdef MASSIVE_PARTICLES
   Makefile.MassiveParticles       = 1;
#  else
   Makefile.MassiveParticles       = 0;
#  endif
#  ifdef TRACER
   Makefile.Tracer                 = 1;
#  else
   Makefile.Tracer                 = 0;
#  endif
#  ifdef STORE_PAR_ACC
   Makefile.StoreParAcc            = 1;
#  else
   Makefile.StoreParAcc            = 0;
#  endif

#  ifdef STAR_FORMATION
   Makefile.StarFormation          = 1;
#  else
   Makefile.StarFormation          = 0;
#  endif

#  ifdef FEEDBACK
   Makefile.Feedback               = 1;
#  else
   Makefile.Feedback               = 0;
#  endif

   Makefile.Par_NAttFltUser        = PAR_NATT_FLT_USER;
   Makefile.Par_NAttIntUser        = PAR_NATT_INT_USER;

#  ifdef FLOAT8_PAR
   Makefile.Float8_Par             = 1;
#  else
   Makefile.Float8_Par             = 0;
#  endif

#  ifdef INT8_PAR
   Makefile.Int8_Par               = 1;
#  else
   Makefile.Int8_Par               = 0;
#  endif
#  endif // #ifdef PARTICLE

} // FUNCTION : FillIn_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_SymConst
// Description :  Fill in the SymConst_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  SymConst : Pointer to be filled in
//-------------------------------------------------------------------------------------------------------
void FillIn_SymConst( SymConst_t &SymConst )
{

// model-independent variables
   SymConst.NCompFluid           = NCOMP_FLUID;
   SymConst.NCompPassive         = NCOMP_PASSIVE;
   SymConst.PatchSize            = PS1;
   SymConst.Flu_NIn              = FLU_NIN;
   SymConst.Flu_NOut             = FLU_NOUT;
   SymConst.Flu_NIn_T            = FLU_NIN_T;
   SymConst.Flu_NIn_S            = FLU_NIN_S;
   SymConst.Flu_NOut_S           = FLU_NOUT_S;
   SymConst.NFluxFluid           = NFLUX_FLUID;
   SymConst.NFluxPassive         = NFLUX_PASSIVE;
   SymConst.Flu_GhostSize        = FLU_GHOST_SIZE;
   SymConst.Flu_Nxt              = FLU_NXT;
#  ifdef DEBUG_HDF5
   SymConst.Debug_HDF5           = 1;
#  else
   SymConst.Debug_HDF5           = 0;
#  endif
   SymConst.SibOffsetNonperiodic = SIB_OFFSET_NONPERIODIC;
#  ifdef LOAD_BALANCE
   SymConst.SonOffsetLB          = SON_OFFSET_LB;
#  endif

   SymConst.TinyNumber           = TINY_NUMBER;
   SymConst.HugeNumber           = HUGE_NUMBER;
   SymConst.MaxError             = MAX_ERROR;


// model-dependent variables
#  ifdef GRAVITY
   SymConst.Gra_NIn              = GRA_NIN;
   SymConst.Pot_GhostSize        = POT_GHOST_SIZE;
   SymConst.Gra_GhostSize        = GRA_GHOST_SIZE;
   SymConst.Rho_GhostSize        = RHO_GHOST_SIZE;
   SymConst.Pot_Nxt              = POT_NXT;
   SymConst.Gra_Nxt              = GRA_NXT;
   SymConst.Rho_Nxt              = RHO_NXT;

#  ifdef UNSPLIT_GRAVITY
   SymConst.USG_GhostSizeF       = USG_GHOST_SIZE_F;
   SymConst.USG_GhostSizeG       = USG_GHOST_SIZE_G;
   SymConst.USG_NxtF             = USG_NXT_F;
   SymConst.USG_NxtG             = USG_NXT_G;
#  endif

   SymConst.ExtPot_BlockSize     = EXTPOT_BLOCK_SIZE;
   SymConst.Gra_BlockSize        = GRA_BLOCK_SIZE;
   SymConst.ExtPotNAuxMax        = EXT_POT_NAUX_MAX;
   SymConst.ExtAccNAuxMax        = EXT_ACC_NAUX_MAX;
   SymConst.ExtPotNGeneMax       = EXT_POT_NGENE_MAX;

#  if   ( POT_SCHEME == SOR )
   SymConst.Pot_BlockSize_z      = POT_BLOCK_SIZE_Z;
#  ifdef SOR_RHO_SHARED
   SymConst.SOR_RhoShared        = 1;
#  else
   SymConst.SOR_RhoShared        = 0;
#  endif
#  ifdef SOR_CPOT_SHARED
   SymConst.SOR_CPotShared       = 1;
#  else
   SymConst.SOR_CPotShared       = 0;
#  endif
#  ifdef SOR_USE_SHUFFLE
   SymConst.SOR_UseShuffle       = 1;
#  else
   SymConst.SOR_UseShuffle       = 0;
#  endif
#  ifdef SOR_USE_PADDING
   SymConst.SOR_UsePadding       = 1;
#  else
   SymConst.SOR_UsePadding       = 0;
#  endif
   SymConst.SOR_ModReduction     = SOR_MOD_REDUCTION;

#  elif ( POT_SCHEME == MG  )
   SymConst.Pot_BlockSize_x      = POT_BLOCK_SIZE_X;
#  endif // POT_SCHEME
#  endif // #ifdef GRAVITY


#  ifdef PARTICLE
   SymConst.Par_NAttFltStored    = PAR_NATT_FLT_STORED;
   SymConst.Par_NAttIntStored    = PAR_NATT_INT_STORED;
   SymConst.Par_NType            = PAR_NTYPE;
#  ifdef GRAVITY
   SymConst.RhoExt_GhostSize     = RHOEXT_GHOST_SIZE;
#  endif

#  ifdef DEBUG_PARTICLE
   SymConst.Debug_Particle       = 1;
#  else
   SymConst.Debug_Particle       = 0;
#  endif

   SymConst.ParList_GrowthFactor = PARLIST_GROWTH_FACTOR;
   SymConst.ParList_ReduceFactor = PARLIST_REDUCE_FACTOR;
#  endif // #ifdef PARTICLE


#  ifdef BIT_REP_FLUX
   SymConst.BitRep_Flux          = 1;
#  else
   SymConst.BitRep_Flux          = 0;
#  endif

#  ifdef MHD
#  ifdef BIT_REP_ELECTRIC
   SymConst.BitRep_Electric      = 1;
#  else
   SymConst.BitRep_Electric      = 0;
#  endif
#  endif // #ifdef MHD

#  ifdef INTERP_MASK
   SymConst.InterpMask           = 1;
#  else
   SymConst.InterpMask           = 0;
#  endif

#  ifdef FB_SEP_FLUOUT
   SymConst.FB_SepFluOut         = 1;
#  else
   SymConst.FB_SepFluOut         = 0;
#  endif


#  if   ( MODEL == HYDRO )
   SymConst.Flu_BlockSize_x      = FLU_BLOCK_SIZE_X;
   SymConst.Flu_BlockSize_y      = FLU_BLOCK_SIZE_Y;
#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   SymConst.CheckUnphyInFluid = 1;
#  else
   SymConst.CheckUnphyInFluid = 0;
#  endif
#  ifdef CHAR_RECONSTRUCTION
   SymConst.CharReconstruction   = 1;
#  else
   SymConst.CharReconstruction   = 0;
#  endif
#  ifdef LR_EINT
   SymConst.LR_Eint              = 1;
#  else
   SymConst.LR_Eint              = 0;
#  endif
#  ifdef CHECK_INTERMEDIATE
   SymConst.CheckIntermediate    = CHECK_INTERMEDIATE;
#  else
   SymConst.CheckIntermediate    = 0;
#  endif
#  ifdef RSOLVER_RESCUE
   SymConst.RSolverRescue        = RSOLVER_RESCUE;
#  else
   SymConst.RSolverRescue        = 0;
#  endif
#  ifdef HLL_NO_REF_STATE
   SymConst.HLL_NoRefState       = 1;
#  else
   SymConst.HLL_NoRefState       = 0;
#  endif
#  ifdef HLL_INCLUDE_ALL_WAVES
   SymConst.HLL_IncludeAllWaves  = 1;
#  else
   SymConst.HLL_IncludeAllWaves  = 0;
#  endif
   SymConst.HLLC_WaveSpeed       = HLLC_WAVESPEED;
   SymConst.HLLE_WaveSpeed       = HLLE_WAVESPEED;
#  ifdef MHD
   SymConst.HLLD_WaveSpeed       = HLLD_WAVESPEED;
#  endif
#  ifdef N_FC_VAR
   SymConst.N_FC_Var             = N_FC_VAR;
#  endif
#  ifdef N_SLOPE_PPM
   SymConst.N_Slope_PPM          = N_SLOPE_PPM;
#  endif
#  ifdef MHD
#  ifdef EULERY
   SymConst.EulerY               = 1;
#  else
   SymConst.EulerY               = 0;
#  endif
#  endif // MHD
#  ifdef MHM_CHECK_PREDICT
   SymConst.MHM_CheckPredict     = 1;
#  else
   SymConst.MHM_CheckPredict     = 0;
#  endif
   SymConst.EoSNAuxMax           = EOS_NAUX_MAX;
   SymConst.EoSNTableMax         = EOS_NTABLE_MAX;

#  elif  ( MODEL == ELBDM )
   SymConst.Flu_BlockSize_x      = FLU_BLOCK_SIZE_X;
   SymConst.Flu_BlockSize_y      = FLU_BLOCK_SIZE_Y;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   SymConst.Flu_HJ_BlockSize_y   = FLU_HJ_BLOCK_SIZE_Y;
#  endif
#  if ( WAVE_SCHEME == WAVE_GRAMFE )
   SymConst.GramFEScheme         = GRAMFE_SCHEME;
   SymConst.GramFEGamma          = GRAMFE_GAMMA;
   SymConst.GramFEG              = GRAMFE_G;
   SymConst.GramFENDelta         = GRAMFE_NDELTA;
   SymConst.GramFEOrder          = GRAMFE_ORDER;
   SymConst.GramFEND             = GRAMFE_ND;
   SymConst.GramFEFluNxt         = GRAMFE_FLU_NXT;
#  endif
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   SymConst.dt_Flu_BlockSize     = DT_FLU_BLOCK_SIZE;
#  ifdef DT_FLU_USE_SHUFFLE
   SymConst.dt_Flu_UseShuffle    = 1;
#  else
   SymConst.dt_Flu_UseShuffle    = 0;
#  endif
#  ifdef GRAVITY
   SymConst.dt_Gra_BlockSize     = DT_GRA_BLOCK_SIZE;
#  ifdef DT_GRA_USE_SHUFFLE
   SymConst.dt_Gra_UseShuffle    = 1;
#  else
   SymConst.dt_Gra_UseShuffle    = 0;
#  endif
#  endif

   SymConst.Src_BlockSize        = SRC_BLOCK_SIZE;
   SymConst.Src_GhostSize        = SRC_GHOST_SIZE;
   SymConst.Src_Nxt              = SRC_NXT;
#  if ( MODEL == HYDRO )
   SymConst.Src_NAuxDlep         = SRC_NAUX_DLEP;
   SymConst.Src_DlepProfNVar     = SRC_DLEP_PROF_NVAR;
   SymConst.Src_DlepProfNBinMax  = SRC_DLEP_PROF_NBINMAX;
#  endif
   SymConst.Src_NAuxUser         = SRC_NAUX_USER;

   SymConst.Der_GhostSize        = DER_GHOST_SIZE;
   SymConst.Der_Nxt              = DER_NXT;
   SymConst.Der_NOut_Max         = DER_NOUT_MAX;

#  ifdef FEEDBACK
   SymConst.FB_GhostSize         = FB_GHOST_SIZE;
   SymConst.FB_Nxt               = FB_NXT;
#  endif

   SymConst.NFieldStoredMax      = NFIELD_STORED_MAX;

   SymConst.NConRefMax           = NCONREF_MAX;

} // FUNCTION : FillIn_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_InputPara
// Description :  Fill in the InputPara_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  InputPara     : Pointer to be filled in
//                NFieldStored  : Number of grid fields to be stored on disk
//                FieldLabelOut : Field labels
//-------------------------------------------------------------------------------------------------------
void FillIn_InputPara( InputPara_t &InputPara, const int NFieldStored, char FieldLabelOut[][MAX_STRING] )
{

// simulation scale
   InputPara.BoxSize                 = BOX_SIZE;
   for (int d=0; d<3; d++)
   InputPara.NX0_Tot[d]              = NX0_TOT[d];
   InputPara.MPI_NRank               = MPI_NRank;
   for (int d=0; d<3; d++)
   InputPara.MPI_NRank_X[d]          = MPI_NRank_X[d];
   InputPara.OMP_NThread             = OMP_NTHREAD;
   InputPara.EndT                    = END_T;
   InputPara.EndStep                 = END_STEP;

// test problems
   InputPara.TestProb_ID             = TESTPROB_ID;

// code units
   InputPara.Opt__Unit               = OPT__UNIT;
   InputPara.Unit_L                  = UNIT_L;
   InputPara.Unit_M                  = UNIT_M;
   InputPara.Unit_T                  = UNIT_T;
   InputPara.Unit_V                  = UNIT_V;
   InputPara.Unit_D                  = UNIT_D;
   InputPara.Unit_E                  = UNIT_E;
   InputPara.Unit_P                  = UNIT_P;
#  ifdef MHD
   InputPara.Unit_B                  = UNIT_B;
#  endif

// boundary condition
   for (int t=0; t<6; t++)
   InputPara.Opt__BC_Flu[t]          = OPT__BC_FLU[t];
#  ifdef GRAVITY
   InputPara.Opt__BC_Pot             = OPT__BC_POT;
   InputPara.GFunc_Coeff0            = GFUNC_COEFF0;
#  endif

// particle
#  ifdef PARTICLE
   InputPara.Par_Init                = amr->Par->Init;
   InputPara.Par_ICFormat            = amr->Par->ParICFormat;
   InputPara.Par_ICMass              = amr->Par->ParICMass;
   InputPara.Par_ICType              = amr->Par->ParICType;
   InputPara.Par_ICFloat8            = PAR_IC_FLOAT8;
   InputPara.Par_ICInt8              = PAR_IC_INT8;
   InputPara.Par_Interp              = amr->Par->Interp;
   InputPara.Par_InterpTracer        = amr->Par->InterpTracer;
   InputPara.Par_Integ               = amr->Par->Integ;
   InputPara.Par_IntegTracer         = amr->Par->IntegTracer;
   InputPara.Par_ImproveAcc          = amr->Par->ImproveAcc;
   InputPara.Par_PredictPos          = amr->Par->PredictPos;
   InputPara.Par_TracerVelCorr       = amr->Par->TracerVelCorr;
   InputPara.Par_RemoveCell          = amr->Par->RemoveCell;
   InputPara.Opt__FreezePar          = OPT__FREEZE_PAR;
   InputPara.Par_GhostSize           = amr->Par->GhostSize;
   InputPara.Par_GhostSizeTracer     = amr->Par->GhostSizeTracer;
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
   InputPara.ParAttFltLabel[v]       = ParAttFltLabel[v];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
   InputPara.ParAttIntLabel[v]       = ParAttIntLabel[v];
#  endif

// cosmology
#  ifdef COMOVING
   InputPara.A_Init                  = A_INIT;
   InputPara.OmegaM0                 = OMEGA_M0;
   InputPara.Hubble0                 = HUBBLE0;
#  endif

// time-step determination
   InputPara.Dt__Max                 = DT__MAX;
   InputPara.Dt__Fluid               = DT__FLUID;
   InputPara.Dt__FluidInit           = DT__FLUID_INIT;
#  ifdef GRAVITY
   InputPara.Dt__Gravity             = DT__GRAVITY;
#  endif
#  if ( MODEL == ELBDM )
   InputPara.Dt__Phase               = DT__PHASE;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   InputPara.Dt__HybridCFL           = DT__HYBRID_CFL;
   InputPara.Dt__HybridCFLInit       = DT__HYBRID_CFL_INIT;
   InputPara.Dt__HybridVelocity      = DT__HYBRID_VELOCITY;
   InputPara.Dt__HybridVelocityInit  = DT__HYBRID_VELOCITY_INIT;
#  endif
#  endif // #if ( MODEL == ELBDM )
#  ifdef PARTICLE
   InputPara.Dt__ParVel              = DT__PARVEL;
   InputPara.Dt__ParVelMax           = DT__PARVEL_MAX;
   InputPara.Dt__ParAcc              = DT__PARACC;
#  endif
#  ifdef CR_DIFFUSION
   InputPara.Dt__CR_Diffusion        = DT__CR_DIFFUSION;
#  endif
#  ifdef COMOVING
   InputPara.Dt__MaxDeltaA           = DT__MAX_DELTA_A;
#  endif
   InputPara.Dt__SyncParentLv        = DT__SYNC_PARENT_LV;
   InputPara.Dt__SyncChildrenLv      = DT__SYNC_CHILDREN_LV;
   InputPara.Opt__DtUser             = OPT__DT_USER;
   InputPara.Opt__DtLevel            = OPT__DT_LEVEL;
   InputPara.Opt__RecordDt           = OPT__RECORD_DT;
   InputPara.AutoReduceDt            = AUTO_REDUCE_DT;
   InputPara.AutoReduceDtFactor      = AUTO_REDUCE_DT_FACTOR;
   InputPara.AutoReduceDtFactorMin   = AUTO_REDUCE_DT_FACTOR_MIN;
#  if ( MODEL == HYDRO )
   InputPara.AutoReduceMinModFactor  = AUTO_REDUCE_MINMOD_FACTOR;
   InputPara.AutoReduceMinModMin     = AUTO_REDUCE_MINMOD_MIN;
#  endif
   InputPara.AutoReduceIntMonoFactor = AUTO_REDUCE_INT_MONO_FACTOR;
   InputPara.AutoReduceIntMonoMin    = AUTO_REDUCE_INT_MONO_MIN;

// domain refinement
   InputPara.RegridCount             = REGRID_COUNT;
   InputPara.RefineNLevel            = REFINE_NLEVEL;
   InputPara.FlagBufferSize          = FLAG_BUFFER_SIZE;
   InputPara.FlagBufferSizeMaxM1Lv   = FLAG_BUFFER_SIZE_MAXM1_LV;
   InputPara.FlagBufferSizeMaxM2Lv   = FLAG_BUFFER_SIZE_MAXM2_LV;
   InputPara.MaxLevel                = MAX_LEVEL;
   InputPara.Opt__Flag_Rho           = OPT__FLAG_RHO;
   InputPara.Opt__Flag_RhoGradient   = OPT__FLAG_RHO_GRADIENT;
   InputPara.Opt__Flag_Angular       = OPT__FLAG_ANGULAR;
   InputPara.FlagAngular_CenX        = FLAG_ANGULAR_CEN_X;
   InputPara.FlagAngular_CenY        = FLAG_ANGULAR_CEN_Y;
   InputPara.FlagAngular_CenZ        = FLAG_ANGULAR_CEN_Z;
   InputPara.Opt__Flag_Radial        = OPT__FLAG_RADIAL;
   InputPara.FlagRadial_CenX         = FLAG_RADIAL_CEN_X;
   InputPara.FlagRadial_CenY         = FLAG_RADIAL_CEN_Y;
   InputPara.FlagRadial_CenZ         = FLAG_RADIAL_CEN_Z;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Flag_PresGradient  = OPT__FLAG_PRES_GRADIENT;
   InputPara.Opt__Flag_Vorticity     = OPT__FLAG_VORTICITY;
   InputPara.Opt__Flag_Jeans         = OPT__FLAG_JEANS;
#  ifdef MHD
   InputPara.Opt__Flag_Current       = OPT__FLAG_CURRENT;
#  endif
#  ifdef SRHD
   InputPara.Dt__SpeedOfLight        = DT__SPEED_OF_LIGHT;
   InputPara.Opt__Flag_LrtzGradient  = OPT__FLAG_LRTZ_GRADIENT;
#  endif
#  ifdef COSMIC_RAY
   InputPara.Opt__Flag_CRay          = OPT__FLAG_CRAY;
#  endif
#  endif
#  if ( MODEL == ELBDM )
   InputPara.Opt__Flag_EngyDensity   = OPT__FLAG_ENGY_DENSITY;
   InputPara.Opt__Flag_Spectral      = OPT__FLAG_SPECTRAL;
   InputPara.Opt__Flag_Spectral_N    = OPT__FLAG_SPECTRAL_N;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   InputPara.Opt__Flag_Interference  = OPT__FLAG_INTERFERENCE;
#  endif
#  endif // #if ( MODEL == ELBDM )
   InputPara.Opt__Flag_LohnerDens    = OPT__FLAG_LOHNER_DENS;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Flag_LohnerEngy    = OPT__FLAG_LOHNER_ENGY;
   InputPara.Opt__Flag_LohnerPres    = OPT__FLAG_LOHNER_PRES;
   InputPara.Opt__Flag_LohnerTemp    = OPT__FLAG_LOHNER_TEMP;
   InputPara.Opt__Flag_LohnerEntr    = OPT__FLAG_LOHNER_ENTR;
#  ifdef COSMIC_RAY
   InputPara.Opt__Flag_LohnerCRay    = OPT__FLAG_LOHNER_CRAY;
#  endif
#  endif
   InputPara.Opt__Flag_LohnerForm    = OPT__FLAG_LOHNER_FORM;
   InputPara.Opt__Flag_User          = OPT__FLAG_USER;
   InputPara.Opt__Flag_User_Num      = OPT__FLAG_USER_NUM;
   InputPara.Opt__Flag_Region        = OPT__FLAG_REGION;
#  ifdef PARTICLE
   InputPara.Opt__Flag_NParPatch     = OPT__FLAG_NPAR_PATCH;
   InputPara.Opt__Flag_NParCell      = OPT__FLAG_NPAR_CELL;
   InputPara.Opt__Flag_ParMassCell   = OPT__FLAG_PAR_MASS_CELL;
#  endif
   InputPara.Opt__NoFlagNearBoundary = OPT__NO_FLAG_NEAR_BOUNDARY;
   InputPara.Opt__PatchCount         = OPT__PATCH_COUNT;
#  ifdef PARTICLE
   InputPara.Opt__ParticleCount      = OPT__PARTICLE_COUNT;
#  endif
   InputPara.Opt__ReuseMemory        = OPT__REUSE_MEMORY;
   InputPara.Opt__MemoryPool         = OPT__MEMORY_POOL;

// load balance
#  ifdef LOAD_BALANCE
   InputPara.LB_WLI_Max              = amr->LB->WLI_Max;
#  ifdef PARTICLE
   InputPara.LB_Par_Weight           = amr->LB->Par_Weight;
#  endif
   InputPara.Opt__RecordLoadBalance  = OPT__RECORD_LOAD_BALANCE;
   InputPara.Opt__LB_ExchangeFather  = OPT__LB_EXCHANGE_FATHER;
#  endif
   InputPara.Opt__MinimizeMPIBarrier = OPT__MINIMIZE_MPI_BARRIER;

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   InputPara.Gamma                   = GAMMA;
   InputPara.MolecularWeight         = MOLECULAR_WEIGHT;
   InputPara.MuNorm                  = MU_NORM;
   InputPara.IsoTemp                 = ISO_TEMP;
   InputPara.MinMod_Coeff            = MINMOD_COEFF;
   InputPara.MinMod_MaxIter          = MINMOD_MAX_ITER;
   InputPara.Opt__LR_Limiter         = OPT__LR_LIMITER;
   InputPara.Opt__1stFluxCorr        = OPT__1ST_FLUX_CORR;
   InputPara.Opt__1stFluxCorrScheme  = OPT__1ST_FLUX_CORR_SCHEME;
#  ifdef DUAL_ENERGY
   InputPara.DualEnergySwitch        = DUAL_ENERGY_SWITCH;
#  endif
#  ifdef MHD
   InputPara.Opt__SameInterfaceB     = OPT__SAME_INTERFACE_B;
#  endif
#  endif // HYDRO

// ELBDM solvers
#  if ( MODEL == ELBDM )
   InputPara.ELBDM_Mass              = ELBDM_MASS;
   InputPara.ELBDM_PlanckConst       = ELBDM_PLANCK_CONST;
#  ifdef QUARTIC_SELF_INTERACTION
   InputPara.ELBDM_Lambda            = ELBDM_LAMBDA;
#  endif
   InputPara.ELBDM_Taylor3_Coeff     = ELBDM_TAYLOR3_COEFF;
   InputPara.ELBDM_Taylor3_Auto      = ELBDM_TAYLOR3_AUTO;
   InputPara.ELBDM_RemoveMotionCM    = ELBDM_REMOVE_MOTION_CM;
   InputPara.ELBDM_BaseSpectral      = ELBDM_BASE_SPECTRAL;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   InputPara.ELBDM_FirstWaveLevel    = ELBDM_FIRST_WAVE_LEVEL;
#  endif
#  endif // ELBDM

// fluid solvers in different models
   InputPara.Flu_GPU_NPGroup         = FLU_GPU_NPGROUP;
   InputPara.GPU_NStream             = GPU_NSTREAM;
   InputPara.Opt__FixUp_Flux         = OPT__FIXUP_FLUX;
   InputPara.FixUpFlux_Var           = FixUpVar_Flux;
#  ifdef MHD
   InputPara.Opt__FixUp_Electric     = OPT__FIXUP_ELECTRIC;
#  endif
   InputPara.Opt__FixUp_Restrict     = OPT__FIXUP_RESTRICT;
   InputPara.FixUpRestrict_Var       = FixUpVar_Restrict;
   InputPara.Opt__CorrAfterAllSync   = OPT__CORR_AFTER_ALL_SYNC;
   InputPara.Opt__NormalizePassive   = OPT__NORMALIZE_PASSIVE;

   InputPara.NormalizePassive_NVar   = PassiveNorm_NVar;

   for (int v=0; v<NCOMP_PASSIVE; v++)
   InputPara.NormalizePassive_VarIdx[v] = PassiveNorm_VarIdx[v];

   InputPara.Opt__IntFracPassive_LR  = OPT__INT_FRAC_PASSIVE_LR;

   InputPara.IntFracPassive_NVar     = PassiveIntFrac_NVar;

   for (int v=0; v<NCOMP_PASSIVE; v++)
   InputPara.IntFracPassive_VarIdx[v] = PassiveIntFrac_VarIdx[v];

   for (int v=0; v<NFieldStored; v++)
   InputPara.FieldLabel[v]           = FieldLabelOut[v];

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   InputPara.MagLabel[v]             = MagLabel[v];
#  endif

   InputPara.Opt__OverlapMPI         = OPT__OVERLAP_MPI;
   InputPara.Opt__ResetFluid         = OPT__RESET_FLUID;
   InputPara.Opt__ResetFluidInit     = OPT__RESET_FLUID_INIT;
   InputPara.Opt__FreezeFluid        = OPT__FREEZE_FLUID;
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   InputPara.MinDens                 = MIN_DENS;
#  endif
#  if ( MODEL == HYDRO )
   InputPara.MinPres                 = MIN_PRES;
   InputPara.MinEint                 = MIN_EINT;
   InputPara.MinTemp                 = MIN_TEMP;
   InputPara.MinEntr                 = MIN_ENTR;
   InputPara.Opt__CheckPresAfterFlu  = OPT__CHECK_PRES_AFTER_FLU,
   InputPara.Opt__LastResortFloor    = OPT__LAST_RESORT_FLOOR;
   InputPara.JeansMinPres            = JEANS_MIN_PRES;
   InputPara.JeansMinPres_Level      = JEANS_MIN_PRES_LEVEL;
   InputPara.JeansMinPres_NCell      = JEANS_MIN_PRES_NCELL;
#  endif

// self-gravity
#  ifdef GRAVITY
   InputPara.NewtonG                 = NEWTON_G;
#  if   ( POT_SCHEME == SOR )
   InputPara.SOR_Omega               = SOR_OMEGA;
   InputPara.SOR_MaxIter             = SOR_MAX_ITER;
   InputPara.SOR_MinIter             = SOR_MIN_ITER;
#  elif ( POT_SCHEME == MG )
   InputPara.MG_MaxIter              = MG_MAX_ITER;
   InputPara.MG_NPreSmooth           = MG_NPRE_SMOOTH;
   InputPara.MG_NPostSmooth          = MG_NPOST_SMOOTH;
   InputPara.MG_ToleratedError       = MG_TOLERATED_ERROR;
#  endif
   InputPara.Pot_GPU_NPGroup         = POT_GPU_NPGROUP;
   InputPara.Opt__GraP5Gradient      = OPT__GRA_P5_GRADIENT;
   InputPara.Opt__SelfGravity        = OPT__SELF_GRAVITY;
   InputPara.Opt__ExtAcc             = OPT__EXT_ACC;
   InputPara.Opt__ExtPot             = OPT__EXT_POT;
   InputPara.ExtPotTable_Name        = EXT_POT_TABLE_NAME;
   for (int d=0; d<3; d++)
   InputPara.ExtPotTable_NPoint[d]   = EXT_POT_TABLE_NPOINT[d];
   for (int d=0; d<3; d++)
   InputPara.ExtPotTable_dh[d]       = EXT_POT_TABLE_DH[d];
   for (int d=0; d<3; d++)
   InputPara.ExtPotTable_EdgeL[d]    = EXT_POT_TABLE_EDGEL[d];
   InputPara.ExtPotTable_Float8      = EXT_POT_TABLE_FLOAT8;
   InputPara.Opt__GravityExtraMass   = OPT__GRAVITY_EXTRA_MASS;
#  endif

// source terms
   InputPara.Src_Deleptonization     = SrcTerms.Deleptonization;
   InputPara.Src_User                = SrcTerms.User;
   InputPara.Src_GPU_NPGroup         = SRC_GPU_NPGROUP;

// Grackle
#  ifdef SUPPORT_GRACKLE
   InputPara.Grackle_Activate        = GRACKLE_ACTIVATE;
   InputPara.Grackle_Verbose         = GRACKLE_VERBOSE;
   InputPara.Grackle_Cooling         = GRACKLE_COOLING;
   InputPara.Grackle_Primordial      = GRACKLE_PRIMORDIAL;
   InputPara.Grackle_Metal           = GRACKLE_METAL;
   InputPara.Grackle_UV              = GRACKLE_UV;
   InputPara.Grackle_CMB_Floor       = GRACKLE_CMB_FLOOR;
   InputPara.Grackle_PE_Heating      = GRACKLE_PE_HEATING;
   InputPara.Grackle_PE_HeatingRate  = GRACKLE_PE_HEATING_RATE;
   InputPara.Grackle_CloudyTable     = GRACKLE_CLOUDY_TABLE;
   InputPara.Grackle_ThreeBodyRate   = GRACKLE_THREE_BODY_RATE;
   InputPara.Grackle_CIE_Cooling     = GRACKLE_CIE_COOLING;
   InputPara.Grackle_H2_OpaApprox    = GRACKLE_H2_OPA_APPROX;
   InputPara.Che_GPU_NPGroup         = CHE_GPU_NPGROUP;
#  endif

// star formation
#  ifdef STAR_FORMATION
   InputPara.SF_CreateStar_Scheme       = SF_CREATE_STAR_SCHEME;
   InputPara.SF_CreateStar_RSeed        = SF_CREATE_STAR_RSEED;
   InputPara.SF_CreateStar_DetRandom    = SF_CREATE_STAR_DET_RANDOM;
   InputPara.SF_CreateStar_MinLevel     = SF_CREATE_STAR_MIN_LEVEL;
   InputPara.SF_CreateStar_MinGasDens   = SF_CREATE_STAR_MIN_GAS_DENS;
   InputPara.SF_CreateStar_MassEff      = SF_CREATE_STAR_MASS_EFF;
   InputPara.SF_CreateStar_MinStarMass  = SF_CREATE_STAR_MIN_STAR_MASS;
   InputPara.SF_CreateStar_MaxStarMFrac = SF_CREATE_STAR_MAX_STAR_MFRAC;
#  endif

// feedback
#  ifdef FEEDBACK
   InputPara.FB_Level                = FB_LEVEL;
   InputPara.FB_RSeed                = FB_RSEED;
   InputPara.FB_SNe                  = FB_SNE;
   InputPara.FB_User                 = FB_USER;
#  endif

// cosmic ray
#  ifdef COSMIC_RAY
   InputPara.CR_Gamma                = GAMMA_CR;
#  ifdef CR_DIFFUSION
   InputPara.CR_Diffusion_ParaCoeff  = CR_DIFF_PARA;
   InputPara.CR_Diffusion_PerpCoeff  = CR_DIFF_PERP;
   InputPara.CR_Diffusion_MinB       = CR_DIFF_MIN_B;
#  endif
#  endif // #ifdef COSMIC_RAY

// initialization
   InputPara.Opt__Init               = OPT__INIT;
   InputPara.RestartLoadNRank        = RESTART_LOAD_NRANK;
   InputPara.Opt__RestartReset       = OPT__RESTART_RESET;
   InputPara.Opt__UM_IC_Level        = OPT__UM_IC_LEVEL;
   InputPara.Opt__UM_IC_NLevel       = OPT__UM_IC_NLEVEL;
   InputPara.Opt__UM_IC_NVar         = OPT__UM_IC_NVAR;
   InputPara.Opt__UM_IC_Format       = OPT__UM_IC_FORMAT;
   InputPara.Opt__UM_IC_Float8       = OPT__UM_IC_FLOAT8;
   InputPara.Opt__UM_IC_Downgrade    = OPT__UM_IC_DOWNGRADE;
   InputPara.Opt__UM_IC_Refine       = OPT__UM_IC_REFINE;
   InputPara.Opt__UM_IC_LoadNRank    = OPT__UM_IC_LOAD_NRANK;

   if ( OPT__INIT == INIT_BY_FILE  &&  OPT__UM_IC_NLEVEL > 1  &&  UM_IC_RefineRegion != NULL )
   {
      const int (*RefineRegion)[6] = ( int(*)[6] )UM_IC_RefineRegion;

      for (int t=0; t<NLEVEL-1; t++)
         if ( t < OPT__UM_IC_NLEVEL - 1 )
            for (int s=0; s<6; s++)
               InputPara.UM_IC_RefineRegion[t][s] = RefineRegion[t][s];
         else
            for (int s=0; s<6; s++)
               InputPara.UM_IC_RefineRegion[t][s] = -1;
   }
   else
   {
      for (int t=0; t<NLEVEL-1; t++)
      for (int s=0; s<6; s++)
         InputPara.UM_IC_RefineRegion[t][s] = -1;
   }

   InputPara.Opt__InitRestrict       = OPT__INIT_RESTRICT;
   InputPara.Opt__InitGridWithOMP    = OPT__INIT_GRID_WITH_OMP;
   InputPara.Opt__GPUID_Select       = OPT__GPUID_SELECT;
   InputPara.Init_Subsampling_NCell  = INIT_SUBSAMPLING_NCELL;
#  ifdef MHD
   InputPara.Opt__InitBFieldByVecPot = OPT__INIT_BFIELD_BYVECPOT;
#  endif
#  ifdef SUPPORT_FFTW
   InputPara.Opt__FFTW_Startup       = OPT__FFTW_STARTUP;
#  endif

// interpolation schemes
   InputPara.Opt__Int_Time               = OPT__INT_TIME;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Int_Prim               = OPT__INT_PRIM;
#  endif
#  if ( MODEL == ELBDM )
   InputPara.Opt__Int_Phase              = OPT__INT_PHASE;
   InputPara.Opt__Res_Phase              = OPT__RES_PHASE;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   InputPara.Opt__Hybrid_Match_Phase     = ELBDM_MATCH_PHASE;
#  endif
#  endif // ELBDM
   InputPara.Opt__Flu_IntScheme          = OPT__FLU_INT_SCHEME;
   InputPara.Opt__RefFlu_IntScheme       = OPT__REF_FLU_INT_SCHEME;
#  ifdef MHD
   InputPara.Opt__Mag_IntScheme          = OPT__MAG_INT_SCHEME;
   InputPara.Opt__RefMag_IntScheme       = OPT__REF_MAG_INT_SCHEME;
#  endif
#  ifdef GRAVITY
   InputPara.Opt__Pot_IntScheme          = OPT__POT_INT_SCHEME;
   InputPara.Opt__Rho_IntScheme          = OPT__RHO_INT_SCHEME;
   InputPara.Opt__Gra_IntScheme          = OPT__GRA_INT_SCHEME;
   InputPara.Opt__RefPot_IntScheme       = OPT__REF_POT_INT_SCHEME;
#  endif
   InputPara.IntMonoCoeff                = INT_MONO_COEFF;
#  ifdef MHD
   InputPara.IntMonoCoeffB               = INT_MONO_COEFF_B;
#  endif
   InputPara.Mono_MaxIter                = MONO_MAX_ITER;
   InputPara.IntOppSign0thOrder          = INT_OPP_SIGN_0TH_ORDER;
#  ifdef SUPPORT_SPECTRAL_INT
   InputPara.SpecInt_TablePath           = SPEC_INT_TABLE_PATH;
   InputPara.SpecInt_GhostBoundary       = SPEC_INT_GHOST_BOUNDARY;
#  if ( MODEL == ELBDM )
   InputPara.SpecInt_XY_Instead_DePha    = SPEC_INT_XY_INSTEAD_DEPHA;
   InputPara.SpecInt_VortexThreshold     = SPEC_INT_VORTEX_THRESHOLD;
#  endif
#  endif // #ifdef SUPPORT_SPECTRAL_INT

// data dump
   InputPara.Opt__Output_Total           = OPT__OUTPUT_TOTAL;
   InputPara.Opt__Output_Part            = OPT__OUTPUT_PART;
   InputPara.Opt__Output_User            = OPT__OUTPUT_USER;
#  ifdef PARTICLE
   InputPara.Opt__Output_Par_Mode        = OPT__OUTPUT_PAR_MODE;
   InputPara.Opt__Output_Par_Mesh        = OPT__OUTPUT_PAR_MESH;
#  endif
   InputPara.Opt__Output_BasePS          = OPT__OUTPUT_BASEPS;
   InputPara.Opt__Output_Base            = OPT__OUTPUT_BASE;
#  ifdef GRAVITY
   InputPara.Opt__Output_Pot             = OPT__OUTPUT_POT;
#  endif
#  ifdef PARTICLE
   InputPara.Opt__Output_ParDens         = OPT__OUTPUT_PAR_DENS;
#  endif
#  ifdef MHD
   InputPara.Opt__Output_CC_Mag          = OPT__OUTPUT_CC_MAG;
#  endif
#  if ( MODEL == HYDRO )
   InputPara.Opt__Output_Pres            = OPT__OUTPUT_PRES;
   InputPara.Opt__Output_Temp            = OPT__OUTPUT_TEMP;
   InputPara.Opt__Output_Entr            = OPT__OUTPUT_ENTR;
   InputPara.Opt__Output_Cs              = OPT__OUTPUT_CS;
   InputPara.Opt__Output_DivVel          = OPT__OUTPUT_DIVVEL;
   InputPara.Opt__Output_Mach            = OPT__OUTPUT_MACH;
#  ifdef MHD
   InputPara.Opt__Output_DivMag          = OPT__OUTPUT_DIVMAG;
#  endif
#  ifdef SRHD
   InputPara.Opt__Output_3Velocity       = OPT__OUTPUT_3VELOCITY;
   InputPara.Opt__Output_Lorentz         = OPT__OUTPUT_LORENTZ;
   InputPara.Opt__Output_Enthalpy        = OPT__OUTPUT_ENTHALPY;
#  endif
#  endif // #if ( MODEL == HYDRO )
   InputPara.Opt__Output_UserField       = OPT__OUTPUT_USER_FIELD;
   InputPara.Opt__Output_Mode            = OPT__OUTPUT_MODE;
   InputPara.Opt__Output_Restart         = OPT__OUTPUT_RESTART;
   InputPara.Opt__Output_Step            = OUTPUT_STEP;
   InputPara.Opt__Output_Dt              = OUTPUT_DT;
   InputPara.Opt__Output_Text_Format_Flt = OPT__OUTPUT_TEXT_FORMAT_FLT;
   InputPara.Opt__Output_Text_Length_Int = OPT__OUTPUT_TEXT_LENGTH_INT;
   InputPara.Output_PartX                = OUTPUT_PART_X;
   InputPara.Output_PartY                = OUTPUT_PART_Y;
   InputPara.Output_PartZ                = OUTPUT_PART_Z;
   InputPara.InitDumpID                  = INIT_DUMPID;

// libyt jupyter
#  if ( defined(SUPPORT_LIBYT) && defined(LIBYT_JUPYTER) )
   InputPara.Yt_JupyterUseConnectionFile = YT_JUPYTER_USE_CONNECTION_FILE;
#  endif

// miscellaneous
   InputPara.Opt__Verbose            = OPT__VERBOSE;
   InputPara.Opt__TimingBarrier      = OPT__TIMING_BARRIER;
   InputPara.Opt__TimingBalance      = OPT__TIMING_BALANCE;
   InputPara.Opt__TimingMPI          = OPT__TIMING_MPI;
   InputPara.Opt__RecordNote         = OPT__RECORD_NOTE;
   InputPara.Opt__RecordUnphy        = OPT__RECORD_UNPHY;
   InputPara.Opt__RecordMemory       = OPT__RECORD_MEMORY;
   InputPara.Opt__RecordPerformance  = OPT__RECORD_PERFORMANCE;
   InputPara.Opt__ManualControl      = OPT__MANUAL_CONTROL;
   InputPara.Opt__RecordCenter       = OPT__RECORD_CENTER;
   InputPara.COM_CenX                = COM_CEN_X;
   InputPara.COM_CenY                = COM_CEN_Y;
   InputPara.COM_CenZ                = COM_CEN_Z;
   InputPara.COM_MaxR                = COM_MAX_R;
   InputPara.COM_MinRho              = COM_MIN_RHO;
   InputPara.COM_TolErrR             = COM_TOLERR_R;
   InputPara.COM_MaxIter             = COM_MAX_ITER;
   InputPara.Opt__RecordUser         = OPT__RECORD_USER;
   InputPara.Opt__OptimizeAggressive = OPT__OPTIMIZE_AGGRESSIVE;
   InputPara.Opt__SortPatchByLBIdx   = OPT__SORT_PATCH_BY_LBIDX;

// simulation checks
   InputPara.Opt__Ck_Refine          = OPT__CK_REFINE;
   InputPara.Opt__Ck_ProperNesting   = OPT__CK_PROPER_NESTING;
   InputPara.Opt__Ck_Conservation    = OPT__CK_CONSERVATION;
   InputPara.AngMom_OriginX          = ANGMOM_ORIGIN_X;
   InputPara.AngMom_OriginY          = ANGMOM_ORIGIN_Y;
   InputPara.AngMom_OriginZ          = ANGMOM_ORIGIN_Z;
   InputPara.Opt__Ck_NormPassive     = OPT__CK_NORMALIZE_PASSIVE;
   InputPara.Opt__Ck_Restrict        = OPT__CK_RESTRICT;
   InputPara.Opt__Ck_Finite          = OPT__CK_FINITE;
   InputPara.Opt__Ck_PatchAllocate   = OPT__CK_PATCH_ALLOCATE;
   InputPara.Opt__Ck_FluxAllocate    = OPT__CK_FLUX_ALLOCATE;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Ck_Negative        = OPT__CK_NEGATIVE;
#  endif
   InputPara.Opt__Ck_MemFree         = OPT__CK_MEMFREE;
#  ifdef PARTICLE
   InputPara.Opt__Ck_Particle        = OPT__CK_PARTICLE;
#  endif
#  ifdef MHD
   InputPara.Opt__Ck_InterfaceB      = OPT__CK_INTERFACE_B;
   InputPara.Opt__Ck_DivergenceB     = OPT__CK_DIVERGENCE_B;
#  endif
   InputPara.Opt__Ck_InputFluid      = OPT__CK_INPUT_FLUID;

// flag tables
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      InputPara.FlagTable_Rho         [lv]    = FlagTable_Rho         [lv];
      InputPara.FlagTable_RhoGradient [lv]    = FlagTable_RhoGradient [lv];

      for (int t=0; t<5; t++)
      InputPara.FlagTable_Lohner      [lv][t] = FlagTable_Lohner      [lv][t];

      for (int t=0; t<3; t++)
      InputPara.FlagTable_Angular     [lv][t] = FlagTable_Angular     [lv][t];

      InputPara.FlagTable_Radial      [lv]    = FlagTable_Radial      [lv];

      InputPara.FlagTable_User        [lv].p   = malloc( OPT__FLAG_USER_NUM*sizeof(double) );
      InputPara.FlagTable_User        [lv].len = OPT__FLAG_USER_NUM;
      for (int t=0; t<OPT__FLAG_USER_NUM; t++)
      ( (double *) InputPara.FlagTable_User[lv].p )[t] = FlagTable_User[lv][t];

#     if   ( MODEL == HYDRO )
      InputPara.FlagTable_PresGradient[lv]    = FlagTable_PresGradient[lv];
      InputPara.FlagTable_Vorticity   [lv]    = FlagTable_Vorticity   [lv];
      InputPara.FlagTable_Jeans       [lv]    = FlagTable_Jeans       [lv];
#     ifdef MHD
      InputPara.FlagTable_Current     [lv]    = FlagTable_Current     [lv];
#     endif
#     ifdef SRHD
      InputPara.FlagTable_LrtzGradient[lv]    = FlagTable_LrtzGradient[lv];
#     endif
#     ifdef COSMIC_RAY
      InputPara.FlagTable_CRay        [lv]    = FlagTable_CRay        [lv];
#     endif

#     elif ( MODEL == ELBDM )
      for (int t=0; t<2; t++) {
      InputPara.FlagTable_EngyDensity [lv][t] = FlagTable_EngyDensity [lv][t];
      }
      for (int t=0; t<2; t++) {
      InputPara.FlagTable_Spectral    [lv][t] = FlagTable_Spectral    [lv][t];
      }
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      for (int t=0; t<4; t++) {
      InputPara.FlagTable_Interference[lv][t] = FlagTable_Interference[lv][t];
      }
#     endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
#     endif // #elif ( MODEL == ELBDM )

#     ifdef PARTICLE
      InputPara.FlagTable_NParPatch   [lv]    = FlagTable_NParPatch   [lv];
      InputPara.FlagTable_NParCell    [lv]    = FlagTable_NParCell    [lv];
      InputPara.FlagTable_ParMassCell [lv]    = FlagTable_ParMassCell [lv];
#     endif
   } // for (int lv=0; lv<NLEVEL-1; lv++)

   InputPara.UserDerField_Num = UserDerField_Num;
   for (int v=0; v<UserDerField_Num; v++)
   {
      InputPara.UserDerField_Label[v] = UserDerField_Label[v];
      InputPara.UserDerField_Unit [v] = UserDerField_Unit [v];
   }

} // FUNCTION : FillIn_InputPara



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_KeyInfo
// Description :  Create the HDF5 compound datatype for KeyInfo
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_KeyInfo( hid_t &H5_TypeID )
{

// create the array type
   const hsize_t H5_ArrDims_3Var         = 3;             // array size of [3]
   const hsize_t H5_ArrDims_NLv          = NLEVEL;        // array size of [NLEVEL]
   const hsize_t H5_ArrDims_NConRef      = 1+NCONREF_MAX; // array size of [1+NCONREF_MAX]

   const hid_t   H5_TypeID_Arr_3Double   = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_3Var    );
   const hid_t   H5_TypeID_Arr_3Int      = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_3Var    );
   const hid_t   H5_TypeID_Arr_NLvInt    = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_NLv     );
   const hid_t   H5_TypeID_Arr_NLvLong   = H5Tarray_create( H5T_NATIVE_LONG,   1, &H5_ArrDims_NLv     );
   const hid_t   H5_TypeID_Arr_NLvDouble = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_NLv     );
   const hid_t   H5_TypeID_Arr_NConRef   = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_NConRef );


// create the "variable-length string" datatype
   hid_t  H5_TypeID_VarStr;
   herr_t H5_Status;

   H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
   H5_Status        = H5Tset_size( H5_TypeID_VarStr, H5T_VARIABLE );


// get the compound type
   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(KeyInfo_t) );

   H5Tinsert( H5_TypeID, "FormatVersion",        HOFFSET(KeyInfo_t,FormatVersion       ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Model",                HOFFSET(KeyInfo_t,Model               ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Float8",               HOFFSET(KeyInfo_t,Float8              ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Gravity",              HOFFSET(KeyInfo_t,Gravity             ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Particle",             HOFFSET(KeyInfo_t,Particle            ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NLevel",               HOFFSET(KeyInfo_t,NLevel              ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NCompFluid",           HOFFSET(KeyInfo_t,NCompFluid          ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NCompPassive",         HOFFSET(KeyInfo_t,NCompPassive        ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "PatchSize",            HOFFSET(KeyInfo_t,PatchSize           ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "DumpID",               HOFFSET(KeyInfo_t,DumpID              ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NX0",                  HOFFSET(KeyInfo_t,NX0                 ), H5_TypeID_Arr_3Int      );
   H5Tinsert( H5_TypeID, "BoxScale",             HOFFSET(KeyInfo_t,BoxScale            ), H5_TypeID_Arr_3Int      );
   H5Tinsert( H5_TypeID, "NPatch",               HOFFSET(KeyInfo_t,NPatch              ), H5_TypeID_Arr_NLvInt    );
   H5Tinsert( H5_TypeID, "CellScale",            HOFFSET(KeyInfo_t,CellScale           ), H5_TypeID_Arr_NLvInt    );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Magnetohydrodynamics", HOFFSET(KeyInfo_t,Magnetohydrodynamics), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "SRHydrodynamics",      HOFFSET(KeyInfo_t,SRHydrodynamics     ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "CosmicRay",            HOFFSET(KeyInfo_t,CosmicRay           ), H5T_NATIVE_INT          );
#  endif

   H5Tinsert( H5_TypeID, "Step",                 HOFFSET(KeyInfo_t,Step                ), H5T_NATIVE_LONG         );
   H5Tinsert( H5_TypeID, "AdvanceCounter",       HOFFSET(KeyInfo_t,AdvanceCounter      ), H5_TypeID_Arr_NLvLong   );
   H5Tinsert( H5_TypeID, "NFieldStored",         HOFFSET(KeyInfo_t,NFieldStored        ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NMagStored",           HOFFSET(KeyInfo_t,NMagStored          ), H5T_NATIVE_INT          );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Par_NPar",             HOFFSET(KeyInfo_t,Par_NPar            ), H5T_NATIVE_LONG         );
   H5Tinsert( H5_TypeID, "Par_NAttFltStored",    HOFFSET(KeyInfo_t,Par_NAttFltStored   ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Par_NAttIntStored",    HOFFSET(KeyInfo_t,Par_NAttIntStored   ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Float8_Par",           HOFFSET(KeyInfo_t,Float8_Par          ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Int8_Par",             HOFFSET(KeyInfo_t,Int8_Par            ), H5T_NATIVE_INT          );
#  endif

#  ifdef COSMIC_RAY
   H5Tinsert( H5_TypeID, "CR_Diffusion",         HOFFSET(KeyInfo_t,CR_Diffusion        ), H5T_NATIVE_INT          );
#  endif

   H5Tinsert( H5_TypeID, "BoxSize",              HOFFSET(KeyInfo_t,BoxSize             ), H5_TypeID_Arr_3Double   );
   H5Tinsert( H5_TypeID, "Time",                 HOFFSET(KeyInfo_t,Time                ), H5_TypeID_Arr_NLvDouble );
   H5Tinsert( H5_TypeID, "CellSize",             HOFFSET(KeyInfo_t,CellSize            ), H5_TypeID_Arr_NLvDouble );
   H5Tinsert( H5_TypeID, "dTime_AllLv",          HOFFSET(KeyInfo_t,dTime_AllLv         ), H5_TypeID_Arr_NLvDouble );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "AveDens_Init",         HOFFSET(KeyInfo_t,AveDens_Init        ), H5T_NATIVE_DOUBLE       );
#  endif
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "UseWaveScheme",        HOFFSET(KeyInfo_t,UseWaveScheme       ), H5_TypeID_Arr_NLvInt    );
#  endif

   H5Tinsert( H5_TypeID, "CodeVersion",          HOFFSET(KeyInfo_t,CodeVersion         ), H5_TypeID_VarStr        );
   H5Tinsert( H5_TypeID, "DumpWallTime",         HOFFSET(KeyInfo_t,DumpWallTime        ), H5_TypeID_VarStr        );
   H5Tinsert( H5_TypeID, "GitBranch",            HOFFSET(KeyInfo_t,GitBranch           ), H5_TypeID_VarStr        );
   H5Tinsert( H5_TypeID, "GitCommit",            HOFFSET(KeyInfo_t,GitCommit           ), H5_TypeID_VarStr        );
   H5Tinsert( H5_TypeID, "UniqueDataID",         HOFFSET(KeyInfo_t,UniqueDataID        ), H5T_NATIVE_LONG         );

   H5Tinsert( H5_TypeID, "ConRef",               HOFFSET(KeyInfo_t,ConRef              ), H5_TypeID_Arr_NConRef   );

// free memory
   H5_Status = H5Tclose( H5_TypeID_Arr_3Double   );
   H5_Status = H5Tclose( H5_TypeID_Arr_3Int      );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvInt    );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvLong   );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvDouble );
   H5_Status = H5Tclose( H5_TypeID_Arr_NConRef   );
   H5_Status = H5Tclose( H5_TypeID_VarStr        );

} // FUNCTION : GetCompound_KeyInfo



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_Makefile
// Description :  Create the HDF5 compound datatype for Makefile
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_Makefile( hid_t &H5_TypeID )
{

   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(Makefile_t) );

   H5Tinsert( H5_TypeID, "Model",                  HOFFSET(Makefile_t,Model                  ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Gravity",                HOFFSET(Makefile_t,Gravity                ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Comoving",               HOFFSET(Makefile_t,Comoving               ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Particle",               HOFFSET(Makefile_t,Particle               ), H5T_NATIVE_INT );

   H5Tinsert( H5_TypeID, "UseGPU",                 HOFFSET(Makefile_t,UseGPU                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "GAMER_Debug",            HOFFSET(Makefile_t,GAMER_Debug            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "BitwiseReproducibility", HOFFSET(Makefile_t,BitwiseReproducibility ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Timing",                 HOFFSET(Makefile_t,Timing                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "TimingSolver",           HOFFSET(Makefile_t,TimingSolver           ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Float8",                 HOFFSET(Makefile_t,Float8                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Serial",                 HOFFSET(Makefile_t,Serial                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "LoadBalance",            HOFFSET(Makefile_t,LoadBalance            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "OverlapMPI",             HOFFSET(Makefile_t,OverlapMPI             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "OpenMP",                 HOFFSET(Makefile_t,OpenMP                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "GPU_Arch",               HOFFSET(Makefile_t,GPU_Arch               ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Laohu",                  HOFFSET(Makefile_t,Laohu                  ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportHDF5",            HOFFSET(Makefile_t,SupportHDF5            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportGSL",             HOFFSET(Makefile_t,SupportGSL             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportSpectralInt",     HOFFSET(Makefile_t,SupportSpectralInt     ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportFFTW",            HOFFSET(Makefile_t,SupportFFTW            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportLibYT",           HOFFSET(Makefile_t,SupportLibYT           ), H5T_NATIVE_INT );
#  ifdef SUPPORT_LIBYT
   H5Tinsert( H5_TypeID, "LibYTUsePatchGroup",     HOFFSET(Makefile_t,LibYTUsePatchGroup     ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "LibYTInteractive",       HOFFSET(Makefile_t,LibYTInteractive       ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "LibYTReload",            HOFFSET(Makefile_t,LibYTReload            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "LibYTJupyter",           HOFFSET(Makefile_t,LibYTJupyter           ), H5T_NATIVE_INT );
#  endif
   H5Tinsert( H5_TypeID, "SupportGrackle",         HOFFSET(Makefile_t,SupportGrackle         ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "RandomNumber",           HOFFSET(Makefile_t,RandomNumber           ), H5T_NATIVE_INT );

   H5Tinsert( H5_TypeID, "NLevel",                 HOFFSET(Makefile_t,NLevel                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "MaxPatch",               HOFFSET(Makefile_t,MaxPatch               ), H5T_NATIVE_INT );

#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "PotScheme",              HOFFSET(Makefile_t,PotScheme              ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "StorePotGhost",          HOFFSET(Makefile_t,StorePotGhost          ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "UnsplitGravity",         HOFFSET(Makefile_t,UnsplitGravity         ), H5T_NATIVE_INT );
#  endif

#  if   ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "FluScheme",              HOFFSET(Makefile_t,FluScheme              ), H5T_NATIVE_INT );
#  ifdef LR_SCHEME
   H5Tinsert( H5_TypeID, "LRScheme",               HOFFSET(Makefile_t,LRScheme               ), H5T_NATIVE_INT );
#  endif
#  ifdef RSOLVER
   H5Tinsert( H5_TypeID, "RSolver",                HOFFSET(Makefile_t,RSolver                ), H5T_NATIVE_INT );
#  endif
   H5Tinsert( H5_TypeID, "DualEnergy",             HOFFSET(Makefile_t,DualEnergy             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Magnetohydrodynamics",   HOFFSET(Makefile_t,Magnetohydrodynamics   ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SRHydrodynamics",        HOFFSET(Makefile_t,SRHydrodynamics        ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "CosmicRay",              HOFFSET(Makefile_t,CosmicRay              ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "EoS",                    HOFFSET(Makefile_t,EoS                    ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "BarotropicEoS",          HOFFSET(Makefile_t,BarotropicEoS          ), H5T_NATIVE_INT );

#  elif ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "ELBDMScheme",            HOFFSET(Makefile_t,ELBDMScheme            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "WaveScheme",             HOFFSET(Makefile_t,WaveScheme             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "ConserveMass",           HOFFSET(Makefile_t,ConserveMass           ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Laplacian4th",           HOFFSET(Makefile_t,Laplacian4th           ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SelfInteraction4",       HOFFSET(Makefile_t,SelfInteraction4       ), H5T_NATIVE_INT );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "MassiveParticles",       HOFFSET(Makefile_t,MassiveParticles       ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Tracer",                 HOFFSET(Makefile_t,Tracer                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "StoreParAcc",            HOFFSET(Makefile_t,StoreParAcc            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "StarFormation",          HOFFSET(Makefile_t,StarFormation          ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Feedback",               HOFFSET(Makefile_t,Feedback               ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Par_NAttFltUser",        HOFFSET(Makefile_t,Par_NAttFltUser        ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Par_NAttIntUser",        HOFFSET(Makefile_t,Par_NAttIntUser        ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Float8_Par",             HOFFSET(Makefile_t,Float8_Par             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Int8_Par",               HOFFSET(Makefile_t,Int8_Par               ), H5T_NATIVE_INT );
#  endif

#  ifdef COSMIC_RAY
   H5Tinsert( H5_TypeID, "CR_Diffusion",           HOFFSET(Makefile_t,CR_Diffusion           ), H5T_NATIVE_INT );
#  endif

} // FUNCTION : GetCompound_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_SymConst
// Description :  Create the HDF5 compound datatype for SymConst
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_SymConst( hid_t &H5_TypeID )
{

   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(SymConst_t) );

   H5Tinsert( H5_TypeID, "NCompFluid",           HOFFSET(SymConst_t,NCompFluid          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "NCompPassive",         HOFFSET(SymConst_t,NCompPassive        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "PatchSize",            HOFFSET(SymConst_t,PatchSize           ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NIn",              HOFFSET(SymConst_t,Flu_NIn             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NOut",             HOFFSET(SymConst_t,Flu_NOut            ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NIn_T",            HOFFSET(SymConst_t,Flu_NIn_T           ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NIn_S",            HOFFSET(SymConst_t,Flu_NIn_S           ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NOut_S",           HOFFSET(SymConst_t,Flu_NOut_S          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "NFluxFluid",           HOFFSET(SymConst_t,NFluxFluid          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "NFluxPassive",         HOFFSET(SymConst_t,NFluxPassive        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_GhostSize",        HOFFSET(SymConst_t,Flu_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_Nxt",              HOFFSET(SymConst_t,Flu_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Debug_HDF5",           HOFFSET(SymConst_t,Debug_HDF5          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SibOffsetNonperiodic", HOFFSET(SymConst_t,SibOffsetNonperiodic), H5T_NATIVE_INT    );
#  ifdef LOAD_BALANCE
   H5Tinsert( H5_TypeID, "SonOffsetLB",          HOFFSET(SymConst_t,SonOffsetLB         ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "TinyNumber",           HOFFSET(SymConst_t,TinyNumber          ), H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_TypeID, "HugeNumber",           HOFFSET(SymConst_t,HugeNumber          ), H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_TypeID, "MaxError",             HOFFSET(SymConst_t,MaxError            ), H5T_NATIVE_DOUBLE );

#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Gra_NIn",              HOFFSET(SymConst_t,Gra_NIn             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Pot_GhostSize",        HOFFSET(SymConst_t,Pot_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Gra_GhostSize",        HOFFSET(SymConst_t,Gra_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Rho_GhostSize",        HOFFSET(SymConst_t,Rho_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Pot_Nxt",              HOFFSET(SymConst_t,Pot_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Gra_Nxt",              HOFFSET(SymConst_t,Gra_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Rho_Nxt",              HOFFSET(SymConst_t,Rho_Nxt             ), H5T_NATIVE_INT    );
#  ifdef UNSPLIT_GRAVITY
   H5Tinsert( H5_TypeID, "USG_GhostSizeF",       HOFFSET(SymConst_t,USG_GhostSizeF      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "USG_GhostSizeG",       HOFFSET(SymConst_t,USG_GhostSizeG      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "USG_NxtF",             HOFFSET(SymConst_t,USG_NxtF            ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "USG_NxtG",             HOFFSET(SymConst_t,USG_NxtG            ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "ExtPot_BlockSize",     HOFFSET(SymConst_t,ExtPot_BlockSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Gra_BlockSize",        HOFFSET(SymConst_t,Gra_BlockSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ExtPotNAuxMax",        HOFFSET(SymConst_t,ExtPotNAuxMax       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ExtAccNAuxMax",        HOFFSET(SymConst_t,ExtAccNAuxMax       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ExtPotNGeneMax",       HOFFSET(SymConst_t,ExtPotNGeneMax      ), H5T_NATIVE_INT    );
#  if   ( POT_SCHEME == SOR )
   H5Tinsert( H5_TypeID, "Pot_BlockSize_z",      HOFFSET(SymConst_t,Pot_BlockSize_z     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_RhoShared",        HOFFSET(SymConst_t,SOR_RhoShared       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_CPotShared",       HOFFSET(SymConst_t,SOR_CPotShared      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_UseShuffle",       HOFFSET(SymConst_t,SOR_UseShuffle      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_UsePadding",       HOFFSET(SymConst_t,SOR_UsePadding      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_ModReduction",     HOFFSET(SymConst_t,SOR_ModReduction    ), H5T_NATIVE_INT    );
#  elif ( POT_SCHEME == MG  )
   H5Tinsert( H5_TypeID, "Pot_BlockSize_x",      HOFFSET(SymConst_t,Pot_BlockSize_x     ), H5T_NATIVE_INT    );
#  endif
#  endif // #ifdef GRAVITY

#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Par_NAttFltStored",    HOFFSET(SymConst_t,Par_NAttFltStored   ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Par_NAttIntStored",    HOFFSET(SymConst_t,Par_NAttIntStored   ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Par_NType",            HOFFSET(SymConst_t,Par_NType           ), H5T_NATIVE_INT    );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "RhoExt_GhostSize",     HOFFSET(SymConst_t,RhoExt_GhostSize    ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "Debug_Particle",       HOFFSET(SymConst_t,Debug_Particle      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ParList_GrowthFactor", HOFFSET(SymConst_t,ParList_GrowthFactor), H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_TypeID, "ParList_ReduceFactor", HOFFSET(SymConst_t,ParList_ReduceFactor), H5T_NATIVE_DOUBLE );
#  endif // #ifdef PARTICLE

   H5Tinsert( H5_TypeID, "BitRep_Flux",          HOFFSET(SymConst_t,BitRep_Flux         ), H5T_NATIVE_INT    );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "BitRep_Electric",      HOFFSET(SymConst_t,BitRep_Electric     ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "InterpMask",           HOFFSET(SymConst_t,InterpMask          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "FB_SepFluOut",         HOFFSET(SymConst_t,FB_SepFluOut        ), H5T_NATIVE_INT    );

#  if   ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Flu_BlockSize_x",      HOFFSET(SymConst_t,Flu_BlockSize_x     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_BlockSize_y",      HOFFSET(SymConst_t,Flu_BlockSize_y     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "CheckUnphyInFluid",    HOFFSET(SymConst_t,CheckUnphyInFluid   ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "CharReconstruction",   HOFFSET(SymConst_t,CharReconstruction  ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "LR_Eint",              HOFFSET(SymConst_t,LR_Eint             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "CheckIntermediate",    HOFFSET(SymConst_t,CheckIntermediate   ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "RSolverRescue",        HOFFSET(SymConst_t,RSolverRescue       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLL_NoRefState",       HOFFSET(SymConst_t,HLL_NoRefState      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLL_IncludeAllWaves",  HOFFSET(SymConst_t,HLL_IncludeAllWaves ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLLC_WaveSpeed",       HOFFSET(SymConst_t,HLLC_WaveSpeed      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLLE_WaveSpeed",       HOFFSET(SymConst_t,HLLE_WaveSpeed      ), H5T_NATIVE_INT    );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "HLLD_WaveSpeed",       HOFFSET(SymConst_t,HLLD_WaveSpeed      ), H5T_NATIVE_INT    );
#  endif
#  ifdef N_FC_VAR
   H5Tinsert( H5_TypeID, "N_FC_Var",             HOFFSET(SymConst_t,N_FC_Var            ), H5T_NATIVE_INT    );
#  endif
#  ifdef N_SLOPE_PPM
   H5Tinsert( H5_TypeID, "N_Slope_PPM",          HOFFSET(SymConst_t,N_Slope_PPM         ), H5T_NATIVE_INT    );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "EulerY",               HOFFSET(SymConst_t,EulerY              ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "MHM_CheckPredict",     HOFFSET(SymConst_t,MHM_CheckPredict    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "EoSNAuxMax",           HOFFSET(SymConst_t,EoSNAuxMax          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "EoSNTableMax",         HOFFSET(SymConst_t,EoSNTableMax        ), H5T_NATIVE_INT    );

#  elif  ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Flu_BlockSize_x",      HOFFSET(SymConst_t,Flu_BlockSize_x     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_BlockSize_y",      HOFFSET(SymConst_t,Flu_BlockSize_y     ), H5T_NATIVE_INT    );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "Flu_HJ_BlockSize_y",   HOFFSET(SymConst_t,Flu_HJ_BlockSize_y  ), H5T_NATIVE_INT    );
#  endif

#  if ( WAVE_SCHEME == WAVE_GRAMFE )
   H5Tinsert( H5_TypeID, "GramFEScheme",         HOFFSET(SymConst_t,GramFEScheme        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "GramFEGamma",          HOFFSET(SymConst_t,GramFEGamma         ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "GramFEG",              HOFFSET(SymConst_t,GramFEG             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "GramFENDelta",         HOFFSET(SymConst_t,GramFENDelta        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "GramFEOrder",          HOFFSET(SymConst_t,GramFEOrder         ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "GramFEND",             HOFFSET(SymConst_t,GramFEND            ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "GramFEFluNxt",         HOFFSET(SymConst_t,GramFEFluNxt        ), H5T_NATIVE_INT    );
#  endif

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   H5Tinsert( H5_TypeID, "dt_Flu_BlockSize",     HOFFSET(SymConst_t,dt_Flu_BlockSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "dt_Flu_UseShuffle",    HOFFSET(SymConst_t,dt_Flu_UseShuffle   ), H5T_NATIVE_INT    );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "dt_Gra_BlockSize",     HOFFSET(SymConst_t,dt_Gra_BlockSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "dt_Gra_UseShuffle",    HOFFSET(SymConst_t,dt_Gra_UseShuffle   ), H5T_NATIVE_INT    );
#  endif

   H5Tinsert( H5_TypeID, "Src_BlockSize",        HOFFSET(SymConst_t,Src_BlockSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_GhostSize",        HOFFSET(SymConst_t,Src_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_Nxt",              HOFFSET(SymConst_t,Src_Nxt             ), H5T_NATIVE_INT    );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Src_NAuxDlep",         HOFFSET(SymConst_t,Src_NAuxDlep        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_DlepProfNVar",     HOFFSET(SymConst_t,Src_DlepProfNVar    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_DlepProfNBinMax",  HOFFSET(SymConst_t,Src_DlepProfNBinMax ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "Src_NAuxUser",         HOFFSET(SymConst_t,Src_NAuxUser        ), H5T_NATIVE_INT    );

   H5Tinsert( H5_TypeID, "Der_GhostSize",        HOFFSET(SymConst_t,Der_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Der_Nxt",              HOFFSET(SymConst_t,Der_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Der_NOut_Max",         HOFFSET(SymConst_t,Der_NOut_Max        ), H5T_NATIVE_INT    );

#  ifdef FEEDBACK
   H5Tinsert( H5_TypeID, "FB_GhostSize",         HOFFSET(SymConst_t,FB_GhostSize        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "FB_Nxt",               HOFFSET(SymConst_t,FB_Nxt              ), H5T_NATIVE_INT    );
#  endif

   H5Tinsert( H5_TypeID, "NFieldStoredMax",      HOFFSET(SymConst_t,NFieldStoredMax     ), H5T_NATIVE_INT    );

   H5Tinsert( H5_TypeID, "NConRefMax",           HOFFSET(SymConst_t,NConRefMax          ), H5T_NATIVE_INT    );

} // FUNCTION : GetCompound_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_InputPara
// Description :  Create the HDF5 compound datatype for InputPara
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID    : HDF5 type ID for storing the compound datatype
//                NFieldStored : Number of grid fields to be stored on disk
//-------------------------------------------------------------------------------------------------------
void GetCompound_InputPara( hid_t &H5_TypeID, const int NFieldStored )
{

// create the array type
   const hsize_t H5_ArrDims_3Var              = 3;                    // array size of [3]
   const hsize_t H5_ArrDims_6Var              = 6;                    // array size of [6]
#  if ( NCOMP_PASSIVE > 0 )
   const hsize_t H5_ArrDims_NPassive          = NCOMP_PASSIVE;        // array size of [NCOMP_PASSIVE]
#  endif
#  if ( NLEVEL > 1 )
   const hsize_t H5_ArrDims_NLvM1             = NLEVEL-1;             // array size of [NLEVEL-1]
   const hsize_t H5_ArrDims_NLvM1_2[2]        = { NLEVEL-1, 2 };      // array size of [NLEVEL-1][2]
   const hsize_t H5_ArrDims_NLvM1_3[2]        = { NLEVEL-1, 3 };      // array size of [NLEVEL-1][3]
   const hsize_t H5_ArrDims_NLvM1_4[2]        = { NLEVEL-1, 4 };      // array size of [NLEVEL-1][4]
   const hsize_t H5_ArrDims_NLvM1_5[2]        = { NLEVEL-1, 5 };      // array size of [NLEVEL-1][5]
   const hsize_t H5_ArrDims_NLvM1_6[2]        = { NLEVEL-1, 6 };      // array size of [NLEVEL-1][6]
#  endif

   const hid_t   H5_TypeID_Arr_3Int           = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_3Var      );
   const hid_t   H5_TypeID_Arr_6Int           = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_6Var      );
#  if ( NCOMP_PASSIVE > 0 )
   const hid_t   H5_TypeID_Arr_NPassive       = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_NPassive  );
#  endif
#  if ( NLEVEL > 1 )
   const hid_t   H5_TypeID_Arr_NLvM1Int       = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_NLvM1     );
   const hid_t   H5_TypeID_Arr_NLvM1_6Int     = H5Tarray_create( H5T_NATIVE_INT,    2,  H5_ArrDims_NLvM1_6   );
   const hid_t   H5_TypeID_Arr_NLvM1Double    = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_NLvM1     );
   const hid_t   H5_TypeID_Arr_NLvM1_2Double  = H5Tarray_create( H5T_NATIVE_DOUBLE, 2,  H5_ArrDims_NLvM1_2   );
   const hid_t   H5_TypeID_Arr_NLvM1_3Double  = H5Tarray_create( H5T_NATIVE_DOUBLE, 2,  H5_ArrDims_NLvM1_3   );
   const hid_t   H5_TypeID_Arr_NLvM1_4Double  = H5Tarray_create( H5T_NATIVE_DOUBLE, 2,  H5_ArrDims_NLvM1_4   );
   const hid_t   H5_TypeID_Arr_NLvM1_5Double  = H5Tarray_create( H5T_NATIVE_DOUBLE, 2,  H5_ArrDims_NLvM1_5   );
   const hid_t   H5_TypeID_Arr_NLvM1_VLDouble = H5Tvlen_create ( H5T_NATIVE_DOUBLE );
#  endif
   const hid_t   H5_TypeID_Arr_3Double        = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_3Var      );


// create the "variable-length string" datatype
   hid_t  H5_TypeID_VarStr;
   herr_t H5_Status;

   H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
   H5_Status        = H5Tset_size( H5_TypeID_VarStr, H5T_VARIABLE );


// get the size of a single pointer, which is used for storing the array of variable-length strings
// --> FieldLabel[], MagLabel[], ParAttFltLabel[], ParAttIntLabel[]
   const int PtrSize     = sizeof( char* );
   const int PtrSize_hvl = sizeof( hvl_t );
   char Key[MAX_STRING];


// get the compound type
   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(InputPara_t) );

// simulation scale
   H5Tinsert( H5_TypeID, "BoxSize",                 HOFFSET(InputPara_t,BoxSize                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "NX0_Tot",                 HOFFSET(InputPara_t,NX0_Tot                ), H5_TypeID_Arr_3Int );
   H5Tinsert( H5_TypeID, "MPI_NRank",               HOFFSET(InputPara_t,MPI_NRank              ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "MPI_NRank_X",             HOFFSET(InputPara_t,MPI_NRank_X            ), H5_TypeID_Arr_3Int );
   H5Tinsert( H5_TypeID, "OMP_NThread",             HOFFSET(InputPara_t,OMP_NThread            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "EndT",                    HOFFSET(InputPara_t,EndT                   ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "EndStep",                 HOFFSET(InputPara_t,EndStep                ), H5T_NATIVE_LONG    );

// test problems
   H5Tinsert( H5_TypeID, "TestProb_ID",             HOFFSET(InputPara_t,TestProb_ID            ), H5T_NATIVE_INT     );

// code units
   H5Tinsert( H5_TypeID, "Opt__Unit",               HOFFSET(InputPara_t,Opt__Unit              ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Unit_L",                  HOFFSET(InputPara_t,Unit_L                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_M",                  HOFFSET(InputPara_t,Unit_M                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_T",                  HOFFSET(InputPara_t,Unit_T                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_V",                  HOFFSET(InputPara_t,Unit_V                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_D",                  HOFFSET(InputPara_t,Unit_D                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_E",                  HOFFSET(InputPara_t,Unit_E                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_P",                  HOFFSET(InputPara_t,Unit_P                 ), H5T_NATIVE_DOUBLE  );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Unit_B",                  HOFFSET(InputPara_t,Unit_B                 ), H5T_NATIVE_DOUBLE  );
#  endif

// boundary condition
   H5Tinsert( H5_TypeID, "Opt__BC_Flu",             HOFFSET(InputPara_t,Opt__BC_Flu            ), H5_TypeID_Arr_6Int );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Opt__BC_Pot",             HOFFSET(InputPara_t,Opt__BC_Pot            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "GFunc_Coeff0",            HOFFSET(InputPara_t,GFunc_Coeff0           ), H5T_NATIVE_DOUBLE  );
#  endif

// particle
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Par_Init",                HOFFSET(InputPara_t,Par_Init               ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ICFormat",            HOFFSET(InputPara_t,Par_ICFormat           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ICMass",              HOFFSET(InputPara_t,Par_ICMass             ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Par_ICType",              HOFFSET(InputPara_t,Par_ICType             ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ICFloat8",            HOFFSET(InputPara_t,Par_ICFloat8           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ICInt8",              HOFFSET(InputPara_t,Par_ICInt8             ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_Interp",              HOFFSET(InputPara_t,Par_Interp             ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_InterpTracer",        HOFFSET(InputPara_t,Par_InterpTracer       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_Integ",               HOFFSET(InputPara_t,Par_Integ              ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_IntegTracer",         HOFFSET(InputPara_t,Par_IntegTracer        ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ImproveAcc",          HOFFSET(InputPara_t,Par_ImproveAcc         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_PredictPos",          HOFFSET(InputPara_t,Par_PredictPos         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_TracerVelCorr",       HOFFSET(InputPara_t,Par_TracerVelCorr      ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_RemoveCell",          HOFFSET(InputPara_t,Par_RemoveCell         ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Opt__FreezePar",          HOFFSET(InputPara_t,Opt__FreezePar         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_GhostSize",           HOFFSET(InputPara_t,Par_GhostSize          ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_GhostSizeTracer",     HOFFSET(InputPara_t,Par_GhostSizeTracer    ), H5T_NATIVE_INT     );

// store the name of all particle attributes
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
   {
//    key for each particle attribute
      sprintf( Key, "ParAttFltLabel%02d", v );

//    assuming the offset between successive ParAttFltLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,ParAttFltLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }

   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
   {
//    key for each particle attribute
      sprintf( Key, "ParAttIntLabel%02d", v );

//    assuming the offset between successive ParAttIntLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,ParAttIntLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }
#  endif

// cosmology
#  ifdef COMOVING
   H5Tinsert( H5_TypeID, "A_Init",                  HOFFSET(InputPara_t,A_Init                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "OmegaM0",                 HOFFSET(InputPara_t,OmegaM0                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Hubble0",                 HOFFSET(InputPara_t,Hubble0                ), H5T_NATIVE_DOUBLE  );
#  endif

// time-step determination
   H5Tinsert( H5_TypeID, "Dt__Max",                 HOFFSET(InputPara_t,Dt__Max                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__Fluid",               HOFFSET(InputPara_t,Dt__Fluid              ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__FluidInit",           HOFFSET(InputPara_t,Dt__FluidInit          ), H5T_NATIVE_DOUBLE  );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Dt__Gravity",             HOFFSET(InputPara_t,Dt__Gravity            ), H5T_NATIVE_DOUBLE  );
#  endif
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Dt__Phase",               HOFFSET(InputPara_t,Dt__Phase              ), H5T_NATIVE_DOUBLE  );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "Dt__HybridCFL",           HOFFSET(InputPara_t,Dt__HybridCFL          ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__HybridCFLInit",       HOFFSET(InputPara_t,Dt__HybridCFLInit      ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__HybridVelocity",      HOFFSET(InputPara_t,Dt__HybridVelocity     ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__HybridVelocityInit",  HOFFSET(InputPara_t,Dt__HybridVelocityInit ), H5T_NATIVE_DOUBLE  );
#  endif
#  endif // ELBDM
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Dt__ParVel",              HOFFSET(InputPara_t,Dt__ParVel             ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__ParVelMax",           HOFFSET(InputPara_t,Dt__ParVelMax          ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__ParAcc",              HOFFSET(InputPara_t,Dt__ParAcc             ), H5T_NATIVE_DOUBLE  );
#  endif
#  ifdef SRHD
   H5Tinsert( H5_TypeID, "Dt__SpeedOfLight",        HOFFSET(InputPara_t,Dt__SpeedOfLight       ), H5T_NATIVE_INT     );
#  endif
#  ifdef CR_DIFFUSION
   H5Tinsert( H5_TypeID, "Dt__CR_Diffusion",        HOFFSET(InputPara_t,Dt__CR_Diffusion       ), H5T_NATIVE_DOUBLE  );
#  endif
#  ifdef COMOVING
   H5Tinsert( H5_TypeID, "Dt__MaxDeltaA",           HOFFSET(InputPara_t,Dt__MaxDeltaA          ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "Dt__SyncParentLv",        HOFFSET(InputPara_t,Dt__SyncParentLv       ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__SyncChildrenLv",      HOFFSET(InputPara_t,Dt__SyncChildrenLv     ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Opt__DtUser",             HOFFSET(InputPara_t,Opt__DtUser            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__DtLevel",            HOFFSET(InputPara_t,Opt__DtLevel           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__RecordDt",           HOFFSET(InputPara_t,Opt__RecordDt          ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "AutoReduceDt",            HOFFSET(InputPara_t,AutoReduceDt           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "AutoReduceDtFactor",      HOFFSET(InputPara_t,AutoReduceDtFactor     ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "AutoReduceDtFactorMin",   HOFFSET(InputPara_t,AutoReduceDtFactorMin  ), H5T_NATIVE_DOUBLE  );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "AutoReduceMinModFactor",  HOFFSET(InputPara_t,AutoReduceMinModFactor ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "AutoReduceMinModMin",     HOFFSET(InputPara_t,AutoReduceMinModMin    ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "AutoReduceIntMonoFactor", HOFFSET(InputPara_t,AutoReduceIntMonoFactor), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "AutoReduceIntMonoMin",    HOFFSET(InputPara_t,AutoReduceIntMonoMin   ), H5T_NATIVE_DOUBLE  );


// domain refinement
   H5Tinsert( H5_TypeID, "RegridCount",             HOFFSET(InputPara_t,RegridCount            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "RefineNLevel",            HOFFSET(InputPara_t,RefineNLevel           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagBufferSize",          HOFFSET(InputPara_t,FlagBufferSize         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagBufferSizeMaxM1Lv",   HOFFSET(InputPara_t,FlagBufferSizeMaxM1Lv  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagBufferSizeMaxM2Lv",   HOFFSET(InputPara_t,FlagBufferSizeMaxM2Lv  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "MaxLevel",                HOFFSET(InputPara_t,MaxLevel               ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Rho",           HOFFSET(InputPara_t,Opt__Flag_Rho          ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_RhoGradient",   HOFFSET(InputPara_t,Opt__Flag_RhoGradient  ), H5T_NATIVE_INT     );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Flag_PresGradient",  HOFFSET(InputPara_t,Opt__Flag_PresGradient ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Vorticity",     HOFFSET(InputPara_t,Opt__Flag_Vorticity    ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Jeans",         HOFFSET(InputPara_t,Opt__Flag_Jeans        ), H5T_NATIVE_INT     );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Flag_Current",       HOFFSET(InputPara_t,Opt__Flag_Current      ), H5T_NATIVE_INT     );
#  endif
#  ifdef SRHD
   H5Tinsert( H5_TypeID, "Opt__Flag_LrtzGradient",  HOFFSET(InputPara_t,Opt__Flag_LrtzGradient ), H5T_NATIVE_INT     );
#  endif
#  ifdef COSMIC_RAY
   H5Tinsert( H5_TypeID, "Opt__Flag_CRay",          HOFFSET(InputPara_t,Opt__Flag_CRay         ), H5T_NATIVE_INT     );
#  endif
#  endif
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Opt__Flag_EngyDensity",   HOFFSET(InputPara_t,Opt__Flag_EngyDensity  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Spectral",      HOFFSET(InputPara_t,Opt__Flag_Spectral     ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Spectral_N",    HOFFSET(InputPara_t,Opt__Flag_Spectral_N   ), H5T_NATIVE_INT     );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "Opt__Flag_Interference",  HOFFSET(InputPara_t,Opt__Flag_Interference ), H5T_NATIVE_INT     );
#  endif
#  endif // ELBDM
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerDens",    HOFFSET(InputPara_t,Opt__Flag_LohnerDens   ), H5T_NATIVE_INT     );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerEngy",    HOFFSET(InputPara_t,Opt__Flag_LohnerEngy   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerPres",    HOFFSET(InputPara_t,Opt__Flag_LohnerPres   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerTemp",    HOFFSET(InputPara_t,Opt__Flag_LohnerTemp   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerEntr",    HOFFSET(InputPara_t,Opt__Flag_LohnerEntr   ), H5T_NATIVE_INT     );
#  ifdef COSMIC_RAY
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerCRay",    HOFFSET(InputPara_t,Opt__Flag_LohnerCRay   ), H5T_NATIVE_INT     );
#  endif
#  endif
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerForm",    HOFFSET(InputPara_t,Opt__Flag_LohnerForm   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_User",          HOFFSET(InputPara_t,Opt__Flag_User         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_User_Num",      HOFFSET(InputPara_t,Opt__Flag_User_Num     ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Region",        HOFFSET(InputPara_t,Opt__Flag_Region       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Angular",       HOFFSET(InputPara_t,Opt__Flag_Angular      ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagAngular_CenX",        HOFFSET(InputPara_t,FlagAngular_CenX       ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "FlagAngular_CenY",        HOFFSET(InputPara_t,FlagAngular_CenY       ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "FlagAngular_CenZ",        HOFFSET(InputPara_t,FlagAngular_CenZ       ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Opt__Flag_Radial",        HOFFSET(InputPara_t,Opt__Flag_Radial       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagRadial_CenX",         HOFFSET(InputPara_t,FlagRadial_CenX        ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "FlagRadial_CenY",         HOFFSET(InputPara_t,FlagRadial_CenY        ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "FlagRadial_CenZ",         HOFFSET(InputPara_t,FlagRadial_CenZ        ), H5T_NATIVE_DOUBLE  );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Flag_NParPatch",     HOFFSET(InputPara_t,Opt__Flag_NParPatch    ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_NParCell",      HOFFSET(InputPara_t,Opt__Flag_NParCell     ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_ParMassCell",   HOFFSET(InputPara_t,Opt__Flag_ParMassCell  ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__NoFlagNearBoundary", HOFFSET(InputPara_t,Opt__NoFlagNearBoundary), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__PatchCount",         HOFFSET(InputPara_t,Opt__PatchCount        ), H5T_NATIVE_INT     );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__ParticleCount",      HOFFSET(InputPara_t,Opt__ParticleCount     ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__ReuseMemory",        HOFFSET(InputPara_t,Opt__ReuseMemory       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__MemoryPool",         HOFFSET(InputPara_t,Opt__MemoryPool        ), H5T_NATIVE_INT     );

// load balance
#  ifdef LOAD_BALANCE
   H5Tinsert( H5_TypeID, "LB_WLI_Max",              HOFFSET(InputPara_t,LB_WLI_Max             ), H5T_NATIVE_DOUBLE  );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "LB_Par_Weight",           HOFFSET(InputPara_t,LB_Par_Weight          ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "Opt__RecordLoadBalance",  HOFFSET(InputPara_t,Opt__RecordLoadBalance ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__LB_ExchangeFather",  HOFFSET(InputPara_t,Opt__LB_ExchangeFather ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__MinimizeMPIBarrier", HOFFSET(InputPara_t,Opt__MinimizeMPIBarrier), H5T_NATIVE_INT     );

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Gamma",                   HOFFSET(InputPara_t,Gamma                  ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "MolecularWeight",         HOFFSET(InputPara_t,MolecularWeight        ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "MuNorm",                  HOFFSET(InputPara_t,MuNorm                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "IsoTemp",                 HOFFSET(InputPara_t,IsoTemp                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "MinMod_Coeff",            HOFFSET(InputPara_t,MinMod_Coeff           ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "MinMod_MaxIter",          HOFFSET(InputPara_t,MinMod_MaxIter         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__LR_Limiter",         HOFFSET(InputPara_t,Opt__LR_Limiter        ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__1stFluxCorr",        HOFFSET(InputPara_t,Opt__1stFluxCorr       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__1stFluxCorrScheme",  HOFFSET(InputPara_t,Opt__1stFluxCorrScheme ), H5T_NATIVE_INT     );
#  ifdef DUAL_ENERGY
   H5Tinsert( H5_TypeID, "DualEnergySwitch",        HOFFSET(InputPara_t,DualEnergySwitch       ), H5T_NATIVE_DOUBLE  );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__SameInterfaceB",     HOFFSET(InputPara_t,Opt__SameInterfaceB    ), H5T_NATIVE_INT     );
#  endif
#  endif // HYDRO

// ELBDM solvers
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "ELBDM_Mass",              HOFFSET(InputPara_t,ELBDM_Mass             ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "ELBDM_PlanckConst",       HOFFSET(InputPara_t,ELBDM_PlanckConst      ), H5T_NATIVE_DOUBLE  );
#  ifdef QUARTIC_SELF_INTERACTION
   H5Tinsert( H5_TypeID, "ELBDM_Lambda",            HOFFSET(InputPara_t,ELBDM_Lambda           ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "ELBDM_Taylor3_Coeff",     HOFFSET(InputPara_t,ELBDM_Taylor3_Coeff    ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "ELBDM_Taylor3_Auto",      HOFFSET(InputPara_t,ELBDM_Taylor3_Auto     ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "ELBDM_RemoveMotionCM",    HOFFSET(InputPara_t,ELBDM_RemoveMotionCM   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "ELBDM_BaseSpectral",      HOFFSET(InputPara_t,ELBDM_BaseSpectral     ), H5T_NATIVE_INT     );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "ELBDM_FirstWaveLevel",    HOFFSET(InputPara_t,ELBDM_FirstWaveLevel   ), H5T_NATIVE_INT     );
#  endif
#  endif // ELBDM

// fluid solvers in different models
   H5Tinsert( H5_TypeID, "Flu_GPU_NPGroup",         HOFFSET(InputPara_t,Flu_GPU_NPGroup        ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "GPU_NStream",             HOFFSET(InputPara_t,GPU_NStream            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__FixUp_Flux",         HOFFSET(InputPara_t,Opt__FixUp_Flux        ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FixUpFlux_Var",           HOFFSET(InputPara_t,FixUpFlux_Var          ), H5T_NATIVE_LONG    );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__FixUp_Electric",     HOFFSET(InputPara_t,Opt__FixUp_Electric    ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__FixUp_Restrict",     HOFFSET(InputPara_t,Opt__FixUp_Restrict    ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FixUpRestrict_Var",       HOFFSET(InputPara_t,FixUpRestrict_Var      ), H5T_NATIVE_LONG    );
   H5Tinsert( H5_TypeID, "Opt__CorrAfterAllSync",   HOFFSET(InputPara_t,Opt__CorrAfterAllSync  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__NormalizePassive",   HOFFSET(InputPara_t,Opt__NormalizePassive  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "NormalizePassive_NVar",   HOFFSET(InputPara_t,NormalizePassive_NVar  ), H5T_NATIVE_INT     );
#  if ( NCOMP_PASSIVE > 0 )
   H5Tinsert( H5_TypeID, "NormalizePassive_VarIdx", HOFFSET(InputPara_t,NormalizePassive_VarIdx), H5_TypeID_Arr_NPassive );
#  endif
   H5Tinsert( H5_TypeID, "Opt__IntFracPassive_LR",  HOFFSET(InputPara_t,Opt__IntFracPassive_LR ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "IntFracPassive_NVar",     HOFFSET(InputPara_t,IntFracPassive_NVar    ), H5T_NATIVE_INT     );
#  if ( NCOMP_PASSIVE > 0 )
   H5Tinsert( H5_TypeID, "IntFracPassive_VarIdx",   HOFFSET(InputPara_t,IntFracPassive_VarIdx  ), H5_TypeID_Arr_NPassive );
#  endif

// store the name of all fields
   for (int v=0; v<NFieldStored; v++)
   {
//    key for each field
      sprintf( Key, "FieldLabel%02d", v );

//    assuming the offset between successive FieldLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,FieldLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
//    key for each field
      sprintf( Key, "MagLabel%02d", v );

//    assuming the offset between successive MagLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,MagLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }
#  endif

   H5Tinsert( H5_TypeID, "Opt__OverlapMPI",         HOFFSET(InputPara_t,Opt__OverlapMPI        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ResetFluid",         HOFFSET(InputPara_t,Opt__ResetFluid        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ResetFluidInit",     HOFFSET(InputPara_t,Opt__ResetFluidInit    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__FreezeFluid",        HOFFSET(InputPara_t,Opt__FreezeFluid       ), H5T_NATIVE_INT              );
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "MinDens",                 HOFFSET(InputPara_t,MinDens                ), H5T_NATIVE_DOUBLE           );
#  endif
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "MinPres",                 HOFFSET(InputPara_t,MinPres                ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "MinEint",                 HOFFSET(InputPara_t,MinEint                ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "MinTemp",                 HOFFSET(InputPara_t,MinTemp                ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "MinEntr",                 HOFFSET(InputPara_t,MinEntr                ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Opt__CheckPresAfterFlu",  HOFFSET(InputPara_t,Opt__CheckPresAfterFlu ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__LastResortFloor",    HOFFSET(InputPara_t,Opt__LastResortFloor   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "JeansMinPres",            HOFFSET(InputPara_t,JeansMinPres           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "JeansMinPres_Level",      HOFFSET(InputPara_t,JeansMinPres_Level     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "JeansMinPres_NCell",      HOFFSET(InputPara_t,JeansMinPres_NCell     ), H5T_NATIVE_INT              );
#  endif

// self-gravity
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "NewtonG",                 HOFFSET(InputPara_t,NewtonG                ), H5T_NATIVE_DOUBLE           );
#  if   ( POT_SCHEME == SOR )
   H5Tinsert( H5_TypeID, "SOR_Omega",               HOFFSET(InputPara_t,SOR_Omega              ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "SOR_MaxIter",             HOFFSET(InputPara_t,SOR_MaxIter            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "SOR_MinIter",             HOFFSET(InputPara_t,SOR_MinIter            ), H5T_NATIVE_INT              );
#  elif ( POT_SCHEME == MG )
   H5Tinsert( H5_TypeID, "MG_MaxIter",              HOFFSET(InputPara_t,MG_MaxIter             ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "MG_NPreSmooth",           HOFFSET(InputPara_t,MG_NPreSmooth          ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "MG_NPostSmooth",          HOFFSET(InputPara_t,MG_NPostSmooth         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "MG_ToleratedError",       HOFFSET(InputPara_t,MG_ToleratedError      ), H5T_NATIVE_DOUBLE           );
#  endif
   H5Tinsert( H5_TypeID, "Pot_GPU_NPGroup",         HOFFSET(InputPara_t,Pot_GPU_NPGroup        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__GraP5Gradient",      HOFFSET(InputPara_t,Opt__GraP5Gradient     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__SelfGravity",        HOFFSET(InputPara_t,Opt__SelfGravity       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ExtAcc",             HOFFSET(InputPara_t,Opt__ExtAcc            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ExtPot",             HOFFSET(InputPara_t,Opt__ExtPot            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "ExtPotTable_Name",        HOFFSET(InputPara_t,ExtPotTable_Name       ), H5_TypeID_VarStr            );
   H5Tinsert( H5_TypeID, "ExtPotTable_NPoint",      HOFFSET(InputPara_t,ExtPotTable_NPoint     ), H5_TypeID_Arr_3Int          );
   H5Tinsert( H5_TypeID, "ExtPotTable_dh",          HOFFSET(InputPara_t,ExtPotTable_dh         ), H5_TypeID_Arr_3Double       );
   H5Tinsert( H5_TypeID, "ExtPotTable_EdgeL",       HOFFSET(InputPara_t,ExtPotTable_EdgeL      ), H5_TypeID_Arr_3Double       );
   H5Tinsert( H5_TypeID, "ExtPotTable_Float8",      HOFFSET(InputPara_t,ExtPotTable_Float8     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__GravityExtraMass",   HOFFSET(InputPara_t,Opt__GravityExtraMass  ), H5T_NATIVE_INT              );
#  endif // #ifdef GRAVITY

// source terms
   H5Tinsert( H5_TypeID, "Src_Deleptonization",     HOFFSET(InputPara_t,Src_Deleptonization    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Src_User",                HOFFSET(InputPara_t,Src_User               ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Src_GPU_NPGroup",         HOFFSET(InputPara_t,Src_GPU_NPGroup        ), H5T_NATIVE_INT              );

// Grackle
#  ifdef SUPPORT_GRACKLE
   H5Tinsert( H5_TypeID, "Grackle_Activate",        HOFFSET(InputPara_t,Grackle_Activate       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Verbose",         HOFFSET(InputPara_t,Grackle_Verbose        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Cooling",         HOFFSET(InputPara_t,Grackle_Cooling        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Primordial",      HOFFSET(InputPara_t,Grackle_Primordial     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Metal",           HOFFSET(InputPara_t,Grackle_Metal          ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_UV",              HOFFSET(InputPara_t,Grackle_UV             ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_CMB_Floor",       HOFFSET(InputPara_t,Grackle_CMB_Floor      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_PE_Heating",      HOFFSET(InputPara_t,Grackle_PE_Heating     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_PE_HeatingRate",  HOFFSET(InputPara_t,Grackle_PE_HeatingRate ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Grackle_CloudyTable",     HOFFSET(InputPara_t,Grackle_CloudyTable    ), H5_TypeID_VarStr            );
   H5Tinsert( H5_TypeID, "Grackle_ThreeBodyRate",   HOFFSET(InputPara_t,Grackle_ThreeBodyRate  ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_CIE_Cooling",     HOFFSET(InputPara_t,Grackle_CIE_Cooling    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_H2_OpaApprox",    HOFFSET(InputPara_t,Grackle_H2_OpaApprox   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Che_GPU_NPGroup",         HOFFSET(InputPara_t,Che_GPU_NPGroup        ), H5T_NATIVE_INT              );
#  endif

// star formation
#  ifdef STAR_FORMATION
   H5Tinsert( H5_TypeID, "SF_CreateStar_Scheme",       HOFFSET(InputPara_t,SF_CreateStar_Scheme       ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_RSeed",        HOFFSET(InputPara_t,SF_CreateStar_RSeed        ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_DetRandom",    HOFFSET(InputPara_t,SF_CreateStar_DetRandom    ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MinLevel",     HOFFSET(InputPara_t,SF_CreateStar_MinLevel     ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MinGasDens",   HOFFSET(InputPara_t,SF_CreateStar_MinGasDens   ), H5T_NATIVE_DOUBLE    );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MassEff",      HOFFSET(InputPara_t,SF_CreateStar_MassEff      ), H5T_NATIVE_DOUBLE    );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MinStarMass",  HOFFSET(InputPara_t,SF_CreateStar_MinStarMass  ), H5T_NATIVE_DOUBLE    );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MaxStarMFrac", HOFFSET(InputPara_t,SF_CreateStar_MaxStarMFrac ), H5T_NATIVE_DOUBLE    );
#  endif

// feedback
#  ifdef FEEDBACK
   H5Tinsert( H5_TypeID, "FB_Level",                HOFFSET(InputPara_t,FB_Level               ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "FB_RSeed",                HOFFSET(InputPara_t,FB_RSeed               ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "FB_SNe",                  HOFFSET(InputPara_t,FB_SNe                 ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "FB_User",                 HOFFSET(InputPara_t,FB_User                ), H5T_NATIVE_INT              );
#  endif

// cosmic ray
#  ifdef COSMIC_RAY
   H5Tinsert( H5_TypeID, "CR_Gamma",               HOFFSET(InputPara_t,CR_Gamma               ), H5T_NATIVE_DOUBLE            );
#  ifdef CR_DIFFUSION
   H5Tinsert( H5_TypeID, "CR_Diffusion_ParaCoeff", HOFFSET(InputPara_t,CR_Diffusion_ParaCoeff ), H5T_NATIVE_DOUBLE            );
   H5Tinsert( H5_TypeID, "CR_Diffusion_PerpCoeff", HOFFSET(InputPara_t,CR_Diffusion_PerpCoeff ), H5T_NATIVE_DOUBLE            );
   H5Tinsert( H5_TypeID, "CR_Diffusion_MinB",      HOFFSET(InputPara_t,CR_Diffusion_MinB      ), H5T_NATIVE_DOUBLE            );
#  endif
#  endif // #ifdef COSMIC_RAY

// initialization
   H5Tinsert( H5_TypeID, "Opt__Init",               HOFFSET(InputPara_t,Opt__Init               ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "RestartLoadNRank",        HOFFSET(InputPara_t,RestartLoadNRank        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RestartReset",       HOFFSET(InputPara_t,Opt__RestartReset       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Level",        HOFFSET(InputPara_t,Opt__UM_IC_Level        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_NLevel",       HOFFSET(InputPara_t,Opt__UM_IC_NLevel       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_NVar",         HOFFSET(InputPara_t,Opt__UM_IC_NVar         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Format",       HOFFSET(InputPara_t,Opt__UM_IC_Format       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Float8",       HOFFSET(InputPara_t,Opt__UM_IC_Float8       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Downgrade",    HOFFSET(InputPara_t,Opt__UM_IC_Downgrade    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Refine",       HOFFSET(InputPara_t,Opt__UM_IC_Refine       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_LoadNRank",    HOFFSET(InputPara_t,Opt__UM_IC_LoadNRank    ), H5T_NATIVE_INT              );
#  if ( NLEVEL > 1 )
   H5Tinsert( H5_TypeID, "UM_IC_RefineRegion",      HOFFSET(InputPara_t,UM_IC_RefineRegion      ), H5_TypeID_Arr_NLvM1_6Int    );
#  endif
   H5Tinsert( H5_TypeID, "Opt__InitRestrict",       HOFFSET(InputPara_t,Opt__InitRestrict       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__InitGridWithOMP",    HOFFSET(InputPara_t,Opt__InitGridWithOMP    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__GPUID_Select",       HOFFSET(InputPara_t,Opt__GPUID_Select       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Init_Subsampling_NCell",  HOFFSET(InputPara_t,Init_Subsampling_NCell  ), H5T_NATIVE_INT              );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__InitBFieldByVecPot", HOFFSET(InputPara_t,Opt__InitBFieldByVecPot ), H5T_NATIVE_INT              );
#  endif
#  ifdef SUPPORT_FFTW
   H5Tinsert( H5_TypeID, "Opt__FFTW_Startup",       HOFFSET(InputPara_t,Opt__FFTW_Startup       ), H5T_NATIVE_INT              );
#  endif

// interpolation schemes
   H5Tinsert( H5_TypeID, "Opt__Int_Time",           HOFFSET(InputPara_t,Opt__Int_Time          ), H5T_NATIVE_INT              );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Int_Prim",           HOFFSET(InputPara_t,Opt__Int_Prim          ), H5T_NATIVE_INT              );
#  endif
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Opt__Int_Phase",          HOFFSET(InputPara_t,Opt__Int_Phase         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Res_Phase",          HOFFSET(InputPara_t,Opt__Res_Phase         ), H5T_NATIVE_INT              );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "Opt__Hybrid_Match_Phase", HOFFSET(InputPara_t,Opt__Hybrid_Match_Phase), H5T_NATIVE_INT              );
#  endif
#  endif // ELBDM
   H5Tinsert( H5_TypeID, "Opt__Flu_IntScheme",      HOFFSET(InputPara_t,Opt__Flu_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RefFlu_IntScheme",   HOFFSET(InputPara_t,Opt__RefFlu_IntScheme  ), H5T_NATIVE_INT              );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Mag_IntScheme",      HOFFSET(InputPara_t,Opt__Mag_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RefMag_IntScheme",   HOFFSET(InputPara_t,Opt__RefMag_IntScheme  ), H5T_NATIVE_INT              );
#  endif
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Opt__Pot_IntScheme",      HOFFSET(InputPara_t,Opt__Pot_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Rho_IntScheme",      HOFFSET(InputPara_t,Opt__Rho_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Gra_IntScheme",      HOFFSET(InputPara_t,Opt__Gra_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RefPot_IntScheme",   HOFFSET(InputPara_t,Opt__RefPot_IntScheme  ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "IntMonoCoeff",            HOFFSET(InputPara_t,IntMonoCoeff           ), H5T_NATIVE_DOUBLE           );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "IntMonoCoeffB",           HOFFSET(InputPara_t,IntMonoCoeffB          ), H5T_NATIVE_DOUBLE           );
#  endif
   H5Tinsert( H5_TypeID, "Mono_MaxIter",            HOFFSET(InputPara_t,Mono_MaxIter           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "IntOppSign0thOrder",      HOFFSET(InputPara_t,IntOppSign0thOrder     ), H5T_NATIVE_INT              );
#  ifdef SUPPORT_SPECTRAL_INT
   H5Tinsert( H5_TypeID, "SpecInt_TablePath",       HOFFSET(InputPara_t,SpecInt_TablePath      ), H5_TypeID_VarStr            );
   H5Tinsert( H5_TypeID, "SpecInt_GhostBoundary",   HOFFSET(InputPara_t,SpecInt_GhostBoundary  ), H5T_NATIVE_INT              );
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "SpecInt_XY_Instead_DePha",HOFFSET(InputPara_t,SpecInt_XY_Instead_DePha),H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "SpecInt_VortexThreshold", HOFFSET(InputPara_t,SpecInt_VortexThreshold), H5T_NATIVE_DOUBLE           );
#  endif
#  endif // #ifdef SUPPORT_SPECTRAL_INT

// data dump
   H5Tinsert( H5_TypeID, "Opt__Output_Total",           HOFFSET(InputPara_t,Opt__Output_Total          ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Part",            HOFFSET(InputPara_t,Opt__Output_Part           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_User",            HOFFSET(InputPara_t,Opt__Output_User           ), H5T_NATIVE_INT              );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Output_Par_Mode",        HOFFSET(InputPara_t,Opt__Output_Par_Mode       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Par_Mesh",        HOFFSET(InputPara_t,Opt__Output_Par_Mesh       ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Output_BasePS",          HOFFSET(InputPara_t,Opt__Output_BasePS         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Base",            HOFFSET(InputPara_t,Opt__Output_Base           ), H5T_NATIVE_INT              );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Opt__Output_Pot",             HOFFSET(InputPara_t,Opt__Output_Pot            ), H5T_NATIVE_INT              );
#  endif
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Output_ParDens",         HOFFSET(InputPara_t,Opt__Output_ParDens        ), H5T_NATIVE_INT              );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Output_CC_Mag",          HOFFSET(InputPara_t,Opt__Output_CC_Mag         ), H5T_NATIVE_INT              );
#  endif
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Output_Pres",            HOFFSET(InputPara_t,Opt__Output_Pres           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Temp",            HOFFSET(InputPara_t,Opt__Output_Temp           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Entr",            HOFFSET(InputPara_t,Opt__Output_Entr           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Cs",              HOFFSET(InputPara_t,Opt__Output_Cs             ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_DivVel",          HOFFSET(InputPara_t,Opt__Output_DivVel         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Mach",            HOFFSET(InputPara_t,Opt__Output_Mach           ), H5T_NATIVE_INT              );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Output_DivMag",          HOFFSET(InputPara_t,Opt__Output_DivMag         ), H5T_NATIVE_INT              );
#  endif
#  ifdef SRHD
   H5Tinsert( H5_TypeID, "Opt__Output_3Velocity",       HOFFSET(InputPara_t,Opt__Output_3Velocity      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Lorentz",         HOFFSET(InputPara_t,Opt__Output_Lorentz        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Enthalpy",        HOFFSET(InputPara_t,Opt__Output_Enthalpy       ), H5T_NATIVE_INT              );
#  endif
#  endif // #if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Output_UserField",       HOFFSET(InputPara_t,Opt__Output_UserField      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Mode",            HOFFSET(InputPara_t,Opt__Output_Mode           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Restart",         HOFFSET(InputPara_t,Opt__Output_Restart        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Step",            HOFFSET(InputPara_t,Opt__Output_Step           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Dt",              HOFFSET(InputPara_t,Opt__Output_Dt             ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Opt__Output_Text_Format_Flt", HOFFSET(InputPara_t,Opt__Output_Text_Format_Flt), H5_TypeID_VarStr            );
   H5Tinsert( H5_TypeID, "Opt__Output_Text_Length_Int", HOFFSET(InputPara_t,Opt__Output_Text_Length_Int), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Output_PartX",                HOFFSET(InputPara_t,Output_PartX               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Output_PartY",                HOFFSET(InputPara_t,Output_PartY               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Output_PartZ",                HOFFSET(InputPara_t,Output_PartZ               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "InitDumpID",                  HOFFSET(InputPara_t,InitDumpID                 ), H5T_NATIVE_INT              );

// libyt jupyter
#  if ( defined(SUPPORT_LIBYT) && defined(LIBYT_JUPYTER) )
   H5Tinsert( H5_TypeID, "Yt_JupyterUseConnectionFile", HOFFSET(InputPara_t,Yt_JupyterUseConnectionFile), H5T_NATIVE_INT              );
#  endif

// miscellaneous
   H5Tinsert( H5_TypeID, "Opt__Verbose",            HOFFSET(InputPara_t,Opt__Verbose           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__TimingBarrier",      HOFFSET(InputPara_t,Opt__TimingBarrier     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__TimingBalance",      HOFFSET(InputPara_t,Opt__TimingBalance     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__TimingMPI",          HOFFSET(InputPara_t,Opt__TimingMPI         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordNote",         HOFFSET(InputPara_t,Opt__RecordNote        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordUnphy",        HOFFSET(InputPara_t,Opt__RecordUnphy       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordMemory",       HOFFSET(InputPara_t,Opt__RecordMemory      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordPerformance",  HOFFSET(InputPara_t,Opt__RecordPerformance ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ManualControl",      HOFFSET(InputPara_t,Opt__ManualControl     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordCenter",       HOFFSET(InputPara_t,Opt__RecordCenter      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "COM_CenX",                HOFFSET(InputPara_t,COM_CenX               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "COM_CenY",                HOFFSET(InputPara_t,COM_CenY               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "COM_CenZ",                HOFFSET(InputPara_t,COM_CenZ               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "COM_MaxR",                HOFFSET(InputPara_t,COM_MaxR               ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "COM_MinRho",              HOFFSET(InputPara_t,COM_MinRho             ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "COM_TolErrR",             HOFFSET(InputPara_t,COM_TolErrR            ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "COM_MaxIter",             HOFFSET(InputPara_t,COM_MaxIter            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordUser",         HOFFSET(InputPara_t,Opt__RecordUser        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__OptimizeAggressive", HOFFSET(InputPara_t,Opt__OptimizeAggressive), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__SortPatchByLBIdx",   HOFFSET(InputPara_t,Opt__SortPatchByLBIdx  ), H5T_NATIVE_INT              );

// simulation checks
   H5Tinsert( H5_TypeID, "Opt__Ck_Refine",          HOFFSET(InputPara_t,Opt__Ck_Refine         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_ProperNesting",   HOFFSET(InputPara_t,Opt__Ck_ProperNesting  ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_Conservation",    HOFFSET(InputPara_t,Opt__Ck_Conservation   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "AngMom_OriginX",          HOFFSET(InputPara_t,AngMom_OriginX         ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "AngMom_OriginY",          HOFFSET(InputPara_t,AngMom_OriginY         ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "AngMom_OriginZ",          HOFFSET(InputPara_t,AngMom_OriginZ         ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Opt__Ck_NormPassive",     HOFFSET(InputPara_t,Opt__Ck_NormPassive    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_Restrict",        HOFFSET(InputPara_t,Opt__Ck_Restrict       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_Finite",          HOFFSET(InputPara_t,Opt__Ck_Finite         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_PatchAllocate",   HOFFSET(InputPara_t,Opt__Ck_PatchAllocate  ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_FluxAllocate",    HOFFSET(InputPara_t,Opt__Ck_FluxAllocate   ), H5T_NATIVE_INT              );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Ck_Negative",        HOFFSET(InputPara_t,Opt__Ck_Negative       ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Ck_MemFree",         HOFFSET(InputPara_t,Opt__Ck_MemFree        ), H5T_NATIVE_DOUBLE           );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Ck_Particle",        HOFFSET(InputPara_t,Opt__Ck_Particle       ), H5T_NATIVE_INT              );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Ck_InterfaceB",      HOFFSET(InputPara_t,Opt__Ck_InterfaceB     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_DivergenceB",     HOFFSET(InputPara_t,Opt__Ck_DivergenceB    ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Ck_InputFluid",      HOFFSET(InputPara_t,Opt__Ck_InputFluid     ), H5T_NATIVE_INT              );

// flag tables
#  if ( NLEVEL > 1 )
   H5Tinsert( H5_TypeID, "FlagTable_Rho",          HOFFSET(InputPara_t,FlagTable_Rho           ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_RhoGradient",  HOFFSET(InputPara_t,FlagTable_RhoGradient   ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_Lohner",       HOFFSET(InputPara_t,FlagTable_Lohner        ), H5_TypeID_Arr_NLvM1_5Double );
   H5Tinsert( H5_TypeID, "FlagTable_Angular",      HOFFSET(InputPara_t,FlagTable_Angular       ), H5_TypeID_Arr_NLvM1_3Double );
   H5Tinsert( H5_TypeID, "FlagTable_Radial",       HOFFSET(InputPara_t,FlagTable_Radial        ), H5_TypeID_Arr_NLvM1Double   );

// store the user-defined thresholds at all levels
   for (int lv=0; lv<MAX_LEVEL; lv++)
   {
//    key for each level
      sprintf( Key, "FlagTable_User_Lv%02d", lv );

//    assuming the offset between successive FlagTable_User pointers is "PtrSize_hvl", which is equal to "sizeof( hvl_t )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,FlagTable_User)+lv*PtrSize_hvl, H5_TypeID_Arr_NLvM1_VLDouble );
   }

#  if   ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "FlagTable_PresGradient", HOFFSET(InputPara_t,FlagTable_PresGradient  ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_Vorticity",    HOFFSET(InputPara_t,FlagTable_Vorticity     ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_Jeans",        HOFFSET(InputPara_t,FlagTable_Jeans         ), H5_TypeID_Arr_NLvM1Double   );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "FlagTable_Current",      HOFFSET(InputPara_t,FlagTable_Current       ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  ifdef SRHD
   H5Tinsert( H5_TypeID, "FlagTable_LrtzGradient", HOFFSET(InputPara_t,FlagTable_LrtzGradient  ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  ifdef COSMIC_RAY
   H5Tinsert( H5_TypeID, "FlagTable_CRay",         HOFFSET(InputPara_t,FlagTable_CRay          ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  elif ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "FlagTable_EngyDensity",  HOFFSET(InputPara_t,FlagTable_EngyDensity   ), H5_TypeID_Arr_NLvM1_2Double );
   H5Tinsert( H5_TypeID, "FlagTable_Spectral",     HOFFSET(InputPara_t,FlagTable_Spectral      ), H5_TypeID_Arr_NLvM1_2Double );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   H5Tinsert( H5_TypeID, "FlagTable_Interference", HOFFSET(InputPara_t,FlagTable_Interference  ), H5_TypeID_Arr_NLvM1_4Double );
#  endif
#  endif // MODEL
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "FlagTable_NParPatch",    HOFFSET(InputPara_t,FlagTable_NParPatch     ), H5_TypeID_Arr_NLvM1Int      );
   H5Tinsert( H5_TypeID, "FlagTable_NParCell",     HOFFSET(InputPara_t,FlagTable_NParCell      ), H5_TypeID_Arr_NLvM1Int      );
   H5Tinsert( H5_TypeID, "FlagTable_ParMassCell",  HOFFSET(InputPara_t,FlagTable_ParMassCell   ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  endif // #if ( NLEVEL > 1 )

// user-defined derived fields
   H5Tinsert( H5_TypeID, "UserDerField_Num",        HOFFSET(InputPara_t,UserDerField_Num       ), H5T_NATIVE_INT              );

// --> only need to insert UserDerField_Num strings even though *UserDerField_Label/Unit are pointer arrays with
//     a size DER_NOUT_MAX in InputPara_t
// --> it should be fine as long as the "offset" (i.e., HOFFSET(InputPara_t,UserDerField_Label[0])+v*PtrSize) is correct
   for (int v=0; v<UserDerField_Num; v++)
   {
//    key for each field label
      sprintf( Key, "UserDerField_Label%02d", v );

//    assuming the offset between successive UserDerField_Label pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,UserDerField_Label[0])+v*PtrSize, H5_TypeID_VarStr );
   }

   for (int v=0; v<UserDerField_Num; v++)
   {
//    key for each field unit
      sprintf( Key, "UserDerField_Unit%02d", v );

//    assuming the offset between successive UserDerField_Unit pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,UserDerField_Unit [0])+v*PtrSize, H5_TypeID_VarStr );
   }


// free memory
   H5_Status = H5Tclose( H5_TypeID_Arr_3Int           );
   H5_Status = H5Tclose( H5_TypeID_Arr_6Int           );
#  if ( NCOMP_PASSIVE > 0 )
   H5_Status = H5Tclose( H5_TypeID_Arr_NPassive       );
#  endif
#  if ( NLEVEL > 1 )
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1Int       );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1Double    );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_2Double  );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_3Double  );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_4Double  );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_5Double  );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_VLDouble );
#  endif
   H5_Status = H5Tclose( H5_TypeID_Arr_3Double        );
   H5_Status = H5Tclose( H5_TypeID_VarStr             );

} // FUNCTION : GetCompound_InputPara



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_General
// Description :  Create the HDF5 compound datatype for HDF5_Output_t
//
// Note        :  1. HDF5_Output_t is defined in HDF5_Typedef.h
//                2. Support int, long, uint, ulong, float, double, and string datatypes
//
// Parameter   :  H5_TypeID   : HDF5 type ID for storing the compound datatype
//                HDF5_Output : Structure storing all parameters to be written
//-------------------------------------------------------------------------------------------------------
void GetCompound_General( hid_t &H5_TypeID, const HDF5_Output_t *HDF5_Output )
{

   if ( HDF5_Output->TotalSize == 0 )   Aux_Error( ERROR_INFO, "HDF5_Output_t structure must not be empty !!\n" );

   herr_t H5_Status;
   hid_t  H5_Type;
   hid_t  H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );

   H5_Status = H5Tset_size( H5_TypeID_VarStr, MAX_STRING ); // H5T_VARIABLE will cause segmentation fault
   H5_TypeID = H5Tcreate( H5T_COMPOUND, HDF5_Output->TotalSize );

   size_t offset = 0;

   for (int i=0; i<HDF5_Output->NPara; i++)
   {
      const int type = HDF5_Output->Type[i];
      switch ( type )
      {
//       must match TYPE_* defined in HDF5_Typedef.h
         case 1: H5_Type = H5T_NATIVE_INT;    break;
         case 2: H5_Type = H5T_NATIVE_LONG;   break;
         case 3: H5_Type = H5T_NATIVE_UINT;   break;
         case 4: H5_Type = H5T_NATIVE_ULONG;  break;
         case 5: H5_Type = H5T_NATIVE_INT;    break; // bool is stored as int
         case 6: H5_Type = H5T_NATIVE_FLOAT;  break;
         case 7: H5_Type = H5T_NATIVE_DOUBLE; break;
         case 8: H5_Type = H5_TypeID_VarStr;  break;
         default: Aux_Error( ERROR_INFO, "Unrecognized type: %d !!\n", type ); break;
      } // switch ( type )

      H5Tinsert( H5_TypeID, HDF5_Output->Key[i], offset, H5_Type );
      offset += HDF5_Output->TypeSize[i];

   } // for (int i=0; i<HDF5_Output->NPara; i++)
   H5_Status = H5Tclose( H5_TypeID_VarStr );

} // FUNCTION : GetCompound_General



//-------------------------------------------------------------------------------------------------------
// Function    :  H5_write_compound
// Description :  Write all parameters in HDF5_Output to an HDF5 compound dataset
//
// Note        :  1. HDF5_Output_t is defined in HDF5_Typedef.h
//                2. Support int, long, uint, ulong, float, double, and string datatypes
//
// Parameter   :  H5_SetID    : HDF5 dataset ID
//                H5_TypeID   : HDF5 compound type ID associated with HDF5_Output
//                HDF5_Output : Structure storing all parameters to be written
//
// Return      :  H5_Status_write : Status of the write operation
//-------------------------------------------------------------------------------------------------------
herr_t H5_write_compound( const hid_t H5_SetID, const hid_t H5_TypeID, const HDF5_Output_t *HDF5_Output )
{

   herr_t H5_Status, H5_Status_write;
   hid_t  H5_TypeID_VarStr;
   hid_t  H5_Type;

   H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
   H5_Status        = H5Tset_size( H5_TypeID_VarStr, MAX_STRING ); // H5T_VARIABLE will cause segmentation fault

   char   *data = new char [HDF5_Output->TotalSize];
   size_t offset = 0;

   for (int i=0; i<HDF5_Output->NPara; i++)
   {
      const int    type      = HDF5_Output->Type    [i];
      const size_t type_size = HDF5_Output->TypeSize[i];
      switch ( type )
      {
//       must match TYPE_* defined in HDF5_Typedef.h
         case 1: H5_Type = H5T_NATIVE_INT;    break;
         case 2: H5_Type = H5T_NATIVE_LONG;   break;
         case 3: H5_Type = H5T_NATIVE_UINT;   break;
         case 4: H5_Type = H5T_NATIVE_ULONG;  break;
         case 5: H5_Type = H5T_NATIVE_INT;    break; // bool is stored as int
         case 6: H5_Type = H5T_NATIVE_FLOAT;  break;
         case 7: H5_Type = H5T_NATIVE_DOUBLE; break;
         case 8: H5_Type = H5_TypeID_VarStr;  break;
         default: Aux_Error( ERROR_INFO, "Unrecognized type: %d !!\n", type ); break;
      } // switch ( type )

      if ( type == 5 )
      {
//       convert int to bool
         const int temp = ( *(bool *)(HDF5_Output->Ptr[i]) ) ? 1 : 0;

         memcpy( data+offset, &temp, type_size );
      }
      else
         memcpy( data+offset, HDF5_Output->Ptr[i], type_size );

      offset += type_size;
   } // for (int i=0; i<HDF5_Output->NPara; i++)

   H5_Status_write = H5Dwrite( H5_SetID, H5_TypeID, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );

   H5_Status = H5Tclose( H5_TypeID_VarStr );

   delete [] data;

   return H5_Status_write;

} // FUNCTION : H5_write_compound



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_HDF5_UserPara_Template
// Description :  Template for storing user-specified parameters in an HDF5 snapshot at User/UserPara
//
// Note         : 1. This function is only called by the root MPI rank
//                2. Support int, uint, long, ulong, bool, float, double, and string datatypes
//                3. HDF5_UserPara MUST store at least one parameter
//                4. The data pointer (i.e., the second argument passed to HDF5_UserPara->Add()) MUST persist outside this function (e.g., global variables)
//                5. Linked to the function pointer Output_HDF5_UserPara_Ptr
//
// Parameter   :  HDF5_UserPara : Structure storing all parameters to be written
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_HDF5_UserPara_Template( HDF5_Output_t *HDF5_UserPara )
{

// HDF5_UserPara->Add( "Your_Data_Label", &Your_Data_Pointer );

} // FUNCTION : Output_HDF5_UserPara_Template



#endif // #ifdef SUPPORT_HDF5
