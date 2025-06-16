#ifdef SUPPORT_HDF5

#include "GAMER.h"
#include "HDF5_Typedef.h"
#include <ctime>

void FillIn_HDF5_KeyInfo  ( HDF5_Output_t *HDF5_KeyInfo, int NFieldStored, const bool Load_RS,
                            const int RS_FormatVersion, const bool ReenablePar );
void FillIn_HDF5_Makefile ( HDF5_Output_t *HDF5_Makefile, const int RS_FormatVersion );
void FillIn_HDF5_SymConst ( HDF5_Output_t *HDF5_SymConst, const int RS_FormatVersion );
void FillIn_HDF5_InputPara( HDF5_Output_t *HDF5_InputPara, const int NFieldStored, char FieldLabelOut[][MAX_STRING],
                            const bool Load_RS );

void (*Output_HDF5_InputTest_Ptr)( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest ) = NULL;
void (*Output_HDF5_UserPara_Ptr)( HDF5_Output_t *HDF5_UserPara ) = NULL;

static void Output_HDF5_UserPara_Template( HDF5_Output_t *HDF5_UserPara );

static int H5_Arr_3[1] = { 3 };
static int H5_Arr_6[1] = { 6 };



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
1. Edit "FillIn_XXX" to fill in the new variables
2. Update FormatVersion
======================================================================================================*/




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Total_HDF5 (FormatVersion = 2504)
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
//                2504 : 2025/04/29 --> output OPT__PAR_INIT_CHECK
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
      const bool Load_RS_No = false; // not loading the restart file

//    3-1. collect all information to be recorded
      HDF5_Output_t HDF5_KeyInfo, HDF5_Makefile, HDF5_SymConst, HDF5_InputPara, HDF5_InputTest, HDF5_UserPara;

      FillIn_HDF5_KeyInfo  ( &HDF5_KeyInfo, NFieldStored, Load_RS_No, NULL_INT, NULL_BOOL );
      FillIn_HDF5_Makefile ( &HDF5_Makefile, NULL_INT );
      FillIn_HDF5_SymConst ( &HDF5_SymConst, NULL_INT );
      FillIn_HDF5_InputPara( &HDF5_InputPara, NFieldStored, FieldLabelOut, Load_RS_No );
      if ( Output_HDF5_InputTest_Ptr != NULL )  Output_HDF5_InputTest_Ptr( LOAD_HDF5_OUTPUT, NULL, &HDF5_InputTest );
      if ( Output_HDF5_UserPara_Ptr  != NULL )  Output_HDF5_UserPara_Ptr( &HDF5_UserPara );


//    3-2. create the "compound" datatype
      HDF5_KeyInfo.GetCompound( H5_TypeID_Com_KeyInfo );
      HDF5_Makefile.GetCompound( H5_TypeID_Com_Makefile );
      HDF5_SymConst.GetCompound( H5_TypeID_Com_SymConst );
      HDF5_InputPara.GetCompound( H5_TypeID_Com_InputPara );
      if ( Output_HDF5_InputTest_Ptr != NULL )  HDF5_InputTest.GetCompound( H5_TypeID_Com_InputTest );
      if ( Output_HDF5_UserPara_Ptr  != NULL )  HDF5_UserPara.GetCompound( H5_TypeID_Com_UserPara );


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
      H5_Status          = HDF5_KeyInfo.Write2CompoundDataset( H5_SetID_KeyInfo, H5_TypeID_Com_KeyInfo );
      H5_Status          = H5Dclose( H5_SetID_KeyInfo );

//    3-4-2. Makefile
      H5_SetID_Makefile  = H5Dcreate( H5_GroupID_Info, "Makefile", H5_TypeID_Com_Makefile, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_Makefile < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Makefile" );
      H5_Status          = HDF5_Makefile.Write2CompoundDataset( H5_SetID_Makefile, H5_TypeID_Com_Makefile );
      H5_Status          = H5Dclose( H5_SetID_Makefile );

//    3-4-3. SymConst
      H5_SetID_SymConst  = H5Dcreate( H5_GroupID_Info, "SymConst", H5_TypeID_Com_SymConst, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_SymConst < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "SymConst" );
      H5_Status          = HDF5_SymConst.Write2CompoundDataset( H5_SetID_SymConst, H5_TypeID_Com_SymConst );
      H5_Status          = H5Dclose( H5_SetID_SymConst );

//    3-4-4. InputPara
      H5_SetID_InputPara = H5Dcreate( H5_GroupID_Info, "InputPara", H5_TypeID_Com_InputPara, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_InputPara < 0 ) Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "InputPara" );
      H5_Status          = HDF5_InputPara.Write2CompoundDataset( H5_SetID_InputPara, H5_TypeID_Com_InputPara );
      H5_Status          = H5Dclose( H5_SetID_InputPara );

//    3-4-5. InputTest
      if ( Output_HDF5_InputTest_Ptr != NULL )
      {
      H5_SetID_InputTest = H5Dcreate( H5_GroupID_Info, "InputTest", H5_TypeID_Com_InputTest, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_InputTest < 0 ) Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "InputTest" );
      H5_Status          = HDF5_InputTest.Write2CompoundDataset( H5_SetID_InputTest, H5_TypeID_Com_InputTest );
      H5_Status          = H5Dclose( H5_SetID_InputTest );
      }

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
         H5_Status         = HDF5_UserPara.Write2CompoundDataset( H5_SetID_UserPara, H5_TypeID_Com_UserPara );
         H5_Status         = H5Dclose( H5_SetID_UserPara );
      } // if ( Output_HDF5_UserPara_Ptr != NULL )

      H5_Status = H5Gclose( H5_GroupID_User );
      H5_Status = H5Fclose( H5_FileID );
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
// Function    :  FillIn_HDF5_KeyInfo
// Description :  Fill in the HDF5_Output_t structure with the important information
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  HDF5_KeyInfo     : Pointer to be filled in
//                NFieldStored     : Number of grid fields to be stored on disk
//                Load_RS          : Whether the structure is used for checking the restart data
//                RS_FormatVersion : The format version of the restart data
//                ReenablePar      : Whether re-enableing PARTICLE (for restart only)
//-------------------------------------------------------------------------------------------------------
void FillIn_HDF5_KeyInfo( HDF5_Output_t *HDF5_KeyInfo, int NFieldStored, const bool Load_RS,
                          const int RS_FormatVersion, const bool ReenablePar )
{

   const time_t CalTime = time( NULL );   // calendar time

   const bool Compare_Yes = true, Compare_No = false;
   const bool   Fatal_Yes = true,   Fatal_No = false;


// format version
// ----------------------------------------------------------------------------------------------------
   HDF5_KeyInfo->Add( "FormatVersion", 2504,            0, NULL, Compare_No, NULL_BOOL, NULL_BOOL );


// physics modules
// ----------------------------------------------------------------------------------------------------
   HDF5_KeyInfo->Add( "Model",         MODEL,           0, NULL, Compare_Yes, Fatal_Yes, Fatal_Yes );
   HDF5_KeyInfo->Add( "NLevel",        NLEVEL,          0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "NCompFluid",    NCOMP_FLUID,     0, NULL, Compare_Yes, Fatal_No,  Fatal_Yes );
   HDF5_KeyInfo->Add( "NCompPassive",  NCOMP_PASSIVE,   0, NULL, Compare_Yes, Fatal_No,  Fatal_Yes );
   HDF5_KeyInfo->Add( "PatchSize",     PS1,             0, NULL, Compare_Yes, Fatal_Yes, Fatal_Yes );
   HDF5_KeyInfo->Add( "DumpID",       &DumpID,          0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "Step",         &Step,            0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
#  ifdef GRAVITY
   HDF5_KeyInfo->Add( "AveDens_Init", &AveDensity_Init, 0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "Gravity",       1,               0, NULL, Compare_Yes, Fatal_Yes, Fatal_Yes );
#  else
   HDF5_KeyInfo->Add( "Gravity",       0,               0, NULL, Compare_Yes, Fatal_Yes, Fatal_Yes );
#  endif
#  ifdef PARTICLE
   HDF5_KeyInfo->Add( "Particle",      1,               0, NULL, Compare_Yes, Fatal_Yes, Fatal_No  );
#  else
   HDF5_KeyInfo->Add( "Particle",      0,               0, NULL, Compare_Yes, Fatal_Yes, Fatal_No  );
#  endif
#  ifdef FLOAT8
   HDF5_KeyInfo->Add( "Float8",        1,               0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
#  else
   HDF5_KeyInfo->Add( "Float8",        0,               0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
#  endif
   HDF5_KeyInfo->Add( "NFieldStored",  NFieldStored,    0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "NMagStored",    NCOMP_MAG,       0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );


// particle
// ----------------------------------------------------------------------------------------------------
#  ifdef PARTICLE
   const bool Compare_Par      = ( !ReenablePar ) ? Compare_Yes : Compare_No;
   const bool Fatal_ParNAttFlt = ( !ReenablePar  &&  RS_FormatVersion >= 2300 ) ? Fatal_Yes : Fatal_No;
   const bool Fatal_ParNAttInt = ( !ReenablePar  &&  RS_FormatVersion >= 2500 ) ? Fatal_Yes : Fatal_No;
   HDF5_KeyInfo->Add( "Par_NPar",          &amr->Par->NPar_Active_AllRank, 0, NULL, Compare_Par, Fatal_No,         Fatal_No  );
   HDF5_KeyInfo->Add( "Par_NAttFltStored",  PAR_NATT_FLT_STORED,           0, NULL, Compare_Par, Fatal_ParNAttFlt, Fatal_No  );
   HDF5_KeyInfo->Add( "Par_NAttIntStored",  PAR_NATT_INT_STORED,           0, NULL, Compare_Par, Fatal_ParNAttInt, Fatal_No  );
#  ifdef FLOAT8_PAR
   HDF5_KeyInfo->Add( "Float8_Par",         1,                             0, NULL, Compare_Yes, Fatal_No,         Fatal_No  );
#  else
   HDF5_KeyInfo->Add( "Float8_Par",         0,                             0, NULL, Compare_Yes, Fatal_No,         Fatal_No  );
#  endif
#  ifdef INT8_PAR
   HDF5_KeyInfo->Add( "Int8_Par",           1,                             0, NULL, Compare_Yes, Fatal_No,         Fatal_No  );
#  else
   HDF5_KeyInfo->Add( "Int8_Par",           0,                             0, NULL, Compare_Yes, Fatal_No,         Fatal_No  );
#  endif
#  endif // #ifdef PARTICLE


// HYDRO physics modules
// ----------------------------------------------------------------------------------------------------
#  if ( MODEL == HYDRO )
   const bool Compare_MHD  = ( RS_FormatVersion >= 2400 ) ? Compare_Yes : Compare_No;
   const bool Compare_SRHD = ( RS_FormatVersion >= 2473 ) ? Compare_Yes : Compare_No;
   const bool Compare_CR   = ( RS_FormatVersion >= 2421 ) ? Compare_Yes : Compare_No;
#  ifdef MHD
   HDF5_KeyInfo->Add( "Magnetohydrodynamics", 1, 0, NULL, Compare_MHD,  Fatal_Yes, Fatal_Yes );
#  else
   HDF5_KeyInfo->Add( "Magnetohydrodynamics", 0, 0, NULL, Compare_MHD,  Fatal_Yes, Fatal_Yes );
#  endif
#  ifdef SRHD
   HDF5_KeyInfo->Add( "SRHydrodynamics",      1, 0, NULL, Compare_SRHD, Fatal_Yes, Fatal_Yes );
#  else
   HDF5_KeyInfo->Add( "SRHydrodynamics",      0, 0, NULL, Compare_SRHD, Fatal_Yes, Fatal_Yes );
#  endif
#  ifdef COSMIC_RAY
   HDF5_KeyInfo->Add( "CosmicRay",            1, 0, NULL, Compare_CR,   Fatal_Yes, Fatal_Yes );
#  ifdef CR_DIFFUSION
   HDF5_KeyInfo->Add( "CR_Diffusion",         1, 0, NULL, Compare_Yes,  Fatal_Yes, Fatal_No  );
#  else
   HDF5_KeyInfo->Add( "CR_Diffusion",         0, 0, NULL, Compare_Yes,  Fatal_Yes, Fatal_No  );
#  endif
#  else // #ifdef COSMIC_RAY
   HDF5_KeyInfo->Add( "CosmicRay",            0, 0, NULL, Compare_CR,   Fatal_Yes, Fatal_Yes );
#  endif // #ifdef COSMIC_RAY .. else ...
#  endif // #if ( MODEL == HYDRO )


// simulation
// ----------------------------------------------------------------------------------------------------
   int H5_Arr_NLv[1] = { NLEVEL };

   HDF5_KeyInfo->Add( "NX0",            &NX0_TOT[0],            1, H5_Arr_3,   Compare_Yes, Fatal_Yes, Fatal_Yes );
   HDF5_KeyInfo->Add( "BoxScale",       &amr->BoxScale[0],      1, H5_Arr_3,   Compare_Yes, Fatal_Yes, Fatal_No  );
   HDF5_KeyInfo->Add( "BoxSize",        &amr->BoxSize[0],       1, H5_Arr_3,   Compare_Yes, Fatal_Yes, Fatal_Yes );
   HDF5_KeyInfo->Add( "Time",           &Time[0],               1, H5_Arr_NLv, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "CellSize",       &amr->dh[0],            1, H5_Arr_NLv, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "CellScale",      &amr->scale[0],         1, H5_Arr_NLv, Compare_Yes, Fatal_Yes, Fatal_No  );
   HDF5_KeyInfo->Add( "NPatch",         &NPatchTotal[0],        1, H5_Arr_NLv, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "AdvanceCounter", &AdvanceCounter[0],     1, H5_Arr_NLv, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "dTime_AllLv",    &dTime_AllLv[0],        1, H5_Arr_NLv, Compare_No,  NULL_BOOL, NULL_BOOL );
#  if ( MODEL == ELBDM  &&  ELBDM_SCHEME == ELBDM_HYBRID )
   HDF5_KeyInfo->Add( "UseWaveScheme",  &amr->use_wave_flag[0], 1, H5_Arr_NLv, Compare_No,  NULL_BOOL, NULL_BOOL );
#  endif


// code information
// ----------------------------------------------------------------------------------------------------
   char *temp;
   char CodeVersion[MAX_STRING] = VERSION;
   char DumpWalltime[MAX_STRING], GitBranch[MAX_STRING], GitCommit[MAX_STRING];

   temp = ctime( &CalTime );
   temp[ strlen(temp)-1 ] = '\0';  // remove the last character '\n'
   strncpy( DumpWalltime, temp, MAX_STRING );

   temp = (char*)EXPAND_AND_QUOTE( GIT_BRANCH );
   strncpy( GitBranch, temp, MAX_STRING );

   temp = (char*)EXPAND_AND_QUOTE( GIT_COMMIT );
   strncpy( GitCommit, temp, MAX_STRING );

   //###REVISE: replace rand() by UUID
   srand( time(NULL) );
   long UniqueDataID = rand();

   HDF5_KeyInfo->Add( "CodeVersion",  &CodeVersion[0],  0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_KeyInfo->Add( "DumpWallTime", &DumpWalltime[0], 0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_KeyInfo->Add( "GitBranch",    &GitBranch[0],    0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_KeyInfo->Add( "GitCommit",    &GitCommit[0],    0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_KeyInfo->Add( "UniqueDataID", &UniqueDataID,    0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );


// conserved variables
// ----------------------------------------------------------------------------------------------------
   int H5_Arr_NCM_P1[1] = { 1+NCONREF_MAX };

   if ( !ConRefInitialized  &&  !Load_RS )   Aux_Error( ERROR_INFO, "Reference values for conserved variables have not been assigned yet !!\n" );
   HDF5_KeyInfo->Add( "ConRef", &ConRef[0], 1, H5_Arr_NCM_P1, Compare_No, NULL_BOOL, NULL_BOOL );

} // FUNCTION : FillIn_HDF5_KeyInfo



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_HDF5_Makefile
// Description :  Fill in the HDF5_Output_t structure with the makefile information
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  HDF5_Makefile    : Pointer to be filled in
//                RS_FormatVersion : The format version of the restart data
//-------------------------------------------------------------------------------------------------------
void FillIn_HDF5_Makefile( HDF5_Output_t *HDF5_Makefile, const int RS_FormatVersion )
{

   const bool Compare_Yes = true, Compare_No = false;
   const bool   Fatal_Yes = true,   Fatal_No = false;


// physics modules options
// ----------------------------------------------------------------------------------------------------
   HDF5_Makefile->Add( "Model",    MODEL, 0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );

#  ifdef GRAVITY
   HDF5_Makefile->Add( "Gravity",  1,     0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  else
   HDF5_Makefile->Add( "Gravity",  0,     0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  endif

#  ifdef COMOVING
   HDF5_Makefile->Add( "Comoving", 1,     0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  else
   HDF5_Makefile->Add( "Comoving", 0,     0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  endif

#  ifdef PARTICLE
   HDF5_Makefile->Add( "Particle", 1,     0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  else
   HDF5_Makefile->Add( "Particle", 0,     0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif


// HYDRO options
// ----------------------------------------------------------------------------------------------------
#  if ( MODEL == HYDRO )
   HDF5_Makefile->Add( "FluScheme",            FLU_SCHEME,  0, NULL, Compare_Yes, Fatal_No, Fatal_No  );

#  ifdef LR_SCHEME
   HDF5_Makefile->Add( "LRScheme",             LR_SCHEME,   0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif

#  ifdef RSOLVER
   HDF5_Makefile->Add( "RSolver",              RSOLVER,     0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif

#  ifdef DUAL_ENERGY
   HDF5_Makefile->Add( "DualEnergy",           DUAL_ENERGY, 0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  else
   HDF5_Makefile->Add( "DualEnergy",           0,           0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif

#  ifdef MHD
   HDF5_Makefile->Add( "Magnetohydrodynamics", 1,           0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  else
   HDF5_Makefile->Add( "Magnetohydrodynamics", 0,           0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  endif

#  ifdef SRHD
   HDF5_Makefile->Add( "SRHydrodynamics",      1,           0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  else
   HDF5_Makefile->Add( "SRHydrodynamics",      0,           0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  endif

#  ifdef COSMIC_RAY
   HDF5_Makefile->Add( "CosmicRay",            1,           0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  ifdef CR_DIFFUSION
   HDF5_Makefile->Add( "CR_Diffusion",         1,           0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  else
   HDF5_Makefile->Add( "CR_Diffusion",         0,           0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif
#  else // #ifdef COSMIC_RAY
   HDF5_Makefile->Add( "CosmicRay",            0,           0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  endif // #ifdef COSMIC_RAY .. else ...

   HDF5_Makefile->Add( "EoS",                  EOS,         0, NULL, Compare_Yes, Fatal_No, Fatal_No  );

#  ifdef BAROTROPIC_EOS
   HDF5_Makefile->Add( "BarotropicEoS",        1,           0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  else
   HDF5_Makefile->Add( "BarotropicEoS",        0,           0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif

#  endif // #if ( MODEL == HYDRO )


// ELBDM options
// ----------------------------------------------------------------------------------------------------
#  if ( MODEL == ELBDM )
   HDF5_Makefile->Add( "ELBDMScheme",      ELBDM_SCHEME, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_Makefile->Add( "WaveScheme",       WAVE_SCHEME,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  ifdef CONSERVE_MASS
   HDF5_Makefile->Add( "ConserveMass",     1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "ConserveMass",     0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef LAPLACIAN_4TH
   HDF5_Makefile->Add( "Laplacian4th",     1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "Laplacian4th",     0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef QUARTIC_SELF_INTERACTION
   HDF5_Makefile->Add( "SelfInteraction4", 1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "SelfInteraction4", 0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // #if ( MODEL == ELBDM )


// gravity options
// ----------------------------------------------------------------------------------------------------
#  ifdef GRAVITY
   HDF5_Makefile->Add( "PotScheme",      POT_SCHEME, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  ifdef STORE_POT_GHOST
   HDF5_Makefile->Add( "StorePotGhost",  1,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "StorePotGhost",  0,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef UNSPLIT_GRAVITY
   HDF5_Makefile->Add( "UnsplitGravity", 1,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "UnsplitGravity", 0,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  endif // #ifdef GRAVITY


// particle options
// ----------------------------------------------------------------------------------------------------
#  ifdef PARTICLE
#  ifdef MASSIVE_PARTICLES
   HDF5_Makefile->Add( "MassiveParticles", 1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "MassiveParticles", 0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif
#  ifdef TRACER
   HDF5_Makefile->Add( "Tracer",           1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "Tracer",           0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif
#  ifdef STORE_PAR_ACC
   HDF5_Makefile->Add( "StoreParAcc",      1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "StoreParAcc",      0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif

#  ifdef STAR_FORMATION
   HDF5_Makefile->Add( "StarFormation",    1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "StarFormation",    0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif

#  ifdef FEEDBACK
   HDF5_Makefile->Add( "Feedback",         1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "Feedback",         0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif

   const bool Fatal_ParNAttFlt = ( RS_FormatVersion > 2300 ) ? Fatal_Yes : Fatal_No;
   const bool Fatal_ParNAttInt = ( RS_FormatVersion > 2500 ) ? Fatal_Yes : Fatal_No;
   HDF5_Makefile->Add( "Par_NAttFltUser",  PAR_NATT_FLT_USER, 0, NULL, Compare_Yes, Fatal_No, Fatal_ParNAttFlt );
   HDF5_Makefile->Add( "Par_NAttIntUser",  PAR_NATT_INT_USER, 0, NULL, Compare_Yes, Fatal_No, Fatal_ParNAttInt );

#  ifdef FLOAT8_PAR
   HDF5_Makefile->Add( "Float8_Par",       1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "Float8_Par",       0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif

#  ifdef INT8_PAR
   HDF5_Makefile->Add( "Int8_Par",         1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_Makefile->Add( "Int8_Par",         0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif
#  endif // #ifdef PARTICLE


// GRACKLE options
// ----------------------------------------------------------------------------------------------------
#  ifdef SUPPORT_GRACKLE
   HDF5_Makefile->Add( "SupportGrackle", 1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "SupportGrackle", 0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// miscellaneous options
// ----------------------------------------------------------------------------------------------------
   HDF5_Makefile->Add( "NLevel",                 NLEVEL,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_Makefile->Add( "MaxPatch",               MAX_PATCH,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  ifdef GAMER_DEBUG
   HDF5_Makefile->Add( "GAMER_Debug",            1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "GAMER_Debug",            0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef BITWISE_REPRODUCIBILITY
   HDF5_Makefile->Add( "BitwiseReproducibility", 1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "BitwiseReproducibility", 0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef TIMING
   HDF5_Makefile->Add( "Timing",                 1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "Timing",                 0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef TIMING_SOLVER
   HDF5_Makefile->Add( "TimingSolver",           1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "TimingSolver",           0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef FLOAT8
   HDF5_Makefile->Add( "Float8",                 1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "Float8",                 0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef LAOHU
   HDF5_Makefile->Add( "Laohu",                  1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "Laohu",                  0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef SUPPORT_HDF5
   HDF5_Makefile->Add( "SupportHDF5",            1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "SupportHDF5",            0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef SUPPORT_GSL
   HDF5_Makefile->Add( "SupportGSL",             1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "SupportGSL",             0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef SUPPORT_SPECTRAL_INT
   HDF5_Makefile->Add( "SupportSpectralInt",     1,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "SupportSpectralInt",     0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef SUPPORT_FFTW
   HDF5_Makefile->Add( "SupportFFTW",            SUPPORT_FFTW,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "SupportFFTW",            0,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

   HDF5_Makefile->Add( "RandomNumber",           RANDOM_NUMBER, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// LibYT options
// ----------------------------------------------------------------------------------------------------
#  ifdef SUPPORT_LIBYT
   HDF5_Makefile->Add( "SupportLibYT",       1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  ifdef LIBYT_USE_PATCH_GROUP
   HDF5_Makefile->Add( "LibYTUsePatchGroup", 1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "LibYTUsePatchGroup", 0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef LIBYT_INTERACTIVE
   HDF5_Makefile->Add( "LibYTInteractive",   1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "LibYTInteractive",   0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef LIBYT_RELOAD
   HDF5_Makefile->Add( "LibYTReload",        1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "LibYTReload",        0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef LIBYT_JUPYTER
   HDF5_Makefile->Add( "LibYTJupyter",       1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "LibYTJupyter",       0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  else  // #ifdef SUPPORT_LIBYT
   HDF5_Makefile->Add( "SupportLibYT",       0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif // #ifdef SUPPORT_LIBYT ... else ...


// parallel options
// ----------------------------------------------------------------------------------------------------
#  ifdef SERIAL
   HDF5_Makefile->Add( "Serial",      1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "Serial",      0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef LOAD_BALANCE
   HDF5_Makefile->Add( "LoadBalance", LOAD_BALANCE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "LoadBalance", 0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef OVERLAP_MPI
   HDF5_Makefile->Add( "OverlapMPI",  1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "OverlapMPI",  0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef OPENMP
   HDF5_Makefile->Add( "OpenMP",      1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "OpenMP",      0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef GPU
   HDF5_Makefile->Add( "UseGPU",      1,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "UseGPU",      0,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef GPU
   HDF5_Makefile->Add( "GPU_Arch",    GPU_ARCH,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_Makefile->Add( "GPU_Arch",    NULL_INT,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

} // FUNCTION : FillIn_HDF5_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_HDF5_SymConst
// Description :  Fill in the HDF5_Output_t structure with symbolic constants
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  HDF5_SymConst    : Pointer to be filled in
//                RS_FormatVersion : The format version of the restart data
//-------------------------------------------------------------------------------------------------------
void FillIn_HDF5_SymConst( HDF5_Output_t *HDF5_SymConst, const int RS_FormatVersion )
{

   const bool Compare_Yes = true, Compare_No = false;
   const bool   Fatal_Yes = true,   Fatal_No = false;


// model-independent variables
// ----------------------------------------------------------------------------------------------------
   HDF5_SymConst->Add( "NCompFluid",           NCOMP_FLUID,            0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
   HDF5_SymConst->Add( "NCompPassive",         NCOMP_PASSIVE,          0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
   HDF5_SymConst->Add( "PatchSize",            PS1,                    0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
   HDF5_SymConst->Add( "Flu_NIn",              FLU_NIN,                0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
   HDF5_SymConst->Add( "Flu_NOut",             FLU_NOUT,               0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "Flu_NIn_T",            FLU_NIN_T,              0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "Flu_NIn_S",            FLU_NIN_S,              0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "Flu_NOut_S",           FLU_NOUT_S,             0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "NFluxFluid",           NFLUX_FLUID,            0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "NFluxPassive",         NFLUX_PASSIVE,          0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "Flu_GhostSize",        FLU_GHOST_SIZE,         0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "Flu_Nxt",              FLU_NXT,                0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  ifdef DEBUG_HDF5
   HDF5_SymConst->Add( "Debug_HDF5",           1,                      0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  else
   HDF5_SymConst->Add( "Debug_HDF5",           0,                      0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif
   HDF5_SymConst->Add( "SibOffsetNonperiodic", SIB_OFFSET_NONPERIODIC, 0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  ifdef LOAD_BALANCE
   HDF5_SymConst->Add( "SonOffsetLB",          SON_OFFSET_LB,          0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif

   HDF5_SymConst->Add( "TinyNumber",           TINY_NUMBER,            0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "HugeNumber",           HUGE_NUMBER,            0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_SymConst->Add( "MaxError",             MAX_ERROR,              0, NULL, Compare_Yes, Fatal_No, Fatal_No  );


// gravity variables
// ----------------------------------------------------------------------------------------------------
#  ifdef GRAVITY
   HDF5_SymConst->Add( "Gra_NIn",          GRA_NIN,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Pot_GhostSize",    POT_GHOST_SIZE,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Gra_GhostSize",    GRA_GHOST_SIZE,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Rho_GhostSize",    RHO_GHOST_SIZE,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Pot_Nxt",          POT_NXT,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Gra_Nxt",          GRA_NXT,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Rho_Nxt",          RHO_NXT,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  ifdef UNSPLIT_GRAVITY
   HDF5_SymConst->Add( "USG_GhostSizeF",   USG_GHOST_SIZE_F,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "USG_GhostSizeG",   USG_GHOST_SIZE_G,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "USG_NxtF",         USG_NXT_F,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "USG_NxtG",         USG_NXT_G,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

   HDF5_SymConst->Add( "ExtPot_BlockSize", EXTPOT_BLOCK_SIZE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Gra_BlockSize",    GRA_BLOCK_SIZE,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "ExtPotNAuxMax",    EXT_POT_NAUX_MAX,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "ExtAccNAuxMax",    EXT_ACC_NAUX_MAX,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "ExtPotNGeneMax",   EXT_POT_NGENE_MAX, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  if   ( POT_SCHEME == SOR )
   HDF5_SymConst->Add( "Pot_BlockSize_z",  POT_BLOCK_SIZE_Z,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef SOR_RHO_SHARED
   HDF5_SymConst->Add( "SOR_RhoShared",    1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "SOR_RhoShared",    0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef SOR_CPOT_SHARED
   HDF5_SymConst->Add( "SOR_CPotShared",   1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "SOR_CPotShared",   0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef SOR_USE_SHUFFLE
   HDF5_SymConst->Add( "SOR_UseShuffle",   1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "SOR_UseShuffle",   0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef SOR_USE_PADDING
   HDF5_SymConst->Add( "SOR_UsePadding",   1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "SOR_UsePadding",   0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_SymConst->Add( "SOR_ModReduction", SOR_MOD_REDUCTION, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );

#  elif ( POT_SCHEME == MG  )
   HDF5_SymConst->Add( "Pot_BlockSize_x",  POT_BLOCK_SIZE_X,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif // POT_SCHEME
#  endif // #ifdef GRAVITY


// particle variables
// ----------------------------------------------------------------------------------------------------
#  ifdef PARTICLE
   const bool Fatal_ParNAttFlt = ( RS_FormatVersion > 2300 ) ? Fatal_Yes : Fatal_No;
   const bool Fatal_ParNAttInt = ( RS_FormatVersion > 2500 ) ? Fatal_Yes : Fatal_No;
   HDF5_SymConst->Add( "Par_NAttFltStored",    PAR_NATT_FLT_STORED,   0, NULL, Compare_Yes, Fatal_No, Fatal_ParNAttFlt );
   HDF5_SymConst->Add( "Par_NAttIntStored",    PAR_NATT_INT_STORED,   0, NULL, Compare_Yes, Fatal_No, Fatal_ParNAttInt );
   HDF5_SymConst->Add( "Par_NType",            PAR_NTYPE,             0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  ifdef GRAVITY
   HDF5_SymConst->Add( "RhoExt_GhostSize",     RHOEXT_GHOST_SIZE,     0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif

#  ifdef DEBUG_PARTICLE
   HDF5_SymConst->Add( "Debug_Particle",       1,                     0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  else
   HDF5_SymConst->Add( "Debug_Particle",       0,                     0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif

   HDF5_SymConst->Add( "ParList_GrowthFactor", PARLIST_GROWTH_FACTOR, 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
   HDF5_SymConst->Add( "ParList_ReduceFactor", PARLIST_REDUCE_FACTOR, 0, NULL, Compare_Yes, Fatal_No, Fatal_No         );
#  endif // #ifdef PARTICLE


// bitwise reproducibility variables
// ----------------------------------------------------------------------------------------------------
#  ifdef BIT_REP_FLUX
   HDF5_SymConst->Add( "BitRep_Flux",     1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "BitRep_Flux",     0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef MHD
#  ifdef BIT_REP_ELECTRIC
   HDF5_SymConst->Add( "BitRep_Electric", 1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "BitRep_Electric", 0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // #ifdef MHD


// HYDRO variables
// ----------------------------------------------------------------------------------------------------
#  if   ( MODEL == HYDRO )
   HDF5_SymConst->Add( "Flu_BlockSize_x",     FLU_BLOCK_SIZE_X,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Flu_BlockSize_y",     FLU_BLOCK_SIZE_Y,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   HDF5_SymConst->Add( "CheckUnphyInFluid",   1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "CheckUnphyInFluid",   0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef CHAR_RECONSTRUCTION
   HDF5_SymConst->Add( "CharReconstruction",  1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "CharReconstruction",  0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef LR_EINT
   HDF5_SymConst->Add( "LR_Eint",             1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "LR_Eint",             0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef CHECK_INTERMEDIATE
   HDF5_SymConst->Add( "CheckIntermediate",   CHECK_INTERMEDIATE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "CheckIntermediate",   0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef RSOLVER_RESCUE
   HDF5_SymConst->Add( "RSolverRescue",       RSOLVER_RESCUE,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "RSolverRescue",       0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef HLL_NO_REF_STATE
   HDF5_SymConst->Add( "HLL_NoRefState",      1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "HLL_NoRefState",      0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef HLL_INCLUDE_ALL_WAVES
   HDF5_SymConst->Add( "HLL_IncludeAllWaves", 1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "HLL_IncludeAllWaves", 0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_SymConst->Add( "HLLC_WaveSpeed",      HLLC_WAVESPEED,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "HLLE_WaveSpeed",      HLLE_WAVESPEED,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef MHD
   HDF5_SymConst->Add( "HLLD_WaveSpeed",      HLLD_WAVESPEED,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef N_FC_VAR
   HDF5_SymConst->Add( "N_FC_Var",            N_FC_VAR,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef N_SLOPE_PPM
   HDF5_SymConst->Add( "N_Slope_PPM",         N_SLOPE_PPM,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef MHD
#  ifdef EULERY
   HDF5_SymConst->Add( "EulerY",              1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "EulerY",              0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // MHD
#  ifdef MHM_CHECK_PREDICT
   HDF5_SymConst->Add( "MHM_CheckPredict",    1,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "MHM_CheckPredict",    0,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_SymConst->Add( "EoSNAuxMax",          EOS_NAUX_MAX,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "EoSNTableMax",        EOS_NTABLE_MAX,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// ELBDM variables
// ----------------------------------------------------------------------------------------------------
#  elif  ( MODEL == ELBDM )
   HDF5_SymConst->Add( "Flu_BlockSize_x",    FLU_BLOCK_SIZE_X,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Flu_BlockSize_y",    FLU_BLOCK_SIZE_Y,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   HDF5_SymConst->Add( "Flu_HJ_BlockSize_y", FLU_HJ_BLOCK_SIZE_Y, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  if ( WAVE_SCHEME == WAVE_GRAMFE )
   HDF5_SymConst->Add( "GramFEScheme",       GRAMFE_SCHEME,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "GramFEGamma",        GRAMFE_GAMMA,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "GramFEG",            GRAMFE_G,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "GramFENDelta",       GRAMFE_NDELTA,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "GramFEOrder",        GRAMFE_ORDER,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "GramFEND",           GRAMFE_ND,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "GramFEFluNxt",       GRAMFE_FLU_NXT,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// dt variables
// ----------------------------------------------------------------------------------------------------
   HDF5_SymConst->Add( "dt_Flu_BlockSize",  DT_FLU_BLOCK_SIZE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef DT_FLU_USE_SHUFFLE
   HDF5_SymConst->Add( "dt_Flu_UseShuffle", 1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "dt_Flu_UseShuffle", 0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef GRAVITY
   HDF5_SymConst->Add( "dt_Gra_BlockSize",  DT_GRA_BLOCK_SIZE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef DT_GRA_USE_SHUFFLE
   HDF5_SymConst->Add( "dt_Gra_UseShuffle", 1,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "dt_Gra_UseShuffle", 0,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif


// source term variables
// ----------------------------------------------------------------------------------------------------
   HDF5_SymConst->Add( "Src_BlockSize",       SRC_BLOCK_SIZE,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Src_GhostSize",       SRC_GHOST_SIZE,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Src_Nxt",             SRC_NXT,               0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == HYDRO )
   HDF5_SymConst->Add( "Src_NAuxDlep",        SRC_NAUX_DLEP,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Src_DlepProfNVar",    SRC_DLEP_PROF_NVAR,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Src_DlepProfNBinMax", SRC_DLEP_PROF_NBINMAX, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_SymConst->Add( "Src_NAuxUser",        SRC_NAUX_USER,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );

   HDF5_SymConst->Add( "Der_GhostSize",       DER_GHOST_SIZE,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Der_Nxt",             DER_NXT,               0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "Der_NOut_Max",        DER_NOUT_MAX,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// feedback variables
// ----------------------------------------------------------------------------------------------------
#  ifdef FEEDBACK
   HDF5_SymConst->Add( "FB_GhostSize", FB_GHOST_SIZE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_SymConst->Add( "FB_Nxt",       FB_NXT,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif

#  ifdef FB_SEP_FLUOUT
   HDF5_SymConst->Add( "FB_SepFluOut", 1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "FB_SepFluOut", 0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// interpolate variables
// ----------------------------------------------------------------------------------------------------
#  ifdef INTERP_MASK
   HDF5_SymConst->Add( "InterpMask",   1, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  else
   HDF5_SymConst->Add( "InterpMask",   0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// maximum number variables
// ----------------------------------------------------------------------------------------------------
   HDF5_SymConst->Add( "NFieldStoredMax", NFIELD_STORED_MAX, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );

   HDF5_SymConst->Add( "NConRefMax",      NCONREF_MAX,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );

} // FUNCTION : FillIn_HDF5_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_HDF5_InputPara
// Description :  Fill in the HDF5_Output_t structure with input parameters
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  HDF5_InputPara : Pointer to be filled in
//                NFieldStored   : Number of grid fields to be stored on disk
//                FieldLabelOut  : Field labels
//                Load_RS        : Whether the structure is used for checking the restart data
//-------------------------------------------------------------------------------------------------------
void FillIn_HDF5_InputPara( HDF5_Output_t *HDF5_InputPara, const int NFieldStored,
                            char FieldLabelOut[][MAX_STRING], const bool Load_RS )
{

   const bool Compare_Yes = true, Compare_No = false;
   const bool   Fatal_Yes = true,   Fatal_No = false;

   char Key[MAX_STRING];


// simulation scale
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "BoxSize",     &BOX_SIZE,       0, NULL,     Compare_Yes, Fatal_Yes, Fatal_No );
   HDF5_InputPara->Add( "NX0_Tot",     &NX0_TOT[0],     1, H5_Arr_3, Compare_Yes, Fatal_Yes, Fatal_No );
   HDF5_InputPara->Add( "MPI_NRank",   &MPI_NRank,      0, NULL,     Compare_Yes, Fatal_No,  Fatal_No );
   HDF5_InputPara->Add( "MPI_NRank_X", &MPI_NRank_X[0], 1, H5_Arr_3, Compare_Yes, Fatal_No,  Fatal_No );
   HDF5_InputPara->Add( "OMP_NThread", &OMP_NTHREAD,    0, NULL,     Compare_Yes, Fatal_No,  Fatal_No );
   HDF5_InputPara->Add( "EndT",        &END_T,          0, NULL,     Compare_Yes, Fatal_No,  Fatal_No );
   HDF5_InputPara->Add( "EndStep",     &END_STEP,       0, NULL,     Compare_Yes, Fatal_No,  Fatal_No );


// test problems
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "TestProb_ID", &TESTPROB_ID, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// code units
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Opt__Unit",  &OPT__UNIT, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_L",     &UNIT_L,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_M",     &UNIT_M,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_T",     &UNIT_T,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_V",     &UNIT_V,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_D",     &UNIT_D,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_E",     &UNIT_E,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Unit_P",     &UNIT_P,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef MHD
   HDF5_InputPara->Add( "Unit_B",     &UNIT_B,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// boundary condition
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Opt__BC_Flu",  &OPT__BC_FLU[0], 1, H5_Arr_6, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef GRAVITY
   HDF5_InputPara->Add( "Opt__BC_Pot",  &OPT__BC_POT,    0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "GFunc_Coeff0", &GFUNC_COEFF0,   0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
#  endif


// particle
// ----------------------------------------------------------------------------------------------------
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Par_Init",            &amr->Par->Init,            0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_ICFormat",        &amr->Par->ParICFormat,     0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_ICMass",          &amr->Par->ParICMass,       0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_ICType",          &amr->Par->ParICType,       0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_ICFloat8",        &PAR_IC_FLOAT8,             0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_ICInt8",          &PAR_IC_INT8,               0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_Interp",          &amr->Par->Interp,          0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_InterpTracer",    &amr->Par->InterpTracer,    0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_Integ",           &amr->Par->Integ,           0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_IntegTracer",     &amr->Par->IntegTracer,     0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_ImproveAcc",      &amr->Par->ImproveAcc,      0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_PredictPos",      &amr->Par->PredictPos,      0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_RemoveCell",      &amr->Par->RemoveCell,      0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__FreezePar",      &OPT__FREEZE_PAR,           0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_GhostSize",       &amr->Par->GhostSize,       0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_GhostSizeTracer", &amr->Par->GhostSizeTracer, 0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Par_TracerVelCorr",   &amr->Par->TracerVelCorr,   0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__ParInitCheck",   &OPT__PAR_INIT_CHECK,       0, NULL, Compare_Yes, Fatal_No,  Fatal_No  );
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
   {
   sprintf( Key, "ParAttFltLabel%02d", v );
   HDF5_InputPara->Add( Key,                   &ParAttFltLabel[v][0],      0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   }
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
   {
   sprintf( Key, "ParAttIntLabel%02d", v );
   HDF5_InputPara->Add( Key,                   &ParAttIntLabel[v][0],      0, NULL, Compare_No,  NULL_BOOL, NULL_BOOL );
   }
#  endif


// cosmology
// ----------------------------------------------------------------------------------------------------
#  ifdef COMOVING
   HDF5_InputPara->Add( "A_Init",  &A_INIT,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "OmegaM0", &OMEGA_M0, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Hubble0", &HUBBLE0,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// time-step determination
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Dt__Max",                 &DT__MAX,                     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__Fluid",               &DT__FLUID,                   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__FluidInit",           &DT__FLUID_INIT,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef GRAVITY
   HDF5_InputPara->Add( "Dt__Gravity",             &DT__GRAVITY,                 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  if ( MODEL == ELBDM )
   HDF5_InputPara->Add( "Dt__Phase",               &DT__PHASE,                   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   HDF5_InputPara->Add( "Dt__HybridCFL",           &DT__HYBRID_CFL,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__HybridCFLInit",       &DT__HYBRID_CFL_INIT,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__HybridVelocity",      &DT__HYBRID_VELOCITY,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__HybridVelocityInit",  &DT__HYBRID_VELOCITY_INIT,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // #if ( MODEL == ELBDM )
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Dt__ParVel",              &DT__PARVEL,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__ParVelMax",           &DT__PARVEL_MAX,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__ParAcc",              &DT__PARACC,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef CR_DIFFUSION
   HDF5_InputPara->Add( "Dt__CR_Diffusion",        &DT__CR_DIFFUSION,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef COMOVING
   HDF5_InputPara->Add( "Dt__MaxDeltaA",           &DT__MAX_DELTA_A,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Dt__SyncParentLv",        &DT__SYNC_PARENT_LV,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Dt__SyncChildrenLv",      &DT__SYNC_CHILDREN_LV,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__DtUser",             &OPT__DT_USER,                0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__DtLevel",            &OPT__DT_LEVEL,               0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordDt",           &OPT__RECORD_DT,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AutoReduceDt",            &AUTO_REDUCE_DT,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AutoReduceDtFactor",      &AUTO_REDUCE_DT_FACTOR,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AutoReduceDtFactorMin",   &AUTO_REDUCE_DT_FACTOR_MIN,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "AutoReduceMinModFactor",  &AUTO_REDUCE_MINMOD_FACTOR,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AutoReduceMinModMin",     &AUTO_REDUCE_MINMOD_MIN,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "AutoReduceIntMonoFactor", &AUTO_REDUCE_INT_MONO_FACTOR, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AutoReduceIntMonoMin",    &AUTO_REDUCE_INT_MONO_MIN,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// domain refinement
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "RegridCount",             &REGRID_COUNT,               0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "RefineNLevel",            &REFINE_NLEVEL,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagBufferSize",          &FLAG_BUFFER_SIZE,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagBufferSizeMaxM1Lv",   &FLAG_BUFFER_SIZE_MAXM1_LV,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagBufferSizeMaxM2Lv",   &FLAG_BUFFER_SIZE_MAXM2_LV,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MaxLevel",                &MAX_LEVEL,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Rho",           &OPT__FLAG_RHO,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_RhoGradient",   &OPT__FLAG_RHO_GRADIENT,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Angular",       &OPT__FLAG_ANGULAR,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Radial",        &OPT__FLAG_RADIAL,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Opt__Flag_PresGradient",  &OPT__FLAG_PRES_GRADIENT,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Vorticity",     &OPT__FLAG_VORTICITY,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Jeans",         &OPT__FLAG_JEANS,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__Flag_Current",       &OPT__FLAG_CURRENT,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef SRHD
   HDF5_InputPara->Add( "Dt__SpeedOfLight",        &DT__SPEED_OF_LIGHT,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_LrtzGradient",  &OPT__FLAG_LRTZ_GRADIENT,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef COSMIC_RAY
   HDF5_InputPara->Add( "Opt__Flag_CRay",          &OPT__FLAG_CRAY,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif
#  if ( MODEL == ELBDM )
   HDF5_InputPara->Add( "Opt__Flag_EngyDensity",   &OPT__FLAG_ENGY_DENSITY,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Spectral",      &OPT__FLAG_SPECTRAL,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Spectral_N",    &OPT__FLAG_SPECTRAL_N,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   HDF5_InputPara->Add( "Opt__Flag_Interference",  &OPT__FLAG_INTERFERENCE,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // #if ( MODEL == ELBDM )
   HDF5_InputPara->Add( "Opt__Flag_LohnerDens",    &OPT__FLAG_LOHNER_DENS,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Opt__Flag_LohnerEngy",    &OPT__FLAG_LOHNER_ENGY,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_LohnerPres",    &OPT__FLAG_LOHNER_PRES,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_LohnerTemp",    &OPT__FLAG_LOHNER_TEMP,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_LohnerEntr",    &OPT__FLAG_LOHNER_ENTR,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef COSMIC_RAY
   HDF5_InputPara->Add( "Opt__Flag_LohnerCRay",    &OPT__FLAG_LOHNER_CRAY,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif
   HDF5_InputPara->Add( "Opt__Flag_LohnerForm",    &OPT__FLAG_LOHNER_FORM,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_User",          &OPT__FLAG_USER,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_User_Num",      &OPT__FLAG_USER_NUM,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_Region",        &OPT__FLAG_REGION,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Opt__Flag_NParPatch",     &OPT__FLAG_NPAR_PATCH,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_NParCell",      &OPT__FLAG_NPAR_CELL,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Flag_ParMassCell",   &OPT__FLAG_PAR_MASS_CELL,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__NoFlagNearBoundary", &OPT__NO_FLAG_NEAR_BOUNDARY, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__PatchCount",         &OPT__PATCH_COUNT,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Opt__ParticleCount",      &OPT__PARTICLE_COUNT,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__ReuseMemory",        &OPT__REUSE_MEMORY,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__MemoryPool",         &OPT__MEMORY_POOL,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// load balance
// ----------------------------------------------------------------------------------------------------
#  ifdef LOAD_BALANCE
   HDF5_InputPara->Add( "LB_WLI_Max",              &amr->LB->WLI_Max,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef PARTICLE
   HDF5_InputPara->Add( "LB_Par_Weight",           &amr->LB->Par_Weight,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__RecordLoadBalance",  &OPT__RECORD_LOAD_BALANCE,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__LB_ExchangeFather",  &OPT__LB_EXCHANGE_FATHER,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__MinimizeMPIBarrier", &OPT__MINIMIZE_MPI_BARRIER, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// fluid solvers in HYDRO
// ----------------------------------------------------------------------------------------------------
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Gamma",                  &GAMMA,                     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MolecularWeight",        &MOLECULAR_WEIGHT,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MuNorm",                 &MU_NORM,                   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "IsoTemp",                &ISO_TEMP,                  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MinMod_Coeff",           &MINMOD_COEFF,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MinMod_MaxIter",         &MINMOD_MAX_ITER,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__LR_Limiter",        &OPT__LR_LIMITER,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__1stFluxCorr",       &OPT__1ST_FLUX_CORR,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__1stFluxCorrScheme", &OPT__1ST_FLUX_CORR_SCHEME, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef DUAL_ENERGY
   HDF5_InputPara->Add( "DualEnergySwitch",       &DUAL_ENERGY_SWITCH,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__SameInterfaceB",    &OPT__SAME_INTERFACE_B,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // HYDRO


// ELBDM solvers
// ----------------------------------------------------------------------------------------------------
#  if ( MODEL == ELBDM )
   HDF5_InputPara->Add( "ELBDM_Mass",           &ELBDM_MASS,             0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_InputPara->Add( "ELBDM_PlanckConst",    &ELBDM_PLANCK_CONST,     0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  ifdef QUARTIC_SELF_INTERACTION
   HDF5_InputPara->Add( "ELBDM_Lambda",         &ELBDM_LAMBDA,           0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  endif
   HDF5_InputPara->Add( "ELBDM_Taylor3_Coeff",  &ELBDM_TAYLOR3_COEFF,    0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_InputPara->Add( "ELBDM_Taylor3_Auto",   &ELBDM_TAYLOR3_AUTO,     0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_InputPara->Add( "ELBDM_RemoveMotionCM", &ELBDM_REMOVE_MOTION_CM, 0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
   HDF5_InputPara->Add( "ELBDM_BaseSpectral",   &ELBDM_BASE_SPECTRAL,    0, NULL, Compare_Yes, Fatal_No, Fatal_No  );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
// ELBDM_FIRST_WAVE_LEVEL currently cannot be changed upon restart because the code cannot robustly handle the conversion
// from Re/Im to Dens/Phase due to the phase ambiguity introduced by vortices
   HDF5_InputPara->Add( "ELBDM_FirstWaveLevel", &ELBDM_FIRST_WAVE_LEVEL, 0, NULL, Compare_Yes, Fatal_No, Fatal_Yes );
#  endif
#  endif // ELBDM


// fluid solvers in different models
// ----------------------------------------------------------------------------------------------------
#  if ( NCOMP_PASSIVE > 0 )
   int H5_Arr_NP[1] = { NCOMP_PASSIVE };
#  endif

   HDF5_InputPara->Add( "Flu_GPU_NPGroup",         &FLU_GPU_NPGROUP,           0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "GPU_NStream",             &GPU_NSTREAM,               0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__FixUp_Flux",         &OPT__FIXUP_FLUX,           0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "FixUpFlux_Var",           &FixUpVar_Flux,             0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__FixUp_Electric",     &OPT__FIXUP_ELECTRIC,       0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  endif
   HDF5_InputPara->Add( "Opt__FixUp_Restrict",     &OPT__FIXUP_RESTRICT,       0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "FixUpRestrict_Var",       &FixUpVar_Restrict,         0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__CorrAfterAllSync",   &OPT__CORR_AFTER_ALL_SYNC,  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__NormalizePassive",   &OPT__NORMALIZE_PASSIVE,    0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "NormalizePassive_NVar",   &PassiveNorm_NVar,          0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  if ( NCOMP_PASSIVE > 0 )
   HDF5_InputPara->Add( "NormalizePassive_VarIdx", &PassiveNorm_VarIdx[0],     1, H5_Arr_NP, Compare_Yes, Fatal_No,  Fatal_No  );
#  endif
   HDF5_InputPara->Add( "Opt__IntFracPassive_LR",  &OPT__INT_FRAC_PASSIVE_LR,  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "IntFracPassive_NVar",     &PassiveIntFrac_NVar,       0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  if ( NCOMP_PASSIVE > 0 )
   HDF5_InputPara->Add( "IntFracPassive_VarIdx",   &PassiveIntFrac_VarIdx[0],  1, H5_Arr_NP, Compare_Yes, Fatal_No,  Fatal_No  );
#  endif

   for (int v=0; v<NFieldStored; v++)
   {
   sprintf( Key, "FieldLabel%02d", v );
   HDF5_InputPara->Add( Key,                       &FieldLabelOut[v][0],       0, NULL,      Compare_No,  NULL_BOOL, NULL_BOOL );
   }

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
   sprintf( Key, "MagLabel%02d", v );
   HDF5_InputPara->Add( Key,                       &MagLabel[v][0],            0, NULL,      Compare_No,  NULL_BOOL, NULL_BOOL );
   }
#  endif

   HDF5_InputPara->Add( "Opt__OverlapMPI",         &OPT__OVERLAP_MPI,          0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__ResetFluid",         &OPT__RESET_FLUID,          0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__ResetFluidInit",     &OPT__RESET_FLUID_INIT,     0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__FreezeFluid",        &OPT__FREEZE_FLUID,         0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   HDF5_InputPara->Add( "MinDens",                 &MIN_DENS,                  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  endif
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "MinPres",                 &MIN_PRES,                  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "MinEint",                 &MIN_EINT,                  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "MinTemp",                 &MIN_TEMP,                  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "MinEntr",                 &MIN_ENTR,                  0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__CheckPresAfterFlu",  &OPT__CHECK_PRES_AFTER_FLU, 0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__LastResortFloor",    &OPT__LAST_RESORT_FLOOR,    0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "JeansMinPres",            &JEANS_MIN_PRES,            0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "JeansMinPres_Level",      &JEANS_MIN_PRES_LEVEL,      0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "JeansMinPres_NCell",      &JEANS_MIN_PRES_NCELL,      0, NULL,      Compare_Yes, Fatal_No,  Fatal_No  );
#  endif


// self-gravity
// ----------------------------------------------------------------------------------------------------
#  ifdef GRAVITY
   HDF5_InputPara->Add( "NewtonG",               &NEWTON_G,                0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
#  if   ( POT_SCHEME == SOR )
   HDF5_InputPara->Add( "SOR_Omega",             &SOR_OMEGA,               0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SOR_MaxIter",           &SOR_MAX_ITER,            0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SOR_MinIter",           &SOR_MIN_ITER,            0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
#  elif ( POT_SCHEME == MG )
   HDF5_InputPara->Add( "MG_MaxIter",            &MG_MAX_ITER,             0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MG_NPreSmooth",         &MG_NPRE_SMOOTH,          0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MG_NPostSmooth",        &MG_NPOST_SMOOTH,         0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "MG_ToleratedError",     &MG_TOLERATED_ERROR,      0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Pot_GPU_NPGroup",       &POT_GPU_NPGROUP,         0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__GraP5Gradient",    &OPT__GRA_P5_GRADIENT,    0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__SelfGravity",      &OPT__SELF_GRAVITY,       0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__ExtAcc",           &OPT__EXT_ACC,            0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__ExtPot",           &OPT__EXT_POT,            0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "ExtPotTable_Name",       EXT_POT_TABLE_NAME,      0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "ExtPotTable_NPoint",    &EXT_POT_TABLE_NPOINT[0], 1, H5_Arr_3, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "ExtPotTable_dh",        &EXT_POT_TABLE_DH[0],     1, H5_Arr_3, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "ExtPotTable_EdgeL",     &EXT_POT_TABLE_EDGEL[0],  1, H5_Arr_3, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "ExtPotTable_Float8",    &EXT_POT_TABLE_FLOAT8,    0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__GravityExtraMass", &OPT__GRAVITY_EXTRA_MASS, 0, NULL,     Compare_Yes, Fatal_No, Fatal_No );
#  endif


// source terms
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Src_Deleptonization", &SrcTerms.Deleptonization, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Src_User",            &SrcTerms.User,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Src_GPU_NPGroup",     &SRC_GPU_NPGROUP,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// Grackle
// ----------------------------------------------------------------------------------------------------
#  ifdef SUPPORT_GRACKLE
   HDF5_InputPara->Add( "Grackle_Activate",       &GRACKLE_ACTIVATE,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_Verbose",        &GRACKLE_VERBOSE,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_Cooling",        &GRACKLE_COOLING,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_Primordial",     &GRACKLE_PRIMORDIAL,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_Metal",          &GRACKLE_METAL,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_UV",             &GRACKLE_UV,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_CMB_Floor",      &GRACKLE_CMB_FLOOR,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_PE_Heating",     &GRACKLE_PE_HEATING,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_PE_HeatingRate", &GRACKLE_PE_HEATING_RATE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_CloudyTable",    &GRACKLE_CLOUDY_TABLE,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_ThreeBodyRate",  &GRACKLE_THREE_BODY_RATE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_CIE_Cooling",    &GRACKLE_CIE_COOLING,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Grackle_H2_OpaApprox",   &GRACKLE_H2_OPA_APPROX,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Che_GPU_NPGroup",        &CHE_GPU_NPGROUP,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// star formation
// ----------------------------------------------------------------------------------------------------
#  ifdef STAR_FORMATION
   HDF5_InputPara->Add( "SF_CreateStar_Scheme",       &SF_CREATE_STAR_SCHEME,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_RSeed",        &SF_CREATE_STAR_RSEED,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_DetRandom",    &SF_CREATE_STAR_DET_RANDOM,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_MinLevel",     &SF_CREATE_STAR_MIN_LEVEL,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_MinGasDens",   &SF_CREATE_STAR_MIN_GAS_DENS,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_MassEff",      &SF_CREATE_STAR_MASS_EFF,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_MinStarMass",  &SF_CREATE_STAR_MIN_STAR_MASS,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SF_CreateStar_MaxStarMFrac", &SF_CREATE_STAR_MAX_STAR_MFRAC, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// feedback
// ----------------------------------------------------------------------------------------------------
#  ifdef FEEDBACK
   HDF5_InputPara->Add( "FB_Level", &FB_LEVEL, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FB_RSeed", &FB_RSEED, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FB_SNe",   &FB_SNE,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FB_User",  &FB_USER,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// cosmic ray
// ----------------------------------------------------------------------------------------------------
#  ifdef COSMIC_RAY
   HDF5_InputPara->Add( "CR_Gamma",               &GAMMA_CR,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef CR_DIFFUSION
   HDF5_InputPara->Add( "CR_Diffusion_ParaCoeff", &CR_DIFF_PARA,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "CR_Diffusion_PerpCoeff", &CR_DIFF_PERP,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "CR_Diffusion_MinB",      &CR_DIFF_MIN_B, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // #ifdef COSMIC_RAY


// initialization
// ----------------------------------------------------------------------------------------------------
   int H5_Arr_NLv_M1_6[2] = { NLEVEL-1, 6 };
   int RefineRegion[NLEVEL-1][6];
   for (int t=0; t<NLEVEL-1; t++)   for (int s=0; s<6; s++)   RefineRegion[t][s] = -1;
   if ( OPT__INIT == INIT_BY_FILE  &&  OPT__UM_IC_NLEVEL > 1  &&  UM_IC_RefineRegion != NULL )
   {
      const int (*temp)[6] = ( int(*)[6] )UM_IC_RefineRegion;
      for (int t=0; t<NLEVEL-1; t++)
         if ( t < OPT__UM_IC_NLEVEL - 1 )
            for (int s=0; s<6; s++)
               RefineRegion[t][s] = temp[t][s];
   }

   HDF5_InputPara->Add( "Opt__Init",               &OPT__INIT,                 0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "RestartLoadNRank",        &RESTART_LOAD_NRANK,        0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__RestartReset",       &OPT__RESTART_RESET,        0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_Level",        &OPT__UM_IC_LEVEL,          0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_NLevel",       &OPT__UM_IC_NLEVEL,         0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_NVar",         &OPT__UM_IC_NVAR,           0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_Format",       &OPT__UM_IC_FORMAT,         0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_Float8",       &OPT__UM_IC_FLOAT8,         0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_Downgrade",    &OPT__UM_IC_DOWNGRADE,      0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_Refine",       &OPT__UM_IC_REFINE,         0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__UM_IC_LoadNRank",    &OPT__UM_IC_LOAD_NRANK,     0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "UM_IC_RefineRegion",      &RefineRegion[0][0],        2, H5_Arr_NLv_M1_6, Compare_No,  NULL_BOOL, NULL_BOOL );
   HDF5_InputPara->Add( "Opt__InitRestrict",       &OPT__INIT_RESTRICT,        0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__InitGridWithOMP",    &OPT__INIT_GRID_WITH_OMP,   0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Opt__GPUID_Select",       &OPT__GPUID_SELECT,         0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
   HDF5_InputPara->Add( "Init_Subsampling_NCell",  &INIT_SUBSAMPLING_NCELL,    0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__InitBFieldByVecPot", &OPT__INIT_BFIELD_BYVECPOT, 0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
#  endif
#  ifdef SUPPORT_FFTW
   HDF5_InputPara->Add( "Opt__FFTW_Startup",       &OPT__FFTW_STARTUP,         0, NULL,            Compare_Yes, Fatal_No,  Fatal_No  );
#  endif


// interpolation schemes
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Opt__Int_Time",            &OPT__INT_TIME,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Opt__Int_Prim",            &OPT__INT_PRIM,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  if ( MODEL == ELBDM )
   HDF5_InputPara->Add( "Opt__Int_Phase",           &OPT__INT_PHASE,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Res_Phase",           &OPT__RES_PHASE,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   HDF5_InputPara->Add( "Opt__Hybrid_Match_Phase",  &ELBDM_MATCH_PHASE,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // ELBDM
   HDF5_InputPara->Add( "Opt__Flu_IntScheme",       &OPT__FLU_INT_SCHEME,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RefFlu_IntScheme",    &OPT__REF_FLU_INT_SCHEME,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__Mag_IntScheme",       &OPT__MAG_INT_SCHEME,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RefMag_IntScheme",    &OPT__REF_MAG_INT_SCHEME,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef GRAVITY
   HDF5_InputPara->Add( "Opt__Pot_IntScheme",       &OPT__POT_INT_SCHEME,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Rho_IntScheme",       &OPT__RHO_INT_SCHEME,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Gra_IntScheme",       &OPT__GRA_INT_SCHEME,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RefPot_IntScheme",    &OPT__REF_POT_INT_SCHEME,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "IntMonoCoeff",             &INT_MONO_COEFF,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef MHD
   HDF5_InputPara->Add( "IntMonoCoeffB",            &INT_MONO_COEFF_B,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Mono_MaxIter",             &MONO_MAX_ITER,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "IntOppSign0thOrder",       &INT_OPP_SIGN_0TH_ORDER,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef SUPPORT_SPECTRAL_INT
   HDF5_InputPara->Add( "SpecInt_TablePath",         SPEC_INT_TABLE_PATH,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SpecInt_GhostBoundary",    &SPEC_INT_GHOST_BOUNDARY,   0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == ELBDM )
   HDF5_InputPara->Add( "SpecInt_XY_Instead_DePha", &SPEC_INT_XY_INSTEAD_DEPHA, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "SpecInt_VortexThreshold",  &SPEC_INT_VORTEX_THRESHOLD, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  endif // #ifdef SUPPORT_SPECTRAL_INT


// data dump
// ----------------------------------------------------------------------------------------------------
   const bool Compare_Part = ( OPT__OUTPUT_PART ) ? Compare_Yes : Compare_No;
#  ifdef PARTICLE
   const bool Compare_Out  = ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PAR_MODE ) ? Compare_Yes : Compare_No;
#  else
   const bool Compare_Out  = ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS ) ? Compare_Yes : Compare_No;
#  endif
   HDF5_InputPara->Add( "Opt__Output_Total",           &OPT__OUTPUT_TOTAL,           0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Part",            &OPT__OUTPUT_PART,            0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_User",            &OPT__OUTPUT_USER,            0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Opt__Output_Par_Mode",        &OPT__OUTPUT_PAR_MODE,        0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Par_Mesh",        &OPT__OUTPUT_PAR_MESH,        0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__Output_BasePS",          &OPT__OUTPUT_BASEPS,          0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Base",            &OPT__OUTPUT_BASE,            0, NULL, Compare_Part, Fatal_No, Fatal_No );
#  ifdef GRAVITY
   HDF5_InputPara->Add( "Opt__Output_Pot",             &OPT__OUTPUT_POT,             0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  endif
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Opt__Output_ParDens",         &OPT__OUTPUT_PAR_DENS,        0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  endif
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__Output_CC_Mag",          &OPT__OUTPUT_CC_MAG,          0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  endif
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Opt__Output_Pres",            &OPT__OUTPUT_PRES,            0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Temp",            &OPT__OUTPUT_TEMP,            0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Entr",            &OPT__OUTPUT_ENTR,            0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Cs",              &OPT__OUTPUT_CS,              0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_DivVel",          &OPT__OUTPUT_DIVVEL,          0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Mach",            &OPT__OUTPUT_MACH,            0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__Output_DivMag",          &OPT__OUTPUT_DIVMAG,          0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  endif
#  ifdef SRHD
   HDF5_InputPara->Add( "Opt__Output_3Velocity",       &OPT__OUTPUT_3VELOCITY,       0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Lorentz",         &OPT__OUTPUT_LORENTZ,         0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Enthalpy",        &OPT__OUTPUT_ENTHALPY,        0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
#  endif
#  endif // #if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Opt__Output_UserField",       &OPT__OUTPUT_USER_FIELD,      0, NULL, Compare_Yes,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Mode",            &OPT__OUTPUT_MODE,            0, NULL, Compare_Out,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Restart",         &OPT__OUTPUT_RESTART,         0, NULL, Compare_Out,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Step",            &OUTPUT_STEP,                 0, NULL, Compare_Out,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Dt",              &OUTPUT_DT,                   0, NULL, Compare_Out,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Text_Format_Flt",  OPT__OUTPUT_TEXT_FORMAT_FLT, 0, NULL, Compare_Out,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Output_Text_Length_Int", &OPT__OUTPUT_TEXT_LENGTH_INT, 0, NULL, Compare_Out,  Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Output_PartX",                &OUTPUT_PART_X,               0, NULL, Compare_Part, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Output_PartY",                &OUTPUT_PART_Y,               0, NULL, Compare_Part, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Output_PartZ",                &OUTPUT_PART_Z,               0, NULL, Compare_Part, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "InitDumpID",                  &INIT_DUMPID,                 0, NULL, Compare_Yes,  Fatal_No, Fatal_No );


// libyt jupyter
// ----------------------------------------------------------------------------------------------------
#  if (  defined( SUPPORT_LIBYT )  &&  defined( LIBYT_JUPYTER )  )
   HDF5_InputPara->Add( "Yt_JupyterUseConnectionFile", &YT_JUPYTER_USE_CONNECTION_FILE,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif


// miscellaneous
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Opt__Verbose",            &OPT__VERBOSE,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__TimingBarrier",      &OPT__TIMING_BARRIER,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__TimingBalance",      &OPT__TIMING_BALANCE,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__TimingMPI",          &OPT__TIMING_MPI,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordNote",         &OPT__RECORD_NOTE,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordUnphy",        &OPT__RECORD_UNPHY,        0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordMemory",       &OPT__RECORD_MEMORY,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordPerformance",  &OPT__RECORD_PERFORMANCE,  0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__ManualControl",      &OPT__MANUAL_CONTROL,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordCenter",       &OPT__RECORD_CENTER,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_CenX",                &COM_CEN_X,                0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_CenY",                &COM_CEN_Y,                0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_CenZ",                &COM_CEN_Z,                0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_MaxR",                &COM_MAX_R,                0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_MinRho",              &COM_MIN_RHO,              0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_TolErrR",             &COM_TOLERR_R,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "COM_MaxIter",             &COM_MAX_ITER,             0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__RecordUser",         &OPT__RECORD_USER,         0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__OptimizeAggressive", &OPT__OPTIMIZE_AGGRESSIVE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__SortPatchByLBIdx",   &OPT__SORT_PATCH_BY_LBIDX, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// simulation checks
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "Opt__Ck_Refine",        &OPT__CK_REFINE,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_ProperNesting", &OPT__CK_PROPER_NESTING,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_Conservation",  &OPT__CK_CONSERVATION,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AngMom_OriginX",        &ANGMOM_ORIGIN_X,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AngMom_OriginY",        &ANGMOM_ORIGIN_Y,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "AngMom_OriginZ",        &ANGMOM_ORIGIN_Z,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_NormPassive",   &OPT__CK_NORMALIZE_PASSIVE, 0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_Restrict",      &OPT__CK_RESTRICT,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_Finite",        &OPT__CK_FINITE,            0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_PatchAllocate", &OPT__CK_PATCH_ALLOCATE,    0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_FluxAllocate",  &OPT__CK_FLUX_ALLOCATE,     0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  if ( MODEL == HYDRO )
   HDF5_InputPara->Add( "Opt__Ck_Negative",      &OPT__CK_NEGATIVE,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__Ck_MemFree",       &OPT__CK_MEMFREE,           0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  ifdef PARTICLE
   HDF5_InputPara->Add( "Opt__Ck_Particle",      &OPT__CK_PARTICLE,          0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
#  ifdef MHD
   HDF5_InputPara->Add( "Opt__Ck_InterfaceB",    &OPT__CK_INTERFACE_B,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "Opt__Ck_DivergenceB",   &OPT__CK_DIVERGENCE_B,      0, NULL, Compare_Yes, Fatal_No, Fatal_No );
#  endif
   HDF5_InputPara->Add( "Opt__Ck_InputFluid",    &OPT__CK_INPUT_FLUID,       0, NULL, Compare_Yes, Fatal_No, Fatal_No );


// flag tables
// ----------------------------------------------------------------------------------------------------
#  if ( NLEVEL > 1 )
#  if   ( MODEL == HYDRO )
#  ifndef COSMIC_RAY
   const bool OPT__FLAG_LOHNER_CRAY = false;
#  endif
   const bool Opt__FlagLohner = ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES ||
                                  OPT__FLAG_LOHNER_TEMP || OPT__FLAG_LOHNER_ENTR || OPT__FLAG_LOHNER_CRAY );
#  elif ( MODEL == ELBDM )
   const bool Opt__FlagLohner = OPT__FLAG_LOHNER_DENS;
#  else
#  error : ERROR : unsupported MODEL !!
#  endif
   const bool Compare_Flag_Rho     = ( OPT__FLAG_RHO          ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_RhoGrad = ( OPT__FLAG_RHO_GRADIENT ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_Lohner  = ( Opt__FlagLohner        ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_Angular = ( OPT__FLAG_ANGULAR      ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_Radial  = ( OPT__FLAG_RADIAL       ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_User    = ( OPT__FLAG_USER         ) ? Compare_Yes : Compare_No;

   int H5_Arr_User[1] = { OPT__FLAG_USER_NUM };
   int H5_Arr_Lv  [1] = { (Load_RS) ? MAX_LEVEL : NLEVEL-1    };
   int H5_Arr_Lv_2[2] = { (Load_RS) ? MAX_LEVEL : NLEVEL-1, 2 };
   int H5_Arr_Lv_3[2] = { (Load_RS) ? MAX_LEVEL : NLEVEL-1, 3 };
   int H5_Arr_Lv_4[2] = { (Load_RS) ? MAX_LEVEL : NLEVEL-1, 4 };
   int H5_Arr_Lv_5[2] = { (Load_RS) ? MAX_LEVEL : NLEVEL-1, 5 };
   int H5_Arr_Lv_6[2] = { (Load_RS) ? MAX_LEVEL : NLEVEL-1, 6 };

   HDF5_InputPara->Add( "FlagAngular_CenX",       &FLAG_ANGULAR_CEN_X,            0, NULL,        Compare_Flag_Angular,     Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagAngular_CenY",       &FLAG_ANGULAR_CEN_Y,            0, NULL,        Compare_Flag_Angular,     Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagAngular_CenZ",       &FLAG_ANGULAR_CEN_Z,            0, NULL,        Compare_Flag_Angular,     Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagRadial_CenX",        &FLAG_RADIAL_CEN_X,             0, NULL,        Compare_Flag_Radial,      Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagRadial_CenY",        &FLAG_RADIAL_CEN_Y,             0, NULL,        Compare_Flag_Radial,      Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagRadial_CenZ",        &FLAG_RADIAL_CEN_Z,             0, NULL,        Compare_Flag_Radial,      Fatal_No, Fatal_No );
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   {
   HDF5_InputPara->Add( "FlagTable_Rho",          &FlagTable_Rho         [0],     1, H5_Arr_Lv,   Compare_Flag_Rho,         Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_RhoGradient",  &FlagTable_RhoGradient [0],     1, H5_Arr_Lv,   Compare_Flag_RhoGrad,     Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_Lohner",       &FlagTable_Lohner      [0][0],  2, H5_Arr_Lv_5, Compare_Flag_Lohner,      Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_Angular",      &FlagTable_Angular     [0][0],  2, H5_Arr_Lv_3, Compare_Flag_Angular,     Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_Radial",       &FlagTable_Radial      [0],     1, H5_Arr_Lv,   Compare_Flag_Radial,      Fatal_No, Fatal_No );
   }
   for (int lv=0; lv<MAX_LEVEL; lv++)
   {
   sprintf( Key, "FlagTable_User_Lv%02d", lv );
   HDF5_InputPara->Add( Key,                      &FlagTable_User        [lv][0], 1, H5_Arr_User, Compare_Flag_User,        Fatal_No, Fatal_No );
   }

#  if   ( MODEL == HYDRO )
   const bool Compare_Flag_PresGrad  = ( OPT__FLAG_PRES_GRADIENT ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_Vorticity = ( OPT__FLAG_VORTICITY     ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_Jeans     = ( OPT__FLAG_JEANS         ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   {
   HDF5_InputPara->Add( "FlagTable_PresGradient", &FlagTable_PresGradient[0],     1, H5_Arr_Lv,   Compare_Flag_PresGrad,    Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_Vorticity",    &FlagTable_Vorticity   [0],     1, H5_Arr_Lv,   Compare_Flag_Vorticity,   Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_Jeans",        &FlagTable_Jeans       [0],     1, H5_Arr_Lv,   Compare_Flag_Jeans,       Fatal_No, Fatal_No );
   }
#  ifdef MHD
   const bool Compare_Flag_Current   = ( OPT__FLAG_CURRENT       ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   HDF5_InputPara->Add( "FlagTable_Current",      &FlagTable_Current     [0],     1, H5_Arr_Lv,   Compare_Flag_Current,     Fatal_No, Fatal_No );
#  endif
#  ifdef SRHD
   const bool Compare_Flag_LrtzGrad  = ( OPT__FLAG_LRTZ_GRADIENT ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   HDF5_InputPara->Add( "FlagTable_LrtzGradient", &FlagTable_LrtzGradient[0],     1, H5_Arr_Lv,   Compare_Flag_LrtzGrad,    Fatal_No, Fatal_No );
#  endif
#  ifdef COSMIC_RAY
   const bool Compare_Flag_CRay      = ( OPT__FLAG_CRAY          ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   HDF5_InputPara->Add( "FlagTable_CRay",         &FlagTable_CRay        [0],     1, H5_Arr_Lv,   Compare_Flag_CRay,        Fatal_No, Fatal_No );
#  endif

#  elif ( MODEL == ELBDM )
   const bool Compare_Flag_EngyDens  = ( OPT__FLAG_ENGY_DENSITY  ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_Spectral  = ( OPT__FLAG_SPECTRAL      ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   {
   HDF5_InputPara->Add( "FlagTable_EngyDensity",  &FlagTable_EngyDensity [0][0],  2, H5_Arr_Lv_2, Compare_Flag_EngyDens,    Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_Spectral",     &FlagTable_Spectral    [0][0],  2, H5_Arr_Lv_2, Compare_Flag_Spectral,    Fatal_No, Fatal_No );
   }
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   const bool Compare_Flag_Interfer = ( OPT__FLAG_INTERFERENCE  ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   HDF5_InputPara->Add( "FlagTable_Interference", &FlagTable_Interference[0][0],  2, H5_Arr_Lv_4, Compare_Flag_Interfer,    Fatal_No, Fatal_No );
#  endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
#  endif // #elif ( MODEL == ELBDM )

#  ifdef PARTICLE
   const bool Compare_Flag_NParPatch   = ( OPT__FLAG_NPAR_PATCH    ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_NParCell    = ( OPT__FLAG_NPAR_CELL     ) ? Compare_Yes : Compare_No;
   const bool Compare_Flag_ParMassCell = ( OPT__FLAG_PAR_MASS_CELL ) ? Compare_Yes : Compare_No;
   if ( (Load_RS  &&  MAX_LEVEL > 0)  ||  !Load_RS )
   {
   HDF5_InputPara->Add( "FlagTable_NParPatch",    &FlagTable_NParPatch   [0],     1, H5_Arr_Lv,   Compare_Flag_NParPatch,   Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_NParCell",     &FlagTable_NParCell    [0],     1, H5_Arr_Lv,   Compare_Flag_NParCell,    Fatal_No, Fatal_No );
   HDF5_InputPara->Add( "FlagTable_ParMassCell",  &FlagTable_ParMassCell [0],     1, H5_Arr_Lv,   Compare_Flag_ParMassCell, Fatal_No, Fatal_No );
   }
#  endif
#  endif // #if ( NLEVEL > 1 )


// user derived fields
// ----------------------------------------------------------------------------------------------------
   HDF5_InputPara->Add( "UserDerField_Num", &UserDerField_Num,         0, NULL, Compare_No, NULL_BOOL, NULL_BOOL );
   for (int v=0; v<UserDerField_Num; v++)
   {
   sprintf( Key, "UserDerField_Label%02d", v );
   HDF5_InputPara->Add( Key,                &UserDerField_Label[v][0], 0, NULL, Compare_No, NULL_BOOL, NULL_BOOL );
   }
   for (int v=0; v<UserDerField_Num; v++)
   {
   sprintf( Key, "UserDerField_Unit%02d", v );
   HDF5_InputPara->Add( Key,                &UserDerField_Unit [v][0], 0, NULL, Compare_No, NULL_BOOL, NULL_BOOL );
   }

} // FUNCTION : FillIn_HDF5_InputPara



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_HDF5_UserPara_Template
// Description :  Template for storing user-specified parameters in an HDF5 snapshot at User/UserPara
//
// Note        :  1. This function is only called by the root MPI rank
//                2. Support int, uint, long, ulong, bool, float, double, string, and their corresponding array types
//                3. HDF5_UserPara MUST store at least one parameter
//                4. The second argument passed to HDF5_UserPara->Add() can be:
//                   a. pointer  : Support int, uint, long, ulong, bool, float, double, and string datatypes
//                   b. variable : Support int, uint, long, ulong, bool, float, and double datatypes
//                5. Linked to the function pointer Output_HDF5_UserPara_Ptr
//                6. The string size MUST be `MAX_STRING`
//                7. The memeory space of the array MUST be continuous
//
// Parameter   :  HDF5_UserPara : Structure storing all parameters to be written
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_HDF5_UserPara_Template( HDF5_Output_t *HDF5_UserPara )
{

   const bool Compare_Yes = true, Compare_No = false;
   const bool   Fatal_Yes = true,   Fatal_No = false;

   int     ArrLen_1D[1] = {3};
   int     ArrLen_2D[2] = {2, 3};

   long    long_data                   = 1L;
   int     IntArr_1D[4]                = {1, 2, 3, 4};
   double  DblArr_2D[2][3]             = { {1., 2., 3.}, {4., 5., 6.}};
   char    StrArr_2D[2][3][MAX_STRING] = { {"1", "2", "3"}, {"4", "5", "6"}};

// example
// HDF5_UserPara->Add( "Your_Data1_Label", &Your_Data1_Pointer, arr_dim, arr_len, compare_or_not, fatal_exist, fatal_comapre );
// HDF5_UserPara->Add( "Your_Data2_Label", const_var, 0, NULL, Compare_No, NULL_BOOL, NULL_BOOL );

   HDF5_UserPara->Add( "const_long_data", &long_data,          0, NULL,      Compare_No, NULL_BOOL, NULL_BOOL );
   HDF5_UserPara->Add( "int_1D_array",    &IntArr_1D[0],       1, ArrLen_1D, Compare_No, NULL_BOOL, NULL_BOOL );
   HDF5_UserPara->Add( "double_2D_array", &DblArr_2D[0][0],    2, ArrLen_2D, Compare_No, NULL_BOOL, NULL_BOOL );
   HDF5_UserPara->Add( "string_2D_array", &StrArr_2D[0][0][0], 2, ArrLen_2D, Compare_No, NULL_BOOL, NULL_BOOL );

} // FUNCTION : Output_HDF5_UserPara_Template



#endif // #ifdef SUPPORT_HDF5
