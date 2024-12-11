#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_ByFile_Default( real fluid_out[], const real fluid_in[], const int nvar_in,
                                 const double x, const double y, const double z, const double Time,
                                 const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*Init_ByFile_User_Ptr)( real fluid_out[], const real fluid_in[], const int nvar_in,
                              const double x, const double y, const double z, const double Time,
                              const int lv, double AuxArray[] ) = Init_ByFile_Default;

static void Init_ByFile_AssignData( const char UM_Filename[], const int UM_lv, const int UM_lv0, const int UM_NVar,
                                    const int UM_LoadNRank, const UM_IC_Format_t UM_Format, const long UM_Size3D[][3],
                                    const int FlagPatch[][6] );
static void Load_RefineRegion( const char Filename[] );
static void Flag_RefineRegion( const int lv, const int FlagPatch[6] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile
// Description :  Set up initial condition from an input uniform-mesh array
//
// Note        :  1. Create levels from 0 to OPT__UM_IC_LEVEL
//                   --> No levels above OPT__UM_IC_LEVEL will be created unless OPT__UM_IC_REFINE is enabled,
//                       for which levels OPT__UM_IC_LEVEL+1 ~ MAX_LEVEL may also be generated based on the
//                       given refinement criteria
//                   --> ALL patches on levels 0 to OPT__UM_IC_LEVEL will be created. In other words,
//                       the simulation domain will be fully refined to level OPT__UM_IC_LEVEL.
//                       --> But if OPT__UM_IC_DOWNGRADE is enabled, patches on levels 1~OPT__UM_IC_LEVEL
//                           may be removed if not satisfying the refinement criteria
//                2. The uniform-mesh input file must be named as "UM_IC"
//                3. This function can load any number of input variables per cell (from 1 to NCOMP_TOTAL)
//                   --> Determined by the input parameter "OPT__UM_IC_NVAR"
//                   --> If "OPT__UM_IC_NVAR < NCOMP_TOTAL", one must specify the way to assign values to all
//                       variables using the function pointer Init_ByFile_User_Ptr()
//                       --> Two exceptions:
//                           (1) When enabling DUAL_ENERGY, we will calculate the dual-energy field
//                               (e.g., entropy) directly instead of loading it from the disk
//                               --> Only NCOMP_TOTAL-1 (which equals 5+NCOMP_PASSIVE_USER currently) fields
//                                   should be stored in "UM_IC"
//                           (2) For ELBDM:
//                               --> Only NCOMP_TOTAL-1 (which equals 2+NCOMP_PASSIVE_USER currently) fields
//                                   should be stored in "UM_IC"
//
//                               ELBDM_SCHEME == ELBDM_WAVE:
//                               --> We will load the real and imaginary parts from the disk on all levels
//                                   and calculate the density field from the input wave function
//                                   directly instead of loading it from the disk.
//
//                               ELBDM_SCHEME == ELBDM_HYBRID:
//                               --> We will load the density and phase fields from the disk on all levels.
//                                   There is no need to separately calculate the density field.
//                4. The data format of the UM_IC file is controlled by the runtime parameter OPT__UM_IC_FORMAT
//                5. Does not work with rectangular domain decomposition anymore
//                   --> Must enable either SERIAL or LOAD_BALANCE
//                6. OpenMP is not supported yet
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ByFile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const char   UM_Filename[] = "UM_IC";
   const char   RR_Filename[] = "Input__UM_IC_RefineRegion";   // RR = Refinement Region
#  if ( defined PARTICLE  &&  defined LOAD_BALANCE )
   const double Par_Weight    = amr->LB->Par_Weight;
#  else
   const double Par_Weight    = 0.0;
#  endif
#  ifdef LOAD_BALANCE
   const UseLBFunc_t UseLB    = USELB_YES;
#  else
   const UseLBFunc_t UseLB    = USELB_NO;
#  endif


// 1. determine the refinement regions for OPT__UM_IC_NLEVEL > 1
//    FlagPatch: the range of patches to be flagged on each side along each direction on levels
//               OPT__UM_IC_LEVEL ~ OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-2
//    --> the first patch along each direction is labeled as 0
//    --> [lv][0/1]: leftmost/rightmost patches along the x direction to be flagged on level lv
   int FlagPatch[NLEVEL][6];

   if ( OPT__UM_IC_NLEVEL > 1 )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting refinement level ...\n" );

//    1-1. load the information of refinement regions for OPT__UM_IC_NLEVEL > 1
      Load_RefineRegion( RR_Filename );


//    1-2. set FlagPatch[]
      const int (*RefineRegion)[6] = ( int(*)[6] )UM_IC_RefineRegion;

//    count the number of patches to be skipped
      for (int t=0; t<OPT__UM_IC_NLEVEL-1; t++)
      for (int s=0; s<6; s++)
      {
         FlagPatch[t][s] = RefineRegion[t][s];

//       include the number of patches to be skipped on the parent level
         if ( t > 0 )   FlagPatch[t][s] += FlagPatch[t-1][s]*2;
      }

//    set the rightmost patches to be flagged
      for (int t=0; t<OPT__UM_IC_NLEVEL-1; t++)
      for (int d=0; d<3; d++)
      {
         const int lv = OPT__UM_IC_LEVEL + t;
         const int s  = 2*d + 1;
         const int NP = Mis_Scale2Cell( amr->BoxScale[d], lv ) / PS1;

         FlagPatch[t][s] = NP - 1 - FlagPatch[t][s];

//       check if rightmost patch > leftmost patch
         if ( FlagPatch[t][s] <= FlagPatch[t][s-1] )
            Aux_Error( ERROR_INFO, "FlagPatch: [%d][%d] (%d) <= [%d][%d] (%d) !!\n",
                       t, s, FlagPatch[t][s], t, s-1, FlagPatch[t][s-1] );
      }

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
      for (int t=0; t<OPT__UM_IC_NLEVEL-1; t++)
         Aux_Message( stdout, "   FlagPatch[%d]: [%4d, %4d, %4d, %4d, %4d, %4d]\n",
                   t, FlagPatch[t][0], FlagPatch[t][1], FlagPatch[t][2], FlagPatch[t][3], FlagPatch[t][4], FlagPatch[t][5] );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting refinement level ... done\n" );

   } // if ( OPT__UM_IC_NLEVEL > 1 )


// 1-3. determine the file size on each level
   long UM_Size3D[NLEVEL][3];

// level OPT__UM_IC_LEVEL
   for (int d=0; d<3; d++)
      UM_Size3D[0][d] = NX0_TOT[d]*(1<<OPT__UM_IC_LEVEL);

// levels OPT__UM_IC_LEVEL+1 ~ OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1
   for (int t=1; t<OPT__UM_IC_NLEVEL; t++)
   for (int d=0; d<3; d++)
      UM_Size3D[t][d] = ( FlagPatch[t-1][2*d+1] - FlagPatch[t-1][2*d] + 1 )*PS2;



// check
#  if ( !defined SERIAL  &&  !defined LOAD_BALANCE )
      Aux_Error( ERROR_INFO, "must enable either SERIAL or LOAD_BALANCE for %s !!\n", __FUNCTION__ );
#  endif

   if ( OPT__UM_IC_LEVEL < 0  ||  OPT__UM_IC_LEVEL > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "OPT__UM_IC_LEVEL (%d) > TOP_LEVEL (%d) !!\n", OPT__UM_IC_LEVEL, TOP_LEVEL );

   if ( OPT__UM_IC_NVAR < 1  ||  OPT__UM_IC_NVAR > NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "invalid OPT__UM_IC_NVAR = %d (accepeted range: %d ~ %d) !!\n",
                 OPT__UM_IC_NVAR, 1, NCOMP_TOTAL );

   if ( !Aux_CheckFileExist(UM_Filename) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", UM_Filename );

   if ( OPT__UM_IC_LEVEL + OPT__UM_IC_NLEVEL - 1 > MAX_LEVEL )
      Aux_Error( ERROR_INFO, "OPT__UM_IC_LEVEL (%d) + OPT__UM_IC_NLEVEL (%d) - 1 = %d > MAX_LEVEL (%d) !!\n",
                 OPT__UM_IC_LEVEL, OPT__UM_IC_NLEVEL, OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1, MAX_LEVEL );

   if ( OPT__RESET_FLUID_INIT  &&  Flu_ResetByUser_API_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_API_Ptr == NULL for OPT__RESET_FLUID_INIT !!\n" );

// check file size
   long FileSize, ExpectSize;

   FILE *FileTemp = fopen( UM_Filename, "rb" );
   fseek( FileTemp, 0, SEEK_END );
   FileSize = ftell( FileTemp );
   fclose( FileTemp );

   ExpectSize = 0;
   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
      ExpectSize += long(OPT__UM_IC_NVAR)*UM_Size3D[t][0]*UM_Size3D[t][1]*UM_Size3D[t][2]*( (OPT__UM_IC_FLOAT8)?sizeof(double):sizeof(float) );

   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> (%ld) != expected (%ld) !!\n", UM_Filename, FileSize, ExpectSize );


   MPI_Barrier( MPI_COMM_WORLD );



// 2. allocate all real patches on levels 0 ~ OPT__UM_IC_LEVEL
   const bool FindHomePatchForPar_Yes = true;
   const bool FindHomePatchForPar_No  = false;

   for (int lv=0; lv<=OPT__UM_IC_LEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Allocating level %d ...\n", lv );

//    associate particles with home patches on the level **OPT__UM_IC_LEVEL** only
      Init_UniformGrid( lv, (lv==OPT__UM_IC_LEVEL)?FindHomePatchForPar_Yes:FindHomePatchForPar_No );

//    construct IdxList_Real[]
//    --> necessary for calling LB_Init_LoadBalance() later with Redistribute_No
#     ifdef LOAD_BALANCE
      if ( amr->LB->IdxList_Real         [lv] != NULL )   delete [] amr->LB->IdxList_Real         [lv];
      if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )   delete [] amr->LB->IdxList_Real_IdxTable[lv];

      amr->LB->IdxList_Real         [lv] = new long [ amr->NPatchComma[lv][1] ];
      amr->LB->IdxList_Real_IdxTable[lv] = new int  [ amr->NPatchComma[lv][1] ];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

      Mis_Heapsort( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );
#     endif

//    get the total number of real patches
      Mis_GetTotalPatchNumber( lv );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Allocating level %d ... done\n", lv );
   } // for (int lv=0; lv<=OPT__UM_IC_LEVEL; lv++)



// 3. initialize load-balancing (or construct patch relation for SERIAL)
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing patch relation ...\n" );

#  ifdef LOAD_BALANCE
// no need to redistribute patches since Init_UniformGrid() already takes into account load balancing
// --> but particle weighting is not considered yet
// --> we will invoke LB_Init_LoadBalance() again after Flu_FixUp_Restrict() for that

// must not reset load-balance variables (i.e., must adopt ResetLB_No) to avoid overwritting IdxList_Real[]
// and IdxList_Real_IdxList[] already set above
   const double ParWeight_Zero   = 0.0;
   const bool   Redistribute_Yes = true;
   const bool   Redistribute_No  = false;
   const bool   SendGridData_Yes = true;
   const bool   SendGridData_No  = false;
   const bool   ResetLB_Yes      = true;
   const bool   ResetLB_No       = false;
   const bool   SortRealPatch_No = false;
   const int    AllLv            = -1;

   LB_Init_LoadBalance( Redistribute_No, SendGridData_No, ParWeight_Zero, ResetLB_No, SortRealPatch_No, AllLv );

#  else // for SERIAL

   for (int lv=0; lv<=OPT__UM_IC_LEVEL; lv++)
   {
//    set up BaseP[]
      if ( lv == 0 )
      Init_RecordBasePatch();

//    construct father-son relation
      if ( lv > 0 )
      FindFather( lv, 1 );

//    construct sibling relation
      SiblingSearch( lv );

//    allocate flux and electric field arrays on "lv-1"
//    --> must do this after constructing the patch relation on lv-1 and lv
      if ( lv > 0 )
      {
         if ( amr->WithFlux )       Flu_AllocateFluxArray( lv-1 );
#        ifdef MHD
         if ( amr->WithElectric )   MHD_AllocateElectricArray( lv-1 );
#        endif
      }
   }
#  endif // #ifdef LOAD_BALANCE ... else ...

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing patch relation ... done\n" );



// 4. assign data on level OPT__UM_IC_LEVEL by the input file UM_IC
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Assigning data on level %d ...\n", OPT__UM_IC_LEVEL );

   Init_ByFile_AssignData( UM_Filename, OPT__UM_IC_LEVEL, OPT__UM_IC_LEVEL, OPT__UM_IC_NVAR, OPT__UM_IC_LOAD_NRANK,
                           OPT__UM_IC_FORMAT, UM_Size3D, FlagPatch );

   if ( OPT__RESET_FLUID_INIT )
      Flu_ResetByUser_API_Ptr( OPT__UM_IC_LEVEL, amr->FluSg[OPT__UM_IC_LEVEL], amr->MagSg[OPT__UM_IC_LEVEL],
                               Time[OPT__UM_IC_LEVEL], 0.0 );

#  ifdef LOAD_BALANCE
   Buf_GetBufferData( OPT__UM_IC_LEVEL, amr->FluSg[OPT__UM_IC_LEVEL], amr->MagSg[OPT__UM_IC_LEVEL], NULL_INT,
                      DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Assigning data on level %d ... done\n", OPT__UM_IC_LEVEL );



// 5. construct levels OPT__UM_IC_LEVEL+1 to OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1 by the input file UM_IC
   for (int t=0; t<OPT__UM_IC_NLEVEL-1; t++)
   {
      const int  FaLv         = OPT__UM_IC_LEVEL + t;
      const int  SonLv        = FaLv + 1;
      const bool AllocData_No = false;

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level %d ...\n", SonLv );

//    flag FaLv
      Flag_RefineRegion( FaLv, FlagPatch[t] );

//    create SonLv
//    --> do not allocate grid data here to reduce memory consumption
#     ifdef LOAD_BALANCE
      LB_Init_Refine( FaLv, AllocData_No );
#     else // for SERIAL
      Init_Refine( FaLv );
#     endif

//    get and check the total number of patches on SonLv
      Mis_GetTotalPatchNumber( SonLv );
      const int NPatchExpect = 8*( FlagPatch[t][1] - FlagPatch[t][0] + 1 )
                                *( FlagPatch[t][3] - FlagPatch[t][2] + 1 )
                                *( FlagPatch[t][5] - FlagPatch[t][4] + 1 );
      if ( NPatchTotal[SonLv] != NPatchExpect )
         Aux_Error( ERROR_INFO, "Total number of patches on level %d (%d) != expected (%d) !!\n",
                    SonLv, NPatchTotal[SonLv], NPatchExpect );

//    redistribute patches for load balancing
//    --> no need to send grid data since it hasn't been assigned yet
#     ifdef LOAD_BALANCE
      LB_Init_LoadBalance( Redistribute_Yes, SendGridData_No, Par_Weight, ResetLB_Yes, SortRealPatch_No, SonLv );
#     endif

//    assign data on SonLv
      Init_ByFile_AssignData( UM_Filename, SonLv, OPT__UM_IC_LEVEL, OPT__UM_IC_NVAR, OPT__UM_IC_LOAD_NRANK,
                              OPT__UM_IC_FORMAT, UM_Size3D, FlagPatch );

      if ( OPT__RESET_FLUID_INIT )
         Flu_ResetByUser_API_Ptr( SonLv, amr->FluSg[SonLv], amr->MagSg[SonLv], Time[SonLv], 0.0 );

//    fill the buffer patches on SonLv
#     ifdef LOAD_BALANCE
      Buf_GetBufferData( SonLv, amr->FluSg[SonLv], amr->MagSg[SonLv], NULL_INT,
                         DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level %d ... done\n", SonLv );
   } // for (int t=0; t<OPT__UM_IC_NLEVEL-1; t++)



// 6. restrict data on levels 0 ~ MAX_LEVEL-1
//    --> assign data on levels 0 ~ OPT__UM_IC_LEVEL-1
//    --> ensure data consistency on levels OPT__UM_IC_LEVEL ~ MAX_LEVEL when OPT__UM_IC_NLEVEL>1
   for (int lv=MAX_LEVEL-1; lv>=0; lv--)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Restricting level %d ... ", lv );

      const long ResVar = ( lv < OPT__UM_IC_LEVEL ) ? _TOTAL : FixUpVar_Restrict;

      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT,
                          ResVar, _MAG );

#     ifdef LOAD_BALANCE
      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_RESTRICT, ResVar, _MAG, NULL_INT,    USELB_YES );

//    use DATA_GENERAL instead of DATA_AFTER_FIXUP since we haven't call Buf_GetBufferData() on levels 0 ~ OPT__UM_IC_LEVEL-1
      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL,  ResVar, _MAG, Flu_ParaBuf, USELB_YES );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }



// 7. optimize load-balancing to take into account particle weighting
#  if ( defined PARTICLE  &&  defined LOAD_BALANCE )
   if ( Par_Weight > 0.0 )
      LB_Init_LoadBalance( Redistribute_Yes, SendGridData_Yes, Par_Weight, ResetLB_Yes, SortRealPatch_No, AllLv );
#  endif



// 8. derefine the uniform-mesh data from levels OPT__UM_IC_LEVEL to 1
   if ( OPT__UM_IC_DOWNGRADE )
   for (int lv=OPT__UM_IC_LEVEL-1; lv>=0; lv--)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Downgrading level %d ...\n", lv+1 );

      Flag_Real( lv, UseLB );

      Refine( lv, UseLB );

#     ifdef LOAD_BALANCE
      LB_Init_LoadBalance( Redistribute_Yes, SendGridData_Yes, Par_Weight, ResetLB_Yes, SortRealPatch_No, lv+1 );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Downgrading level %d ... done\n", lv+1 );
   } // for (int lv=OPT__UM_IC_LEVEL-1; lv>=0; lv--)



// 9. refine the uniform-mesh data from levels OPT__UM_IC_LEVEL to MAX_LEVEL-1
//    --> we currently do not apply OPT__RESET_FLUID_INIT after Refine() because
//        (1) the existing patches from level OPT__UM_IC_LEVEL to level OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1
//            already include the reset results (by Flu_ResetByUser_API_Ptr() directly),
//        (2) the refined patches on levels higher than OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1
//            also include the reset results (by interpolation), and
//        (3) to avoid double-counting the reset data when Flu_ResetByUser_API_Ptr() uses fluid[]+=... instead of fluid[]=...
   if ( OPT__UM_IC_REFINE )
   for (int lv=OPT__UM_IC_LEVEL; lv<MAX_LEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ...\n", lv );

      Flag_Real( lv, UseLB );

      Refine( lv, UseLB );

#     ifdef LOAD_BALANCE
      LB_Init_LoadBalance( Redistribute_Yes, SendGridData_Yes, Par_Weight, ResetLB_Yes, SortRealPatch_No, lv+1 );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ... done\n", lv );
   } // for (int lv=OPT__UM_IC_LEVEL; lv<MAX_LEVEL; lv++)



// 10. restrict data again
//    --> for bitwise reproducibility only
//    --> strictly speaking, it is only necessary for C-binary output (i.e., OPT__OUTPUT_TOTAL=2)
//        since that output format does not store non-leaf patch data
#  ifdef BITWISE_REPRODUCIBILITY
   for (int lv=MAX_LEVEL-1; lv>=0; lv--)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Restricting level %d ... ", lv );

      if ( NPatchTotal[lv+1] == 0 )    continue;

//    no need to restrict potential since it will be recalculated later
      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT,
                          FixUpVar_Restrict, _MAG );

#     ifdef LOAD_BALANCE
      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_RESTRICT,    FixUpVar_Restrict, _MAG, NULL_INT,    USELB_YES );

      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_AFTER_FIXUP, FixUpVar_Restrict, _MAG, Flu_ParaBuf, USELB_YES );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }
#  endif // # ifdef BITWISE_REPRODUCIBILITY


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ByFile



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile_AssignData
// Description :  Use the input uniform-mesh array stored in the file "UM_Filename" to assign data to all
//                real patches on level "UM_lv"
//
// Note        :  1. The function pointer Init_ByFile_User_Ptr() points to Init_ByFile_Default() by default
//                   but may be overwritten by various test problem initializers
//                2. Can be applied to levels OPT__UM_IC_LEVEL ~ OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1
//
// Parameter   :  UM_Filename  : Target file name
//                UM_lv        : Target AMR level --> OPT__UM_IC_LEVEL ~ OPT__UM_IC_LEVEL+OPT__UM_IC_NLEVEL-1
//                UM_lv0       : First AMR level in the file UM_Filename (i.e., OPT__UM_IC_LEVEL)
//                UM_NVar      : Number of variables
//                UM_LoadNRank : Number of parallel I/O
//                UM_Format    : Data format of the file UM_Filename
//                               --> UM_IC_FORMAT_VZYX: [field][z][y][x] in a row-major order
//                                   UM_IC_FORMAT_ZYXV: [z][y][x][field] in a row-major order
//                UM_Size3D    : Size of the file UM_Filename on each level
//                FlagPatch    : Range of patches to be flagged on each level
//                               --> For details, see the description under its declaration at the beginning
//                                   of Init_ByFile()
//
// Return      :  amr->patch->fluid
//-------------------------------------------------------------------------------------------------------
void Init_ByFile_AssignData( const char UM_Filename[], const int UM_lv, const int UM_lv0, const int UM_NVar,
                             const int UM_LoadNRank, const UM_IC_Format_t UM_Format, const long UM_Size3D[][3],
                             const int FlagPatch[][6] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "      Loading data from the input file on level %d ...\n", UM_lv );


// check
   if ( Init_ByFile_User_Ptr == NULL )  Aux_Error( ERROR_INFO, "Init_ByFile_User_Ptr == NULL !!\n" );


   const int    dlv         = UM_lv - UM_lv0;
   const long   UM_Size1v   = UM_Size3D[dlv][0]*UM_Size3D[dlv][1]*UM_Size3D[dlv][2];
   const int    NVarPerLoad = ( UM_Format == UM_IC_FORMAT_ZYXV ) ? UM_NVar : 1;
   const int    scale       = amr->scale[UM_lv];
   const double dh          = amr->dh[UM_lv];

   long   Offset3D_File0[3], Offset_File0, Offset_File, Offset_PG, Offset_lv;
   real   fluid_in[UM_NVar], fluid_out[NCOMP_TOTAL];
   double x, y, z;


// determine the load_data_size and allocate buffer for loading UM_IC
   size_t load_data_size = ( OPT__UM_IC_FLOAT8 ) ? sizeof(double) : sizeof(float);
   char*  PG_Data        = new char [ CUBE(PS2)*UM_NVar*load_data_size ];


// calculate the file offset of the target level
   Offset_lv = 0;
   for (int t=0; t<dlv; t++)
      Offset_lv += long(UM_NVar)*UM_Size3D[t][0]*UM_Size3D[t][1]*UM_Size3D[t][2]*load_data_size;


// load data with UM_LoadNRank ranks at a time
   for (int TRank0=0; TRank0<MPI_NRank; TRank0+=UM_LoadNRank)
   {
      if ( MPI_Rank >= TRank0  &&  MPI_Rank < TRank0+UM_LoadNRank )
      {
         if ( MPI_Rank == TRank0 )  Aux_Message( stdout, "      Loading ranks %4d -- %4d ... ",
                                                 TRank0, MIN(TRank0+UM_LoadNRank-1, MPI_NRank-1) );

         FILE *File = fopen( UM_Filename, "rb" );

//       load one patch group at a time
         for (int PID0=0; PID0<amr->NPatchComma[UM_lv][1]; PID0+=8)
         {
//          calculate Offset_File0, which is the file offset of the target patch group relative to Offset_lv
            for (int d=0; d<3; d++)    Offset3D_File0[d] = amr->patch[0][UM_lv][PID0]->corner[d] / scale;

            if ( dlv > 0 )
            for (int d=0; d<3; d++)
            {
               Offset3D_File0[d] -= FlagPatch[dlv-1][2*d]*PS2;

               if ( Offset3D_File0[d] < 0 )
                  Aux_Error( ERROR_INFO, "Offset3D_File0[%d] = %ld < 0 !!\n", d, Offset3D_File0[d] );
            }

            Offset_File0  = IDX321( Offset3D_File0[0], Offset3D_File0[1], Offset3D_File0[2],
                                    UM_Size3D[dlv][0], UM_Size3D[dlv][1] );
            Offset_File0 *= (long)NVarPerLoad*load_data_size;


//          load data from the disk (one row at a time)
            Offset_PG = 0;

            for (int v=0; v<UM_NVar; v+=NVarPerLoad )
            {
               for (int k=0; k<PS2; k++)
               for (int j=0; j<PS2; j++)
               {
                  Offset_File = Offset_lv + Offset_File0
                                + (long)NVarPerLoad*load_data_size*( ((long)k*UM_Size3D[dlv][1] + j)*UM_Size3D[dlv][0] )
                                + v*UM_Size1v*load_data_size;

                  fseek( File, Offset_File, SEEK_SET );
                  fread( PG_Data+Offset_PG, load_data_size, NVarPerLoad*PS2, File );

//                verify that the file size is not exceeded
                  if ( feof(File) )   Aux_Error( ERROR_INFO, "reaching the end of the file \"%s\" !!\n", UM_Filename );

                  Offset_PG += NVarPerLoad*PS2*load_data_size;
               }
            }


//          copy data to each patch
            for (int LocalID=0; LocalID<8; LocalID++)
            {
               const int PID    = PID0 + LocalID;
               const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1 );
               const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1 );
               const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1 );

               for (int k=0; k<PS1; k++)  {  z = amr->patch[0][UM_lv][PID]->EdgeL[2] + (k+0.5)*dh;
               for (int j=0; j<PS1; j++)  {  y = amr->patch[0][UM_lv][PID]->EdgeL[1] + (j+0.5)*dh;
               for (int i=0; i<PS1; i++)  {  x = amr->patch[0][UM_lv][PID]->EdgeL[0] + (i+0.5)*dh;

                  Offset_PG = (long)NVarPerLoad*IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS2, PS2 )*load_data_size;

                  if ( UM_Format == UM_IC_FORMAT_ZYXV )
                  {
                     if ( OPT__UM_IC_FLOAT8 )
                        for (int v=0; v<UM_NVar; v++) fluid_in[v] = (real)( *((double*)( PG_Data + Offset_PG + v*load_data_size )) ) ;
                     else
                        for (int v=0; v<UM_NVar; v++) fluid_in[v] = (real)( *((float* )( PG_Data + Offset_PG + v*load_data_size )) ) ;
                  }

                  else
                  {
                     if ( OPT__UM_IC_FLOAT8 )
                        for (int v=0; v<UM_NVar; v++) fluid_in[v] = (real)( *((double*)( PG_Data + Offset_PG + v*CUBE(PS2)*load_data_size )) );
                     else
                        for (int v=0; v<UM_NVar; v++) fluid_in[v] = (real)( *((float* )( PG_Data + Offset_PG + v*CUBE(PS2)*load_data_size )) );
                  }

                  Init_ByFile_User_Ptr( fluid_out, fluid_in, UM_NVar, x, y, z, Time[UM_lv], UM_lv, NULL );

                  for (int v=0; v<NCOMP_TOTAL; v++)
                     amr->patch[ amr->FluSg[UM_lv] ][UM_lv][PID]->fluid[v][k][j][i] = fluid_out[v];
               }}}
            } // for (int LocalID=0; LocalID<8; LocalID++)
         } // for (int PID0=0; PID0<amr->NPatchComma[UM_lv][1]; PID0+=8)

         fclose( File );

         if ( MPI_Rank == TRank0 )  Aux_Message( stdout, "done\n" );
      } // if ( MPI_Rank >= TRank0  &&  MPI_Rank < TRank0+UM_LoadNRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TRank0=0; TRank0<MPI_NRank; TRank0+=UM_LoadNRank)

   delete [] PG_Data;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "      Loading data from the input file on level %d ... done\n", UM_lv );

} // FUNCTION : Init_ByFile_AssignData



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile_Default
// Description :  Function to actually set the fluid field from the input uniform-mesh array
//
// Note        :  1. Invoked by Init_ByFile_AssignData() using the function pointer Init_ByFile_User_Ptr()
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Does not floor and normalize passive scalars
//                3. Calculate the dual-energy variable automatically instead of load it from the disk
//                   --> When adopting DUAL_ENERGY, the input uniform-mesh array must NOT include the dual-energy
//                       variable
//                4. ELBDM:
//                   ELBDM_SCHEME == ELBDM_WAVE:
//                   --> Calculate the density field automatically instead of loading it from the disk for ELBDM
//                   --> The input uniform-mesh array must NOT include the density field
//                   ELBDM_SCHEME == ELBDM_HYBRID:
//                   --> We will load the density and phase fields from the disk on all levels
//                   --> There is no need to separately calculate the density field
//                5. Assuming nvar_in (i.e., OPT__UM_IC_NVAR) == NCOMP_TOTAL
//                   --> Unless either DUAL_ENERGY or ELBDM is adopted, for which it assumes nvar_in == NCOMP_TOTAL-1
//
// Parameter   :  fluid_out : Fluid field to be set
//                fluid_in  : Fluid field loaded from the uniform-mesh array (UM_IC)
//                nvar_in   : Number of variables in fluid_in
//                x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target AMR level
//                AuxArray  : Auxiliary array
//
// Return      :  fluid_out
//-------------------------------------------------------------------------------------------------------
void Init_ByFile_Default( real fluid_out[], const real fluid_in[], const int nvar_in,
                          const double x, const double y, const double z, const double Time,
                          const int lv, double AuxArray[] )
{

#  ifdef GAMER_DEBUG
#  if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
   if ( nvar_in != NCOMP_TOTAL-1 )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != NCOMP_TOTAL-1 (%d) when enabling DUAL_ENERGY !!\n", nvar_in, NCOMP_TOTAL-1 );

#  elif ( MODEL == ELBDM )
   if ( nvar_in != NCOMP_TOTAL-1 )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != NCOMP_TOTAL-1 (%d) for ELBDM !!\n", nvar_in, NCOMP_TOTAL-1 );

#  else
   if ( nvar_in != NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != NCOMP_TOTAL (%d) !!\n", nvar_in, NCOMP_TOTAL );
#  endif
#  endif // #ifdef GAMER_DEBUG

   for (int v_in=0, v_out=0; v_in<nvar_in; v_in++, v_out++)
   {
//    skip the dual-energy field for HYDRO
#     if   ( MODEL == HYDRO )
#     ifdef DUAL_ENERGY
      if ( v_out == DUAL )    v_out ++;
#     endif

//    skip the density field for ELBDM_SCHEME == ELBDM_WAVE
#     elif ( MODEL == ELBDM )
#     if ( ELBDM_SCHEME == ELBDM_WAVE )
      if ( v_out == DENS )    v_out ++;
#     endif
#     endif // MODEL

      fluid_out[v_out] = fluid_in[v_in];
   }

// calculate the dual-energy field for HYDRO
#  if   ( MODEL == HYDRO )

#  ifdef MHD
   Aux_Error( ERROR_INFO, "MHD is NOT supported yet !!\n" );
   const real Emag = NULL_REAL;
#  else
   const real Emag = NULL_REAL;
#  endif

#  ifdef DUAL_ENERGY
   fluid_out[DUAL] = Hydro_Con2Dual( fluid_in[DENS], fluid_in[MOMX], fluid_in[MOMY], fluid_in[MOMZ], fluid_in[ENGY], Emag,
                                     EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#  endif

#  elif ( MODEL == ELBDM )

#  if   ( ELBDM_SCHEME == ELBDM_WAVE )
// calculate the density field for ELBDM wave scheme
   fluid_out[DENS] = SQR( fluid_out[REAL] ) + SQR( fluid_out[IMAG] );

#  elif ( ELBDM_SCHEME == ELBDM_HYBRID )
// convert density and phase to real and imaginary part on wave levels for hybrid scheme
   if ( amr->use_wave_flag[lv] ) {
      const real Phase = fluid_out[PHAS];
      const real Amp   = SQRT( fluid_out[DENS] );
      fluid_out[REAL] = Amp * COS( Phase );
      fluid_out[IMAG] = Amp * SIN( Phase );
   } else {
      fluid_out[STUB] = (real)0.0;
   }
#  endif // ELBDM_SCHEME
#  endif // MODEL

} // FUNCTION : Init_ByFile_Default



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_RefineRegion
// Description :  Load the information of refinement regions from the table "Filename"
//
// Note        :  1. Invoked by Init_ByFile()
//                2. Only useful when OPT__UM_IC_NLEVEL>1
//                3. UM_IC_RefineRegion[] will be allocated here
//                   --> Deallocated by End_MemFree()
//
// Parameter   :  Filename : Target file name
//
// Return      :  UM_IC_RefineRegion[]
//-------------------------------------------------------------------------------------------------------
void Load_RefineRegion( const char Filename[] )
{

// check
   if ( ! Aux_CheckFileExist(Filename) )  Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Filename );


   const bool RowMajor_Yes = true;                 // load data into the row-major order
   const bool AllocMem_Yes = true;                 // allocate memory for UM_IC_RefineRegion[]
   const int  NCol        = 6;                     // total number of columns to load
   const int  Col[NCol]   = {1,2,3,4,5,6};         // target columns (skip the first column)
   const int  NLv         = OPT__UM_IC_NLEVEL - 1; // number of AMR levels to load

   int NRow, (*RefineRegion)[6];

   NRow         = Aux_LoadTable( UM_IC_RefineRegion, Filename, NCol, Col, RowMajor_Yes, AllocMem_Yes );
   RefineRegion = ( int(*)[6] )UM_IC_RefineRegion;


// check
   if ( NRow < NLv )
      Aux_Error( ERROR_INFO, "Number of rows in %s = %d < OPT__UM_IC_NLEVEL-1 = %d !!\n",
                 Filename, NRow, NLv );

   for (int t=0; t<NLv; t++)
   for (int s=0; s<6; s++)
      if ( RefineRegion[t][s] < 1 )
         Aux_Error( ERROR_INFO, "RefineRegion[%d][%d] = %d < 1 !!\n", t, s, RefineRegion[t][s] );


// log
   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
   {
      Aux_Message( stdout, "   %3s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                   "dLv", "NP_Skip_xL", "NP_Skip_xR", "NP_Skip_yL", "NP_Skip_yR", "NP_Skip_zL", "NP_Skip_zR" );

      for (int t=0; t<NLv; t++)
      Aux_Message( stdout, "   %3d  %10d  %10d  %10d  %10d  %10d  %10d\n",
                   t+1, RefineRegion[t][0], RefineRegion[t][1], RefineRegion[t][2],
                        RefineRegion[t][3], RefineRegion[t][4], RefineRegion[t][5] );
   }

} // FUNCTION : Load_RefineRegion



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_RefineRegion
// Description :  Flag patches within the refinement region specified by FlagPatch[]
//
// Note        :  1. Invoked by Init_ByFile()
//                2. Only useful when OPT__UM_IC_NLEVEL>1
//
// Parameter   :  lv        : Target AMR level
//                FlagPatch : Range of patches to be flagged
//                            --> For details, see the description under its declaration at the beginning
//                                of Init_ByFile()
//
// Return      :  amr->patch[0][lv][*]->flag
//-------------------------------------------------------------------------------------------------------
void Flag_RefineRegion( const int lv, const int FlagPatch[6] )
{

// initialize all flags as false
   for (int PID=0; PID<amr->num[lv]; PID++)  amr->patch[0][lv][PID]->flag = false;


// flag patches
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      int order[3]; // order of the target patch along each direction
      for (int d=0; d<3; d++)    order[d] = Mis_Scale2Cell( amr->patch[0][lv][PID]->corner[d], lv )/PS1;

      if ( order[0] >= FlagPatch[0]  &&  order[0] <= FlagPatch[1]  &&
           order[1] >= FlagPatch[2]  &&  order[1] <= FlagPatch[3]  &&
           order[2] >= FlagPatch[4]  &&  order[2] <= FlagPatch[5] )
         amr->patch[0][lv][PID]->flag = true;
   }

} // FUNCTION : Flag_RefineRegion
