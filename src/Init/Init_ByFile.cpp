#include "GAMER.h"

static void UM_Downgrade( const real *Input, real *Output, const int NX, const int NY, const int NZ,
                          const int NVar );
static void UM_AssignData( const int lv, real *UM_Data, const int NVar );
static void UM_CreateLevel( real *UM_Data, const int lv, int **FlagMap, const int Buffer, const int NVar );
static void UM_RecordFlagMap( const int lv, int *FlagMap, const int ip, const int jp, const int kp );
static void UM_FindAncestor( const int Son_lv, const int Son_ip, const int Son_jp, const int Son_kp,
                             int **FlagMap );
static void UM_Flag( const int lv, const int *FlagMap );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile
// Description :  Set up the initial condition from an input uniform-mesh array
//
// Note        :  1. Create levels from 0 to OPT__UM_IC_LEVEL
//                   --> No levels above OPT__UM_IC_LEVEL will be created
//                       --> Unless OPT__UM_IC_REFINE is enabled, for which levels at OPT__UM_IC_LEVEL+1 ~ MAX_LEVEL
//                           will also be generated based on the given refinement criteria
//                   --> ALL patches at levels 0 to OPT__UM_IC_LEVEL will be created. In other words,
//                       the simulation domain will be fully refined to level OPT__UM_IC_LEVEL.
//                       --> But if OPT__UM_IC_DOWNGRADE, patches at levels 1~OPT__UM_IC_LEVEL
//                           may be removed if not satisfying the refinement criteria
//                2. The uniform-mesh input file should be named as "UM_IC"
//                3. This function can load any number of input values per cell (from 1 to NCOMP_TOTAL)
//                   --> Determined by the input parameter "OPT__UM_IC_NVAR"
//                   --> If "OPT__UM_IC_NVAR < NCOMP_TOTAL", one must specify the way to assign values to all
//                       variables in Model_Init_ByFile_AssignData()
//                4. The data format in the UM_IC file should be [k][j][i][v] instead of [v][k][j][i]
//                   --> Different from the data layout of patch->fluid
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ByFile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


   const char  FileName[]     = "UM_IC";
   const int   UM_lv          = OPT__UM_IC_LEVEL;
   const int   UM_NVar        = OPT__UM_IC_NVAR;
   const int   UM_Size_Tot[3] = {  NX0_TOT[0]*(1<<UM_lv),         // size of the input data
                                   NX0_TOT[1]*(1<<UM_lv),
                                   NX0_TOT[2]*(1<<UM_lv) };
   const int   UM_Size[3]     = {  UM_Size_Tot[0]/MPI_NRank_X[0], // size of the input data loaded by each rank
                                   UM_Size_Tot[1]/MPI_NRank_X[1],
                                   UM_Size_Tot[2]/MPI_NRank_X[2] };

   real *Input;               // store the input data
   real *UM_Data[UM_lv+1];    // downgraded variables at each level
   int  *FlagMap[UM_lv];      // record the positions of patches to be flagged at each level ( 1-->flag )



// initial check
// ===========================================================================================================
   if ( UM_lv < 0  ||  UM_lv > NLEVEL-1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d > NLEVEL-1 !!\n", "UM_lv", UM_lv );

   if ( UM_NVar < 1  ||  UM_NVar > NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "UM_NVar", UM_NVar );

   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *FileTemp = fopen( FileName, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long ExpectSize = long(UM_NVar)*UM_Size_Tot[0]*UM_Size_Tot[1]*UM_Size_Tot[2]*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                 FileName, FileSize, ExpectSize );

   fclose( FileTemp );

   MPI_Barrier( MPI_COMM_WORLD );



// set the file offset for this rank
// ===========================================================================================================
   const long offset0 = (long)UM_NVar * sizeof(real) * ( (long)MPI_Rank_X[2]*UM_Size[2]*UM_Size_Tot[1]*UM_Size_Tot[0]
                                                             + MPI_Rank_X[1]*UM_Size[1]*UM_Size_Tot[0]
                                                             + MPI_Rank_X[0]*UM_Size[0] );
   long offset, ID;



// load the initial condition from the input uniform-mesh data
// ===========================================================================================================
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading data ...\n" );

   Input = new real [ (long)UM_NVar*UM_Size[0]*UM_Size[1]*UM_Size[2] ];

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

         FILE *File = fopen( FileName, "rb" );

//       begin to load data
         for (int k=0; k<UM_Size[2]; k++)
         for (int j=0; j<UM_Size[1]; j++)
         {
            offset = offset0 + (long)UM_NVar*sizeof(real)*( (long)k*UM_Size_Tot[1]*UM_Size_Tot[0] + j*UM_Size_Tot[0] );
            ID     = (long)UM_NVar * ( (long)k*UM_Size[1]*UM_Size[0] + j*UM_Size[0] );

            fseek( File, offset, SEEK_SET );
            fread( &Input[ID], sizeof(real), UM_NVar*UM_Size[0], File );

//          verify that we do NOT exceed the file size
            if ( feof(File) )
               Aux_Error( ERROR_INFO, "rank %d reach the end of the file <%s> in the function <%s> !!\n",
                          MPI_Rank, FileName, __FUNCTION__ );
         }

         fclose( File );
      }

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)



// allocate and initialize the FlagMap array
// ===========================================================================================================
   for (int lv=0; lv<UM_lv; lv++)
   {
      const int NPatch1D[3] = { (NX0[0]/PATCH_SIZE)*(1<<lv) + 4,
                                (NX0[1]/PATCH_SIZE)*(1<<lv) + 4,
                                (NX0[2]/PATCH_SIZE)*(1<<lv) + 4  };
      const int NPatch = NPatch1D[0]*NPatch1D[1]*NPatch1D[2];

      FlagMap[lv] = new int [NPatch];

      for (int P=0; P<NPatch; P++)  FlagMap[lv][P] = -1;
   }



// compute the downgraded data at levels 0 -> UM_lv
// ===========================================================================================================
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Downgrading data ...\n" );

   UM_Data[UM_lv] = Input;

   for (int lv=UM_lv-1; lv>=0; lv--)
   {
      const int NX = NX0[0] * (1<<lv);
      const int NY = NX0[1] * (1<<lv);
      const int NZ = NX0[2] * (1<<lv);

      UM_Data[lv] = new real [ (long)UM_NVar*NX*NY*NZ ];

      UM_Downgrade( UM_Data[lv+1], UM_Data[lv], NX, NY, NZ, UM_NVar );
   }



// create level : 0
// ===========================================================================================================
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level 0 ...\n" );

   Init_BaseLevel();

#  ifdef PARTICLE
   Par_FindHomePatch_Base( BaseP );
#  endif

// assign data for the base level
   UM_AssignData( 0, UM_Data[0], UM_NVar );

// get the buffer data for the base level
   Buf_GetBufferData( 0, amr->FluSg[0], NULL_INT, DATA_GENERAL, _TOTAL, Flu_ParaBuf, USELB_NO );



// create level : UM_lv -> 1
// ===========================================================================================================
// construct the FlagMap for levels UM_lv-1 -> 0
   for (int lv=UM_lv; lv>0; lv--)
      UM_CreateLevel( UM_Data[lv], lv, FlagMap, (lv==MAX_LEVEL-1)?FLAG_BUFFER_SIZE_MAXM1_LV:
                                                (lv==MAX_LEVEL-2)?FLAG_BUFFER_SIZE_MAXM2_LV:
                                                                  FLAG_BUFFER_SIZE,
                      UM_NVar );

// flag levels 0 -> UM_lv-1 and construct levels 1-> UM_lv accordingly
// (we still have to loop over ALL levels since several lists needed to be initialized)
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level %d ...\n", lv+1 );

      if ( lv < UM_lv )  UM_Flag( lv, FlagMap[lv] );

      Buf_RecordBoundaryFlag( lv );

      MPI_ExchangeBoundaryFlag( lv );

      Flag_Buffer( lv );

      Init_Refine( lv );

      UM_AssignData( lv+1, UM_Data[lv+1], UM_NVar );

      Buf_GetBufferData( lv+1, amr->FluSg[lv+1], NULL_INT, DATA_GENERAL, _TOTAL, Flu_ParaBuf, USELB_NO );
   }

// get the total number of patches in all ranks
   for (int lv=0; lv<NLEVEL; lv++)     Mis_GetTotalPatchNumber( lv );

   delete [] Input;
   for (int lv=0; lv<UM_lv; lv++)   delete [] FlagMap[lv];
   for (int lv=0; lv<UM_lv; lv++)   delete [] UM_Data[lv];


// downgrade the uniform-mesh data from level=OPT__UM_IC_LEVEL to level=0
// ===========================================================================================================
   if ( OPT__UM_IC_DOWNGRADE )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Derefining the uniform-mesh data ...\n" );

      for (int lv=OPT__UM_IC_LEVEL-1; lv>=0; lv--)
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "      Lv %2d ... ", lv+1 );

         Flag_Real( lv, USELB_NO );

         MPI_ExchangeBoundaryFlag( lv );

         Flag_Buffer( lv );

         Refine( lv, USELB_NO );

         Buf_GetBufferData( lv+1, amr->FluSg[lv+1], NULL_INT, DATA_AFTER_REFINE, _TOTAL, Flu_ParaBuf, USELB_NO );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=OPT__UM_IC_LEVEL-1; lv>=0; lv--)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Downgrading the uniform-mesh data ... done\n" );
   } // if ( OPT__UM_IC_DOWNGRADE )


// refine the uniform-mesh data from level=OPT__UM_IC_LEVEL to level=MAX_LEVEL
// ===========================================================================================================
   if ( OPT__UM_IC_REFINE )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining the uniform-mesh data ...\n" );

      for (int lv=OPT__UM_IC_LEVEL; lv<MAX_LEVEL; lv++)
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "      Lv %2d ... ", lv );

         Flag_Real( lv, USELB_NO );

         MPI_ExchangeBoundaryFlag( lv );

         Flag_Buffer( lv );

         Refine( lv, USELB_NO );

         Buf_GetBufferData( lv+1, amr->FluSg[lv+1], NULL_INT, DATA_AFTER_REFINE, _TOTAL, Flu_ParaBuf, USELB_NO );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=OPT__UM_IC_LEVEL; lv<MAX_LEVEL; lv++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining the uniform-mesh data ... done\n" );
   } // if ( OPT__UM_IC_REFINE )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ByFile



//-------------------------------------------------------------------------------------------------------
// Function    :  UM_CreateLevel
// Description :  Create all patches at level lv
//
// Note        :  The flag buffer zones are also included
//
// Parameter   :  UM_Data : Input uniform-mesh array at level "lv"
//                lv      : Target refinement level to be constructed
//                FlagMap : Map recording the refinement flag of each patch
//                Buffer  : Size of the flag buffer
//                NVar    : Number of variables
//-------------------------------------------------------------------------------------------------------
void UM_CreateLevel( real *UM_Data, const int lv, int **FlagMap, const int Buffer, const int NVar )
{

   const int NPatch1D[3]    = { (NX0[0]/PATCH_SIZE)*(1<<(lv  )) + 4,
                                (NX0[1]/PATCH_SIZE)*(1<<(lv  )) + 4,
                                (NX0[2]/PATCH_SIZE)*(1<<(lv  )) + 4  };
   const int Fa_NPatch1D[3] = { (NX0[0]/PATCH_SIZE)*(1<<(lv-1)) + 4,
                                (NX0[1]/PATCH_SIZE)*(1<<(lv-1)) + 4,
                                (NX0[2]/PATCH_SIZE)*(1<<(lv-1)) + 4  };

   const int NX             = PATCH_SIZE * ( NPatch1D[0]-4 );
   const int NY             = PATCH_SIZE * ( NPatch1D[1]-4 );


// skip buffer patches
   for (int kp=2; kp<NPatch1D[2]-2; kp++)         // kp : (kp)th patch in the z direction
   for (int jp=2; jp<NPatch1D[1]-2; jp++)
   for (int ip=2; ip<NPatch1D[0]-2; ip++)
   {
//    const int RhoID0  = NVar * ( (kp-2)*PATCH_SIZE*NY*NX + (jp-2)*PATCH_SIZE*NX + (ip-2)*PATCH_SIZE );
      bool mark         = false;

      for (int k=0; k<PATCH_SIZE; k++)   {  if (mark) break;
      for (int j=0; j<PATCH_SIZE; j++)   {  if (mark) break;
      for (int i=0; i<PATCH_SIZE; i++)   {  if (mark) break;

//       const int RhoID = RhoID0 + NVar*(k*NY*NX + j*NX + i);

//       if ( UM_Data[RhoID] > MinRho )
         if ( true )    // always flag all patches at level <= OPT__UM_IC_LEVEL
         {
//          get the right index and corner of the patch 0 in the local patch group (8 patches = 1 patch group)
            const int ip2           = ip - ip%2;
            const int jp2           = jp - jp%2;
            const int kp2           = kp - kp%2;

            const int Fa_Corner[3]  = { ip2*PATCH_SIZE/2, jp2*PATCH_SIZE/2, kp2*PATCH_SIZE/2 };


//          evaluate the corner of the grid at level lv-1 enclosing the grid (i,j,k) at level lv
            const int i2            = i - i%2;
            const int j2            = j - j%2;
            const int k2            = k - k%2;
            const int Fa_FlagPos[3] = { (ip*PATCH_SIZE+i2)/2, (jp*PATCH_SIZE+j2)/2, (kp*PATCH_SIZE+k2)/2 };


//          evaluate the patch id in each direction including the buffer zone
            const int ip_start = ip2 - (  Fa_FlagPos[0]-Buffer <  Fa_Corner[0]  ?  2:0  );
            const int jp_start = jp2 - (  Fa_FlagPos[1]-Buffer <  Fa_Corner[1]  ?  2:0  );
            const int kp_start = kp2 - (  Fa_FlagPos[2]-Buffer <  Fa_Corner[2]  ?  2:0  );
            const int ip_end   = ip2 + (  Fa_FlagPos[0]+Buffer >= Fa_Corner[0]+PATCH_SIZE  ?  3:1  );
            const int jp_end   = jp2 + (  Fa_FlagPos[1]+Buffer >= Fa_Corner[1]+PATCH_SIZE  ?  3:1  );
            const int kp_end   = kp2 + (  Fa_FlagPos[2]+Buffer >= Fa_Corner[2]+PATCH_SIZE  ?  3:1  );


//          create patches including the buffer zone
            for ( int kp3=kp_start; kp3<kp_end; kp3+=2 )
            for ( int jp3=jp_start; jp3<jp_end; jp3+=2 )
            for ( int ip3=ip_start; ip3<ip_end; ip3+=2 )
            {
               const int Fa_ip = 1 + ip3/2;
               const int Fa_jp = 1 + jp3/2;
               const int Fa_kp = 1 + kp3/2;

               const int Fa_FlagID = Fa_kp*Fa_NPatch1D[1]*Fa_NPatch1D[0] + Fa_jp*Fa_NPatch1D[0] + Fa_ip;

               if ( FlagMap[lv-1][Fa_FlagID] == -1 )
               {
//                create new patches in level 0 ~ lv-1 (to ensure the proper-nesting condition)
                  UM_FindAncestor( lv, ip3, jp3, kp3, FlagMap );

//                record the position of patches being flagged at level "lv-1"
                  UM_RecordFlagMap( lv-1, FlagMap[lv-1], Fa_ip, Fa_jp, Fa_kp );
               }

            } // kp3, jp3, ip3

//          set mark = true to stop looping the root patch
            mark = true;

         } // if ( true )
      }}} // k, j, i
   } // kp, jp, ip
} // FUNCTION : UM_CreateLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  UM_FindAncestor
// Description :  For a given Son at level "Son_lv" with origin equal to "Son_ip, Son_jp, Son_kp",
//                find all its ancestors according to the proper-nesting condition
//
// Parameter   :  Son_lv  : Tefinement level of the son patch
//                Son_ip  : "Son_ip"-th son patch in the x direction
//                Son_jp  : "Son_jp"-th son patch in the y direction
//                Son_kp  : "Son_kp"-th son patch in the z direction
//                FlagMap : Map recording the refinement flag of each patch
//-------------------------------------------------------------------------------------------------------
void UM_FindAncestor( const int Son_lv, const int Son_ip, const int Son_jp, const int Son_kp, int **FlagMap )
{

   const int Fa_lv               = Son_lv - 1;
   const int GrandPa_NPatch1D[3] = { (NX0[0]/PATCH_SIZE)*(1<<(Fa_lv-1)) + 4,
                                     (NX0[1]/PATCH_SIZE)*(1<<(Fa_lv-1)) + 4,
                                     (NX0[2]/PATCH_SIZE)*(1<<(Fa_lv-1)) + 4  };

   const int Fa_ip               = 1 + Son_ip/2;
   const int Fa_jp               = 1 + Son_jp/2;
   const int Fa_kp               = 1 + Son_kp/2;

   const int Fa_ip2              = Fa_ip - Fa_ip%2;         // the ip of the patch 0 in the local patch group
   const int Fa_jp2              = Fa_jp - Fa_jp%2;
   const int Fa_kp2              = Fa_kp - Fa_kp%2;

   const int dip                 = 2*( -1 + 2*(Fa_ip%2) );  // the interval of ip of the sibling patch group
   const int djp                 = 2*( -1 + 2*(Fa_jp%2) );  // (even -> -2 ; odd -> +2)
   const int dkp                 = 2*( -1 + 2*(Fa_kp%2) );


// base-level patches will always exist
   if ( Fa_lv == 0 )    return;


// allcoate/find eight patch groups (each with eight patches) at level Fa_lv to surround the target son patch
   for (int Fa_kp3=Fa_kp2, kc=0; kc<2; Fa_kp3+=dkp, kc++)
   for (int Fa_jp3=Fa_jp2, jc=0; jc<2; Fa_jp3+=djp, jc++)
   for (int Fa_ip3=Fa_ip2, ic=0; ic<2; Fa_ip3+=dip, ic++)
   {
      const int GrandPa_ip = 1 + Fa_ip3/2;
      const int GrandPa_jp = 1 + Fa_jp3/2;
      const int GrandPa_kp = 1 + Fa_kp3/2;

      const int GrandPa_FlagID =   GrandPa_kp*GrandPa_NPatch1D[1]*GrandPa_NPatch1D[0]
                                 + GrandPa_jp*GrandPa_NPatch1D[0] + GrandPa_ip;


      if ( FlagMap[Fa_lv-1][GrandPa_FlagID] == -1 )
      {
//       create new patches in level 0 ~ Fa_lv-1 (to ensure the proper-nesting condition)
         UM_FindAncestor( Fa_lv, Fa_ip3, Fa_jp3, Fa_kp3, FlagMap );

//       record the position of patches being flagged in the level "Fa_lv-1"
         UM_RecordFlagMap( Fa_lv-1, FlagMap[Fa_lv-1], GrandPa_ip, GrandPa_jp, GrandPa_kp );
      }
   }

} // FUNCTION : UM_FindAncestor



//-------------------------------------------------------------------------------------------------------
// Function    :  UM_Downgrade
// Description :  Downgrade the data of Input array by a factor of two and store in the Output array
//
// Parameter   :  Input  : Input array
//                Output : Output array
//                NX     : Size of the Output array in the x direction
//                NY     : Size of the Output array in the y direction
//                NZ     : Size of the Output array in the z direction
//                NVar   : Number of variables
//-------------------------------------------------------------------------------------------------------
void UM_Downgrade( const real *Input, real *Output, const int NX, const int NY, const int NZ, const int NVar )
{

   const int NX2 = 2*NX;
   const int NY2 = 2*NY;
   int ii, jj, kk;

   for (int k=0; k<NZ; k++)   {  kk = 2*k;
   for (int j=0; j<NY; j++)   {  jj = 2*j;
   for (int i=0; i<NX; i++)   {  ii = 2*i;
   for (int v=0; v<NVar; v++) {

      const long OutID   =    (long)NVar * ( (long)k*NY*NX+ j*NX+ i ) + v;
      const long InID[8] = {  (long)NVar * ( (long)(kk+0)*NY2*NX2 + (jj+0)*NX2 + (ii+0) ) + v,
                              (long)NVar * ( (long)(kk+0)*NY2*NX2 + (jj+0)*NX2 + (ii+1) ) + v,
                              (long)NVar * ( (long)(kk+0)*NY2*NX2 + (jj+1)*NX2 + (ii+0) ) + v,
                              (long)NVar * ( (long)(kk+1)*NY2*NX2 + (jj+0)*NX2 + (ii+0) ) + v,
                              (long)NVar * ( (long)(kk+0)*NY2*NX2 + (jj+1)*NX2 + (ii+1) ) + v,
                              (long)NVar * ( (long)(kk+1)*NY2*NX2 + (jj+1)*NX2 + (ii+0) ) + v,
                              (long)NVar * ( (long)(kk+1)*NY2*NX2 + (jj+0)*NX2 + (ii+1) ) + v,
                              (long)NVar * ( (long)(kk+1)*NY2*NX2 + (jj+1)*NX2 + (ii+1) ) + v  };

      Output[OutID] = 0.125*(   Input[InID[0]] + Input[InID[1]] + Input[InID[2]] + Input[InID[3]]
                              + Input[InID[4]] + Input[InID[5]] + Input[InID[6]] + Input[InID[7]] );
   }}}}

} // FUNCTION : UM_Downgrade



//-------------------------------------------------------------------------------------------------------
// Function    :  UM_RecordFlagMap
// Description :  Record the postion of patches being flagged at level "lv"
//
// Note        :  +1 : Flagged
//                -1 : Unflagged
//
// Parameter   :  lv       : Target refinement level to be flagged
//                FlagMap  : Map recording the refinement flag of each patch
//                ip,jp,kp : Sequence numbers of patch in each direction
//-------------------------------------------------------------------------------------------------------
void UM_RecordFlagMap( const int lv, int *FlagMap, const int ip, const int jp, const int kp )
{

   const int NPatch1D[3] = { (NX0[0]/PATCH_SIZE)*(1<<lv) + 4,
                             (NX0[1]/PATCH_SIZE)*(1<<lv) + 4,
                             (NX0[2]/PATCH_SIZE)*(1<<lv) + 4  };

   const int FlagID      = kp*NPatch1D[1]*NPatch1D[0] + jp*NPatch1D[0] + ip;

   FlagMap[FlagID] = 1;

} // FUNCTION : UM_RecordFlagMap



//-------------------------------------------------------------------------------------------------------
// Function    :  UM_Flag
// Description :  Flag patches at level "lv" according to the FlagMap[lv]
//
// Parameter   :  lv      : Refinement level to be flagged
//                FlagMap : Map recording the refinement flag of each patch
//-------------------------------------------------------------------------------------------------------
void UM_Flag( const int lv, const int *FlagMap )
{

   const int scale0      = amr->scale[ 0];
   const int scale       = amr->scale[lv];
   const int NPatch1D[3] = { (NX0[0]/PATCH_SIZE)*(1<<lv) + 4,
                             (NX0[1]/PATCH_SIZE)*(1<<lv) + 4,
                             (NX0[2]/PATCH_SIZE)*(1<<lv) + 4  };

   int ip, jp, kp, FlagID;

   for (int PID=0; PID<amr->num[lv]; PID++)
   {
      ip       = ( amr->patch[0][lv][PID]->corner[0] - MPI_Rank_X[0]*NX0[0]*scale0 ) / (PATCH_SIZE*scale) + 2;
      jp       = ( amr->patch[0][lv][PID]->corner[1] - MPI_Rank_X[1]*NX0[1]*scale0 ) / (PATCH_SIZE*scale) + 2;
      kp       = ( amr->patch[0][lv][PID]->corner[2] - MPI_Rank_X[2]*NX0[2]*scale0 ) / (PATCH_SIZE*scale) + 2;

      FlagID   = kp*NPatch1D[1]*NPatch1D[0] + jp*NPatch1D[0] + ip;

      amr->patch[0][lv][PID]->flag = ( FlagMap[FlagID] == 1 ) ? true : false;
   }

} // FUNCTION : UM_Flag



//-------------------------------------------------------------------------------------------------------
// Function    :  UM_AssignData
// Description :  Use the input uniform-mesh array to assign data to all patches at level "lv"
//
// Note        :  If "NVar == NCOMP_TOTAL", we just copy the values recorded in UM_Data to all patches.
//                Otherwise, the model-dependent function XXX_Init_ByFile_AssignData() must be provided to
//                specify the way to assign data.
//
// Parameter   :  lv      : Target refinement level to assign data
//                UM_Data : Input uniform-mesh array
//                NVar    : Number of variables stored in UM_Data
//-------------------------------------------------------------------------------------------------------
void UM_AssignData( const int lv, real *UM_Data, const int NVar )
{

   if ( NVar == NCOMP_TOTAL )
   {
      const int NX[3]  = { NX0[0]*(1<<lv), NX0[1]*(1<<lv), NX0[2]*(1<<lv) };
      const int scale0 = amr->scale[ 0];
      const int scale  = amr->scale[lv];

      int *Corner;
      int ii, jj, kk;
      long Idx;

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         Corner = amr->patch[0][lv][PID]->corner;

         for (int k=0; k<PATCH_SIZE; k++)    { kk = ( Corner[2] - MPI_Rank_X[2]*NX0[2]*scale0 ) / scale + k;
         for (int j=0; j<PATCH_SIZE; j++)    { jj = ( Corner[1] - MPI_Rank_X[1]*NX0[1]*scale0 ) / scale + j;
         for (int i=0; i<PATCH_SIZE; i++)    { ii = ( Corner[0] - MPI_Rank_X[0]*NX0[0]*scale0 ) / scale + i;

            Idx = (long)NVar*( (long)kk*NX[1]*NX[0] + jj*NX[0]+ ii );

//          load all variables
            for (int v=0; v<NCOMP_TOTAL; v++)
               amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = UM_Data[ Idx + v ];

//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
               amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i]
                  = FMAX( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
            {
               real Passive[NCOMP_PASSIVE];

               for (int v=0; v<NCOMP_PASSIVE; v++)
                  Passive[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ NCOMP_FLUID + v ][k][j][i];

               CPU_NormalizePassive( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i],
                                     Passive, PassiveNorm_NVar, PassiveNorm_VarIdx );

               for (int v=0; v<NCOMP_PASSIVE; v++)
                  amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ NCOMP_FLUID + v ][k][j][i] = Passive[v];
            }
#           endif
         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // if ( NVar == NCOMP_TOTAL )

   else
   {
#     if   ( MODEL == HYDRO )
      Hydro_Init_ByFile_AssignData( lv, UM_Data, NVar );

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      ELBDM_Init_ByFile_AssignData( lv, UM_Data, NVar );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   } // if ( NVar == NCOMP_TOTAL ) ... else ...

} // FUNCTION : UM_AssignData
