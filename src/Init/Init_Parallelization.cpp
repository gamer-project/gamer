#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Parallelization
// Description :  Initialize parameters for parallelization
//
// Note        :  This function should be invoked even in the SERIAL mode
//-------------------------------------------------------------------------------------------------------
void Init_Parallelization()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// 0. validate the number of MPI ranks
//    --> check here first mainly because the particle routine in this function (step 6) needs to know "MPI_Rank"
#  ifdef SERIAL
   int NRank = 1;
#  else
   int NRank;
   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
#  endif

   if ( MPI_NRank != NRank )
      Aux_Error( ERROR_INFO, "MPI_NRank (%d) != MPI_Comm_size (%d) --> Is the runtime parameter MPI_NRANK consistent with mpirun?\n",
                 MPI_NRank, NRank );


// 1. number of parallel buffer cells for setting boundary conditions by either data copy or interpolation
// 1.1 get the number of ghost zones required in different interpolations
   int NSide_Useless, NGhost_Flu, NGhost_RefFlu;
   Int_Table( OPT__FLU_INT_SCHEME,     NSide_Useless, NGhost_Flu );
   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Useless, NGhost_RefFlu );

   int IntGhostSize_Flu;
   IntGhostSize_Flu = ( (FLU_GHOST_SIZE == 0) ? 0 : (FLU_GHOST_SIZE+1)/2 + NGhost_Flu );
   IntGhostSize_Flu = MAX( IntGhostSize_Flu, NGhost_RefFlu );

#  ifdef GRAVITY
   int NGhost_Pot, NGhost_Rho, NGhost_Gra, NGhost_RefPot;
   Int_Table( OPT__POT_INT_SCHEME,     NSide_Useless, NGhost_Pot );
   Int_Table( OPT__RHO_INT_SCHEME,     NSide_Useless, NGhost_Rho );
   Int_Table( OPT__GRA_INT_SCHEME,     NSide_Useless, NGhost_Gra );
   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Useless, NGhost_RefPot );

   int IntGhostSize_Pot, IntGhostSize_Rho;
   IntGhostSize_Pot = ( (POT_GHOST_SIZE == 0) ? 0 : (POT_GHOST_SIZE+1)/2 + NGhost_Pot );
// the following line may result in an unnecessarily large IntGhostSize_Pot if GRA_GHOST_SIZE>0 is used for STORE_POT_GHOST
   IntGhostSize_Pot = MAX( IntGhostSize_Pot, ( (GRA_GHOST_SIZE == 0) ? 0 : (GRA_GHOST_SIZE+1)/2 + NGhost_Gra ) );
   IntGhostSize_Pot = MAX( IntGhostSize_Pot, NGhost_RefPot );
   IntGhostSize_Rho = ( (RHO_GHOST_SIZE == 0) ? 0 : (RHO_GHOST_SIZE+1)/2 + NGhost_Rho );
#  endif

// 1.2 set the sizes of parallel buffers
   Flu_ParaBuf = MAX( FLU_GHOST_SIZE, IntGhostSize_Flu );
#  ifdef GRAVITY
   Pot_ParaBuf = MAX( GRA_GHOST_SIZE, IntGhostSize_Pot );
   Rho_ParaBuf = MAX( RHO_GHOST_SIZE, IntGhostSize_Rho );

   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   Flu_ParaBuf = MAX( Flu_ParaBuf,    IntGhostSize_Rho );   // ensure that the fluid ghost zone is large enough for Poisson
#  endif


// the following operations are useless for LOAD_BALANCE during restart
#  ifdef LOAD_BALANCE
   if ( OPT__INIT == INIT_BY_RESTART )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

      return;
   }
#  endif


// check
   if ( MPI_Rank == 0 )
   {
      if ( MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank )
         Aux_Error( ERROR_INFO, "MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank !!\n" );

      for (int d=0; d<3; d++)
      if ( MPI_NRank_X[d] <= 0 )    Aux_Error( ERROR_INFO, "MPI_NRank_X[%d] = %d <= 0 !!\n", d,  MPI_NRank_X[d] );
   }


// 2. number of coarse-grid cells in one rank for the rectangular domain decomposition
   for (int d=0; d<3; d++)    NX0[d] = NX0_TOT[d] / MPI_NRank_X[d];


// 3. MPI ranks in different spatial directions
   MPI_Rank_X[0] =  MPI_Rank%MPI_NRank_X[0];
   MPI_Rank_X[1] = (MPI_Rank/MPI_NRank_X[0]) % MPI_NRank_X[1];
   MPI_Rank_X[2] = (MPI_Rank/MPI_NRank_X[0]) / MPI_NRank_X[1];


// 4. sibling MPI ranks
   const int Buf     = 1;
   const int Size[3] = { MPI_NRank_X[0]+2*Buf, MPI_NRank_X[1]+2*Buf, MPI_NRank_X[2]+2*Buf };

   int ii, jj, kk, Width[3], Disp[3], ID1, ID2;
   int *RankMap = new int [ Size[0]*Size[1]*Size[2] ];

// interior ranks
   for (int k=Buf; k<Buf+MPI_NRank_X[2]; k++)   {  kk = k - Buf;
   for (int j=Buf; j<Buf+MPI_NRank_X[1]; j++)   {  jj = j - Buf;
   for (int i=Buf; i<Buf+MPI_NRank_X[0]; i++)   {  ii = i - Buf;

      ID1          = ( k *       Size[1] + j  )*       Size[0] + i;
      ID2          = ( kk*MPI_NRank_X[1] + jj )*MPI_NRank_X[0] + ii;
      RankMap[ID1] = ID2;
   }}}

// exterior ranks
   for (int s=0; s<26; s++)
   {
      for (int d=0; d<3; d++)
      {
         Width[d] = TABLE_01( s, 'x'+d, Buf, MPI_NRank_X[d], Buf                );
         Disp [d] = TABLE_01( s, 'x'+d,   0,            Buf, Buf+MPI_NRank_X[d] );
      }

      for (int k=Disp[2]; k<Disp[2]+Width[2]; k++) {  kk = k - Buf;
      for (int j=Disp[1]; j<Disp[1]+Width[1]; j++) {  jj = j - Buf;
      for (int i=Disp[0]; i<Disp[0]+Width[0]; i++) {  ii = i - Buf;

         ID1 = ( k*Size[1] + j )*Size[0] + i;

//       check the periodicity
         bool Periodic = true;
         for (int d=0; d<3; d++)
         {
            if (  ( TABLE_01(s,'x'+d, true,false,false) && OPT__BC_FLU[2*d+0] != BC_FLU_PERIODIC ) ||
                  ( TABLE_01(s,'x'+d,false,false, true) && OPT__BC_FLU[2*d+1] != BC_FLU_PERIODIC )   )
            {
               Periodic = false;
               break;
            }
         }

//       periodic BC
         if ( Periodic )
            RankMap[ID1] = ( kk + MPI_NRank_X[2] )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] +
                           ( jj + MPI_NRank_X[1] )%MPI_NRank_X[1]*MPI_NRank_X[0]                +
                           ( ii + MPI_NRank_X[0] )%MPI_NRank_X[0];

//       non-periodic BC
//       --> set to "SIB_OFFSET_NONPERIODIC-s"
//       --> similar to the sibling patches outside the simulation domain
         else
            RankMap[ID1] = SIB_OFFSET_NONPERIODIC - s;
      }}} // for k,j,i
   } // for (int s=0; s<26; s++)

// record the sibling ranks
   for (int s=0; s<26; s++)
   {
      for (int d=0; d<3; d++)    Disp[d] = TABLE_01( s, 'x'+d, -1, 0, +1 );

      ID1 = ( (Buf+MPI_Rank_X[2]+Disp[2])*Size[1] + (Buf+MPI_Rank_X[1]+Disp[1]) )*Size[0] + (Buf+MPI_Rank_X[0]+Disp[0]);
      MPI_SibRank[s] = RankMap[ID1];
   }

   delete [] RankMap;


// 5. left/right edges of each subdomain
#  ifndef SERIAL
   const double dh_min             = amr->dh[TOP_LEVEL];
   const int    SubDomain_Scale[3] = { amr->BoxScale[0]/MPI_NRank_X[0],
                                       amr->BoxScale[1]/MPI_NRank_X[1],
                                       amr->BoxScale[2]/MPI_NRank_X[2] };

   double (*SubDomain_EdgeL)[3] = new double [MPI_NRank][3];
   double (*SubDomain_EdgeR)[3] = new double [MPI_NRank][3];

   double (*SubDomain_EdgeL3D)[ MPI_NRank_X[1] ][ MPI_NRank_X[0] ][3]
                = ( double (*)[ MPI_NRank_X[1] ][ MPI_NRank_X[0] ][3] )SubDomain_EdgeL;
   double (*SubDomain_EdgeR3D)[ MPI_NRank_X[1] ][ MPI_NRank_X[0] ][3]
                = ( double (*)[ MPI_NRank_X[1] ][ MPI_NRank_X[0] ][3] )SubDomain_EdgeR;

// calculate the left/right edges by rank 0 only to ensure that all ranks see EXACTLY the same values (no round-off errors)
   if ( MPI_Rank == 0 )
   {
      int Rank_X[3];

//    set left edge
      for (int k=0; k<MPI_NRank_X[2]; k++)
      for (int j=0; j<MPI_NRank_X[1]; j++)
      for (int i=0; i<MPI_NRank_X[0]; i++)
      {
         SubDomain_EdgeL3D[k][j][i][0] = (double)(i*SubDomain_Scale[0])*dh_min;
         SubDomain_EdgeL3D[k][j][i][1] = (double)(j*SubDomain_Scale[1])*dh_min;
         SubDomain_EdgeL3D[k][j][i][2] = (double)(k*SubDomain_Scale[2])*dh_min;
      }

//    force EdgeR == EdgeL at the same sub-domain interfaces
      for (int k=0; k<MPI_NRank_X[2]; k++)
      for (int j=0; j<MPI_NRank_X[1]; j++)
      for (int i=0; i<MPI_NRank_X[0]; i++)
      {
         SubDomain_EdgeR3D[k][j][i][0] = ( i == MPI_NRank_X[0]-1 ) ? amr->BoxSize[0] : SubDomain_EdgeL3D[k  ][j  ][i+1][0];
         SubDomain_EdgeR3D[k][j][i][1] = ( j == MPI_NRank_X[1]-1 ) ? amr->BoxSize[1] : SubDomain_EdgeL3D[k  ][j+1][i  ][1];
         SubDomain_EdgeR3D[k][j][i][2] = ( k == MPI_NRank_X[2]-1 ) ? amr->BoxSize[2] : SubDomain_EdgeL3D[k+1][j  ][i  ][2];
      }
   } // if ( MPI_Rank == 0 )

// send results to all ranks
   MPI_Scatter( SubDomain_EdgeL, 3, MPI_DOUBLE, amr->ParaVar->SubDomain_EdgeL, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Scatter( SubDomain_EdgeR, 3, MPI_DOUBLE, amr->ParaVar->SubDomain_EdgeR, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   delete [] SubDomain_EdgeL;
   delete [] SubDomain_EdgeR;
#  endif // ifndef SERIAL


// 6. number of particles for each rank (only during the initialization)
#  ifdef PARTICLE
   if ( amr->Par->Init != PAR_INIT_BY_RESTART )
   {
      if ( amr->Par->NPar_Active_AllRank < 0 )
         Aux_Error( ERROR_INFO, "NPar_Active_AllRank = %ld < 0 !!\n", amr->Par->NPar_Active_AllRank );

      const int NPar_per_Rank  = amr->Par->NPar_Active_AllRank / MPI_NRank;
      const int Rank_with_more = amr->Par->NPar_Active_AllRank % MPI_NRank;

//    assuming all particles are active initially (some of them may be marked as inactive when calling Par_Aux_InitCheck)
      amr->Par->NPar_AcPlusInac = NPar_per_Rank;

      if ( MPI_Rank < Rank_with_more )    amr->Par->NPar_AcPlusInac ++;

//    check
#     ifdef DEBUG_PARTICLE
      long NPar_Sum = 0;
      MPI_Reduce( &amr->Par->NPar_AcPlusInac, &NPar_Sum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0  &&  NPar_Sum != amr->Par->NPar_Active_AllRank )
         Aux_Error( ERROR_INFO, "Total number of active particles in all ranks (%ld) != expected (%ld) !!\n",
                    NPar_Sum, amr->Par->NPar_Active_AllRank );
#     endif
   } // if ( OPT__INIT != INIT_BY_RESTART )
#  endif // #ifdef PARTICLE


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_Parallelization
