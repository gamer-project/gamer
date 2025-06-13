#include "GAMER.h"


#ifdef SUPPORT_HYPRE
void Hypre_FillArrays( const int lv, const int NExtend, const double TimeNew, const real Poi_Coeff );
void Hypre_InitialGuess( const int lv );
void Hypre_SetA( const int lv );
void Hypre_SetBC( const int entry, int *stencil_indices, const int var, const int lv, const int *cornerL, const int *cornerR );
void Hypre_Solve( const Hypre_Solver_t Solver, int *N_iter, real *final_res_norm );



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_PrepareSingleLevel
// Description :  Prepare the single level Hypre arrays
//
// Parameter   :  lv      : Target level
//                NExtend : Number of cells to extend, mainly for Dirichlet boundary condition (NExtend = 1)
//-------------------------------------------------------------------------------------------------------
void Hypre_PrepareSingleLevel( const int lv, const int NExtend )
{

   const int  NParts      = 1; // number of AMR levels
   const int  NVars       = 1; // one var, potential
   const int  part        = 0; // one part only, no need to iterate
   const int  var         = 0; // one variable only, no need to iterate
   const int  NDim        = 3;
   const int  NEntries    = 2*NDim + 1;
   const int  object_type = HYPRE_SSTRUCT; // or this HYPRE_PARCSR // TODO: check which one should I use
   int Offsets[NEntries][NDim] = { {  0,  0,  0 },
                                   { -1,  0,  0 }, // -x
                                   {  1,  0,  0 }, // +x
                                   {  0, -1,  0 }, // -y
                                   {  0,  1,  0 }, // +y
                                   {  0,  0, -1 }, // -z
                                   {  0,  0,  1 }, // +z
                                 };

// create grid object
   HYPRE_CHECK_FUNC(   HYPRE_SStructGridCreate( HYPRE_MPI_COMM, NDim, NParts, &Hypre_grid )   );

// set each patch
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
#     ifdef DEBUG_HYPRE
      for (int d=0; d<3; d++)
      if ( amr->patch[0][lv][PID]->cornerR[d] - amr->patch[0][lv][PID]->cornerL[d] != PS1-1 )   Aux_Error( ERROR_INFO, "cornerR[%d]-cornerL[%d] != PS1-1 !!\n", d, d );
#     endif

      // TODO: Do not include the patch has son (not this level) ? If yes, need to update how to fill the array boundary
      HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetExtents( Hypre_grid, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR )   );
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// solve one cell-centered variables
   HYPRE_SStructVariable vartypes[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

   HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetVariables( Hypre_grid, part, NVars, vartypes )   );

// set periodic
   if ( OPT__BC_POT == BC_POT_PERIODIC )
   {
      int periodicity[3] = { NX0_TOT[0]*(int)(1L<<lv), NX0_TOT[1]*(int)(1L<<lv), NX0_TOT[2]*(int)(1L<<lv) };
      HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetPeriodic( Hypre_grid, part, periodicity )   );
      // Aux_Message( stdout, "HYPRE: periodic lv %d (%d %d %d)\n", lv, periodicity[0], periodicity[1], periodicity[2] );
   }

// assamble grid
   HYPRE_CHECK_FUNC(   HYPRE_SStructGridAssemble( Hypre_grid )   );

// create stencils
   HYPRE_CHECK_FUNC(   HYPRE_SStructStencilCreate( NDim, NEntries, &Hypre_stencil )   );

// set entries
   for (int e=0; e<NEntries; e++)
      HYPRE_CHECK_FUNC(   HYPRE_SStructStencilSetEntry( Hypre_stencil, e, Offsets[e], var )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructGraphCreate( HYPRE_MPI_COMM, Hypre_grid, &Hypre_graph )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructGraphSetObjectType( Hypre_graph, object_type )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructGraphSetStencil( Hypre_graph, part, var, Hypre_stencil )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructGraphAssemble( Hypre_graph )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixCreate( HYPRE_MPI_COMM, Hypre_graph, &Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorCreate( HYPRE_MPI_COMM, Hypre_grid,  &Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorCreate( HYPRE_MPI_COMM, Hypre_grid,  &Hypre_b )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetObjectType( Hypre_A, object_type )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetObjectType( Hypre_x, object_type )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetObjectType( Hypre_b, object_type )   );

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixInitialize( Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorInitialize( Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorInitialize( Hypre_b )   );

} // FUNCTION : Hypre_PrepareSingleLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_Free
// Description :  Destory the Hypre arrays
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Hypre_Free()
{

   HYPRE_CHECK_FUNC(   HYPRE_SStructGridDestroy( Hypre_grid )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructGraphDestroy( Hypre_graph )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructStencilDestroy( Hypre_stencil )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixDestroy( Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorDestroy( Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorDestroy( Hypre_b )   );

} // FUNCITON : Hypre_Free



#if ( defined GRAVITY  &&  POT_SCHEME == HYPRE_POI )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_SolvePoisson
// Description :  Solve the poisson equation
//
// Note        :  Use GetBoxValues more efficient then GetValues
//
// Parameter   :  SaveSg_Pot :
//                lv         : Target level
//-------------------------------------------------------------------------------------------------------
void Hypre_SolvePoisson( const int SaveSg_Pot, const int lv, const double TimeNew, const real Poi_Coeff )
{

   const int NExtend = 1;
   const int part    = 0; // single level only, no need to iterate parts
   const int var     = 0; // single variable only, no need to iterate variables
   int  N_iter;
   real final_residual;

   if ( NPatchTotal[lv] == 0 )   return;

   Hypre_PrepareSingleLevel( lv, NExtend );

// set Hypre arrays
   const bool   In_time = ( amr->PotSgTime[lv][SaveSg_Pot] == TimeNew  ||  amr->PotSgTime[lv][1-SaveSg_Pot] == TimeNew );
   const double Time_tmp = amr->PotSgTime[lv][SaveSg_Pot];
   if ( !In_time )  amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;
   Hypre_FillArrays( lv, NExtend, TimeNew, Poi_Coeff );
   if ( !In_time )  amr->PotSgTime[lv][SaveSg_Pot] = Time_tmp;

// setup solver and solve
   Hypre_Solve( HYPRE_SOLVER, &N_iter, &final_residual );

// collect the potential
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGather( Hypre_x )   );

// update GAMER array
   real *pote = new real [CUBE(PS1)];
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = IDX321( i, j, k, PS1, PS1 );
         amr->patch[ SaveSg_Pot ][lv][PID]->pot[k][j][i] = pote[idx];
      } // i, j, k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   delete [] pote;

   Hypre_Free();

// record
   Hypre_Aux_Record( "Poisson", lv, N_iter, final_residual );

// update Sg
   amr->PotSg    [lv]             = SaveSg_Pot;
   amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;

} // FUNCTION : Hypre_SolvePoisson
#endif // #if ( defined GRAVITY  &&  POT_SCHEME == HYPRE_POI )



void Hypre_FillArrays( const int lv, const int NExtend, const double TimeNew, const real Poi_Coeff )
{

   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const int    NDim              = 3;
   const int    NEntries          = 2*NDim + 1;
   const int    var               = 0;
   const int    part              = 0;
   const double dh                = amr->dh[lv];
   const double dh2               = SQR( dh );
   const double coeff             = - Poi_Coeff * dh2;

   int stencil_indices[NEntries];

   for (int i=0; i<NEntries; i++)  stencil_indices[i] = i;

// determin if the patch group is at boundary or course-fine edge
   int *PID0_List    = new int [ amr->NPatchComma[lv][1] / 8 ];
   int *PID0_BC_List = new int [ amr->NPatchComma[lv][1] / 8 ];
   int  NPG_BC       = 0;
   int  NPG          = 0;

   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
      bool atBC = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         for (int s=0; s<6; s++)
         {
            const int d     = s/2;
            const int TDir1 = (d+1) % 3;
            const int TDir2 = (d+2) % 3;
            const int SibPID = amr->patch[0][lv][PID]->sibling[s];

            if ( SibPID >= 0 )   continue;
            if ( SibPID < SIB_OFFSET_NONPERIODIC  &&  OPT__BC_POT == BC_POT_PERIODIC )   continue;

            atBC = true;
            break;
         }
         if ( atBC )   break;
      } // for (int PID=PID0; PID<PID0+8; PID++)

      if ( atBC )
      {
         PID0_BC_List[NPG_BC] = PID0;
         NPG_BC++;
      }
      PID0_List[NPG] = PID0;
      NPG++;
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1] / 8; PID0+=8)

   real   *Matrix_Laplace = new real [NEntries*CUBE(PS1)];
   real   *pote           = new real [CUBE(PS1)];
   real   *dens           = new real [CUBE(PS1)];
   real   *bcBox          = new real [1*SQR(PS1)];
   real   *bcVal          = new real [SQR(PS1)];
   real  (*Pot_Array)[CUBE(PS1+2)] = NULL;
   real  (*Dens_Array)[CUBE(PS1)]  = NULL;

// fill common Laplace matrix
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
      const int idx = IDX321( i, j, k, PS1, PS1 ) * NEntries;

      Matrix_Laplace[idx+0] =  6; // i,   j,   k
      Matrix_Laplace[idx+1] = -1; // i-1, j,   k
      Matrix_Laplace[idx+2] = -1; // i+1, j,   k
      Matrix_Laplace[idx+3] = -1; // i,   j-1, k
      Matrix_Laplace[idx+4] = -1; // i,   j+1, k
      Matrix_Laplace[idx+5] = -1; // i,   j,   k-1
      Matrix_Laplace[idx+6] = -1; // i,   j,   k+1
   }

// prepare density with particle density
   Dens_Array = new real [8*NPG][ CUBE(PS1) ];

   Prepare_PatchData( lv, TimeNew, Dens_Array[0], NULL,
                      0, NPG, PID0_List, _TOTAL_DENS, _NONE, OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00,
                      IntPhase_No, OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

// fill patches
#  ifdef COMOVING
   const bool Comoving = true;
#  else
   const bool Comoving = false;
#  endif
   real RhoSubtract;
   if ( OPT__BC_POT == BC_POT_PERIODIC ) RhoSubtract = (real)AveDensity_Init;
   else if ( Comoving )                  RhoSubtract = (real)1.0;
   else                                  RhoSubtract = (real)0.0;

   for (int PG=0; PG<NPG; PG++)
   for (int LocalID=0; LocalID<8; LocalID++)
   {
      const int PID = PID0_List[PG] + LocalID;

      for (int k=0; k<PS1; k++) { const double z = amr->patch[0][lv][PID]->EdgeL[2] + (0.5+k)*dh;
      for (int j=0; j<PS1; j++) { const double y = amr->patch[0][lv][PID]->EdgeL[1] + (0.5+j)*dh;
      for (int i=0; i<PS1; i++) { const double x = amr->patch[0][lv][PID]->EdgeL[0] + (0.5+i)*dh;
         const int idx = IDX321( i, j, k, PS1, PS1 );

         real Dens = Dens_Array[8*PG+LocalID][idx] - RhoSubtract;

//       add extra mass source for gravity if required
         if ( OPT__GRAVITY_EXTRA_MASS )
            Dens += Poi_AddExtraMassForGravity_Ptr( x, y, z, Time[lv], lv, NULL );

         dens[idx] = coeff * Dens;
         pote[idx] = 0.0;
      }}} // i, j, k

      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, dens )   );
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );
      HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, NEntries, stencil_indices, Matrix_Laplace )   );
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// prepare boundary potential
   if ( NPG_BC > 0 )
   {
      Pot_Array = new real [8*NPG_BC][ CUBE(PS1+2) ];

      Prepare_PatchData( lv, TimeNew, Pot_Array[0], NULL,
                         NExtend, NPG_BC, PID0_BC_List, _POTE, _NONE, OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_06,
                         IntPhase_No, OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
   }

// fill boundary condition
   for (int PG=0; PG<NPG_BC; PG++)
   for (int LocalID=0; LocalID<8; LocalID++)
   {
      const int PID = PID0_BC_List[PG] + LocalID;
      int cornerL[3] = { amr->patch[0][lv][PID]->cornerL[0], amr->patch[0][lv][PID]->cornerL[1], amr->patch[0][lv][PID]->cornerL[2] };
      int cornerR[3] = { amr->patch[0][lv][PID]->cornerR[0], amr->patch[0][lv][PID]->cornerR[1], amr->patch[0][lv][PID]->cornerR[2] };

//    remove the connection and update A due to Dirichlet boundary condition
      for (int s=0; s<6; s++)
      {
         const int SibPID = amr->patch[0][lv][PID]->sibling[s];

         if ( SibPID >= 0 )   continue;
         if ( SibPID < SIB_OFFSET_NONPERIODIC  &&  OPT__BC_POT == BC_POT_PERIODIC )   continue;

         const int d     = s/2;
         const int TDir1 = (d+1) % 3;
         const int TDir2 = (d+2) % 3;
         const int lr    = s%2;

         int stencil_indices_bc[1] = { s+1 };
         int cornerL_bc        [3] = { cornerL[0], cornerL[1], cornerL[2] };
         int cornerR_bc        [3] = { cornerR[0], cornerR[1], cornerR[2] };

         if ( lr )   cornerL_bc[d] = cornerR_bc[d];
         else        cornerR_bc[d] = cornerL_bc[d];

         const int Nk = cornerR_bc[2] - cornerL_bc[2] + 1;
         const int Nj = cornerR_bc[1] - cornerL_bc[1] + 1;
         const int Ni = cornerR_bc[0] - cornerL_bc[0] + 1;

         HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_b, part, cornerL_bc, cornerR_bc, var, bcVal )   );

         for (int k=0; k<Nk; k++)
         for (int j=0; j<Nj; j++)
         for (int i=0; i<Ni; i++)
         {
            const int idx = IDX321( i, j, k, Ni, Nj );
            int idx_i, idx_j, idx_k;
            switch ( d )
            {
               case 0: idx_i = ( lr ) ? PS1+2*NExtend-1 : 0;  idx_j = j + NExtend;                   idx_k = k + NExtend;                   break;
               case 1: idx_i = i + NExtend;                   idx_j = ( lr ) ? PS1+2*NExtend-1 : 0;  idx_k = k + NExtend;                   break;
               case 2: idx_i = i + NExtend;                   idx_j = j + NExtend;                   idx_k = ( lr ) ? PS1+2*NExtend-1 : 0;  break;
            } // switch ( d )
            const int idx_arr = idx_k*SQR(PS1+2*NExtend) + idx_j*(PS1+2*NExtend) + idx_i;

            const double temp = bcVal[idx];
            bcBox[idx]  = 0.0;
            bcVal[idx] += Pot_Array[8*PG+LocalID][idx_arr];
         } // i, j, k

         HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, cornerL_bc, cornerR_bc, var, bcVal )   );
         HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, cornerL_bc, cornerR_bc, var, 1, stencil_indices_bc, bcBox )   );
      } // for (int s=0; s<6; s++)
   } // for (int PG=0; PG<NPG_BC; PG++); for (int LocalID=0; LocalID<8; LocalID++)

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixAssemble( Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_b )   );

   delete [] Matrix_Laplace;
   delete [] pote;
   delete [] dens;
   delete [] bcBox;
   delete [] bcVal;
   delete [] PID0_BC_List;
   delete [] PID0_List;
   delete [] Pot_Array;
   delete [] Dens_Array;

} // FUNCITON : Hypre_FillArrays



#endif // #ifdef SUPPORT_HYPRE
