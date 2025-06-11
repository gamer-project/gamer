#include "GAMER.h"


// static bool hypre_debug = true;
static bool hypre_debug = false;

extern bool in_flu_corr;


#ifdef SUPPORT_HYPRE
void Hypre_FillArrays( const int lv, const int NExtend, const double TimeNew, const real Poi_Coeff );
void Hypre_InitialGuess( const int lv );
void Hypre_SetA( const int lv );
void Hypre_SetBC( const int entry, int *stencil_indices, const int var, const int lv, const int *cornerL, const int *cornerR );
void Hypre_Solve( const Hypre_Solver_t Solver, int *N_iter, real *final_res_norm );


//      Aux_Message( stdout, "DEBUG: %s %d\n", __FILE__, __LINE__ );
   // Aux_Message( stdout, "\n" );
   // Aux_Message( stdout, "Track flu prepare %d %24.16e %d %24.16e %s %d\n", amr->FluSg[lv], amr->FluSgTime[lv][amr->FluSg[lv]], 1-amr->FluSg[lv], amr->FluSgTime[lv][1-amr->FluSg[lv]], __FUNCTION__, __LINE__ );
   // Aux_Message( stdout, "Track pot prepare %d %24.16e %d %24.16e %s %d\n", amr->PotSg[lv], amr->PotSgTime[lv][amr->PotSg[lv]], 1-amr->PotSg[lv], amr->PotSgTime[lv][1-amr->PotSg[lv]], __FUNCTION__, __LINE__ );

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
                                 }; // TODO: need to be checked

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
// Parameter   :  SaveSg_Pot :
//                lv         : Target level
//-------------------------------------------------------------------------------------------------------
void Hypre_SolvePoisson( const int SaveSg_Pot, const int lv, const double TimeNew, const real Poi_Coeff )
{

   const int NExtend = 1;
   const int part    = 0; // single level only, no need to iterate parts
   const int var     = 0; // single variable only, no need to iterate variables

   if ( NPatchTotal[lv] == 0 )   return;

   Hypre_PrepareSingleLevel( lv, NExtend );

// set Hypre arrays
   // Aux_Message( stdout, "%s lv: %d, TimeNew: %24.16e\n", __FUNCTION__, lv, TimeNew );

   const bool   In_time = ( amr->PotSgTime[lv][SaveSg_Pot] == TimeNew  ||  amr->PotSgTime[lv][1-SaveSg_Pot] == TimeNew );
   const double Time_tmp = amr->PotSgTime[lv][SaveSg_Pot];
   if ( !In_time )  amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;
   Hypre_FillArrays( lv, NExtend, TimeNew, Poi_Coeff );
   if ( !In_time )  amr->PotSgTime[lv][SaveSg_Pot] = Time_tmp;

// setup solver and solve
   int N_iter;
   real final_residual;
   Hypre_Solve( HYPRE_SOLVER, &N_iter, &final_residual );
   Hypre_Aux_Record( "Poisson", lv, N_iter, final_residual );

// collect the potential
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGather( Hypre_x )   );

// update GAMER array
   real *pote = new real [CUBE(PS1)];
   FILE *f;
   if (hypre_debug && ! in_flu_corr)
   {
      char filename[MAX_STRING];
      sprintf( filename, "Hypre_R%d", MPI_Rank );
      f = fopen( filename, "a" );
   }
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      // NOTE: By FLASH: Use GetBoxValues more efficient then GetValues.
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = IDX321( i, j, k, PS1, PS1 );
         amr->patch[ SaveSg_Pot ][lv][PID]->pot[k][j][i] = pote[idx];
         if (hypre_debug && ! in_flu_corr)
         {
            Aux_Message( stdout, "FILL X Rank: %d, PID: %6d, cell: (%4d %4d %4d), val: %24.16e\n",
                                  MPI_Rank, PID,
                                  amr->patch[0][lv][PID]->cornerL[0] + i,
                                  amr->patch[0][lv][PID]->cornerL[1] + j,
                                  amr->patch[0][lv][PID]->cornerL[2] + k,
                                  pote[idx]
                                 );
            fprintf( f, "FILL X LV: %d, Step: %6ld, PID: %6d, cell: (%4d %4d %4d), val: %24.16e\n",
                         lv, AdvanceCounter[lv], PID,
                         amr->patch[0][lv][PID]->cornerL[0] + i,
                         amr->patch[0][lv][PID]->cornerL[1] + j,
                         amr->patch[0][lv][PID]->cornerL[2] + k,
                         pote[idx]
                        );
         }
         // TODO: apply the floor value here
      } // i, j, k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   {
      const int  NDim              = 3;
      const int  NEntries          = 2*NDim + 1;

      int stencil_indices[NEntries];

      for (int i=0; i<NEntries; i++)  stencil_indices[i] = i;

      real   *Matrix_Laplace = new real [NEntries*CUBE(PS1)];
      real   *dens           = new real [CUBE(PS1)];
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_b, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, dens )   );
         HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );
         HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixGetBoxValues( Hypre_A, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, NEntries, stencil_indices, Matrix_Laplace )   );
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            const int idx = IDX321( i, j, k, PS1, PS1 );

            // if (hypre_debug)
            if (hypre_debug && ! in_flu_corr)
            {
               Aux_Message( stdout, "OUTHYPRE Rank: %d, PID: %6d, cell: (%4d %4d %4d), dens: %24.16e, pote: %24.16e, A: %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                                     MPI_Rank, PID,
                                     amr->patch[0][lv][PID]->cornerL[0] + i,
                                     amr->patch[0][lv][PID]->cornerL[1] + j,
                                     amr->patch[0][lv][PID]->cornerL[2] + k,
                                     dens[idx], pote[idx],
                                     Matrix_Laplace[7*idx+0],
                                     Matrix_Laplace[7*idx+1],
                                     Matrix_Laplace[7*idx+2],
                                     Matrix_Laplace[7*idx+3],
                                     Matrix_Laplace[7*idx+4],
                                     Matrix_Laplace[7*idx+5],
                                     Matrix_Laplace[7*idx+6]
                                    );
               fprintf( f, "OUTHYPRE Rank: %d, PID: %6d, cell: (%4d %4d %4d), dens: %24.16e, pote: %24.16e, A: %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                            MPI_Rank, PID,
                            amr->patch[0][lv][PID]->cornerL[0] + i,
                            amr->patch[0][lv][PID]->cornerL[1] + j,
                            amr->patch[0][lv][PID]->cornerL[2] + k,
                            dens[idx], pote[idx],
                            Matrix_Laplace[7*idx+0],
                            Matrix_Laplace[7*idx+1],
                            Matrix_Laplace[7*idx+2],
                            Matrix_Laplace[7*idx+3],
                            Matrix_Laplace[7*idx+4],
                            Matrix_Laplace[7*idx+5],
                            Matrix_Laplace[7*idx+6]
                           );
            }
         }

      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] Matrix_Laplace;
      delete [] dens;
   }
   if (hypre_debug && ! in_flu_corr)   fclose(f);

   delete [] pote;

   Hypre_Free();

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

// prepare boundary potential
   if ( NPG_BC > 0 )
   {
      Pot_Array = new real [8*NPG_BC][ CUBE(PS1+2) ];

      Prepare_PatchData( lv, TimeNew, Pot_Array[0], NULL,
                         NExtend, NPG_BC, PID0_BC_List, _POTE, _NONE, OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_06,
                         IntPhase_No, OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
   }

// fill patches
   FILE *f;
   if (hypre_debug && ! in_flu_corr)
   {
      char filename[MAX_STRING];
      sprintf( filename, "Hypre_R%d", MPI_Rank );
      f = fopen( filename, "a" );
   }

   real RhoSubtract;
#  ifdef COMOVING
   const bool Comoving = true;
#  else
   const bool Comoving = false;
#  endif
   if ( OPT__BC_POT == BC_POT_PERIODIC ) RhoSubtract = (real)AveDensity_Init;
   else if ( Comoving )                  RhoSubtract = (real)1.0;
   else                                  RhoSubtract = (real)0.0;

   for (int PG=0; PG<NPG; PG++)
   for (int LocalID=0; LocalID<8; LocalID++)
   {
      const int PID = PID0_List[PG] + LocalID;
      // Aux_Message( stdout, "cornerL: (%6d %6d %6d)\n", amr->patch[0][lv][PID]->cornerL[0], amr->patch[0][lv][PID]->cornerL[1], amr->patch[0][lv][PID]->cornerL[2] );
      for (int k=0; k<PS1; k++) { const double z = amr->patch[0][lv][PID]->EdgeL[2] + (0.5+k)*dh;
      for (int j=0; j<PS1; j++) { const double y = amr->patch[0][lv][PID]->EdgeL[1] + (0.5+j)*dh;
      for (int i=0; i<PS1; i++) { const double x = amr->patch[0][lv][PID]->EdgeL[0] + (0.5+i)*dh;
         const int idx = IDX321( i, j, k, PS1, PS1 );

         // real Dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] - RhoSubtract;
         real Dens = Dens_Array[8*PG+LocalID][idx] - RhoSubtract;

//       add extra mass source for gravity if required
         if ( OPT__GRAVITY_EXTRA_MASS )
            Dens += Poi_AddExtraMassForGravity_Ptr( x, y, z, Time[lv], lv, NULL );

         dens[idx] = coeff * Dens;
         pote[idx] = 0.0;

         // if (hypre_debug)
         if (hypre_debug && ! in_flu_corr)
         {
            Aux_Message( stdout, "FILL B Rank: %d, PID: %6d, cell: (%4d %4d %4d), val: %24.16e\n",
                                  MPI_Rank, PID,
                                  amr->patch[0][lv][PID]->cornerL[0] + i,
                                  amr->patch[0][lv][PID]->cornerL[1] + j,
                                  amr->patch[0][lv][PID]->cornerL[2] + k,
                                  dens[idx]
                                 );
            fprintf( f, "FILL B LV: %d, Step: %6ld, PID: %6d, cell: (%4d %4d %4d), val: %24.16e\n",
                         lv, AdvanceCounter[lv], PID,
                         amr->patch[0][lv][PID]->cornerL[0] + i,
                         amr->patch[0][lv][PID]->cornerL[1] + j,
                         amr->patch[0][lv][PID]->cornerL[2] + k,
                         dens[idx]
                        );
         }
      }}}

      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, dens )   );
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );
      HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, NEntries, stencil_indices, Matrix_Laplace )   );
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   if (hypre_debug && ! in_flu_corr)   fclose(f);

// fill boundary condition
   if (hypre_debug && ! in_flu_corr)
   {
      char filename[MAX_STRING];
      sprintf( filename, "Hypre_R%d", MPI_Rank );
      f = fopen( filename, "a" );
   }

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

            // if (hypre_debug)
            if (hypre_debug && ! in_flu_corr)
            {
               Aux_Message( stdout, "FILL A Rank: %d, PID: %6d, s: %d, cell: (%4d %4d %4d), val: %24.16e -> %24.16e\n",
                                     MPI_Rank, PID, s+1,
                                     cornerL_bc[0] + i,
                                     cornerL_bc[1] + j,
                                     cornerL_bc[2] + k,
                                     temp, bcVal[idx]
                                    );

               fprintf( f, "FILL A LV: %d, Step: %6ld, PID: %6d, s: %d, cell: (%4d %4d %4d), val: %24.16e -> %24.16e\n",
                            lv, AdvanceCounter[lv], PID, s+1,
                            cornerL_bc[0] + i,
                            cornerL_bc[1] + j,
                            cornerL_bc[2] + k,
                            temp, bcVal[idx]
                           );
            }
         } // i, j, k

         HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, cornerL_bc, cornerR_bc, var, bcVal )   );
         HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, cornerL_bc, cornerR_bc, var, 1, stencil_indices_bc, bcBox )   );
      } // for (int s=0; s<6; s++)
   } // for (int PG=0; PG<NPG_BC; PG++); for (int LocalID=0; LocalID<8; LocalID++)
   if (hypre_debug && ! in_flu_corr)   fclose(f);

   // delete [] Matrix_Laplace;
   // delete [] pote;
   // delete [] dens;
   // delete [] bcBox;
   // delete [] bcVal;
   // delete [] PID0_BC_List;
   // delete [] Pot_Array;

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixAssemble( Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_b )   );

// check
   // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   // {
   //    HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_b, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, dens )   );
   //    HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );
   //    HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixGetBoxValues( Hypre_A, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, NEntries, stencil_indices, Matrix_Laplace )   );
   //    for (int k=0; k<PS1; k++)
   //    for (int j=0; j<PS1; j++)
   //    for (int i=0; i<PS1; i++)
   //    {
   //       const int idx = IDX321( i, j, k, PS1, PS1 );

   //       // if (hypre_debug)
   //       if (hypre_debug && ! in_flu_corr)
   //       {
   //          Aux_Message( stdout, "OUTHYPRE Rank: %d, PID: %6d, cell: (%4d %4d %4d), dens: %24.16e, pote: %24.16e, A: %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
   //                                MPI_Rank, PID,
   //                                amr->patch[0][lv][PID]->cornerL[0] + i,
   //                                amr->patch[0][lv][PID]->cornerL[1] + j,
   //                                amr->patch[0][lv][PID]->cornerL[2] + k,
   //                                dens[idx], pote[idx],
   //                                Matrix_Laplace[7*idx+0],
   //                                Matrix_Laplace[7*idx+1],
   //                                Matrix_Laplace[7*idx+2],
   //                                Matrix_Laplace[7*idx+3],
   //                                Matrix_Laplace[7*idx+4],
   //                                Matrix_Laplace[7*idx+5],
   //                                Matrix_Laplace[7*idx+6]
   //                               );
   //       }
   //    }

   // } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

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



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_InitialGuess
// Description :  Set the initial guess
//
// Parameter   :  lv : Target level
//-------------------------------------------------------------------------------------------------------
void Hypre_InitialGuess( const int lv )
{

   const int var  = 0;
   const int part = 0;

   real *pote = new real [CUBE(PS1)];
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         // TODO: check the data structure
         const int idx = k*SQR(PS1) + j*PS1 + i;
         pote[idx] = 0.0;
      }

      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );

   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   delete [] pote;

   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_x )   );

} // FUNCTION : Hypre_InitialGuess



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_SetA
// Description :  Set the array A
//
// Parameter   :  lv : Target level
//-------------------------------------------------------------------------------------------------------
void Hypre_SetA( const int lv )
{

   const int NDim     = 3;
   const int NEntries = 2*NDim + 1;
   const int var      = 0;
   const int part     = 0;
   const bool fix_one_cell_sol = ( OPT__BC_POT != BC_POT_ISOLATED );

   int stencil_indices[NEntries];

   for (int i=0; i<NEntries; i++)  stencil_indices[i] = i;

   real *boxVal = new real [NEntries*CUBE(PS1)];
   real *bcVal  = new real [SQR(PS1)];

   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = (k*SQR(PS1) + j*PS1 + i) * NEntries;
//       \nabla matrix
         boxVal[idx+0] =  6; // i,   j,   k
         boxVal[idx+1] = -1; // i-1, j,   k
         boxVal[idx+2] = -1; // i+1, j,   k
         boxVal[idx+3] = -1; // i,   j-1, k
         boxVal[idx+4] = -1; // i,   j+1, k
         boxVal[idx+5] = -1; // i,   j,   k-1
         boxVal[idx+6] = -1; // i,   j,   k+1
      }

      HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR,
                                                           var, NEntries, stencil_indices, boxVal )   );

//    set BC
//       if ( amr->patch[0][lv][PID]->EdgeL[0] == 0.0             )   Hypre_SetBC( 1, stencil_indices, var, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR );
//       if ( amr->patch[0][lv][PID]->EdgeR[0] == amr->BoxSize[0] )   Hypre_SetBC( 2, stencil_indices, var, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR );
//       if ( amr->patch[0][lv][PID]->EdgeL[1] == 0.0             )   Hypre_SetBC( 3, stencil_indices, var, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR );
//       if ( amr->patch[0][lv][PID]->EdgeR[1] == amr->BoxSize[1] )   Hypre_SetBC( 4, stencil_indices, var, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR );
//       if ( amr->patch[0][lv][PID]->EdgeL[2] == 0.0             )   Hypre_SetBC( 5, stencil_indices, var, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR );
//       if ( amr->patch[0][lv][PID]->EdgeR[2] == amr->BoxSize[2] )   Hypre_SetBC( 6, stencil_indices, var, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR );
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   delete [] boxVal;
   delete [] bcVal;

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixAssemble( Hypre_A )   );

} // FUNCTION : Hypre_SetA



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_SetBC
// Description :  Set boundary condition
//
// Parameter   :  entry           :
//                stencil_indices :
//                var             :
//                lv              :
//                cornerL         :
//                cornerR         :
//-------------------------------------------------------------------------------------------------------
// dir: 1~6 -> 1:-x, 2:+x, 3:-y, 4:+y, 5:-z, 6:+z
void Hypre_SetBC( const int entry, int *stencil_indices, const int var, const int lv, const int *cornerL, const int *cornerR )
{

   int hypre_ierr;
   int bcCornerL[3], bcCornerR[3];

   for (int d=0; d<3; d++)
   {
      bcCornerL[d] = cornerL[d];
      bcCornerR[d] = cornerR[d];
   }

   switch ( entry )
   {
      case 1: bcCornerR[0] = cornerL[0];  break;
      case 2: bcCornerL[0] = cornerR[0];  break;
      case 3: bcCornerR[1] = cornerL[1];  break;
      case 4: bcCornerL[1] = cornerR[1];  break;
      case 5: bcCornerR[2] = cornerL[2];  break;
      case 6: bcCornerL[2] = cornerR[2];  break;
      default:  Aux_Error( ERROR_INFO, "Incorrect entry (1~6)\n" ); break;
   }

   real *bcVal  = new real [SQR(PS1)];

   if ( OPT__BC_POT == BC_POT_ISOLATED )
   {
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = j*PS1 + i;
         bcVal[idx] = 0;
      }

      HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, lv, bcCornerL, bcCornerR, var, 1, stencil_indices+entry, bcVal )   );
   }

   delete [] bcVal;

} // FUNCITON : Hypre_SetBC



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_FillVector
// Description :  Fill Hypre vector
//
// Parameter   :  Hypre_verctor :
//                field_idx     :
//                factor        :
//-------------------------------------------------------------------------------------------------------
void Hypre_FillVector( HYPRE_SStructVector Hypre_vector, const int field_idx, const real factor )
{

   int var = 0;

   real *field = new real [CUBE(PS1)];
   for (int lv=0; lv<MAX_LEVEL+1; lv++)
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            // TODO: check the data structure
            const int idx = k*SQR(PS1) + j*PS1 + i;
            field[idx] = factor * amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[field_idx][k][j][i];
         }

         HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_vector, lv, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, field )   );

      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<MAX_LEVEL+1; lv++)

   delete [] field;

} // FUNCTION : Hypre_FillVector
#endif // #ifdef SUPPORT_HYPRE
