#include "GAMER.h"




#ifdef SUPPORT_HYPRE
void Hypre_FillArrays( const int lv, const int NExtend, const double TimeNew );
void Hypre_InitialGuess( const int lv );
void Hypre_SetA( const int lv );
void Hypre_SetBC( const int entry, int *stencil_indices, const int var, const int lv, const int *cornerL, const int *cornerR );
void Hypre_Solve();



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
   const int  NVars       = 1; // One var, potential?
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
      if ( amr->patch[0][lv][PID]->cornerR[0] - amr->patch[0][lv][PID]->cornerL[0] != PS1-1 )   Aux_Error( ERROR_INFO, "cornerR[0]-cornerL[0] != PS1-1 !!\n" );
      if ( amr->patch[0][lv][PID]->cornerR[1] - amr->patch[0][lv][PID]->cornerL[1] != PS1-1 )   Aux_Error( ERROR_INFO, "cornerR[1]-cornerL[1] != PS1-1 !!\n" );
      if ( amr->patch[0][lv][PID]->cornerR[2] - amr->patch[0][lv][PID]->cornerL[2] != PS1-1 )   Aux_Error( ERROR_INFO, "cornerR[2]-cornerL[2] != PS1-1 !!\n" );
#     endif

      // TODO: Do not include the patch has son (not this level) ? If yes, need to update how to fill the array boundary

      HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetExtents( Hypre_grid, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR )   );
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// solve one cell-centered variables
   HYPRE_SStructVariable vartypes[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

   HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetVariables( Hypre_grid, part, NVars, vartypes )   );

// set periodic
   if ( OPT__BC_POT != BC_POT_ISOLATED )
   {
      int periodicity[3] = { NX0_TOT[0]*(1L<<lv), NX0_TOT[1]*(1L<<lv), NX0_TOT[2]*(1L<<lv) };
      HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetPeriodic( Hypre_grid, part, periodicity )   );
   }

// assamble grid
   HYPRE_CHECK_FUNC(   HYPRE_SStructGridAssemble( Hypre_grid )   );

// create stencils
   HYPRE_CHECK_FUNC(   HYPRE_SStructStencilCreate( NDim, NEntries, &Hypre_stencil )   );

// set entries
   for (int e=0; e<NEntries; e++)
   {
      HYPRE_CHECK_FUNC(   HYPRE_SStructStencilSetEntry( Hypre_stencil, e, Offsets[e], var )   );
   }

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



//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_SolvePoisson
// Description :  Solve the poisson equation
//
// Parameter   :  SaveSg_Pot :
//                lv         : Target level
//-------------------------------------------------------------------------------------------------------
void Hypre_SolvePoisson( const int SaveSg_Pot, const int lv, const double TimeNew )
{

   const int NExtend = 1;
   const int part    = 0; // single level only, no need to iterate parts
   const int var     = 0; // single variable only, no need to iterate variables

   // Hypre_Free();
   Hypre_PrepareSingleLevel( lv, NExtend );

// set Hypre arrays
   Hypre_FillArrays( lv, NExtend, TimeNew );

// setup solver and solve
   Hypre_Solve();

// collect the potential
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGather( Hypre_x )   );

// update GAMER array
   real *pote = new real [CUBE(PS1)];
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      // NOTE: By FLASH: Use GetBoxValues more efficient then GetValues.
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorGetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = k*SQR(PS1) + j*PS1 + i;
         amr->patch[ SaveSg_Pot ][lv][PID]->pot[k][j][i] = pote[idx];
         // TODO: apply the floor value here
      } // i, j, k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   delete [] pote;

// clean b array
   real *dens = new real [CUBE(PS1)];
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         // TODO: check the data structure
         const int idx = i*SQR(PS1) + j*PS1 + k;
         dens[idx] = 0.0;
      }

      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, dens )   );

   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   delete [] dens;

} // FUNCTION : Hypre_SolvePoisson



void Hypre_FillArrays( const int lv, const int NExtend, const double TimeNew )
{

   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;
   const int  NDim              = 3;
   const int  NEntries          = 2*NDim + 1;
   const int  var               = 0;
   const int  part              = 0;
   const real dh2               = SQR( amr->dh[lv] );
   const real coeff             = -4.0 * M_PI * NEWTON_G * dh2;

   int stencil_indices[NEntries];

   for (int i=0; i<NEntries; i++)  stencil_indices[i] = i;

// determin if the patch group is at boundary or course-fine edge
   int *PID0_BC_List = new int [ amr->NPatchComma[lv][1] / 8 ];
   int  NPG_BC       = 0;

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

            if ( amr->patch[0][lv][PID]->sibling[s] >= 0 )   continue;
            if ( amr->patch[0][lv][PID]->sibling[s] < -100  &&  OPT__BC_POT == BC_POT_PERIODIC )   continue;

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
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1] / 8; PID0+=8)

   real   *Matrix_Laplace = new real [NEntries*CUBE(PS1)];
   real   *pote           = new real [CUBE(PS1)];
   real   *dens           = new real [CUBE(PS1)];
   real   *bcBox          = new real [1*SQR(PS1)];
   real   *bcVal          = new real [SQR(PS1)];
   real  (*Pot_Array)[CUBE(PS1+2)] = NULL;

// fill common Laplace matrix
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
      const int idx = (k*SQR(PS1) + j*PS1 + i)*NEntries;

      Matrix_Laplace[idx+0] =  6; // i,   j,   k
      Matrix_Laplace[idx+1] = -1; // i-1, j,   k
      Matrix_Laplace[idx+2] = -1; // i+1, j,   k
      Matrix_Laplace[idx+3] = -1; // i,   j-1, k
      Matrix_Laplace[idx+4] = -1; // i,   j+1, k
      Matrix_Laplace[idx+5] = -1; // i,   j,   k-1
      Matrix_Laplace[idx+6] = -1; // i,   j,   k+1
   }

// prepare boundary potential
   if ( NPG_BC > 0 )
   {
      Pot_Array = new real [8*NPG_BC][ CUBE(PS1+2) ];

      Prepare_PatchData( lv, TimeNew, Pot_Array[0], NULL,
                         NExtend, NPG_BC, PID0_BC_List, _POTE, _NONE, OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_06,
                         IntPhase_No, OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
   }

// fill patches
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const int idx = k*SQR(PS1) + j*PS1 + i;

         dens[idx] = coeff * amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
         pote[idx] = 0.0;
      }

      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, dens )   );
      HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_x, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, pote )   );
      HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, amr->patch[0][lv][PID]->cornerL, amr->patch[0][lv][PID]->cornerR, var, NEntries, stencil_indices, Matrix_Laplace )   );
   } // for (int PID=PID0; PID<PID0+8; PID++)

// fill boundary condition
   for (int PG=0; PG<NPG_BC; PG++)
   {
      const int PID0 = PID0_BC_List[PG];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;
         int cornerL[3] = { amr->patch[0][lv][PID]->cornerL[0], amr->patch[0][lv][PID]->cornerL[1], amr->patch[0][lv][PID]->cornerL[2] };
         int cornerR[3] = { amr->patch[0][lv][PID]->cornerR[0], amr->patch[0][lv][PID]->cornerR[1], amr->patch[0][lv][PID]->cornerR[2] };

//       remove the connection and update A due to Dirichlet boundary condition
         for (int s=0; s<6; s++)
         {
            if ( amr->patch[0][lv][PID]->sibling[s] >= 0 )   continue;
            if ( amr->patch[0][lv][PID]->sibling[s] < -100  &&  OPT__BC_POT == BC_POT_PERIODIC )   continue;

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
               const int idx = k * Nj * Ni + j * Ni + i;
               int idx_i, idx_j, idx_k;
               switch ( d )
               {
                  case 0: idx_i = ( lr ) ? PS1+2*NExtend-1 : 0;  idx_j = j + NExtend;                   idx_k = k + NExtend;                   break;
                  case 1: idx_i = i + NExtend;                   idx_j = ( lr ) ? PS1+2*NExtend-1 : 0;  idx_k = k + NExtend;                   break;
                  case 2: idx_i = i + NExtend;                   idx_j = j + NExtend;                   idx_k = ( lr ) ? PS1+2*NExtend-1 : 0;  break;
               } // switch ( d )
               const int idx_arr = idx_k*SQR(PS1+2*NExtend) + idx_j*(PS1+2*NExtend) + idx_i;

               bcBox[idx]  = 0.0;
               bcVal[idx] += Pot_Array[8*PG+LocalID][idx_arr];
            } // i, j, k

            HYPRE_CHECK_FUNC(   HYPRE_SStructVectorSetBoxValues( Hypre_b, part, cornerL_bc, cornerR_bc, var, bcVal )   );
            HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixSetBoxValues( Hypre_A, part, cornerL_bc, cornerR_bc, var, 1, stencil_indices_bc, bcBox )   );
         } // for (int s=0; s<6; s++)
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int PG=0; PG<NPG_BC; PG++)

   delete [] Matrix_Laplace;
   delete [] pote;
   delete [] dens;
   delete [] bcBox;
   delete [] bcVal;
   delete [] PID0_BC_List;
   delete [] Pot_Array;

   HYPRE_CHECK_FUNC(   HYPRE_SStructMatrixAssemble( Hypre_A )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_x )   );
   HYPRE_CHECK_FUNC(   HYPRE_SStructVectorAssemble( Hypre_b )   );

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
