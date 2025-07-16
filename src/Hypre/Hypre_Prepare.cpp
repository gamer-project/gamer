#include "GAMER.h"


#ifdef SUPPORT_HYPRE
//-------------------------------------------------------------------------------------------------------
// Function    :  Hypre_PrepareSingleLevel
// Description :  Prepare the single level Hypre arrays
//
// Parameter   :  lv       : Target level
//                NExtend  : Number of cells to extend, mainly for Dirichlet boundary condition (NExtend = 1)
//                Periodic : Whether the grid is periodic or not
//-------------------------------------------------------------------------------------------------------
// TODO: rethink if we need NExtend
void Hypre_PrepareSingleLevel( const int lv, const int NExtend, const bool Periodic[] )
{

   const int NParts      = 1; // number of AMR levels
   const int NVars       = 1; // one var, potential
   const int part        = 0; // one part only, no need to iterate
   const int var         = 0; // one variable only, no need to iterate
   const int NDim        = 3;
   const int NEntries    = 2*NDim + 1;
   const int object_type = HYPRE_SSTRUCT; // or this HYPRE_PARCSR // TODO: check which one should I use
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
   bool SetPeriodic = false;
   int  periodicity[3];
   for (int d=0; d<3; d++)
   {
      periodicity[d]  = ( Periodic[d] ) ? NX0_TOT[d]*(int)(1L<<lv) : 0;
      SetPeriodic    |= Periodic[d];
   }

   if ( SetPeriodic )
      HYPRE_CHECK_FUNC(   HYPRE_SStructGridSetPeriodic( Hypre_grid, part, periodicity )   );

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
#endif // #ifdef SUPPORT_HYPRE
