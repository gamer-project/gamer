#include "GAMER.h"

#if ( MODEL == HYDRO )

static void BC_Reflecting_xm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                              const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                              const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Reflecting_xp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                              const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                              const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Reflecting_ym( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                              const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                              const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Reflecting_yp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                              const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                              const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Reflecting_zm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                              const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                              const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Reflecting_zp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                              const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                              const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_BoundaryCondition_Reflecting
// Description :  Fill up the ghost-zone values by the reflecting B.C.
//
// Note        :  1. Work for the functions "Prepare_PatchData, InterpolateGhostZone, Refine, LB_Refine_AllocateNewPatch"
//                2. Similar to the outflow (i.e., zero-gradient) B.C. except that the normal vecotor components change sign
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                BC_Face        : Boundary face (0~5) --> (-x,+x,-y,+y,-z,+z)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void Hydro_BoundaryCondition_Reflecting( real *Array, const int BC_Face, const int NVar_Flu, const int GhostSize,
                                         const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                         const int Idx_Start[], const int Idx_End[], const int TFluVarIdxList[],
                                         const int NVar_Der, const int TDerVarList[] )
{

// check the index range
#  ifdef GAMER_DEBUG
   switch ( BC_Face )
   {
      case 0:  if ( Idx_Start[0] != 0  ||  Idx_End[0] != GhostSize-1 )
                  Aux_Error( ERROR_INFO, "incorrect index range (Start %d, End %d, GhostSize %d, Face %d) !!\n",
                             Idx_Start[0], Idx_End[0], GhostSize, BC_Face );  break;

      case 1:  if ( Idx_Start[0] != ArraySizeX-GhostSize  ||  Idx_End[0] != ArraySizeX-1 )
                  Aux_Error( ERROR_INFO, "incorrect index range (Start %d, End %d, GhostSize %d, Face %d) !!\n",
                             Idx_Start[0], Idx_End[0], GhostSize, BC_Face );  break;

      case 2:  if ( Idx_Start[1] != 0  ||  Idx_End[1] != GhostSize-1 )
                  Aux_Error( ERROR_INFO, "incorrect index range (Start %d, End %d, GhostSize %d, Face %d) !!\n",
                             Idx_Start[1], Idx_End[1], GhostSize, BC_Face );  break;

      case 3:  if ( Idx_Start[1] != ArraySizeY-GhostSize  ||  Idx_End[1] != ArraySizeY-1 )
                  Aux_Error( ERROR_INFO, "incorrect index range (Start %d, End %d, GhostSize %d, Face %d) !!\n",
                             Idx_Start[1], Idx_End[1], GhostSize, BC_Face );  break;

      case 4:  if ( Idx_Start[2] != 0  ||  Idx_End[2] != GhostSize-1 )
                  Aux_Error( ERROR_INFO, "incorrect index range (Start %d, End %d, GhostSize %d, Face %d) !!\n",
                             Idx_Start[2], Idx_End[2], GhostSize, BC_Face );  break;

      case 5:  if ( Idx_Start[2] != ArraySizeZ-GhostSize  ||  Idx_End[2] != ArraySizeZ-1 )
                  Aux_Error( ERROR_INFO, "incorrect index range (Start %d, End %d, GhostSize %d, Face %d) !!\n",
                             Idx_Start[2], Idx_End[2], GhostSize, BC_Face );  break;

      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   } // switch ( BC_Face )

   if ( NVar_Flu != 0  &&  TFluVarIdxList == NULL )
      Aux_Error( ERROR_INFO, "NVar_Flu = %d != 0, TFluVarIdxList == NULL !!\n", NVar_Flu );

   if ( NVar_Der != 0  &&  TDerVarList == NULL )
      Aux_Error( ERROR_INFO, "NVar_Der = %d != 0, TDerVarList == NULL !!\n", NVar_Der );
#  endif // #ifdef GAMER_DEBUG


// set the boundary values at different boundary faces
   switch ( BC_Face )
   {
      case 0:  BC_Reflecting_xm( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                                 ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 1:  BC_Reflecting_xp( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                                 ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 2:  BC_Reflecting_ym( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                                 ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 3:  BC_Reflecting_yp( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                                 ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 4:  BC_Reflecting_zm( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                                 ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 5:  BC_Reflecting_zp( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                                 ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }

} // FUNCTION : Hydro_BoundaryCondition_Reflecting



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Reflecting_xm
// Description :  Set the reflecting B.C. at the -x boundary
//
// Note        :  1. Work for the function "Hydro_BoundaryCondition_Reflecting"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Reflecting_xm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                       const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                       const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = 2*GhostSize-1;    // reference i index
   int TFluVarIdx, ii;

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMX )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = -Array3D[v][k][j][ii];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = Array3D[v][k][j][ii];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELX )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = -Array3D[v][k][j][ii];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = Array3D[v][k][j][ii];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Reflecting_xm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Reflecting_xp
// Description :  Set the reflecting B.C. at the +x boundary
//
// Note        :  1. Work for the function "Hydro_BoundaryCondition_Reflecting"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Reflecting_xp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                       const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                       const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = 2*( ArraySizeX - GhostSize ) - 1;    // reference i index
   int TFluVarIdx, ii;

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMX )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = -Array3D[v][k][j][ii];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = Array3D[v][k][j][ii];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELX )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = -Array3D[v][k][j][ii];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {  ii = i_ref - i;

            Array3D[v][k][j][i] = Array3D[v][k][j][ii];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Reflecting_xp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Reflecting_ym
// Description :  Set the reflecting B.C. at the -y boundary
//
// Note        :  1. Work for the function "Hydro_BoundaryCondition_Reflecting"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Reflecting_ym( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                       const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                       const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = 2*GhostSize-1;    // reference j index
   int TFluVarIdx, jj;

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][k][jj][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][jj][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][k][jj][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][jj][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Reflecting_ym



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Reflecting_yp
// Description :  Set the reflecting B.C. at the +y boundary
//
// Note        :  1. Work for the function "Hydro_BoundaryCondition_Reflecting"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Reflecting_yp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                       const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                       const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = 2*( ArraySizeY - GhostSize ) - 1;  // reference j index
   int TFluVarIdx, jj;

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][k][jj][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][jj][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][k][jj][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {  jj = j_ref - j;
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][jj][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Reflecting_yp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Reflecting_zm
// Description :  Set the reflecting B.C. at the -z boundary
//
// Note        :  1. Work for the function "Hydro_BoundaryCondition_Reflecting"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Reflecting_zm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                       const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                       const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = 2*GhostSize-1;    // reference k index
   int TFluVarIdx, kk;

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][kk][j][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][kk][j][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][kk][j][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][kk][j][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Reflecting_zm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Reflecting_zp
// Description :  Set the reflecting B.C. at the +z boundary
//
// Note        :  1. Work for the function "Hydro_BoundaryCondition_Reflecting"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables is NOT included)
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Reflecting_zp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                       const int TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                       const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = 2*( ArraySizeZ - GhostSize ) - 1;  // reference k index
   int TFluVarIdx, kk;

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][kk][j][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][kk][j][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = -Array3D[v][kk][j][i];

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {  kk = k_ref - k;
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][kk][j][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Reflecting_zp



#endif // if ( MODEL == HYDRO )
