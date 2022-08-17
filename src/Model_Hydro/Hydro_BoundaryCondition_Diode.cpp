#include "GAMER.h"

#if ( MODEL == HYDRO )

static void BC_Diode_xm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                         const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                         const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Diode_xp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                         const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                         const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Diode_ym( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                         const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                         const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Diode_yp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                         const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                         const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Diode_zm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                         const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                         const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Diode_zp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                         const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                         const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_BoundaryCondition_Diode
// Description :  Fill up the ghost-zone values by the Diode B.C.
//
// Note        :  1. Work for Prepare_PatchData(), InterpolateGhostZone(), Refine(), and LB_Refine_GetNewRealPatchList()
//                2. Similar to the outflow (i.e., zero-gradient) B.C. except that quantities only "flow" outward
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
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
void Hydro_BoundaryCondition_Diode( real *Array, const int BC_Face, const int NVar_Flu, const int GhostSize,
                                    const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                    const int Idx_Start[], const int Idx_End[], const int TFluVarIdxList[],
                                    const int NVar_Der, const long TDerVarList[] )
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
      case 0:  BC_Diode_xm( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                            ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 1:  BC_Diode_xp( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                            ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 2:  BC_Diode_ym( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                            ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 3:  BC_Diode_yp( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                            ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 4:  BC_Diode_zm( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                            ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      case 5:  BC_Diode_zp( Array, NVar_Flu, TFluVarIdxList, NVar_Der, TDerVarList, GhostSize,
                            ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );   break;

      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }

} // FUNCTION : Hydro_BoundaryCondition_Diode



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Diode_xm
// Description :  Set the diode B.C. at the -x boundary
//
// Note        :  1. Work for Hydro_BoundaryCondition_Diode()
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  See Hydro_BoundaryCondition_Diode()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Diode_xm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                  const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                  const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = GhostSize;    // reference i index
   int TFluVarIdx;

// 1D array -> 3D array
   typedef real (*vla)[ArraySizeZ][ArraySizeY][ArraySizeX];
   vla Array3D = ( vla )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMX )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MIN(Array3D[v][k][j][i_ref], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j][i_ref];

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
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MIN(Array3D[v][k][j][i_ref], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j][i_ref];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Diode_xm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Diode_xp
// Description :  Set the diode B.C. at the +x boundary
//
// Note        :  1. Work for Hydro_BoundaryCondition_Diode()
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  See Hydro_BoundaryCondition_Diode()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Diode_xp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                  const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                  const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = ArraySizeX - GhostSize - 1;    // reference i index
   int TFluVarIdx;

// 1D array -> 3D array
   typedef real (*vla)[ArraySizeZ][ArraySizeY][ArraySizeX];
   vla Array3D = ( vla )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMX )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MAX(Array3D[v][k][j][i_ref], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j][i_ref];

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
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MAX(Array3D[v][k][j][i_ref], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j][i_ref];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Diode_xp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Diode_ym
// Description :  Set the diode B.C. at the -y boundary
//
// Note        :  1. Work for Hydro_BoundaryCondition_Diode()
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  See Hydro_BoundaryCondition_Diode()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Diode_ym( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                  const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                  const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = GhostSize;  // reference j index
   int TFluVarIdx;

// 1D array -> 3D array
   typedef real (*vla)[ArraySizeZ][ArraySizeY][ArraySizeX];
   vla Array3D = ( vla )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MIN(Array3D[v][k][j_ref][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j_ref][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MIN(Array3D[v][k][j_ref][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j_ref][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Diode_ym



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Diode_yp
// Description :  Set the diode B.C. at the +y boundary
//
// Note        :  1. Work for Hydro_BoundaryCondition_Diode()
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  See Hydro_BoundaryCondition_Diode()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Diode_yp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                  const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                  const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = ArraySizeY - GhostSize - 1;   // reference j index
   int TFluVarIdx;

// 1D array -> 3D array
   typedef real (*vla)[ArraySizeZ][ArraySizeY][ArraySizeX];
   vla Array3D = ( vla )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MAX(Array3D[v][k][j_ref][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j_ref][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELY )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MAX(Array3D[v][k][j_ref][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k][j_ref][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Diode_yp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Diode_zm
// Description :  Set the diode B.C. at the -z boundary
//
// Note        :  1. Work for Hydro_BoundaryCondition_Diode()
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  See Hydro_BoundaryCondition_Diode()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Diode_zm( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                  const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                  const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = GhostSize;  // reference k index
   int TFluVarIdx;

// 1D array -> 3D array
   typedef real (*vla)[ArraySizeZ][ArraySizeY][ArraySizeX];
   vla Array3D = ( vla )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MIN(Array3D[v][k_ref][j][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k_ref][j][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MIN(Array3D[v][k_ref][j][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k_ref][j][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Diode_zm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Diode_zp
// Description :  Set the diode B.C. at the +z boundary
//
// Note        :  1. Work for Hydro_BoundaryCondition_Diode()
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target fluid variables
//
// Parameter   :  See Hydro_BoundaryCondition_Diode()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Diode_zp( real *Array, const int NVar_Flu, const int TFluVarIdxList[], const int NVar_Der,
                  const long TDerVarList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                  const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = ArraySizeZ - GhostSize - 1;   // reference k index
   int TFluVarIdx;

// 1D array -> 3D array
   typedef real (*vla)[ArraySizeZ][ArraySizeY][ArraySizeX];
   vla Array3D = ( vla )Array;


// set the boundary values
   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      if ( TFluVarIdx == MOMZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MAX(Array3D[v][k_ref][j][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k_ref][j][i];

         }}}
      }
   } // for (int v=0; v<NVar_Flu; v++)


// derived variables
   for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)
   {
      if ( TDerVarList[v2] == _VELZ )
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = MAX(Array3D[v][k_ref][j][i], 0.0);

         }}}
      }

      else
      {
         for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)    {
         for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)    {
         for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)    {

            Array3D[v][k][j][i] = Array3D[v][k_ref][j][i];

         }}}
      }
   } // for (int v2=0, v=NVar_Flu; v2<NVar_Der; v2++, v++)

} // FUNCTION : BC_Diode_zp



#endif // if ( MODEL == HYDRO )
