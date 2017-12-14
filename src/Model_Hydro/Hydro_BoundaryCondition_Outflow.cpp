#include "GAMER.h"

#if ( MODEL == HYDRO )

static void BC_Outflow_xm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_xp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_ym( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_yp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_zm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_zp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_BoundaryCondition_Outflow
// Description :  Fill up the ghost-zone values by the outflow B.C.
//
// Note        :  Work for the functions "Prepare_PatchData, InterpolateGhostZone, Refine, LB_Refine_AllocateNewPatch"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                BC_Face        : Boundary face (0~5) --> (-x,+x,-y,+y,-z,+z)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void Hydro_BoundaryCondition_Outflow( real *Array, const int BC_Face, const int NVar, const int GhostSize, 
                                      const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ, 
                                      const int Idx_Start[], const int Idx_End[] )
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
#  endif // #ifdef GAMER_DEBUG


// set the boundary values at different boundary faces
   switch ( BC_Face )
   {
      case 0:  BC_Outflow_xm( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break; 
      case 1:  BC_Outflow_xp( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break; 
      case 2:  BC_Outflow_ym( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break; 
      case 3:  BC_Outflow_yp( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break; 
      case 4:  BC_Outflow_zm( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break; 
      case 5:  BC_Outflow_zp( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break; 
      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }

} // FUNCTION : Hydro_BoundaryCondition_Outflow



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_xm
// Description :  Set the outflow B.C. at the -x boundary
//
// Note        :  Work for the function "Hydro_BoundaryCondition_Outflow"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_xm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = GhostSize;  // reference i index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Array3D[v][k][j][i_ref];

} // FUNCTION : BC_Outflow_xm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_xp
// Description :  Set the outflow B.C. at the +x boundary
//
// Note        :  Work for the function "Hydro_BoundaryCondition_Outflow"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_xp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = ArraySizeX - GhostSize - 1;    // reference i index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Array3D[v][k][j][i_ref];

} // FUNCTION : BC_Outflow_xp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_ym
// Description :  Set the outflow B.C. at the -y boundary
//
// Note        :  Work for the function "Hydro_BoundaryCondition_Outflow"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_ym( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = GhostSize;  // reference j index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Array3D[v][k][j_ref][i];

} // FUNCTION : BC_Outflow_ym



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_yp
// Description :  Set the outflow B.C. at the +y boundary
//
// Note        :  Work for the function "Hydro_BoundaryCondition_Outflow"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_yp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = ArraySizeY - GhostSize - 1;    // reference j index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Array3D[v][k][j_ref][i];

} // FUNCTION : BC_Outflow_yp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_zm
// Description :  Set the outflow B.C. at the -z boundary
//
// Note        :  Work for the function "Hydro_BoundaryCondition_Outflow"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_zm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = GhostSize;  // reference k index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Array3D[v][k_ref][j][i];

} // FUNCTION : BC_Outflow_zm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_zp
// Description :  Set the outflow B.C. at the +z boundary
//
// Note        :  Work for the function "Hydro_BoundaryCondition_Outflow"
// 
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of fluid and derived variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side 
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_zp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY, 
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = ArraySizeZ - GhostSize - 1;    // reference k index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Array3D[v][k_ref][j][i];

} // FUNCTION : BC_Outflow_zp



#endif // if ( MODEL == HYDRO )
