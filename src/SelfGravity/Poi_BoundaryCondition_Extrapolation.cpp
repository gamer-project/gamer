#include "GAMER.h"

static void BC_Extrapolation_xm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] );
static void BC_Extrapolation_xp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] );
static void BC_Extrapolation_ym( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] );
static void BC_Extrapolation_yp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] );
static void BC_Extrapolation_zm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] );
static void BC_Extrapolation_zp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_BoundaryCondition_Extrapolation
// Description :  Fill up the ghost-zone potential by extrapolation
//
// Note        :  1. Work for the function "Prepare_PatchData"
//                2. Quadratic polynomial is adopted for extrapolation
//                3. Used only for the isolated boundary condition
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                BC_Face        : Boundary face (0~5) --> (-x,+x,-y,+y,-z,+z)
//                NVar           : Number of variables to be prepared
//                                 --> Currently it's always 1 for potential but it works for NVar>1
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void Poi_BoundaryCondition_Extrapolation( real *Array, const int BC_Face, const int NVar, const int GhostSize,
                                          const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                          const int Idx_Start[], const int Idx_End[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( GhostSize < 0 )
      Aux_Error( ERROR_INFO, "GhostSize (%d) < 0 !!\n", GhostSize );

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


// determine the extrapolation coefficients with the quadratic polynomial
// --> Coeff[0/GhostSize-1][]: closest/furthest ghost    cell to the ghost-interior interface
//     Coeff[][0/2]          : closest/furthest interior cell to the ghost-interior interface
   real (*Coeff)[3] = new real [GhostSize][3];

   for (int t=0; t<GhostSize; t++)
   {
      Coeff[t][0] = +0.5*(t+2.0)*(t+3.0);
      Coeff[t][1] = -1.0*(t+1.0)*(t+3.0);
      Coeff[t][2] = +0.5*(t+1.0)*(t+2.0);
   }


// set the boundary values at different boundary faces
   switch ( BC_Face )
   {
      case 0:  BC_Extrapolation_xm( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Coeff );  break;
      case 1:  BC_Extrapolation_xp( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Coeff );  break;
      case 2:  BC_Extrapolation_ym( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Coeff );  break;
      case 3:  BC_Extrapolation_yp( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Coeff );  break;
      case 4:  BC_Extrapolation_zm( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Coeff );  break;
      case 5:  BC_Extrapolation_zp( Array, NVar, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Coeff );  break;
      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }


// free memory
   delete [] Coeff;

} // FUNCTION : Poi_BoundaryCondition_Extrapolation



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Extrapolation_xm
// Description :  Set the BC at the -x boundary
//
// Note        :  Work for the function "Poi_BoundaryCondition_Extrapolation"
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                Coeff          : Extrapolation coefficients
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Extrapolation_xm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                          const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] )
{

   const int i_ref = GhostSize;  // reference i index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_End[0], ii=0; i>=Idx_Start[0]; i--, ii++)
      Array3D[v][k][j][i] = Coeff[ii][0]*Array3D[v][k][j][i_ref  ] +
                            Coeff[ii][1]*Array3D[v][k][j][i_ref+1] +
                            Coeff[ii][2]*Array3D[v][k][j][i_ref+2];

} // FUNCTION : BC_Extrapolation_xm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Extrapolation_xp
// Description :  Set the BC at the +x boundary
//
// Note        :  Work for the function "Poi_BoundaryCondition_Extrapolation"
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                Coeff          : Extrapolation coefficients
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Extrapolation_xp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                          const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] )
{

   const int i_ref = ArraySizeX - GhostSize - 1;   // reference i index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0], ii=0; i<=Idx_End[0]; i++, ii++)
      Array3D[v][k][j][i] = Coeff[ii][0]*Array3D[v][k][j][i_ref  ] +
                            Coeff[ii][1]*Array3D[v][k][j][i_ref-1] +
                            Coeff[ii][2]*Array3D[v][k][j][i_ref-2];

} // FUNCTION : BC_Extrapolation_xp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Extrapolation_ym
// Description :  Set the BC at the -y boundary
//
// Note        :  Work for the function "Poi_BoundaryCondition_Extrapolation"
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                Coeff          : Extrapolation coefficients
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Extrapolation_ym( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                          const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] )
{

   const int j_ref = GhostSize;  // reference j index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_End[1], jj=0; j>=Idx_Start[1]; j--, jj++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Coeff[jj][0]*Array3D[v][k][j_ref  ][i] +
                            Coeff[jj][1]*Array3D[v][k][j_ref+1][i] +
                            Coeff[jj][2]*Array3D[v][k][j_ref+2][i];

} // FUNCTION : BC_Extrapolation_ym



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Extrapolation_yp
// Description :  Set the BC at the +y boundary
//
// Note        :  Work for the function "Poi_BoundaryCondition_Extrapolation"
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                Coeff          : Extrapolation coefficients
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Extrapolation_yp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                          const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] )
{

   const int j_ref = ArraySizeY - GhostSize - 1;   // reference j index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
   for (int j=Idx_Start[1], jj=0; j<=Idx_End[1]; j++, jj++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Coeff[jj][0]*Array3D[v][k][j_ref  ][i] +
                            Coeff[jj][1]*Array3D[v][k][j_ref-1][i] +
                            Coeff[jj][2]*Array3D[v][k][j_ref-2][i];

} // FUNCTION : BC_Extrapolation_yp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Extrapolation_zm
// Description :  Set the BC at the -z boundary
//
// Note        :  Work for the function "Poi_BoundaryCondition_Extrapolation"
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                Coeff          : Extrapolation coefficients
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Extrapolation_zm( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                          const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] )
{

   const int k_ref = GhostSize;  // reference k index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_End[2], kk=0; k>=Idx_Start[2]; k--, kk++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Coeff[kk][0]*Array3D[v][k_ref  ][j][i] +
                            Coeff[kk][1]*Array3D[v][k_ref+1][j][i] +
                            Coeff[kk][2]*Array3D[v][k_ref+2][j][i];

} // FUNCTION : BC_Extrapolation_zm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Extrapolation_zp
// Description :  Set the BC at the +z boundary
//
// Note        :  Work for the function "Poi_BoundaryCondition_Extrapolation"
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar           : Number of variables to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                Coeff          : Extrapolation coefficients
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Extrapolation_zp( real *Array, const int NVar, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                          const int ArraySizeZ, const int Idx_Start[], const int Idx_End[], const real Coeff[][3] )
{

   const int k_ref = ArraySizeZ - GhostSize - 1;   // reference k index

// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

// set the boundary values
   for (int v=0; v<NVar; v++)
   for (int k=Idx_Start[2], kk=0; k<=Idx_End[2]; k++, kk++)
   for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
   for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
      Array3D[v][k][j][i] = Coeff[kk][0]*Array3D[v][k_ref  ][j][i] +
                            Coeff[kk][1]*Array3D[v][k_ref-1][j][i] +
                            Coeff[kk][2]*Array3D[v][k_ref-2][j][i];

} // FUNCTION : BC_Extrapolation_zp
