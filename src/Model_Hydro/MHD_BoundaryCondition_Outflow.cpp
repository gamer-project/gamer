#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )


static void BC_Outflow_xm( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_xp( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_ym( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_yp( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_zm( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );
static void BC_Outflow_zp( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                           const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_BoundaryCondition_Outflow
// Description :  Fill up the ghost-zone **face-centered magnetic field** by the outflow B.C.
//
// Note        :  1. Work for Prepare_PatchData(), InterpolateGhostZone(), Refine(), and LB_Refine_GetNewRealPatchList()
//                2. Specifically, the so-called outflow B.C. is actually the **zero-gradient** B.C.
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                                 --> Must be a pointer array with dimension NCOMP_MAG (3):
//                                 --> Array[MAGX]: Bx array with dimension [ArraySizeZ  ][ArraySizeY  ][ArraySizeX+1]
//                                     Array[MAGY]: By array with dimension [ArraySizeZ  ][ArraySizeY+1][ArraySizeX  ]
//                                     Array[MAGZ]: Bz array with dimension [ArraySizeZ+1][ArraySizeY  ][ArraySizeX  ]
//                                 --> Array[MAG?] can be NULL if MAG? is not specified in TVarIdxList[]
//                BC_Face        : Boundary face (0~5) --> (-x,+x,-y,+y,-z,+z)
//                NVar           : Number of magnetic fields to be prepared
//                GhostSize      : Number of ghost zones
//                ArraySizeX/Y/Z : Size of the corresponding cell-centered array including the ghost zones on each side
//                                 --> See the description of "Array" above
//                Idx_Start      : Minimum array indices (referred to the corresponding cell-centered array)
//                Idx_End        : Maximum array indices (referred to the corresponding cell-centered array)
//                TVarIdxList    : List recording the target magnetic field indices ( = [0 ... NCOMP_MAG-1] )
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void MHD_BoundaryCondition_Outflow( real **Array, const int BC_Face, const int NVar, const int GhostSize,
                                    const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                    const int Idx_Start[], const int Idx_End[], const int TVarIdxList[] )
{

// check
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

   if ( NVar != 0  &&  TVarIdxList == NULL )
      Aux_Error( ERROR_INFO, "NVar = %d != 0, TVarIdxList == NULL !!\n", NVar );

   for (int v=0; v<NVar; v++)
      if ( Array[ TVarIdxList[v] ] == NULL )    Aux_Error( ERROR_INFO, "Array[%d] is NULL !!\n", TVarIdxList[v] );
#  endif // #ifdef GAMER_DEBUG


// set the boundary values at different boundary faces
   switch ( BC_Face )
   {
      case 0:  BC_Outflow_xm( Array, NVar, TVarIdxList, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break;
      case 1:  BC_Outflow_xp( Array, NVar, TVarIdxList, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break;
      case 2:  BC_Outflow_ym( Array, NVar, TVarIdxList, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break;
      case 3:  BC_Outflow_yp( Array, NVar, TVarIdxList, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break;
      case 4:  BC_Outflow_zm( Array, NVar, TVarIdxList, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break;
      case 5:  BC_Outflow_zp( Array, NVar, TVarIdxList, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End );  break;
      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }

} // FUNCTION : MHD_BoundaryCondition_Outflow



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_xm
// Description :  Set the outflow B.C. at the -x boundary
//
// Note        :  Work for MHD_BoundaryCondition_Outflow()
//
// Parameter   :  See MHD_BoundaryCondition_Outflow()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_xm( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref = GhostSize;  // reference i index for both the longitudinal and transverse components

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
               MagX[k][j][i] = MagX[k][j][i_ref];

            break;
         }

         case MAGY:
         {
            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]+1; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagY[k][j][i] = MagY[k][j][i_ref];

            break;
         }

         case MAGZ:
         {
            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (int k=Idx_Start[2]; k<=Idx_End[2]+1; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagZ[k][j][i] = MagZ[k][j][i_ref];

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_Outflow_xm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_xp
// Description :  Set the outflow B.C. at the +x boundary
//
// Note        :  Work for MHD_BoundaryCondition_Outflow()
//
// Parameter   :  See MHD_BoundaryCondition_Outflow()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_xp( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int i_ref_n = ArraySizeX - GhostSize;  // reference i index for the longitudinal component
   const int i_ref_t = i_ref_n - 1;             // reference i index for the transverse   component

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (int k=Idx_Start[2];   k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1];   j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]+1; i<=Idx_End[0]+1; i++)
               MagX[k][j][i] = MagX[k][j][i_ref_n];

            break;
         }

         case MAGY:
         {
            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]+1; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagY[k][j][i] = MagY[k][j][i_ref_t];

            break;
         }

         case MAGZ:
         {
            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (int k=Idx_Start[2]; k<=Idx_End[2]+1; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagZ[k][j][i] = MagZ[k][j][i_ref_t];

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_Outflow_xp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_ym
// Description :  Set the outflow B.C. at the -y boundary
//
// Note        :  Work for MHD_BoundaryCondition_Outflow()
//
// Parameter   :  See MHD_BoundaryCondition_Outflow()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_ym( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref = GhostSize;  // reference j index for both the longitudinal and transverse components

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]+1; i++)
               MagX[k][j][i] = MagX[k][j_ref][i];

            break;
         }

         case MAGY:
         {
            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
               MagY[k][j][i] = MagY[k][j_ref][i];

            break;
         }

         case MAGZ:
         {
            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (int k=Idx_Start[2]; k<=Idx_End[2]+1; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagZ[k][j][i] = MagZ[k][j_ref][i];

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_Outflow_ym



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_yp
// Description :  Set the outflow B.C. at the +y boundary
//
// Note        :  Work for MHD_BoundaryCondition_Outflow()
//
// Parameter   :  See MHD_BoundaryCondition_Outflow()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_yp( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int j_ref_n = ArraySizeY - GhostSize;  // reference j index for the longitudinal component
   const int j_ref_t = j_ref_n - 1;             // reference j index for the transverse   component

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]+1; i++)
               MagX[k][j][i] = MagX[k][j_ref_t][i];

            break;
         }

         case MAGY:
         {
            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (int k=Idx_Start[2];   k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]+1; j<=Idx_End[1]+1; j++)
            for (int i=Idx_Start[0];   i<=Idx_End[0];   i++)
               MagY[k][j][i] = MagY[k][j_ref_n][i];

            break;
         }

         case MAGZ:
         {
            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (int k=Idx_Start[2]; k<=Idx_End[2]+1; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagZ[k][j][i] = MagZ[k][j_ref_t][i];

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_Outflow_yp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_zm
// Description :  Set the outflow B.C. at the -z boundary
//
// Note        :  Work for MHD_BoundaryCondition_Outflow()
//
// Parameter   :  See MHD_BoundaryCondition_Outflow()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_zm( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref = GhostSize;  // reference k index for both the longitudinal and transverse components

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]+1; i++)
               MagX[k][j][i] = MagX[k_ref][j][i];

            break;
         }

         case MAGY:
         {
            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]+1; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagY[k][j][i] = MagY[k_ref][j][i];

            break;
         }

         case MAGZ:
         {
            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (int k=Idx_Start[2]; k<=Idx_End[2]; k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]; i++)
               MagZ[k][j][i] = MagZ[k_ref][j][i];

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_Outflow_zm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Outflow_zp
// Description :  Set the outflow B.C. at the +z boundary
//
// Note        :  Work for MHD_BoundaryCondition_Outflow()
//
// Parameter   :  See MHD_BoundaryCondition_Outflow()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_Outflow_zp( real **Array, const int NVar, const int TVarIdxList[], const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                    const int ArraySizeZ, const int Idx_Start[], const int Idx_End[] )
{

   const int k_ref_n = ArraySizeZ - GhostSize;  // reference k index for the longitudinal component
   const int k_ref_t = k_ref_n - 1;             // reference k index for the transverse   component

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0]+1; i++)
               MagX[k][j][i] = MagX[k_ref_t][j][i];

            break;
         }

         case MAGY:
         {
            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (int k=Idx_Start[2]; k<=Idx_End[2];   k++)
            for (int j=Idx_Start[1]; j<=Idx_End[1]+1; j++)
            for (int i=Idx_Start[0]; i<=Idx_End[0];   i++)
               MagY[k][j][i] = MagY[k_ref_t][j][i];

            break;
         }

         case MAGZ:
         {
            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (int k=Idx_Start[2]+1; k<=Idx_End[2]+1; k++)
            for (int j=Idx_Start[1];   j<=Idx_End[1];   j++)
            for (int i=Idx_Start[0];   i<=Idx_End[0];   i++)
               MagZ[k][j][i] = MagZ[k_ref_n][j][i];

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_Outflow_zp



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
