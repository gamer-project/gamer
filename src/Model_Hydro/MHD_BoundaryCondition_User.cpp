#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void BC_BField_User( real magnetic[], const double x, const double y, const double z, const double Time,
                            const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*BC_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                            const int lv, double AuxArray[] ) = BC_BField_User;

static void BC_User_xm( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                        const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const double Time, const double dh, const double *Corner, const int lv );
static void BC_User_xp( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                        const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const double Time, const double dh, const double *Corner, const int lv );
static void BC_User_ym( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                        const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const double Time, const double dh, const double *Corner, const int lv );
static void BC_User_yp( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                        const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const double Time, const double dh, const double *Corner, const int lv );
static void BC_User_zm( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                        const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const double Time, const double dh, const double *Corner, const int lv );
static void BC_User_zp( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                        const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const double Time, const double dh, const double *Corner, const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  BC_BField_User
// Description :  User-specified boundary condition for the magnetic field
//
// Note        :  1. Invoked by MHD_BoundaryCondition_User() using the function pointer "BC_BField_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Always return NCOMP_MAG magentic field components
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void BC_BField_User( real magnetic[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] )
{

// put your B.C. here
// ##########################################################################################################
   /*
   magnetic[MAGX] = 1.0;
   magnetic[MAGY] = 2.0;
   magnetic[MAGZ] = 3.0;
   */
// ##########################################################################################################

} // FUNCTION : BC_BField_User



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_BoundaryCondition_User
// Description :  Fill up the ghost-zone **face-centered magnetic field** by the user-specified B.C.
//
// Note        :  1. Work for Prepare_PatchData(), InterpolateGhostZone(), Refine(), and LB_Refine_GetNewRealPatchList()
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                                 --> Must be a pointer array with dimension NCOMP_MAG (3):
//                                 --> Array[MAGX]: Bx array with dimension [ArraySizeZ  ][ArraySizeY  ][ArraySizeX+1]
//                                     Array[MAGY]: By array with dimension [ArraySizeZ  ][ArraySizeY+1][ArraySizeX  ]
//                                     Array[MAGZ]: Bz array with dimension [ArraySizeZ+1][ArraySizeY  ][ArraySizeX  ]
//                                 --> Array[MAG?] can be NULL if MAG? is not specified in TVarIdxList[]
//                BC_Face        : Boundary face (0~5) --> (-x,+x,-y,+y,-z,+z)
//                NVar           : Number of magnetic fields to be prepared
//                ArraySizeX/Y/Z : Size of the corresponding cell-centered array including the ghost zones on each side
//                                 --> See the description of "Array" above
//                Idx_Start      : Minimum array indices (referred to the corresponding cell-centered array)
//                Idx_End        : Maximum array indices (referred to the corresponding cell-centered array)
//                TVarIdxList    : List recording the target magnetic field indices ( = [0 ... NCOMP_MAG-1] )
//                Time           : Current physical time
//                dh             : Cell size
//                Corner         : Cell-centered physcial coordinates of the cell (0,0,0)
//                lv             : Refinement level
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void MHD_BoundaryCondition_User( real **Array, const int BC_Face, const int NVar,
                                 const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                 const int Idx_Start[], const int Idx_End[], const int TVarIdxList[],
                                 const double Time, const double dh, const double *Corner, const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( NVar != 0  &&  TVarIdxList == NULL )
      Aux_Error( ERROR_INFO, "NVar = %d != 0, TVarIdxList == NULL !!\n", NVar );

   for (int v=0; v<NVar; v++)
      if ( Array[ TVarIdxList[v] ] == NULL )    Aux_Error( ERROR_INFO, "Array[%d] is NULL !!\n", TVarIdxList[v] );
#  endif


// set the boundary values at different boundary faces
   switch ( BC_Face )
   {
      case 0:  BC_User_xm( Array, NVar, TVarIdxList, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Time, dh, Corner, lv );  break;
      case 1:  BC_User_xp( Array, NVar, TVarIdxList, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Time, dh, Corner, lv );  break;
      case 2:  BC_User_ym( Array, NVar, TVarIdxList, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Time, dh, Corner, lv );  break;
      case 3:  BC_User_yp( Array, NVar, TVarIdxList, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Time, dh, Corner, lv );  break;
      case 4:  BC_User_zm( Array, NVar, TVarIdxList, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Time, dh, Corner, lv );  break;
      case 5:  BC_User_zp( Array, NVar, TVarIdxList, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, Time, dh, Corner, lv );  break;
      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }

} // FUNCTION : MHD_BoundaryCondition_User



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_xm
// Description :  Set the user-specified B.C. at the -x boundary
//
// Note        :  Work for MHD_BoundaryCondition_User()
//
// Parameter   :  See MHD_BoundaryCondition_User()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_User_xm( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const double Time, const double dh, const double *Corner, const int lv )
{

// cell-centered coordinates at Idx_Start
   const double xyz0_CC[3] = { Corner[0] + (double)Idx_Start[0]*dh,
                               Corner[1] + (double)Idx_Start[1]*dh,
                               Corner[2] + (double)Idx_Start[2]*dh };
   double xyz0_FC[3], x, y, z;
   real   BMag[NCOMP_MAG];
   int    i, j, k;

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            xyz0_FC[0] = xyz0_CC[0] - 0.5*dh;
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagX[k][j][i] = BMag[MAGX];
            }

            break;
         }

         case MAGY:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1] - 0.5*dh;
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]+1; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagY[k][j][i] = BMag[MAGY];
            }

            break;
         }

         case MAGZ:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2] - 0.5*dh;

            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]+1; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagZ[k][j][i] = BMag[MAGZ];
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_User_xm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_xp
// Description :  Set the user-specified B.C. at the +x boundary
//
// Note        :  Work for MHD_BoundaryCondition_User()
//
// Parameter   :  See MHD_BoundaryCondition_User()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_User_xp( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const double Time, const double dh, const double *Corner, const int lv )
{

// cell-centered coordinates at Idx_Start
   const double xyz0_CC[3] = { Corner[0] + (double)Idx_Start[0]*dh,
                               Corner[1] + (double)Idx_Start[1]*dh,
                               Corner[2] + (double)Idx_Start[2]*dh };
   double xyz0_FC[3], x, y, z;
   real   BMag[NCOMP_MAG];
   int    i, j, k;

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            xyz0_FC[0] = xyz0_CC[0] + 0.5*dh;
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (k=Idx_Start[2],   z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1],   y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0]+1, x=xyz0_FC[0]; i<=Idx_End[0]+1; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagX[k][j][i] = BMag[MAGX];
            }

            break;
         }

         case MAGY:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1] - 0.5*dh;
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]+1; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagY[k][j][i] = BMag[MAGY];
            }

            break;
         }

         case MAGZ:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2] - 0.5*dh;

            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]+1; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagZ[k][j][i] = BMag[MAGZ];
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_User_xp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_ym
// Description :  Set the user-specified B.C. at the -y boundary
//
// Note        :  Work for MHD_BoundaryCondition_User()
//
// Parameter   :  See MHD_BoundaryCondition_User()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_User_ym( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const double Time, const double dh, const double *Corner, const int lv )
{

// cell-centered coordinates at Idx_Start
   const double xyz0_CC[3] = { Corner[0] + (double)Idx_Start[0]*dh,
                               Corner[1] + (double)Idx_Start[1]*dh,
                               Corner[2] + (double)Idx_Start[2]*dh };
   double xyz0_FC[3], x, y, z;
   real   BMag[NCOMP_MAG];
   int    i, j, k;

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            xyz0_FC[0] = xyz0_CC[0] - 0.5*dh;
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]+1; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagX[k][j][i] = BMag[MAGX];
            }

            break;
         }

         case MAGY:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1] - 0.5*dh;
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagY[k][j][i] = BMag[MAGY];
            }

            break;
         }

         case MAGZ:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2] - 0.5*dh;

            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]+1; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagZ[k][j][i] = BMag[MAGZ];
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_User_ym



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_yp
// Description :  Set the user-specified B.C. at the +y boundary
//
// Note        :  Work for MHD_BoundaryCondition_User()
//
// Parameter   :  See MHD_BoundaryCondition_User()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_User_yp( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const double Time, const double dh, const double *Corner, const int lv )
{

// cell-centered coordinates at Idx_Start
   const double xyz0_CC[3] = { Corner[0] + (double)Idx_Start[0]*dh,
                               Corner[1] + (double)Idx_Start[1]*dh,
                               Corner[2] + (double)Idx_Start[2]*dh };
   double xyz0_FC[3], x, y, z;
   real   BMag[NCOMP_MAG];
   int    i, j, k;

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            xyz0_FC[0] = xyz0_CC[0] - 0.5*dh;
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]+1; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagX[k][j][i] = BMag[MAGX];
            }

            break;
         }

         case MAGY:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1] + 0.5*dh;
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (k=Idx_Start[2],   z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1]+1, y=xyz0_FC[1]; j<=Idx_End[1]+1; j++, y+=dh)
            for (i=Idx_Start[0],   x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagY[k][j][i] = BMag[MAGY];
            }

            break;
         }

         case MAGZ:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2] - 0.5*dh;

            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]+1; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagZ[k][j][i] = BMag[MAGZ];
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_User_yp



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_zm
// Description :  Set the user-specified B.C. at the -z boundary
//
// Note        :  Work for MHD_BoundaryCondition_User()
//
// Parameter   :  See MHD_BoundaryCondition_User()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_User_zm( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const double Time, const double dh, const double *Corner, const int lv )
{

// cell-centered coordinates at Idx_Start
   const double xyz0_CC[3] = { Corner[0] + (double)Idx_Start[0]*dh,
                               Corner[1] + (double)Idx_Start[1]*dh,
                               Corner[2] + (double)Idx_Start[2]*dh };
   double xyz0_FC[3], x, y, z;
   real   BMag[NCOMP_MAG];
   int    i, j, k;

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            xyz0_FC[0] = xyz0_CC[0] - 0.5*dh;
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]+1; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagX[k][j][i] = BMag[MAGX];
            }

            break;
         }

         case MAGY:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1] - 0.5*dh;
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]+1; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagY[k][j][i] = BMag[MAGY];
            }

            break;
         }

         case MAGZ:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2] - 0.5*dh;

            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2]; k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagZ[k][j][i] = BMag[MAGZ];
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_User_zm



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_zp
// Description :  Set the user-specified B.C. at the +z boundary
//
// Note        :  Work for MHD_BoundaryCondition_User()
//
// Parameter   :  See MHD_BoundaryCondition_User()
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void BC_User_zp( real **Array, const int NVar, const int TVarIdxList[], const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const double Time, const double dh, const double *Corner, const int lv )
{

// cell-centered coordinates at Idx_Start
   const double xyz0_CC[3] = { Corner[0] + (double)Idx_Start[0]*dh,
                               Corner[1] + (double)Idx_Start[1]*dh,
                               Corner[2] + (double)Idx_Start[2]*dh };
   double xyz0_FC[3], x, y, z;
   real   BMag[NCOMP_MAG];
   int    i, j, k;

   for (int v=0; v<NVar; v++)
   {
      switch ( TVarIdxList[v] )
      {
         case MAGX:
         {
            xyz0_FC[0] = xyz0_CC[0] - 0.5*dh;
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagX)[ArraySizeY][ArraySizeX+1] = ( real (*)[ArraySizeY][ArraySizeX+1] )Array[MAGX];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0]+1; i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagX[k][j][i] = BMag[MAGX];
            }

            break;
         }

         case MAGY:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1] - 0.5*dh;
            xyz0_FC[2] = xyz0_CC[2];

            real (*MagY)[ArraySizeY+1][ArraySizeX] = ( real (*)[ArraySizeY+1][ArraySizeX] )Array[MAGY];

            for (k=Idx_Start[2], z=xyz0_FC[2]; k<=Idx_End[2];   k++, z+=dh)
            for (j=Idx_Start[1], y=xyz0_FC[1]; j<=Idx_End[1]+1; j++, y+=dh)
            for (i=Idx_Start[0], x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagY[k][j][i] = BMag[MAGY];
            }

            break;
         }

         case MAGZ:
         {
            xyz0_FC[0] = xyz0_CC[0];
            xyz0_FC[1] = xyz0_CC[1];
            xyz0_FC[2] = xyz0_CC[2] + 0.5*dh;

            real (*MagZ)[ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeY][ArraySizeX] )Array[MAGZ];

            for (k=Idx_Start[2]+1, z=xyz0_FC[2]; k<=Idx_End[2]+1; k++, z+=dh)
            for (j=Idx_Start[1],   y=xyz0_FC[1]; j<=Idx_End[1];   j++, y+=dh)
            for (i=Idx_Start[0],   x=xyz0_FC[0]; i<=Idx_End[0];   i++, x+=dh)
            {
               BC_BField_User_Ptr( BMag, x, y, z, Time, lv, NULL );
               MagZ[k][j][i] = BMag[MAGZ];
            }

            break;
         }

         default:
            Aux_Error( ERROR_INFO, "incorrect B field index (%d) !!\n", TVarIdxList[v] );
      } // switch ( TVarIdxList[v] )
   } // for (int v=0; v<NVar; v++)

} // FUNCTION : BC_User_zp



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
