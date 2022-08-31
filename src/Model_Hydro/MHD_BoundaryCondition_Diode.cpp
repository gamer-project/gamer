#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )

//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_BoundaryCondition_Diode
// Description :  Fill up the ghost-zone **face-centered magnetic field** by the diode B.C.
//
// Note        :  1. Work for Prepare_PatchData(), InterpolateGhostZone(), Refine(), and
//                   LB_Refine_GetNewRealPatchList()
//                2. Similar to the outflow (i.e., zero-gradient) B.C. except that quantities
//                   only "flow" outward
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
void MHD_BoundaryCondition_Diode( real **Array, const int BC_Face, const int NVar, const int GhostSize,
                                  const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                  const int Idx_Start[], const int Idx_End[], const int TVarIdxList[] )
{

// Simply call MHD_BoundaryCondition_Outflow

   MHD_BoundaryCondition_Outflow( Array, BC_Face, NVar, GhostSize, ArraySizeX, ArraySizeY,
                                  ArraySizeZ, Idx_Start, Idx_End, TVarIdxList );

} // FUNCTION : MHD_BoundaryCondition_Diode

#endif // #if ( MODEL == HYDRO  &&  defined MHD )
