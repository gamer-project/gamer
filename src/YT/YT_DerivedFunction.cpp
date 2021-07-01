#include "GAMER.h"


#ifdef SUPPORT_LIBYT

#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  MagX/Y/Z_DerivedFunc
// Description :  Derived function for MagX/Y/Z
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. yt_getGridInfo_Dimensions() gets the grid_dimensions[0][1][2] in [x][y][z] coordinate.
//
// Parameter   :  GID            : Grid GID
//                Converted_MagX : Store the converted field data here.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void MagX_DerivedFunc(long gid, double *Converted_MagX){
    // Get the dimension of the grid, and the data array pointer of the grid
    int Dimensions[3];
    yt_data DataRaw;
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    yt_getGridInfo_FieldData( gid, "MagX", &DataRaw );

    // Cast the DataRaw
    real *Data = (real *) DataRaw.data_ptr;

    // Compute converted field via data. [z, y, x] direction.
    for (int k=0; k<Dimensions[2]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[0]; i++){
                int idx_Bx = IDX321_BX(i, j, k, Dimensions[0], Dimensions[1]);
                int idx_cc = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];
                Converted_MagX[idx_cc] = (double) 0.5 * ( Data[ idx_Bx ] + Data[ idx_Bx + 1 ] );
            }
        }
    }
}

void MagY_DerivedFunc(long gid, double *Converted_MagY){
    // Get the dimension of the grid, and the data array pointer of the grid
    int Dimensions[3];
    yt_data DataRaw;
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    yt_getGridInfo_FieldData( gid, "MagY", &DataRaw );

    // Cast the DataRaw
    real *Data = (real *) DataRaw.data_ptr;

    // Compute converted field via data. [z, y, x] direction.
    for (int k=0; k<Dimensions[2]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[0]; i++){
                int idx_By = IDX321_BY(i, j, k, Dimensions[0], Dimensions[1]);
                int idx_cc = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];
                Converted_MagY[idx_cc] = (double) 0.5 * ( Data[ idx_By ] + Data[ idx_By + Dimensions[0] ] );
            }
        }
    }
}

void MagZ_DerivedFunc(long gid, double *Converted_MagZ){
    // Get the dimension of the grid, and the data array pointer of the grid
    int Dimensions[3];
    yt_data DataRaw;
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    yt_getGridInfo_FieldData( gid, "MagZ", &DataRaw );

    // Cast the DataRaw
    real *Data = (real *) DataRaw.data_ptr;

    // Compute converted field via data. [z, y, x] direction.
    for (int k=0; k<Dimensions[2]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[0]; i++){
                int idx_Bz = IDX321_BZ(i, j, k, Dimensions[0], Dimensions[1]);
                int idx_cc = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];
                Converted_MagZ[idx_cc] = (double) 0.5 * ( Data[ idx_Bz ] + Data[ idx_Bz + Dimensions[0] * Dimensions[1] ] );
            }
        }
    }
}

#endif

#endif