#include "GAMER.h"

#ifdef SUPPORT_LIBYT




#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  MagX/Y/Z_DerivedFunc
// Description :  Derived function for CCMagX/Y/Z
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
    yt_getGridInfo_FieldData( gid, "CCMagX", &DataRaw );

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
    yt_getGridInfo_FieldData( gid, "CCMagY", &DataRaw );

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
    yt_getGridInfo_FieldData( gid, "CCMagZ", &DataRaw );

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
#endif // #ifdef MHD



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  Temperature_DerivedFunc
// Description :  Derived function for fluid temperature invoke by EoS routine
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. yt_getGridInfo_Dimensions() gets the grid_dimensions[0][1][2] in [x][y][z] coordinate.
//
// Parameter   :  GID        : Grid GID
//                TempData   : Store the derived field data here.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Temperature_DerivedFunc(long gid, double *TempData){
    // Get dim of the grid to be return, and all the other NCOMP_TOTAL fields.
    int Dimensions[3];
    yt_data   DataRaw[NCOMP_TOTAL];
    real     *Data[NCOMP_TOTAL];
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    for (int v=0; v<NCOMP_TOTAL; v++){
        yt_getGridInfo_FieldData( gid, FieldLabel[v], &(DataRaw[v]));
        Data[v] = (real *) DataRaw[v].data_ptr;
    }

    // Preparation for getting Passive and Emag
    real Passive[NCOMP_PASSIVE];
    real Emag;
#ifdef MHD
    yt_data   MagDataRaw[NCOMP_MAG];
    real     *MagData[NCOMP_MAG];
    char     *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};
    for (int v=0; v<NCOMP_MAG; v++){
        yt_getGridInfo_FieldData( gid, CCMagLabel[v], &(MagDataRaw[v]));
        MagData[v] = (real *) MagDataRaw[v].data_ptr;
    }
#endif

    // Get temperature cell-by-cell
    for (int k=0; k<Dimensions[2]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[0]; i++){

                int idx = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];

                // Get Passive
                for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++){
                    Passive[v-NCOMP_FLUID] = Data[v][idx];
                }

                // Get Emag
                #ifdef MHD
                Emag = MHD_GetCellCenteredBEnergy(MagData[MAGX], MagData[MAGY], MagData[MAGZ],
                                                  Dimensions[0], Dimensions[1], Dimensions[2], i, j, k);
                #else
                Emag = NULL_REAL;
                #endif

                // Get temperature
                TempData[idx] = (double) Hydro_Con2Temp(Data[DENS][idx], Data[MOMX][idx], Data[MOMY][idx], Data[MOMZ][idx], Data[ENGY][idx],
                                                        Passive, false, NULL_REAL, Emag,
                                                        EoS_DensEint2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
            }
        }
    }

}
#endif // #if ( MODEL == HYDRO )



#endif // #ifdef SUPPORT_LIBYT
