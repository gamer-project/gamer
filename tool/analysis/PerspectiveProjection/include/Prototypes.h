#ifndef __PROTOTYPES__
#define __PROTOTYPES__

#include "Macro.h"


// global variables
extern int MPI_Rank, NRank, RootRank;



// utilities
void checkInt32Overflow( int32_t a, int32_t b, int operation, int line );
bool checkNAN( real input, const char FunctionName[], const int line );
bool checkMinus ( real input, const char FunctionName[], const int line );
bool checkmemoryContiguous( real *Ptr, int sizeOf, int Length );
void OutputBinary( void *Ptr, int size, int count, char Name [] );

// x_ray_map.c
void X_ray_3D( real ***Density, real ***Temperature, real ***Emissivity,
               int numCellXYZ[], int numRow, real *tempTable, real *lambdaTable );
real Xray_emissivity( real dens, real lambda );

void Synchrotron_3D( real ***CREngy, real spectral_index, real ***Synchrotron, real observedFreq,
                     int numCellXYZ[], real BoxSize[], real gamma_max, real gamma_min );
void Leptonic_3D( real scatteredEnergy, real ***CREngy, real gamma_max, real gamma_min, real spectral_index,
                  real ***Emissivity, int numCellXYZ[], real BoxSize[], char **argv );
void Hadronic_3D( real scatteredEnergy, real ***Density, real ***CREngy,
                  real gamma_max, real gamma_min, real spectral_index,
                  real ***Emissivity, int numCellXYZ[], real BoxSize[] );
real GammaRay_Hadronic_Emissivity( real number_density, real CREngy, real scatteredEnergy,
                                   real gamma_max, real gamma_min, real spectral_index );
void Index2Position( const int Index[], double Pos[], const int numCellXYZ[], const real BoxSize[] );

// make_projection.c
real PerspectiveProject( real b, real l, real dt, real *XYZ[], int numCellXYZ[], real dxyz[], real azimuthalAngle,
                         real BoxSize[], real ***TargetQuantity, int numRow, real *tempTable, real *lambdaTable );
void rotationMatrix( real x, real y, real z, real *xp, real *yp, real *zp, real angle );

// make_slice.c
real Slice( int Idx1, int Idx2, int numCellXYZ[], real ***TargetQuantity, char CuttingPlane[] );

// Interpolation.c
real TrilinearInterpolation( real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz );
real LinearInterpolation( real x1, real y1, real x2, real y2, real x );

// Array.c
void ***calloc_3d_array( size_t nt, size_t nr, size_t nc, size_t size );
void free_3d_array( void ***array );
void **calloc_2d_array( size_t nr, size_t nc, size_t size );
void free_2d_array( void *array );

// BinarySearch.c
int BinarySearch( const real Array[], int Min, int Max, const real Key );
int double_BinarySearch( const double Array[], int Min, int Max, const real Key );

// cooling.c
real Lambda( real Temp, int numRow, real *tempTable, real *lambdaTable );

// read_table.c
int CountNumLine( FILE *table );
void ReadTable( int TargetColumn, int numRow, FILE *table, float *columnArray );

// haloFun.c
real haloDensFun( const real x, const real y, const real z );
real xRayHalo( const real x, const real y, const real z, int numRow, real *tempTable, real *lambdaTable );


#endif // #ifndef __PROTOTYPES__
