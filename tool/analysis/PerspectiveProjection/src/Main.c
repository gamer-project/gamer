#include <sys/time.h>

#include "../include/General.h"
#include "hdf5.h"

#define GCC_VERSION ( __GNUC__ * 10000 \
                      + __GNUC_MINOR__ * 100 \
                      + __GNUC_PATCHLEVEL__ )

#ifdef XRAY_ROSAT
# warning: XRAY_ROSAT defined!
#endif

#ifdef XRAY_EROSITA
# warning: XRAY_EROSITA defined!
#endif

#ifdef HADRONIC_GAMMARAY
# warning: HADRONIC_GAMMARAY defined!
#endif

#ifdef LEPTONIC_GAMMARAY
# warning: LEPTONIC_GAMMARAY defined!
#endif

#ifdef SYNCHROTRON
# warning: SYNCHROTRON defined!
#endif


herr_t LoadCompoundData( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target, const hid_t H5_TypeID_Target );
herr_t checkH5Status( const herr_t status1, const herr_t status2 );

int MPI_Rank, NRank, RootRank;



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char **argv )
{
// 1. initialize
   MPI_Init( &argc, &argv );
   MPI_Comm_rank( MPI_COMM_WORLD, &MPI_Rank );
   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
   RootRank = 0;

   MASTER_PRINT( "NRank = %d, NTHREAD = %d\n", NRank, NUM_THREADS );

// there is run-time error on tw3 with compiler/gcc/7.5.0
#  if ( GCC_VERSION == 70500 )
#     error : We recommand use gcc 9.4.0
#  endif

// 2. parameters
#  if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )
   FILE *table_08_keV = fopen( "cf_0.8keV.dat", "r" );
   FILE *table_15_keV = fopen( "cf_1.5keV.dat", "r" );

   if ( MPI_Rank == RootRank )
   {
      if ( table_15_keV == NULL )   ERROR_EXIT( -1, "ERROR : Could not open file table_15_keV !!\n" );
      if ( table_08_keV == NULL )   ERROR_EXIT( -1, "ERROR : Could not open file table_08_keV !!\n" );
   } // if ( MPI_Rank == RootRank )
#  endif // #if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )


#  ifdef PERSPECTIVE_PROJECTION
// numAngle should be a factor of 360
   const int numAzimuthalAngle = 1;
   const int numCellB = 360, numCellL = 360;

   if ( numCellB != numCellL )   ERROR_EXIT( 0, "ERROR : numCellB != numCellL !!\n" );

   const real b_max = +90.0, b_min = -90.0;
   const real l_max = +180.0, l_min = -180.0;
#  endif // #ifdef PERSPECTIVE_PROJECTION

#  ifdef SLICE
   const char CuttingPlane[1] = "x";
#  endif

// 3. Load data
   hid_t  FRB_file;
   hid_t  FRB_group1, FRB_group2;
   hid_t  FRB_dset1G1, FRB_dset2G1, FRB_dset3G1, FRB_dset4G1, FRB_dset5G1;
   hid_t  H5_SetID_KeyInfo, H5_TypeID_KeyInfo;
   herr_t H5_Status;

// 3-1. Open FRB_file using the default properties.
   FRB_file          = H5Fopen( argv[1], H5F_ACC_RDONLY, H5P_DEFAULT );
   FRB_group2        = H5Gopen( FRB_file,   "Info", H5P_DEFAULT );
   H5_SetID_KeyInfo  = H5Dopen( FRB_group2, "Keys", H5P_DEFAULT );
   H5_TypeID_KeyInfo = H5Dget_type( H5_SetID_KeyInfo );

   int numCellX, numCellY, numCellZ;
   real BoxSizeX, BoxSizeY, BoxSizeZ;
   real Unit_L, Unit_M, Unit_T, Unit_V, Unit_D, Unit_E, Unit_P;

   LoadCompoundData( "numCellX", &numCellX, H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "numCellY", &numCellY, H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "numCellZ", &numCellZ, H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "BoxSizeX", &BoxSizeX, H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "BoxSizeY", &BoxSizeY, H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "BoxSizeZ", &BoxSizeZ, H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_L",   &Unit_L,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_M",   &Unit_M,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_T",   &Unit_T,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_V",   &Unit_V,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_D",   &Unit_D,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_E",   &Unit_E,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
   LoadCompoundData( "Unit_P",   &Unit_P,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo );

   MASTER_PRINT( "numCellX=%d\nnumCellY=%d\nnumCellZ=%d\n", numCellX, numCellY, numCellZ );
   MASTER_PRINT( "BoxSizeX=%e\nBoxSizeY=%e\nBoxSizeZ=%e\n", BoxSizeX, BoxSizeY, BoxSizeZ );
   MASTER_PRINT( "Unit_L=%e\nUnit_M=%e\nUnit_T=%e\nUnit_V=%e\nUnit_D=%e\nUnit_E=%e\nUnit_P=%e\n",
                 Unit_L, Unit_M, Unit_T, Unit_V, Unit_D, Unit_E, Unit_P );

   if ( numCellX%NRank != 0 )   ERROR_EXIT( 0, "numCellX(%d) %% NRank(%d) != 0 !!\n", numCellX, NRank );

// 3D array to store the field to be projected
   real*** Density          = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** Temperature      = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** Pressure         = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** Xray_08_keV      = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** Xray_15_keV      = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** CREngy           = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** Passive_0001     = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** HadronicGammaRay = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** LeptonicGammaRay = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );
   real*** Synchrotron      = (real***)calloc_3d_array( (size_t)numCellX, (size_t)numCellY, (size_t)numCellZ, sizeof(real) );

   real deltaX   = (real)BoxSizeX / (real)numCellX;
   real deltaY   = (real)BoxSizeY / (real)numCellY;
   real deltaZ   = (real)BoxSizeZ / (real)numCellZ;


   real *X = (real*)calloc( (size_t)numCellX, sizeof(real) );
   real *Y = (real*)calloc( (size_t)numCellY, sizeof(real) );
   real *Z = (real*)calloc( (size_t)numCellZ, sizeof(real) );

#  ifdef PERSPECTIVE_PROJECTION
   real dxyz[3]      = { deltaX, deltaY, deltaZ };
   real *XYZ[3]      = { X, Y, Z };
   int numCellXYZ[3] = { numCellX, numCellY, numCellZ };
#  endif

   for ( int x=0; x<numCellX; x++ )   X[x] = ((real)x+0.5) * deltaX - 0.5 * (real)BoxSizeX;
   for ( int y=0; y<numCellY; y++ )   Y[y] = ((real)y+0.5) * deltaY - 0.5 * (real)BoxSizeY;
   for ( int z=0; z<numCellZ; z++ )   Z[z] = ((real)z+0.5) * deltaZ - 0.5 * (real)BoxSizeZ;


// Open group
   FRB_group1  = H5Gopen( FRB_file, "Fields", H5P_DEFAULT );

// Open dataset
   FRB_dset1G1 = H5Dopen( FRB_group1, "Dens",          H5P_DEFAULT );
   FRB_dset2G1 = H5Dopen( FRB_group1, "Pres",          H5P_DEFAULT );
   FRB_dset3G1 = H5Dopen( FRB_group1, "Temp",          H5P_DEFAULT );
   FRB_dset4G1 = H5Dopen( FRB_group1, "CRay",          H5P_DEFAULT );
   FRB_dset5G1 = H5Dopen( FRB_group1, "Passive_0001",  H5P_DEFAULT );

// Read the data using the default properties.
   MASTER_PRINT( "Loading Dens ...\n" );
   H5_Status = H5Dread( FRB_dset1G1, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Density[0][0][0] );
   MASTER_PRINT( "Loading Dens ... done\n" );

   MASTER_PRINT( "Loading Pres ...\n" );
   H5_Status = checkH5Status( H5Dread( FRB_dset2G1, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Pressure[0][0][0] ),     H5_Status );
   MASTER_PRINT( "Loading Pres ... done\n" );

   MASTER_PRINT( "Loading Temp ... \n" );
   H5_Status = checkH5Status( H5Dread( FRB_dset3G1, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Temperature[0][0][0] ),  H5_Status );
   MASTER_PRINT( "Loading Temp ... done\n" );

   MASTER_PRINT( "Loading CRay ... \n" );
   H5_Status = checkH5Status( H5Dread( FRB_dset4G1, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &CREngy[0][0][0] ),       H5_Status );
   MASTER_PRINT( "Loading CRay ... done\n" );

   MASTER_PRINT( "Loading Passive_0001 ...\n" );
   H5_Status = checkH5Status( H5Dread( FRB_dset5G1, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Passive_0001[0][0][0] ), H5_Status );
   MASTER_PRINT( "Loading Passive_0001 ... done\n" );

   if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Dread() fails !!\n" );

// 3-2. Load the table
#  if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )
   int numRow_08_keV = CountNumLine( table_08_keV );
   int numRow_15_keV = CountNumLine( table_15_keV );

   real *tempTable_08_keV    = (real*)calloc( numRow_08_keV, sizeof(real) );
   real *tempTable_15_keV    = (real*)calloc( numRow_15_keV, sizeof(real) );

   real *lambdaTable_08_keV  = (real*)calloc( numRow_08_keV, sizeof(real) );
   real *lambdaTable_15_keV  = (real*)calloc( numRow_15_keV, sizeof(real) );

   ReadTable( 1, numRow_08_keV, table_08_keV, tempTable_08_keV   );
   ReadTable( 1, numRow_15_keV, table_15_keV, tempTable_15_keV   );

   ReadTable( 2, numRow_08_keV, table_08_keV, lambdaTable_08_keV );
   ReadTable( 2, numRow_15_keV, table_15_keV, lambdaTable_15_keV );
#  endif // #if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )

// 4. Convert unit
// Dens     --> Dens * UNIT_D
// Pres     --> Pres * UNIT_P
// ===========================================================
   for (int i=0; i<numCellXYZ[0]; i++)
   {
      for (int j=0; j<numCellXYZ[1]; j++)
      {
#        pragma omp parallel for num_threads(NUM_THREADS)
         for (int k=0; k<numCellXYZ[2]; k++)
         {
           // Temperature[i][j][k] *=  MASS_PROTON_GRAM * MU_ELECTRON * SQR(SPEED_OF_LIGHT) / BOLTZMANN_CONST_ERG;
           Temperature[i][j][k] *=  MU_ELECTRON; // TODO: I have not idea why I need this
           Density    [i][j][k] *=  Unit_D; // g/cm**3
           Pressure   [i][j][k] *=  Unit_P; // erg/cm**3
           CREngy     [i][j][k] *=  Unit_P; // erg/cm**3
           CREngy     [i][j][k] *=  ERG2EV; // eV/cm**3
         } // for (int k=0; k<numCellXYZ[2]; k++)
      } // for (int j=0; j<numCellXYZ[1]; j++)
   } // for (int i=0; i<numCellXYZ[0]; i++)

// 5. Make synthetic emissivity in 3D space
// ===========================================================
   real BoxSize[3] = { BoxSizeX, BoxSizeY, BoxSizeZ };

   MPI_Barrier( MPI_COMM_WORLD );

#  ifdef XRAY_EROSITA
   if ( MPI_Rank == RootRank )
   {
     MASTER_PRINT( "Computing X-ray emissivity at 0.8 keV in 3D space ...\n" );
     X_ray_3D( Density, Temperature, Xray_08_keV, numCellXYZ, numRow_08_keV, tempTable_08_keV, lambdaTable_08_keV );
     MASTER_PRINT( "Computing X-ray emissivity at 0.8 keV in 3D space ... done\n" );
   } // if ( MPI_Rank == RootRank )
#  endif

   MPI_Barrier( MPI_COMM_WORLD );

#  ifdef XRAY_ROSAT
   if ( MPI_Rank == RootRank )
   {
     MASTER_PRINT( "Computing X-ray emissivity at 1.5 keV in 3D space ...\n" );
     X_ray_3D( Density, Temperature, Xray_15_keV, numCellXYZ, numRow_15_keV, tempTable_15_keV, lambdaTable_15_keV );
     MASTER_PRINT( "Computing X-ray emissivity at 1.5 keV in 3D space ... done\n" );
   } // if ( MPI_Rank == RootRank )
#  endif

   MPI_Barrier( MPI_COMM_WORLD );

#  if   ( defined HADRONIC_GAMMARAY  ||  defined LEPTONIC_GAMMARAY )
   real scatteredEnergy = atof(argv[2]);  //  (eV)  100 GeV
#  elif ( defined SYNCHROTRON )
   real observedFreq    = atof(argv[2]);  //  GHz
#  endif

#  if ( defined HADRONIC_GAMMARAY  ||  defined LEPTONIC_GAMMARAY  ||  defined SYNCHROTRON )
   real gamma_max       = atof(argv[3]);
   real gamma_min       = 1.0;
   real spectral_index  = atof(argv[4]);
#  endif

#  if   ( defined HADRONIC_GAMMARAY  ||  defined LEPTONIC_GAMMARAY )
   MASTER_PRINT( "gamma_max=%e, gamma_min=%e, spectral_index=%e, scatteredEnergy=%e GHz\n",
                  gamma_max,    gamma_min,    spectral_index,    scatteredEnergy );
#  elif ( defined SYNCHROTRON )
   MASTER_PRINT( "gamma_max=%e, gamma_min=%e, spectral_index=%e, observedFreq=%e GHz\n",
                  gamma_max,    gamma_min,    spectral_index,    observedFreq );
#  endif


#  ifdef HADRONIC_GAMMARAY
   MASTER_PRINT( "Computing hadronic gamma-ray emissivity in 3D space ...\n" );
   Hadronic_3D( scatteredEnergy, Density, CREngy, gamma_max, gamma_min, spectral_index, HadronicGammaRay, numCellXYZ, BoxSize );
   MASTER_PRINT( "Computing hadronic gamma-ray emissivity in 3D space ... done\n" );
#  endif

   MPI_Barrier( MPI_COMM_WORLD );

#  ifdef LEPTONIC_GAMMARAY
   struct timeval begin, end;

// Start measuring time
   if ( MPI_Rank == RootRank )   gettimeofday( &begin, 0 );

   MASTER_PRINT( "Computing leptonic gamma-ray emissivity in 3D space ...\n" );
   Leptonic_3D( scatteredEnergy, CREngy, gamma_max, gamma_min, spectral_index, LeptonicGammaRay, numCellXYZ, BoxSize, argv );
   MASTER_PRINT( "Computing leptonic gamma-ray emissivity in 3D space ... done\n" );

// Stop measuring time and calculate the elapsed time
   if ( MPI_Rank == RootRank )   gettimeofday( &end, 0 );

   MASTER_PRINT( "Leptonic_3D() takes %.3f seconds.\n", end.tv_sec - begin.tv_sec + 1.e-6*(end.tv_usec - begin.tv_usec) );
#  endif // #ifdef LEPTONIC_GAMMARAY

   MPI_Barrier( MPI_COMM_WORLD );

#  ifdef SYNCHROTRON
   if ( MPI_Rank == RootRank )
   {
      MASTER_PRINT( "Computing synchrotron emissivity in 3D space ...\n" );
      Synchrotron_3D( CREngy, spectral_index, Synchrotron, observedFreq, numCellXYZ, BoxSize, gamma_max, gamma_min );
      MASTER_PRINT( "Computing synchrotron emissivity in 3D space ... done\n" );
   } // if ( MPI_Rank == RootRank )
#  endif

   MPI_Barrier( MPI_COMM_WORLD );

// Make a map for b-l plane
// ===========================================================
#  ifdef PERSPECTIVE_PROJECTION
   real delta_b = (b_max-b_min) / (real)numCellB;
   real delta_l = (l_max-l_min) / (real)numCellL;
   real dt      = 0.0;

   dt  = MAX( dt, BoxSizeX/(real)numCellX );
   dt  = MAX( dt, BoxSizeY/(real)numCellY );
   dt  = MAX( dt, BoxSizeZ/(real)numCellZ );
   dt *= 0.5;

   MASTER_PRINT( "Perform perspective projection ...\n" );
#  endif

   int UpperBound_ll, UpperBound_bb;

#  ifdef PERSPECTIVE_PROJECTION
   UpperBound_ll = numCellL;
   UpperBound_bb = numCellB;
#  endif

#  ifdef SLICE
   int numCellXYZ[3] = { numCellX, numCellY, numCellZ };

   if ( atoi(CuttingPlane) == atoi("x") )
   {
      UpperBound_ll = numCellY;
      UpperBound_bb = numCellZ;
   }
   else if ( atoi(CuttingPlane) == atoi("y") )
   {
      UpperBound_ll = numCellX;
      UpperBound_bb = numCellZ;
   }
   else if ( atoi(CuttingPlane) == atoi("z") )
   {
      UpperBound_ll = numCellX;
      UpperBound_bb = numCellY;
   } // if ( atoi(CuttingPlane)== atoi("x") ) ... else if ...
#  endif // #ifdef SLICE


   real ***ProjectedSliceDensity          = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSlicePressure         = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceTemperature      = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceXray_08_keV      = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceXray_15_keV      = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceCRay             = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSlicePassive_0001     = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceHadronicGammaRay = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceLeptonicGammaRay = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );
   real ***ProjectedSliceSynchrotron      = (real***) calloc_3d_array( (size_t)numAzimuthalAngle, (size_t)UpperBound_ll, (size_t)UpperBound_bb, sizeof(real) );

   real azimuthalAngle = 0.0;

// Why
// #   pragma omp parallel for collapse(2)  num_threads(NUM_THREADS)
// fails?
   for (int i=0; i<numAzimuthalAngle; i++)
   {
      MASTER_PRINT( "Project along the azimuthal angle: %4.3f ...\n", azimuthalAngle*180.0/M_PI );
      for (int ll=0; ll<UpperBound_ll; ll++)
      {
#        pragma omp parallel for num_threads( NUM_THREADS )
         for (int bb=0; bb<UpperBound_bb; bb++)
         {
#           ifdef PERSPECTIVE_PROJECTION
            real l = ((real)ll+0.5) * delta_l + l_min;
            real b = ((real)bb+0.5) * delta_b + b_min;
            ProjectedSliceDensity         [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Density         ,            -1,             NULL,              NULL  );
            ProjectedSlicePressure        [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Pressure        ,            -1,             NULL,              NULL  );
            ProjectedSliceTemperature     [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Temperature     ,            -1,             NULL,              NULL  );
#           if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )
            ProjectedSliceXray_08_keV     [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Xray_08_keV     , numRow_08_keV, tempTable_08_keV, lambdaTable_08_keV );
            ProjectedSliceXray_15_keV     [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Xray_15_keV     , numRow_15_keV, tempTable_15_keV, lambdaTable_15_keV );
#           endif
            ProjectedSliceCRay            [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, CREngy          ,            -1,             NULL,              NULL  );
            ProjectedSlicePassive_0001    [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Passive_0001    ,            -1,             NULL,              NULL  );
            ProjectedSliceHadronicGammaRay[i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, HadronicGammaRay,            -1,             NULL,              NULL  );
            ProjectedSliceLeptonicGammaRay[i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, LeptonicGammaRay,            -1,             NULL,              NULL  );
            ProjectedSliceSynchrotron     [i][bb][ll] = PerspectiveProject( b, l, dt, XYZ, numCellXYZ, dxyz, azimuthalAngle, BoxSize, Synchrotron     ,            -1,             NULL,              NULL  );
#           endif // #ifdef PERSPECTIVE_PROJECTION
#           ifdef SLICE
            const int Idx1 = ll, Idx2 = bb;
            ProjectedSliceDensity         [i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, Density,          CuttingPlane );
            ProjectedSlicePressure        [i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, Pressure,         CuttingPlane );
            ProjectedSliceTemperature     [i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, Temperature,      CuttingPlane );
            ProjectedSliceXray_08_keV     [i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, Xray_08_keV,      CuttingPlane );
            ProjectedSliceXray_15_keV     [i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, Xray_15_keV,      CuttingPlane );
            ProjectedSliceHadronicGammaRay[i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, HadronicGammaRay, CuttingPlane );
            ProjectedSliceLeptonicGammaRay[i][Idx1][Idx2] = Slice( Idx1, Idx2, numCellXYZ, LeptonicGammaRay, CuttingPlane );
#           endif
         } // for (int bb=0; bb<UpperBound_bb; bb++)
      } // for (int ll=0; ll<UpperBound_ll; ll++)

      MASTER_PRINT( "Project along the azimuthal angle: %4.3f ... done\n", azimuthalAngle*180.0/M_PI );

      azimuthalAngle += (real)2.0 * M_PI / (real)numAzimuthalAngle;
   } // for (int i=0; i<numAzimuthalAngle; i++)

   MPI_Barrier( MPI_COMM_WORLD );

   free_3d_array( (void***)Density      );
   free_3d_array( (void***)Pressure     );
   free_3d_array( (void***)Temperature  );
   free_3d_array( (void***)Xray_08_keV  );
   free_3d_array( (void***)Xray_15_keV  );
   free_3d_array( (void***)CREngy       );
   free_3d_array( (void***)Passive_0001 );

// Close dataset
   H5_Status = H5Dclose( FRB_dset1G1 );
   H5_Status = checkH5Status( H5Dclose( FRB_dset2G1 ), H5_Status );
   H5_Status = checkH5Status( H5Dclose( FRB_dset3G1 ), H5_Status );
   H5_Status = checkH5Status( H5Dclose( FRB_dset4G1 ), H5_Status );
   H5_Status = checkH5Status( H5Dclose( FRB_dset5G1 ), H5_Status );

   if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Dclose() fails !!\n" );

// Close group
   H5_Status = H5Gclose( FRB_group1 );
   H5_Status = checkH5Status( H5Gclose( FRB_group2 ), H5_Status );

   if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Gclose() fails !!\n" );

// Close file
   H5_Status = H5Fclose( FRB_file );
   if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : H5Fclose() fails !!\n" );


// 6. Storing perspective projection (pp) map on disk
   if ( MPI_Rank == RootRank )
   {
      hid_t pp_file;
      hid_t pp_group1, pp_group2;
      hid_t pp_DS1G1, pp_DS2G1, pp_DS3G1, pp_DS4G1, pp_DS5G1, pp_DS6G1, pp_DS7G1, pp_DS8G1, pp_DS9G1, pp_DS10G1;
      hid_t pp_dsetKey, typeID;
      hid_t pp_spaceDS1G1, pp_spaceDS2G1, pp_spaceDS3G1, pp_spaceDS4G1, pp_spaceDS5G1, pp_spaceDS6G1, pp_spaceDS7G1, pp_spaceDS8G1, pp_spaceDS9G1, pp_spaceDS10G1;
      hid_t pp_spaceKey;

      Keys_t *Keys = (Keys_t *)calloc( 1, sizeof(Keys_t) );

      strcpy( Keys->FRBDataName, argv[1] );

      Keys->UpperBound_ll     = UpperBound_ll;
      Keys->UpperBound_bb     = UpperBound_bb;
      Keys->numAzimuthalAngle = numAzimuthalAngle;
#     if ( defined HADRONIC_GAMMARAY || defined LEPTONIC_GAMMARAY )
      Keys->scatteredEnergy   = scatteredEnergy;
#     elif ( defined SYNCHROTRON )
      Keys->observedFreq      = observedFreq;
#     endif
#     if ( defined HADRONIC_GAMMARAY || defined LEPTONIC_GAMMARAY || defined SYNCHROTRON )
      Keys->gamma_max         = gamma_max;
      Keys->gamma_min         = gamma_min;
      Keys->spectral_index    = spectral_index;
#     endif
#     ifdef PERSPECTIVE_PROJECTION
      Keys->b_max             = b_max;
      Keys->b_min             = b_min;
      Keys->l_max             = l_max;
      Keys->l_min             = l_min;
      Keys->dt                = dt;
#     endif

      LoadCompoundData( "AMRDataName",    &(Keys->AMRDataName),             H5_SetID_KeyInfo, H5_TypeID_KeyInfo );
      LoadCompoundData( "EpochTimeStamp", &(Keys->EpochTimeStampInFRBData), H5_SetID_KeyInfo, H5_TypeID_KeyInfo );

//    Create a new file using the default properties.
#     ifdef PERSPECTIVE_PROJECTION
      char MapName[STRING_LENGTH] = "Projected_FRB_";
#     endif
#     ifdef SLICE
      char MapName[STRING_LENGTH] = "Slice_FRB_";
#     endif

      strcat( MapName, Keys->AMRDataName );

#     if ( defined SYNCHROTRON )
      strcat( MapName, "_Synchrotron" );
#     endif

#     if ( defined LEPTONIC_GAMMARAY  ||  defined HADRONIC_GAMMARAY  ||  SYNCHROTRON )
      strcat( MapName, "_" );
      strcat( MapName, argv[2] );
      strcat( MapName, "_" );
      strcat( MapName, argv[3] );
      strcat( MapName, "_" );
      strcat( MapName, argv[4] );
#     endif

#     if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )
      strcat( MapName, "_XRay" );
#     endif

      pp_file   = H5Fcreate( strcat(MapName, ".h5"), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

//    Create groups in the file.
      pp_group1 = H5Gcreate( pp_file, "/Map",        H5P_DEFAULT,   H5P_DEFAULT, H5P_DEFAULT );
      pp_group2 = H5Gcreate( pp_file, "/Info",       H5P_DEFAULT,   H5P_DEFAULT, H5P_DEFAULT );

      const hsize_t dims[3] = { numAzimuthalAngle, UpperBound_ll, UpperBound_bb };

//    Create dataspace.  Setting maximum size to NULL sets the maximum size to be the current size.
      pp_spaceDS1G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS2G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS3G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS4G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS5G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS6G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS7G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS8G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS9G1  = H5Screate_simple( 3, dims, NULL );
      pp_spaceDS10G1 = H5Screate_simple( 3, dims, NULL );
      pp_spaceKey    = H5Screate( H5S_SCALAR );

      typeID         = H5Tcreate( H5T_COMPOUND, sizeof(Keys_t) );

      H5Tinsert( typeID, "UpperBound_bb",           HOFFSET(Keys_t,           UpperBound_bb ),    H5T_NATIVE_INT );
      H5Tinsert( typeID, "UpperBound_ll",           HOFFSET(Keys_t,           UpperBound_ll ),    H5T_NATIVE_INT );
      H5Tinsert( typeID, "numAzimuthalAngle",       HOFFSET(Keys_t,        numAzimuthalAngle),    H5T_NATIVE_INT );
#     if ( defined HADRONIC_GAMMARAY || defined LEPTONIC_GAMMARAY )
      H5Tinsert( typeID, "scatteredEnergy",         HOFFSET(Keys_t,          scatteredEnergy),    H5T_IEEE_F32LE );
#     elif ( defined SYNCHROTRON )
      H5Tinsert( typeID, "observedFreq",            HOFFSET(Keys_t,             observedFreq),    H5T_IEEE_F32LE );
#     endif
#     if ( defined HADRONIC_GAMMARAY || defined LEPTONIC_GAMMARAY || defined SYNCHROTRON )
      H5Tinsert( typeID, "gamma_max",               HOFFSET(Keys_t,                gamma_max),    H5T_IEEE_F32LE );
      H5Tinsert( typeID, "gamma_min",               HOFFSET(Keys_t,                gamma_min),    H5T_IEEE_F32LE );
      H5Tinsert( typeID, "spectral_index",          HOFFSET(Keys_t,           spectral_index),    H5T_IEEE_F32LE );
#     endif
#     ifdef PERSPECTIVE_PROJECTION
      H5Tinsert( typeID, "b_max",                   HOFFSET(Keys_t,                   b_max ),    H5T_IEEE_F32LE );
      H5Tinsert( typeID, "b_min",                   HOFFSET(Keys_t,                   b_min ),    H5T_IEEE_F32LE );
      H5Tinsert( typeID, "l_max",                   HOFFSET(Keys_t,                   l_max ),    H5T_IEEE_F32LE );
      H5Tinsert( typeID, "l_min",                   HOFFSET(Keys_t,                   l_min ),    H5T_IEEE_F32LE );
      H5Tinsert( typeID, "dt",                      HOFFSET(Keys_t,                      dt ),    H5T_IEEE_F32LE );
#     endif
      H5Tinsert( typeID, "EpochTimeStampInFRBData", HOFFSET(Keys_t, EpochTimeStampInFRBData ), H5T_NATIVE_UINT64 );

//    create the "variable-length string" datatype
      hid_t  H5_TypeID_VarStr;

      H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
      H5_Status        = H5Tset_size( H5_TypeID_VarStr, STRING_LENGTH );
      if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : H5Tset_size() fails: %d !!\n", __LINE__ );

      H5Tinsert( typeID, "AMRDataName",             HOFFSET(Keys_t,             AMRDataName ),  H5_TypeID_VarStr );
      H5Tinsert( typeID, "FRBDataName",             HOFFSET(Keys_t,             FRBDataName ),  H5_TypeID_VarStr );
#     ifdef SLICE
      H5Tinsert( typeID, "CuttingPlane",            HOFFSET(Keys_t,             CuttingPlane ), H5_TypeID_VarStr );
#     endif
      H5Tclose( H5_TypeID_VarStr );


//    Create the dataset.  We will use all default properties for this example.
#     ifdef PERSPECTIVE_PROJECTION
      pp_DS1G1   = H5Dcreate( pp_group1, "ProjectedDensity" ,         H5T_IEEE_F32LE, pp_spaceDS1G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS2G1   = H5Dcreate( pp_group1, "ProjectedPressure",         H5T_IEEE_F32LE, pp_spaceDS2G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS3G1   = H5Dcreate( pp_group1, "ProjectedTemperature",      H5T_IEEE_F32LE, pp_spaceDS3G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS4G1   = H5Dcreate( pp_group1, "ProjectedXray_08_keV",      H5T_IEEE_F32LE, pp_spaceDS4G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS5G1   = H5Dcreate( pp_group1, "ProjectedXray_15_keV",      H5T_IEEE_F32LE, pp_spaceDS5G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS6G1   = H5Dcreate( pp_group1, "ProjectedCRay",             H5T_IEEE_F32LE, pp_spaceDS6G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS7G1   = H5Dcreate( pp_group1, "ProjectedPassive_0001",     H5T_IEEE_F32LE, pp_spaceDS7G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS8G1   = H5Dcreate( pp_group1, "ProjectedHadronicGammaRay", H5T_IEEE_F32LE, pp_spaceDS8G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS9G1   = H5Dcreate( pp_group1, "ProjectedLeptonicGammaRay", H5T_IEEE_F32LE, pp_spaceDS9G1,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS10G1  = H5Dcreate( pp_group1, "ProjectedSynchrotron",      H5T_IEEE_F32LE, pp_spaceDS10G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#     endif
#     ifdef SLICE
      pp_DS1G1   = H5Dcreate( pp_group1, "SliceDensity" ,             H5T_IEEE_F32LE, pp_spaceDS1G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS2G1   = H5Dcreate( pp_group1, "SlicePressure",             H5T_IEEE_F32LE, pp_spaceDS2G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS3G1   = H5Dcreate( pp_group1, "SliceTemperature",          H5T_IEEE_F32LE, pp_spaceDS3G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS4G1   = H5Dcreate( pp_group1, "SliceXray_08_keV",          H5T_IEEE_F32LE, pp_spaceDS4G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS5G1   = H5Dcreate( pp_group1, "SliceXray_15_keV",          H5T_IEEE_F32LE, pp_spaceDS5G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      pp_DS6G1   = H5Dcreate( pp_group1, "SliceGammaRay",             H5T_IEEE_F32LE, pp_spaceDS5G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#     endif
      pp_dsetKey = H5Dcreate( pp_group2, "Keys",                      typeID,         pp_spaceKey,   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

//    Write the data to the dataset.
      H5_Status = H5Dwrite( pp_DS1G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS1G1,  H5P_DEFAULT, ProjectedSliceDensity         [0][0] );
      H5_Status = checkH5Status( H5Dwrite( pp_DS2G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS2G1,  H5P_DEFAULT, ProjectedSlicePressure        [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS3G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS3G1,  H5P_DEFAULT, ProjectedSliceTemperature     [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS4G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS4G1,  H5P_DEFAULT, ProjectedSliceXray_08_keV     [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS5G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS5G1,  H5P_DEFAULT, ProjectedSliceXray_15_keV     [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS6G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS6G1,  H5P_DEFAULT, ProjectedSliceCRay            [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS7G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS7G1,  H5P_DEFAULT, ProjectedSlicePassive_0001    [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS8G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS8G1,  H5P_DEFAULT, ProjectedSliceHadronicGammaRay[0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS9G1,   H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS9G1,  H5P_DEFAULT, ProjectedSliceLeptonicGammaRay[0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_DS10G1,  H5T_IEEE_F32LE, H5S_ALL,  pp_spaceDS10G1, H5P_DEFAULT, ProjectedSliceSynchrotron     [0][0] ), H5_Status );
      H5_Status = checkH5Status( H5Dwrite( pp_dsetKey, typeID,         H5S_ALL,  H5S_ALL,        H5P_DEFAULT, Keys                                 ), H5_Status );

      if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Dwrite() fails: %d !!\n", __LINE__ );

//    Close dataset
      H5_Status = H5Dclose( pp_dsetKey );
      H5_Status = checkH5Status( H5Dclose( pp_DS1G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS2G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS3G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS4G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS5G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS6G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS7G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS8G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS9G1  ), H5_Status );
      H5_Status = checkH5Status( H5Dclose( pp_DS10G1 ), H5_Status );

      if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Dclose() fails: %d !!\n", __LINE__ );

//    Close space
      H5_Status = H5Sclose( pp_spaceKey );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS1G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS2G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS3G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS4G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS5G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS6G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS7G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS8G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS9G1  ), H5_Status );
      H5_Status = checkH5Status( H5Sclose( pp_spaceDS10G1 ), H5_Status );

      if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Sclose() fails: %d !!\n", __LINE__ );

//    Close group
      H5_Status = H5Gclose( pp_group1 );
      H5_Status = checkH5Status( H5Gclose( pp_group2 ), H5_Status );
      if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : One of the H5Gclose() fails: %d !!\n", __LINE__ );

//    Close file
      H5_Status = H5Fclose( pp_file );
      if ( H5_Status < 0 )   ERROR_EXIT( 0, "ERROR : H5Fclose() fails: %d !!\n", __LINE__ );

      free(Keys);
   } // if ( MPI_Rank == RootRank )

   MPI_Barrier( MPI_COMM_WORLD );

   free_3d_array( (void***)ProjectedSliceDensity          );
   free_3d_array( (void***)ProjectedSlicePressure         );
   free_3d_array( (void***)ProjectedSliceTemperature      );
   free_3d_array( (void***)ProjectedSliceXray_08_keV      );
   free_3d_array( (void***)ProjectedSliceXray_15_keV      );
   free_3d_array( (void***)ProjectedSliceCRay             );
   free_3d_array( (void***)ProjectedSlicePassive_0001     );
   free_3d_array( (void***)ProjectedSliceHadronicGammaRay );
   free_3d_array( (void***)ProjectedSliceLeptonicGammaRay );
   free_3d_array( (void***)ProjectedSliceSynchrotron      );

#  if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )
   free( tempTable_08_keV   );
   free( tempTable_15_keV   );
   free( lambdaTable_08_keV );
   free( lambdaTable_15_keV );
   fclose( table_08_keV );
   fclose( table_15_keV );
#  endif

   free( X );
   free( Y );
   free( Z );

   MASTER_PRINT( "Done!!\n" );

   MPI_Finalize();

   return 0;
} // FUNCTION : main



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadCompoundData
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
herr_t LoadCompoundData( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target, const hid_t H5_TypeID_Target )
{
   int    H5_FieldIdx;
   size_t H5_FieldSize;
   hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
   hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
   herr_t H5_Status;

// load
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldName );

   if ( H5_FieldIdx < 0 )   return -1;

   H5_TypeID_Field = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
   H5_FieldSize    = H5Tget_size( H5_TypeID_Field );
   H5_TypeID_Load  = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
   H5_Status       = H5Tinsert( H5_TypeID_Load, FieldName, 0, H5_TypeID_Field );

   H5_Status       = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, FieldPtr );
   if ( H5_Status < 0 )    printf( "failed to load the field!!\n" );

   H5_Status       = H5Tclose( H5_TypeID_Field );
   H5_Status       = H5Tclose( H5_TypeID_Load  );

   return 0;
} // FUNCTION : LoadCompoundData



//-------------------------------------------------------------------------------------------------------
// Function    :  checkH5Status
// Description :  check if any of the H5_Status is fail
// Parameter   :  status1/2 : HDF5 function returned status
//
// Return      :  status1
//-------------------------------------------------------------------------------------------------------
herr_t checkH5Status( const herr_t status1, const herr_t status2 )
{
   if ( status1 < 0 )   return status1;
   if ( status2 < 0 )   return status2;
   return status1;
} // FUNCTION : checkH5Status
