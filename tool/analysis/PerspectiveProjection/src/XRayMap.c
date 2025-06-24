#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  X_ray_3D
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void X_ray_3D( real ***Density, real ***Temperature, real ***Emissivity,
               int numCellXYZ[], int numRow, real *tempTable, real *lambdaTable )
{
   for (int i=0; i<numCellXYZ[0]; i++)
   {
      for (int j=0; j<numCellXYZ[1]; j++)
      {
#        pragma omp parallel for num_threads( NUM_THREADS )
         for (int k=0; k<numCellXYZ[2]; k++)
         {
            real lambda         = Lambda( Temperature[i][j][k], numRow, tempTable, lambdaTable );
            Emissivity[i][j][k] = Xray_emissivity( Density[i][j][k], lambda ); // * ERG2EV;
         } // for (int k=0; k<numCellXYZ[2]; k++)
      } // for (int j=0; j<numCellXYZ[1]; j++)
   } // for (int i=0; i<numCellXYZ[0]; i++)
} // FUNCTION : X_ray_3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Xray_emissivity
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real Xray_emissivity( real dens, real lambda )
{
   return SQR( dens / (MU_ELECTRON * MASS_PROTON_GRAM) ) * lambda;
} // FUNCTION : Xray_emissivity
