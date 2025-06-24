#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "Typedef.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  calloc_3d_array
// Description :  Allocate 3D array with size nt*nr*nc
//
// Parameter   :  nt/r/c :  size of array
//                size    : size of allocated type
//
// Return      :  c
//-------------------------------------------------------------------------------------------------------
void ***calloc_3d_array( size_t nt, size_t nr, size_t nc, size_t size )
{

   void ***array;
   size_t i, j;

   if ( (array = (void ***)calloc( nt, sizeof(void **) )) == NULL )
   {
      printf( "[calloc_3d] failed to allocate memory for %d 1st-pointers\n", (int)nt );
      return NULL;
   }

   if ( (array[0] = (void **)calloc( nt*nr, sizeof(void *) )) == NULL )
   {
      printf( "[calloc_3d] failed to allocate memory for %d 2nd-pointers\n", (int)(nt*nr) );
      free( (void *) array );
      return NULL;
   }

   for (i=1; i<nt; i++)
   {
      array[i] = (void **)( (unsigned char *)array[0] + i*nr*sizeof(void *) );
   }

   if ( (array[0][0] = (void *)calloc( nt*nr*nc, size )) == NULL )
   {
      printf( "[calloc_3d] failed to alloc. memory (%d X %d X %d of size %d)\n",
              (int)nt, (int)nr, (int)nc, (int)size );
      free( (void *)array[0] );
      free( (void *)array    );
      return NULL;
   }

   for (j=1; j<nr; j++)
   {
      array[0][j] = (void **)( (unsigned char *)array[0][j-1] + nc*size );
   }

   for (i=1; i<nt; i++)
   {
      array[i][0] = (void **)( (unsigned char *)array[i-1][0] + nr*nc*size );
      for (j=1; j<nr; j++)
      {
         array[i][j] = (void **)( (unsigned char *)array[i][j-1] + nc*size );
      }
   }

   return array;

} // FUNCTION : ***calloc_3d_array



//-------------------------------------------------------------------------------------------------------
// Function    :  free_3d_array
// Description :  Free the allocated 3D array by calloc_3d_array
//
// Parameter   :  array : pointer of the array
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void free_3d_array( void ***array )
{
   free( array[0][0] );
   free( array[0]    );
   free( array       );
} // FUNCTION : free_3d_array



//-------------------------------------------------------------------------------------------------------
// Function    :  TrilinearInterpolation
// Description :  Linear interpolation in 3-D spaces
//
// Parameter   :  FieldAtVertices : the pointer of the field values at near-by coordinates [2*2*2]
//                                  --> Coordinate order: [ (0,0,0), (0,0,1), (0,1,0), (1,0,0), (0,1,1), (1,0,1), (1,1,0), (1,1,1) ]
//                xyz000          : the left coordinates of x/y/z (0,0,0)
//                dxyz            : the width/distance of x/y/z (dx,dy,dz)
//                xyz             : the coordinates of x/y/z need to be interpolation
//
// Return      :  c
//-------------------------------------------------------------------------------------------------------
real TrilinearInterpolation( real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz )
{

   real c = 0.0, weight[8];

// weight of the left / right
   const real w_xR = (xyz[0]-xyz000[0]) / dxyz[0];
   const real w_yR = (xyz[1]-xyz000[1]) / dxyz[1];
   const real w_zR = (xyz[2]-xyz000[2]) / dxyz[2];

   const real w_xL = 1.0 - w_xR;
   const real w_yL = 1.0 - w_yR;
   const real w_zL = 1.0 - w_zR;

// total weight
   weight[0] = w_xL * w_yL * w_zL;
   weight[1] = w_xL * w_yL * w_zR;
   weight[2] = w_xL * w_yR * w_zL;
   weight[3] = w_xR * w_yL * w_zL;
   weight[4] = w_xL * w_yR * w_zR;
   weight[5] = w_xR * w_yL * w_zR;
   weight[6] = w_xR * w_yR * w_zL;
   weight[7] = w_xR * w_yR * w_zR;

   for (int i=0; i<8; i++)   c += FieldAtVertices[i] * weight[i];

   return c;

} // FUNCTION : TrilinearInterpolation
