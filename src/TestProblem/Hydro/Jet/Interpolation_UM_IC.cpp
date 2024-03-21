#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include"Typedef.h"

/*----------------------------------------------------------------------------*/
/*! \fn void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size)
 *  *  *  \brief Construct 3D array = array[nt][nr][nc]  */
void ***calloc_3d_array (size_t nt, size_t nr, size_t nc, size_t size)
{
  void ***array;
  size_t i, j;

  if ((array = (void ***) calloc (nt, sizeof (void **))) == NULL)
    {
      printf ("[calloc_3d] failed to allocate memory for %d 1st-pointers\n", (int) nt);
      return NULL;
    }

  if ((array[0] = (void **) calloc (nt * nr, sizeof (void *))) == NULL)
    {
      printf ("[calloc_3d] failed to allocate memory for %d 2nd-pointers\n", (int) (nt * nr));
      free ((void *) array);
      return NULL;
    }

  for (i = 1; i < nt; i++)
    {
      array[i] = (void **) ((unsigned char *) array[0] + i * nr * sizeof (void *));
    }

  if ((array[0][0] = (void *) calloc (nt * nr * nc, size)) == NULL)
    {
      printf ("[calloc_3d] failed to alloc. memory (%d X %d X %d of size %d)\n",
      (int) nt, (int) nr, (int) nc, (int) size);
      free ((void *) array[0]);
      free ((void *) array);
      return NULL;
    }
  for (j = 1; j < nr; j++)
    {
      array[0][j] = (void **) ((unsigned char *) array[0][j - 1] + nc * size);
    }
for (i = 1; i < nt; i++)
    {
      array[i][0] = (void **) ((unsigned char *) array[i - 1][0] + nr * nc * size);
      for (j = 1; j < nr; j++)
        {
          array[i][j] = (void **) ((unsigned char *) array[i][j - 1] + nc * size);
        }
    }

  return array;
}

void free_3d_array(void ***array){

  free(array[0][0]);
  free(array[0]);
  free(array);
}

real TrilinearInterpolation(real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz)
{
  real x1, y1, z1, x0, y0, z0, xd, yd, zd, x, y, z;
  real c000, c001, c010, c100, c011, c101, c110, c111, c00, c01, c10, c11, c0, c1, c;

  x0 = xyz000[0];
  y0 = xyz000[1];
  z0 = xyz000[2];

  x1 = xyz000[0] + dxyz[0];
  y1 = xyz000[1] + dxyz[1];
  z1 = xyz000[2] + dxyz[2];

  x = xyz[0];
  y = xyz[1];
  z = xyz[2];

  c000 = FieldAtVertices[0];
  c001 = FieldAtVertices[1];
  c010 = FieldAtVertices[2];
  c100 = FieldAtVertices[3];
  c011 = FieldAtVertices[4];
  c101 = FieldAtVertices[5];
  c110 = FieldAtVertices[6];
  c111 = FieldAtVertices[7];

  xd = (x-x0)/(x1-x0);
  yd = (y-y0)/(y1-y0);
  zd = (z-z0)/(z1-z0);

  c00 = c000*(1.0-xd) + c100*xd;
  c01 = c001*(1.0-xd) + c101*xd;
  c10 = c010*(1.0-xd) + c110*xd;
  c11 = c011*(1.0-xd) + c111*xd;

  c0  = c00*(1.0-yd) + c10*yd;
  c1  = c01*(1.0-yd) + c11*yd;

  c = c0*(1.0-zd) + c1*zd;

  return c;
}
