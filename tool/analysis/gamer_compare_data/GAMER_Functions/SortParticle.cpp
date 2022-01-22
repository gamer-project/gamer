#include "CompareData.h"

static void Mis_Heapsort( const long N, real Array[], long IdxTable[] );
static void Heapsort_SiftDown( const long L, const long R, real Array[], long IdxTable[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  SortParticle
// Description :  Sort particles by their x, y, z, and vx
//-------------------------------------------------------------------------------------------------------
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, const real *VelX, long *IdxTable )
{

   long NParSameX, NParSameY, NParSameZ;

// 0. back up the PosX array since we don't want to modify it during the sorting
   real *PosX_Sorted = new real [NPar];

   memcpy( PosX_Sorted, PosX, sizeof(real)*NPar );


// 1. sort by x
   Mis_Heapsort( NPar, PosX_Sorted, IdxTable );


// 2. sort by y
   for (long x=0; x<NPar-1; x++)
   {
      NParSameX = 1;    // 1 --> itself

//    find all particles with the same x coordinate
      while ( x+NParSameX<NPar  &&  PosX_Sorted[x] == PosX_Sorted[x+NParSameX] )    NParSameX ++;

      if ( NParSameX > 1 )
      {
         long *SortByX_IdxTable = new long [NParSameX];
         long *SortByY_IdxTable = new long [NParSameX];
         real *PosY_Sorted      = new real [NParSameX];

         for (long y=0; y<NParSameX; y++)
         {
            SortByX_IdxTable[y] = IdxTable[ x + y ];
            PosY_Sorted     [y] = PosY[ SortByX_IdxTable[y] ];
         }

         Mis_Heapsort( NParSameX, PosY_Sorted, SortByY_IdxTable );

         for (long y=0; y<NParSameX; y++)    IdxTable[ x + y ] = SortByX_IdxTable[ SortByY_IdxTable[y] ];


//       3. sort by z
         for (long y=0; y<NParSameX-1; y++)
         {
            NParSameY = 1;    // 1 --> itself

//          find all particles with the same x and y coordinates
            while ( y+NParSameY<NParSameX  &&  PosY_Sorted[y] == PosY_Sorted[y+NParSameY] )  NParSameY ++;

            if ( NParSameY > 1 )
            {
               long *SortByZ_IdxTable = new long [NParSameY];
               real *PosZ_Sorted      = new real [NParSameY];

               for (long z=0; z<NParSameY; z++)
               {
                  SortByY_IdxTable[z] = IdxTable[ x + y + z ];
                  PosZ_Sorted     [z] = PosZ[ SortByY_IdxTable[z] ];
               }

               Mis_Heapsort( NParSameY, PosZ_Sorted, SortByZ_IdxTable );

               for (long z=0; z<NParSameY; z++)    IdxTable[ x + y + z ] = SortByY_IdxTable[ SortByZ_IdxTable[z] ];


//             4. sort by vx
               for (long z=0; z<NParSameY-1; z++)
               {
                  NParSameZ = 1;    // 1 --> itself

//                find all particles with the same x, y, and z coordinates
                  while ( z+NParSameZ<NParSameY  &&  PosZ_Sorted[z] == PosZ_Sorted[z+NParSameZ] )  NParSameZ ++;

                  if ( NParSameZ > 1 )
                  {
                     long *SortByVx_IdxTable = new long [NParSameZ];
                     real *VelX_Sorted       = new real [NParSameZ];

                     for (long w=0; w<NParSameZ; w++)
                     {
                        SortByZ_IdxTable[w] = IdxTable[ x + y + z + w ];
                        VelX_Sorted     [w] = VelX[ SortByZ_IdxTable[w] ];
                     }

                     Mis_Heapsort( NParSameZ, VelX_Sorted, SortByVx_IdxTable );

                     for (long w=0; w<NParSameZ; w++)    IdxTable[ x + y + z + w ] = SortByZ_IdxTable[ SortByVx_IdxTable[w] ];

                     delete [] SortByVx_IdxTable;
                     delete [] VelX_Sorted;

                     z += NParSameZ - 1;
                  } // if ( NParSameZ > 1 )
               } // for (long z=0; z<NParSameY-1; z++)

               delete [] SortByZ_IdxTable;
               delete [] PosZ_Sorted;

               y += NParSameY - 1;
            } // if ( NParSameY > 1 )
         } // for (long y=0; y<NParSameX-1; y++)

         delete [] SortByX_IdxTable;
         delete [] SortByY_IdxTable;
         delete [] PosY_Sorted;

         x += NParSameX - 1;
      } // if ( NParSameX > 1 )
   } // for (long x=0; x<NPar-1; x++)

   delete [] PosX_Sorted;

} // FUNCTION : SortParticle



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Heapsort
// Description :  Use the Heapsort algorithm to sort the input array into ascending numerical order
//                --> An index table will also be constructed if "IdxTable != NULL"
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//
// Parameter   :  N        :  Size of Array
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
void Mis_Heapsort( const long N, real Array[], long IdxTable[] )
{

// initialize the IdxTable
   if ( IdxTable != NULL )
      for (long t=0; t<N; t++)   IdxTable[t] = t;

// heap creation
   for (long L=N/2-1; L>=0; L--)    Heapsort_SiftDown( L, N-1, Array, IdxTable );

// retirement-and-promotion
   real Buf;
   for (long R=N-1; R>0; R--)
   {
      Buf      = Array[R];
      Array[R] = Array[0];
      Array[0] = Buf;

      if ( IdxTable != NULL )
      {
         Buf         = IdxTable[R];
         IdxTable[R] = IdxTable[0];
         IdxTable[0] = Buf;
      }

      Heapsort_SiftDown( 0, R-1, Array, IdxTable );
   }

} // FUNCTION : Mis_Heapsort



//-------------------------------------------------------------------------------------------------------
// Function    :  Heapsort_SiftDown
// Description :  Sift-down process for the Heapsort algorithm
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//
// Parameter   :  L        :  Left  range of the sift-down
//                R        :  Right range of the sift-down
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
void Heapsort_SiftDown( const long L, const long R, real Array[], long IdxTable[] )
{

   long Idx_up    = L;
   long Idx_down  = 2*Idx_up + 1;
   real Target    = Array[Idx_up];
   long TargetIdx = ( IdxTable == NULL ) ? -1 : IdxTable[Idx_up];

   while ( Idx_down <= R )
   {
//    find the better employee
      if ( Idx_down < R  &&  Array[Idx_down+1] > Array[Idx_down] )   Idx_down ++;

//    terminate the sift-down process if the target (supervisor) is better than both its employees
      if ( Target >= Array[Idx_down] )    break;

//    otherwise, promote the better employee
      Array[Idx_up] = Array[Idx_down];
      if ( IdxTable != NULL )    IdxTable[Idx_up] = IdxTable[Idx_down];

//    prepare the next sift-down operation
      Idx_up   = Idx_down;
      Idx_down = 2*Idx_up + 1;
   }

// put target at its best position
   Array[Idx_up] = Target;
   if ( IdxTable != NULL )    IdxTable[Idx_up] = TargetIdx;

} // FUNCTION : Heapsort_SiftDown
