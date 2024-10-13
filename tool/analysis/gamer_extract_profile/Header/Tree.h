#ifndef __TREE_H__
#define __TREE_H__



// block size (number of cells of a single block in the x/y/z directions)
#ifndef BLOCK_SIZE
#  define BLOCK_SIZE       PATCH_SIZE
#endif



//-------------------------------------------------------------------------------------------------------
// Structure   :  block_t
// Description :  Data structure of a single block
//
// Data Member :
//
// Method      :  block_t     : Constructor
//                ~block_t    : Destructor
//-------------------------------------------------------------------------------------------------------
struct block_t
{

// data members
// ===================================================================================
   real (*dens)[BLOCK_SIZE][BLOCK_SIZE];  // density array

   int  father, son;    // ID of the father and son block
   long data_pos;       // data offset in the data file
   int  corner[3];      // bottom left "scale" of the block
   real edge[3];        // bottom left coordinates of the block



   //===================================================================================
   // Constructor :  block_t
   // Description :  Constructor of the structure "block_t"
   //
   // Note        :  Initialize data members
   //
   // Parameter   :
   //===================================================================================
   block_t()
   {
//    dens     = NULL;
      data_pos = -1;
   }



   //===================================================================================
   // Destructor  :  ~block_t
   // Description :  Destructor of the structure "block_t"
   //
   // Note        :  Deallocate data members
   //===================================================================================
   ~block_t()
   {
//    mfree();
   }



   /*
   //===================================================================================
   // Constructor :  mnew
   // Description :  Allocate memory for the data array
   //
   // Note        :
   //
   // Parameter   :
   //===================================================================================
   void mnew()
   {
      if ( dens != NULL )
      {
         fprintf( stderr, "WARNING : allocating an existing dens array !!\n" );
         mfree();
      }

      dens = new real [BLOCK_SIZE][BLOCK_SIZE][BLOCK_SIZE];
   }



   //===================================================================================
   // Constructor :  mfree
   // Description :  Deallocate memory for the data array
   //
   // Note        :
   //
   // Parameter   :
   //===================================================================================
   void mfree()
   {
      if ( dens != NULL )  delete [] dens;
      dens = NULL;
   }
   */


}; // struct block_t




//-------------------------------------------------------------------------------------------------------
// Structure   :  tree_t
// Description :  Octree structure
//
// Data Member :
//
// Method      :  pnew  : Allocate all blocks
//                pfree : Deallocate all blocks
//-------------------------------------------------------------------------------------------------------
struct tree_t
{

// data members
// ===================================================================================
   int    nlv, data_id, nblock_tot;    // number of AMR levels / data ID / total number of blocks at all levels
   long   step;                        // simulation step
   double time, box_size[3];           // simulation time and box size
   int    *nblock;                     // number of blocks at each level
   double *dh;                         // cell size at each level

   block_t **block;                    // block pointers



   //===================================================================================
   // Constructor :  tree_t
   // Description :  Constructor of the structure "tree_t"
   //
   // Note        :  Initialize data members
   //===================================================================================
   tree_t( const int nlv_in, const int data_id_in, const double time_in, const long step_in,
           const int nblock_in[], const double dh_in[], const double box_size_in[] )
   {
      nlv        = nlv_in;
      data_id    = data_id_in;
      time       = time_in;
      step       = step_in;
      nblock_tot = 0;
      nblock     = new int    [nlv];
      dh         = new double [nlv];
      block      = NULL;

      for (int lv=0; lv<nlv; lv++)
      {
         nblock[lv] = nblock_in[lv];
         dh    [lv] = dh_in[lv];

         nblock_tot += nblock[lv];
      }

      for (int d=0; d<3; d++)    box_size[d] = box_size_in[d];

      pnew();
   }



   //===================================================================================
   // Denstructor :  ~tree_t
   // Description :  Destructor of the structure "tree_t"
   //
   // Note        :  Free data members
   //===================================================================================
   ~tree_t()
   {
      delete [] nblock;
      delete [] dh;
      pfree();
   }



   //===================================================================================
   // Method      :  pnew
   // Description :  Allocate all blocks
   //
   // Parameter   :
   //===================================================================================
   void pnew()
   {
      int nblock_sum = 0;

      for (int lv=0; lv<nlv; lv++)  nblock_sum += nblock[lv];

      if ( block != NULL )
      {
         fprintf( stderr, "ERROR : allocating existing blocks !!\n" );
         pfree();
      }

      block    = new block_t*[nlv];
      block[0] = new block_t[nblock_sum];

      for (int lv=1; lv<nlv; lv++)  block[lv] = block[lv-1] + nblock[lv-1];
   }



   //===================================================================================
   // Method      :  pfree
   // Description :  Deallocate all blocks
   //
   // Parameter   :
   //===================================================================================
   void pfree()
   {
      if ( block[0] != NULL )
      {
         delete [] block[0];
         block[0] = NULL;
      }

      if ( block != NULL )
      {
         delete [] block;
         block = NULL;
      }
   }

}; // struct tree_t



#endif // #ifndef __TREE_H__


