#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif


// models
#define HYDRO              1
#define MHD                2
#define ELBDM              3
#define SR_HYDRO           5


#ifndef NULL
#define NULL               0
#endif


#ifndef __FLT_MIN__
#  define __FLT_MIN__      1.17549435e-38F
#endif


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE         8


// number of components in each cell and the variable indices in the array "fluid"
#if   ( MODEL == HYDRO || MODEL == SR_HYDRO )
#  define NCOMP_FLUID      5
#  define DENS             0
#  define MOMX             1
#  define MOMY             2
#  define MOMZ             3
#  define ENGY             4

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  define NCOMP_FLUID      5

#elif ( MODEL == ELBDM )
#  define NCOMP_FLUID      3
#  define DENS             0
#  define REAL             1
#  define IMAG             2

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL

#  define NCOMP_TOTAL      ( NCOMP_FLUID + NCOMP_PASSIVE )



//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t
// Description :  data structure of a single patch
//
// Data Member :  fluid	      : fluid variables (mass density, momentum density x, y ,z, energy density)
//		  pot	      : potential
//		  par_dens    : particle density deposited onto grids
//		  corner[3]   : physical coordinates of the patch corner
//		  father      : patch ID of the father patch
//		  son	      : patch ID of the child patch
//
// Method      :  patch_t     : constructor
//		  ~patch_t    : destructor
//-------------------------------------------------------------------------------------------------------
struct patch_t
{

// data members
// ===================================================================================
   real (*fluid   )[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
   real (*pot     )[PATCH_SIZE][PATCH_SIZE];
   real (*par_dens)[PATCH_SIZE][PATCH_SIZE];

   int	corner[3];
   int 	father;
   int  son;
   bool check;



   //===================================================================================
   // Constructor :  patch_t
   // Description :  constructor of the structure "patch_t"
   //
   // Note	  :  initialize data members
   //
   // Parameter   :  x,y,z : physical coordinates of the patch corner
   //		     FaPID : patch ID of the father patch
   //		     Data  : true --> allocate physical data (fluid + pot + par_dens)
   //===================================================================================
   patch_t( const int x, const int y, const int z, const int FaPID, const bool Data )
   {
      corner[0] = x;
      corner[1] = y;
      corner[2] = z;
      father    = FaPID;
      son       = -1;
      check     = false;

      fluid    = NULL;
      pot      = NULL;
      par_dens = NULL;

      if ( Data )
      {
	 fluid    = new real [NCOMP_TOTAL][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
	 pot      = new real              [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
	 par_dens = new real              [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];

	 fluid[0][0][0][0] = -1;
      }
   }

}; // struct patch_t





//-------------------------------------------------------------------------------------------------------
// Structure   :  GAMER_t
// Description :  data structure of GAMER
//
// Data Member :  ptr		 : pointer of each patch
//		  num		 : number of patches (real patch + buffer patch) at each level
//		  scale 	 : grid size at each level (normalize to the grid size at the finest level
//
// Method      :  pnew		 : allocate one patch
//		  pdelete	 : deallocate one patch
//-------------------------------------------------------------------------------------------------------
struct GAMER_t
{

// data members
// ===================================================================================
   patch_t *ptr[NLEVEL][MAX_PATCH];

   int num[NLEVEL];
   int scale[NLEVEL];
   int nx0_tot[3];
   int ngpu_x[3];



   //===================================================================================
   // Constructor :  GAMER_t
   // Description :  constructor of the structure "GAMER_t"
   //
   // Note	  :  initialize data members
   //===================================================================================
   GAMER_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
         num[lv]	    = 0;
         scale[lv]	    = 1<<(NLEVEL-1-lv);
      }

      for (int d=0; d<3; d++)
      {
	 nx0_tot[d] = 0;
	 ngpu_x [d] = 0;
      }

      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<MAX_PATCH; PID++)
	 ptr[lv][PID] = NULL;
   }



   //===================================================================================
   // Method	  :  pnew
   // Description :  allocate a single patch
   //
   // Parameter   :  lv	   : the targeted refinement level
   //		     x,y,z : physical coordinates of the patch corner
   //		     FaPID : the patch ID of the parent patch at level "lv-1"
   //		     Data  : true --> allocate physical data (fluid + pot + par_dens)
   //===================================================================================
   void pnew( const int lv, const int x, const int y, const int z, const int FaPID, const bool Data )
   {
      if ( ptr[lv][num[lv]] != NULL )
      {
	 fprintf( stderr, "ERROR : \"allocate an existing patch !!\" : Lv %d, PID %d, FaPID %d\n",
		  lv, num[lv], FaPID );
	 exit(-1);
      }

      ptr[lv][ num[lv] ] = new patch_t( x, y, z, FaPID, Data );

      num[lv] ++;

      if ( num[lv] > MAX_PATCH )
      {
         fprintf( stderr, "ERROR : exceed MAX_PATCH !!\n" );
         exit(-1);
      }
   }

}; // struct GAMER_t



#endif // #ifndef __TYPEDEF_H__


