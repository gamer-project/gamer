#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif

#ifdef FLOAT8_PAR
typedef double real_par;
#else
typedef float  real_par;
#endif

#ifdef INT8_PAR
typedef long long_par;
#else
typedef int  long_par;
#endif


// models
#define HYDRO              1
#define MHD                2
#define ELBDM              3


#ifndef NULL
#define NULL               0
#endif


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE         8
#define PS1                (1*PATCH_SIZE)
#define PS2                (2*PATCH_SIZE)


// number of components in each cell and the variable indices in the array "fluid"
#if   ( MODEL == HYDRO )
#  define NCOMP_FLUID      5
#  define DENS             0
#  define MOMX             1
#  define MOMY             2
#  define MOMZ             3
#  define ENGY             4

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  define NCOMP_FLUID      8

#elif ( MODEL == ELBDM )
#  define NCOMP_FLUID      3
#  define DENS             0
#  define REAL             1
#  define IMAG             2

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL

#  define NCOMP_TOTAL      ( NCOMP_FLUID + NCOMP_PASSIVE )


// sibling index offset for the non-periodic B.C.
#define SIB_OFFSET_NONPERIODIC   ( -100 )


// macro for the function "Aux_Error"
#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__


// interpolation scheme
enum IntScheme_t { INT_DEFAULT=-1, INT_MINMOD3D=1, INT_MINMOD1D=2, INT_VANLEER=3, INT_CQUAD=4, INT_QUAD=5,
                   INT_CQUAR=6, INT_QUAR=7 };


// sign function
#define SIGN( a )       (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )


// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define  FMAX( a, b )     fmax( a, b )
#  define  FMIN( a, b )     fmin( a, b )
#  define   POW( a, b )      pow( a, b )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define ATAN2( a, b )    atan2( a, b )
#else
#  define  FABS( a )        fabsf( a )
#  define  SQRT( a )        sqrtf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define   POW( a, b )      powf( a, b )
#  define   SIN( a )         sinf( a )
#  define   COS( a )         cosf( a )
#  define ATAN2( a, b )    atan2f( a, b )
#endif


// square/cube function
#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// initial value for some variables (must be negative)
#define WRONG              -9999999


// maximum length of strings
#define MAX_STRING         512




//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t
// Description :  data structure of a single patch
//
// Data Member :  fluid       : fluid variables (mass density, momentum density x, y ,z, energy density)
//                pot         : potential
//                par_dens    : particle density on grids
//                corner[3]   : physical coordinates of the patch corner
//                sibling[26] : patch IDs of the 26 sibling patches
//                father      : patch ID of the father patch
//                son         : patch ID of the child patch
//
// Method      :  patch_t     : constructor
//                ~patch_t    : destructor
//-------------------------------------------------------------------------------------------------------
struct patch_t
{

// data members
// ===================================================================================
   real (*fluid   )[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
   real (*pot     )[PATCH_SIZE][PATCH_SIZE];
   real (*par_dens)[PATCH_SIZE][PATCH_SIZE];

   int  corner[3];
   int  sibling[26];
   int  father;
   int  son;



   //===================================================================================
   // Constructor :  patch_t
   // Description :  constructor of the structure "patch_t"
   //
   // Note        :  initialize data members
   //
   // Parameter   :  x,y,z : physical coordinates of the patch corner
   //                FaPID : patch ID of the father patch
   //                Data  : true --> allocate physical data (fluid + pot + par_dens)
   //===================================================================================
   patch_t( const int x, const int y, const int z, const int FaPID, const bool Data )
   {
      corner[0] = x;
      corner[1] = y;
      corner[2] = z;
      father    = FaPID;
      son       = -1;

      for (int s=0; s<26; s++ )     sibling[s] = -1;     // -1 <--> NO sibling

      fluid    = NULL;
      pot      = NULL;
      par_dens = NULL;

      if ( Data )
      {
         fluid    = new real [NCOMP_TOTAL][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
         pot      = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
         par_dens = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];

         fluid[0][0][0][0] = -1;
      }
   }



   //===================================================================================
   // Destructor  :  ~patch_t
   // Description :  destructor of the structure "patch_t"
   //
   // Note        :  deallocate flux and data arrays
   //===================================================================================
   ~patch_t()
   {
      if ( fluid != NULL )
      {
         delete [] fluid;
         fluid = NULL;
      }

      if ( pot != NULL )
      {
         delete [] pot;
         pot = NULL;
      }

      if ( par_dens != NULL )
      {
         delete [] par_dens;
         par_dens = NULL;
      }
   }


}; // struct patch_t





//-------------------------------------------------------------------------------------------------------
// Structure   :  AMR_t
// Description :  data structure of GAMER
//
// Data Member :  patch    : Pointer of each patch
//                num      : Number of patches (real patch + buffer patch) at each level
//                scale    : Grid size at each level (normalize to the grid size at the finest level
//                dh       : Grid size at each level
//                BoxSize  : Simulation box size
//                BoxScale : Simulation box scale
//
// Method      :  pnew    : Allocate one patch
//                pdelete : Deallocate one patch
//-------------------------------------------------------------------------------------------------------
struct AMR_t
{

// data members
// ===================================================================================
   patch_t *patch[NLEVEL][MAX_PATCH];

   int    num  [NLEVEL];
   int    scale[NLEVEL];
   double dh   [NLEVEL];
   double BoxSize [3];
   int    BoxScale[3];



   //===================================================================================
   // Constructor :  AMR_t
   // Description :  constructor of the structure "AMR_t"
   //
   // Note        :  initialize data members
   //===================================================================================
   AMR_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
         num  [lv] = 0;
         scale[lv] = 1<<(NLEVEL-1-lv);
      }

      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<MAX_PATCH; PID++)
         patch[lv][PID] = NULL;

      for (int d=0; d<3; d++)
      {
         BoxSize [d] = -1.0;
         BoxScale[d] = -1;
      }
   }



   //===================================================================================
   // Method      :  pnew
   // Description :  allocate a single patch
   //
   // Note        :  a. each patch contains two patch pointers --> SANDGLASS (Sg) = 0 / 1
   //                b. Sg = 0 : store both data and relation (father,son.sibling,corner,flag,flux)
   //                   Sg = 1 : store only data
   //
   // Parameter   :  lv    : the targeted refinement level
   //                x,y,z : physical coordinates of the patch corner
   //                FaPID : the patch ID of the parent patch at level "lv-1"
   //                Data  : true --> allocate physical data (fluid + pot + par_dens)
   //===================================================================================
   void pnew( const int lv, const int x, const int y, const int z, const int FaPID, const bool Data )
   {
      if ( num[lv] == MAX_PATCH )
      {
         fprintf( stderr, "ERROR : \"exceeding the maximum number of patches [%d] !!\n", MAX_PATCH );
         exit(-1);
      }

      if ( patch[lv][num[lv]] != NULL )
      {
         fprintf( stderr, "ERROR : \"allocate an existing patch !!\" : Lv %d, PID %d, FaPID %d\n",
                  lv, num[lv], FaPID );
         exit(-1);
      }

      patch[lv][ num[lv] ] = new patch_t( x, y, z, FaPID, Data );

      num[lv] ++;

      if ( num[lv] > MAX_PATCH )
      {
         fprintf( stderr, "ERROR : exceed MAX_PATCH !!\n" );
         exit(-1);
      }
   }



   //===================================================================================
   // Method      :  pdelete
   // Description :  deallocate a single patch
   //
   // Note        :  a. this function should NOT be applied to the base-level patches
   //                b. this function will also deallocate the flux arrays of the targeted patch
   //                c. delete a patch with son is forbidden
   //
   // Parameter   :  lv  : the targeted refinement level
   //                PID : the patch ID to be removed
   //===================================================================================
   void pdelete( const int lv, const int PID )
   {
      if ( lv == 0 )
      {
         fprintf( stderr, "ERROR : delete a base-level patch !!\n" );
         exit(-1);
      }

      if ( patch[lv][PID]->son != -1 )
      {
         fprintf( stderr, "ERROR : \"delete a patch with son !!\" : Lv %d, PID %d, SonPID %d\n",
                  lv, PID, patch[lv][PID]->son );
         exit(-1);
      }

      delete patch[lv][PID];

      patch[lv][PID] = NULL;

      num[lv] --;

      if ( num[lv] < 0 )
      {
         fprintf( stderr, "ERROR : num[%d] = %d < 0 !!\n", lv, num[lv] );
         exit(-1);
      }
   }


}; // struct AMR_t





//-------------------------------------------------------------------------------------------------------
// Structure   :  ParaVar_t
// Description :  data structure collecting variables related to GAMER parallelization
//
// Data Member :  BounP_NList    : the number of boundary patches in 26 sibling directions
//                BounP_IDList   : the IDs of boundary patches in 26 sibling directions
//                BounP_PosList  : the positions of boundary patches recorded in "BounP_IDList"
//
// Method      :
//-------------------------------------------------------------------------------------------------------
struct ParaVar_t
{

// data members
// ===================================================================================
   int  BounP_NList      [NLEVEL  ][26];
   int *BounP_IDList     [NLEVEL  ][26];
   int *BounP_PosList    [NLEVEL  ][26];

   int  SendP_NList      [NLEVEL  ][26];
   int *SendP_IDList     [NLEVEL  ][26];
   int  RecvP_NList      [NLEVEL  ][26];
   int *RecvP_IDList     [NLEVEL  ][26];



   //===================================================================================
   // Constructor :  ParaVar_t
   // Description :  constructor of the structure "ParaVar_t"
   //
   // Note        :  initialize all pointers as NULL and all counters as 0
   //===================================================================================
   ParaVar_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      for (int s=0; s<26; s++)
      {
         BounP_NList       [lv][s] = 0;
         BounP_IDList      [lv][s] = NULL;
         BounP_PosList     [lv][s] = NULL;

         SendP_NList       [lv][s] = 0;
         SendP_IDList      [lv][s] = NULL;
         RecvP_NList       [lv][s] = 0;
         RecvP_IDList      [lv][s] = NULL;
      }
   }


}; // struct ParaVar_t



#endif // #ifndef __TYPEDEF_H__
