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


#define WRONG              -999999


// always turn on gravity
#define GRAVITY


// models
#define HYDRO              1
#define MHD                2
#define ELBDM              3


#ifndef NULL
#define NULL               0
#endif


#ifndef NULL_INT
#  define NULL_INT         __INT_MAX__
#endif


#ifndef NULL_REAL
#  define NULL_REAL        __FLT_MAX__
#endif


#ifndef __FLT_MIN__
#  define __FLT_MIN__      1.17549435e-38F
#endif


#ifndef __FLT_MAX__
#  define __FLT_MAX__      3.40282347e+38F
#endif


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE         8
#define PS1                PATCH_SIZE


// number of components in each cell and the variable indices in the array "fluid"
#if   ( MODEL == HYDRO )
#  define NCOMP_FLUID      5
#  define DENS             0
#  define MOMX             1
#  define MOMY             2
#  define MOMZ             3
#  define ENGY             4

#  define _DENS            ( 1L << (DENS) )
#  define _MOMX            ( 1L << (MOMX) )
#  define _MOMY            ( 1L << (MOMY) )
#  define _MOMZ            ( 1L << (MOMZ) )
#  define _ENGY            ( 1L << (ENGY) )

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  define NCOMP_FLUID      8

#elif ( MODEL == ELBDM )
#  define NCOMP_FLUID      3
#  define DENS             0
#  define REAL             1
#  define IMAG             2

#  define _DENS            ( 1L << (DENS) )
#  define _REAL            ( 1L << (REAL) )
#  define _IMAG            ( 1L << (IMAG) )

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL

#  define NCOMP_TOTAL      ( NCOMP_FLUID + NCOMP_PASSIVE )


// symbolic constants used by all models
#  define _FLUID           (  ( 1L << NCOMP_FLUID ) - 1L           )
#  define _PASSIVE         (  ( 1L << NCOMP_TOTAL ) - 1L - _FLUID  )
#  define _TOTAL           (  ( 1L << NCOMP_TOTAL ) - 1L           )
#  define _POTE            ( 1L << (NCOMP_TOTAL  ) )
#  define _PAR_DENS        ( 1L << (NCOMP_TOTAL+1) )


// 3D to 1D array indices transformation
#define IDX321( i, j, k, Ni, Nj )   (  ( (k)*(Nj) + (j) )*(Ni) + (i)  )


// macro for the function "Aux_Error"
#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__


// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define   LOG( a )         log( a )
#  define  ATAN( a )        atan( a )
#  define  FMAX( a, b )     fmax( a, b )
#  define  FMIN( a, b )     fmin( a, b )
#  define   POW( a, b )      pow( a, b )
#  define  FMOD( a, b )     fmod( a, b )
#  define ATAN2( a, b )    atan2( a, b )
#else
#  define  FABS( a )        fabsf( a )
#  define  SQRT( a )        sqrtf( a )
#  define   SIN( a )         sinf( a )
#  define   COS( a )         cosf( a )
#  define   LOG( a )         logf( a )
#  define  ATAN( a )        atanf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define   POW( a, b )      powf( a, b )
#  define  FMOD( a, b )     fmodf( a, b )
#  define ATAN2( a, b )    atan2f( a, b )
#endif


// sign function
#define SIGN( a )       (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )


// square/cube function
#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// interpolation scheme
enum IntScheme_t { INT_DEFAULT=-1, INT_MINMOD3D=1, INT_MINMOD1D=2, INT_VANLEER=3, INT_CQUAD=4, INT_QUAD=5,
                   INT_CQUAR=6, INT_QUAR=7 };

// options in "Prepare_PatchGroupData"
enum NSide_t    { NSIDE_06=6, NSIDE_26=26 };


void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t
// Description :  data structure of a single patch
//
// Data Member :  fluid       : fluid variables (mass density, momentum density x, y ,z, energy density)
//                pot         : potential
//                par_dens    : particle (or total) density on grids
//                corner[3]   : physical coordinates of the patch corner
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

         fluid   [0][0][0][0] = -1;
         pot     [0][0][0]    = -1;
         par_dens[0][0][0]    = -1;
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



   //===================================================================================
   // Method      :  dnew
   // Description :  Allocate the data array
   //
   // Note        :  An error message will be displayed if array has already been allocated
   //===================================================================================
   void dnew()
   {
#     ifdef GAMER_DEBUG
      if ( fluid != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing fluid array !!\n" );

      if ( pot != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing pot array !!\n" );

      if ( par_dens != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing par_dens array !!\n" );
#     endif

      fluid    = new real [NCOMP_TOTAL][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
      pot      = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
      par_dens = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];

      fluid   [0][0][0][0] = -1;
      pot     [0][0][0]    = -1;
      par_dens[0][0][0]    = -1;
   } // METHOD : dnew


}; // struct patch_t





//-------------------------------------------------------------------------------------------------------
// Structure   :  AMR_t
// Description :  data structure of GAMER
//
// Data Member :  patch    : pointer of each patch
//                num      : number of patches (real patch + buffer patch) at each level
//                scale    : grid size at each level (normalize to the grid size at the finest level
//                dh       : Grid size at each level
//                BoxSize  : Simulation box size
//                BoxScale : Simulation box scale
//
// Method      :  pnew     : allocate one patch
//                pdelete  : deallocate one patch
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
   }



   //===================================================================================
   // Method      :  pnew
   // Description :  allocate a single patch
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
         fprintf( stderr, "ERROR : \"allocate an existing patch !!\" : Lv %d, PID %d\n", lv, num[lv] );
         exit(-1);
      }

      patch[lv][ num[lv] ] = new patch_t( x, y, z, FaPID, Data );

      num[lv] ++;
   }



   //===================================================================================
   // Method      :  pdelete
   // Description :  deallocate a single patch
   //
   // Parameter   :  lv  : the targeted refinement level
   //                PID : the patch ID to be removed
   //===================================================================================
   void pdelete( const int lv, const int PID )
   {
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



#endif // #ifndef __TYPEDEF_H__


