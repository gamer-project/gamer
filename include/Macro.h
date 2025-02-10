#ifndef __MACRO_H__
#define __MACRO_H__



// ****************************************************************************
// ** This header defines the symbolic constants and macros used in GAMER.   **
// ** For clarity, useless options defined in the makefile will be "undef"   **
// ** in the end of this file.                                               **
// ****************************************************************************


// ########################
// ## Symbolic Constants ##
// ########################

// current version
#define VERSION      "gamer-2.1.1.dev"


// option == NONE --> the option is turned off
#define NONE         0


// GPU architecture
#define FERMI        1
#define KEPLER       2
#define MAXWELL      3
#define PASCAL       4
#define VOLTA        5
#define TURING       6
#define AMPERE       7
#define ADA_LOVELACE 8
#define HOPPER       9

#ifdef GPU
#if   ( GPU_COMPUTE_CAPABILITY >= 200  &&  GPU_COMPUTE_CAPABILITY < 300 )
# define GPU_ARCH FERMI
#elif ( GPU_COMPUTE_CAPABILITY >= 300  &&  GPU_COMPUTE_CAPABILITY < 500 )
# define GPU_ARCH KEPLER
#elif ( GPU_COMPUTE_CAPABILITY >= 500  &&  GPU_COMPUTE_CAPABILITY < 600 )
# define GPU_ARCH MAXWELL
#elif ( GPU_COMPUTE_CAPABILITY >= 600  &&  GPU_COMPUTE_CAPABILITY < 700 )
# define GPU_ARCH PASCAL
#elif ( GPU_COMPUTE_CAPABILITY >= 700  &&  GPU_COMPUTE_CAPABILITY < 750 )
# define GPU_ARCH VOLTA
#elif ( GPU_COMPUTE_CAPABILITY >= 750  &&  GPU_COMPUTE_CAPABILITY < 800 )
# define GPU_ARCH TURING
#elif ( GPU_COMPUTE_CAPABILITY >= 800  &&  GPU_COMPUTE_CAPABILITY < 890 )
# define GPU_ARCH AMPERE
#elif ( GPU_COMPUTE_CAPABILITY >= 890  &&  GPU_COMPUTE_CAPABILITY < 900 )
# define GPU_ARCH ADA_LOVELACE
#elif ( GPU_COMPUTE_CAPABILITY >= 900  &&  GPU_COMPUTE_CAPABILITY < 1000 )
# define GPU_ARCH HOPPER
#else
# error : ERROR : Unknown GPU_COMPUTE_CAPABILITY !!
#endif // GPU_COMPUTE_CAPABILITY
#endif // #ifdef GPU


// models
#define HYDRO        1
//#define MHD        2     // MHD is now regarded as an option of HYDRO
#define ELBDM        3
#define PAR_ONLY     4


// hydrodynamic schemes
#define RTVD         1
#define MHM          3
#define MHM_RP       4
#define CTU          5


// data reconstruction schemes
#define PLM          1
#define PPM          2


// Riemann solvers
#define EXACT        1
#define ROE          2
#define HLLE         3
#define HLLC         4
#define HLLD         5


// dual-energy variables
#define DE_ENPY      1
#define DE_EINT      2

#ifdef DUAL_ENERGY
#define DE_UPDATED_BY_ETOT       ('0')
#define DE_UPDATED_BY_DUAL       ('1')
#define DE_UPDATED_BY_MIN_PRES   ('2')
#define DE_UPDATED_BY_ETOT_GRA   ('3')
#endif


// equation of states
#define EOS_GAMMA       1
#define EOS_ISOTHERMAL  2
#define EOS_NUCLEAR     3
#define EOS_TABULAR     4
#define EOS_COSMIC_RAY  5
#define EOS_TAUBMATHEWS 6
#define EOS_USER        7


// ELBDM schemes
#define ELBDM_WAVE      1
#define ELBDM_HYBRID    2


// ELBDM wave schemes
#define WAVE_FD         1
#define WAVE_GRAMFE     2


// ELBDM hybrid schemes
#define HYBRID_UPWIND   1
#define HYBRID_MUSCL    2
#define HYBRID_FROMM    3


// ELBDM gramfe schemes
#define GRAMFE_FFT      1
#define GRAMFE_MATMUL   2


// Poisson solvers
#define SOR          1
#define MG           2


// load-balance parallelization
#define HILBERT      1

// random number implementation
#define RNG_GNU_EXT  1
#define RNG_CPP11    2


// NCOMP_FLUID : number of active components in each cell (for patch->fluid[])
//               --> do not include passive components here, which is set by NCOMP_PASSIVE
// NFLUX_FLUID : number of active components in patch->flux[]
//               --> do not include passive components here, which is set by NFLUX_PASSIVE
// NCOMP_MAG   : number of magnetic field components (for patch->magnetic[])
// NCOMP_ELE   : number of electric field components on each cell face (for patch->electric[])
#if   ( MODEL == HYDRO )
#  define NCOMP_FLUID         5
#  define NFLUX_FLUID         NCOMP_FLUID
# ifdef MHD
#  define NCOMP_MAG           3
#  define NCOMP_ELE           2
# else
#  define NCOMP_MAG           0
#  define NCOMP_ELE           0
# endif

// for ELBDM, we only need the density flux
#elif ( MODEL == ELBDM )
#  define NCOMP_FLUID         3
#  define NFLUX_FLUID         1
#  define NCOMP_MAG           0
#  define NCOMP_ELE           0

#elif ( MODEL == PAR_ONLY )
#  define NCOMP_FLUID         0
#  define NFLUX_FLUID         0
#  define NCOMP_MAG           0
#  define NCOMP_ELE           0

#else
#  error : ERROR : unsupported MODEL (please edit NCOMP_FLUID and NFLUX_FLUID for the new MODEL) !!
#endif // MODEL


// number of passively advected components in each cell

// define NCOMP_PASSIVE_USER if not set in the Makefile
#ifndef NCOMP_PASSIVE_USER
#  define NCOMP_PASSIVE_USER  0
#endif

// add built-in scalars
#if ( MODEL == HYDRO )

// dual-energy variable
# ifdef DUAL_ENERGY
#  define NCOMP_PASSIVE_BUILTIN0    1
# else
#  define NCOMP_PASSIVE_BUILTIN0    0
# endif

// cosmic rays
# ifdef COSMIC_RAY
#  define NCOMP_PASSIVE_BUILTIN1    1
# else
#  define NCOMP_PASSIVE_BUILTIN1    0
# endif

// total number of built-in scalars
#  define NCOMP_PASSIVE_BUILTIN     ( NCOMP_PASSIVE_BUILTIN0 + NCOMP_PASSIVE_BUILTIN1 )

#endif // #if ( MODEL == HYDRO )

// define NCOMP_PASSIVE_BUILTIN if not set yet
#ifndef NCOMP_PASSIVE_BUILTIN
#  define NCOMP_PASSIVE_BUILTIN     0
#endif

// total number of passive scalars
#  define NCOMP_PASSIVE       ( NCOMP_PASSIVE_USER + NCOMP_PASSIVE_BUILTIN )

// assuming all passive scalars have the corresponding fluxes
#  define NFLUX_PASSIVE       NCOMP_PASSIVE


// total number of variables in each cell and in the flux array including both active and passive variables
#  define NCOMP_TOTAL         ( NCOMP_FLUID + NCOMP_PASSIVE )
#  define NFLUX_TOTAL         ( NFLUX_FLUID + NFLUX_PASSIVE )


// maximum number of reference values stored in ConRef[]
#  define NCONREF_MAX         60


// number of input/output fluid variables in the fluid solver
#if   ( MODEL == HYDRO )
#  define FLU_NIN             NCOMP_TOTAL
#  define FLU_NOUT            NCOMP_TOTAL


// for ELBDM, we do not need to transfer the density component into GPU;
// also exclude passive scalars for now since it is not supported yet
// --> consistent with excluding _PASSIVE when calling Prepare_PatchData() in Flu_Prepare.cpp
#elif ( MODEL == ELBDM )
//#  define FLU_NIN             ( NCOMP_TOTAL - 1 )
//#  define FLU_NOUT            ( NCOMP_TOTAL - 0 )
#  define FLU_NIN             ( NCOMP_FLUID - 1 )
#  define FLU_NOUT            ( NCOMP_FLUID - 0 )

#elif ( MODEL == PAR_ONLY )
#  define FLU_NIN             0
#  define FLU_NOUT            0

#else
#  error : ERROR : unsupported MODEL (please edit FLU_NIN and FLU_NOUT for the new MODEL) !!
#endif // MODEL


// number of input fluid variables in the dt solver
// --> EOS_GAMMA/EOS_ISOTHERMAL do not require passive scalars
#if (  MODEL == HYDRO  &&  !defined SRHD  &&  ( EOS == EOS_GAMMA || EOS == EOS_ISOTHERMAL )  )
#  define FLU_NIN_T           NCOMP_FLUID
#else
#  define FLU_NIN_T           NCOMP_TOTAL
#endif


// number of input/output fluid variables in the source-term solver
// --> fixed to NCOMP_TOTAL for now
#  define FLU_NIN_S           NCOMP_TOTAL
#  define FLU_NOUT_S          NCOMP_TOTAL


// maximum number of output derived fields
#  define DER_NOUT_MAX        10


// maximum number of fields to be stored in HDF5 snapshots
#  define NFIELD_STORED_MAX   50


// built-in fields in different models
#if   ( MODEL == HYDRO )
// field indices of fluid[] --> element of [0 ... NCOMP_FLUID-1]
// --> must NOT modify their values
// --> in addition, they must be consistent with the order these fields are declared in Init_Field()
#  define DENS                0
#  define MOMX                1
#  define MOMY                2
#  define MOMZ                3
#  define ENGY                4

// field indices of passive[] --> element of [NCOMP_FLUID ... NCOMP_TOTAL-1]
#if ( NCOMP_PASSIVE > 0 )

// always put the built-in variables at the END of the field list
// --> so that their indices (e.g., DUAL/CRAY) can be determined during compilation
// --> convenient (and probably also more efficient) for the fluid solver
#  define PASSIVE_NEXT_IDX0   ( NCOMP_TOTAL - 1   )

# ifdef DUAL_ENERGY
#  define DUAL                ( PASSIVE_NEXT_IDX0 )
#  define PASSIVE_NEXT_IDX1   ( DUAL - 1          )
# else
#  define PASSIVE_NEXT_IDX1   ( PASSIVE_NEXT_IDX0 )
# endif

# ifdef COSMIC_RAY
#  define CRAY                ( PASSIVE_NEXT_IDX1 )
#  define PASSIVE_NEXT_IDX2   ( CRAY - 1          )
# else
#  define PASSIVE_NEXT_IDX2   ( PASSIVE_NEXT_IDX1 )
# endif

#endif // #if ( NCOMP_PASSIVE > 0 )

// field indices of magnetic --> element of [0 ... NCOMP_MAG-1]
# ifdef MHD
#  define MAGX                0
#  define MAGY                1
#  define MAGZ                2
# endif

// flux indices of flux[] --> element of [0 ... NFLUX_FLUID-1]
#  define FLUX_DENS           0
#  define FLUX_MOMX           1
#  define FLUX_MOMY           2
#  define FLUX_MOMZ           3
#  define FLUX_ENGY           4

// flux indices of flux_passive[] --> element of [NFLUX_FLUID ... NFLUX_TOTAL-1]
#if ( NCOMP_PASSIVE > 0 )

// always put the built-in variables at the END of the list
#  define FLUX_NEXT_IDX0   ( NFLUX_TOTAL - 1 )

# ifdef DUAL_ENERGY
#  define FLUX_DUAL        ( FLUX_NEXT_IDX0  )
#  define FLUX_NEXT_IDX1   ( FLUX_DUAL - 1   )
# else
#  define FLUX_NEXT_IDX1   ( FLUX_NEXT_IDX0  )
# endif

# ifdef COSMIC_RAY
#  define FLUX_CRAY        ( FLUX_NEXT_IDX1  )
#  define FLUX_NEXT_IDX2   ( FLUX_CRAY - 1   )
# else
#  define FLUX_NEXT_IDX2   ( FLUX_NEXT_IDX1  )
# endif

#endif // #if ( NCOMP_PASSIVE > 0 )

// bitwise field indices
// --> must have "_VAR_NAME = 1L<<VAR_NAME" (e.g., _DENS == 1L<<DENS)
// --> convenient for determining subsets of fields (e.g., _DENS|_ENGY)
// --> used as function parameters (e.g., Prepare_PatchData(), Flu_FixUp(), Flu_FixUp_Restrict(), Buf_GetBufferData())
#  define _DENS               ( 1L << DENS )
#  define _MOMX               ( 1L << MOMX )
#  define _MOMY               ( 1L << MOMY )
#  define _MOMZ               ( 1L << MOMZ )
#  define _ENGY               ( 1L << ENGY )

#if ( NCOMP_PASSIVE > 0 )

# ifdef DUAL_ENERGY
#  define _DUAL               ( 1L << DUAL )
# endif

# ifdef COSMIC_RAY
#  define _CRAY               ( 1L << CRAY )
# endif

#endif // #if ( NCOMP_PASSIVE > 0 )

// magnetic field
# ifdef MHD
#  define _MAGX               ( 1L << MAGX )
#  define _MAGY               ( 1L << MAGY )
#  define _MAGZ               ( 1L << MAGZ )
#  define _MAG                ( _MAGX | _MAGY | _MAGZ )
# else
#  define _MAG                0
# endif

// bitwise flux indices
#  define _FLUX_DENS          ( 1L << FLUX_DENS )
#  define _FLUX_MOMX          ( 1L << FLUX_MOMX )
#  define _FLUX_MOMY          ( 1L << FLUX_MOMY )
#  define _FLUX_MOMZ          ( 1L << FLUX_MOMZ )
#  define _FLUX_ENGY          ( 1L << FLUX_ENGY )

#if ( NFLUX_PASSIVE > 0 )

# ifdef DUAL_ENERGY
#  define _FLUX_DUAL          ( 1L << FLUX_DUAL )
# endif

# ifdef COSMIC_RAY
#  define _FLUX_CRAY          ( 1L << FLUX_CRAY )
# endif

#endif // #if ( NFLUX_PASSIVE > 0 )

// bitwise indices of derived fields
// --> start from (1L<<NCOMP_TOTAL) to distinguish from the intrinsic fields
// --> remember to define NDERIVE = total number of derived fields
#  define _VELX               ( 1L << (NCOMP_TOTAL+ 0) )
#  define _VELY               ( 1L << (NCOMP_TOTAL+ 1) )
#  define _VELZ               ( 1L << (NCOMP_TOTAL+ 2) )
#  define _VELR               ( 1L << (NCOMP_TOTAL+ 3) )
#  define _PRES               ( 1L << (NCOMP_TOTAL+ 4) )
#  define _TEMP               ( 1L << (NCOMP_TOTAL+ 5) )
#  define _ENTR               ( 1L << (NCOMP_TOTAL+ 6) )
#  define _EINT               ( 1L << (NCOMP_TOTAL+ 7) )
# ifdef MHD
#  define _MAGX_CC            ( 1L << (NCOMP_TOTAL+ 8) )
#  define _MAGY_CC            ( 1L << (NCOMP_TOTAL+ 9) )
#  define _MAGZ_CC            ( 1L << (NCOMP_TOTAL+10) )
#  define _MAGE_CC            ( 1L << (NCOMP_TOTAL+11) )
# else
#  define _MAGX_CC            0
#  define _MAGY_CC            0
#  define _MAGZ_CC            0
#  define _MAGE_CC            0
# endif // #ifdef MHD ... else ...
#  define _DERIVED            ( _VELX | _VELY | _VELZ | _VELR | _PRES | _TEMP | _ENTR | _EINT | _MAGX_CC | _MAGY_CC | _MAGZ_CC | _MAGE_CC )
#  define NDERIVE             12


#elif ( MODEL == ELBDM )
// field indices of fluid[] --> element of [0 ... NCOMP_FLUID-1]
// --> must NOT modify their values
#  define  DENS               0
#  define  REAL               1
#  define  IMAG               2

# if ( ELBDM_SCHEME == ELBDM_HYBRID )
#  define  PHAS               1
#  define  STUB               2
# endif

// field indices of passive[] --> element of [NCOMP_FLUID ... NCOMP_TOTAL-1]
// none for ELBDM

// flux indices of flux[] --> element of [0 ... NFLUX_FLUID-1]
#  define  FLUX_DENS          0

// bitwise field indices
#  define _DENS               ( 1L << DENS )
#  define _REAL               ( 1L << REAL )
#  define _IMAG               ( 1L << IMAG )
#  define _MAG                0
# if ( ELBDM_SCHEME == ELBDM_HYBRID )
#  define _PHAS               ( 1L << PHAS )
#  define _STUB               ( 1L << STUB )
# endif


// bitwise flux indices
// for the hybrid scheme, we also only need the density flux
#  define _FLUX_DENS          ( 1L << FLUX_DENS )

// bitwise indices of derived fields
#  define _DERIVED            0
#  define NDERIVE             0


#elif ( MODEL == PAR_ONLY )
#  define _MAG                0
#  define _DERIVED            0
#  define NDERIVE             0


#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// bitwise field indices used by all models
#  define _NONE               0
# ifdef GRAVITY
#  define _POTE               ( 1L << (NCOMP_TOTAL+NDERIVE) )
# endif
#  define _FLUID              (  ( 1L << NCOMP_FLUID ) - 1L           )
#  define _PASSIVE            (  ( 1L << NCOMP_TOTAL ) - 1L - _FLUID  )
#  define _TOTAL              (  ( 1L << NCOMP_TOTAL ) - 1L           )

#  define _FLUX_FLUID         (  ( 1L << NFLUX_FLUID ) - 1L                )
#  define _FLUX_PASSIVE       (  ( 1L << NFLUX_TOTAL ) - 1L - _FLUX_FLUID  )
#  define _FLUX_TOTAL         (  ( 1L << NFLUX_TOTAL ) - 1L                )



// symbolic constants for particles
#ifdef PARTICLE

// number of built-in particle attributes
// floating-point: mass, position*3, velocity*3, and time
// integer: type
#  define PAR_NATT_FLT_BUILTIN0   8
#  define PAR_NATT_INT_BUILTIN0   1

// acceleration*3 when STORE_PAR_ACC is adopted
# if ( defined STORE_PAR_ACC  &&  defined GRAVITY )
#  define PAR_NATT_FLT_BUILTIN1   3
# else
#  define PAR_NATT_FLT_BUILTIN1   0
# endif

// particle creation time when STAR_FORMATION is adopted
# ifdef STAR_FORMATION
#  define PAR_NATT_FLT_BUILTIN2   1
# else
#  define PAR_NATT_FLT_BUILTIN2   0
# endif

// **total** number of built-in particle attributes
#  define PAR_NATT_FLT_BUILTIN    ( PAR_NATT_FLT_BUILTIN0 + PAR_NATT_FLT_BUILTIN1 + PAR_NATT_FLT_BUILTIN2 )
#  define PAR_NATT_INT_BUILTIN    ( PAR_NATT_INT_BUILTIN0 )


// number of particle attributes that we do not want to store on disk (currently time + acceleration*3)
#  define PAR_NATT_FLT_UNSTORED   ( 1 + PAR_NATT_FLT_BUILTIN1 )
#  define PAR_NATT_FLT_STORED     ( PAR_NATT_FLT_TOTAL - PAR_NATT_FLT_UNSTORED )
#  define PAR_NATT_INT_UNSTORED   ( 0 )
#  define PAR_NATT_INT_STORED     ( PAR_NATT_INT_TOTAL - PAR_NATT_INT_UNSTORED )


// define PAR_NATT_FLT/INT_USER if not set in the Makefile
# ifndef PAR_NATT_FLT_USER
#  define PAR_NATT_FLT_USER       0
# endif
# ifndef PAR_NATT_INT_USER
#  define PAR_NATT_INT_USER       0
# endif


// total number of particle attributes (built-in + user-defined)
#  define PAR_NATT_FLT_TOTAL      ( PAR_NATT_FLT_BUILTIN + PAR_NATT_FLT_USER )
#  define PAR_NATT_INT_TOTAL      ( PAR_NATT_INT_BUILTIN + PAR_NATT_INT_USER )


// indices of built-in particle floating-point attributes in Par->AttributeFlt[]
// --> must NOT modify their values
#  define  PAR_MASS           0
#  define  PAR_POSX           1
#  define  PAR_POSY           2
#  define  PAR_POSZ           3
#  define  PAR_VELX           4
#  define  PAR_VELY           5
#  define  PAR_VELZ           6

// indices of built-in particle integer attributes in Par->AttributeInt[]
// --> must NOT modify their values
#  define  PAR_TYPE           0

// always put acceleration and time at the END of the particle attribute list
// --> make it easier to discard them when storing data on disk (see Output_DumpData_Total(_HDF5).cpp)
# if ( defined STORE_PAR_ACC  &&  defined GRAVITY )
#  define  PAR_ACCX           ( PAR_NATT_FLT_TOTAL - 4 )
#  define  PAR_ACCY           ( PAR_NATT_FLT_TOTAL - 3 )
#  define  PAR_ACCZ           ( PAR_NATT_FLT_TOTAL - 2 )
# endif
#  define  PAR_TIME           ( PAR_NATT_FLT_TOTAL - 1 )


// bitwise indices of particles
// particle attributes
#  define _PAR_MASS           ( 1L << PAR_MASS )
#  define _PAR_POSX           ( 1L << PAR_POSX )
#  define _PAR_POSY           ( 1L << PAR_POSY )
#  define _PAR_POSZ           ( 1L << PAR_POSZ )
#  define _PAR_VELX           ( 1L << PAR_VELX )
#  define _PAR_VELY           ( 1L << PAR_VELY )
#  define _PAR_VELZ           ( 1L << PAR_VELZ )
# if ( defined STORE_PAR_ACC  &&  defined GRAVITY )
#  define _PAR_ACCX           ( 1L << PAR_ACCX )
#  define _PAR_ACCY           ( 1L << PAR_ACCY )
#  define _PAR_ACCZ           ( 1L << PAR_ACCZ )
# endif
#  define _PAR_TIME           ( 1L << PAR_TIME )
#  define _PAR_POS            ( _PAR_POSX | _PAR_POSY | _PAR_POSZ )
#  define _PAR_VEL            ( _PAR_VELX | _PAR_VELY | _PAR_VELZ )
# if ( defined STORE_PAR_ACC  &&  defined GRAVITY )
#  define _PAR_ACC            ( _PAR_ACCX | _PAR_ACCY | _PAR_ACCZ )
# endif
#  define _PAR_FLT_TOTAL      (  ( 1L << PAR_NATT_FLT_TOTAL ) - 1L )

#  define _PAR_TYPE           ( 1L << PAR_TYPE )
#  define _PAR_INT_TOTAL      (  ( 1L << PAR_NATT_INT_TOTAL ) - 1L )

// grid fields related to particles
// --> note that _POTE = ( 1L << (NCOMP_TOTAL+NDERIVE) )
#  define _PAR_DENS           ( 1L << (NCOMP_TOTAL+NDERIVE+1) )

# if ( MODEL == PAR_ONLY )
#  define _TOTAL_DENS         ( _PAR_DENS )
# else
#  define _TOTAL_DENS         ( 1L << (NCOMP_TOTAL+NDERIVE+2) )
# endif

// particle type macros

// number of particle types (default: 4)
#  define  PAR_NTYPE                4

// particle type indices (must be in the range 0<=index<PAR_NTYPE)
#  define  PTYPE_TRACER             (long_par)0
#  define  PTYPE_GENERIC_MASSIVE    (long_par)1
#  define  PTYPE_DARK_MATTER        (long_par)2
#  define  PTYPE_STAR               (long_par)3

# ifdef GRAVITY
#  define MASSIVE_PARTICLES
# endif

#else // #ifdef PARTICLE

// total density equals gas density if there is no particle
#  define _TOTAL_DENS         ( _DENS )

#endif // #ifdef PARTICLE ... else ...



// number of fluid ghost zones for the fluid solver
#if   ( MODEL == HYDRO )   // hydro

#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
#    if   ( LR_SCHEME == PLM )
#     define LR_GHOST_SIZE          1
#    elif ( LR_SCHEME == PPM )
#     define LR_GHOST_SIZE          2
#    else
#     error : ERROR : unsupported LR_SCHEME !!
#    endif
#  endif // MHM/MHM_RP/CTU

#  if   ( FLU_SCHEME == RTVD )
#     define FLU_GHOST_SIZE         3
#  elif ( FLU_SCHEME == MHM )
#     define FLU_GHOST_SIZE         ( 1 + LR_GHOST_SIZE )
#  elif ( FLU_SCHEME == MHM_RP )
#     define FLU_GHOST_SIZE         ( 2 + LR_GHOST_SIZE )
#  elif ( FLU_SCHEME == CTU )
#    ifdef MHD
#     define FLU_GHOST_SIZE         ( 2 + LR_GHOST_SIZE )
#    else
#     define FLU_GHOST_SIZE         ( 1 + LR_GHOST_SIZE )
#    endif // MHD
#  endif // FLU_SCHEME


#elif ( MODEL == ELBDM )   // ELBDM

#  if ( WAVE_SCHEME == WAVE_FD )
#     ifdef LAPLACIAN_4TH
#        define FLU_GHOST_SIZE         6
#     else
//     hybrid scheme requires FLU_GHOST_SIZE >= HYB_GHOST_SIZE (6)
#      if ( ELBDM_SCHEME == ELBDM_HYBRID )
#        define FLU_GHOST_SIZE         6
#      else
#        define FLU_GHOST_SIZE         3
#      endif
#     endif // LAPLACIAN_4TH
#  elif ( WAVE_SCHEME == WAVE_GRAMFE )
// the accuracy of the local spectral method increases with larger FLU_GHOST_SIZE.
// a minimum of FLU_GHOST_SIZE 6 has been found to be stable with the filter options alpha = 100 and beta = 32 * log(10)
// larger ghost zones should increase stability and accuracy and allow for larger timesteps, but have not extensively tested
// for smaller ghost zones, GRAMFE_ORDER should be decreased to values between 6 and 12 and the filter parameters should be adapted
#        define FLU_GHOST_SIZE         8
#  else // WAVE_SCHEME
#     error : ERROR : unsupported WAVE_SCHEME !!
#  endif // WAVE_SCHEME

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL

// define fluid ghost boundary size for hybrid scheme
// --> it must be smaller than or equal to FLU_GHOST_SIZE because the same fluid arrays are used for both the wave and fluid solvers
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
#        define HYB_GHOST_SIZE         6
#  endif

# if ( WAVE_SCHEME == WAVE_GRAMFE  ||  SUPPORT_SPECTRAL_INT )
//  number of evaluation points of Gram polynomials for computing FC(SVD) continuation
#   define GRAMFE_GAMMA  150
//  number of Fourier modes used in the FC(SVD) continuation
//  roughly GRAMFE_G = GRAMFE_GAMMA/2
#   define GRAMFE_G      63
# endif

// set default parameters of gram extension scheme
# if ( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE )
//  number of boundary points used for Gram polynomial space on boundary
#   define GRAMFE_NDELTA 14
//  maximum order of Gram polynomials on boundary
//  for GRAMFE_ORDER < GRAMFE_NDELTA, the boundary information is projected
//  onto a lower-dimensional polynomial space
//  this increases the stability but decreases the accuracy of the algorithm
#   define GRAMFE_ORDER  14

//  a boundary of size GRAMFE_NDELTA can only support polynomials of degree up to GRAMFE_ORDER
#   if ( GRAMFE_ORDER > GRAMFE_NDELTA )
#     error : ERROR : Gram Fourier extension order must not be higher than NDELTA !!
#   endif

//  size of the extension region
//  --> total size of extended region = GRAMFE_FLU_NXT = FLU_NXT + GRAMFE_ND
#   if   ( GRAMFE_SCHEME == GRAMFE_FFT )
//  default values in order for GRAMFE_FLU_NXT to have small prime factorisations
//  --> this is important for the FFT to be fast
#    if   ( PATCH_SIZE == 8 )
#     define GRAMFE_ND        32 // GRAMFE_FLU_NXT = 2^6
#    elif ( PATCH_SIZE == 16 )
#     define GRAMFE_ND        24 // GRAMFE_FLU_NXT = 2^3 * 3^2
#    elif ( PATCH_SIZE == 32 )
#     define GRAMFE_ND        28 // GRAMFE_FLU_NXT = 2^2 * 3^3
#    elif ( PATCH_SIZE == 64 )
#     define GRAMFE_ND        24 // GRAMFE_FLU_NXT = 2^3 * 3 * 7
#    elif ( PATCH_SIZE == 128 )
#     define GRAMFE_ND        28 // GRAMFE_FLU_NXT = 2^2 * 3 * 5^2
#    else
#     error : ERROR : Unsupported PATCH_SIZE for GRAMFE_FFT!!
#    endif // PATCH_SIZE

#   elif ( GRAMFE_SCHEME == GRAMFE_MATMUL )
//  for GRAMFE_MATMUL extension size is irrelevant since matrix multiplication works for all sizes
#     define GRAMFE_ND        32 // GRAMFE_FLU_NXT = 2^6

#     if ( PATCH_SIZE != 8  &&  PATCH_SIZE != 16 )
#       error : ERROR : Unsupported PATCH_SIZE for GRAMFE_MATMUL (only support 8 and 16) !! Consider switching to GRAMFE_FFT.
#     endif

#   else
#     error : ERROR : Unsupported GRAMFE_SCHEME!!
#   endif

//  total size of extended region
#   define GRAMFE_FLU_NXT     ( FLU_NXT + GRAMFE_ND )

# endif // # if ( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE )


// self-gravity constants
#ifdef GRAVITY

// number of input and output variables in the gravity solver
#  if   ( MODEL == HYDRO )
#        define GRA_NIN             ( NCOMP_FLUID )

// for ELBDM, we do not need to transfer the density component
// --> this remains valid for hybrid solver that also has 2 components
#  elif ( MODEL == ELBDM )
#        define GRA_NIN             ( NCOMP_FLUID - 1 )

#  else
#     error Error : unsupported MODEL (please edit GRA_NIN in the new MODEL) !!
#  endif // MODEL


// number of potential ghost zones for evaluating potential (maximum=5) ~ Poisson solver
#        define POT_GHOST_SIZE      5


// number of potential ghost zones for advancing fluid by gravity ~ Gravity solver
#  if   ( MODEL == HYDRO )
#     ifdef STORE_POT_GHOST
#        define GRA_GHOST_SIZE      2
#     else
#        define GRA_GHOST_SIZE      1
//#      define GRA_GHOST_SIZE      2
#     endif

#  elif ( MODEL == ELBDM )
#     ifdef STORE_POT_GHOST
#        define GRA_GHOST_SIZE      2
#     else
#        define GRA_GHOST_SIZE      0
#     endif

#  elif ( MODEL == PAR_ONLY )
#     ifdef STORE_POT_GHOST
#        define GRA_GHOST_SIZE      2
#     else
#        define GRA_GHOST_SIZE      0
#     endif

#  else
#     error : ERROR : unsupported MODEL !!
#  endif // MODEL


// number of potential ghost zones for correcting the half-step velocity if UNSPLIT_GRAVITY is on
// _F/_G: fluid/gravity solvers
#  ifdef UNSPLIT_GRAVITY
#     if   ( MODEL == HYDRO )
#       ifdef MHD
#        define USG_GHOST_SIZE_F    2
#        define USG_GHOST_SIZE_G    1
#       else
#        define USG_GHOST_SIZE_F    1
#        define USG_GHOST_SIZE_G    1
#       endif
#     elif ( MODEL == ELBDM )
#        define USG_GHOST_SIZE_F    0
#        define USG_GHOST_SIZE_G    0
#     else
#        error : ERROR : unsupported MODEL !!
#     endif // MODEL
#  endif // #ifdef UNSPLIT_GRAVITY


// number of density ghost zones for storing the temporary particle mass density in rho_ext[]
#  ifdef PARTICLE
#        define RHOEXT_GHOST_SIZE   2
#  endif


// number of density ghost zones for the Poisson solver
#        define RHO_GHOST_SIZE      ( POT_GHOST_SIZE-1 )

#endif // #ifdef GRAVITY


// number of ghost zones for the source-term solver
// --> fixed to zero for now since ghost zones in source terms are not supported yet
#        define SRC_GHOST_SIZE      0


// number of ghost zones for computing derived fields
#        define DER_GHOST_SIZE      1


// number of ghost zones for feedback
// --> can be changed manually
// --> set to 0 if applicable to improve performance
#ifdef FEEDBACK
#        define FB_GHOST_SIZE       3
#endif



// patch size (number of cells of a single patch in the x/y/z directions)
#define PS1             ( 1*PATCH_SIZE )
#define PS2             ( 2*PATCH_SIZE )
#define PS2P1           ( PS2 + 1 )
#define PS1M1           ( PS1 - 1 )
#define PS1P1           ( PS1 + 1 )


// size of GPU arrays (in one dimension)
//###REVISE: support interpolation schemes requiring 2 ghost cells on each side for POT_NXT
#  define FLU_NXT       ( PS2 + 2*FLU_GHOST_SIZE )                // use patch group as the unit
#  define FLU_NXT_P1    ( FLU_NXT + 1 )
#ifdef GRAVITY
#  define POT_NXT       ( PS1/2 + 2*( (POT_GHOST_SIZE+3)/2 ) )    // assuming interpolation ghost zone == 1
#  define RHO_NXT       ( PS1 + 2*RHO_GHOST_SIZE )                // POT/RHO/GRA_NXT use patch as the unit
#  define GRA_NXT       ( PS1 + 2*GRA_GHOST_SIZE )
#  ifdef UNSPLIT_GRAVITY
#  define USG_NXT_F     ( PS2 + 2*USG_GHOST_SIZE_F )              // we use patch group as unit for the fluid   solver
#  define USG_NXT_G     ( PS1 + 2*USG_GHOST_SIZE_G )              // we use patch       as unit for the gravity solver
#  else
#  define USG_NXT_F     ( 1 )                                     // still define USG_NXT_F/G since many function prototypes
#  define USG_NXT_G     ( 1 )                                     // require it
#  endif
#  ifdef PARTICLE
#  define RHOEXT_NXT    ( PS1 + 2*RHOEXT_GHOST_SIZE )             // array rho_ext of each patch
#  endif
#else
#  define GRA_NXT       ( 1 )                                     // still define GRA_NXT   ...
#  define USG_NXT_F     ( 1 )                                     // still define USG_NXT_F ...
#endif // GRAVITY
#  define SRC_NXT       ( PS1 + 2*SRC_GHOST_SIZE )                // use patch as the unit
#  define SRC_NXT_P1    ( SRC_NXT + 1 )
#  define DER_NXT       ( PS1 + 2*DER_GHOST_SIZE )                // use patch as the unit
#ifdef FEEDBACK
#  define FB_NXT        ( PS2 + 2*FB_GHOST_SIZE )                 // use patch group as the unit
#endif
#if ( ELBDM_SCHEME == ELBDM_HYBRID )
#  define HYB_NXT       ( PS2 + 2*HYB_GHOST_SIZE )
#else
#  define HYB_NXT       ( 1 )
#endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )


// size of auxiliary arrays and EoS tables
#if ( MODEL == HYDRO )
#  define EOS_NAUX_MAX           20    // EoS_AuxArray_Flt/Int[]
#  define EOS_NTABLE_MAX         20    // *_EoS_Table[]
#else
#  define EOS_NAUX_MAX           0
#  define EOS_NTABLE_MAX         0
#endif

#ifdef GRAVITY
#  define EXT_POT_NAUX_MAX       20    // ExtPot_AuxArray[]
#  define EXT_ACC_NAUX_MAX       20    // ExtAcc_AuxArray[]
#  define EXT_POT_NGENE_MAX       6    // h/d_ExtPotGenePtr
#endif

#if ( MODEL == HYDRO )
#  define SRC_NAUX_DLEP          5     // SrcTerms.Dlep_AuxArray_Flt/Int[]
#  define SRC_DLEP_PROF_NVAR     6     // SrcTerms.Dlep_Profile_DataDevPtr[]/RadiusDevPtr[]
#  define SRC_DLEP_PROF_NBINMAX  4000
#else
#  define SRC_NAUX_DLEP          0
#endif
#  define SRC_NAUX_USER          10    // SrcTerms.User_AuxArray_Flt/Int[]


// bitwise reproducibility in flux and electric field fix-up operations
# ifdef BITWISE_REPRODUCIBILITY
#  define BIT_REP_FLUX
# endif

// enable BIT_REP_ELECTRIC by default even when BITWISE_REPRODUCIBILITY is off
// --> ensures that the B field on the common interface between two nearby patches are fully
//     consistent with each other (even the round-off errors are the same)
// --> reduces the div(B) errors significantly
# ifdef MHD
//# ifdef BITWISE_REPRODUCIBILITY
#  define BIT_REP_ELECTRIC
//# endif
# endif // MHD


// only apply iterations to broken cells in Interpolate_Iterate()
#define INTERP_MASK

// used by INTERP_MASK for now but can be applied to other places in the future
#define MASKED                   true
#define UNMASKED                 false


// in FB_AdvanceDt(), store the updated fluid data in a separate array to avoid data racing among different patch groups
#if ( defined FEEDBACK  &&  FB_GHOST_SIZE > 0 )
#  define FB_SEP_FLUOUT
#endif


// precision for FFT in GRAMFE_FFT and matrix multiplication in GRAMFE_MATMUL
// --> enable double precision for GRAMFE_FFT by default since it is less stable compared to GRAMFE_MATMUL
#if ( GRAMFE_SCHEME == GRAMFE_FFT )
#   define GRAMFE_FFT_FLOAT8
#endif

#if ( ( GRAMFE_SCHEME == GRAMFE_MATMUL ) && defined( FLOAT8 ) )
#   define GRAMFE_MATMUL_FLOAT8
#endif


// extreme values
#ifndef __INT_MAX__
#  define __INT_MAX__            2147483647
#endif

#ifndef __LONG_MAX__
#  define __LONG_MAX__           9223372036854775807L
#endif

#ifndef __UINT_MAX__
#  define __UINT_MAX__           ( __INT_MAX__*2U + 1U )
#endif

#ifndef __ULONG_MAX__
#  define __ULONG_MAX__          18446744073709551615UL     // 2^64-1
#endif

#ifndef __FLT_MAX__
#  define __FLT_MAX__            3.40282347e+38F
#endif

#ifndef __FLT_MIN__
#  define __FLT_MIN__            1.17549435e-38F
#endif

#ifndef __DBL_MAX__
#  define __DBL_MAX__            1.79769313e+308
#endif

#ifndef __DBL_MIN__
#  define __DBL_MIN__            2.22507386e-308
#endif

#ifndef __FLT_EPSILON__
#  define __FLT_EPSILON__        1.19209290e-07F
#endif

#ifndef __DBL_EPSILON__
#  define __DBL_EPSILON__        2.2204460492503131e-16
#endif

// extreme value used for various purposes (e.g., floor value for passive scalars)
#ifdef FLOAT8
#  define TINY_NUMBER            __DBL_MIN__
#  define HUGE_NUMBER            __DBL_MAX__
#else
#  define TINY_NUMBER            __FLT_MIN__
#  define HUGE_NUMBER            __FLT_MAX__
#endif


// maximum allowed error for various purposes (e.g., exact Riemann solver, MHD routines, Mis_CompareRealValue())
#define MAX_ERROR_DBL            1.0e-14
#define MAX_ERROR_FLT            1.0e-06f

#ifdef FLOAT8
#  define MACHINE_EPSILON        __DBL_EPSILON__
#  define MAX_ERROR              MAX_ERROR_DBL
#else
#  define MACHINE_EPSILON        __FLT_EPSILON__
#  define MAX_ERROR              MAX_ERROR_FLT
#endif


// sibling index offset for the non-periodic B.C.
#define SIB_OFFSET_NONPERIODIC   ( -100 )


// son index offset for LOAD_BALANCE
#ifdef LOAD_BALANCE
#  define SON_OFFSET_LB          ( -1000 )
#endif


// flag used in "Buf_RecordBoundaryFlag" and "Flag_Buffer" (must be negative)
#ifndef SERIAL
#  define BUFFER_IS_FLAGGED      ( -999 )
#endif


// marker indicating that the array "pot_ext" has NOT been properly set
#if ( defined GRAVITY  &&  defined STORE_POT_GHOST )
#  define POT_EXT_NEED_INIT      __FLT_MAX__
#endif


// marker indicating that the array "rho_ext" has NOT been properly set
#ifdef PARTICLE
#  define RHO_EXT_NEED_INIT      __FLT_MAX__
#endif


// markers for inactive particles
#ifdef PARTICLE
#  define PAR_INACTIVE_OUTSIDE   ( -1.0 )
#  define PAR_INACTIVE_MPI       ( -2.0 )
#endif


// OpenMP scheduling for particle routines
#if ( defined PARTICLE  &&  defined OPENMP )
#  define PAR_OMP_SCHED          dynamic
#  define PAR_OMP_SCHED_CHUNK    1
#endif


// NULL values
#ifndef NULL
#  define NULL             0
#endif

#ifndef NULL_INT
#  define NULL_INT         __INT_MAX__
#endif

#ifndef NULL_REAL
#  define NULL_REAL        __FLT_MAX__
#endif

#ifndef NULL_BOOL
#  define NULL_BOOL        false
#endif


// GAMER status
// --> if we ever want to swap the following values, must check all MPI functions using MPI_BAND or MPI_BOR
#define GAMER_SUCCESS      1
#define GAMER_FAILED       0


// timer switch
#define TIMER_ON           1
#define TIMER_OFF          0


// symbolic constant for Aux_Error()
#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__


// miscellaneous
#define TOP_LEVEL          ( NLEVEL - 1 )


// maximum length for strings
#define MAX_STRING         512


// MPI floating-point data type
#ifdef FLOAT8
#  define MPI_GAMER_REAL MPI_DOUBLE
#else
#  define MPI_GAMER_REAL MPI_FLOAT
#endif

#ifdef FLOAT8_PAR
#  define MPI_GAMER_REAL_PAR MPI_DOUBLE
#else
#  define MPI_GAMER_REAL_PAR MPI_FLOAT
#endif

#ifdef INT8_PAR
#  define MPI_GAMER_LONG_PAR MPI_LONG
#else
#  define MPI_GAMER_LONG_PAR MPI_INT
#endif



// ############
// ## Macros ##
// ############

// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define   LOG( a )         log( a )
#  define   EXP( a )         exp( a )
#  define  ATAN( a )        atan( a )
#  define FLOOR( a )       floor( a )
#  define ROUND( a )       round( a )
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
#  define   EXP( a )         expf( a )
#  define  ATAN( a )        atanf( a )
#  define FLOOR( a )       floorf( a )
#  define ROUND( a )       roundf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define   POW( a, b )      powf( a, b )
#  define  FMOD( a, b )     fmodf( a, b )
#  define ATAN2( a, b )    atan2f( a, b )
#endif


// sign function
#define SIGN( a )       (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// safe ATAN2 that does not return nan when a = b = 0
#define SATAN2( a, b )   (  ( (a) == (real)0.0  &&  (b) == (real)0.0 ) ? (real)0.0 : ATAN2( (a), (b) )  )


// power functions
#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )
#define POW4( a )       ( (a)*(a)*(a)*(a) )


// 3D to 1D array indices transformation
#define IDX321( i, j, k, Ni, Nj )   (  ( (k)*(Nj) + (j) )*(Ni) + (i)  )


// 3D to 1D array indices transformation for patch->magnetic[]
#ifdef MHD
#define IDX321_BX( i, j, k, Ni, Nj    )   (  ( (k)*((Nj)  ) + (j) )*((Ni)+1) + (i)  )
#define IDX321_BY( i, j, k, Ni, Nj    )   (  ( (k)*((Nj)+1) + (j) )*((Ni)  ) + (i)  )
#define IDX321_BZ( i, j, k, Ni, Nj    )   (  ( (k)*((Nj)  ) + (j) )*((Ni)  ) + (i)  )
#define IDX321_B(  i, j, k, Ni, Nj, d )   (  ( (d) == 0 ) ? IDX321_BX( i, j, k, Ni, Nj ) : \
                                             ( (d) == 1 ) ? IDX321_BY( i, j, k, Ni, Nj ) : \
                                             ( (d) == 2 ) ? IDX321_BZ( i, j, k, Ni, Nj ) : NULL_INT  )
#endif


// helper macros for printing symbolic constants in macros
// ref: https://stackoverflow.com/questions/3419332/c-preprocessor-stringify-the-result-of-a-macro
#  define QUOTE( str )              #str
#  define EXPAND_AND_QUOTE( str )   QUOTE( str )


// convenient macros for defining and declaring global variables
// ==> predefine DEFINE_GLOBAL in the file actually **defines** these global variables
// ==> there should be one and only one file that defines DEFINE_GLOBAL

// SET_GLOBAL will invoke either SET_GLOBAL_INIT or SET_GLOBAL_NOINIT depending on the number of arguments
// ==> http://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments
#define GET_MACRO( _1, _2, TARGET_MACRO, ... )  TARGET_MACRO
#define SET_GLOBAL( ... )   GET_MACRO( __VA_ARGS__, SET_GLOBAL_INIT, SET_GLOBAL_NOINIT ) ( __VA_ARGS__ )

// SET_GLOBAL_INIT/NOINIT are for global variables with/without initialization
#ifdef DEFINE_GLOBAL
# define SET_GLOBAL_INIT( declaration, init_value )   declaration = init_value
# define SET_GLOBAL_NOINIT( declaration )             declaration
#else
# define SET_GLOBAL_INIT( declaration, init_value )   extern declaration
# define SET_GLOBAL_NOINIT( declaration )             extern declaration
#endif


// macro converting an array index (e.g., DENS) to bitwise index (e.g., _DENS=(1L<<DENS))
#define BIDX( idx )     ( 1L << (idx) )


// helper macro for printing warning messages when resetting parameters
#  define FORMAT_INT       %- 21d
#  define FORMAT_LONG      %- 21ld
#  define FORMAT_UINT      %- 21u
#  define FORMAT_ULONG     %- 21lu
#  define FORMAT_BOOL      %- 21d
#  define FORMAT_REAL      %- 21.14e
#  define PRINT_RESET_PARA( name, format, reason )                                                       \
   {                                                                                                     \
      if ( MPI_Rank == 0 )                                                                               \
         Aux_Message( stderr, "WARNING : parameter [%-30s] is reset to [" EXPAND_AND_QUOTE(format) "] "  \
                              "%s\n", #name, name, reason );                                             \
   }


// ################################
// ## Remove useless definitions ##
// ################################
#if ( MODEL == HYDRO )
#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#  undef LR_SCHEME
#  endif

#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#  undef RSOLVER
#  endif
#endif

#if ( MODEL == PAR_ONLY )
#  undef UNSPLIT_GRAVITY
#endif

// currently we always set GPU_ARCH == NONE when GPU is off
#ifndef GPU
#  undef  GPU_ARCH
#  define GPU_ARCH NONE
#endif

// currently we assume that particle acceleration is solely due to gravity
// --> if we ever want to enable STORE_PAR_ACC without GRAVITY, remember to check all STORE_PAR_ACC in this header
#ifndef GRAVITY
#  undef STORE_PAR_ACC
#endif



#endif  // #ifndef __MACRO_H__
