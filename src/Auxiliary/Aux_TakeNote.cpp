#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif
#include <sched.h>
#ifdef __APPLE__
#include <cpuid.h>
#endif
#include "time.h"

static int get_cpuid();




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_TakeNote
// Description :  Record simulation parameters and the content in the file "Input__Note" to the
//                note file "Record__Note"
//-------------------------------------------------------------------------------------------------------
void Aux_TakeNote()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_TakeNote ...\n" );


   const char FileName[] = "Record__Note";
   FILE *Note;

   if ( MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

//    copy the content in the file "Input__Note"
      Note = fopen( FileName, "a" );
      fprintf( Note, "\n\n\nSimulation Notes\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fclose( Note );

      system( "cat ./Input__Note >> Record__Note" );

      Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the simulation options in the Makefile (numerical schemes)
      fprintf( Note, "Makefile Options (numerical schemes)\n" );
      fprintf( Note, "***********************************************************************************\n" );

//    a. options for all models
#     if   ( MODEL == HYDRO )
      fprintf( Note, "MODEL                           HYDRO\n" );
#     elif ( MODEL == ELBDM )
      fprintf( Note, "MODEL                           ELBDM\n" );
#     elif ( MODEL == PAR_ONLY )
      fprintf( Note, "MODEL                           PAR_ONLY\n" );
#     else
      fprintf( Note, "MODEL                           UNKNOWN\n" );
#     endif // MODEL

#     ifdef GRAVITY
      fprintf( Note, "GRAVITY                         ON\n" );
#     else
      fprintf( Note, "GRAVITY                         OFF\n" );
#     endif

#     ifdef GRAVITY
#     if   ( POT_SCHEME == SOR )
      fprintf( Note, "POT_SCHEME                      SOR\n" );
#     elif ( POT_SCHEME == MG )
      fprintf( Note, "POT_SCHEME                      MG\n" );
#     elif ( POT_SCHEME == NONE )
      fprintf( Note, "POT_SCHEME                      NONE\n" );
#     else
      fprintf( Note, "POT_SCHEME                      UNKNOWN\n" );
#     endif

#     ifdef STORE_POT_GHOST
      fprintf( Note, "STORE_POT_GHOST                 ON\n" );
#     else
      fprintf( Note, "STORE_POT_GHOST                 OFF\n" );
#     endif

#     ifdef UNSPLIT_GRAVITY
      fprintf( Note, "UNSPLIT_GRAVITY                 ON\n" );
#     else
      fprintf( Note, "UNSPLIT_GRAVITY                 OFF\n" );
#     endif
#     endif // #ifdef GRAVITY

#     ifdef COMOVING
      fprintf( Note, "COMOVING                        ON\n" );
#     else
      fprintf( Note, "COMOVING                        OFF\n" );
#     endif

#     ifdef PARTICLE
      fprintf( Note, "PARTICLE                        ON\n" );
#     else
      fprintf( Note, "PARTICLE                        OFF\n" );
#     endif

#     ifdef SUPPORT_GRACKLE
      fprintf( Note, "SUPPORT_GRACKLE                 ON\n" );
#     else
      fprintf( Note, "SUPPORT_GRACKLE                 OFF\n" );
#     endif

//    b. options in HYDRO
#     if   ( MODEL == HYDRO )

#     if   ( FLU_SCHEME == RTVD )
      fprintf( Note, "FLU_SCHEME                      RTVD\n" );
#     elif ( FLU_SCHEME == MHM )
      fprintf( Note, "FLU_SCHEME                      MHM\n" );
#     elif ( FLU_SCHEME == MHM_RP )
      fprintf( Note, "FLU_SCHEME                      MHM with Riemann prediction\n" );
#     elif ( FLU_SCHEME == CTU )
      fprintf( Note, "FLU_SCHEME                      CTU\n" );
#     elif ( FLU_SCHEME == NONE )
      fprintf( Note, "FLU_SCHEME                      NONE\n" );
#     else
      fprintf( Note, "FLU_SCHEME                      UNKNOWN\n" );
#     endif

#     if   ( LR_SCHEME == PLM )
      fprintf( Note, "LR_SCHEME                       PLM\n" );
#     elif ( LR_SCHEME == PPM )
      fprintf( Note, "LR_SCHEME                       PPM\n" );
#     elif ( LR_SCHEME == NONE )
      fprintf( Note, "LR_SCHEME                       NONE\n" );
#     else
      fprintf( Note, "LR_SCHEME                       UNKNOWN\n" );
#     endif

#     if   ( RSOLVER == EXACT )
      fprintf( Note, "RSOLVER                         EXACT\n" );
#     elif ( RSOLVER == ROE )
      fprintf( Note, "RSOLVER                         ROE\n" );
#     elif ( RSOLVER == HLLE )
      fprintf( Note, "RSOLVER                         HLLE\n" );
#     elif ( RSOLVER == HLLC )
      fprintf( Note, "RSOLVER                         HLLC\n" );
#     elif ( RSOLVER == HLLD )
      fprintf( Note, "RSOLVER                         HLLD\n" );
#     elif ( RSOLVER == NONE )
      fprintf( Note, "RSOLVER                         NONE\n" );
#     else
      fprintf( Note, "RSOLVER                         UNKNOWN\n" );
#     endif

#     if   ( DUAL_ENERGY == DE_ENPY )
      fprintf( Note, "DUAL_ENERGY                     DE_ENPY\n" );
#     elif ( DUAL_ENERGY == DE_EINT )
      fprintf( Note, "DUAL_ENERGY                     DE_EINT\n" );
#     elif ( DUAL_ENERGY == NONE )
      fprintf( Note, "DUAL_ENERGY                     NONE\n" );
#     else
      fprintf( Note, "DUAL_ENERGY                     UNKNOWN\n" );
#     endif

#     ifdef MHD
      fprintf( Note, "MHD                             ON\n" );
#     else
      fprintf( Note, "MHD                             OFF\n" );
#     endif

#     ifdef SRHD
      fprintf( Note, "SRHD                            ON\n" );
#     else
      fprintf( Note, "SRHD                            OFF\n" );
#     endif

#     ifdef COSMIC_RAY
      fprintf( Note, "COSMIC_RAY                      ON\n" );
#     ifdef CR_DIFFUSION
      fprintf( Note, "CR_DIFFUSION                    ON\n" );
#     else
      fprintf( Note, "CR_DIFFUSION                    OFF\n" );
#     endif
#     else // #ifdef COSMIC_RAY
      fprintf( Note, "COSMIC_RAY                      OFF\n" );
#     endif // #ifdef COSMIC_RAY ... else ...

#     if   ( EOS == EOS_GAMMA )
      fprintf( Note, "EOS                             EOS_GAMMA\n" );
#     elif ( EOS == EOS_ISOTHERMAL )
      fprintf( Note, "EOS                             EOS_ISOTHERMAL\n" );
#     elif ( EOS == EOS_NUCLEAR )
      fprintf( Note, "EOS                             EOS_NUCLEAR\n" );
#     elif ( EOS == EOS_TAUBMATHEWS )
      fprintf( Note, "EOS                             EOS_TAUBMATHEWS\n" );
#     elif ( EOS == EOS_TABULAR )
      fprintf( Note, "EOS                             EOS_TABULAR\n" );
#     elif ( EOS == EOS_COSMIC_RAY )
      fprintf( Note, "EOS                             EOS_COSMIC_RAY\n" );
#     elif ( EOS == EOS_USER )
      fprintf( Note, "EOS                             EOS_USER\n" );
#     else
      fprintf( Note, "EOS                             UNKNOWN\n" );
#     endif

#     ifdef BAROTROPIC_EOS
      fprintf( Note, "BAROTROPIC_EOS                  ON\n" );
#     else
      fprintf( Note, "BAROTROPIC_EOS                  OFF\n" );
#     endif

//    c. options in ELBDM
#     elif ( MODEL == ELBDM )

#     ifdef CONSERVE_MASS
      fprintf( Note, "CONSERVE_MASS                   ON\n" );
#     else
      fprintf( Note, "CONSERVE_MASS                   OFF\n" );
#     endif

#     ifdef LAPLACIAN_4TH
      fprintf( Note, "LAPLACIAN_4TH                   ON\n" );
#     else
      fprintf( Note, "LAPLACIAN_4TH                   OFF\n" );
#     endif

#     ifdef QUARTIC_SELF_INTERACTION
      fprintf( Note, "QUARTIC_SELF_INTERACTION        ON\n" );
#     else
      fprintf( Note, "QUARTIC_SELF_INTERACTION        OFF\n" );
#     endif

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

//    d. options in PARTICLE
#     ifdef PARTICLE

#     ifdef MASSIVE_PARTICLES
      fprintf( Note, "MASSIVE_PARTICLES               ON\n" );
#     else
      fprintf( Note, "MASSIVE_PARTICLES               OFF\n" );
#     endif

#     ifdef TRACER
      fprintf( Note, "TRACER                          ON\n" );
#     else
      fprintf( Note, "TRACER                          OFF\n" );
#     endif

#     ifdef STORE_PAR_ACC
      fprintf( Note, "STORE_PAR_ACC                   ON\n" );
#     else
      fprintf( Note, "STORE_PAR_ACC                   OFF\n" );
#     endif

#     ifdef STAR_FORMATION
      fprintf( Note, "STAR_FORMATION                  ON\n" );
#     else
      fprintf( Note, "STAR_FORMATION                  OFF\n" );
#     endif

#     ifdef FEEDBACK
      fprintf( Note, "FEEDBACK                        ON\n" );
#     else
      fprintf( Note, "FEEDBACK                        OFF\n" );
#     endif

#     endif // #ifdef PARTICLE

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the simulation options in the Makefile (optimization and compilation)
      fprintf( Note, "Makefile Options (optimization and compilation)\n" );
      fprintf( Note, "***********************************************************************************\n" );

#     ifdef GPU
      fprintf( Note, "GPU                             ON\n" );
#     else
      fprintf( Note, "GPU                             OFF\n" );
#     endif

#     ifdef GAMER_DEBUG
      fprintf( Note, "GAMER_DEBUG                     ON\n" );
#     else
      fprintf( Note, "GAMER_DEBUG                     OFF\n" );
#     endif

#     ifdef BITWISE_REPRODUCIBILITY
      fprintf( Note, "BITWISE_REPRODUCIBILITY         ON\n" );
#     else
      fprintf( Note, "BITWISE_REPRODUCIBILITY         OFF\n" );
#     endif

#     ifdef TIMING
      fprintf( Note, "TIMING                          ON\n" );
#     else
      fprintf( Note, "TIMING                          OFF\n" );
#     endif

#     ifdef TIMING_SOLVER
      fprintf( Note, "TIMING_SOLVER                   ON\n" );
#     else
      fprintf( Note, "TIMING_SOLVER                   OFF\n" );
#     endif

#     ifdef FLOAT8
      fprintf( Note, "FLOAT8                          ON\n" );
#     else
      fprintf( Note, "FLOAT8                          OFF\n" );
#     endif

#     ifdef FLOAT8_PAR
      fprintf( Note, "FLOAT8_PAR                      ON\n" );
#     else
      fprintf( Note, "FLOAT8_PAR                      OFF\n" );
#     endif

#     ifdef SERIAL
      fprintf( Note, "SERIAL                          ON\n" );
#     else
      fprintf( Note, "SERIAL                          OFF\n" );
#     endif

#     ifdef LOAD_BALANCE
#     if   ( LOAD_BALANCE == HILBERT )
      fprintf( Note, "LOAD_BALANCE                    HILBERT\n" );
#     else
      fprintf( Note, "LOAD_BALANCE                    UNKNOWN\n" );
#     endif
#     else // #ifdef LOAD_BALANCE
      fprintf( Note, "LOAD_BALANCE                    OFF\n" );
#     endif // #ifdef LOAD_BALANCE ... else ...

#     ifdef OVERLAP_MPI
      fprintf( Note, "OVERLAP_MPI                     ON\n" );
#     else
      fprintf( Note, "OVERLAP_MPI                     OFF\n" );
#     endif

#     ifdef OPENMP
      fprintf( Note, "OPENMP                          ON\n" );
#     else
      fprintf( Note, "OPENMP                          OFF\n" );
#     endif

#     ifdef GPU
#     if   ( GPU_ARCH == FERMI )
      fprintf( Note, "GPU_ARCH                        FERMI\n" );
#     elif ( GPU_ARCH == KEPLER )
      fprintf( Note, "GPU_ARCH                        KEPLER\n" );
#     elif ( GPU_ARCH == MAXWELL )
      fprintf( Note, "GPU_ARCH                        MAXWELL\n" );
#     elif ( GPU_ARCH == PASCAL )
      fprintf( Note, "GPU_ARCH                        PASCAL\n" );
#     elif ( GPU_ARCH == VOLTA )
      fprintf( Note, "GPU_ARCH                        VOLTA\n" );
#     elif ( GPU_ARCH == TURING )
      fprintf( Note, "GPU_ARCH                        TURING\n" );
#     elif ( GPU_ARCH == AMPERE )
      fprintf( Note, "GPU_ARCH                        AMPERE\n" );
#     else
      fprintf( Note, "GPU_ARCH                        UNKNOWN\n" );
#     endif
      fprintf( Note, "GPU_COMPUTE_CAPABILITY          %d\n", GPU_COMPUTE_CAPABILITY );
#     endif // #ifdef GPU

#     ifdef LAOHU
      fprintf( Note, "LAOHU                           ON\n" );
#     else
      fprintf( Note, "LAOHU                           OFF\n" );
#     endif

#     ifdef SUPPORT_HDF5
      fprintf( Note, "SUPPORT_HDF5                    ON\n" );
#     else
      fprintf( Note, "SUPPORT_HDF5                    OFF\n" );
#     endif

#     ifdef SUPPORT_GSL
      fprintf( Note, "SUPPORT_GSL                     ON\n" );
#     else
      fprintf( Note, "SUPPORT_GSL                     OFF\n" );
#     endif

#     ifdef SUPPORT_FFTW
      fprintf( Note, "SUPPORT_FFTW                    ON\n" );
#     else
      fprintf( Note, "SUPPORT_FFTW                    OFF\n" );
#     endif

#     ifdef SUPPORT_LIBYT
      fprintf( Note, "SUPPORT_LIBYT                   ON\n" );
#     else
      fprintf( Note, "SUPPORT_FFTW                    OFF\n" );
#     endif

#     ifdef SUPPORT_LIBYT
      fprintf( Note, "SUPPORT_LIBYT                   ON\n" );

#     ifdef LIBYT_USE_PATCH_GROUP
      fprintf( Note, "LIBYT_USE_PATCH_GROUP           ON\n" );
#     else
      fprintf( Note, "LIBYT_USE_PATCH_GROUP           OFF\n" );
#     endif

#     ifdef LIBYT_INTERACTIVE
      fprintf( Note, "LIBYT_INTERACTIVE               ON\n" );
#     else
      fprintf( Note, "LIBYT_INTERACTIVE               OFF\n" );
#     endif

#     ifdef LIBYT_RELOAD
      fprintf( Note, "LIBYT_RELOAD                    ON\n" );
#     else
      fprintf( Note, "LIBYT_RELOAD                    OFF\n" );
#     endif

#     ifdef LIBYT_JUPYTER
      fprintf( Note, "LIBYT_JUPYTER                   ON\n" );
#     else
      fprintf( Note, "LIBYT_JUPYTER                   OFF\n" );
#     endif

#     else  // #ifdef SUPPORT_LIBYT
      fprintf( Note, "SUPPORT_LIBYT                   OFF\n" );
#     endif // #ifdef SUPPORT_LIBYT ... else ...

#     if   ( RANDOM_NUMBER == RNG_GNU_EXT )
      fprintf( Note, "RANDOM_NUMBER                   RNG_GNU_EXT\n" );
#     elif ( RANDOM_NUMBER == RNG_CPP11 )
      fprintf( Note, "RANDOM_NUMBER                   RNG_CPP11\n" );
#     else
      fprintf( Note, "RANDOM_NUMBER                   UNKNOWN\n" );
#     endif

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the simulation options in Macro.h, CUFLU.h and CUPOT.h
      fprintf( Note, "Other Options (in Macro.h, CUFLU.h and CUPOT.h)\n" );
      fprintf( Note, "***********************************************************************************\n" );

#     ifdef BIT_REP_FLUX
      fprintf( Note, "BIT_REP_FLUX                    ON\n" );
#     else
      fprintf( Note, "BIT_REP_FLUX                    OFF\n" );
#     endif

#     ifdef MHD
#     ifdef BIT_REP_ELECTRIC
      fprintf( Note, "BIT_REP_ELECTRIC                ON\n" );
#     else
      fprintf( Note, "BIT_REP_ELECTRIC                OFF\n" );
#     endif
#     endif

#     ifdef INTERP_MASK
      fprintf( Note, "INTERP_MASK                     ON\n" );
#     else
      fprintf( Note, "INTERP_MASK                     OFF\n" );
#     endif

#     ifdef FB_SEP_FLUOUT
      fprintf( Note, "FB_SEP_FLUOUT                   ON\n" );
#     else
      fprintf( Note, "FB_SEP_FLUOUT                   OFF\n" );
#     endif

#     if   ( MODEL == HYDRO )
#     ifdef CHECK_UNPHYSICAL_IN_FLUID
      fprintf( Note, "CHECK_UNPHYSICAL_IN_FLUID       ON\n" );
#     else
      fprintf( Note, "CHECK_UNPHYSICAL_IN_FLUID       OFF\n" );
#     endif

#     ifdef CHAR_RECONSTRUCTION
      fprintf( Note, "CHAR_RECONSTRUCTION             ON\n" );
#     else
      fprintf( Note, "CHAR_RECONSTRUCTION             OFF\n" );
#     endif

#     ifdef LR_EINT
      fprintf( Note, "LR_EINT                         ON\n" );
#     else
      fprintf( Note, "LR_EINT                         OFF\n" );
#     endif

#     if   ( CHECK_INTERMEDIATE == EXACT )
      fprintf( Note, "CHECK_INTERMEDIATE              EXACT\n" );
#     elif ( CHECK_INTERMEDIATE == HLLE )
      fprintf( Note, "CHECK_INTERMEDIATE              HLLE\n" );
#     elif ( CHECK_INTERMEDIATE == HLLC )
      fprintf( Note, "CHECK_INTERMEDIATE              HLLC\n" );
#     elif ( CHECK_INTERMEDIATE == HLLD )
      fprintf( Note, "CHECK_INTERMEDIATE              HLLD\n" );
#     elif ( CHECK_INTERMEDIATE == NONE )
      fprintf( Note, "CHECK_INTERMEDIATE              OFF\n" );
#     else
      fprintf( Note, "CHECK_INTERMEDIATE              UNKNOWN\n" );
#     endif

#     if   ( RSOLVER_RESCUE == EXACT )
      fprintf( Note, "RSOLVER_RESCUE                  EXACT\n" );
#     elif ( RSOLVER_RESCUE == HLLE )
      fprintf( Note, "RSOLVER_RESCUE                  HLLE\n" );
#     elif ( RSOLVER_RESCUE == HLLC )
      fprintf( Note, "RSOLVER_RESCUE                  HLLC\n" );
#     elif ( RSOLVER_RESCUE == HLLD )
      fprintf( Note, "RSOLVER_RESCUE                  HLLD\n" );
#     elif ( RSOLVER_RESCUE == NONE )
      fprintf( Note, "RSOLVER_RESCUE                  OFF\n" );
#     else
      fprintf( Note, "RSOLVER_RESCUE                  UNKNOWN\n" );
#     endif

#     ifdef HLL_NO_REF_STATE
      fprintf( Note, "HLL_NO_REF_STATE                ON\n" );
#     else
      fprintf( Note, "HLL_NO_REF_STATE                OFF\n" );
#     endif

#     ifdef HLL_INCLUDE_ALL_WAVES
      fprintf( Note, "HLL_INCLUDE_ALL_WAVES           ON\n" );
#     else
      fprintf( Note, "HLL_INCLUDE_ALL_WAVES           OFF\n" );
#     endif

      fprintf( Note, "HLLC_WAVESPEED                 % d\n",      HLLC_WAVESPEED );
      fprintf( Note, "HLLE_WAVESPEED                 % d\n",      HLLE_WAVESPEED );
#     ifdef MHD
      fprintf( Note, "HLLD_WAVESPEED                 % d\n",      HLLD_WAVESPEED );
#     endif

#     ifdef MHD
#     ifdef EULERY
      fprintf( Note, "EULERY                          ON\n" );
#     else
      fprintf( Note, "EULERY                          OFF\n" );
#     endif
#     endif // #ifdef MHD

#     ifdef MHM_CHECK_PREDICT
      fprintf( Note, "MHM_CHECK_PREDICT               ON\n" );
#     else
      fprintf( Note, "MHM_CHECK_PREDICT               OFF\n" );
#     endif

#     elif ( MODEL == ELBDM )

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

#     if ( defined GRAVITY  &&  POT_SCHEME == SOR  &&  defined GPU )
#     ifdef SOR_RHO_SHARED
      fprintf( Note, "SOR_RHO_SHARED                  ON\n" );
#     else
      fprintf( Note, "SOR_RHO_SHARED                  OFF\n" );
#     endif

#     ifdef SOR_CPOT_SHARED
      fprintf( Note, "SOR_CPOT_SHARED                 ON\n" );
#     else
      fprintf( Note, "SOR_CPOT_SHARED                 OFF\n" );
#     endif

#     ifdef SOR_USE_SHUFFLE
      fprintf( Note, "SOR_USE_SHUFFLE                 ON\n" );
#     else
      fprintf( Note, "SOR_USE_SHUFFLE                 OFF\n" );
#     endif

#     ifdef SOR_USE_PADDING
      fprintf( Note, "SOR_USE_PADDING                 ON\n" );
#     else
      fprintf( Note, "SOR_USE_PADDING                 OFF\n" );
#     endif

      fprintf( Note, "SOR_MOD_REDUCTION              % d\n",      SOR_MOD_REDUCTION       );
#     endif // #if ( defined GRAVITY  &&  POT_SCHEME == SOR  &&  defined GPU )

#     ifdef GPU
#     ifdef DT_FLU_USE_SHUFFLE
      fprintf( Note, "DT_FLU_USE_SHUFFLE              ON\n" );
#     else
      fprintf( Note, "DT_FLU_USE_SHUFFLE              OFF\n" );
#     endif
#     ifdef GRAVITY
#     ifdef DT_GRA_USE_SHUFFLE
      fprintf( Note, "DT_GRA_USE_SHUFFLE              ON\n" );
#     else
      fprintf( Note, "DT_GRA_USE_SHUFFLE              OFF\n" );
#     endif
#     endif // #ifdef GRAVITY
#     endif // #ifdef GPU

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the symbolic constants
      fprintf( Note, "Symbolic Constants\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "#define VERSION                 %s\n",      VERSION               );
      fprintf( Note, "#define NCOMP_FLUID            % d\n",      NCOMP_FLUID           );
      fprintf( Note, "#define NCOMP_PASSIVE          % d\n",      NCOMP_PASSIVE         );
      fprintf( Note, "#define FLU_NIN                % d\n",      FLU_NIN               );
      fprintf( Note, "#define FLU_NOUT               % d\n",      FLU_NOUT              );
      fprintf( Note, "#define FLU_NIN_T              % d\n",      FLU_NIN_T             );
      fprintf( Note, "#define FLU_NIN_S              % d\n",      FLU_NIN_S             );
      fprintf( Note, "#define FLU_NOUT_S             % d\n",      FLU_NOUT_S            );
      fprintf( Note, "#define DER_NOUT_MAX           % d\n",      DER_NOUT_MAX          );
      fprintf( Note, "#define NFIELD_STORED_MAX      % d\n",      NFIELD_STORED_MAX     );
      fprintf( Note, "#define NFLUX_FLUID            % d\n",      NFLUX_FLUID           );
      fprintf( Note, "#define NFLUX_PASSIVE          % d\n",      NFLUX_PASSIVE         );
#     ifdef GRAVITY
      fprintf( Note, "#define GRA_NIN                % d\n",      GRA_NIN               );
#     endif
#     ifdef MHD
      fprintf( Note, "#define NCOMP_MAG              % d\n",      NCOMP_MAG             );
      fprintf( Note, "#define NCOMP_ELE              % d\n",      NCOMP_ELE             );
#     endif
      fprintf( Note, "#define PATCH_SIZE             % d\n",      PATCH_SIZE            );
      fprintf( Note, "#define MAX_PATCH              % d\n",      MAX_PATCH             );
      fprintf( Note, "#define NLEVEL                 % d\n",      NLEVEL                );
      fprintf( Note, "\n" );
      fprintf( Note, "#define FLU_GHOST_SIZE         % d\n",      FLU_GHOST_SIZE        );
#     if ( MODEL == HYDRO  &&  defined LR_GHOST_SIZE )
      fprintf( Note, "#define LR_GHOST_SIZE          % d\n",      LR_GHOST_SIZE         );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "#define POT_GHOST_SIZE         % d\n",      POT_GHOST_SIZE        );
      fprintf( Note, "#define RHO_GHOST_SIZE         % d\n",      RHO_GHOST_SIZE        );
      fprintf( Note, "#define GRA_GHOST_SIZE         % d\n",      GRA_GHOST_SIZE        );
#     ifdef UNSPLIT_GRAVITY
      fprintf( Note, "#define USG_GHOST_SIZE_F       % d\n",      USG_GHOST_SIZE_F      );
      fprintf( Note, "#define USG_GHOST_SIZE_G       % d\n",      USG_GHOST_SIZE_G      );
#     endif
#     ifdef PARTICLE
      fprintf( Note, "#define RHOEXT_GHOST_SIZE      % d\n",      RHOEXT_GHOST_SIZE     );
#     endif
#     endif // #ifdef GRAVITY
      fprintf( Note, "#define SRC_GHOST_SIZE          %d\n",      SRC_GHOST_SIZE        );
      fprintf( Note, "#define DER_GHOST_SIZE          %d\n",      DER_GHOST_SIZE        );
#     ifdef FEEDBACK
      fprintf( Note, "#define FB_GHOST_SIZE           %d\n",      FB_GHOST_SIZE         );
#     endif
      fprintf( Note, "#define FLU_NXT                 %d\n",      FLU_NXT               );
#     ifdef GRAVITY
      fprintf( Note, "#define POT_NXT                % d\n",      POT_NXT               );
      fprintf( Note, "#define RHO_NXT                % d\n",      RHO_NXT               );
      fprintf( Note, "#define GRA_NXT                % d\n",      GRA_NXT               );
#     ifdef UNSPLIT_GRAVITY
      fprintf( Note, "#define USG_NXT_F              % d\n",      USG_NXT_F             );
      fprintf( Note, "#define USG_NXT_G              % d\n",      USG_NXT_G             );
#     endif
#     endif // #ifdef GRAVITY
#     ifdef MASSIVE_PARTICLES
      fprintf( Note, "#define RHOEXT_NXT             % d\n",      RHOEXT_NXT          );
#     endif
      fprintf( Note, "#define SRC_NXT                % d\n",      SRC_NXT               );
      fprintf( Note, "#define DER_NXT                % d\n",      DER_NXT               );
#     ifdef FEEDBACK
      fprintf( Note, "#define FB_NXT                 % d\n",      FB_NXT                );
#     endif
      fprintf( Note, "#define SRC_NXT                 %d\n",      SRC_NXT               );
      fprintf( Note, "#define DER_NXT                 %d\n",      DER_NXT               );
#     ifdef FEEDBACK
      fprintf( Note, "#define FB_NXT                  %d\n",      FB_NXT                );
#     endif
#     if ( MODEL == HYDRO )
      fprintf( Note, "#define EOS_NAUX_MAX           % d\n",      EOS_NAUX_MAX          );
      fprintf( Note, "#define EOS_NTABLE_MAX         % d\n",      EOS_NTABLE_MAX        );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "#define EXT_POT_NAUX_MAX       % d\n",      EXT_POT_NAUX_MAX      );
      fprintf( Note, "#define EXT_ACC_NAUX_MAX       % d\n",      EXT_ACC_NAUX_MAX      );
      fprintf( Note, "#define EXT_POT_NGENE_MAX      % d\n",      EXT_POT_NGENE_MAX     );
#     endif
      fprintf( Note, "#define SRC_NAUX_DLEP          % d\n",      SRC_NAUX_DLEP         );
      fprintf( Note, "#define SRC_DLEP_PROF_NVAR     % d\n",      SRC_DLEP_PROF_NVAR    );
      fprintf( Note, "#define SRC_DLEP_PROF_NBINMAX  % d\n",      SRC_DLEP_PROF_NBINMAX );
      fprintf( Note, "#define SRC_NAUX_USER          % d\n",      SRC_NAUX_USER         );
#     ifdef GPU
      fprintf( Note, "#define FLU_BLOCK_SIZE_X       % d\n",      FLU_BLOCK_SIZE_X      );
      fprintf( Note, "#define FLU_BLOCK_SIZE_Y       % d\n",      FLU_BLOCK_SIZE_Y      );
#     ifdef GRAVITY
#     if   ( POT_SCHEME == SOR )
      fprintf( Note, "#define POT_BLOCK_SIZE_Z       % d\n",      POT_BLOCK_SIZE_Z      );
#     elif ( POT_SCHEME == MG )
      fprintf( Note, "#define POT_BLOCK_SIZE_X       % d\n",      POT_BLOCK_SIZE_X      );
#     endif
      fprintf( Note, "#define EXTPOT_BLOCK_SIZE      % d\n",      EXTPOT_BLOCK_SIZE     );
      fprintf( Note, "#define GRA_BLOCK_SIZE         % d\n",      GRA_BLOCK_SIZE        );
#     endif // #ifdef GRAVITY
      fprintf( Note, "#define DT_FLU_BLOCK_SIZE      % d\n",      DT_FLU_BLOCK_SIZE     );
#     ifdef GRAVITY
      fprintf( Note, "#define DT_GRA_BLOCK_SIZE      % d\n",      DT_GRA_BLOCK_SIZE     );
#     endif
      fprintf( Note, "#define SRC_BLOCK_SIZE         % d\n",      SRC_BLOCK_SIZE        );
#     endif // #ifdef GPU
#     ifdef PARTICLE
      fprintf( Note, "#define PAR_NATT_TOTAL         % d\n",      PAR_NATT_TOTAL        );
      fprintf( Note, "#define PAR_NATT_USER          % d\n",      PAR_NATT_USER         );
      fprintf( Note, "#define PAR_NATT_STORED        % d\n",      PAR_NATT_STORED       );
      fprintf( Note, "#define PAR_NTYPE              % d\n",      PAR_NTYPE             );
#     endif
      fprintf( Note, "#define MAX_STRING             % d\n",      MAX_STRING            );
      fprintf( Note, "#define TINY_NUMBER            % 21.14e\n", TINY_NUMBER           );
      fprintf( Note, "#define HUGE_NUMBER            % 21.14e\n", HUGE_NUMBER           );
      fprintf( Note, "#define MAX_ERROR              % 21.14e\n", MAX_ERROR             );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of simulation scale
      fprintf( Note, "Parameters of Simulation Scale\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "BOX_SIZE (input)               % 21.14e\n", BOX_SIZE         );
      fprintf( Note, "BOX_SIZE_X                     % 21.14e\n", amr->BoxSize[0]  );
      fprintf( Note, "BOX_SIZE_Y                     % 21.14e\n", amr->BoxSize[1]  );
      fprintf( Note, "BOX_SIZE_Z                     % 21.14e\n", amr->BoxSize[2]  );
      fprintf( Note, "BOX_SCALE_X                    % d\n",      amr->BoxScale[0] );
      fprintf( Note, "BOX_SCALE_Y                    % d\n",      amr->BoxScale[1] );
      fprintf( Note, "BOX_SCALE_Z                    % d\n",      amr->BoxScale[2] );
      fprintf( Note, "NX0_TOT[0]                     % d\n",      NX0_TOT[0]       );
      fprintf( Note, "NX0_TOT[1]                     % d\n",      NX0_TOT[1]       );
      fprintf( Note, "NX0_TOT[2]                     % d\n",      NX0_TOT[2]       );
      fprintf( Note, "MPI_NRank                      % d\n",      MPI_NRank        );
      fprintf( Note, "MPI_NRank_X[0]                 % d\n",      MPI_NRank_X[0]   );
      fprintf( Note, "MPI_NRank_X[1]                 % d\n",      MPI_NRank_X[1]   );
      fprintf( Note, "MPI_NRank_X[2]                 % d\n",      MPI_NRank_X[2]   );
      fprintf( Note, "OMP_NTHREAD                    % d\n",      OMP_NTHREAD      );
      fprintf( Note, "END_T                          % 21.14e\n", END_T            );
      fprintf( Note, "END_STEP                       % ld\n",     END_STEP         );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of test problems
      fprintf( Note, "Parameters of Test Problems\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "TESTPROB_ID                    % d\n",      TESTPROB_ID      );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of code units
      fprintf( Note, "Parameters of Code Units\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__UNIT                      % d\n",      OPT__UNIT        );
      if ( OPT__UNIT ) {
#     ifdef COMOVING
      fprintf( Note, "\n### All units marked with (*) assume h =% 14.7e ###\n\n", HUBBLE0 );

      const double Current_Matter_Density = OMEGA_M0*3*SQR( 100.0*HUBBLE0*Const_km/Const_Mpc/Const_s )/( 8.0*M_PI*Const_NewtonG );
      fprintf( Note, "rho_bg = current matter density =% 21.14e Msun/kpc^3 (*)\n\n",
               Current_Matter_Density/(Const_Msun/CUBE(Const_kpc)) );

      fprintf( Note, "UNIT_L (length)                % 21.14e Mpc/h\n",          UNIT_L/(Const_Mpc/HUBBLE0)    );
      fprintf( Note, "                              =% 21.14e cm         (*)\n", UNIT_L                        );
      fprintf( Note, "UNIT_M (mass)                  % 21.14e Msun/h\n",         UNIT_M/(Const_Msun/HUBBLE0)   );
      fprintf( Note, "                              =% 21.14e g          (*)\n", UNIT_M                        );
      fprintf( Note, "UNIT_T (time)                  % 21.14e Gyr        (*)\n", UNIT_T/Const_Gyr              );
      fprintf( Note, "                              =% 21.14e s          (*)\n", UNIT_T                        );
      fprintf( Note, "UNIT_V (velocity)              % 21.14e km/s\n",           UNIT_V/(Const_km/Const_s)     );
      fprintf( Note, "                              =% 21.14e cm/s\n",           UNIT_V                        );
      fprintf( Note, "UNIT_D (mass density)          % 21.14e rho_bg     (*)\n", UNIT_D/Current_Matter_Density );
      fprintf( Note, "                              =% 21.14e g/cm^3     (*)\n", UNIT_D                        );
      fprintf( Note, "UNIT_E (energy)                % 21.14e g*cm^2/s^2 (*)\n", UNIT_E                        );
      fprintf( Note, "UNIT_P (energy density)        % 21.14e g/cm/s^2   (*)\n", UNIT_P                        );
#     ifdef MHD
#     warning : ERROR : MHD is not supported here !!!
#     endif

#     else

      fprintf( Note, "UNIT_L                         % 21.14e cm\n",             UNIT_L                        );
      fprintf( Note, "UNIT_M                         % 21.14e g\n",              UNIT_M                        );
      fprintf( Note, "UNIT_T                         % 21.14e s\n",              UNIT_T                        );
      fprintf( Note, "UNIT_V                         % 21.14e cm/s\n",           UNIT_V                        );
      fprintf( Note, "UNIT_D                         % 21.14e g/cm^3\n",         UNIT_D                        );
      fprintf( Note, "UNIT_E (energy)                % 21.14e g*cm^2/s^2\n",     UNIT_E                        );
      fprintf( Note, "UNIT_P (energy density)        % 21.14e g/cm/s^2\n",       UNIT_P                        );
#     ifdef MHD
      fprintf( Note, "UNIT_B (magnetic field)        % 21.14e gauss\n",          UNIT_B                        );
#     endif
#     endif // #ifdef COMOVING ... else ...
      }

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of boundary condition
      fprintf( Note, "Parameters of Boundary Condition\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__BC_FLU[0] (-x)            % d\n",      OPT__BC_FLU[0] );
      fprintf( Note, "OPT__BC_FLU[1] (+x)            % d\n",      OPT__BC_FLU[1] );
      fprintf( Note, "OPT__BC_FLU[2] (-y)            % d\n",      OPT__BC_FLU[2] );
      fprintf( Note, "OPT__BC_FLU[3] (+y)            % d\n",      OPT__BC_FLU[3] );
      fprintf( Note, "OPT__BC_FLU[4] (-z)            % d\n",      OPT__BC_FLU[4] );
      fprintf( Note, "OPT__BC_FLU[5] (+z)            % d\n",      OPT__BC_FLU[5] );
#     ifdef GRAVITY
      fprintf( Note, "OPT__BC_POT                    % d\n",      OPT__BC_POT    );
      fprintf( Note, "GFUNC_COEFF0                   % 14.7e\n",  GFUNC_COEFF0   );
#     endif
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of particle
#     ifdef PARTICLE
      fprintf( Note, "Parameters of Particle\n" );
      fprintf( Note, "***********************************************************************************\n" );
#     ifdef DEBUG_PARTICLE
      fprintf( Note, "DEBUG_PARTICLE                  ON\n" );
#     else
      fprintf( Note, "DEBUG_PARTICLE                  OFF\n" );
#     endif
      fprintf( Note, "Par->NPar_Active_AllRank       % ld\n",     amr->Par->NPar_Active_AllRank );
      fprintf( Note, "Par->Init                      % d\n",      amr->Par->Init                );
      fprintf( Note, "Par->ParICFormat               % d\n",      amr->Par->ParICFormat         );
      fprintf( Note, "PAR_IC_FLOAT8                  % d\n",      PAR_IC_FLOAT8                 );
      fprintf( Note, "Par->ParICMass                 % 14.7e\n",  amr->Par->ParICMass           );
      fprintf( Note, "Par->ParICType                 % d\n",      amr->Par->ParICType           );
      fprintf( Note, "Par->Interp                    % d\n",      amr->Par->Interp              );
      fprintf( Note, "Par->Integ                     % d\n",      amr->Par->Integ               );
      fprintf( Note, "Par->GhostSize                 % d\n",      amr->Par->GhostSize           );
      fprintf( Note, "Par->ImproveAcc                % d\n",      amr->Par->ImproveAcc          );
      fprintf( Note, "Par->PredictPos                % d\n",      amr->Par->PredictPos          );
      fprintf( Note, "Par->RemoveCell                % 14.7e\n",  amr->Par->RemoveCell          );
      fprintf( Note, "Par->InterpTracer              % d\n",      amr->Par->InterpTracer        );
      fprintf( Note, "Par->IntegTracer               % d\n",      amr->Par->IntegTracer         );
      fprintf( Note, "Par->GhostSizeTracer           % d\n",      amr->Par->GhostSizeTracer     );
      fprintf( Note, "Par->TracerVelCorr             % d\n",      amr->Par->TracerVelCorr       );
      fprintf( Note, "OPT__FREEZE_PAR                % d\n",      OPT__FREEZE_PAR               );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif


//    record the parameters of cosmological simulations (comoving frame)
#     ifdef COMOVING
      fprintf( Note, "Parameters of Cosmological Simulation\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "A_INIT                         % 14.7e\n",  A_INIT   );
      fprintf( Note, "OMEGA_M0                       % 14.7e\n",  OMEGA_M0 );
      fprintf( Note, "HUBBLE0 (h)                    % 14.7e\n",  HUBBLE0  );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif


//    record the parameters of time-step determination
      fprintf( Note, "Parameters of Time-step Determination\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "DT__MAX                        % 14.7e\n",  DT__MAX                     );
      fprintf( Note, "DT__FLUID                      % 14.7e\n",  DT__FLUID                   );
      fprintf( Note, "DT__FLUID_INIT                 % 14.7e\n",  DT__FLUID_INIT              );
#     ifdef SRHD
      fprintf( Note, "DT__SPEED_OF_LIGHT             % d\n",      DT__SPEED_OF_LIGHT          );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "DT__GRAVITY                    % 14.7e\n",  DT__GRAVITY                 );
#     endif
#     if ( MODEL == ELBDM )
      fprintf( Note, "DT__PHASE                      % 14.7e\n",  DT__PHASE                   );
#     endif
#     ifdef PARTICLE
      fprintf( Note, "DT__PARVEL                     % 14.7e\n",  DT__PARVEL                  );
      fprintf( Note, "DT__PARVEL_MAX                 % 14.7e\n",  DT__PARVEL_MAX              );
      fprintf( Note, "DT__PARACC                     % 14.7e\n",  DT__PARACC                  );
#     endif
#     ifdef CR_DIFFUSION
      fprintf( Note, "DT__CR_DIFFUSION               % 14.7e\n",  DT__CR_DIFFUSION            );
#     endif
#     ifdef COMOVING
      fprintf( Note, "DT__MAX_DELTA_A                % 14.7e\n",  DT__MAX_DELTA_A             );
#     endif
      fprintf( Note, "DT__SYNC_PARENT_LV             % 14.7e\n",  DT__SYNC_PARENT_LV          );
      fprintf( Note, "DT__SYNC_CHILDREN_LV           % 14.7e\n",  DT__SYNC_CHILDREN_LV        );
      fprintf( Note, "OPT__DT_USER                   % d\n",      OPT__DT_USER                );
      fprintf( Note, "OPT__DT_LEVEL                  % d\n",      OPT__DT_LEVEL               );
      fprintf( Note, "AUTO_REDUCE_DT                 % d\n",      AUTO_REDUCE_DT              );
      fprintf( Note, "AUTO_REDUCE_DT_FACTOR          % 14.7e\n",  AUTO_REDUCE_DT_FACTOR       );
      fprintf( Note, "AUTO_REDUCE_DT_FACTOR_MIN      % 14.7e\n",  AUTO_REDUCE_DT_FACTOR_MIN   );
#     if ( MODEL == HYDRO )
      fprintf( Note, "AUTO_REDUCE_MINMOD_FACTOR      % 14.7e\n",  AUTO_REDUCE_MINMOD_FACTOR   );
      fprintf( Note, "AUTO_REDUCE_MINMOD_MIN         % 14.7e\n",  AUTO_REDUCE_MINMOD_MIN      );
#     endif
      fprintf( Note, "AUTO_REDUCE_INT_MONO_FACTOR    % 14.7e\n",  AUTO_REDUCE_INT_MONO_FACTOR );
      fprintf( Note, "AUTO_REDUCE_INT_MONO_MIN       % 14.7e\n",  AUTO_REDUCE_INT_MONO_MIN    );
      fprintf( Note, "OPT__RECORD_DT                 % d\n",      OPT__RECORD_DT              );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of domain refinement
      fprintf( Note, "Parameters of Domain Refinement\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of yt inline analysis
#     ifdef SUPPORT_LIBYT
      fprintf( Note, "Parameters of YT Inline Analysis\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "YT_SCRIPT                           %s\n",      YT_SCRIPT  );
      fprintf( Note, "YT_VERBOSE                         % d\n",      YT_VERBOSE );
      fprintf( Note, "YT_FIG_BASENAME                     %s\n",      YT_FIG_BASENAME );
#     ifdef LIBYT_JUPYTER
      fprintf( Note, "YT_JUPYTER_USE_CONNECTION_FILE     % d\n",      YT_JUPYTER_USE_CONNECTION_FILE );
#     endif
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif


//    record the parameters of miscellaneous purposes
      fprintf( Note, "Parameters of Miscellaneous Purposes\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__VERBOSE                   % d\n",      OPT__VERBOSE             );
      fprintf( Note, "OPT__TIMING_BARRIER            % d\n",      OPT__TIMING_BARRIER      );
      fprintf( Note, "OPT__TIMING_BALANCE            % d\n",      OPT__TIMING_BALANCE      );
      fprintf( Note, "OPT__TIMING_MPI                % d\n",      OPT__TIMING_MPI          );
      fprintf( Note, "OPT__RECORD_NOTE               % d\n",      OPT__RECORD_NOTE         );
      fprintf( Note, "OPT__RECORD_UNPHY              % d\n",      OPT__RECORD_UNPHY        );
      fprintf( Note, "OPT__RECORD_MEMORY             % d\n",      OPT__RECORD_MEMORY       );
      fprintf( Note, "OPT__RECORD_PERFORMANCE        % d\n",      OPT__RECORD_PERFORMANCE  );
      fprintf( Note, "OPT__RECORD_CENTER             % d\n",      OPT__RECORD_CENTER       );
      if ( OPT__RECORD_CENTER )
      {
      fprintf( Note, "   COM_CEN_X                   % 14.7e\n",  COM_CEN_X                );
      fprintf( Note, "   COM_CEN_Y                   % 14.7e\n",  COM_CEN_Y                );
      fprintf( Note, "   COM_CEN_Z                   % 14.7e\n",  COM_CEN_Z                );
      fprintf( Note, "   COM_MAX_R                   % 14.7e\n",  COM_MAX_R                );
      fprintf( Note, "   COM_MIN_RHO                 % 14.7e\n",  COM_MIN_RHO              );
      fprintf( Note, "   COM_TOLERR_R                % 14.7e\n",  COM_TOLERR_R             );
      fprintf( Note, "   COM_MAX_ITER                % d\n",      COM_MAX_ITER             );
      }
      fprintf( Note, "OPT__MANUAL_CONTROL            % d\n",      OPT__MANUAL_CONTROL      );
      fprintf( Note, "OPT__RECORD_USER               % d\n",      OPT__RECORD_USER         );
      fprintf( Note, "OPT__OPTIMIZE_AGGRESSIVE       % d\n",      OPT__OPTIMIZE_AGGRESSIVE );
      fprintf( Note, "OPT__SORT_PATCH_BY_LBIDX       % d\n",      OPT__SORT_PATCH_BY_LBIDX );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of simulation checks
      fprintf( Note, "Parameters of Simulation Checks\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__CK_REFINE                 % d\n",      OPT__CK_REFINE            );
      fprintf( Note, "OPT__CK_PROPER_NESTING         % d\n",      OPT__CK_PROPER_NESTING    );
      fprintf( Note, "OPT__CK_CONSERVATION           % d\n",      OPT__CK_CONSERVATION      );
      if ( OPT__CK_CONSERVATION )
      {
      fprintf( Note, "   ANGMOM_ORIGIN_X             % 14.7e\n",  ANGMOM_ORIGIN_X           );
      fprintf( Note, "   ANGMOM_ORIGIN_Y             % 14.7e\n",  ANGMOM_ORIGIN_Y           );
      fprintf( Note, "   ANGMOM_ORIGIN_Z             % 14.7e\n",  ANGMOM_ORIGIN_Z           );
      }
      fprintf( Note, "OPT__CK_NORMALIZE_PASSIVE      % d\n",      OPT__CK_NORMALIZE_PASSIVE );
      fprintf( Note, "OPT__CK_RESTRICT               % d\n",      OPT__CK_RESTRICT          );
      fprintf( Note, "OPT__CK_FINITE                 % d\n",      OPT__CK_FINITE            );
      fprintf( Note, "OPT__CK_PATCH_ALLOCATE         % d\n",      OPT__CK_PATCH_ALLOCATE    );
      fprintf( Note, "OPT__CK_FLUX_ALLOCATE          % d\n",      OPT__CK_FLUX_ALLOCATE     );
#     if ( MODEL == HYDRO )
      fprintf( Note, "OPT__CK_NEGATIVE               % d\n",      OPT__CK_NEGATIVE          );
#     endif
      fprintf( Note, "OPT__CK_MEMFREE                % 14.7e\n",  OPT__CK_MEMFREE           );
#     ifdef PARTICLE
      fprintf( Note, "OPT__CK_PARTICLE               % d\n",      OPT__CK_PARTICLE          );
#     endif
#     ifdef MHD
      fprintf( Note, "OPT__CK_INTERFACE_B            % d\n",      OPT__CK_INTERFACE_B       );
      fprintf( Note, "OPT__CK_DIVERGENCE_B           % d\n",      OPT__CK_DIVERGENCE_B      );
#     endif
      fprintf( Note, "OPT__CK_INPUT_FLUID            % d\n",      OPT__CK_INPUT_FLUID       );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the flag criterion (density/density gradient/pressure gradient/user-defined)
      if ( OPT__FLAG_RHO )
      {
         fprintf( Note, "Flag Criterion (Density)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level             Density\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Rho[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_RHO_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Density Gradient)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level    Density Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_RhoGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

#     if   ( MODEL == HYDRO )
      if ( OPT__FLAG_PRES_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Pressure Gradient in HYDRO)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level   Pressure Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_PresGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_VORTICITY )
      {
         fprintf( Note, "Flag Criterion (Vorticity in HYDRO)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level           Vorticity\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Vorticity[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_JEANS )
      {
         fprintf( Note, "Flag Criterion (Jeans Length over Cell Size in HYDRO)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level       lambda_J / dh\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Jeans[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

#     ifdef MHD
      if ( OPT__FLAG_CURRENT )
      {
         fprintf( Note, "Flag Criterion (Current Density in MHD)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level             Current\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Current[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif

#     ifdef SRHD
      if ( OPT__FLAG_LRTZ_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Lorentz Factor Gradient in SRHD)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level   Lorentz Factor Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%26.7e\n", lv, FlagTable_LrtzGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif

#     ifdef COSMIC_RAY
      if ( OPT__FLAG_CRAY )
      {
         fprintf( Note, "Flag Criterion (Cosmic Ray Energy)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level             Cosmic Ray Energy\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_CRay[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif
#     endif // #if ( MODEL == HYDRO )

#     if ( MODEL == ELBDM )
      if ( OPT__FLAG_ENGY_DENSITY )
      {
         fprintf( Note, "Flag Criterion (Energy Density in ELBDM)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level     Angle_over_2*PI              Soften\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)
            fprintf( Note, "%7d%20.7e%20.7e\n", lv, FlagTable_EngyDensity[lv][0], FlagTable_EngyDensity[lv][1] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif

#     if   ( MODEL == HYDRO )
#     ifndef COSMIC_RAY
      const bool OPT__FLAG_LOHNER_CRAY = false;
#     endif
      if ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES || OPT__FLAG_LOHNER_TEMP ||
           OPT__FLAG_LOHNER_ENTR || OPT__FLAG_LOHNER_CRAY )
#     elif ( MODEL == ELBDM )
      if ( OPT__FLAG_LOHNER_DENS )
#     endif
      {
         fprintf( Note, "Flag Criterion (Lohner Error Estimator)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level    Threshold_Refine  Threshold_Derefine              Filter              Soften      MinimumDensity\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)
            fprintf( Note, "%7d%20.7e%20.7e%20.7e%20.7e%20.7e\n", lv, FlagTable_Lohner[lv][0], FlagTable_Lohner[lv][1],
                     FlagTable_Lohner[lv][2], FlagTable_Lohner[lv][3], FlagTable_Lohner[lv][4] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_USER )
      {
         fprintf( Note, "Flag Criterion (User-defined)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level           Threshold\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)
         {
            fprintf( Note, "%7d",    lv );
            for (int t=0; t<OPT__FLAG_USER_NUM; t++)   fprintf( Note, "%20.7e", FlagTable_User[lv][t] );
            fprintf( Note, "\n" );
         }
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

#     ifdef PARTICLE
      if ( OPT__FLAG_NPAR_PATCH )
      {
         fprintf( Note, "Flag Criterion (# of Particles per Patch)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level      # of Particles\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20d\n", lv, FlagTable_NParPatch[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_NPAR_CELL )
      {
         fprintf( Note, "Flag Criterion (# of Particles per Cell)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level      # of Particles\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20d\n", lv, FlagTable_NParCell[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_PAR_MASS_CELL )
      {
         fprintf( Note, "Flag Criterion (Particle Mass per Cell)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level       Particle Mass\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_ParMassCell[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif // #ifdef PARTICLE


//    record the grid size in different refinement level
      fprintf( Note, "Cell Size and Scale (scale = number of cells at the finest level)\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "%7s%*c%26s%*c%16s\n", "Level", 10, ' ', "Cell Size", 10, ' ', "Cell Scale" );
      for (int lv=0; lv<NLEVEL; lv++)
      fprintf( Note, "%7d%*c%26.20lf%*c%16d\n", lv, 10, ' ', amr->dh[lv], 10, ' ', amr->scale[lv] );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the compilation time of the file Aux_TakeNote.cpp
      fprintf( Note, "Compilation Time\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "%s %s\n", __DATE__, __TIME__ );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the current time when running GAMER
      time_t t = time( NULL );
      fprintf( Note, "Current Time\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "%s", ctime( &t ) );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the git information
      fprintf( Note, "Git information\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "Branch : %s\n", EXPAND_AND_QUOTE(GIT_BRANCH) );
      fprintf( Note, "Commit : %s\n", EXPAND_AND_QUOTE(GIT_COMMIT) );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n" );


      fclose( Note );
   } // if ( MPI_Rank == 0 )


// record the hostname and PID of each MPI process (the function "CUAPI_DiagnoseDevice" will also record them)
#  ifndef GPU
   const int PID = getpid();
   char Host[1024];
   gethostname( Host, 1024 );

   if ( MPI_Rank == 0 )
   {
       Note = fopen( FileName, "a" );
       fprintf( Note, "Device Diagnosis\n" );
       fprintf( Note, "***********************************************************************************\n" );
       fclose( Note );
   }

   for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)
   {
      if ( MPI_Rank == YourTurn )
      {
         Note = fopen( FileName, "a" );
         if ( MPI_Rank != 0 )    fprintf( Note, "\n" );
         fprintf( Note, "MPI_Rank = %3d, hostname = %10s, PID = %5d\n", MPI_Rank, Host, PID );
         fprintf( Note, "CPU Info :\n" );
         fflush( Note );

         Aux_GetCPUInfo( FileName );

         fclose( Note );
      }

      MPI_Barrier( MPI_COMM_WORLD );
   }

   if ( MPI_Rank == 0 )
   {
      Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
      fclose( Note );
   }
#  endif // #ifndef GPU


// record the OpenMP status
#  ifdef OPENMP
   int omp_nthread, omp_chunk_size, omp_nested;
   omp_sched_t omp_schedule;

   omp_nested = omp_get_max_active_levels();
   omp_get_schedule( &omp_schedule, &omp_chunk_size );

#  pragma omp parallel
#  pragma omp master
   { omp_nthread = omp_get_num_threads(); }

   if ( MPI_Rank == 0 )
   {
      Note = fopen( FileName, "a" );
      fprintf( Note, "OpenMP Diagnosis\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "Schedule                        %s\n",      ( omp_schedule == omp_sched_static  ) ? "STATIC"  :
                                                                  ( omp_schedule == omp_sched_dynamic ) ? "DYNAMIC" :
                                                                  ( omp_schedule == omp_sched_guided  ) ? "GUIDED"  :
                                                                  ( omp_schedule == omp_sched_auto    ) ? "AUTO"    : "UNKNOWN" );
      fprintf( Note, "Chunk size                     % d\n",      omp_chunk_size );
      fprintf( Note, "Max number of nested levels    % d\n",      omp_nested );
      fprintf( Note, "\n" );
      fprintf( Note, "CPU core IDs of all OpenMP threads (tid == thread ID):\n" );
      fprintf( Note, "------------------------------------------------------------------------\n" );
      fprintf( Note, "%5s  %10s  %7s", "Rank", "Host", "NThread" );
      for (int t=0; t<omp_nthread; t++)
      fprintf( Note, "  tid-%02d", t );
      fprintf( Note, "\n" );
      fclose( Note );
   }

// record the CPU core id of each OpenMP thread
   int *omp_core_id = new int [omp_nthread];

#  pragma omp parallel
   { omp_core_id[ omp_get_thread_num() ] = get_cpuid(); }

   for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)
   {
      if ( MPI_Rank == YourTurn )
      {
         char MPI_Host[1024];
         gethostname( MPI_Host, 1024 );

         Note = fopen( FileName, "a" );
         fprintf( Note, "%5d  %10s  %7d", MPI_Rank, MPI_Host, omp_nthread );
         for (int t=0; t<omp_nthread; t++)   fprintf( Note, "  %6d", omp_core_id[t] );
         fprintf( Note, "\n" );
         fflush( Note );
         fclose( Note );
      }

      MPI_Barrier( MPI_COMM_WORLD );
   }

   delete [] omp_core_id;

   if ( MPI_Rank == 0 )
   {
      Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
      fclose( Note );
   }
#  endif // #ifdef OPENMP


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_TakeNote ... done\n" );

} // FUNCTION : Aux_TakeNote



//-------------------------------------------------------------------------------------------------------
// Function    :  get_cpuid
// Description :  Get the CPU ID
//
// Note        :  Work on both macOS and Linux systems
//-------------------------------------------------------------------------------------------------------
int get_cpuid()
{

// See https://stackoverflow.com/questions/33745364/sched-getcpu-equivalent-for-os-x
   int CPU;

#  ifdef __APPLE__
   uint32_t CPUInfo[4];
   __cpuid_count(1, 0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
   if ((CPUInfo[3] & (1 << 9)) == 0) {
      CPU = -1;  /* no APIC on chip */
   } else {
      CPU = (unsigned)CPUInfo[1] >> 24;
   }
   if (CPU < 0) CPU = 0;
#  else
   CPU = sched_getcpu();
#  endif

   return CPU;

} // FUNCTION : get_cpuid
