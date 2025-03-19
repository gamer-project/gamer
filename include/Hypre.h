#ifndef __HYPRE_H__
#define __HYPRE_H__



#ifdef GAMER_DEBUG
#  define DEBUG_HYPRE
#endif


#include "HYPRE_config.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"


// macro for checking Hypre function return
#ifdef DEBUG_HYPRE

#  define HYPRE_CHECK_FUNC( call )                                         \
   {                                                                       \
      int hypre_ierr = call;                                               \
                                                                           \
      if ( hypre_ierr )                                                    \
         Aux_Error( ERROR_INFO, "hypre error code: %d !!\n", hypre_ierr ); \
   }

#else

#  define HYPRE_CHECK_FUNC( call ) call

#endif



#endif // #ifndef __HYPRE_H__
