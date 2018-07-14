#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_PassiveVariable
// Description :  Initialize settings for the passive variables
//
// Note        :  1. Set field names for passively advected scalars and passive particle attributes
//                   --> Default names are "PassiveXX", where XX = [0 ... NCOMP_PASSIVE-1]
//                       --> Can be modified by users
//                   --> These names are only used in the output files
//                       --> To set symbolic constants used internally for different passive variables,
//                           please edit the header "Macro.h" and follow the example "METAL"
//                2. Set passive scalars to be normalized (for OPT__NORMALIZE_PASSIVE)
//                   --> Default is to normalize all passive scalars except for the internal energy
//                       used by the dual-energy formalism
//                       --> Can be modified by users
//
// Return      :  PassiveFieldName_Grid, PassiveFieldName_Par, PassiveNorm_NVar, PassiveNorm_VarIdx
//-------------------------------------------------------------------------------------------------------
void Init_PassiveVariable()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. Set field names
// ==============================================================================================
// 1-1. allocate memory
// passive scalars on grids
   for (int v=0; v<NCOMP_PASSIVE; v++)    PassiveFieldName_Grid[v] = new char [MAX_STRING];

// passive attributes of particles
#  ifdef PARTICLE
   for (int v=0; v<PAR_NPASSIVE; v++)     PassiveFieldName_Par [v] = new char [MAX_STRING];
#  endif


// 1-2-a. default names
// passive scalars on grids
   for (int v=0; v<NCOMP_PASSIVE; v++)    sprintf( PassiveFieldName_Grid[v], "Passive%02d", v );

// passive attributes of particles
#  ifdef PARTICLE
   for (int v=0; v<PAR_NPASSIVE; v++)     sprintf( PassiveFieldName_Par[v],  "Passive%02d", v );
#  endif


// 1-2-b. example of user-specified field names
// (1) please comment out the sprintf() above
// (2) set field names as follows (assuming three passive fields for both grids and particles)
/*
#  if ( NCOMP_PASSIVE > 0 )
   sprintf( PassiveFieldName_Grid[0], "Metal" );
   sprintf( PassiveFieldName_Grid[1], "HI" );
   sprintf( PassiveFieldName_Grid[2], "HII" );
#  endif

#  if ( defined PARTICLE  &&  PAR_NPASSIVE > 0 )
   sprintf( PassiveFieldName_Par [0], "ParMetalFrac" );
   sprintf( PassiveFieldName_Par [1], "ParHIFrac" );
   sprintf( PassiveFieldName_Par [2], "ParHIIFrac" );
#  endif
*/


// 1-3. entropy (or internal energy) in the dual-energy formalism
// --> we always store it as the last passive variable
#  if   ( DUAL_ENERGY == DE_ENPY )
   sprintf( PassiveFieldName_Grid[NCOMP_PASSIVE-1], "Entropy" );
#  elif ( DUAL_ENERGY == DE_EINT )
   sprintf( PassiveFieldName_Grid[NCOMP_PASSIVE-1], "Eint" );
#  endif


// 1-4. particle creation time when STAR_FORMATION is enabled
// --> we always store it as the last passive attribute
#  if ( defined PARTICLE  &&  defined STAR_FORMATION )
   sprintf( PassiveFieldName_Par[PAR_NPASSIVE-1], "ParCreTime" );
#  endif



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_PassiveVariable



//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_PassiveFieldName
// Description :  Free memory allocated by Init_PassiveVariable()
//
// Note        :  1. Called by End_MemFree()
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_PassiveFieldName()
{

   for (int v=0; v<NCOMP_PASSIVE; v++)
   {
      delete [] PassiveFieldName_Grid[v];
      PassiveFieldName_Grid[v] = NULL;
   }

#  ifdef PARTICLE
   for (int v=0; v<PAR_NPASSIVE; v++)
   {
      delete [] PassiveFieldName_Par[v];
      PassiveFieldName_Par[v] = NULL;
   }
#  endif

} // FUNCTION : End_MemFree_PassiveFieldName

