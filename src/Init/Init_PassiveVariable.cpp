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



// 2. Set target passive scalars on grids to be normalized
// --> ensure "sum(passive_scalar_density) == gas_density"
// ==============================================================================================
// 2-1. number of passive scalars to be normalized
   if ( OPT__NORMALIZE_PASSIVE )
   {
//    2-1-a. default --> all passive scalars except for the internal energy used by the dual-energy formalism
#     ifdef DUAL_ENERGY
      PassiveNorm_NVar = NCOMP_PASSIVE - 1;
#     else
      PassiveNorm_NVar = NCOMP_PASSIVE;
#     endif

//    2-1-b. example of user-specified variables
/*
      PassiveNorm_NVar = 2;
*/
   }

   else
      PassiveNorm_NVar = 0;


// 2-2. indices of passive scalars to be normalized
// --> note that here the indices of passive scalars start from 0 instead of NCOMP_FLUID
// --> also be careful not to include the internal energy used by the dual-energy formalism

// 2-2-1. initialize as -1 so that later on we can check whether it's set properly
   for (int v=0; v<NCOMP_PASSIVE; v++)    PassiveNorm_VarIdx[v] = -1;

   if ( PassiveNorm_NVar > 0 )
   {
//    2-2-a. default --> all passive scalars except for the internal energy used by the dual-energy formalism
      for (int v=0; v<PassiveNorm_NVar; v++)    PassiveNorm_VarIdx[v] = v;

//    2-2-b. example of user-specified variables
/*
      PassiveNorm_VarIdx[0] = 0;
      PassiveNorm_VarIdx[1] = 2;
*/
   }


// 2-3. check
   if ( OPT__NORMALIZE_PASSIVE  &&  PassiveNorm_NVar <= 0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : found no target passive scalars to be normalized even though OPT__NORMALIZE_PASSIVE is on !!\n" );

   if ( !OPT__NORMALIZE_PASSIVE  &&  PassiveNorm_NVar != 0 )
      Aux_Error( ERROR_INFO, "ERROR : PassiveNorm_NVar (%d) != 0 when OPT__NORMALIZE_PASSIVE is off !!\n", PassiveNorm_NVar );

   if ( OPT__NORMALIZE_PASSIVE )
   {
      const int MinIdx = 0;
#     ifdef DUAL_ENERGY
      const int MaxIdx = NCOMP_PASSIVE - 2;
#     else
      const int MaxIdx = NCOMP_PASSIVE - 1;
#     endif

      for (int v=0; v<PassiveNorm_NVar; v++)
      {
         if ( PassiveNorm_VarIdx[v] < MinIdx  ||  PassiveNorm_VarIdx[v] > MaxIdx )
            Aux_Error( ERROR_INFO, "PassiveNorm_VarIdx[%d] = %d is not within the correct range ([%d <= idx <= %d]) !!\n",
                       v, PassiveNorm_VarIdx[v], MinIdx, MaxIdx );
      }
   }


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

