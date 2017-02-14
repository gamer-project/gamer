#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_SetPassiveFieldName
// Description :  Set field names for passively advected scalars and passive particle attributes
//
// Note        :  1. Default names are "PassiveXX", where XX = [0 ... NCOMP_PASSIVE-1]
//                2. These names are only used in the output files
//                   --> To set symbolic constants used internally for different passive variables,
//                       please edit the header "Macro.h" and follow the example "METAL"
//
// Return      :  PassiveFieldName_Grid, PassiveFieldName_Par
//-------------------------------------------------------------------------------------------------------
void Init_SetPassiveFieldName()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// allocate memory
// ==============================================================================================
// passive scalars on grids
   for (int v=0; v<NCOMP_PASSIVE; v++)    PassiveFieldName_Grid[v] = new char [MAX_STRING];

// passive attributes of particles
#  ifdef PARTICLE
   for (int v=0; v<PAR_NPASSIVE; v++)     PassiveFieldName_Par [v] = new char [MAX_STRING];
#  endif


// default names
// ==============================================================================================
// passive scalars on grids
   for (int v=0; v<NCOMP_PASSIVE; v++)    sprintf( PassiveFieldName_Grid[v], "Passive%02d", v );

// passive attributes of particles
#  ifdef PARTICLE
   for (int v=0; v<PAR_NPASSIVE; v++)     sprintf( PassiveFieldName_Par[v],  "Passive%02d", v );
#  endif


// example of user-specified field names
// ==============================================================================================
// (1) please comment out the sprintf() above
// (2) set field names as follows (assuming three passive fields for both grids and particles)
/*
   sprintf( PassiveFieldName_Grid[0], "Metal" );
   sprintf( PassiveFieldName_Grid[1], "Oxygen" );
   sprintf( PassiveFieldName_Grid[2], "Fe" );

#  ifdef PARTICLE
   sprintf( PassiveFieldName_Par [0], "Par_Metal" );
   sprintf( PassiveFieldName_Par [1], "Par_Oxygen" );
   sprintf( PassiveFieldName_Par [2], "Par_Fe" );
#  endif
*/


// internal energy in the dual-energy formalism
// --> we always store it as the last passive variable
// ==============================================================================================
#  ifdef DUAL_ENERGY
   sprintf( PassiveFieldName_Grid[NCOMP_PASSIVE-1], "Eint" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_SetPassiveFieldName



//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_PassiveFieldName
// Description :  Free memory allocated by Init_SetPassiveFieldName()
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

