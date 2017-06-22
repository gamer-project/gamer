#include "Copyright.h"
#include "GAMER.h"
#include <typeinfo>

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Init
// Description :  Initialize the chemistry and radiative cooling library Grackle
//
// Note        :  1. Must be called AFTER Init_Load_Parameter(), Init_Unit(), and Init_OpenMP()
//                2. COMOVING is not supported yet
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Grackle_Init()
{

// nothing to do if Grackle is disabled
   if ( GRACKLE_MODE == GRACKLE_MODE_NONE )  return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
// floating-point type (don't know how to validate it yet...)
   /*
   if ( typeid(real) != typeid(gr_float) )
      Aux_Error( ERROR_INFO, "inconsistent floating-point type (GAMER: %d, Grackle: %d) !!\n",
                 sizeof(real), sizeof(gr_float) );
                 */

// comoving frame is not supported yet
#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "SUPPORT_GRACKLE does not work with COMOVING yet !!\n" );
#  endif

   if ( GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_CLOUDY )
      Aux_Error( ERROR_INFO, "only support \"GRACKLE_PRIMORDIAL\" = %d for now !!\n", GRACKLE_PRI_CHE_CLOUDY );

   if ( GRACKLE_METAL )
      Aux_Error( ERROR_INFO, "\"GRACKLE_METAL\" is not supported yet !!\n" );

   if ( GRACKLE_UV )
      Aux_Error( ERROR_INFO, "\"GRACKLE_UV\" is not supported yet !!\n" );

   if (  ( GRACKLE_PRIMORDIAL == GRACKLE_PRI_CHE_CLOUDY || GRACKLE_METAL || GRACKLE_UV )  &&
         !Aux_CheckFileExist(GRACKLE_CLOUDY_TABLE)  )
      Aux_Error( ERROR_INFO, "Grackle data file \"%s\" does not exist !!\n", GRACKLE_CLOUDY_TABLE );


// enable output (only for the root rank)
   grackle_verbose = ( MPI_Rank == 0 ) ? GRACKLE_VERBOSE : 0;


// units in cgs
   code_units my_units;

#  ifdef COMOVING
   my_units.comoving_coordinates = 1;
   my_units.density_units        = NULL_REAL;   // not sure how to set the units in the comoving coordinates yet...
   my_units.length_units         = NULL_REAL;   // --> see http://grackle.readthedocs.io/en/latest/Integration.html
   my_units.time_units           = NULL_REAL;
   my_units.velocity_units       = NULL_REAL;
   my_units.a_units              = 1.0;
   my_units.a_value              = Time[0];

#  else
   my_units.comoving_coordinates = 0;
   my_units.density_units        = UNIT_D;
   my_units.length_units         = UNIT_L;
   my_units.time_units           = UNIT_T;
   my_units.velocity_units       = UNIT_V;
   my_units.a_units              = 1.0;
   my_units.a_value              = 1.0;
#  endif


// set the default chemsitry
// --> note that "my_data" will be attached to the Grackle internal pointer "grackle_data"
//     after calling set_default_chemistry_parameters()
// --> we must NOT deallocate "my_data" during the simulation
// --> currently it's deallocated by Grackle_End()
  chemistry_data *my_data = new chemistry_data;

  if ( set_default_chemistry_parameters(my_data) == 0 )
    Aux_Error( ERROR_INFO, "set_default_chemistry_parameters() failed !!\n" );


// set chemistry by accessing "grackle_data"
   grackle_data->use_grackle            = ( GRACKLE_MODE == GRACKLE_MODE_ORI );
   grackle_data->with_radiative_cooling = GRACKLE_COOLING;
   grackle_data->primordial_chemistry   = GRACKLE_PRIMORDIAL;
   grackle_data->metal_cooling          = GRACKLE_METAL;
   grackle_data->UVbackground           = GRACKLE_UV;
   grackle_data->grackle_data_file      = GRACKLE_CLOUDY_TABLE;
#  ifdef OPENMP
// currently we apply the OpenMP parallelization to **multiple patches** instead of **cells inside a patch**
// --> therefore, we should disable OpenMP for Grackle (i.e., no nested parallelization)
// grackle_data->omp_nthreads           = OMP_NTHREAD;
   grackle_data->omp_nthreads           = 1;
#  endif

#  if ( MODEL == HYDRO )
   grackle_data->Gamma                  = GAMMA;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif


// initialize the chemistry object
   if ( initialize_chemistry_data(&my_units) == 0 )
     Aux_Error( ERROR_INFO, "initialize_chemistry_data() failed !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Grackle_Init



#endif // #ifdef SUPPORT_GRACKLE
