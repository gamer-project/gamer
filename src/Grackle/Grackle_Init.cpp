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
   if ( !GRACKLE_ACTIVATE )   return;


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

   if (  ( GRACKLE_PRIMORDIAL == GRACKLE_PRI_CHE_CLOUDY || GRACKLE_METAL || GRACKLE_UV )  &&
         !Aux_CheckFileExist(GRACKLE_CLOUDY_TABLE)  )
      Aux_Error( ERROR_INFO, "Grackle data file \"%s\" does not exist !!\n", GRACKLE_CLOUDY_TABLE );


// enable output (only for the root rank)
   grackle_verbose = ( MPI_Rank == 0 ) ? GRACKLE_VERBOSE : 0;


// units in cgs
// --> Che_Units is declared as a global variable since all Grackle solvers require that as well
#  ifdef COMOVING
   Che_Units.comoving_coordinates = 1;
   Che_Units.density_units        = NULL_REAL;  // not sure how to set the units in the comoving coordinates yet...
   Che_Units.length_units         = NULL_REAL;  // --> see http://grackle.readthedocs.io/en/latest/Integration.html
   Che_Units.time_units           = NULL_REAL;
   Che_Units.velocity_units       = NULL_REAL;
   Che_Units.a_units              = 1.0;
   Che_Units.a_value              = Time[0];

#  else
   Che_Units.comoving_coordinates = 0;
   Che_Units.density_units        = UNIT_D;
   Che_Units.length_units         = UNIT_L;
   Che_Units.time_units           = UNIT_T;
   Che_Units.velocity_units       = UNIT_V;
   Che_Units.a_units              = 1.0;
   Che_Units.a_value              = 1.0;
#  endif


// set the default chemsitry
// --> note that "Che_Data" will be attached to the Grackle internal pointer "grackle_data"
//     after calling set_default_chemistry_parameters()
// --> we must NOT deallocate "Che_Data" during the simulation
// --> currently it's deallocated by Grackle_End()
  chemistry_data *Che_Data = new chemistry_data;

  if ( set_default_chemistry_parameters(Che_Data) == 0 )
    Aux_Error( ERROR_INFO, "set_default_chemistry_parameters() failed !!\n" );


// set chemistry by accessing "grackle_data"
   grackle_data->use_grackle                = GRACKLE_ACTIVATE;
   grackle_data->with_radiative_cooling     = GRACKLE_COOLING;
   grackle_data->primordial_chemistry       = GRACKLE_PRIMORDIAL;
   grackle_data->metal_cooling              = GRACKLE_METAL;
   grackle_data->UVbackground               = GRACKLE_UV;
   grackle_data->cmb_temperature_floor      = GRACKLE_CMB_FLOOR;
   grackle_data->photoelectric_heating      = GRACKLE_PE_HEATING;
   grackle_data->photoelectric_heating_rate = GRACKLE_PE_HEATING_RATE;
   grackle_data->grackle_data_file          = GRACKLE_CLOUDY_TABLE;

#  ifdef OPENMP
// currently we adopt the OpenMP implementation in Grackle directly, which applies the parallelization to
// **different cells inside a patch group** instead of **different patch groups**
// --> this approach is found to be more efficient
// --> therefore, we should enable OpenMP for Grackle and disable OpenMP in CPU_GrackleSolver()
//     to avoid the nested parallelization
   grackle_data->omp_nthreads               = OMP_NTHREAD;
#  endif

#  if ( MODEL == HYDRO )
   grackle_data->Gamma                      = GAMMA;
#  endif


// initialize the chemistry object
   if ( initialize_chemistry_data(&Che_Units) == 0 )
     Aux_Error( ERROR_INFO, "initialize_chemistry_data() failed !!\n" );


// initialize the "grackle_field_data" object of Grackle
   Grackle_Init_FieldData();


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Grackle_Init



#endif // #ifdef SUPPORT_GRACKLE
