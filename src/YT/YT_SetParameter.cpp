#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_SetParameter
// Description :  Set YT-specific parameters for the inline analysis
//
// Note        :  1. This function must be called in advance **every time** we invoke the inline analysis
//                2. Invoked by YT_Inline().
//                3. Set up num_species, species_list for supporting PARTICLE.
//
// Parameter   :  NPatchAllLv : Total number of patches at all levels
//                NField      : Total number of fields
//                NPatchLocal : Number of local patches at all levels
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalAllLv )
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. prepare the simulation information for libyt
   yt_param_yt param_yt;

   param_yt.frontend                = "gamer";           // simulation frontend
   if ( strcmp(YT_FIG_BASENAME, "") != 0 )
       param_yt.fig_basename = YT_FIG_BASENAME;          // figure base name, use default if not set (default=Fig%09d)

   param_yt.length_unit             = UNIT_L;            // units are in cgs
   param_yt.mass_unit               = UNIT_M;
   param_yt.time_unit               = UNIT_T;

#  ifdef MHD
   param_yt.magnetic_unit           = UNIT_B;
#  endif

   param_yt.current_time            = Time[0];
   param_yt.dimensionality          = 3;
   param_yt.refine_by               = 2;
   param_yt.num_fields              = NField;

#  ifdef LIBYT_USE_PATCH_GROUP
   if ( NPatchAllLv % 8 != 0 || NPatchLocalAllLv % 8 != 0 ) Aux_Error( ERROR_INFO, "Using patch group in libyt failed !!\n" );
   param_yt.num_grids               = NPatchAllLv / 8;
   param_yt.num_grids_local         = NPatchLocalAllLv / 8;
#  else
   param_yt.num_grids               = NPatchAllLv;
   param_yt.num_grids_local         = NPatchLocalAllLv;
#  endif

#  ifdef PARTICLE
   yt_species *species_list         = new yt_species [1];
   species_list[0].species_name     = "io";
   species_list[0].num_attr         = PAR_NATT_TOTAL;

   param_yt.num_species             = 1;
   param_yt.species_list            = species_list;
#  endif

   for (int d=0; d<3; d++)
   {
      param_yt.domain_dimensions[d] = NX0_TOT[d];
      param_yt.domain_left_edge [d] = 0.0;
      param_yt.domain_right_edge[d] = amr->BoxSize[d];
      param_yt.periodicity      [d] = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC ) ? 1 : 0;
   }

#  ifdef COMOVING
   param_yt.cosmological_simulation = 1;
   param_yt.current_redshift        = 1.0/Time[0] - 1.0;
   param_yt.omega_lambda            = 1.0 - OMEGA_M0;
   param_yt.omega_matter            = OMEGA_M0;
   param_yt.hubble_constant         = HUBBLE0;
#  else
   param_yt.cosmological_simulation = 0;
   param_yt.current_redshift        = 0.0;
   param_yt.omega_lambda            = 0.0;
   param_yt.omega_matter            = 0.0;
   param_yt.hubble_constant         = 0.0;
#  endif


// 2. transfer simulation information to libyt
   if ( yt_set_parameter( &param_yt ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_set_parameter() failed !!\n" );

// 2-1. free no longer used resource
#  ifdef PARTICLE
   delete [] species_list;
#  endif

// 3. set code specific parameter
#  ifdef MHD
   const int mhd = 1;
#  else
   const int mhd = 0;
#  endif
   if (yt_add_user_parameter_int("mhd", 1, &mhd) != YT_SUCCESS)  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add mhd failed !!\n" );

#  if ( MODEL == HYDRO )
   const double gamma = (double) GAMMA;
   const double mu = (double) MOLECULAR_WEIGHT;
#  ifdef SRHD
   const int srhd = 1;
#  else
   const int srhd = 0;
#  endif
   if (yt_add_user_parameter_double("gamma", 1, &gamma) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add GAMMA failed !!\n" );
   if (yt_add_user_parameter_double("mu", 1, &mu) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add MOLECULAR_WEIGHT failed !!\n" );
   if (yt_add_user_parameter_int("srhd", 1, &srhd) != YT_SUCCESS ) Aux_Error( ERROR_INFO, "yt_add_user_parameter() add srhd failed !!\n" );
#  endif

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );


} // FUNCTION : YT_SetParameter



#endif // #ifdef SUPPORT_LIBYT
