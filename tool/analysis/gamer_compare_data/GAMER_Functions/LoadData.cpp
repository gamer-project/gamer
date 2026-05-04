#include "CompareData.h"
#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif


static void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal );
static void Load_Parameter_After_2000( FILE *File, const int FormatVersion, const long HeaderOffset_Makefile,
                                       const long HeaderOffset_Constant, const long HeaderOffset_Parameter,
                                       int *NX0_Tot, bool &WithPar, int &NField, int &NMag, int &NParAttFlt,
                                       int &NParAttInt );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load the input data from the binary file "FileName"
//
// Note        :  1. *Label[MAX_STRING], ParFltData[], and ParIntData[] will be allocated here and must be deallocated manually
//
// Parameter   :  FileName       : Name of the input file
//                amr            : Target AMR_t pointer
//                Format         : 1/2 --> C-binary/HDF5
//                NField         : Number of cell-centered fields stored in the file
//                NMag           : Number of face-centered magnetic components stored in the file
//                NParAttFlt     : Number of particle floating-point attributes stored in the file
//                NParAttInt     : Number of particle integer        attributes stored in the file
//                NPar           : NUmber of particles
//                ParFltData     : Particle floating-point data array (allocated here --> must be dallocated manually later)
//                ParIntData     : Particle integer        data array (allocated here --> must be dallocated manually later)
//                FieldLabel     : Labels of cell-centered fields
//                MagLabel       : Labels of face-centered magnetic components
//                ParAttFltLabel : Labels of particle floating-point attributes
//                ParAttIntLabel : Labels of particle integer        attributes
//
// Return      :  amr, Format, NField, NMag, NParAttFlt, NParAttInt, NPar, ParFltData, ParIntData
//                FieldLabel, MagLabel, ParAttFltLabel, ParAttIntLabel (for HDF5 only)
//-------------------------------------------------------------------------------------------------------
void LoadData( const char *FileName, AMR_t &amr, int &Format, int &NField, int &NMag, int &NParAttFlt, int &NParAttInt,
               long &NPar, real_par **&ParFltData, long_par **&ParIntData, char (*&FieldLabel)[MAX_STRING],
               char (*&MagLabel)[MAX_STRING], char (*&ParAttFltLabel)[MAX_STRING], char (*&ParAttIntLabel)[MAX_STRING] )
{

// load the HDF5 data
#  ifdef SUPPORT_HDF5
   if (  Aux_CheckFileExist(FileName)  &&  H5Fis_hdf5(FileName)  )
   {
      LoadData_HDF5( FileName, amr, Format, NField, NMag, NParAttFlt, NParAttInt, NPar,
                     ParFltData, ParIntData, FieldLabel, MagLabel, ParAttFltLabel, ParAttIntLabel );
      return;
   }
#  endif


   Aux_Message( stdout, "Loading binary data %s ...\n", FileName );


// check if the target file exists
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "rb" );


// record the size of different data types
   const int size_bool   = sizeof( bool   );
   const int size_int    = sizeof( int    );
   const int size_uint   = sizeof( uint   );
   const int size_long   = sizeof( long   );
   const int size_ulong  = sizeof( ulong  );
   const int size_real   = sizeof( real   );
   const int size_double = sizeof( double );

   if ( size_int != size_uint )
      Aux_Error( ERROR_INFO, "sizeof(int) = %d != sizeof(uint) = %d !!\n", size_int, size_uint );

   if ( size_long != size_ulong )
      Aux_Error( ERROR_INFO, "sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );


// a. load the information of data format
// =================================================================================================
   long FormatVersion, CheckCode;
   long HeaderSize_Format, HeaderSize_Makefile, HeaderSize_Constant, HeaderSize_Parameter, HeaderSize_SimuInfo;
   long HeaderOffset_Format, HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter, HeaderOffset_SimuInfo;
   long HeaderSize_Total;


// load the file format
   fread( &FormatVersion, sizeof(long), 1, File );


// load the size of each header (in bytes)
   fread( &HeaderSize_Format,    sizeof(long), 1, File );
   fread( &HeaderSize_Makefile,  sizeof(long), 1, File );
   fread( &HeaderSize_Constant,  sizeof(long), 1, File );
   fread( &HeaderSize_Parameter, sizeof(long), 1, File );
   fread( &HeaderSize_SimuInfo,  sizeof(long), 1, File );
   fread( &CheckCode,            sizeof(long), 1, File );


// determine the offset of each header (in bytes)
   HeaderOffset_Format    = 0;    // it must be zero
   HeaderOffset_Makefile  = HeaderOffset_Format    + HeaderSize_Format;
   HeaderOffset_Constant  = HeaderOffset_Makefile  + HeaderSize_Makefile;
   HeaderOffset_Parameter = HeaderOffset_Constant  + HeaderSize_Constant;
   HeaderOffset_SimuInfo  = HeaderOffset_Parameter + HeaderSize_Parameter;

   HeaderSize_Total       = HeaderOffset_SimuInfo  + HeaderSize_SimuInfo;


// verify the input data format version
   Aux_Message( stdout, "   The format version of the input file = %ld\n", FormatVersion );

   if ( FormatVersion < 2000 )
      Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2000) !!\n" );


// check if the size of different data types are consistent
   int size_bool_restart, size_int_restart, size_long_restart, size_real_restart, size_double_restart;

   fread( &size_bool_restart,   sizeof(int), 1, File );
   fread( &size_int_restart,    sizeof(int), 1, File );
   fread( &size_long_restart,   sizeof(int), 1, File );
   fread( &size_real_restart,   sizeof(int), 1, File );
   fread( &size_double_restart, sizeof(int), 1, File );

   if ( size_bool_restart != size_bool )
      Aux_Error( ERROR_INFO, "sizeof(bool) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_bool_restart, size_bool );

   if ( size_int_restart != size_int )
      Aux_Error( ERROR_INFO, "sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_int_restart, size_int );

   if ( size_long_restart != size_long )
      Aux_Error( ERROR_INFO, "sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_long_restart, size_long );

   if ( size_real_restart != size_real )
      Aux_Error( ERROR_INFO, "sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_real_restart, size_real );

   if ( size_double_restart != size_double )
      Aux_Error( ERROR_INFO, "sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_double_restart, size_double );


// b. load all simulation parameters
// =================================================================================================
   bool WithPar;
   Load_Parameter_After_2000( File, FormatVersion, HeaderOffset_Makefile, HeaderOffset_Constant,
                              HeaderOffset_Parameter, amr.nx0_tot, WithPar, NField, NMag, NParAttFlt,
                              NParAttInt );


// c. load the simulation information
// =================================================================================================
   Aux_Message( stdout, "   Loading simulation information ...\n" );

// verify the check code
   long checkcode;

   fseek( File, HeaderOffset_SimuInfo, SEEK_SET );

   fread( &checkcode, sizeof(long), 1, File );

   if ( checkcode != CheckCode )
      Aux_Error( ERROR_INFO, "incorrect check code in the RESTART file (input %ld <-> expeect %ld) !!\n",
                 checkcode, CheckCode );


// load information necessary for restart
   long   AdvanceCounter[NLEVEL];
   int    NPatchTotal[NLEVEL], NDataPatch_Total[NLEVEL], DumpID;
   long   Step, FileOffset_Particle;
   double Time[NLEVEL];

   fread( &DumpID,          sizeof(int),         1, File );
   fread( Time,             sizeof(double), NLEVEL, File );
   fread( &Step,            sizeof(long),        1, File );
   fread( NPatchTotal,      sizeof(int),    NLEVEL, File );
   fread( NDataPatch_Total, sizeof(int),    NLEVEL, File );
   fread( AdvanceCounter,   sizeof(long),   NLEVEL, File );
   fseek( File, sizeof(double), SEEK_CUR );

   if ( WithPar )
   {
      fread( &NPar,                sizeof(long), 1, File );
      fread( &FileOffset_Particle, sizeof(long), 1, File );
   }
   else
      NPar = 0;


// verify the size of the RESTART file
   Aux_Message( stdout, "      Verifying the size of the RESTART file ...\n" );

   long ExpectSize, InputSize, PatchDataSize, DataSize[NLEVEL];

   PatchDataSize  = CUBE(PS1)*NField*sizeof(real);
   PatchDataSize += PS1P1*SQR(PS1)*NMag*sizeof(real);
   ExpectSize     = HeaderSize_Total;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      DataSize[lv]  = 0;
      DataSize[lv] += (long)NPatchTotal[lv]*4*sizeof(int);        // 4 = corner(3) + son(1)
      DataSize[lv] += (long)NDataPatch_Total[lv]*PatchDataSize;
      ExpectSize   += (long)DataSize[lv];
   }

   if ( WithPar )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
//       2 = NPar + starting particle index stored in each leaf patch
         const long ParInfoSize = (long)NDataPatch_Total[lv]*2*sizeof(long);

         DataSize[lv] += ParInfoSize;
         ExpectSize   += ParInfoSize;
      }

      ExpectSize += (long)NParAttFlt*NPar*sizeof(real_par);
      ExpectSize += (long)NParAttInt*NPar*sizeof(long_par);
   }

   fseek( File, 0, SEEK_END );
   InputSize = ftell( File );

   if ( InputSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
                 FileName, InputSize, ExpectSize );

#  if ( (defined FLOAT8 && !defined FLOAT8_PAR) || (!defined FLOAT8 && defined FLOAT8_PAR) )
   if ( WithPar  &&  FormatVersion < 2300 )   Aux_Error( ERROR_INFO, "Must adopt FLOAT8_PAR=FLOAT8 for C-binary snapshots for versions < 2300 !!\n" );
#  endif

   Aux_Message( stdout, "      Verifying the size of the RESTART file ... passed\n" );
   Aux_Message( stdout, "   Loading simulation information ... done\n" );


// d. load the simulation data
// =================================================================================================
   int  LoadCorner[3], LoadSon, PID;

   fseek( File, HeaderSize_Total, SEEK_SET );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Aux_Message( stdout, "   Loading grid data at level %2d ... ", lv );

      for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
      {
//       d1. load the patch information
         fread(  LoadCorner, sizeof(int), 3, File );
         fread( &LoadSon,    sizeof(int), 1, File );

//       d2. create the patch and load the physical data only if it is a leaf patch
         if ( LoadSon == -1 )
         {
//          skip particle info
            if ( WithPar )    fseek( File, 2*sizeof(long), SEEK_CUR );

            PID = amr.num[lv];

            amr.pnew( lv, LoadCorner[0], LoadCorner[1], LoadCorner[2], -1, true, NField, NMag );
            amr.patch[lv][PID]->son = -1;    // set the SonPID as -1 to indicate that it's a leaf patch

//          cell-centered fields
            fread( amr.patch[lv][PID]->field, sizeof(real), CUBE(PS1)*NField,    File );

//          face-centered B field
            fread( amr.patch[lv][PID]->mag,   sizeof(real), PS1P1*SQR(PS1)*NMag, File );
         }
      } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)

      Aux_Message( stdout, "done\n" );
   } // for (int lv=0; lv<NLEVEL; lv++)

   fclose( File );


// load particles
   if ( NPar > 0 )
   {
      const long ParFltDataSize1v = NPar*sizeof(real_par);
      const long ParIntDataSize1v = NPar*sizeof(long_par);

      Aux_AllocateArray2D<real_par>( ParFltData, NParAttFlt, NPar );
      Aux_AllocateArray2D<long_par>( ParIntData, NParAttInt, NPar );

      File = fopen( FileName, "rb" );

      for (int v=0; v<NParAttFlt; v++)
      {
         fseek( File, FileOffset_Particle + v*ParFltDataSize1v, SEEK_SET );
         fread( ParFltData[v], sizeof(real_par), NPar, File );
      }

      for (int v=0; v<NParAttInt; v++)
      {
         fseek( File, FileOffset_Particle + NParAttFlt*ParFltDataSize1v + v*ParIntDataSize1v, SEEK_SET );
         fread( ParIntData[v], sizeof(long_par), NPar, File );
      }

      fclose( File );
   } // if ( NPar > 0 )


// e. record parameters
// =================================================================================================
   Format = 1;

   Aux_Message( stdout, "\n" );
   Aux_Message( stdout, "   Data information\n" );
   Aux_Message( stdout, "   ==============================\n" );
   Aux_Message( stdout, "   DumpID      = %d\n",  DumpID      );
   Aux_Message( stdout, "   Step        = %ld\n", Step        );
   Aux_Message( stdout, "   Time        = %lf\n", Time[0]     );
   Aux_Message( stdout, "   Format      = %d\n",  Format      );
   Aux_Message( stdout, "   WithPar     = %d\n",  WithPar     );
   Aux_Message( stdout, "   NField      = %d\n",  NField      );
   Aux_Message( stdout, "   NMag        = %d\n",  NMag        );
   if ( WithPar ) {
   Aux_Message( stdout, "   NParAttFlt  = %d\n",  NParAttFlt  );
   Aux_Message( stdout, "   NParAttInt  = %d\n",  NParAttInt  );
   Aux_Message( stdout, "   NPar        = %ld\n", NPar        ); }
   for (int lv=0; lv<NLEVEL; lv++)
   Aux_Message( stdout, "   NPatch[%2d] = %d\n",  lv, amr.num[lv] );
   Aux_Message( stdout, "   ==============================\n" );
   Aux_Message( stdout, "\n" );


   Aux_Message( stdout, "Loading binary data %s ... done\n", FileName );

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_After_2000
// Description :  Load all simulation parameters from the RESTART file with format version >= 2000
//
// Note        :  All floating-point variables are declared as double after version 2000
//
// Parameter   :  File           : RESTART file pointer
//                FormatVersion  : Format version of the RESTART file
//                HeaderOffset_X : Offsets of different headers
//                NX0_Tot        : Total number of base-level cells along each direction
//                WithPar        : Whether or not the RESTART file stores particles
//                NField         : Number of cell-centered fields stored in the file
//                NMag           : Number of face-centered magnetic components stored in the file
//                NParAttFlt     : Number of particle floating-point attributes stored in the file
//                NParAttInt     : Number of particle integer        attributes stored in the file
//
// Return      :  NX0_Tot, WithPar, NField, NMag, NParAttFlt, NParAttInt
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, const long HeaderOffset_Makefile,
                                const long HeaderOffset_Constant, const long HeaderOffset_Parameter,
                                int *NX0_Tot, bool &WithPar, int &NField, int &NMag, int &NParAttFlt,
                                int &NParAttInt )
{

   Aux_Message( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options and parameters defined in the Makefile
// =================================================================================================
   bool gravity, individual_timestep, comoving, gpu, gamer_optimization, gamer_debug, timing, timing_solver;
   bool intel, float8, serial, overlap_mpi, openmp;
   int  model, pot_scheme, flu_scheme, lr_scheme, rsolver, load_balance, nlevel, max_patch, gpu_arch, ncomp_passive;
   bool store_pot_ghost, unsplit_gravity, particle;
   bool conserve_mass, laplacian_4th, self_interaction, laohu, support_hdf5, mhd;

   fseek( File, HeaderOffset_Makefile, SEEK_SET );

   fread( &model,                      sizeof(int),                     1,             File );
   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &pot_scheme,                 sizeof(int),                     1,             File );
   fread( &individual_timestep,        sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &flu_scheme,                 sizeof(int),                     1,             File );
   fread( &lr_scheme,                  sizeof(int),                     1,             File );
   fread( &rsolver,                    sizeof(int),                     1,             File );
   fread( &gpu,                        sizeof(bool),                    1,             File );
   fread( &gamer_optimization,         sizeof(bool),                    1,             File );
   fread( &gamer_debug,                sizeof(bool),                    1,             File );
   fread( &timing,                     sizeof(bool),                    1,             File );
   fread( &timing_solver,              sizeof(bool),                    1,             File );
   fread( &intel,                      sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );
   fread( &serial,                     sizeof(bool),                    1,             File );
   fread( &load_balance,               sizeof(int),                     1,             File );
   fread( &overlap_mpi,                sizeof(bool),                    1,             File );
   fread( &openmp,                     sizeof(bool),                    1,             File );
   fread( &gpu_arch,                   sizeof(int),                     1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );
   fread( &store_pot_ghost,            sizeof(bool),                    1,             File );
   fread( &unsplit_gravity,            sizeof(bool),                    1,             File );
   fread( &particle,                   sizeof(bool),                    1,             File );
   fread( &ncomp_passive,              sizeof(int),                     1,             File );
   fread( &conserve_mass,              sizeof(bool),                    1,             File );
   fread( &laplacian_4th,              sizeof(bool),                    1,             File );
   fread( &self_interaction,           sizeof(bool),                    1,             File );
   fread( &laohu,                      sizeof(bool),                    1,             File );
   fread( &support_hdf5,               sizeof(bool),                    1,             File );
   fread( &mhd,                        sizeof(bool),                    1,             File );

   if ( FormatVersion < 2100 )   particle      = false;
   if ( FormatVersion < 2120 )   ncomp_passive = 0;
   if ( FormatVersion < 2210 )   mhd           = false;


// b. load the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
// =================================================================================================
   bool   enforce_positive, char_reconstruction, hll_no_ref_state, hll_include_all_waves, waf_dissipate;
   bool   use_psolver_10to14;
   int    ncomp_fluid, patch_size, flu_ghost_size, pot_ghost_size, gra_ghost_size, check_intermediate;
   int    flu_block_size_x, flu_block_size_y, pot_block_size_x, pot_block_size_z, gra_block_size_z;
   int    par_nvar_flt, par_npassive_flt;
   int    par_nvar_int, par_npassive_int;
   double min_value, max_error;

   fseek( File, HeaderOffset_Constant, SEEK_SET );

   fread( &ncomp_fluid,                sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &min_value,                  sizeof(double),                  1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );
   fread( &enforce_positive,           sizeof(bool),                    1,             File );
   fread( &char_reconstruction,        sizeof(bool),                    1,             File );
   fread( &check_intermediate,         sizeof(int),                     1,             File );
   fread( &hll_no_ref_state,           sizeof(bool),                    1,             File );
   fread( &hll_include_all_waves,      sizeof(bool),                    1,             File );
   fread( &waf_dissipate,              sizeof(bool),                    1,             File );
   fread( &max_error,                  sizeof(double),                  1,             File );
   fread( &flu_block_size_x,           sizeof(int),                     1,             File );
   fread( &flu_block_size_y,           sizeof(int),                     1,             File );
   fread( &use_psolver_10to14,         sizeof(bool),                    1,             File );
   fread( &pot_block_size_x,           sizeof(int),                     1,             File );
   fread( &pot_block_size_z,           sizeof(int),                     1,             File );
   fread( &gra_block_size_z,           sizeof(int),                     1,             File );
   if ( particle ) {
   fread( &par_nvar_flt,               sizeof(int),                     1,             File );
   fread( &par_npassive_flt,           sizeof(int),                     1,             File );
   fread( &par_nvar_int,               sizeof(int),                     1,             File );
   fread( &par_npassive_int,           sizeof(int),                     1,             File ); }



// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   opt__adaptive_dt, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_test_error, opt__output_base, opt__output_pot;
   bool   opt__output_baseps, opt__timing_balance, opt__int_phase, opt__corr_unphy, opt__unit;
   bool   opt__output_cc_mag;
   int    nx0_tot[3], mpi_nrank, mpi_nrank_x[3], omp_nthread, regrid_count, opt__output_par_dens;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter, flu_gpu_npgroup, gpu_nstream;
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step, opt__corr_unphy_scheme;
   long   end_step;
   double lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, elbdm_planck_const, newton_g, sor_omega;
   double mg_tolerated_error, output_part_x, output_part_y, output_part_z, molecular_weight;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt, hubble0;
   double unit_l, unit_m, unit_t, unit_v, unit_d, unit_e, unit_p, unit_b;

   fseek( File, HeaderOffset_Parameter, SEEK_SET );

   fread( &box_size,                   sizeof(double),                  1,             File );
   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &mpi_nrank,                  sizeof(int),                     1,             File );
   fread(  mpi_nrank_x,                sizeof(int),                     3,             File );
   fread( &omp_nthread,                sizeof(int),                     1,             File );
   fread( &end_t,                      sizeof(double),                  1,             File );
   fread( &end_step,                   sizeof(long),                    1,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__phase,                  sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &opt__adaptive_dt,           sizeof(bool),                    1,             File );
   fread( &opt__dt_user,               sizeof(bool),                    1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &max_level,                  sizeof(int),                     1,             File );
   fread( &opt__flag_rho,              sizeof(bool),                    1,             File );
   fread( &opt__flag_rho_gradient,     sizeof(bool),                    1,             File );
   fread( &opt__flag_pres_gradient,    sizeof(bool),                    1,             File );
   fread( &opt__flag_engy_density,     sizeof(bool),                    1,             File );
   fread( &opt__flag_user,             sizeof(bool),                    1,             File );
   fread( &lb_wli_max,                 sizeof(double),                  1,             File );
   fread( &gamma,                      sizeof(double),                  1,             File );
   fread( &minmod_coeff,               sizeof(double),                  1,             File );
   fread( &ep_coeff,                   sizeof(double),                  1,             File );
   fread( &opt__lr_limiter,            sizeof(int),                     1,             File );
   fread( &opt__waf_limiter,           sizeof(int),                     1,             File );
   fread( &elbdm_mass,                 sizeof(double),                  1,             File );
   fread( &elbdm_planck_const,         sizeof(double),                  1,             File );
   fread( &flu_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &gpu_nstream,                sizeof(int),                     1,             File );
   fread( &opt__fixup_flux,            sizeof(bool),                    1,             File );
   fread( &opt__fixup_restrict,        sizeof(bool),                    1,             File );
   fread( &opt__overlap_mpi,           sizeof(bool),                    1,             File );
   fread( &newton_g,                   sizeof(double),                  1,             File );
   fread( &sor_omega,                  sizeof(double),                  1,             File );
   fread( &sor_max_iter,               sizeof(int),                     1,             File );
   fread( &sor_min_iter,               sizeof(int),                     1,             File );
   fread( &mg_max_iter,                sizeof(int),                     1,             File );
   fread( &mg_npre_smooth,             sizeof(int),                     1,             File );
   fread( &mg_npost_smooth,            sizeof(int),                     1,             File );
   fread( &mg_tolerated_error,         sizeof(double),                  1,             File );
   fread( &pot_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__int_time,              sizeof(bool),                    1,             File );
   fread( &opt__int_phase,             sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &opt__output_part,           sizeof(int),                     1,             File );
   fread( &opt__output_test_error,     sizeof(bool),                    1,             File );
   fread( &opt__output_base,           sizeof(bool),                    1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_mode,           sizeof(int),                     1,             File );
   fread( &output_step,                sizeof(int),                     1,             File );
   fread( &output_dt,                  sizeof(double),                  1,             File );
   fread( &output_part_x,              sizeof(double),                  1,             File );
   fread( &output_part_y,              sizeof(double),                  1,             File );
   fread( &output_part_z,              sizeof(double),                  1,             File );
   fread( &opt__timing_balance,        sizeof(bool),                    1,             File );
   fread( &opt__output_baseps,         sizeof(bool),                    1,             File );
   fread( &opt__corr_unphy,            sizeof(bool),                    1,             File );
   fread( &opt__corr_unphy_scheme,     sizeof(int),                     1,             File );
   fread( &opt__output_par_dens,       sizeof(int),                     1,             File );
   fread( &hubble0,                    sizeof(double),                  1,             File );
   fread( &opt__unit,                  sizeof(bool),                    1,             File );
   fread( &unit_l,                     sizeof(double),                  1,             File );
   fread( &unit_m,                     sizeof(double),                  1,             File );
   fread( &unit_t,                     sizeof(double),                  1,             File );
   fread( &unit_v,                     sizeof(double),                  1,             File );
   fread( &unit_d,                     sizeof(double),                  1,             File );
   fread( &unit_e,                     sizeof(double),                  1,             File );
   fread( &unit_p,                     sizeof(double),                  1,             File );
   fread( &molecular_weight,           sizeof(double),                  1,             File );
   fread( &opt__output_cc_mag,         sizeof(bool),                    1,             File );
   fread( &unit_b,                     sizeof(double),                  1,             File );

   if ( !particle )  opt__output_par_dens = 0;
   if ( !mhd      )  opt__output_cc_mag   = 0;


   Aux_Message( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   Aux_Message( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
      Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
#  else
   if (  float8 )
      Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
#  endif

   CompareVar( "MODEL",      model,      MODEL,      Fatal );
   CompareVar( "NLEVEL",     nlevel,     NLEVEL,     Fatal );
   CompareVar( "PATCH_SIZE", patch_size, PATCH_SIZE, Fatal );


   Aux_Message( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   for (int d=0; d<3; d++)    NX0_Tot[d] = nx0_tot[d];

   WithPar = particle;
   NField  = ncomp_fluid + ncomp_passive;
   NMag    = ( mhd ) ? 3 : 0;

   if ( opt__output_pot )        NField ++;
   if ( opt__output_par_dens )   NField ++;
   if ( opt__output_cc_mag )     NField += 3;

   if ( WithPar )
   {
      if ( FormatVersion >= 2300 )
      {
         NParAttFlt = par_nvar_flt;           // after version 2300, par_nvar_flt = PAR_NATT_FLT_STORED
         NParAttInt = par_nvar_int;           // after version 2300, par_nvar_int = PAR_NATT_INT_STORED
      }
      else if ( FormatVersion > 2131 )
      {
         NParAttFlt = par_nvar_flt;           // after version 2131, par_nvar_flt = PAR_NATT_STORED
         NParAttInt = 0;
      }
      else
      {
         NParAttFlt = 7 + par_npassive_flt;   // mass, position x/y/z, velocity x/y/z, and passive variables
         NParAttInt = 0;
      }
   }
   else
   {
         NParAttFlt = 0;
         NParAttInt = 0;
   }

} // FUNCTION : Load_Parameter_After_2000



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareVar
// Description :  Compare the input variables
//
// Note        :  This function is overloaded to work with different data types
//
// Parameter   :  VarName     : Name of the targeted variable
//                RestartVar  : Variable loaded from the RESTART file
//                RuntimeVar  : Variable loaded from the Input__Parameter
//                Fatal       : Whether or not the difference between RestartVar and RuntimeVar is fatal
//                              --> true  : terminate the program if the input variables are different
//                                  false : display warning message if the input variables are different
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%d) != runtime (%d) !!\n",
                    VarName, RestartVar, RuntimeVar );
      else
         Aux_Message( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n",
                      VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (int)
