#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_FlagCriteria
// Description :  Load the flag criteria from several input files named "Input__Flag_XXX"
//
// Note        :  The loaded data will be used as the refinement thresholds at each refinement level
//-------------------------------------------------------------------------------------------------------
void Init_Load_FlagCriteria()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  if ( MODEL != HYDRO )
   const bool OPT__FLAG_PRES_GRADIENT = false;
   double *FlagTable_PresGradient     = NULL;

   const bool OPT__FLAG_VORTICITY     = false;
   double *FlagTable_Vorticity        = NULL;

   const bool OPT__FLAG_JEANS         = false;
   double *FlagTable_Jeans            = NULL;
#  endif

#  ifndef MHD
   const bool OPT__FLAG_CURRENT       = false;
   double *FlagTable_Current          = NULL;
#  endif

#  ifndef SRHD
   const bool OPT__FLAG_LRTZ_GRADIENT = false;
   double *FlagTable_LrtzGradient     = NULL;
#  endif
#  ifndef COSMIC_RAY
   const bool OPT__FLAG_LOHNER_CRAY   = false;
   const bool OPT__FLAG_CRAY          = false;
   double *FlagTable_CRay             = NULL;
#  endif

#  if ( MODEL != ELBDM )
   const bool OPT__FLAG_ENGY_DENSITY  = false;
   double FlagTable_EngyDensity [NLEVEL-1][2];

   const bool OPT__FLAG_SPECTRAL      = false;
   double FlagTable_Spectral    [NLEVEL-1][2];
#  endif

#  if ( ELBDM_SCHEME != ELBDM_HYBRID )
   const bool OPT__FLAG_INTERFERENCE  = false;
   double FlagTable_Interference[NLEVEL-1][4];
#  endif

#  ifndef PARTICLE
   const bool OPT__FLAG_NPAR_PATCH    = false;
   int *FlagTable_NParPatch           = NULL;

   const bool OPT__FLAG_NPAR_CELL     = false;
   int *FlagTable_NParCell            = NULL;

   const bool OPT__FLAG_PAR_MASS_CELL = false;
   double *FlagTable_ParMassCell      = NULL;
#  endif

#  if   ( MODEL == HYDRO )
   const bool OPT__FLAG_LOHNER = ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES ||
                                   OPT__FLAG_LOHNER_TEMP || OPT__FLAG_LOHNER_ENTR || OPT__FLAG_LOHNER_CRAY );
#  elif ( MODEL == ELBDM )
   const bool OPT__FLAG_LOHNER = OPT__FLAG_LOHNER_DENS;
#  else
#  error : unsupported MODEL !!
#  endif

   const int  NFlagMode         = 18;
   const bool Flag[NFlagMode]   = { OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_PRES_GRADIENT,
                                    OPT__FLAG_ENGY_DENSITY, OPT__FLAG_LOHNER, OPT__FLAG_USER,
                                    (bool)OPT__FLAG_NPAR_PATCH, OPT__FLAG_NPAR_CELL, OPT__FLAG_PAR_MASS_CELL,
                                    OPT__FLAG_VORTICITY, OPT__FLAG_JEANS, OPT__FLAG_CURRENT,
                                    OPT__FLAG_CRAY, OPT__FLAG_LRTZ_GRADIENT, OPT__FLAG_INTERFERENCE,
                                    OPT__FLAG_SPECTRAL, OPT__FLAG_ANGULAR, OPT__FLAG_RADIAL };
   const char ModeName[][100]   = { "OPT__FLAG_RHO", "OPT__FLAG_RHO_GRADIENT", "OPT__FLAG_PRES_GRADIENT",
                                    "OPT__FLAG_ENGY_DENSITY", "OPT__FLAG_LOHNER", "OPT__FLAG_USER",
                                    "OPT__FLAG_NPAR_PATCH", "OPT__FLAG_NPAR_CELL", "OPT__FLAG_PAR_MASS_CELL",
                                    "OPT__FLAG_VORTICITY", "OPT__FLAG_JEANS", "OPT__FLAG_CURRENT",
                                    "OPT__FLAG_CRAY", "OPT__FLAG_LRTZ_GRADIENT", "OPT__FLAG_INTERFERENCE",
                                    "OPT__FLAG_SPECTRAL", "OPT__FLAG_ANGULAR", "OPT__FLAG_RADIAL" };
   const char FileName[][100]   = { "Input__Flag_Rho", "Input__Flag_RhoGradient", "Input__Flag_PresGradient",
                                    "Input__Flag_EngyDensity", "Input__Flag_Lohner", "Input__Flag_User",
                                    "Input__Flag_NParPatch", "Input__Flag_NParCell", "Input__Flag_ParMassCell",
                                    "Input__Flag_Vorticity", "Input__Flag_Jeans", "Input__Flag_Current",
                                    "Input__Flag_CRay", "Input__Flag_LrtzGradient", "Input__Flag_Interference",
                                    "Input__Flag_Spectral", "Input__Flag_AngularResolution",
                                    "Input__Flag_RadialResolution"};
   double *FlagTable[NFlagMode] = { FlagTable_Rho, FlagTable_RhoGradient, FlagTable_PresGradient,
                                    NULL, NULL, NULL,
                                    NULL, NULL, FlagTable_ParMassCell,
                                    FlagTable_Vorticity, FlagTable_Jeans, FlagTable_Current,
                                    FlagTable_CRay, FlagTable_LrtzGradient, NULL, NULL, NULL,
                                    FlagTable_Radial };

   FILE *File;
   char *input_line = NULL, TargetName[100];
   size_t len = 0;
   int Trash, n;


// initialize the flag tables as -1
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      FlagTable_Rho         [lv]    = -1.0;
      FlagTable_RhoGradient [lv]    = -1.0;

      for (int t=0; t<5; t++)
      FlagTable_Lohner      [lv][t] = -1.0;

      for (int t=0; t<3; t++)
      FlagTable_Angular     [lv][t] = -1.0;

      FlagTable_Radial      [lv]    = -1.0;

      FlagTable_User        [lv]    = (double *) malloc( OPT__FLAG_USER_NUM*sizeof(double) );
      for (int t=0; t<OPT__FLAG_USER_NUM; t++)
      FlagTable_User        [lv][t] = -1.0;

#     if   ( MODEL == HYDRO )
      FlagTable_PresGradient[lv]    = -1.0;
      FlagTable_Vorticity   [lv]    = -1.0;
      FlagTable_Jeans       [lv]    = -1.0;
#     ifdef MHD
      FlagTable_Current     [lv]    = -1.0;
#     endif
#     ifdef COSMIC_RAY
      FlagTable_CRay        [lv]    = -1.0;
#     endif
#     ifdef SRHD
      FlagTable_LrtzGradient[lv]    = -1.0;
#     endif

#     elif ( MODEL == ELBDM )
      for (int t=0; t<2; t++) {
      FlagTable_EngyDensity [lv][t] = -1.0;
      }

      for (int t=0; t<2; t++) {
      FlagTable_Spectral    [lv][t] = -1.0;
      }

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      for (int t=0; t<4; t++) {
      FlagTable_Interference[lv][t] = -1.0;
      }
#     endif
#     endif // MODEL

#     ifdef PARTICLE
      FlagTable_NParPatch   [lv]    = -1;
      FlagTable_NParCell    [lv]    = -1;
      FlagTable_ParMassCell [lv]    = -1.0;
#     endif
   } // for (int lv=0; lv<NLEVEL-1; lv++)


// nothing to do if there is no refinement level
   if ( MAX_LEVEL > 0 )
   for (int FlagMode=0; FlagMode<NFlagMode; FlagMode++)
   {
      if ( !Flag[FlagMode] )   continue;

//    open the target file
      strcpy( TargetName, FileName[FlagMode] );

      if ( !Aux_CheckFileExist(TargetName) )
      {
         Aux_Message( stderr, "\n" );
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist for the mode \"%s\" !!\n",
                    TargetName, ModeName[FlagMode] );
      }

      File = fopen( TargetName, "r" );

//    skip the header
      getline( &input_line, &len, File );

//    begin to read
      for (int lv=0; lv<MAX_LEVEL; lv++)
      {
         n = getline( &input_line, &len, File );

//       check
         if ( n <= 1 )
         {
            Aux_Message( stderr, "\n" );
            Aux_Error( ERROR_INFO, "incorrect reading at level %d of the file <%s> !!\n",
                       lv, TargetName );
         }

//       OPT__FLAG_ENGY_DENSITY and OPT__FLAG_LOHNER have two and five columns to be loaded, respectively
         if      ( FlagMode == 3 )  sscanf( input_line, "%d%lf%lf", &Trash, &FlagTable_EngyDensity[lv][0],
                                                                            &FlagTable_EngyDensity[lv][1] );
         else if ( FlagMode == 4 )  sscanf( input_line, "%d%lf%lf%lf%lf%lf", &Trash, &FlagTable_Lohner[lv][0],
                                                                                     &FlagTable_Lohner[lv][1],
                                                                                     &FlagTable_Lohner[lv][2],
                                                                                     &FlagTable_Lohner[lv][3],
                                                                                     &FlagTable_Lohner[lv][4] );
//       OPT__FLAG_USER has OPT__FLAG_USER_NUM columns to be loaded
         else if ( FlagMode == 5 )
         {
            char *input_line0 = input_line;
            int offset;

            sscanf( input_line, "%d%n", &Trash, &offset );
            input_line += offset;

            for (int t=0; t<OPT__FLAG_USER_NUM; t++)
            {
               sscanf( input_line, "%lf%n", &FlagTable_User[lv][t], &offset );
               input_line += offset;
            }

            input_line = input_line0;  // reset input_line
         }
//       OPT__FLAG_NPAR_PATCH/CELL load integers
         else if ( FlagMode == 6 )  sscanf( input_line, "%d%d",  &Trash, &FlagTable_NParPatch[lv] );
         else if ( FlagMode == 7 )  sscanf( input_line, "%d%d",  &Trash, &FlagTable_NParCell [lv] );
//       OPT__FLAG_INTERFERENCE loads four columns
         else if ( FlagMode == 14 ) sscanf( input_line, "%d%lf%lf%lf%lf", &Trash, &FlagTable_Interference[lv][0],
                                                                                  &FlagTable_Interference[lv][1],
                                                                                  &FlagTable_Interference[lv][2],
                                                                                  &FlagTable_Interference[lv][3]
                                                                                  );
//       OPT__FLAG_SPECTRAL loads two columns
         else if ( FlagMode == 15 ) sscanf( input_line, "%d%lf%lf", &Trash,       &FlagTable_Spectral[lv][0],
                                                                                  &FlagTable_Spectral[lv][1]
                                                                                  );
//       OPT__FLAG_ANGULAR loads three columns
         else if ( FlagMode == 16 )
         {
            sscanf( input_line, "%d%lf%lf%lf", &Trash, &FlagTable_Angular[lv][0], &FlagTable_Angular[lv][1], &FlagTable_Angular[lv][2] );
            if ( FlagTable_Angular[lv][0] >= 0.0  &&
                 FlagTable_Angular[lv][1] >= 0.0  &&
                 FlagTable_Angular[lv][0] >  0.5*FlagTable_Angular[lv][1] )
               Aux_Message( stderr, "WARNING : AngRes_Max (%20.14e) > 0.5*AngRes_Min (%20.14e) at lv %d (note that AngRes_Max has higher priority) !!\n",
                            FlagTable_Angular[lv][0], 0.5*FlagTable_Angular[lv][1], lv );
         }

//       others use the default format: (integer, double)
         else                       sscanf( input_line, "%d%lf", &Trash, &FlagTable[FlagMode][lv] );
      } // for (int lv=0; lv<MAX_LEVEL; lv++)

      fclose( File );

   } // for (int FlagMode=0; FlagMode<NFlagMode; FlagMode++)


   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Load_FlagCriteria
