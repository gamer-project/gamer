#include "GAMER.h"


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_SubGrid_Fields
// Description :  Populate SubGridField[], SubGridDerLabel[], and the associated counters from
//                Input__Sub_Grid
//
// Note        :  1. Input__Sub_Grid is REQUIRED when OPT__OUTPUT_SUBDIV_GRID is enabled; a missing
//                   file is a fatal error (no silent fallback)
//                2. Each non-comment line must contain exactly one field label (case-sensitive), which
//                   can be either a native field (FieldLabel[v]) or a derived field recognized by
//                   Output_DumpData_Total_HDF5() (e.g., Pote, ParDens/TotalDens, Pres, Temp, Entr,
//                   CCMagX/Y/Z, user-defined derived fields, ...)
//                3. Sub-dump derived fields are controlled solely by this list (OPT__OUTPUT_* flags
//                   apply to main dumps only); listing ParDens/TotalDens without OPT__OUTPUT_PAR_DENS
//                   (which selects the deposition mode) or a user-defined field without
//                   OPT__OUTPUT_USER_FIELD (which initializes the labels) is a fatal error
//                4. Called by Init_GAMER() when OPT__OUTPUT_SUBDIV_GRID is true
//                   --> after Init_DerivedField_User_Ptr() so that UserDerField_Label[] is available
//-------------------------------------------------------------------------------------------------------
void Init_SubGrid_Fields()
{

// initialize
   for (int v=0; v<NCOMP_TOTAL; v++)   SubGridField[v] = false;
   SubGridField_Num    = 0;
   SubGridDerLabel_Num = 0;

   char FileName[] = "Input__Sub_Grid";
   FILE *File      = fopen( FileName, "r" );

// the field-list file is required when OPT__OUTPUT_SUBDIV_GRID is enabled
   if ( File == NULL )
      Aux_Error( ERROR_INFO, "%s is required when OPT__OUTPUT_SUBDIV_GRID is enabled "
                 "(see the template at example/input/Input__Sub_Grid) !!\n", FileName );

// collect the derived-field labels recognized by Output_DumpData_Total_HDF5()
   const char *DerList[NFIELD_STORED_MAX];
   int NDer = 0;
#  ifdef GRAVITY
   DerList[NDer++] = PotLabel;
#  endif
#  ifdef MASSIVE_PARTICLES
   DerList[NDer++] = "ParDens";
   DerList[NDer++] = "TotalDens";
#  endif
#  ifdef MHD
   DerList[NDer++] = "CCMagX";
   DerList[NDer++] = "CCMagY";
   DerList[NDer++] = "CCMagZ";
#  endif
#  if ( MODEL == HYDRO )
   DerList[NDer++] = "Pres";
   DerList[NDer++] = "Temp";
   DerList[NDer++] = "Entr";
   DerList[NDer++] = "Cs";
   DerList[NDer++] = "DivVel";
   DerList[NDer++] = "Mach";
#  ifdef MHD
   DerList[NDer++] = "DivMag";
#  endif
#  ifdef SRHD
   DerList[NDer++] = "Lrtz";
   DerList[NDer++] = "VelX";
   DerList[NDer++] = "VelY";
   DerList[NDer++] = "VelZ";
   DerList[NDer++] = "Enth";
#  endif
#  ifdef SUPPORT_GRACKLE
   DerList[NDer++] = "GrackleTemp";
   DerList[NDer++] = "GrackleMu";
   DerList[NDer++] = "GrackleTCool";
#  endif
#  endif // #if ( MODEL == HYDRO )

   char *Line      = new char [MAX_STRING];
   char  FirstItem[MAX_STRING];

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      if ( sscanf(Line, "%s", FirstItem) <= 0  ||  FirstItem[0] == '#' )   continue;

//    (a) native fields
      bool found = false;
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         if ( strcmp(FieldLabel[v], FirstItem) == 0 )
         {
            SubGridField[v] = true;
            SubGridField_Num++;
            found = true;
            break;
         }
      }

//    (b) derived fields (built-in and user-defined)
      if ( !found )
      {
         for (int t=0; t<NDer; t++)
         {
            if ( strcmp(DerList[t], FirstItem) == 0 )
            {
               found = true;
               break;
            }
         }

         if ( !found  &&  OPT__OUTPUT_USER_FIELD )
         for (int v=0; v<UserDerField_Num; v++)
         {
            if ( strcmp(UserDerField_Label[v], FirstItem) == 0 )
            {
               found = true;
               break;
            }
         }

         if ( found )
         {
#           ifdef MASSIVE_PARTICLES
            if (  ( strcmp(FirstItem, "ParDens") == 0  ||  strcmp(FirstItem, "TotalDens") == 0 )  &&
                  OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_NONE  )
               Aux_Error( ERROR_INFO, "\"%s\" in %s requires OPT__OUTPUT_PAR_DENS, which selects the deposition mode !!\n",
                          FirstItem, FileName );
#           endif

            if ( SubGridDerLabel_Num >= NFIELD_STORED_MAX )
               Aux_Error( ERROR_INFO, "number of derived fields in %s exceeds NFIELD_STORED_MAX (%d) !!\n",
                          FileName, NFIELD_STORED_MAX );

            sprintf( SubGridDerLabel[ SubGridDerLabel_Num ++ ], "%s", FirstItem );
         }
      }

      if ( !found )
         Aux_Error( ERROR_INFO, "unknown field label \"%s\" in %s "
                    "(user-defined derived fields additionally require OPT__OUTPUT_USER_FIELD) !!\n",
                    FirstItem, FileName );
   }

   fclose( File );
   delete [] Line;

   if ( SubGridField_Num == 0  &&  SubGridDerLabel_Num == 0 )
      Aux_Error( ERROR_INFO, "%s contains no valid field labels !!\n", FileName );

} // FUNCTION : Init_SubGrid_Fields



//-------------------------------------------------------------------------------------------------------
// Function    :  SubGrid_DerFieldSelected
// Description :  Check whether a derived field is selected for grid sub-dumps
//
// Note        :  1. Invoked by Output_DumpData_Total_HDF5() when writing grid sub-dumps (SubGridMode)
//
// Parameter   :  Label : Derived-field label to check (e.g., "Pres", PotLabel, ...)
//
// Return      :  true if the field should be written in grid sub-dumps
//-------------------------------------------------------------------------------------------------------
bool SubGrid_DerFieldSelected( const char *Label )
{

   for (int t=0; t<SubGridDerLabel_Num; t++)
      if ( strcmp( SubGridDerLabel[t], Label ) == 0 )    return true;

   return false;

} // FUNCTION : SubGrid_DerFieldSelected
