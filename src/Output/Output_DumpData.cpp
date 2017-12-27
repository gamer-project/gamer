#include "GAMER.h"

static void Write_DumpRecord();
extern void (*Output_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData
// Description :  Trigger the output functions "Output_DumpData_Total, Output_DumpData_Part, Output_User,
//                Output_BasePowerSpectrum, Par_Output_TextFile"
//
// Note        :  1. The function pointer "Output_User_Ptr" points to "Output_User()" by default
//                   but may be overwritten by various test problem initializers
//
// Parameter   :  Stage : 0 : beginning of the run
//                        1 : during the evolution
//                        2 : end of the run
//-------------------------------------------------------------------------------------------------------
void Output_DumpData( const int Stage )
{

// check
   if ( Stage < 0  ||  Stage > 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Stage", Stage );


// nothing to do if all output options are off
#  ifdef PARTICLE
   if ( !OPT__OUTPUT_TOTAL && !OPT__OUTPUT_PART && !OPT__OUTPUT_USER && !OPT__OUTPUT_BASEPS && !OPT__OUTPUT_PAR_TEXT )
#  else
   if ( !OPT__OUTPUT_TOTAL && !OPT__OUTPUT_PART && !OPT__OUTPUT_USER && !OPT__OUTPUT_BASEPS )
#  endif
      return;


// set the first dump time
   static int DumpTableID;

   if ( Stage == 0 )
   {
      switch ( OPT__OUTPUT_MODE )
      {
         case OUTPUT_CONST_STEP :
         break;   // do nothing

         case OUTPUT_CONST_DT :
         {
            if ( OPT__INIT != INIT_BY_RESTART  ||  OPT__RESTART_RESET )
               DumpTime = Time[0];

            else
            {
               DumpTime = ( int(Time[0]/OUTPUT_DT) + 1.0 )*OUTPUT_DT;

//             be careful about round-off errors
               if (   (  DumpTime <= Time[0]  )                                            ||
                      (  Time[0] != 0.0 && fabs( (Time[0]-DumpTime)/Time[0] ) < 1.0e-8  )  ||
                      (  Time[0] == 0.0 && fabs(  Time[0]-DumpTime          ) < 1.0e-12 )      )   DumpTime += OUTPUT_DT;
            }
         }
         break;

         case OUTPUT_USE_TABLE :
         {
            if ( OPT__INIT != INIT_BY_RESTART  ||  OPT__RESTART_RESET )
            {
               for (DumpTableID=0; DumpTableID<DumpTable_NDump; DumpTableID++)
               {
                  DumpTime = DumpTable[DumpTableID];

                  if (   (  DumpTime >= Time[0]  )                                            ||
                         (  Time[0] != 0.0 && fabs( (Time[0]-DumpTime)/Time[0] ) < 1.0e-8  )  ||
                         (  Time[0] == 0.0 && fabs(  Time[0]-DumpTime          ) < 1.0e-12 )      )   break;
               }
            }

            else
            {
               for (DumpTableID=0; DumpTableID<DumpTable_NDump; DumpTableID++)
               {
                  DumpTime = DumpTable[DumpTableID];

                  if (   (  DumpTime >= Time[0]  )                                            &&
                        !(  Time[0] != 0.0 && fabs( (Time[0]-DumpTime)/Time[0] ) < 1.0e-8  )  &&
                        !(  Time[0] == 0.0 && fabs(  Time[0]-DumpTime          ) < 1.0e-12 )      )   break;
               }
            }

            if ( DumpTableID >= DumpTable_NDump )
               Aux_Error( ERROR_INFO, "no proper data dump time is found in the dump table !!\n" );
         }
         break;

      } // switch ( OPT__OUTPUT_MODE )
   } // if ( Stage == 0 )


// do not output the initial data for the restart run
   if ( OPT__INIT == INIT_BY_RESTART  &&  Stage == 0  &&  !OPT__RESTART_RESET )  return;


// set the file names for all output functions
   char FileName_Total[50], FileName_Part[50], FileName_Temp[50], FileName_PS[50];
#  ifdef PARTICLE
   char FileName_Particle[50];
#  endif

   if ( OPT__OUTPUT_TOTAL )   sprintf( FileName_Total, "Data_%06d", DumpID );

   if ( OPT__OUTPUT_PART )
   {
      switch ( OPT__OUTPUT_PART )
      {
         case OUTPUT_XY :    sprintf( FileName_Temp, "XYslice_z%.3f_%06d", OUTPUT_PART_Z, DumpID );   break;
         case OUTPUT_YZ :    sprintf( FileName_Temp, "YZslice_x%.3f_%06d", OUTPUT_PART_X, DumpID );   break;
         case OUTPUT_XZ :    sprintf( FileName_Temp, "XZslice_y%.3f_%06d", OUTPUT_PART_Y, DumpID );   break;
         case OUTPUT_X  :    sprintf( FileName_Temp, "Xline_y%.3f_z%.3f_%06d", OUTPUT_PART_Y, OUTPUT_PART_Z, DumpID );  break;
         case OUTPUT_Y  :    sprintf( FileName_Temp, "Yline_x%.3f_z%.3f_%06d", OUTPUT_PART_X, OUTPUT_PART_Z, DumpID );  break;
         case OUTPUT_Z  :    sprintf( FileName_Temp, "Zline_x%.3f_y%.3f_%06d", OUTPUT_PART_X, OUTPUT_PART_Y, DumpID );  break;
         case OUTPUT_DIAG :  sprintf( FileName_Temp, "Diag_%06d", DumpID );   break;
         default :           Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__OUTPUT_PART", OPT__OUTPUT_PART );
      } // switch ( OPT__OUTPUT_PART )

      if ( OPT__OUTPUT_BASE )
      {
         sprintf( FileName_Part, "%s", "Base" );
         strcat( FileName_Part, FileName_Temp );
      }
      else
         strcpy( FileName_Part, FileName_Temp );

   } // if ( OPT__OUTPUT_PART )

   if ( OPT__OUTPUT_BASEPS )
      sprintf( FileName_PS, "PowerSpec_%06d", DumpID );

#  ifdef PARTICLE
   if ( OPT__OUTPUT_PAR_TEXT )
      sprintf( FileName_Particle, "Particle_%06d", DumpID );
#  endif


// determine whether or not to output data during the simulation
   static int PreviousDumpStep   = -999;
   bool OutputData               = false;

   switch ( OPT__OUTPUT_MODE )
   {
      case OUTPUT_CONST_STEP :   if ( Step%OUTPUT_STEP == 0 )     OutputData = true;
                                 break;

      case OUTPUT_CONST_DT :
      case OUTPUT_USE_TABLE :    if (   ( Time[0] != 0.0 && fabs( (Time[0]-DumpTime)/Time[0] ) < 1.0e-8  )
                                     || ( Time[0] == 0.0 && fabs(  Time[0]-DumpTime          ) < 1.0e-12 )   )
                                    OutputData = true;
                                 break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__OUTPUT_MODE", OPT__OUTPUT_MODE );
   } // switch ( OPT__OUTPUT_MODE )


// always output data in the end of the simulation (unless it has already been output)
   if ( Stage == 2 )
   {
      if ( Step == PreviousDumpStep )     OutputData = false;
      else                                OutputData = true;
   }


// dump data if any process has detected the file named "DUMP_GAMER_DUMP"
   int OutputData_RunTime = false;

// enable this functionality only if OPT__MANUAL_CONTROL is on
   if ( OPT__MANUAL_CONTROL )    Output_DumpManually( OutputData_RunTime );


// output data
   if ( OutputData || OutputData_RunTime )
   {
//    apply various corrections (e.g., synchronize particles, restrict data, recalculate potential and particle acceleration)
//    before dumpting data --> for bitwise reproducibility
      if ( OPT__CORR_AFTER_ALL_SYNC == CORR_AFTER_SYNC_BEFORE_DUMP  &&  Stage != 0 )  Flu_CorrAfterAllSync();

      if ( OPT__OUTPUT_TOTAL )            Output_DumpData_Total( FileName_Total );
      if ( OPT__OUTPUT_PART  )            Output_DumpData_Part( OPT__OUTPUT_PART, OPT__OUTPUT_BASE, OUTPUT_PART_X,
                                                                OUTPUT_PART_Y, OUTPUT_PART_Z, FileName_Part );
      if ( OPT__OUTPUT_USER  &&
           Output_User_Ptr != NULL )      Output_User_Ptr();
#     ifdef GRAVITY
      if ( OPT__OUTPUT_BASEPS )           Output_BasePowerSpectrum( FileName_PS );
#     endif
#     ifdef PARTICLE
      if ( OPT__OUTPUT_PAR_TEXT )         Par_Output_TextFile( FileName_Particle );
#     endif

      Write_DumpRecord();

      DumpID ++;

      if ( OutputData )
      {
         if ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT  )  DumpTime = Time[0] + OUTPUT_DT;
         if ( OPT__OUTPUT_MODE == OUTPUT_USE_TABLE )  DumpTime = DumpTable[ ++DumpTableID ];
      }

      PreviousDumpStep = Step;
   } // if ( OutputData || OutputData_RunTime )

} // FUNCTION : Output_DumpData



//-------------------------------------------------------------------------------------------------------
// Function    :  Write_DumpRecord
// Description :  Record the information of each data dump in the file "Record__Dump"
//-------------------------------------------------------------------------------------------------------
void Write_DumpRecord()
{

   const char FileName[] = "Record__Dump";


// create the "Record__Dump" file at the first dump
   static bool FirstTime = true;

   if ( MPI_Rank == 0  &&  FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      else
      {
         FILE *File = fopen( FileName, "w" );
         fprintf( File, "%6s\t\t%20s\t\t%9s\n", "DumpID", "Time", "Step" );
         fclose( File );
      }

      FirstTime = false;
   }


// record the information of data dump in the file "Record__Dump"
   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );
      fprintf( File, "%6d\t\t%20.14e\t\t%9ld\n", DumpID, Time[0], Step );
      fclose( File );
   }

} // FUNCTION : Write_DumpRecord
