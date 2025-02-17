#include "GAMER.h"

extern Timer_t Timer_OutputWalltime;

static void Write_DumpRecord();




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData
// Description :  Trigger the output functions Output_DumpData_Total(), Output_DumpData_Part(), Output_User(),
//                Output_BasePowerSpectrum(), Par_Output_TextFile(), Par_Output_BinaryFile()
//
// Note        :  1. For OUTPUT_USER, the function pointer "Output_User_Ptr" must be set by a
//                   test problem initializer
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
   if ( !OPT__OUTPUT_TOTAL && !OPT__OUTPUT_PART && !OPT__OUTPUT_USER && !OPT__OUTPUT_BASEPS && !OPT__OUTPUT_PAR_MODE )
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
            if ( OPT__INIT != INIT_BY_RESTART  ||  OPT__RESTART_RESET  ||  OPT__OUTPUT_RESTART )
               DumpTime = Time[0];

            else
            {
               DumpTime = round( int(Time[0]/OUTPUT_DT) + 1.0 )*OUTPUT_DT;

//             be careful about round-off errors
               if (   (  DumpTime <= Time[0]  )                                            ||
                      (  Time[0] != 0.0 && fabs( (Time[0]-DumpTime)/Time[0] ) < 1.0e-8  )  ||
                      (  Time[0] == 0.0 && fabs(  Time[0]-DumpTime          ) < 1.0e-12 )      )   DumpTime += OUTPUT_DT;
            }
         }
         break;

         case OUTPUT_USE_TABLE :
         {
            if ( OPT__INIT != INIT_BY_RESTART  ||  OPT__RESTART_RESET  ||  OPT__OUTPUT_RESTART )
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


// do not output the initial data for the restart run (unless enabling OPT__OUTPUT_RESTART)
   if ( OPT__INIT == INIT_BY_RESTART  &&  Stage == 0  &&  !OPT__RESTART_RESET  &&  !OPT__OUTPUT_RESTART )   return;


// set the file names for all output functions
   char FileName_Total[2*MAX_STRING], FileName_Part[2*MAX_STRING], FileName_Temp[2*MAX_STRING], FileName_PS[2*MAX_STRING];
#  ifdef PARTICLE
   char FileName_Particle[2*MAX_STRING];
#  endif

   if ( OPT__OUTPUT_TOTAL )
   {
      sprintf( FileName_Total, "%s/Data_%06d", OUTPUT_DIR, DumpID );
   }

   if ( OPT__OUTPUT_PART )
   {
      sprintf( FileName_Part, "%s/", OUTPUT_DIR );
      switch ( OPT__OUTPUT_PART )
      {
         case OUTPUT_XY   :  sprintf( FileName_Temp, "XYslice_z%.3f_%06d", OUTPUT_PART_Z, DumpID );   break;
         case OUTPUT_YZ   :  sprintf( FileName_Temp, "YZslice_x%.3f_%06d", OUTPUT_PART_X, DumpID );   break;
         case OUTPUT_XZ   :  sprintf( FileName_Temp, "XZslice_y%.3f_%06d", OUTPUT_PART_Y, DumpID );   break;
         case OUTPUT_X    :  sprintf( FileName_Temp, "Xline_y%.3f_z%.3f_%06d", OUTPUT_PART_Y, OUTPUT_PART_Z, DumpID );  break;
         case OUTPUT_Y    :  sprintf( FileName_Temp, "Yline_x%.3f_z%.3f_%06d", OUTPUT_PART_X, OUTPUT_PART_Z, DumpID );  break;
         case OUTPUT_Z    :  sprintf( FileName_Temp, "Zline_x%.3f_y%.3f_%06d", OUTPUT_PART_X, OUTPUT_PART_Y, DumpID );  break;
         case OUTPUT_DIAG :  sprintf( FileName_Temp, "Diag_%06d", DumpID );   break;
         case OUTPUT_BOX  :  sprintf( FileName_Temp, "Box_%06d", DumpID );   break;
         default :           Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__OUTPUT_PART", OPT__OUTPUT_PART );
      } // switch ( OPT__OUTPUT_PART )

      if ( OPT__OUTPUT_BASE )
      {
         strcat( FileName_Part, "Base" );
         strcat( FileName_Part, FileName_Temp );
      }
      else
         strcat( FileName_Part, FileName_Temp );

   } // if ( OPT__OUTPUT_PART )

   if ( OPT__OUTPUT_BASEPS )
      sprintf( FileName_PS, "%s/PowerSpec_%06d", OUTPUT_DIR, DumpID );

#  ifdef PARTICLE
   if ( OPT__OUTPUT_PAR_MODE == OUTPUT_PAR_TEXT )
      sprintf( FileName_Particle, "%s/Particle_%06d.txt", OUTPUT_DIR, DumpID );
   if ( OPT__OUTPUT_PAR_MODE == OUTPUT_PAR_CBIN )
      sprintf( FileName_Particle, "%s/Particle_%06d.cbin", OUTPUT_DIR, DumpID );
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


// dump data if the elapsed walltime exceeds the user-defined walltime
   int OutputData_Walltime = false;

   if ( OUTPUT_WALLTIME > 0.0 )
   {
      Timer_OutputWalltime.Stop();

      double ElapsedWalltime = Timer_OutputWalltime.GetValue();

#     ifndef SERIAL
      MPI_Allreduce( MPI_IN_PLACE, &ElapsedWalltime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#     endif

      switch ( OUTPUT_WALLTIME_UNIT )
      {
         case 0 :                               break;
         case 1 : ElapsedWalltime /=    60.0;   break;
         case 2 : ElapsedWalltime /=  3600.0;   break;
         case 3 : ElapsedWalltime /= 86400.0;   break;
         default: Aux_Error( ERROR_INFO, "unsupported unit (%d) for output walltime !!\n", OUTPUT_WALLTIME_UNIT );
      }


      if ( ElapsedWalltime >= OUTPUT_WALLTIME )
      {
         OutputData_Walltime = true;
         Timer_OutputWalltime.Reset();
      }

      Timer_OutputWalltime.Start();
   } // if ( OUTPUT_WALLTIME > 0.0 )


// set potential to zero when disabling both self-gravity and external potential
// to make outputs deterministic and more reasonable
#  ifdef GRAVITY
   if ( !OPT__SELF_GRAVITY && !OPT__EXT_POT )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
         amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i] = (real)0.0;
   }
#  endif


// set the acceleration of tracer particles to zero to make outputs deterministic
#  if ( defined TRACER  &&  defined STORE_PAR_ACC )
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
      if ( amr->Par->Type[p] == PTYPE_TRACER )
      {
         amr->Par->AccX[p] = (real)0.0;
         amr->Par->AccY[p] = (real)0.0;
         amr->Par->AccZ[p] = (real)0.0;
      }
   }
#  endif


// output data
   if ( OutputData || OutputData_RunTime || OutputData_Walltime )
   {
//    sort real patches by their load-balance indices to improve bitwise reproducibility
#     ifdef LOAD_BALANCE
      if ( OPT__SORT_PATCH_BY_LBIDX )
      {
#        ifdef PARTICLE
         const double ParWeight         = amr->LB->Par_Weight;
#        else
         const double ParWeight         = 0.0;
#        endif
         const bool   Redistribute_Yes  = true;
         const bool   SendGridData_Yes  = true;
         const bool   ResetLB_Yes       = true;
         const bool   SortRealPatch_Yes = true;
         const int    AllLv             = -1;

         LB_Init_LoadBalance( Redistribute_Yes, SendGridData_Yes, ParWeight, ResetLB_Yes, SortRealPatch_Yes, AllLv );
      }
#     endif // #ifdef LOAD_BALANCE

//    apply various corrections (e.g., synchronize particles, restrict data, recalculate potential and particle acceleration)
//    before dumpting data --> for bitwise reproducibility
      if ( OPT__CORR_AFTER_ALL_SYNC == CORR_AFTER_SYNC_BEFORE_DUMP  &&  Stage != 0 )  Flu_CorrAfterAllSync();

//    perform user-specified work before dumping data
      if ( Output_UserWorkBeforeOutput_Ptr != NULL )  Output_UserWorkBeforeOutput_Ptr();

//    start dumping data
      if ( OPT__OUTPUT_TOTAL )            Output_DumpData_Total( FileName_Total );
      if ( OPT__OUTPUT_PART  )            Output_DumpData_Part( OPT__OUTPUT_PART, OPT__OUTPUT_BASE, OUTPUT_PART_X,
                                                                OUTPUT_PART_Y, OUTPUT_PART_Z, FileName_Part );
      if ( OPT__OUTPUT_USER )
      {
         if ( Output_User_Ptr != NULL )   Output_User_Ptr();
         else
            Aux_Error( ERROR_INFO, "Output_User_Ptr == NULL for OPT__OUTPUT_USER !!\n" );
      }
#     ifdef SUPPORT_FFTW
      if ( OPT__OUTPUT_BASEPS )           Output_BasePowerSpectrum( FileName_PS, _TOTAL_DENS );
#     endif
#     ifdef PARTICLE
      if ( OPT__OUTPUT_PAR_MODE == OUTPUT_PAR_TEXT )  Par_Output_TextFile( FileName_Particle );
      if ( OPT__OUTPUT_PAR_MODE == OUTPUT_PAR_CBIN )  Par_Output_BinaryFile( FileName_Particle );
#     endif

      Write_DumpRecord();

      DumpID ++;

      if ( OutputData )
      {
         if ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT  )  DumpTime = round( Time[0]/OUTPUT_DT + 1.0 )*OUTPUT_DT;
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

   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Dump", OUTPUT_DIR );


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
