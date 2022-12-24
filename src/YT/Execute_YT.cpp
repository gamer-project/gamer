#include "GAMER.h"

static void Write_ExecuteYTRecord();




//-------------------------------------------------------------------------------------------------------
// Function    :  Execute_YT
// Description :  Trigger the yt inline analysis functions YT_Inline()
//
// Note        :  
//
// Parameter   :  Stage : 0 : beginning of the run
//                        1 : during the evolution
//                        2 : end of the run
//-------------------------------------------------------------------------------------------------------
void Execute_YT( const int Stage )
{

// check
   if ( Stage < 0  ||  Stage > 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Stage", Stage );


// set the first execute time
   static int ExecuteYTTableID;

   if ( Stage == 0 )
   {
      switch ( OPT__EXECUTE_YT_MODE )
      {
         case EXECUTE_YT_CONST_STEP :
         break;   // do nothing

         case EXECUTE_YT_CONST_DT :
         {
            if ( OPT__INIT != INIT_BY_RESTART  ||  OPT__RESTART_RESET  ||  OPT__EXECUTE_YT_RESTART )
               ExecuteYTTime = Time[0];

            else
            {
               ExecuteYTTime = ( int(Time[0]/EXECUTE_YT_DT) + 1.0 )*EXECUTE_YT_DT;

//             be careful about round-off errors
               if (   (  ExecuteYTTime <= Time[0]  )                                            ||
                      (  Time[0] != 0.0 && fabs( (Time[0]-ExecuteYTTime)/Time[0] ) < 1.0e-8  )  ||
                      (  Time[0] == 0.0 && fabs(  Time[0]-ExecuteYTTime          ) < 1.0e-12 )      )   ExecuteYTTime += EXECUTE_YT_DT;
            }
         }
         break;

         case EXECUTE_YT_USE_TABLE :
         {
            if ( OPT__INIT != INIT_BY_RESTART  ||  OPT__RESTART_RESET  ||  OPT__EXECUTE_YT_RESTART )
            {
               for (ExecuteYTTableID=0; ExecuteYTTableID<ExecuteYTTable_NExecute; ExecuteYTTableID++)
               {
                  ExecuteYTTime = ExecuteYTTable[ExecuteYTTableID];

                  if (   (  ExecuteYTTime >= Time[0]  )                                            ||
                         (  Time[0] != 0.0 && fabs( (Time[0]-ExecuteYTTime)/Time[0] ) < 1.0e-8  )  ||
                         (  Time[0] == 0.0 && fabs(  Time[0]-ExecuteYTTime          ) < 1.0e-12 )      )   break;
               }
            }

            else
            {
               for (ExecuteYTTableID=0; ExecuteYTTableID<ExecuteYTTable_NExecute; ExecuteYTTableID++)
               {
                  ExecuteYTTime = ExecuteYTTable[ExecuteYTTableID];

                  if (   (  ExecuteYTTime >= Time[0]  )                                            &&
                        !(  Time[0] != 0.0 && fabs( (Time[0]-ExecuteYTTime)/Time[0] ) < 1.0e-8  )  &&
                        !(  Time[0] == 0.0 && fabs(  Time[0]-ExecuteYTTime          ) < 1.0e-12 )      )   break;
               }
            }

            if ( ExecuteYTTableID >= ExecuteYTTable_NExecute )
               Aux_Error( ERROR_INFO, "no proper data dump time is found in the execute YT table !!\n" );
         }
         break;

      } // switch ( OPT__EXECUTE_YT_MODE )
   } // if ( Stage == 0 )


// do not execute yt inline analysis for the initial data for the restart run (unless enabling OPT__EXECUTE_YT_RESTART)
   if ( OPT__INIT == INIT_BY_RESTART  &&  Stage == 0  &&  !OPT__RESTART_RESET  &&  !OPT__EXECUTE_YT_RESTART )   return;


// determine whether or not to execute yt inline analysis during the simulation
   static int PreviousExecuteYTStep   = -999;
   bool ExecuteYT                     = false;

   switch ( OPT__EXECUTE_YT_MODE )
   {
      case EXECUTE_YT_CONST_STEP :   if ( Step%EXECUTE_YT_STEP == 0 )     ExecuteYT = true;
                                     break;

      case EXECUTE_YT_CONST_DT   :   // rely on the same execution statement as EXECUTE_YT_USE_TABLE; order is important (first case EXECUTE_YT_CONST_DT then case EXECUTE_YT_USE_TABLE)!!
      case EXECUTE_YT_USE_TABLE  :   if (   ( Time[0] != 0.0 && fabs( (Time[0]-ExecuteYTTime)/Time[0] ) < 1.0e-8  )
                                         || ( Time[0] == 0.0 && fabs(  Time[0]-ExecuteYTTime          ) < 1.0e-12 )   )
                                        ExecuteYT = true;
                                     break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__EXECUTE_YT_MODE", OPT__EXECUTE_YT_MODE );
   } // switch ( OPT__EXECUTE_YT_MODE )


// always execute yt inline analysis in the end of the simulation (unless it has already been analyzed)
   if ( Stage == 2 )
   {
      if ( Step == PreviousExecuteYTStep )     ExecuteYT = false;
      else                                     ExecuteYT = true;
   }


// Should be done at Output_DumpData.cpp
//// set the acceleration of tracer particles to zero to make the output deterministic
//#  if ( defined TRACER  &&  defined STORE_PAR_ACC )
//   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
//   {
//      if ( amr->Par->Type[p] == PTYPE_TRACER )
//      {
//         amr->Par->AccX[p] = (real)0.0;
//         amr->Par->AccY[p] = (real)0.0;
//         amr->Par->AccZ[p] = (real)0.0;
//      }
//   }
//#  endif


// execute yt inline analysis
   if ( ExecuteYT )
   {
      YT_Inline();

      Write_ExecuteYTRecord();

      ExecuteYTID ++;

      if ( OPT__EXECUTE_YT_MODE == EXECUTE_YT_CONST_DT  )  ExecuteYTTime = Time[0] + EXECUTE_YT_DT;
      if ( OPT__EXECUTE_YT_MODE == EXECUTE_YT_USE_TABLE )  ExecuteYTTime = ExecuteYTTable[ ++ExecuteYTTableID ];

      PreviousExecuteYTStep = Step;
   } // if ( ExecuteYT )

} // FUNCTION : Execute_YT



//-------------------------------------------------------------------------------------------------------
// Function    :  Write_ExecuteYTRecord
// Description :  Record the information of each yt  in the file "Record__ExecuteYT"
//-------------------------------------------------------------------------------------------------------
void Write_ExecuteYTRecord()
{

   const char FileName[] = "Record__ExecuteYT";


// create the "Record__ExecuteYT" file at the first execution
   static bool FirstTime = true;

   if ( MPI_Rank == 0  &&  FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      else
      {
         FILE *File = fopen( FileName, "w" );
         fprintf( File, "%6s\t\t%12s\t\t%9s\n", "ExecuteYTID", "Time", "Step" );
         fclose( File );
      }

      FirstTime = false;
   }


// record the information of data dump in the file "Record__Dump"
   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );
      fprintf( File, "%6d\t\t%20.14e\t\t%9ld\n", ExecuteYTID, Time[0], Step );
      fclose( File );
   }

} // FUNCTION : Write_ExecuteYTRecord
