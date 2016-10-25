#include "Copyright.h"
#include "GAMER.h"

extern double Myr;
extern double Msun;

extern double BH_SinkMass;
extern double BH_SinkMomX;
extern double BH_SinkMomY;
extern double BH_SinkMomZ;
extern double BH_SinkMomXAbs;
extern double BH_SinkMomYAbs;
extern double BH_SinkMomZAbs;
extern double BH_SinkEk;
extern double BH_SinkEt;
extern int    BH_SinkNCell;




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_RecordUser
// Description :  Record user-specified information
//
// Note        :  1. Please turn on the option "OPT__RECORD_USER" 
//                2. This function will be called both during the program initialization and after each full update
// 
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Aux_RecordUser( )
{

   const char FileName[] = "Record__BHAccretionRate";
   static bool FirstTime = true;
   static double Time0, dTime;

   if ( FirstTime )
   {
//    get the total number of cells within the void region
      int SinkNCell_Sum;

      MPI_Reduce( &BH_SinkNCell, &SinkNCell_Sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );


//    header
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )    Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#Number of void cells : %d\n", SinkNCell_Sum );
         fprintf( File_User, "#%9s%16s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n",
                  "Step", "Time [yr]", "Mass [Msun]", "Time [yr]", "dM/dt [Msun/yr]",
                  "MomX [g*cm/s]", "MomY [g*cm/s]", "MomZ [g*cm/s]", "MomXAbs [g*cm/s]", "MomYAbs [g*cm/s]", "MomZAbs [g*cm/s]",
                  "Ek [erg]", "Et [erg]" );
         fclose( File_User );
      }

      FirstTime = false;
      Time0     = Time[0];
   }

   else
   {
      double Mass_Sum, MomX_Sum, MomY_Sum, MomZ_Sum, MomXAbs_Sum, MomYAbs_Sum, MomZAbs_Sum, Ek_Sum, Et_Sum;

      MPI_Reduce( &BH_SinkMass,    &Mass_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkMomX,    &MomX_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkMomY,    &MomY_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkMomZ,    &MomZ_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkMomXAbs, &MomXAbs_Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkMomYAbs, &MomYAbs_Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkMomZAbs, &MomZAbs_Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkEk,      &Ek_Sum,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &BH_SinkEt,      &Et_Sum,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      Mass_Sum    *= UNIT_M/Msun;
      MomX_Sum    *= UNIT_M*UNIT_L/UNIT_T;
      MomY_Sum    *= UNIT_M*UNIT_L/UNIT_T;
      MomZ_Sum    *= UNIT_M*UNIT_L/UNIT_T;
      MomXAbs_Sum *= UNIT_M*UNIT_L/UNIT_T;
      MomYAbs_Sum *= UNIT_M*UNIT_L/UNIT_T;
      MomZAbs_Sum *= UNIT_M*UNIT_L/UNIT_T;
      Ek_Sum      *= UNIT_E;
      Et_Sum      *= UNIT_E;

      if ( MPI_Rank == 0 )
      {
         const double yr = Myr*1.0e-6;

         dTime  = Time[0] - Time0;
         Time0  = Time[0];
         dTime *= UNIT_T/yr;

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "%10ld%16.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e\n",
                  Step, Time[0]*UNIT_T/yr, Mass_Sum, dTime, Mass_Sum/dTime,
                  MomX_Sum, MomY_Sum, MomZ_Sum, MomXAbs_Sum, MomYAbs_Sum, MomZAbs_Sum, Ek_Sum, Et_Sum );
         fclose( File_User );
      }
   } // if ( FirstTime ) ... else ...

// reset the cumulative mass to zero
   BH_SinkMass    = 0.0;
   BH_SinkMomX    = 0.0;
   BH_SinkMomY    = 0.0;
   BH_SinkMomZ    = 0.0;
   BH_SinkMomXAbs = 0.0;
   BH_SinkMomYAbs = 0.0;
   BH_SinkMomZAbs = 0.0;
   BH_SinkEk      = 0.0;
   BH_SinkEt      = 0.0;

} // FUNCTION : Aux_RecordUser


