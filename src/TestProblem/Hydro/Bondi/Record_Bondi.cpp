#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )



// specific global variables declared in Init_TestProb_Hydro_Bondi.cpp
// =======================================================================================
extern double Bondi_MassBH;
extern double Bondi_SinkMass;
extern double Bondi_SinkMomX;
extern double Bondi_SinkMomY;
extern double Bondi_SinkMomZ;
extern double Bondi_SinkMomXAbs;
extern double Bondi_SinkMomYAbs;
extern double Bondi_SinkMomZAbs;
extern double Bondi_SinkEk;
extern double Bondi_SinkEt;
extern int    Bondi_SinkNCell;
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Bondi
// Description :  Record the accretion rates of various quantities
//
// Note        :  1. Linked to the function pointer "Aux_Record_User_Ptr" by Init_TestProb_Hydro_Bondi()
//                2. Please turn on the runtime option "OPT__RECORD_USER"
//                3. This function only outputs the header when first calling it
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Record_Bondi()
{

   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__BondiAccretionRate", OUTPUT_DIR );

   static bool   FirstTime = true;
   static double Time0, dTime;

   if ( FirstTime )
   {
//    header
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )    Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#%9s%16s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n",
                  "Step", "Time [yr]", "NVoidCell", "Mass [Msun]", "Time [yr]", "dM/dt [Msun/yr]",
                  "MomX [g*cm/s]", "MomY [g*cm/s]", "MomZ [g*cm/s]", "MomXAbs [g*cm/s]", "MomYAbs [g*cm/s]", "MomZAbs [g*cm/s]",
                  "Ek [erg]", "Et [erg]", "BHMass [Msun]" );
         fclose( File_User );
      }

      FirstTime = false;
      Time0     = Time[0];
   }

   else
   {
//    get the total number of cells within the void region
      int SinkNCell_Sum;

      MPI_Reduce( &Bondi_SinkNCell, &SinkNCell_Sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );


//    get the total amount of sunk variables
      double Mass_Sum, MomX_Sum, MomY_Sum, MomZ_Sum, MomXAbs_Sum, MomYAbs_Sum, MomZAbs_Sum, Ek_Sum, Et_Sum, Mass_BH;

      MPI_Reduce( &Bondi_SinkMass,    &Mass_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkMomX,    &MomX_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkMomY,    &MomY_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkMomZ,    &MomZ_Sum,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkMomXAbs, &MomXAbs_Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkMomYAbs, &MomYAbs_Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkMomZAbs, &MomZAbs_Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkEk,      &Ek_Sum,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &Bondi_SinkEt,      &Et_Sum,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      Mass_Sum    *= UNIT_M/Const_Msun;
      MomX_Sum    *= UNIT_M*UNIT_L/UNIT_T;
      MomY_Sum    *= UNIT_M*UNIT_L/UNIT_T;
      MomZ_Sum    *= UNIT_M*UNIT_L/UNIT_T;
      MomXAbs_Sum *= UNIT_M*UNIT_L/UNIT_T;
      MomYAbs_Sum *= UNIT_M*UNIT_L/UNIT_T;
      MomZAbs_Sum *= UNIT_M*UNIT_L/UNIT_T;
      Ek_Sum      *= UNIT_E;
      Et_Sum      *= UNIT_E;
      Mass_BH      = Bondi_MassBH*UNIT_M/Const_Msun;

      if ( MPI_Rank == 0 )
      {
         dTime  = Time[0] - Time0;
         Time0  = Time[0];
         dTime *= UNIT_T/Const_yr;

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "%10ld%16.7e%20d%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e\n",
                  Step, Time[0]*UNIT_T/Const_yr, SinkNCell_Sum, Mass_Sum, dTime, Mass_Sum/dTime,
                  MomX_Sum, MomY_Sum, MomZ_Sum, MomXAbs_Sum, MomYAbs_Sum, MomZAbs_Sum, Ek_Sum, Et_Sum, Mass_BH );
         fclose( File_User );
      }
   } // if ( FirstTime ) ... else ...

// reset the cumulative variables to zero
   Bondi_SinkMass    = 0.0;
   Bondi_SinkMomX    = 0.0;
   Bondi_SinkMomY    = 0.0;
   Bondi_SinkMomZ    = 0.0;
   Bondi_SinkMomXAbs = 0.0;
   Bondi_SinkMomYAbs = 0.0;
   Bondi_SinkMomZAbs = 0.0;
   Bondi_SinkEk      = 0.0;
   Bondi_SinkEt      = 0.0;

} // FUNCTION : Record_Bondi



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
