#include "GAMER.h"




extern int      Merger_Coll_NumBHs;

extern double  *CM_BH_Mass;
extern double  *CM_BH_Mdot_tot;
extern double  *CM_BH_Mdot_hot;
extern double  *CM_BH_Mdot_cold;

extern double  *CM_Bondi_SinkMass;
extern double  *CM_Bondi_SinkMomX;
extern double  *CM_Bondi_SinkMomY;
extern double  *CM_Bondi_SinkMomZ;
extern double  *CM_Bondi_SinkMomXAbs;
extern double  *CM_Bondi_SinkMomYAbs;
extern double  *CM_Bondi_SinkMomZAbs;
extern double  *CM_Bondi_SinkE;
extern double  *CM_Bondi_SinkEk;
extern double  *CM_Bondi_SinkEt;
extern int     *CM_Bondi_SinkNCell;

extern double (*CM_RAcc_GasVel)[3];      // gas velocity
extern double  *CM_RAcc_SoundSpeed;
extern double  *CM_RAcc_GasDens;
extern double  *CM_RAcc_RelativeVel;
extern double  *CM_RAcc_ColdGasMass;
extern double (*CM_Jet_Vec)[3];          // jet direction
extern double  *CM_Jet_Mdot;             // the feedback injection rate
extern double  *CM_Jet_Pdot;
extern double  *CM_Jet_Edot;
extern double   E_inj_exp[3];
extern double   M_inj_exp[3];
static double   E_power_inj[3];          // the injection power
extern double   CM_ClusterCen[3][3];
extern double   CM_BH_Vel[3][3];
extern int     *CM_Cluster_NPar_close;   // total number of particles inside the target region of each cluster



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_ClusterMerger
// Description :  Record the cluster centers
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//                4. Must enable Merger_Coll_LabelCenter
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_ClusterMerger()
{

   const char FileName[] = "Record__ClusterCenter";
   static bool FirstTime = true;

// header
   if ( FirstTime )
   {
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist( FileName ) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#%19s%20s",  "Time", "Step" );
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "ClusterCen_x",     c, "ClusterCen_y",     c, "ClusterCen_z",      c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "BHVel_x[km/s]",    c, "BHVel_y",          c, "BHVel_z",           c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "GasVel_x",         c, "GasVel_y",         c, "GasVel_z",          c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "RelativeVel",      c, "SoundSpeed",       c, "GasDens(cgs)",      c );
            fprintf( File_User, " %19s%1d",                 "mass_BH[Msun]",    c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "Mdot_tot_BH(cgs)", c, "Mdot_hot_BH(cgs)", c, "Mdot_cold_BH(cgs)", c );
            fprintf( File_User, " %19s%1d",                 "NVoidCell",        c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "MomXInj(cgs)",     c, "MomYInj",          c, "MomZInj",           c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "MomXInjAbs",       c, "MomYInjAbs",       c, "MomZInjAbs",        c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "EInj_exp[erg]",    c, "E_Inj[erg]",       c, "E_Inj_err",         c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "Ek_Inj[erg]",      c, "Et_Inj[erg]",      c, "PowerInj(cgs)",     c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "MInjexp[Msun]",    c, "MassInj[Msun]",    c, "M_Inj_err",         c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "Mdot(cgs)",        c, "Pdot(cgs)",        c, "Edot(cgs)",         c );
            fprintf( File_User, " %19s%1d %19s%1d %19s%1d", "Jet_Vec_x",        c, "Jet_Vec_y",        c, "Jet_Vec_z",         c );
            fprintf( File_User, " %19s%1d %19s%1d",         "num_par_sum",      c, "ColdGasMass",      c );
         }
         fprintf( File_User, "\n" );
         fclose( File_User );
      } // if ( MPI_Rank == 0 )

      FirstTime = false;
   } // if ( FirstTime )

// sum over the variables and convert units
   int SinkNCell_Sum[3];
   double Mass_Sum[3], MomX_Sum[3], MomY_Sum[3], MomZ_Sum[3], MomXAbs_Sum[3], MomYAbs_Sum[3], MomZAbs_Sum[3], E_Sum[3], Ek_Sum[3], Et_Sum[3];

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      MPI_Reduce( &CM_Bondi_SinkNCell[c],   &SinkNCell_Sum[c], 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMass[c],    &Mass_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomX[c],    &MomX_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomY[c],    &MomY_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomZ[c],    &MomZ_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomXAbs[c], &MomXAbs_Sum[c],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomYAbs[c], &MomYAbs_Sum[c],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomZAbs[c], &MomZAbs_Sum[c],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkE[c],       &E_Sum[c],         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkEk[c],      &Ek_Sum[c],        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkEt[c],      &Et_Sum[c],        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      Mass_Sum[c]    *= UNIT_M / Const_Msun;
      MomX_Sum[c]    *= UNIT_M * UNIT_V;
      MomY_Sum[c]    *= UNIT_M * UNIT_V;
      MomZ_Sum[c]    *= UNIT_M * UNIT_V;
      MomXAbs_Sum[c] *= UNIT_M * UNIT_V;
      MomYAbs_Sum[c] *= UNIT_M * UNIT_V;
      MomZAbs_Sum[c] *= UNIT_M * UNIT_V;
      E_Sum[c]       *= UNIT_E;
      Ek_Sum[c]      *= UNIT_E;
      Et_Sum[c]      *= UNIT_E;
      E_inj_exp[c]   *= UNIT_E;
      M_inj_exp[c]   *= UNIT_M / Const_Msun;
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)

   for (int c=0; c<Merger_Coll_NumBHs; c++)   E_power_inj[c] = E_Sum[c]/(dTime_Base*UNIT_T);

// output the properties of the cluster centers
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%20.7e%20ld", Time[0], Step );
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_ClusterCen[c][0], CM_ClusterCen[c][1], CM_ClusterCen[c][2] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_BH_Vel[c][0]*UNIT_V/(Const_km/Const_s), CM_BH_Vel[c][1]*UNIT_V/(Const_km/Const_s), CM_BH_Vel[c][2]*UNIT_V/(Const_km/Const_s) );
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_RAcc_GasVel[c][0]*UNIT_V/(Const_km/Const_s), CM_RAcc_GasVel[c][1]*UNIT_V/(Const_km/Const_s), CM_RAcc_GasVel[c][2]*UNIT_V/(Const_km/Const_s) );
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_RAcc_RelativeVel[c]*UNIT_V/(Const_km/Const_s), CM_RAcc_SoundSpeed[c]*UNIT_V/(Const_km/Const_s), CM_RAcc_GasDens[c]*UNIT_D );
         fprintf( File_User, " %20.7e",               CM_BH_Mass[c]*UNIT_M/Const_Msun );
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_BH_Mdot_tot[c]*UNIT_M/UNIT_T, CM_BH_Mdot_hot[c]*UNIT_M/UNIT_T, CM_BH_Mdot_cold[c]*UNIT_M/UNIT_T );
         fprintf( File_User, " %20d",                 SinkNCell_Sum[c] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", MomX_Sum[c], MomY_Sum[c], MomZ_Sum[c] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", MomXAbs_Sum[c], MomYAbs_Sum[c], MomZAbs_Sum[c] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", E_inj_exp[c], E_Sum[c], (E_Sum[c]-E_inj_exp[c])/E_inj_exp[c] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", Ek_Sum[c], Et_Sum[c], E_power_inj[c] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", M_inj_exp[c], Mass_Sum[c], (Mass_Sum[c]-M_inj_exp[c])/M_inj_exp[c] );
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_Jet_Mdot[c]*UNIT_M/UNIT_T, CM_Jet_Pdot[c]*UNIT_M*UNIT_V/UNIT_T, CM_Jet_Edot[c]*UNIT_E/UNIT_T );
         fprintf( File_User, " %20.7e %20.7e %20.7e", CM_Jet_Vec[c][0], CM_Jet_Vec[c][1], CM_Jet_Vec[c][2] );
         fprintf( File_User, " %20d %20.7e",          CM_Cluster_NPar_close[c], CM_RAcc_ColdGasMass[c]*UNIT_M/Const_Msun );
      }
      fprintf( File_User, "\n" );
      fclose( File_User );
   } // if ( MPI_Rank == 0 )

// reset the cumulative variables to zero
   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      CM_Bondi_SinkMass[c]    = 0.0;
      CM_Bondi_SinkMomX[c]    = 0.0;
      CM_Bondi_SinkMomY[c]    = 0.0;
      CM_Bondi_SinkMomZ[c]    = 0.0;
      CM_Bondi_SinkMomXAbs[c] = 0.0;
      CM_Bondi_SinkMomYAbs[c] = 0.0;
      CM_Bondi_SinkMomZAbs[c] = 0.0;
      CM_Bondi_SinkE[c]       = 0.0;
      CM_Bondi_SinkEk[c]      = 0.0;
      CM_Bondi_SinkEt[c]      = 0.0;
      E_inj_exp[c]            = 0.0;
      M_inj_exp[c]            = 0.0;
   }

} // FUNCTION : Aux_Record_ClusterMerger
