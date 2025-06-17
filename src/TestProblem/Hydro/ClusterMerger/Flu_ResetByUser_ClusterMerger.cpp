#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )


extern int        Merger_Coll_NumBHs, Accretion_Mode;
extern double     eta, eps_f, eps_m, R_acc, R_dep;         // parameters of jet feedback
extern bool       AGN_feedback;
extern FieldIdx_t Idx_ParHalo;

extern double     CM_Bondi_SinkMass[3];
extern double     CM_Bondi_SinkMomX[3];
extern double     CM_Bondi_SinkMomY[3];
extern double     CM_Bondi_SinkMomZ[3];
extern double     CM_Bondi_SinkMomXAbs[3];
extern double     CM_Bondi_SinkMomYAbs[3];
extern double     CM_Bondi_SinkMomZAbs[3];
extern double     CM_Bondi_SinkE[3];
extern double     CM_Bondi_SinkEk[3];
extern double     CM_Bondi_SinkEt[3];
extern int        CM_Bondi_SinkNCell[3];

extern double     Bondi_MassBH1;
extern double     Bondi_MassBH2;
extern double     Bondi_MassBH3;
extern double     Mdot_tot_BH1;                            // the total accretion rate
extern double     Mdot_tot_BH2;
extern double     Mdot_tot_BH3;
extern double     Mdot_hot_BH1;                            // the hot   accretion rate
extern double     Mdot_hot_BH2;
extern double     Mdot_hot_BH3;
extern double     Mdot_cold_BH1;                           // the cold  accretion rate
extern double     Mdot_cold_BH2;
extern double     Mdot_cold_BH3;
extern double     Jet_HalfHeight1;
extern double     Jet_HalfHeight2;
extern double     Jet_HalfHeight3;
extern double     Jet_Radius1;
extern double     Jet_Radius2;
extern double     Jet_Radius3;
extern double     Jet_Vec[3][3];                           // jet direction
extern double     Mdot[3];                                 // the feedback injection rate
extern double     Pdot[3];
extern double     Edot[3];
extern double     GasVel[3][3];                            // gas velocity
extern double     SoundSpeed[3];
extern double     GasDens[3];
extern double     RelativeVel[3];                          // the relative velocity between BH and gas
extern double     ColdGasMass[3];
extern double     GasMass[3];
extern double     ParMass[3];
extern double     ClusterCen[3][3];
extern double     BH_Pos[3][3];                            // BH position (for updating ClusterCen)
extern double     BH_Vel[3][3];                            // BH velocity

static double     Jet_WaveK[3];                            // jet wavenumber used in the sin() function to have smooth bidirectional jets
static double     Jet_HalfHeight[3];
static double     Jet_Radius[3];
static double     V_cyl[3];                                // the volume of jet source
static double     M_inj[3], P_inj[3], E_inj[3];            // the injected density
static double     normalize_const[3];                      // the exact normalization constant

// the variables that need to be recorded
       double     E_inj_exp[3] = { 0.0, 0.0, 0.0 };        // the expected amount of injected energy
       double     M_inj_exp[3] = { 0.0, 0.0, 0.0 };        // the expected amount of injected gas mass

static bool       FirstTime = true;
extern int        JetDirection_NBin;                       // number of bins of the jet direction table
extern double    *Time_table;                              // the time table of jet direction
extern double    *Theta_table[3];                          // the theta table of jet direction for 3 clusters
extern double    *Phi_table[3];                            // the phi table of jet direction for 3 clusters

extern bool       AdjustBHPos;
extern bool       AdjustBHVel;
extern double     AdjustPeriod;
extern int        AdjustCount;                             // count the number of adjustments
extern int        JetDirection_case;                       // methods for choosing the jet direction
       int        merge_index = 0;                         // record BH 1 merge BH 2 / BH 2 merge BH 1

       double     ang_mom_sum[3][3] = { { 1.0, 0.0, 0.0 },
                                        { 1.0, 0.0, 0.0 },
                                        { 1.0, 0.0, 0.0 } };
       bool       if_overlap = false;

extern void GetClusterCenter( int lv, bool AdjustPos, bool AdjustVel, double Cen_old[][3], double Cen_new[][3], double Cen_Vel[][3] );

#ifdef MHD
extern double (*MHD_ResetByUser_VecPot_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[] );
extern double (*MHD_ResetByUser_BField_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[],
                                             const double B_in, const bool UseVecPot, const real *Ax, const real *Ay,
                                             const real *Az, const int i, const int j, const int k );
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  BH_accretion_rate
// Description :  Calculate the black hole accretion rate
//
// Note        :
//
// Parameter   :  mode         : Accretion mode
//                Mdot_tot     : Pointer to store the total accretion rate
//                Mdot_hot     : Pointer to store the hot   accretion rate
//                Mdot_cold    : Pointer to store the cold  accretion rate
//                r_acc        : Accretion radius of the black hole
//                mass_BH      : Mass of the black hole
//                rho_gas      : Average gas density within the accretion radius
//                cs           : Average sound speed within the accretion radius
//                v_rel        : Relative velocity of the gas within the accretion radius and the black hole
//                mass_coldGas : Mass of the cold gas within the accretion radius
//                mass_gas     : Mass of the gas      within the accretion radius
//                mass_par     : Mass of the particle within the accretion radius
//
// Return      :  Mdot_tot  : Total accretion rate
//                Mdot_hot  : Hot   accretion rate
//                Mdot_cold : Cold  accretion rate
//-------------------------------------------------------------------------------------------------------
void BH_accretion_rate( const int mode, double *Mdot_tot, double *Mdot_hot, double *Mdot_cold, const double r_acc,
                        const double mass_BH, const double rho_gas, const double cs, const double v_rel,
                        const double mass_coldGas, const double mass_gas, const double mass_par )
{

   double acc_cold = 0.0, acc_hot = 0.0;

// hot accretion rate
   if ( (mode == 1  ||  mode == 3)  &&  rho_gas > 0.0 )
   {
      acc_hot = 4.0 * M_PI * SQR(NEWTON_G) * SQR(mass_BH) * rho_gas /
                pow( SQR(cs) + SQR(v_rel), 1.5 );
   }

// cold accretion rate
   if ( ( mode == 2  ||  mode == 3 )  &&  (mass_gas+mass_par) > 0.0 )
   {
      const double t_ff = sqrt( 2*CUBE(r_acc) / NEWTON_G / (mass_gas+mass_par) );
      acc_cold = mass_coldGas / t_ff;
   }

   *Mdot_tot  = acc_hot + acc_cold;
   *Mdot_hot  = acc_hot;
   *Mdot_cold = acc_cold;

} // FUNCTION : BH_accretion_rate



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func_ClusterMerger
// Description :  Function to reset the fluid field in the Bondi accretion problem
//
// Note        :  1. Invoked by Flu_ResetByUser_API_ClusterMerger() and Hydro_Init_ByFunction_AssignData() using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                   --> This function pointer is reset by Init_TestProb_Hydro_ClusterMerger()
//                   --> Hydro_Init_ByFunction_AssignData(): constructing initial condition
//                       Flu_ResetByUser_API_ClusterMerger()       : after each update
//                2. Input fluid[] stores the original values
//                3. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically
//                4. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                dt       : Time interval to advance solution
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
int Flu_ResetByUser_Func_ClusterMerger( real fluid[], const double Emag, const double x, const double y, const double z,
                                        const double Time, const double dt, const int lv, double AuxArray[] )
{

   const double Pos[3]  = { x, y, z };

// (1) SMBH Accretion (note: gas depletion is currently disabled)
   double dr2[3][3], r2[3];
   const double V_dep = 4.0 / 3.0 * M_PI * pow( R_dep, 3.0 ); // the volume to remove gas
// the density need to be removed
   double D_dep[3] = { Mdot_tot_BH1*dt/V_dep, Mdot_tot_BH2*dt/V_dep, Mdot_tot_BH3*dt/V_dep };
   int    reset   = false; // mark whether this cell is reset or not [false/true]

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      for (int d=0; d<3; d++)   dr2[c][d] = SQR( Pos[d] - ClusterCen[c][d] );
      r2[c] = dr2[c][0] + dr2[c][1] + dr2[c][2];
      if ( r2[c] <= SQR(R_dep) )
      {

           double dens_old = fluid[DENS];
           fluid[DENS] -= D_dep[c];

           if ( fluid[DENS] < MIN_DENS )  fluid[DENS] = MIN_DENS;

           fluid[MOMX] *= fluid[DENS]/dens_old;
           fluid[MOMY] *= fluid[DENS]/dens_old;
           fluid[MOMZ] *= fluid[DENS]/dens_old;
           fluid[ENGY] *= fluid[DENS]/dens_old;
           reset = true;

      } // if ( r2[c] <= SQR(R_dep) )
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)


// (2) Jet Feedback
   double Jet_dr, Jet_dh;
   double Dis_c2m, Vec_c2m[3];
   real   EngySin;
   int    status, n_jet;
   int    which_cluster = 0;

   if ( if_overlap )
   {
      if ( Merger_Coll_NumBHs == 1 )  Aux_Error( ERROR_INFO, "Error: Merger_Coll_NumBHs = 1 but if_overlap = true!\n" );
      if ( Edot[0] >= Edot[1] )   status = 0;   // only inject cluster 1
      else                        status = 1;   // only inject cluster 2
      n_jet = Merger_Coll_NumBHs-1;
   }
   else // if ( if_overlap )
   {
      status = 0;
      n_jet = Merger_Coll_NumBHs;
   } // if ( if_overlap ) ...  else ...

   for (int c=status; c<(n_jet+status); c++)
   {

//    distance: jet center to mesh
      for (int d=0; d<3; d++)   Vec_c2m[d] = Pos[d] - ClusterCen[c][d];

      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );

      Jet_dh = fabs( Jet_Vec[c][0]*Vec_c2m[0] + Jet_Vec[c][1]*Vec_c2m[1] + Jet_Vec[c][2]*Vec_c2m[2] );
      Jet_dr = sqrt( SQR(Dis_c2m) - SQR(Jet_dh) );

      if ( Jet_dh <= Jet_HalfHeight[c]  &&  Jet_dr <= Jet_Radius[c] )
      {
         which_cluster += c+1;

//       record the old momentum
         const real momx_old = fluid[MOMX];
         const real momy_old = fluid[MOMY];
         const real momz_old = fluid[MOMZ];
         const real dens_old = fluid[DENS];

//       accrete mass
         fluid[DENS] += M_inj[c];

//       transfer into BH frame
         fluid[MOMX] -= BH_Vel[c][0]*dens_old;
         fluid[MOMY] -= BH_Vel[c][1]*dens_old;
         fluid[MOMZ] -= BH_Vel[c][2]*dens_old;

//       use a sine function to make the velocity smooth within the jet from +Jet_Vec to -Jet_Vec
         EngySin = E_inj[c]*normalize_const[c]*sin( Jet_WaveK[c]*Jet_dh );
         real P_SQR    = SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]);
//       the new momentum is calculated from the old density, new density, old momentum and injected energy
         real Ekin_old  = ( dens_old == 0.0 ) ? (real)0.0 : 0.5*P_SQR/dens_old;
         real P_new     = sqrt( 2*fluid[DENS]*(EngySin + Ekin_old) );
         P_new *= SIGN( Vec_c2m[0]*Jet_Vec[c][0] + Vec_c2m[1]*Jet_Vec[c][1] + Vec_c2m[2]*Jet_Vec[c][2] );
         fluid[MOMX] = P_new * Jet_Vec[c][0];
         fluid[MOMY] = P_new * Jet_Vec[c][1];
         fluid[MOMZ] = P_new * Jet_Vec[c][2];

//       transfer back into the rest frame
         fluid[MOMX] += BH_Vel[c][0] * fluid[DENS];
         fluid[MOMY] += BH_Vel[c][1] * fluid[DENS];
         fluid[MOMZ] += BH_Vel[c][2] * fluid[DENS];

         fluid[ENGY] += 0.5*( (SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ])) / fluid[DENS] -
                              (SQR(momx_old)    + SQR(momy_old)    + SQR(momz_old)   ) / dens_old );
      } // if ( Jet_dh <= Jet_HalfHeight[c]  &&  Jet_dr <= Jet_Radius[c] )
   } // for (int c=status; c<(n_jet+status); c++)

   if ( which_cluster >= 3 )   Aux_Error( ERROR_INFO, "Error: which_cluster >= 3!\n" );
   return which_cluster;

} // FUNCTION : Flu_ResetByUser_Func_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API_ClusterMerger
// Description :  API for resetting the fluid array in the Bondi accretion problem
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked using the function pointer "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer is reset by Init_TestProb_Hydro_ClusterMerger()
//                3. Currently does not work with "OPT__OVERLAP_MPI"
//                4. Invoke Flu_ResetByUser_Func_ClusterMerger() directly
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Target fluid sandglass
//                MagSg   : Target B field sandglass
//                TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_ClusterMerger( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt )
{

   const bool    CurrentMaxLv = ( NPatchTotal[lv] > 0  &&  lv == MAX_LEVEL        ) ? true :
                                ( NPatchTotal[lv] > 0  &&  NPatchTotal[lv+1] == 0 ) ? true : false;
   const double  dh           = amr->dh[lv];
   const real    dv           = CUBE(dh);
#  if ( MODEL == HYDRO  &&  !defined SRHD )
   const real    Gamma_m1     = GAMMA - (real)1.0;
   const real   _Gamma_m1     = (real)1.0 / Gamma_m1;
#  endif

// (1) get the BH position and velocity and adjust them if needed
   bool AdjustPosNow = false;
   bool AdjustVelNow = false;

// only adjust the BHs on the current maximum level
   if ( CurrentMaxLv  &&  AdjustCount < int(TimeNew/AdjustPeriod) )
   {
      if ( AdjustBHPos == true )   AdjustPosNow = true;
      if ( AdjustBHVel == true )   AdjustVelNow = true;
      AdjustCount += 1;
   } // if ( CurrentMaxLv  &&  AdjustCount < int(TimeNew/AdjustPeriod))
   GetClusterCenter( lv, AdjustPosNow, AdjustVelNow, BH_Pos, ClusterCen, BH_Vel );


// (2) decide whether to merge BHs
   double RelativeBHPos[3] = { BH_Pos[0][0]-BH_Pos[1][0], BH_Pos[0][1]-BH_Pos[1][1], BH_Pos[0][2]-BH_Pos[1][2] };
   double RelativeBHVel[3] = { BH_Vel[0][0]-BH_Vel[1][0], BH_Vel[0][1]-BH_Vel[1][1], BH_Vel[0][2]-BH_Vel[1][2] };
   double AbsRelPos        = sqrt( SQR(RelativeBHPos[0]) + SQR(RelativeBHPos[1]) + SQR(RelativeBHPos[2]) );
   double AbsRelVel        = sqrt( SQR(RelativeBHVel[0]) + SQR(RelativeBHVel[1]) + SQR(RelativeBHVel[2]) );
   double escape_vel;
   double soften           = amr->dh[MAX_LEVEL];
   if ( AbsRelPos > soften )
   {
      escape_vel = sqrt( 2 * NEWTON_G * (Bondi_MassBH1+Bondi_MassBH2) / AbsRelPos );
   }
   else
   {
      escape_vel = sqrt( 2 * NEWTON_G * (Bondi_MassBH1+Bondi_MassBH2) / soften );
   } // if ( AbsRelPos > soften ) ... else ...

// merge the two BHs if they are located within R_acc, and the relative velocity is small enough
   if ( Merger_Coll_NumBHs == 2 )
   {
      if ( AbsRelPos < R_acc  &&  AbsRelVel < 3*escape_vel )
      {
         Merger_Coll_NumBHs -= 1;
         if ( Bondi_MassBH1 >= Bondi_MassBH2 )   merge_index = 1;   // record BH 1 merge BH 2 / BH 2 merge BH 1
         else                                    merge_index = 2;
         Bondi_MassBH1 += Bondi_MassBH2;
         Bondi_MassBH2 = 0.0;
//       relabel the BH and DM particles being merged
         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
            if ( amr->Par->AttributeInt[Idx_ParHalo][p] != (long_par)1 )   continue;
            if ( amr->Par->Mass[p] < (real)0.0 )   continue;

            amr->Par->Type[p] = PTYPE_DARK_MATTER;
         } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         Aux_Message( stdout, "BHs Merge! In rank %d, TimeNew = %14.8e; merge_index = %d, "
                              "BHPos1 = %14.8e, %14.8e, %14.8e; BHPos2 = %14.8e, %14.8e, %14.8e; "
                              "BHVel1 = %14.8e, %14.8e, %14.8e; BHVel2 = %14.8e, %14.8e, %14.8e; "
                              "AbsRelPos = %14.8e, AbsRelVel = %14.8e, escape_vel = %14.8e.\n",
                      MPI_Rank, TimeNew, merge_index,
                      BH_Pos[0][0], BH_Pos[0][1], BH_Pos[0][2], BH_Pos[1][0], BH_Pos[1][1], BH_Pos[1][2],
                      BH_Vel[0][0], BH_Vel[0][1], BH_Vel[0][2], BH_Vel[1][0], BH_Vel[1][1], BH_Vel[1][2],
                      AbsRelPos, AbsRelVel, escape_vel );
      } // if ( AbsRelPos < R_acc  &&  AbsRelVel < 3*escape_vel )
   } // if ( Merger_Coll_NumBHs == 2 )


   if ( AGN_feedback )
   {
//    (3) set the injection parameters
      double Mdot_tot_BH [3] = { Mdot_tot_BH1,  Mdot_tot_BH2,  Mdot_tot_BH3  };
      double Mdot_hot_BH [3] = { Mdot_hot_BH1,  Mdot_hot_BH2,  Mdot_hot_BH3  };
      double Mdot_cold_BH[3] = { Mdot_cold_BH1, Mdot_cold_BH2, Mdot_cold_BH3 };
      double Bondi_MassBH[3] = { Bondi_MassBH1, Bondi_MassBH2, Bondi_MassBH3 };

      Jet_HalfHeight[0] = Jet_HalfHeight1;
      Jet_HalfHeight[1] = Jet_HalfHeight2;
      Jet_HalfHeight[2] = Jet_HalfHeight3;
      Jet_Radius[0]     = Jet_Radius1;
      Jet_Radius[1]     = Jet_Radius2;
      Jet_Radius[2]     = Jet_Radius3;

//    set the jet direction vector
      if ( JetDirection_case == 1 ) // fixed at x-axis
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            Jet_Vec[c][0] = 1.0;
            Jet_Vec[c][1] = 0.0;
            Jet_Vec[c][2] = 0.0;
         }
      }
      else if ( JetDirection_case == 2 ) // import from table
      {
         const double Time_period      = Time_table[JetDirection_NBin-1];
         const double Time_interpolate = fmod( TimeNew, Time_period );
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            const double theta = Mis_InterpolateFromTable( JetDirection_NBin, Time_table, Theta_table[c], Time_interpolate );
            const double phi   = Mis_InterpolateFromTable( JetDirection_NBin, Time_table, Phi_table[c],   Time_interpolate );

            Jet_Vec[c][0] = cos(theta);
            Jet_Vec[c][1] = sin(theta)*cos(phi);
            Jet_Vec[c][2] = sin(theta)*sin(phi);
         }
      }
      else if ( JetDirection_case == 3 ) // align with angular momentum
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            const double ang_mom_norm = sqrt( SQR(ang_mom_sum[c][0]) + SQR(ang_mom_sum[c][1]) + SQR(ang_mom_sum[c][2]) );
            for (int d=0; d<3; d++)   Jet_Vec[c][d] = ang_mom_sum[c][d] / ang_mom_norm;
         }
      } // if ( JetDirection_case == 1 ) ... else ...


//    (4) calculate the accretion and feedback
      real   fluid[NCOMP_TOTAL], fluid_old[NCOMP_TOTAL], fluid_acc[NCOMP_TOTAL];
      double x, y, z, x2, y2, z2, x02, y02, z02;

//    reset to 0 since we only want to record the number of void cells **for one sub-step**
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         CM_Bondi_SinkNCell[c] = 0;
         Jet_WaveK[c]          = 0.5 * M_PI / Jet_HalfHeight[c];
      }

//    variables for each rank
      int  num        [3]       =  { 0, 0, 0 };        // the number of cells inside the accretion radius
      double gas_mass [3]       =  { 0.0, 0.0, 0.0 };  // total gas mass inside the accretion radius
      real rho        [3]       =  { 0.0, 0.0, 0.0 };  // the average density inside the accretion radius (hot gas)
      double mass_cold[3]       =  { 0.0, 0.0, 0.0 };  // cold gas mass (T < 5e5 K) inside the accretion radius
      real Cs         [3]       =  { 0.0, 0.0, 0.0 };  // the average sound speed inside the accretion radius
      real gas_mom    [3][3]    = {{ 0.0, 0.0, 0.0 },  // average gas momentum
                                   { 0.0, 0.0, 0.0 },
                                   { 0.0, 0.0, 0.0 }};
      double ang_mom  [3][3]    = {{ 0.0, 0.0, 0.0 },  // total angular momentum inside the accretion radius
                                   { 0.0, 0.0, 0.0 },
                                   { 0.0, 0.0, 0.0 }};
      real V_cyl_exact[3]       =  { 0.0, 0.0, 0.0 };  // exact volume of jet cylinder
      real normalize  [3]       =  { 0.0, 0.0, 0.0 };  // for computing the correct normalization constant
      bool if_overlap_each_rank = false;

//    variables for all ranks
      int    num_sum[3];
      real rho_sum[3], Cs_sum[3], gas_mom_sum[3][3], gas_vel_sum[3][3], V_cyl_exact_sum[3], normalize_sum[3];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         x02 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         y02 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         z02 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)   { z2 = z02 + k*dh;
         for (int j=0; j<PS1; j++)   { y2 = y02 + j*dh;
         for (int i=0; i<PS1; i++)   { x2 = x02 + i*dh;

            for (int v=0; v<NCOMP_TOTAL; v++)   fluid_acc[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#           else
            const real Emag = (real)0.0;
#           endif

            int status_overlap = 0;

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
//             calculate the average density, sound speed and gas velocity inside accretion radius
               if ( SQR(x2-ClusterCen[c][0]) + SQR(y2-ClusterCen[c][1]) + SQR(z2-ClusterCen[c][2]) <= SQR(R_acc) )
               {
                  gas_mass[c] += fluid_acc[0]*dv;
#                 ifdef DUAL_ENERGY
                  const real Pres = Hydro_DensDual2Pres( fluid_acc[0], fluid_acc[DUAL], EoS_AuxArray_Flt[1],
                                                         false, NULL_REAL );
                  const real Eint = EoS_DensPres2Eint_CPUPtr( fluid_acc[0], Pres, NULL, EoS_AuxArray_Flt,
                                                              EoS_AuxArray_Int, h_EoS_Table );
                  const real Temp = EoS_DensEint2Temp_CPUPtr( fluid_acc[0], Eint, NULL, EoS_AuxArray_Flt,
                                                              EoS_AuxArray_Int, h_EoS_Table );
#                 else
                  const real Pres = Hydro_Con2Pres( fluid_acc[0], fluid_acc[1], fluid_acc[2], fluid_acc[3],
                                                    fluid_acc[4], fluid_acc+NCOMP_FLUID, true, MIN_PRES, Emag,
                                                    EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                                    EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                    h_EoS_Table, NULL );
                  const real Temp = Hydro_Con2Temp( fluid_acc[0], fluid_acc[1], fluid_acc[2], fluid_acc[3],
                                                    fluid_acc[4], fluid_acc+NCOMP_FLUID, false, MIN_TEMP, Emag,
                                                    EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                                    EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                    h_EoS_Table );
#                 endif
                  if ( Temp <= 5e5 )   mass_cold[c] += fluid_acc[0]*dv;
                  else
                  {
                     rho[c] += fluid_acc[0]*dv;
                     Cs[c] += sqrt( EoS_DensPres2CSqr_CPUPtr( fluid_acc[0], Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table ) );
                     for (int d=0; d<3; d++)   gas_mom[c][d] += fluid_acc[d+1]*dv;
                     num[c] += 1;
                  } // if ( Temp <= 5e5 )
                  double dr[3] = { x2-ClusterCen[c][0], y2-ClusterCen[c][1], z2-ClusterCen[c][2] };
                  ang_mom[c][0] += dv * ( dr[1]*fluid_acc[3] - dr[2]*fluid_acc[2] );
                  ang_mom[c][1] += dv * ( dr[2]*fluid_acc[1] - dr[0]*fluid_acc[3] );
                  ang_mom[c][2] += dv * ( dr[0]*fluid_acc[2] - dr[1]*fluid_acc[1] );
               } // if ( SQR(x2-ClusterCen[c][0]) + SQR(y2-ClusterCen[c][1]) + SQR(z2-ClusterCen[c][2]) <= SQR(R_acc) )

//             calculate the exact volume of jet cylinder and normalization
               if ( CurrentMaxLv )
               {
                  double Jet_dr_2, Jet_dh_2, Dis_c2m_2, Vec_c2m_2[3];
                  double Pos_2[3] = {x2, y2, z2};

                  for (int d=0; d<3; d++)   Vec_c2m_2[d] = Pos_2[d] - ClusterCen[c][d];
                  Dis_c2m_2 = sqrt( SQR(Vec_c2m_2[0]) + SQR(Vec_c2m_2[1]) + SQR(Vec_c2m_2[2]) );
                  Jet_dh_2 = fabs( Jet_Vec[c][0]*Vec_c2m_2[0] + Jet_Vec[c][1]*Vec_c2m_2[1] + Jet_Vec[c][2]*Vec_c2m_2[2] );
                  Jet_dr_2 = sqrt( SQR(Dis_c2m_2) - SQR(Jet_dh_2) );

                  if ( Jet_dh_2 <= Jet_HalfHeight[c]  &&  Jet_dr_2 <= Jet_Radius[c] )
                  {
                     V_cyl_exact[c] += dv;
                     normalize[c]   += sin( Jet_WaveK[c]*Jet_dh_2 )*dv;
                     status_overlap += c+1;
                  } // if ( Jet_dh_2 <= Jet_HalfHeight[c]  &&  Jet_dr_2 <= Jet_Radius[c] )
               } // if ( CurrentMaxLv )
            } // for (int c=0; c<Merger_Coll_NumBHs; c++)

            if ( status_overlap >= 3 )   if_overlap_each_rank = true;

         }}} // for (int k=0; k<PS1; k++); for (int j=0; j<PS1; j++); for (int i=0; i<PS1; i++)
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         MPI_Allreduce( &num[c],         &num_sum[c],         1, MPI_INT,        MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &gas_mass[c],    &GasMass[c],         1, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &rho[c],         &rho_sum[c],         1, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &mass_cold[c],   &ColdGasMass[c],     1, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &Cs[c],          &Cs_sum[c],          1, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( gas_mom[c],       gas_mom_sum[c],     3, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( ang_mom[c],       ang_mom_sum[c],     3, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &V_cyl_exact[c], &V_cyl_exact_sum[c], 1, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &normalize[c],   &normalize_sum[c],   1, MPI_GAMER_REAL, MPI_SUM, MPI_COMM_WORLD );
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)
      MPI_Allreduce( &if_overlap_each_rank, &if_overlap, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD );

//    calculate the total DM mass inside the accretion radius
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         if ( Accretion_Mode == 2  ||  Accretion_Mode == 3 )
         {
            double par_mass = 0.0;
            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               const double *EdgeL        = amr->patch[0][lv][PID]->EdgeL;
               const double *EdgeR        = amr->patch[0][lv][PID]->EdgeR;
               const double  patch_pos[3] = { (EdgeL[0]+EdgeR[0])*0.5, (EdgeL[1]+EdgeR[1])*0.5, (EdgeL[2]+EdgeR[2])*0.5 };
               const double  patch_d      = sqrt( SQR(EdgeL[0]-EdgeR[0]) + SQR(EdgeL[1]-EdgeR[1]) + SQR(EdgeL[2]-EdgeR[2]) ) * 0.5;
               if ( DIST_SQR_3D( patch_pos, ClusterCen[c] ) <= SQR(2*R_acc+patch_d) )
               {
                  for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
                  {
                     const long ParID = amr->patch[0][lv][PID]->ParList[p];
                     const real_par ParM  = amr->Par->Mass[ParID];
                     const real_par ParPos[3] = { amr->Par->PosX[ParID], amr->Par->PosY[ParID], amr->Par->PosZ[ParID] };
                     if ( DIST_SQR_3D( ParPos, ClusterCen[c] ) <= SQR(R_acc) )
                        par_mass += ParM;
                  } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               } // if ( DIST_SQR_3D( patch_pos, ClusterCen[c] ) <= SQR(2*R_acc+patch_d) )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            MPI_Allreduce( &par_mass, &ParMass[c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         }
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//    calculate the accretion rate

      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         double v = 0.0;  // the relative velocity between BH and gas
         if ( num_sum[c] == 0 )
         {
            GasDens[c]     = 0.0;
            SoundSpeed[c]  = 0.0;
            RelativeVel[c] = 0.0;
            for (int d=0; d<3; d++)  GasVel[c][d] = 0.0;
         }
         else
         {
            for (int d=0; d<3; d++)  gas_vel_sum[c][d] = gas_mom_sum[c][d] / rho_sum[c];
            rho_sum[c] /= ( 4.0 / 3.0 * M_PI * pow(R_acc, 3) );
            Cs_sum[c]  /= (double)num_sum[c];
            for (int d=0; d<3; d++)   v += SQR( BH_Vel[c][d] - gas_vel_sum[c][d] );

            GasDens[c]     = rho_sum[c];
            SoundSpeed[c]  = Cs_sum[c];
            RelativeVel[c] = sqrt(v);
            for (int d=0; d<3; d++)  GasVel[c][d] = gas_vel_sum[c][d];
         } // if ( num_sum[c] == 0 ) ... else ...

         BH_accretion_rate( Accretion_Mode, Mdot_tot_BH+c, Mdot_hot_BH+c, Mdot_cold_BH+c,
                            R_acc, Bondi_MassBH[c], GasDens[c], SoundSpeed[c], RelativeVel[c],
                            ColdGasMass[c], GasMass[c], ParMass[c] );

         if ( V_cyl_exact_sum[c] != 0 )   normalize_const[c] = V_cyl_exact_sum[c] / normalize_sum[c];
         else                             normalize_const[c] = 0.5 * M_PI;

      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

      Mdot_tot_BH1  = Mdot_tot_BH[0];
      Mdot_tot_BH2  = Mdot_tot_BH[1];
      Mdot_tot_BH3  = Mdot_tot_BH[2];

      Mdot_hot_BH1  = Mdot_hot_BH[0];
      Mdot_hot_BH2  = Mdot_hot_BH[1];
      Mdot_hot_BH3  = Mdot_hot_BH[2];

      Mdot_cold_BH1 = Mdot_cold_BH[0];
      Mdot_cold_BH2 = Mdot_cold_BH[1];
      Mdot_cold_BH3 = Mdot_cold_BH[2];

//    update BH mass
      for (int c=0; c<Merger_Coll_NumBHs; c++)   if ( CurrentMaxLv )   Bondi_MassBH[c] += Mdot_tot_BH[c] * dt;
      Bondi_MassBH1 = Bondi_MassBH[0];
      Bondi_MassBH2 = Bondi_MassBH[1];
      Bondi_MassBH3 = Bondi_MassBH[2];

//    (5) calculate the injection rate
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         Mdot[c]  = eta * Mdot_tot_BH[c];
         Pdot[c]  = sqrt(2*eta*eps_f*(1.0-eps_m)) * Mdot_tot_BH[c] * (Const_c/UNIT_V);
         Edot[c]  = eps_f * Mdot_tot_BH[c] * SQR(Const_c/UNIT_V);
         V_cyl[c] = M_PI * SQR(Jet_Radius[c]) * 2 * Jet_HalfHeight[c];

//       calculate the density that need to be injected
         if ( CurrentMaxLv  &&  V_cyl_exact_sum[c] != 0 )
         {
            M_inj[c] = Mdot[c] * dt / V_cyl_exact_sum[c];
            P_inj[c] = Pdot[c] * dt / V_cyl_exact_sum[c];
            E_inj[c] = Edot[c] * dt / V_cyl_exact_sum[c];
         }
         else
         {
            M_inj[c] = Mdot[c] * dt / V_cyl[c];
            P_inj[c] = Pdot[c] * dt / V_cyl[c];
            E_inj[c] = Edot[c] * dt / V_cyl[c];
         } // if ( CurrentMaxLv  &&  V_cyl_exact_sum[c] != 0 ) ... else ...
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

      if ( CurrentMaxLv )
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)   E_inj_exp[c] += Edot[c]*dt;
         for (int c=0; c<Merger_Coll_NumBHs; c++)   M_inj_exp[c] += Mdot[c]*dt;
      } // if ( CurrentMaxLv )


//    (6) perform injection
//    use the "static" schedule for reproducibility
#     pragma omp parallel for private( fluid, fluid_old, x, y, z ) schedule( static ) \
      reduction( +:CM_Bondi_SinkMass, CM_Bondi_SinkMomX, CM_Bondi_SinkMomY, CM_Bondi_SinkMomZ, \
                   CM_Bondi_SinkMomXAbs, CM_Bondi_SinkMomYAbs, CM_Bondi_SinkMomZAbs, \
                   CM_Bondi_SinkE, CM_Bondi_SinkEk, CM_Bondi_SinkEt, CM_Bondi_SinkNCell )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         int Reset;
         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               fluid   [v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//             backup the unmodified values since we want to record the amount of sunk variables removed at the maximum level
               fluid_old[v] = fluid[v];
            }

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#           else
            const real Emag = (real)0.0;
#           endif

//          reset this cell
            if ( CurrentMaxLv )   Reset = Flu_ResetByUser_Func_ClusterMerger( fluid, Emag, x, y, z, TimeNew, dt, lv, NULL );

//          operations necessary only when this cell has been reset
            if ( Reset != 0 )
            {
#              if ( MODEL == HYDRO  &&  !defined SRHD )
//             abort the program instead of applying a density floor
               if ( fluid[DENS] < MIN_DENS )   Aux_Error( ERROR_INFO, "Fluid density has reached the floor!\n" );

//             apply an energy floor
               fluid[ENGY] = Hydro_CheckMinEintInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                       (real)MIN_EINT, Emag );

//             calculate the dual-energy variable
#              ifdef DUAL_ENERGY
               fluid[DUAL] = Hydro_Con2Dual( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Emag,
                                             EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#              endif

//             floor and normalize passive scalars
#              if ( NCOMP_PASSIVE > 0 )
               for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)   fluid[v] = FMAX( fluid[v], TINY_NUMBER );

               if ( OPT__NORMALIZE_PASSIVE )
                  Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#              endif
#              endif // if ( MODEL == HYDRO  ||  MODEL == MHD )

//             record the amount of sunk variables removed at the maximum level
               if ( CurrentMaxLv )
               {
                  const real Ekin_old = 0.5*( SQR(fluid_old[MOMX]) + SQR(fluid_old[MOMY]) + SQR(fluid_old[MOMZ]) ) /
                                             (fluid_old[DENS]);
                  const real Ekin_new = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
                  const real Eint_old = fluid_old[ENGY] - Ekin_old - Emag;
                  const real Eint_new = fluid[ENGY]     - Ekin_new - Emag;

                  CM_Bondi_SinkMass[Reset-1]    += dv *     ( fluid[DENS] - fluid_old[DENS] );
                  CM_Bondi_SinkMomX[Reset-1]    += dv *     ( fluid[MOMX] - fluid_old[MOMX] );
                  CM_Bondi_SinkMomY[Reset-1]    += dv *     ( fluid[MOMY] - fluid_old[MOMY] );
                  CM_Bondi_SinkMomZ[Reset-1]    += dv *     ( fluid[MOMZ] - fluid_old[MOMZ] );
                  CM_Bondi_SinkMomXAbs[Reset-1] += dv * FABS( fluid[MOMX] - fluid_old[MOMX] );
                  CM_Bondi_SinkMomYAbs[Reset-1] += dv * FABS( fluid[MOMY] - fluid_old[MOMY] );
                  CM_Bondi_SinkMomZAbs[Reset-1] += dv * FABS( fluid[MOMZ] - fluid_old[MOMZ] );
                  CM_Bondi_SinkE[Reset-1]       += dv *     ( fluid[ENGY] - fluid_old[ENGY] );
                  CM_Bondi_SinkEk[Reset-1]      += dv *     ( Ekin_new    - Ekin_old        );
                  CM_Bondi_SinkEt[Reset-1]      += dv *     ( Eint_new    - Eint_old        );
                  CM_Bondi_SinkNCell[Reset-1]   ++;
               } // if ( CurrentMaxLv )

//             store the reset values
               for (int v=0; v<NCOMP_TOTAL; v++)
                  amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];

            } // if ( Reset != 0 )

         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // if ( AGN_feedback )

} // FUNCTION : Flu_ResetByUser_API_ClusterMerger



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )
