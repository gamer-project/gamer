#include "GAMER.h"



#if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )
// problem-specific global variables
// =======================================================================================
// (1) read from Input__TestProb
extern int        Merger_Coll_NumHalos;
extern double    *CM_BH_Mass;
extern double    *Jet_HalfHeight;
extern double    *Jet_Radius;

extern bool       AGN_feedback;
extern int        Accretion_Mode, JetDirection_case;
extern double     eta, eps_f, eps_m, R_acc, R_dep;
extern bool       fixBH, AdjustBHPos, AdjustBHVel;
extern double     AdjustPeriod;
// ---------------------------------------------------------------------------------------
// (2) read from files
extern int        JetDirection_NBin;                       // number of bins of the jet direction table
extern double    *CM_Jet_Time_table;                       // the time  table of jet direction
extern double   **CM_Jet_Theta_table;                      // the theta table of jet direction for 3 clusters
extern double   **CM_Jet_Phi_table;                        // the phi   table of jet direction for 3 clusters
// ---------------------------------------------------------------------------------------
// (3) variables to record
extern double   (*CM_ClusterCen)[3];
extern double   (*CM_BH_Pos)[3];                           // BH position (for updating CM_ClusterCen)
extern double   (*CM_BH_Vel)[3];                           // BH velocity
extern double    *CM_BH_Mdot_tot;                          // the total accretion rate
extern double    *CM_BH_Mdot_hot;                          // the hot   accretion rate
extern double    *CM_BH_Mdot_cold;                         // the cold  accretion rate
extern double    *CM_Jet_Mdot;                             // the feedback injection rate
extern double    *CM_Jet_Pdot;
extern double    *CM_Jet_Edot;
extern double   (*CM_Jet_Vec)[3];                          // jet direction
extern double   (*CM_RAcc_GasVel)[3];                      // gas velocity
extern double    *CM_RAcc_SoundSpeed;
extern double    *CM_RAcc_GasDens;
extern double    *CM_RAcc_RelativeVel;                     // the relative velocity between BH and gas
extern double    *CM_RAcc_ColdGasMass;
extern double    *CM_RAcc_GasMass;
extern double    *CM_RAcc_ParMass;
extern int       *CM_Cluster_NPar_close;
extern double    *CM_Bondi_SinkMass, *CM_Bondi_SinkMomX, *CM_Bondi_SinkMomY, *CM_Bondi_SinkMomZ;
extern double    *CM_Bondi_SinkMomXAbs, *CM_Bondi_SinkMomYAbs, *CM_Bondi_SinkMomZAbs;
extern double    *CM_Bondi_SinkE, *CM_Bondi_SinkEk, *CM_Bondi_SinkEt;
extern int       *CM_Bondi_SinkNCell;
// ---------------------------------------------------------------------------------------
// (4) other variables
extern double    *E_inj_exp;                               // the expected amount of injected energy
extern double    *M_inj_exp;                               // the expected amount of injected gas mass
extern double   (*ang_mom_sum)[3];
extern int        AdjustCount;                             // count the number of adjustments
extern int        Merger_Coll_NumBHs;

extern FieldIdx_t Idx_ParHalo;

extern double    *Jet_WaveK;                               // jet wavenumber used in the sin() function to have smooth bidirectional jets
extern double    *V_cyl;                                   // the volume of jet source
extern double    *M_inj, *P_inj, *E_inj;                   // the injected density
extern double    *normalize_const;                         // the exact normalization constant

extern long_par  *CM_ClusterIdx_Cur;

static bool       if_overlap = false;
static int        merge_index = 0;                         // record BH 1 merge BH 2 / BH 2 merge BH 1

static bool       FirstTime = true;
// =======================================================================================


// problem-specific function prototypes
// =======================================================================================
// (1) external functions
#ifdef MHD
extern double (*MHD_ResetByUser_VecPot_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[] );
extern double (*MHD_ResetByUser_BField_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[],
                                             const double B_in, const bool UseVecPot, const real *Ax, const real *Ay,
                                             const real *Az, const int i, const int j, const int k );
#endif
// ---------------------------------------------------------------------------------------
// (2) internal functions
static void GetClusterCenter( int lv, bool AdjustPos, bool AdjustVel, double Cen_old[][3], double Cen_new[][3], double Cen_Vel[][3] );
static void SetJetDirection( const double TimeNew, const int lv, const int FluSg );
// =======================================================================================



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
   if ( ( mode == 1  ||  mode == 3 )  &&  rho_gas > 0.0 )
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

   const double Pos[3] = { x, y, z };

// (1) SMBH Accretion (note: gas depletion is currently disabled)
   double dr2[3][3], r2[3];
   const double V_dep = 4.0 / 3.0 * M_PI * pow( R_dep, 3.0 ); // the volume to remove gas
   double D_dep[Merger_Coll_NumBHs]; // the density need to be removed
   for (int c=0; c<Merger_Coll_NumBHs; c++)   D_dep[c] = CM_BH_Mdot_tot[c]*dt/V_dep;

   int reset = false; // mark whether this cell is reset or not [false/true]

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      for (int d=0; d<3; d++)   dr2[c][d] = SQR( Pos[d] - CM_ClusterCen[c][d] );
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
      if ( Merger_Coll_NumBHs == 1 )   Aux_Error( ERROR_INFO, "Error: Merger_Coll_NumBHs = 1 but if_overlap = true!\n" );
      if ( CM_Jet_Edot[0] >= CM_Jet_Edot[1] )   status = 0;   // only inject cluster 1
      else                                      status = 1;   // only inject cluster 2
      n_jet = Merger_Coll_NumBHs-1;
   }
   else // if ( if_overlap )
   {
      status = 0;
      n_jet  = Merger_Coll_NumBHs;
   } // if ( if_overlap ) ...  else ...

   for (int c=status; c<(n_jet+status); c++)
   {
//    distance: jet center to mesh
      for (int d=0; d<3; d++)   Vec_c2m[d] = Pos[d] - CM_ClusterCen[c][d];

      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );

      Jet_dh = fabs( CM_Jet_Vec[c][0]*Vec_c2m[0] + CM_Jet_Vec[c][1]*Vec_c2m[1] + CM_Jet_Vec[c][2]*Vec_c2m[2] );
      Jet_dr = sqrt( SQR(Dis_c2m) - SQR(Jet_dh) );

      if ( Jet_dh <= Jet_HalfHeight[c]  &&  Jet_dr <= Jet_Radius[c] )
      {
         which_cluster += c+1;

         const bool Check_MinEint_No = false;

//       record the old fluid variables
         const real dens_old = fluid[DENS];
#        ifdef MHD
         const real emag_old = (real)0.5 * ( SQR(fluid[MAG_OFFSET+0]) + SQR(fluid[MAG_OFFSET+1]) + SQR(fluid[MAG_OFFSET+2]) );
#        else
         const real emag_old = (real)0.0;
#        endif
         const real eint_old = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ],
                                               fluid[ENGY], Check_MinEint_No, NULL_REAL, PassiveFloorMask, emag_old,
                                               NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                               h_EoS_Table );

//       accrete mass
         fluid[DENS] += M_inj[c];

//       transfer into BH frame
         fluid[MOMX] -= CM_BH_Vel[c][0]*dens_old;
         fluid[MOMY] -= CM_BH_Vel[c][1]*dens_old;
         fluid[MOMZ] -= CM_BH_Vel[c][2]*dens_old;

//       use a sine function to make the velocity smooth within the jet from +CM_Jet_Vec to -CM_Jet_Vec
         EngySin = E_inj[c]*normalize_const[c]*sin( Jet_WaveK[c]*Jet_dh );

//       the new momentum is calculated from the old density, new density, old momentum and injected energy
//       the momentum injection alters only the component along the jet direction
         real   P_old_sqr = SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]);
         real   Ekin_old  = ( dens_old == 0.0 ) ? (real)0.0 : 0.5*P_old_sqr/dens_old;
         real   P_new     = SQRT( 2*fluid[DENS]*(EngySin + Ekin_old) );
         double JetSign   = SIGN( Vec_c2m[0]*CM_Jet_Vec[c][0] + Vec_c2m[1]*CM_Jet_Vec[c][1] + Vec_c2m[2]*CM_Jet_Vec[c][2] );

         double P_old_perp[3], P_old_para, P_new_para;
         P_old_para    = fluid[MOMX] * CM_Jet_Vec[c][0] + fluid[MOMY] * CM_Jet_Vec[c][1] + fluid[MOMZ] * CM_Jet_Vec[c][2];
         P_old_perp[0] = fluid[MOMX] - P_old_para * CM_Jet_Vec[c][0];
         P_old_perp[1] = fluid[MOMY] - P_old_para * CM_Jet_Vec[c][1];
         P_old_perp[2] = fluid[MOMZ] - P_old_para * CM_Jet_Vec[c][2];
         P_new_para    = SQRT( SQR(P_new) - SQR(P_old_perp[0]) - SQR(P_old_perp[1]) - SQR(P_old_perp[2]) );
         P_new_para   *= JetSign;

         fluid[MOMX] = P_new_para * CM_Jet_Vec[c][0] + P_old_perp[0];
         fluid[MOMY] = P_new_para * CM_Jet_Vec[c][1] + P_old_perp[1];
         fluid[MOMZ] = P_new_para * CM_Jet_Vec[c][2] + P_old_perp[2];

//       transfer back into the rest frame
         fluid[MOMX] += CM_BH_Vel[c][0] * fluid[DENS];
         fluid[MOMY] += CM_BH_Vel[c][1] * fluid[DENS];
         fluid[MOMZ] += CM_BH_Vel[c][2] * fluid[DENS];

         fluid[ENGY]  = 0.5 * (SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ])) / fluid[DENS] +
                              eint_old + emag_old;
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
   const double  dv           = CUBE(dh);
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
   GetClusterCenter( lv, AdjustPosNow, AdjustVelNow, CM_BH_Pos, CM_ClusterCen, CM_BH_Vel );


// (2) decide whether to merge BHs
   const double AbsRelPos = DIST_3D_DBL( CM_BH_Pos[0], CM_BH_Pos[1] );
   const double AbsRelVel = DIST_3D_DBL( CM_BH_Vel[0], CM_BH_Vel[1] );
   const double soften    = amr->dh[MAX_LEVEL];
   double escape_vel;
   if ( AbsRelPos > soften )
   {
      escape_vel = sqrt( 2 * NEWTON_G * (CM_BH_Mass[0]+CM_BH_Mass[1]) / AbsRelPos );
   }
   else
   {
      escape_vel = sqrt( 2 * NEWTON_G * (CM_BH_Mass[0]+CM_BH_Mass[1]) / soften );
   } // if ( AbsRelPos > soften ) ... else ...

// merge the two BHs if they are located within R_acc, and the relative velocity is small enough
   if ( Merger_Coll_NumBHs == 2 )
   {
      if ( AbsRelPos < R_acc  &&  AbsRelVel < 3*escape_vel )
      {
         Merger_Coll_NumBHs -= 1;
         if ( CM_BH_Mass[0] >= CM_BH_Mass[1] )   merge_index = 1;   // record BH 1 merge BH 2 / BH 2 merge BH 1
         else                                    merge_index = 2;
         CM_BH_Mass[0] += CM_BH_Mass[1];
         CM_BH_Mass[1] = 0.0;

//       relabel the BH and DM particles being merged
         CM_ClusterIdx_Cur[1] = 0;
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
                      CM_BH_Pos[0][0], CM_BH_Pos[0][1], CM_BH_Pos[0][2], CM_BH_Pos[1][0], CM_BH_Pos[1][1], CM_BH_Pos[1][2],
                      CM_BH_Vel[0][0], CM_BH_Vel[0][1], CM_BH_Vel[0][2], CM_BH_Vel[1][0], CM_BH_Vel[1][1], CM_BH_Vel[1][2],
                      AbsRelPos, AbsRelVel, escape_vel );
      } // if ( AbsRelPos < R_acc  &&  AbsRelVel < 3*escape_vel )
   } // if ( Merger_Coll_NumBHs == 2 )


   if ( AGN_feedback )
   {
//    (3) set the injection parameters
//    set the jet direction vector
      SetJetDirection( TimeNew, lv, FluSg );


//    (4) calculate the accretion and feedback
//    reset to 0 since we only want to record the number of void cells **for one sub-step**
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         CM_Bondi_SinkNCell[c] = 0;
         Jet_WaveK[c]          = 0.5 * M_PI / Jet_HalfHeight[c];
      }

//    variables for each rank
      int    num        [Merger_Coll_NumBHs];     // the number of cells inside the accretion radius
      double gas_mass   [Merger_Coll_NumBHs];     // total gas mass inside the accretion radius
      double rho        [Merger_Coll_NumBHs];     // the average density inside the accretion radius (hot gas)
      double mass_cold  [Merger_Coll_NumBHs];     // cold gas mass (T < 5e5 K) inside the accretion radius
      double Cs         [Merger_Coll_NumBHs];     // the average sound speed inside the accretion radius
      double gas_mom    [Merger_Coll_NumBHs][3];  // average gas momentum
      double V_cyl_exact[Merger_Coll_NumBHs];     // exact volume of jet cylinder
      double normalize  [Merger_Coll_NumBHs];     // for computing the correct normalization constant
      bool if_overlap_each_rank = false;
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         num        [c] = 0;
         gas_mass   [c] = 0.0;
         rho        [c] = 0.0;
         mass_cold  [c] = 0.0;
         Cs         [c] = 0.0;
         V_cyl_exact[c] = 0.0;
         normalize  [c] = 0.0;
         for (int d=0; d<3; d++)
            gas_mom[c][d] = 0.0;
      }


      real fluid_acc[NCOMP_TOTAL];

//    variables for all ranks
      int    num_sum[Merger_Coll_NumBHs];
      double rho_sum[Merger_Coll_NumBHs], Cs_sum[Merger_Coll_NumBHs], gas_mom_sum[Merger_Coll_NumBHs][3];
      double gas_vel_sum[Merger_Coll_NumBHs][3], V_cyl_exact_sum[Merger_Coll_NumBHs], normalize_sum[Merger_Coll_NumBHs];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         const double x02 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y02 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z02 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)   { const double z2 = z02 + k*dh;
         for (int j=0; j<PS1; j++)   { const double y2 = y02 + j*dh;
         for (int i=0; i<PS1; i++)   { const double x2 = x02 + i*dh;

            double Pos_2[3] = { x2, y2, z2 };

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
               if ( DIST_SQR_3D( Pos_2, CM_ClusterCen[c] ) <= SQR(R_acc) )
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
                                                    fluid_acc[4], fluid_acc+NCOMP_FLUID, true, MIN_PRES, PassiveFloorMask, Emag,
                                                    EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                                    EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                    h_EoS_Table, NULL );
                  const real Temp = Hydro_Con2Temp( fluid_acc[0], fluid_acc[1], fluid_acc[2], fluid_acc[3],
                                                    fluid_acc[4], fluid_acc+NCOMP_FLUID, false, MIN_TEMP, PassiveFloorMask, Emag,
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
               } // if ( DIST_SQR_3D( Pos_2, CM_ClusterCen[c] ) <= SQR(R_acc) )

//             calculate the exact volume of jet cylinder and normalization
               if ( CurrentMaxLv )
               {
                  double Jet_dr_2, Jet_dh_2, Dis_c2m_2, Vec_c2m_2[3];

                  for (int d=0; d<3; d++)   Vec_c2m_2[d] = Pos_2[d] - CM_ClusterCen[c][d];
                  Dis_c2m_2 = sqrt( SQR(Vec_c2m_2[0]) + SQR(Vec_c2m_2[1]) + SQR(Vec_c2m_2[2]) );
                  Jet_dh_2 = fabs( CM_Jet_Vec[c][0]*Vec_c2m_2[0] + CM_Jet_Vec[c][1]*Vec_c2m_2[1] + CM_Jet_Vec[c][2]*Vec_c2m_2[2] );
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
         MPI_Allreduce( &num[c],         &num_sum[c],             1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &gas_mass[c],    &CM_RAcc_GasMass[c],     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &rho[c],         &rho_sum[c],             1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &mass_cold[c],   &CM_RAcc_ColdGasMass[c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &Cs[c],          &Cs_sum[c],              1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( gas_mom[c],       gas_mom_sum[c],         3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &V_cyl_exact[c], &V_cyl_exact_sum[c],     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &normalize[c],   &normalize_sum[c],       1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
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

               if ( DIST_SQR_3D( patch_pos, CM_ClusterCen[c] ) > SQR(2*R_acc+patch_d) )   continue;

               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  const long     ParID     = amr->patch[0][lv][PID]->ParList[p];
                  const real_par ParM      = amr->Par->Mass[ParID];
                  const real_par ParPos[3] = { amr->Par->PosX[ParID], amr->Par->PosY[ParID], amr->Par->PosZ[ParID] };

                  if ( DIST_SQR_3D( ParPos, CM_ClusterCen[c] ) > SQR(R_acc) )   continue;

                  par_mass += ParM;
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            MPI_Allreduce( &par_mass, &CM_RAcc_ParMass[c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         } // if ( Accretion_Mode == 2  ||  Accretion_Mode == 3 )
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//    calculate the accretion rate
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         double v = 0.0;  // the relative velocity between BH and gas
         if ( num_sum[c] == 0 )
         {
            CM_RAcc_GasDens[c]     = 0.0;
            CM_RAcc_SoundSpeed[c]  = 0.0;
            CM_RAcc_RelativeVel[c] = 0.0;
            for (int d=0; d<3; d++)  CM_RAcc_GasVel[c][d] = 0.0;
         }
         else
         {
            for (int d=0; d<3; d++)  gas_vel_sum[c][d] = gas_mom_sum[c][d] / rho_sum[c];
            rho_sum[c] /= ( 4.0 / 3.0 * M_PI * pow(R_acc, 3) );
            Cs_sum[c]  /= (double)num_sum[c];
            for (int d=0; d<3; d++)   v += SQR( CM_BH_Vel[c][d] - gas_vel_sum[c][d] );

            CM_RAcc_GasDens[c]     = rho_sum[c];
            CM_RAcc_SoundSpeed[c]  = Cs_sum[c];
            CM_RAcc_RelativeVel[c] = sqrt(v);
            for (int d=0; d<3; d++)  CM_RAcc_GasVel[c][d] = gas_vel_sum[c][d];
         } // if ( num_sum[c] == 0 ) ... else ...

         BH_accretion_rate( Accretion_Mode, CM_BH_Mdot_tot+c, CM_BH_Mdot_hot+c, CM_BH_Mdot_cold+c,
                            R_acc, CM_BH_Mass[c], CM_RAcc_GasDens[c], CM_RAcc_SoundSpeed[c], CM_RAcc_RelativeVel[c],
                            CM_RAcc_ColdGasMass[c], CM_RAcc_GasMass[c], CM_RAcc_ParMass[c] );

         if ( V_cyl_exact_sum[c] != 0 )   normalize_const[c] = V_cyl_exact_sum[c] / normalize_sum[c];
         else                             normalize_const[c] = 0.5 * M_PI;

      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//    update BH mass
      for (int c=0; c<Merger_Coll_NumBHs; c++)   if ( CurrentMaxLv )   CM_BH_Mass[c] += CM_BH_Mdot_tot[c] * dt;

//    (5) calculate the injection rate
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         CM_Jet_Mdot[c] = eta * CM_BH_Mdot_tot[c];
         CM_Jet_Pdot[c] = sqrt(2*eta*eps_f*(1.0-eps_m)) * CM_BH_Mdot_tot[c] * (Const_c/UNIT_V);
         CM_Jet_Edot[c] = eps_f * CM_BH_Mdot_tot[c] * SQR(Const_c/UNIT_V);
         V_cyl      [c] = M_PI * SQR(Jet_Radius[c]) * 2 * Jet_HalfHeight[c];

//       calculate the density that need to be injected
         if ( CurrentMaxLv  &&  V_cyl_exact_sum[c] != 0.0 )
         {
            M_inj[c] = CM_Jet_Mdot[c] * dt / V_cyl_exact_sum[c];
            P_inj[c] = CM_Jet_Pdot[c] * dt / V_cyl_exact_sum[c];
            E_inj[c] = CM_Jet_Edot[c] * dt / V_cyl_exact_sum[c];
         }
         else
         {
            M_inj[c] = CM_Jet_Mdot[c] * dt / V_cyl[c];
            P_inj[c] = CM_Jet_Pdot[c] * dt / V_cyl[c];
            E_inj[c] = CM_Jet_Edot[c] * dt / V_cyl[c];
         } // if ( CurrentMaxLv  &&  V_cyl_exact_sum[c] != 0.0 ) ... else ...
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

      if ( CurrentMaxLv )
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)   E_inj_exp[c] += CM_Jet_Edot[c]*dt;
         for (int c=0; c<Merger_Coll_NumBHs; c++)   M_inj_exp[c] += CM_Jet_Mdot[c]*dt;
      } // if ( CurrentMaxLv )


//    (6) perform injection
//    use the "static" schedule for reproducibility
#     pragma omp parallel for schedule( static )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         int Reset = 0;
         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)  {  const double z = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  const double y = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  const double x = x0 + i*dh;
            real fluid[NCOMP_TOTAL], fluid_old[NCOMP_TOTAL];

            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

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
                                                       (real)MIN_EINT, PassiveFloorMask, Emag );

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

#                 pragma omp critical
                  {
                  CM_Bondi_SinkMass   [Reset-1] += dv *     ( fluid[DENS] - fluid_old[DENS] );
                  CM_Bondi_SinkMomX   [Reset-1] += dv *     ( fluid[MOMX] - fluid_old[MOMX] );
                  CM_Bondi_SinkMomY   [Reset-1] += dv *     ( fluid[MOMY] - fluid_old[MOMY] );
                  CM_Bondi_SinkMomZ   [Reset-1] += dv *     ( fluid[MOMZ] - fluid_old[MOMZ] );
                  CM_Bondi_SinkMomXAbs[Reset-1] += dv * FABS( fluid[MOMX] - fluid_old[MOMX] );
                  CM_Bondi_SinkMomYAbs[Reset-1] += dv * FABS( fluid[MOMY] - fluid_old[MOMY] );
                  CM_Bondi_SinkMomZAbs[Reset-1] += dv * FABS( fluid[MOMZ] - fluid_old[MOMZ] );
                  CM_Bondi_SinkE      [Reset-1] += dv *     ( fluid[ENGY] - fluid_old[ENGY] );
                  CM_Bondi_SinkEk     [Reset-1] += dv *     ( Ekin_new    - Ekin_old        );
                  CM_Bondi_SinkEt     [Reset-1] += dv *     ( Eint_new    - Eint_old        );
                  CM_Bondi_SinkNCell  [Reset-1] ++;
                  }
               } // if ( CurrentMaxLv )

//             store the reset values
               for (int v=0; v<NCOMP_TOTAL; v++)
                  amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
            } // if ( Reset != 0 )
         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // if ( AGN_feedback )

} // FUNCTION : Flu_ResetByUser_API_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  GetClusterCenter
// Description :  Get the cluster centers
//
// Note        :  1. Must enable Merger_Coll_LabelCenter
//
// Parameter   :  lv        : Refeinement level to search on
//                AdjustPos : If true, adjust the positions of the BHs
//                AdjustVel : If true, adjust the velocities of the BHs
//                Cen_old   : Old BH position array with the shape of [Merger_Coll_NumBHs][3]
//                Cen_new   : New BH position array with the shape of [Merger_Coll_NumBHs][3]
//                Cen_Vel   : BH CoM velocity array with the shape of [Merger_Coll_NumBHs][3]
//
// Return      :  Cen_new, Cen_Vel (passed into and updated by this function)
//-------------------------------------------------------------------------------------------------------
void GetClusterCenter( int lv, bool AdjustPos, bool AdjustVel, double Cen_old[][3], double Cen_new[][3], double Cen_Vel[][3] )
{

// fix the BH position and rest BH
   if ( fixBH )
   {
      for (int d=0; d<3; d++)   Cen_new[0][d] = amr->BoxCenter[d];
      for (int d=0; d<3; d++)   Cen_Vel[0][d] = 0.0;

      for (int c=0; c<Merger_Coll_NumBHs; c++)
         for (int d=0; d<3; d++)  Cen_old[c][d] = Cen_new[c][d];

      return;
   } // if ( fixBH )

   double pos_min[Merger_Coll_NumBHs][3], DM_Vel[Merger_Coll_NumBHs][3];   // the updated BH position and velocity
   const bool CurrentMaxLv = ( NPatchTotal[lv] > 0  &&  lv == MAX_LEVEL        ) ? true :
                             ( NPatchTotal[lv] > 0  &&  NPatchTotal[lv+1] == 0 ) ? true : false;

// initialize pos_min to be the old center
   for (int c=0; c<Merger_Coll_NumBHs; c++)   for (int d=0; d<3; d++)   pos_min[c][d] = Cen_old[c][d];

   if ( CurrentMaxLv  &&  (AdjustPos  ||  AdjustVel) )
   {

      const double dis_exp = 1e-6;  // to check if the output BH positions of each calculaiton are close enough
      bool   converged     = false; // if the BH positions are close enough, then complete the calculation
      int    counter       = 0;     // how many times the calculation is performed (minimum: 2, maximum: 10)
      double Cen_new_pre[Merger_Coll_NumBHs][3];

      while ( converged == false  &&  counter <= 10 )
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)
            for (int d=0; d<3; d++)  Cen_new_pre[c][d] = pos_min[c][d];

         int N_max[Merger_Coll_NumBHs]; // maximum particle numbers (to allocate the array size)
         for (int c=0; c<Merger_Coll_NumBHs; c++)   N_max[c] = 10000;

         int      num_par[Merger_Coll_NumBHs];   // (each rank) number of particles inside the target region of each cluster
         real_par **ParX = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         real_par **ParY = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         real_par **ParZ = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         real_par **ParM = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         real_par **VelX = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         real_par **VelY = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         real_par **VelZ = (real_par**)malloc( Merger_Coll_NumBHs*sizeof(real_par*) );
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            ParX[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
            ParY[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
            ParZ[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
            ParM[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
            VelX[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
            VelY[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
            VelZ[c] = (real_par*)malloc( N_max[c]*sizeof(real_par) );
         }

//       find the particles within 10 times the accretion radius
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            num_par[c] = 0;
            CM_Cluster_NPar_close[c] = 0;
            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               const double *EdgeL        = amr->patch[0][lv][PID]->EdgeL;
               const double *EdgeR        = amr->patch[0][lv][PID]->EdgeR;
               const double  patch_pos[3] = { (EdgeL[0]+EdgeR[0])*0.5, (EdgeL[1]+EdgeR[1])*0.5, (EdgeL[2]+EdgeR[2])*0.5 };
               const double  patch_d      = DIST_3D_DBL( EdgeL, EdgeR ) * 0.5;

               if ( DIST_SQR_3D( patch_pos, Cen_new_pre[c] ) > SQR(10*R_acc+patch_d) )   continue;

               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  const long_par ParID         = amr->patch[0][lv][PID]->ParList[p];
                  const real_par ParX_tmp      = amr->Par->PosX[ParID];
                  const real_par ParY_tmp      = amr->Par->PosY[ParID];
                  const real_par ParZ_tmp      = amr->Par->PosZ[ParID];
                  const real_par ParM_tmp      = amr->Par->Mass[ParID];
                  const real_par VelX_tmp      = amr->Par->VelX[ParID];
                  const real_par VelY_tmp      = amr->Par->VelY[ParID];
                  const real_par VelZ_tmp      = amr->Par->VelZ[ParID];
                  const real_par ParPos_tmp[3] = { ParX_tmp, ParY_tmp, ParZ_tmp };

                  if ( CM_ClusterIdx_Cur[amr->Par->AttributeInt[Idx_ParHalo][ParID]] != (long_par)c )   continue;
                  if ( DIST_SQR_3D( ParPos_tmp, Cen_new_pre[c] ) > SQR(10*R_acc) )   continue;

//                record the mass, position and velocity of this particle
                  ParX[c][num_par[c]] = ParX_tmp;
                  ParY[c][num_par[c]] = ParY_tmp;
                  ParZ[c][num_par[c]] = ParZ_tmp;
                  ParM[c][num_par[c]] = ParM_tmp;
                  VelX[c][num_par[c]] = VelX_tmp;
                  VelY[c][num_par[c]] = VelY_tmp;
                  VelZ[c][num_par[c]] = VelZ_tmp;
                  num_par[c] += 1;

                  if ( num_par[c] >= N_max[c] )
                  {
//                   increase the new maximum size if needed
                     N_max[c] = (int)ceil( PARLIST_GROWTH_FACTOR*(num_par[c]+1) );
                     ParX[c]  = (real_par*)realloc( ParX[c], N_max[c]*sizeof(real_par) );
                     ParY[c]  = (real_par*)realloc( ParY[c], N_max[c]*sizeof(real_par) );
                     ParZ[c]  = (real_par*)realloc( ParZ[c], N_max[c]*sizeof(real_par) );
                     ParM[c]  = (real_par*)realloc( ParM[c], N_max[c]*sizeof(real_par) );
                     VelX[c]  = (real_par*)realloc( VelX[c], N_max[c]*sizeof(real_par) );
                     VelY[c]  = (real_par*)realloc( VelY[c], N_max[c]*sizeof(real_par) );
                     VelZ[c]  = (real_par*)realloc( VelZ[c], N_max[c]*sizeof(real_par) );
                  }
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//       collect the number of target particles from each rank
         MPI_Allreduce( num_par, CM_Cluster_NPar_close, Merger_Coll_NumBHs, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

         int num_par_eachRank[Merger_Coll_NumBHs][MPI_NRank];
         int displs[Merger_Coll_NumBHs][MPI_NRank];
         int CM_Cluster_NPar_close_max = 0;
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            MPI_Allgather( &num_par[c], 1, MPI_INT, num_par_eachRank[c], 1, MPI_INT, MPI_COMM_WORLD );
            displs[c][0] = 0;
            for (int i=1; i<MPI_NRank; i++)  displs[c][i] = displs[c][i-1] + num_par_eachRank[c][i-1];
            CM_Cluster_NPar_close_max = MAX( CM_Cluster_NPar_close_max, CM_Cluster_NPar_close[c] );
         }

//       collect the mass, position and velocity of target particles to the root rank
         real_par **ParX_sum = new real_par* [Merger_Coll_NumBHs];
         real_par **ParY_sum = new real_par* [Merger_Coll_NumBHs];
         real_par **ParZ_sum = new real_par* [Merger_Coll_NumBHs];
         real_par **ParM_sum = new real_par* [Merger_Coll_NumBHs];
         real_par **VelX_sum = new real_par* [Merger_Coll_NumBHs];
         real_par **VelY_sum = new real_par* [Merger_Coll_NumBHs];
         real_par **VelZ_sum = new real_par* [Merger_Coll_NumBHs];
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            ParX_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
            ParY_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
            ParZ_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
            ParM_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
            VelX_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
            VelY_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
            VelZ_sum[c] = new real_par [CM_Cluster_NPar_close[c]];
         }

         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            MPI_Allgatherv( ParX[c], num_par[c], MPI_GAMER_REAL_PAR, ParX_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
            MPI_Allgatherv( ParY[c], num_par[c], MPI_GAMER_REAL_PAR, ParY_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
            MPI_Allgatherv( ParZ[c], num_par[c], MPI_GAMER_REAL_PAR, ParZ_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
            MPI_Allgatherv( ParM[c], num_par[c], MPI_GAMER_REAL_PAR, ParM_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
            MPI_Allgatherv( VelX[c], num_par[c], MPI_GAMER_REAL_PAR, VelX_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
            MPI_Allgatherv( VelY[c], num_par[c], MPI_GAMER_REAL_PAR, VelY_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
            MPI_Allgatherv( VelZ[c], num_par[c], MPI_GAMER_REAL_PAR, VelZ_sum[c], num_par_eachRank[c],
                            displs[c], MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
         }

//       compute potential and find the minimum position, and calculate the average DM velocity on the root rank
         if ( AdjustPos )
         {
            double  soften        = amr->dh[MAX_LEVEL];
            double *pote_AllRank  = new double [CM_Cluster_NPar_close_max];
            double *pote_ThisRank = new double [CM_Cluster_NPar_close_max];
            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
//             distribute MPI jobs
               int par_per_rank = CM_Cluster_NPar_close[c] / MPI_NRank;
               int remainder    = CM_Cluster_NPar_close[c] % MPI_NRank;
               int start        = MPI_Rank*par_per_rank + MIN( MPI_Rank, remainder );
               int end          = start + par_per_rank + (MPI_Rank < remainder ? 1 : 0);

#              pragma omp parallel for schedule( static )
               for (int i=start; i<end; i++)
               {
                  pote_ThisRank[i-start] = 0.0;
                  for (int j=0; j<CM_Cluster_NPar_close[c]; j++)
                  {
                     if ( i == j )   continue;

                     const double rel_pos = sqrt( SQR(ParX_sum[c][i]-ParX_sum[c][j]) + SQR(ParY_sum[c][i]-ParY_sum[c][j]) +
                                                  SQR(ParZ_sum[c][i]-ParZ_sum[c][j]) );
                     if       ( rel_pos >  soften )   pote_ThisRank[i-start] += ParM_sum[c][j] / rel_pos;
                     else if  ( rel_pos <= soften )   pote_ThisRank[i-start] += ParM_sum[c][j] / soften;
                  }
                  pote_ThisRank[i-start] *= -NEWTON_G;
               }

               int N_recv[MPI_NRank], N_disp[MPI_NRank];
               for (int i=0; i<MPI_NRank; i++)  N_recv[i] = (i < remainder ? par_per_rank+1 : par_per_rank);
               N_disp[0] = 0;
               for (int i=1; i<MPI_NRank; i++)  N_disp[i] = N_disp[i-1] + N_recv[i-1];
               MPI_Allgatherv( pote_ThisRank, end-start, MPI_DOUBLE, pote_AllRank, N_recv, N_disp, MPI_DOUBLE, MPI_COMM_WORLD );

               double Pote_min = 0.0;
               for (int i=0; i<CM_Cluster_NPar_close[c]; i++)
               {
                  if ( pote_AllRank[i] >= Pote_min )  continue;
                  Pote_min      = pote_AllRank[i];
                  pos_min[c][0] = ParX_sum[c][i];
                  pos_min[c][1] = ParY_sum[c][i];
                  pos_min[c][2] = ParZ_sum[c][i];
               }
            } // for (int c=0; c<Merger_Coll_NumBHs; c++)
            delete[] pote_AllRank;
            delete[] pote_ThisRank;
         } // if ( AdjustPos )

//       calculate the average DM velocity
         if ( AdjustVel )
         {
            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
               for (int d=0; d<3; d++)   DM_Vel[c][d] = 0.0;
               double ParM_Tot = 0.0;
               for (int i=0; i<CM_Cluster_NPar_close[c]; i++)
               {
                  DM_Vel[c][0] += VelX_sum[c][i]*ParM_sum[c][i];
                  DM_Vel[c][1] += VelY_sum[c][i]*ParM_sum[c][i];
                  DM_Vel[c][2] += VelZ_sum[c][i]*ParM_sum[c][i];
                  ParM_Tot     += ParM_sum[c][i];
               }
               for (int d=0; d<3; d++)   DM_Vel[c][d] /= ParM_Tot;
            }
         } // if ( AdjustVel )

//       iterate the above calculation until the output BH positions become close enough
         counter += 1;
         double dis[Merger_Coll_NumBHs];
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            dis[c] = 0.0;
            for (int d=0; d<3; d++)   dis[c] += SQR( pos_min[c][d] - Cen_new_pre[c][d] );
         }

         if ( counter > 1 )
         {
            bool all_BH_converged = true;
            for (int c=0; c<Merger_Coll_NumBHs; c++)
               if ( sqrt(dis[c]) >= dis_exp )   all_BH_converged = false;

            converged |= all_BH_converged;
         }

         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            delete [] ParX_sum[c];
            delete [] ParY_sum[c];
            delete [] ParZ_sum[c];
            delete [] ParM_sum[c];
            delete [] VelX_sum[c];
            delete [] VelY_sum[c];
            delete [] VelZ_sum[c];
         }
         delete [] ParX_sum;
         delete [] ParY_sum;
         delete [] ParZ_sum;
         delete [] ParM_sum;
         delete [] VelX_sum;
         delete [] VelY_sum;
         delete [] VelZ_sum;

         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            free( ParX[c] );
            free( ParY[c] );
            free( ParZ[c] );
            free( ParM[c] );
            free( VelX[c] );
            free( VelY[c] );
            free( VelZ[c] );
         }
         free( ParX );
         free( ParY );
         free( ParZ );
         free( ParM );
         free( VelX );
         free( VelY );
         free( VelZ );
      } // while ( converged == false )
   } // if ( CurrentMaxLv  &&  (AdjustPos  ||  AdjustVel) )

// find the BH particles and adjust their position and velocity
   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      double Cen_Tmp[3] = { -__FLT_MAX__, -__FLT_MAX__, -__FLT_MAX__ }; // set to -inf
      double Vel_Tmp[3] = { -__FLT_MAX__, -__FLT_MAX__, -__FLT_MAX__ };
      for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
      {
         if ( amr->Par->Mass[p] < (real_par)0.0 )   continue;
         if ( CM_ClusterIdx_Cur[amr->Par->AttributeInt[Idx_ParHalo][p]] != (long_par)c )   continue;
         if ( amr->Par->Type[p] != PTYPE_BLACK_HOLE )   continue;

         if ( CurrentMaxLv  &&  AdjustPos )
         {
            amr->Par->PosX[p] = pos_min[c][0];
            amr->Par->PosY[p] = pos_min[c][1];
            amr->Par->PosZ[p] = pos_min[c][2];
         }
         if ( CurrentMaxLv  &&  AdjustVel )
         {
            amr->Par->VelX[p] = DM_Vel[c][0];
            amr->Par->VelY[p] = DM_Vel[c][1];
            amr->Par->VelZ[p] = DM_Vel[c][2];
         }
         Cen_Tmp[0] = amr->Par->PosX[p];
         Cen_Tmp[1] = amr->Par->PosY[p];
         Cen_Tmp[2] = amr->Par->PosZ[p];
         Vel_Tmp[0] = amr->Par->VelX[p];
         Vel_Tmp[1] = amr->Par->VelY[p];
         Vel_Tmp[2] = amr->Par->VelZ[p];
         break;
      } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)

//    use MPI_MAX since Cen_Tmp[] is initialized as -inf
      MPI_Allreduce( Cen_Tmp, Cen_new[c], 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( Vel_Tmp, Cen_Vel[c], 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)

   if ( CurrentMaxLv  &&  AdjustPos )
   {
      const bool TimingSendPar_No = false;
      Par_PassParticle2Sibling( lv, TimingSendPar_No );
      Par_PassParticle2Son_MultiPatch( lv, PAR_PASS2SON_EVOLVE, TimingSendPar_No, NULL_INT, NULL );
   }

   for (int c=0; c<Merger_Coll_NumBHs; c++)
      for (int d=0; d<3; d++)  Cen_old[c][d] = Cen_new[c][d];

} // FUNCTION : GetClusterCenter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetJetDirection
// Description :  Set the jet direction
//
// Parameter   :  TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//
// Return      :  CM_Jet_Vec (global variable)
//-------------------------------------------------------------------------------------------------------
void SetJetDirection( const double TimeNew, const int lv, const int FluSg )
{

   double ang_mom[Merger_Coll_NumBHs][3]; // total angular momentum inside the accretion radius

   switch ( JetDirection_case )
   {
      case 1: // fixed at x-axis
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            CM_Jet_Vec[c][0] = 1.0;
            CM_Jet_Vec[c][1] = 0.0;
            CM_Jet_Vec[c][2] = 0.0;
         }
         break;
      case 2: // import from table
         const double Time_period      = CM_Jet_Time_table[JetDirection_NBin-1];
         const double Time_interpolate = fmod( TimeNew, Time_period );
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            const double theta = Mis_InterpolateFromTable( JetDirection_NBin, CM_Jet_Time_table, CM_Jet_Theta_table[c], Time_interpolate );
            const double phi   = Mis_InterpolateFromTable( JetDirection_NBin, CM_Jet_Time_table, CM_Jet_Phi_table[c],   Time_interpolate );

            CM_Jet_Vec[c][0] = cos(theta);
            CM_Jet_Vec[c][1] = sin(theta)*cos(phi);
            CM_Jet_Vec[c][2] = sin(theta)*sin(phi);
         }
         break;
      case 3: // align with angular momentum
         const double dh = amr->dh[lv];
         const double dv = CUBE(dh);

         for (int c=0; c<Merger_Coll_NumBHs; c++)
            for (int d=0; d<3; d++)
               ang_mom[c][d] = 0.0;

         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
               const double pos2[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (0.5+i)*dh,
                                        amr->patch[0][lv][PID]->EdgeL[1] + (0.5+j)*dh,
                                        amr->patch[0][lv][PID]->EdgeL[2] + (0.5+k)*dh };

               for (int c=0; c<Merger_Coll_NumBHs; c++)
               {
                  if ( DIST_SQR_3D( pos2, CM_ClusterCen[c] ) > SQR(R_acc) ) continue;

                  double dr[3] = { pos2[0]-CM_ClusterCen[c][0], pos2[1]-CM_ClusterCen[c][1], pos2[2]-CM_ClusterCen[c][2] };
                  ang_mom[c][0] += dv * ( dr[1]*amr->patch[FluSg][lv][PID]->fluid[3][k][j][i] - dr[2]*amr->patch[FluSg][lv][PID]->fluid[2][k][j][i] );
                  ang_mom[c][1] += dv * ( dr[2]*amr->patch[FluSg][lv][PID]->fluid[1][k][j][i] - dr[0]*amr->patch[FluSg][lv][PID]->fluid[3][k][j][i] );
                  ang_mom[c][2] += dv * ( dr[0]*amr->patch[FluSg][lv][PID]->fluid[2][k][j][i] - dr[1]*amr->patch[FluSg][lv][PID]->fluid[1][k][j][i] );
               } // for (int c=0; c<Merger_Coll_NumBHs; c++)
            } // for (int k=0; k<PS1; k++); for (int j=0; j<PS1; j++); for (int i=0; i<PS1; i++)
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            MPI_Allreduce( ang_mom[c], ang_mom_sum[c], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            const double ang_mom_norm = sqrt( SQR(ang_mom_sum[c][0]) + SQR(ang_mom_sum[c][1]) + SQR(ang_mom_sum[c][2]) );

            for (int d=0; d<3; d++)   CM_Jet_Vec[c][d] = ang_mom_sum[c][d] / ang_mom_norm;
         } // for (int c=0; c<Merger_Coll_NumBHs; c++)
         break;
      default:
         Aux_Error( ERROR_INFO, "Unsupported JetDirection_case %d [1/2/3] !!\n", JetDirection_case );
         break;
   } // switch ( JetDirection_case )

} // FUNCTION : SetJetDirection



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )
