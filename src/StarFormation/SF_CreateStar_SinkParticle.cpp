#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )

//-------------------------------------------------------------------------------------------------------
// Function    :  ProximityCheck_SecondDenThres
// Description :  Check whether the gas cell is close to any existing particle and satisfies the second density threshold
//
// Note        :  1. Ref: Clarke et al., 2017, MNRAS, 468, 2489, eqn (5)
//
// Parameter   :  ExistingNPar   : Number of existing particles
//                ParAtt_Local   : Attribute of existing particles
//                PosX/Y/Z       : Position of the gas cell
//                VelX/Y/Z       : Velocity of the gas cell
//                GasDens        : Density of the gas cell
//                AccRadius      : Accretion radius
//                Coeff_FreeFall : Coefficient for free-fall time
//
// Return      :  true : The gas cell passes the check
//                false: The gas cell fails the check
//-------------------------------------------------------------------------------------------------------
bool ProximityCheck_SecondDenThres( const int ExistingNPar, real *ParAtt_Local[],
                                    const real PosX, const real PosY, const real PosZ,
                                    const real VelX, const real VelY, const real VelZ,
                                    const real GasDens, const real AccRadius, const real Coeff_FreeFall )
{
   real Par2Cell[3], Par2CellDist, Par2CellVel[3];
   real NorPar2Cell[3];
   real GasDensFreeFall;

   for (int p=0; p<ExistingNPar; p++) // loop over all nearby existing particles
   {
      Par2Cell[0] = PosX - ParAtt_Local[PAR_POSX][p];
      Par2Cell[1] = PosY - ParAtt_Local[PAR_POSY][p];
      Par2Cell[2] = PosZ - ParAtt_Local[PAR_POSZ][p];
      Par2CellDist = SQRT(SQR(Par2Cell[0])+SQR(Par2Cell[1])+SQR(Par2Cell[2]));

      if ( Par2CellDist <= 2*AccRadius )  return false;

      Par2CellVel[0] = VelX - ParAtt_Local[PAR_VELX][p];
      Par2CellVel[1] = VelY - ParAtt_Local[PAR_VELY][p];
      Par2CellVel[2] = VelZ - ParAtt_Local[PAR_VELZ][p];

      NorPar2Cell[0] = Par2Cell[0]/Par2CellDist;
      NorPar2Cell[1] = Par2Cell[1]/Par2CellDist;
      NorPar2Cell[2] = Par2Cell[2]/Par2CellDist;

      if ( Par2CellVel[0]/NorPar2Cell[0] >= 0 )                       continue; 
      if ( Par2CellVel[1]/NorPar2Cell[1] >= 0 )                       continue;
      if ( Par2CellVel[2]/NorPar2Cell[2] >= 0 )                       continue; // the gas is moving away from the existing particle

      GasDensFreeFall = SQR((1/Coeff_FreeFall)*(NorPar2Cell[0]*Par2CellVel[0] + NorPar2Cell[1]*Par2CellVel[1] + NorPar2Cell[2]*Par2CellVel[2])/Par2CellDist); // Clarke et al. 2017, eqn (5)
      if ( GasDens < GasDensFreeFall )    return false;
   } // for (int p=0; p<ExistingNPar; p++) 

   return true;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  GravMini
// Description :  Check whether the gas cell has the minimum gravitational potential inside the control volume
//
// Parameter   :  pi, pj, pk   : Indices of the current cell
//                AccCellNum   : Accretion radius in cells
//                Pot_Array_F  : Potential array
//                Size_Flu     : Size of the fluid array
//                Phi000       : Potential of the current cell
//
// Return      :  true : The gas cell is the local minimum
//                false: The gas cell is NOT the local minimum
//-------------------------------------------------------------------------------------------------------
bool GravMini( const int pi, const int pj, const int pk, const real AccCellNum, 
               const real *Pot_Array_F, const int Size_Flu, const real Phi000 )
{
   for (int vk=pk-AccCellNum; vk<=pk+AccCellNum; vk++)
   for (int vj=pj-AccCellNum; vj<=pj+AccCellNum; vj++)
   for (int vi=pi-AccCellNum; vi<=pi+AccCellNum; vi++)
   {
      if ( SQRT(SQR(vi - pi)+SQR(vj - pj)+SQR(vk - pk)) > AccCellNum )           continue;

      const int vt = IDX321( vi, vj, vk, Size_Flu, Size_Flu );
      const real Phiijk = Pot_Array_F[vt];

      if ( Phiijk < Phi000 )  return false;
   }

   return true;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  ConvergingFlow
// Description :  Check whether the gas in the nearby six cells converge toward the central one
//
// Parameter   :  pi, pj, pk      : Indices of the current cell
//                Size_Flu        : Size of the fluid array
//                Flu_Array_F_In  : Fluid array (4D: [v][k][j][i])
//
// Return      :  true : The flow is converging
//                false: The flow is NOT converging
//-------------------------------------------------------------------------------------------------------
bool ConvergingFlow( const int pi, const int pj, const int pk, const int Size_Flu, 
                     const real (*Flu_Array_F_In)[Size_Flu][Size_Flu][Size_Flu] )
{
   real VelNeighbor[6];
   
   VelNeighbor[0] = Flu_Array_F_In[MOMX][pk][pj][pi+1] / Flu_Array_F_In[DENS][pk][pj][pi+1]; // x+ (NeighborID 0)
   VelNeighbor[1] = Flu_Array_F_In[MOMX][pk][pj][pi-1] / Flu_Array_F_In[DENS][pk][pj][pi-1]; // x- (NeighborID 1)
   VelNeighbor[2] = Flu_Array_F_In[MOMY][pk][pj+1][pi] / Flu_Array_F_In[DENS][pk][pj+1][pi]; // y+ (NeighborID 2)
   VelNeighbor[3] = Flu_Array_F_In[MOMY][pk][pj-1][pi] / Flu_Array_F_In[DENS][pk][pj-1][pi]; // y- (NeighborID 3)
   VelNeighbor[4] = Flu_Array_F_In[MOMZ][pk+1][pj][pi] / Flu_Array_F_In[DENS][pk+1][pj][pi]; // z+ (NeighborID 4)
   VelNeighbor[5] = Flu_Array_F_In[MOMZ][pk-1][pj][pi] / Flu_Array_F_In[DENS][pk-1][pj][pi]; // z- (NeighborID 5)

   if ( (VelNeighbor[0] - VelNeighbor[1]) >= 0 || 
        (VelNeighbor[2] - VelNeighbor[3]) >= 0 || 
        (VelNeighbor[4] - VelNeighbor[5]) >= 0 )                      
      return false;

   return true;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  JeansInstability_BoundState
// Description :  Check for Jeans instability and bound state
//
// Parameter   :  pi, pj, pk      : Indices of the current cell
//                AccCellNum      : Accretion radius in cells
//                Flu_Array_F_In  : Fluid array (4D: [v][k][j][i])
//                Size_Flu        : Size of the fluid array
//                Mag_Array_F_In  : Magnetic field array (linearized)
//                Mag_Stride      : Stride for magnetic field array
//                dv              : Cell volume
//                dh              : Cell size
//
// Return      :  true : The control volume is unstable and bound (create particle)
//                false: Otherwise
//-------------------------------------------------------------------------------------------------------
bool JeansInstability_BoundState( const int pi, const int pj, const int pk, const real AccCellNum,
                                  const int Size_Flu, const real (*Flu_Array_F_In)[Size_Flu][Size_Flu][Size_Flu],
                                  const real *Mag_Array_F_In, const int Mag_Stride,
                                  const real dv, const real dh )
{
   real ControlFluid[FLU_NIN], ControlFluidj[FLU_NIN];
   real Pres, Cs2, vEmag=NULL_REAL;

   // calculate bulk velocity
   real TotalMass = (real)0.0, MassVel[3] = { (real)0.0, (real)0.0, (real)0.0}, BulkVel[3];
   
   for (int vk=pk-AccCellNum; vk<=pk+AccCellNum; vk++)
   for (int vj=pj-AccCellNum; vj<=pj+AccCellNum; vj++)
   for (int vi=pi-AccCellNum; vi<=pi+AccCellNum; vi++)
   {
      if ( SQRT(SQR(vi - pi)+SQR(vj - pj)+SQR(vk - pk)) > AccCellNum )           continue;

      for (int v=0; v<FLU_NIN; v++)    ControlFluid[v] = Flu_Array_F_In[v][vk][vj][vi];

      MassVel[0] += ControlFluid[MOMX]*dv;
      MassVel[1] += ControlFluid[MOMY]*dv;
      MassVel[2] += ControlFluid[MOMZ]*dv;
      TotalMass  += ControlFluid[DENS]*dv;
   }

   BulkVel[0] = MassVel[0]/TotalMass;
   BulkVel[1] = MassVel[1]/TotalMass;
   BulkVel[2] = MassVel[2]/TotalMass;

   // get the energy
   real Egtot = (real)0.0, Ethtot = (real)0.0, Emagtot = (real)0.0, Ekintot = (real)0.0;
   const bool CheckMinPres_No = false;

   for (int vk=pk-AccCellNum; vk<=pk+AccCellNum; vk++)
   for (int vj=pj-AccCellNum; vj<=pj+AccCellNum; vj++)
   for (int vi=pi-AccCellNum; vi<=pi+AccCellNum; vi++)
   {
      if ( SQRT(SQR(vi - pi)+SQR(vj - pj)+SQR(vk - pk)) > AccCellNum )           continue;

      for (int v=0; v<FLU_NIN; v++)    ControlFluid[v] = Flu_Array_F_In[v][vk][vj][vi];

      // Storing Egtot
      real SelfPhiijk = (real)0.0;
      for (int vkj=pk-AccCellNum; vkj<=pk+AccCellNum; vkj++)
      for (int vjj=pj-AccCellNum; vjj<=pj+AccCellNum; vjj++)
      for (int vij=pi-AccCellNum; vij<=pi+AccCellNum; vij++)
      {
         if ( SQRT(SQR(vij - pi)+SQR(vjj - pj)+SQR(vkj - pk)) > AccCellNum )           continue;

         int rijPix = SQRT(SQR(vi - vij)+SQR(vj - vjj)+SQR(vk - vkj));
         if ( rijPix == 0 )                        continue;

         real rij = rijPix*dh;
         
         for (int v=0; v<FLU_NIN; v++)    ControlFluidj[v] = Flu_Array_F_In[v][vkj][vjj][vij];

         SelfPhiijk += -NEWTON_G*ControlFluidj[DENS]*dv/rij;
      }

      Egtot += 0.5*ControlFluid[DENS]*dv*SelfPhiijk;

      // Storing Emagtot
#     ifdef MHD
      vEmag = MHD_GetCellCenteredBEnergy( Mag_Array_F_In + MAGX*Mag_Stride,
                                          Mag_Array_F_In + MAGY*Mag_Stride,
                                          Mag_Array_F_In + MAGZ*Mag_Stride,
                                          Size_Flu, Size_Flu, Size_Flu, vi, vj, vk );
      Emagtot += vEmag*dv;
#     endif

      // Storing Ethtot and Ekintot
      Pres = Hydro_Con2Pres( ControlFluid[DENS], ControlFluid[MOMX], ControlFluid[MOMY], ControlFluid[MOMZ], ControlFluid[ENGY],
                             ControlFluid+NCOMP_FLUID, CheckMinPres_No, NULL_REAL, PassiveFloorMask, vEmag,
                             EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                             EoS_AuxArray_Flt,
                             EoS_AuxArray_Int, h_EoS_Table, NULL );
      Cs2  = EoS_DensPres2CSqr_CPUPtr( ControlFluid[DENS], Pres, ControlFluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      
      Ethtot += 0.5*ControlFluid[DENS]*dv*Cs2;

      Ekintot += 0.5*ControlFluid[DENS]*dv*( SQR(ControlFluid[MOMX]/ControlFluid[DENS] - BulkVel[0]) + SQR(ControlFluid[MOMY]/ControlFluid[DENS] - BulkVel[1]) + SQR(ControlFluid[MOMZ]/ControlFluid[DENS] - BulkVel[2]));
   }

   if ( FABS(Egtot) <= 2*Ethtot )                      return false;
   if (( Egtot + Ethtot + Ekintot + Emagtot ) >= 0)    return false;

   return true;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_SinkParticle
// Description :  Create sink particles based on FALSH prescription
//
// Note        :  1. Ref: (1) C. Federrath et al., 2010, ApJ, 713, 269
//                        (2) S. D. Clarke et al., 2017, MNRAS, 468, 2489
//                2. One must turn on STORE_POT_GHOST when adopting STORE_PAR_ACC
//                   --> It is because, currently, this function always uses the pot_ext[] array of each patch
//                       to calculate the gravitationally acceleration of the new star particles
//                3. One must invoke Buf_GetBufferData( ..., _TOTAL, ... ) "before" calling this function
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Current physical time (after advancing solution by dt)
//                GasDensThres : Minimum gas density for creating sink particles                (--> "SF_CREATE_SINK_MIN_GAS_DENS"  )
//                AccCellNum   : Accretion radius in cells at highest refinement level          (--> "SF_CREATE_SINK_ACC_RADIUS"    )
//                MaxNewPar    : Maximum number of particle per MPI rank                        (--> "SF_CREATE_SINK_MAX_NPAR_MPI"  )
//
// Return      :  1. Particle repository will be updated
//                2. fluid[] array of gas will be updated
//-------------------------------------------------------------------------------------------------------

void SF_CreateStar_SinkParticle( const int lv, const real TimeNew, const real GasDensThres, const real AccCellNum, 
                                 const int MaxNewPar )
{
// check
#  if ( defined STORE_PAR_ACC  &&  !defined STORE_POT_GHOST )
#     error : STAR_FORMATION + STORE_PAR_ACC must work with STORE_POT_GHOST !!
#  endif

#  ifndef GRAVITY
#     error : must turn on GRAVITY for SF_CreateStar_SinkParticle() !!
#  endif

#  ifdef COMOVING
#     error : SF_CreateStar_SinkParticle() does not support COMOVING yet !!
#  endif

if ( lv != MAX_LEVEL )
   Aux_Error( ERROR_INFO, "Sink particles should be created at MAX_LEVEL, please set SF_CREATE_STAR_MIN_LEVEL to -1.\n" );

#  ifdef GAMER_DEBUG
   if ( Idx_ParCreTime == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParCreTime is undefined !!\n" );

   if ( Idx_ParID == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParID is undefined !!\n" );
#  endif // #ifdef GAMER_DEBUG

// constant parameters
   const double dh             = amr->dh[lv];
   const double dv             = CUBE( dh );

   const int    NGhost         = PS1 / 2; // the number of ghost cell at each side
   const int    Size_Flu       = PS2 + 2*NGhost; // final cube size
   const int    Size_Flu_P1    = Size_Flu + 1; // for face-centered B field
   const int    Size_Pot       = Size_Flu; // for potential
   const int    NPG            = 1;
   
   const real   AccRadius      = AccCellNum*dh;
   const int    FluSg          = amr->FluSg[lv];
   const real   Coeff_FreeFall = SQRT( (32.0*NEWTON_G)/(3.0*M_PI) );
// const real   GraConst       = ( OPT__GRA_P5_GRADIENT ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh);
   const real   GraConst       = ( false                ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // P5 is NOT supported yet

   int          NNewPar = 0; // number of new particle
   real       (*RemovalFlu)[5]                      = new real     [MaxNewPar][5]; // information used to remove gas from the cell
   long       (*RemovalPos)[4]                      = new long     [MaxNewPar][4]; // some IDs used to remove gas from the cell
   real_par   (*NewParAttFlt)[PAR_NATT_FLT_TOTAL]   = new real_par [MaxNewPar][PAR_NATT_FLT_TOTAL]; // attribute of the new particles
   long_par   (*NewParAttInt)[PAR_NATT_INT_TOTAL]   = new long_par [MaxNewPar][PAR_NATT_INT_TOTAL];
   long        *NewParRepo                          = new long     [MaxNewPar]; // particle respository that will be added
   long        *NewParPID                           = new long     [MaxNewPar]; // PID where the particle is formed

// Some check
   if (AccCellNum > 0.5*PS1)
      Aux_Error( ERROR_INFO, "AccCellNum should be smaller than half of PATCH_SIZE !!" );

   if (OPT__MINIMIZE_MPI_BARRIER)
      Aux_Error( ERROR_INFO, "Please disable OPT__MINIMIZE_MPI_BARRIER !!" );

   const bool TimingSendPar_Yes = true;
   const bool JustCountNPar_No  = false;
   const bool PredictPos_No     = false;
   const bool SibBufPatch_Yes   = true;
   const bool FaSibBufPatch_No  = false;
   const long ParAttFltBitIdx_In = _PAR_FLT_TOTAL;
   const long ParAttIntBitIdx_In = _PAR_INT_TOTAL;
   Par_CollectParticle2OneLevel( lv, ParAttFltBitIdx_In, ParAttIntBitIdx_In, PredictPos_No, TimeNew, SibBufPatch_Yes, FaSibBufPatch_No,
                                 JustCountNPar_No, TimingSendPar_Yes );

// start of OpenMP parallel region
#  pragma omp parallel
   {
   real PosX, PosY, PosZ, VelX, VelY, VelZ; // position and velocity of the current test cell
   real GasDens; // gas density of the current test cell
   real GasMFracLeft; // the fraction of gas mass of the cell that will be took out by particle
   real Corner_Array_F[3]; // the corner of the ghost zone
   real fluid[FLU_NIN]; // fluid in the current test cell

   real PotNeighbor[6]; // record the neighboring cell potential [x+, x-, y+, y+, z+, z-]

   real   (*Flu_Array_F_In)[Size_Flu][Size_Flu][Size_Flu] = new real [FLU_NIN][Size_Flu][Size_Flu][Size_Flu];
   real   (*Mag_Array_F_In)[Size_Flu_P1*SQR(Size_Flu)]     = new real [NCOMP_MAG][Size_Flu_P1*SQR(Size_Flu)];
   real   (*Pot_Array_F)                                   = new real [CUBE(Size_Pot)];

   int LocalPID, delta_t, PGi, PGj, PGk;

// get the sibling index differences along different directions
   int NSibPID_Delta[26], *SibPID_Delta[26];

   TABLE_GetSibPID_Delta( NSibPID_Delta, SibPID_Delta );

// loop over all real patch groups
// use static schedule to ensure bitwise reproducibility when running with the same numbers of OpenMP threads and MPI ranks
// --> bitwise reproducibility will still break when running with different numbers of OpenMP threads and/or MPI ranks
//     unless both BITWISE_REPRODUCIBILITY and SF_CREATE_STAR_DET_RANDOM are enabled
#  pragma omp for schedule( static )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
#     ifndef MHD
      const int    OPT__MAG_INT_SCHEME = INT_NONE;
#     endif
      const bool   IntPhase_No         = false;
      const real   MinDens_No          = -1.0;
      const real   MinPres_No          = -1.0;
      const real   MinTemp_No          = -1.0;
      const real   MinEntr_No          = -1.0;
      const bool   DE_Consistency_Yes  = true;
      const bool   DE_Consistency_No   = false;
      const bool   DE_Consistency      = ( OPT__OPTIMIZE_AGGRESSIVE ) ? DE_Consistency_No : DE_Consistency_Yes;
      const real   MinDens             = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;

#     ifdef MHD
      real *Mag_Array = Mag_Array_F_In[0];
#     else
      real *Mag_Array = NULL;
#     endif
      Prepare_PatchData( lv, TimeNew, (real*)Flu_Array_F_In, Mag_Array,
                        NGhost, NPG, &PID0, _TOTAL, _MAG,
                        OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                        OPT__BC_FLU, BC_POT_NONE, MinDens,    MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency );

//    prepare the potential array
      if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
      Prepare_PatchData( lv, TimeNew, Pot_Array_F, NULL,
                        NGhost, NPG, &PID0, _POTE, _NONE,
                        OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                        OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

//    prepare the corner array
      for (int d=0; d<3; d++)    Corner_Array_F[d] = amr->patch[0][lv][PID0]->EdgeL[d] + 0.5*dh - dh*NGhost;

//    prepare the sibling patches id for this patch group
//    collect patch information
      const int NNearbyPatchMax = 64;  // maximum number of neaby patches of a patch group (including 8 local patches)
      int Nearby_PID_List[NNearbyPatchMax], NNearbyPatch, SibPID0_List[26];

//    get nearby patches
      NNearbyPatch = 0;

//    local patches
      for (int PID=PID0; PID<PID0+8; PID++)  Nearby_PID_List[ NNearbyPatch ++ ] = PID;

//    sibling patches
      TABLE_GetSibPID_Based( lv, PID0, SibPID0_List );

//###OPTIMIZATION: skip sibling patches if the maximum feedback radius is zero
      for (int s=0; s<26; s++)
      {
         const int SibPID0 = SibPID0_List[s];   // first target patch in the sibling patch group

//       only consider leaf patches on FB_LEVEL (including both real and buffer patches)
         if ( SibPID0 >= 0 )
         for (int c=0; c<NSibPID_Delta[s]; c++)
         {
            const int SibPID = SibPID0 + SibPID_Delta[s][c];
            Nearby_PID_List[ NNearbyPatch ++ ] = SibPID;
         }
      }

//    prepare local arrays for particle attributes in this patch group
      real  *ParAtt_Local[PAR_NATT_FLT_TOTAL];
      int    NParAllPatch   = 0; // the number of existing particles in all nearby patches

      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   ParAtt_Local[v] = NULL; // array to store nearby particle

      for (int t=0; t<NNearbyPatch; t++)
      {
         const int PPID = Nearby_PID_List[t]; // check the particle number for this patch (PPID)
         int   NParMaxInPatch   = 0;

      // check both NPar and NPar_Copy (NPar_Copy may be -1, which is fine)
         NParMaxInPatch = MAX( NParMaxInPatch, amr->patch[0][lv][PPID]->NPar      );
         NParMaxInPatch = MAX( NParMaxInPatch, amr->patch[0][lv][PPID]->NPar_Copy );

         if (NParMaxInPatch > 0 )   NParAllPatch += NParMaxInPatch;
      }

      int    ExistingNPar = 0;
      if ( NParAllPatch > 0 )
      {
         for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
            if ( ParAttFltBitIdx_In & BIDX(v) )    ParAtt_Local[v] = new real [NParAllPatch];

         for (int t=0; t<NNearbyPatch; t++)
         {
            const int PPID = Nearby_PID_List[t];
            long  *ParList = NULL;
            int    NPar;
            bool   UseParAttCopy;

//          the check "son == -1" is actually useless for now since we fix LEVEL == MAX_LEVEL
            if ( amr->patch[0][lv][PPID]->son == -1  &&  PPID < amr->NPatchComma[lv][1] )
            {
               NPar          = amr->patch[0][lv][PPID]->NPar;
               ParList       = amr->patch[0][lv][PPID]->ParList;
               UseParAttCopy = false;

#              ifdef DEBUG_PARTICLE
               if ( amr->patch[0][lv][PPID]->NPar_Copy != -1 )
                  Aux_Error( ERROR_INFO, "lv %d, PPID %d, NPar_Copy = %d != -1 !!\n",
                           lv, PPID, amr->patch[0][lv][PPID]->NPar_Copy );
#              endif
            }

            else
            {
//             note that amr->patch[0][lv][PPID]->NPar>0 is still possible
               NPar          = amr->patch[0][lv][PPID]->NPar_Copy;
#              ifdef LOAD_BALANCE
               ParList       = NULL;
               UseParAttCopy = true;
#              else
               ParList       = amr->patch[0][lv][PPID]->ParList_Copy;
               UseParAttCopy = false;
#              endif
            } // if ( amr->patch[0][lv][PPID]->son == -1  &&  PPID < amr->NPatchComma[lv][1] ) ... else ...

#           ifdef LOAD_BALANCE
            if ( UseParAttCopy ) {
               for (int v=0; v<PAR_NATT_FLT_TOTAL; v++) {
                  if ( ParAttFltBitIdx_In & BIDX(v) ) {

#                 ifdef DEBUG_PARTICLE
                  if ( NPar > 0  &&  amr->patch[0][lv][PPID]->ParAtt_Copy[v] == NULL )
                     Aux_Error( ERROR_INFO, "ParAtt_Copy == NULL for NPar (%d) > 0 (lv %d, PPID %d, v %d) !!\n",
                                 NPar, lv, PPID, v );
#                 endif

                  for (int p=0; p<NPar; p++)
                     ParAtt_Local[v][ExistingNPar + p] = amr->patch[0][lv][PPID]->ParAttFlt_Copy[v][p];
            }}}

            else
#           endif // #ifdef LOAD_BALANCE
            {
#              ifdef DEBUG_PARTICLE
               if ( NPar > 0  &&  ParList == NULL )
                  Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, PPID %d) !!\n",
                              NPar, lv, PPID );
#              endif

               for (int v=0; v<PAR_NATT_FLT_TOTAL; v++) {
                  if ( ParAttFltBitIdx_In & BIDX(v) )
                     for (int p=0; p<NPar; p++)
                        ParAtt_Local[v][ExistingNPar + p] = amr->Par->AttributeFlt[v][ ParList[p] ];
               }
            } // if ( UseParAttCopy ) ... else ...

            ExistingNPar += NPar;
         } // for (int t=0; t<NNearbyPatch; t++)
      } // if ( NParAllPatch > 0 )

//    loop over all cells in the patch group
      for (int pk=NGhost; pk<PS2 + NGhost; pk++)
      for (int pj=NGhost; pj<PS2 + NGhost; pj++)
      for (int pi=NGhost; pi<PS2 + NGhost; pi++)
      {  
         PGi = pi - NGhost;
         PGj = pj - NGhost;
         PGk = pk - NGhost; // the cell id inside patch group

         LocalPID = 2*2*(PGk/PS1) + 2*(PGj/PS1) + (PGi/PS1);

         const int Disp_i = TABLE_02( LocalPID, 'x', 0, PS1 ); // the cell index within PID
         const int Disp_j = TABLE_02( LocalPID, 'y', 0, PS1 );
         const int Disp_k = TABLE_02( LocalPID, 'z', 0, PS1 );
         const int PID = PID0 + LocalPID; // record the current PID

//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         PosX = Corner_Array_F[0] + pi*dh;
         PosY = Corner_Array_F[1] + pj*dh;
         PosZ = Corner_Array_F[2] + pk*dh;

         const int t = IDX321( pi, pj, pk, Size_Flu, Size_Flu );
         for (int v=0; v<FLU_NIN; v++)    fluid[v] = Flu_Array_F_In[v][pk][pj][pi];
         VelX = fluid[MOMX]/fluid[DENS];
         VelY = fluid[MOMY]/fluid[DENS];
         VelZ = fluid[MOMZ]/fluid[DENS];

//       First density threshold:
//       Check whether the gas density is larger than the threshold
//       ===========================================================================================================
         GasDens = fluid[DENS];
         if ( GasDens <= GasDensThres )    continue;

//       Proximity check + second density threshold:
//       The gas cell should not near any existing particle
//       We can use the distance between the gas cell and nearby particle and the relative velocity to define a 
//       a density threshold (Clarke et al. 2017, eqn (5))
//       ===========================================================================================================
         if ( !ProximityCheck_SecondDenThres( ExistingNPar, ParAtt_Local, PosX, PosY, PosZ, VelX, VelY, VelZ, 
                                              GasDens, AccRadius, Coeff_FreeFall ) )
            continue;

//       Gravitational minimum check inside the control volume:
//       The gas cell should have the minimum gravitational potential inside the control volume
//       ===========================================================================================================
         const real Phi000 = Pot_Array_F[t]; // the potential of the current cell
         if ( !GravMini( pi, pj, pk, AccCellNum, Pot_Array_F, Size_Flu, Phi000 ) )
            continue;
         
//       Converging flow check:
//       The gas in the nearby six cells should converge toward the central one
//       ===========================================================================================================
         if ( !ConvergingFlow( pi, pj, pk, Size_Flu, Flu_Array_F_In ) )
            continue;

//       Jeans instability check + check for bound state
//       The control volume is Jeans unstable if | Egtot | >= 2*Ethtot
//       The control volume is bound if ( Ethtot + Ekintot + Emagtot ) >= | Egtot |
//       ===========================================================================================================
         if ( !JeansInstability_BoundState( pi, pj, pk, AccCellNum, Size_Flu, Flu_Array_F_In, 
                                            (real*)Mag_Array_F_In, Size_Flu_P1*SQR(Size_Flu), dv, dh ) )
            continue;

//       Store the information of new star particles
//       ===========================================================================================================
         if ( NNewPar >= MaxNewPar )
            Aux_Error( ERROR_INFO, "NNewPar (%d) >= MaxNewPar (%d) !! Please try a larger MaxNewPar.\n", NNewPar, MaxNewPar );

#        pragma omp critical
         {
            NewParAttFlt[NNewPar][PAR_MASS]  = (GasDens - GasDensThres)*dv;
            NewParAttFlt[NNewPar][PAR_POSX]  = PosX;
            NewParAttFlt[NNewPar][PAR_POSY]  = PosY;
            NewParAttFlt[NNewPar][PAR_POSZ]  = PosZ;
            NewParAttFlt[NNewPar][PAR_VELX]  = VelX;
            NewParAttFlt[NNewPar][PAR_VELY]  = VelY;
            NewParAttFlt[NNewPar][PAR_VELZ]  = VelZ;
            NewParAttFlt[NNewPar][PAR_TIME]  = TimeNew;
            NewParAttInt[NNewPar][PAR_TYPE]  = PTYPE_STAR;
            NewParAttInt[NNewPar][Idx_ParID] = -1; // initialize the value to be -1

   //       particle acceleration
   #        ifdef STORE_PAR_ACC
            real GasAcc[3] = { (real)0.0, (real)0.0, (real)0.0 };

   //       external acceleration
            if ( OPT__EXT_ACC )  CPUExtAcc_Ptr( GasAcc, PosX, PosY, PosZ, TimeNew, ExtAcc_AuxArray );

   //       self-gravity and external potential
            if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
            {
               for (int NeighborID=0; NeighborID<6; NeighborID++)
               {  
                  if      (NeighborID == 0) delta_t = IDX321(  1,  0,  0, Size_Flu, Size_Flu );
                  else if (NeighborID == 1) delta_t = IDX321( -1,  0,  0, Size_Flu, Size_Flu );
                  else if (NeighborID == 2) delta_t = IDX321(  0,  1,  0, Size_Flu, Size_Flu );
                  else if (NeighborID == 3) delta_t = IDX321(  0, -1,  0, Size_Flu, Size_Flu );
                  else if (NeighborID == 4) delta_t = IDX321(  0,  0,  1, Size_Flu, Size_Flu );
                  else if (NeighborID == 5) delta_t = IDX321(  0,  0, -1, Size_Flu, Size_Flu );

                  const int Neighbort = t + delta_t;
                  PotNeighbor[NeighborID] = Pot_Array_F[Neighbort];
               }
               GasAcc[0] += GraConst*(PotNeighbor[0] - PotNeighbor[1]);
               GasAcc[1] += GraConst*(PotNeighbor[2] - PotNeighbor[3]);
               GasAcc[2] += GraConst*(PotNeighbor[4] - PotNeighbor[5]);
            }

            NewParAttFlt[NNewPar][PAR_ACCX] = GasAcc[0];
            NewParAttFlt[NNewPar][PAR_ACCY] = GasAcc[1];
            NewParAttFlt[NNewPar][PAR_ACCZ] = GasAcc[2];
   #        endif // ifdef STORE_PAR_ACC

            NewParAttFlt[NNewPar][Idx_ParCreTime  ] = TimeNew;
            NewParPID[NNewPar] = PID;

            GasMFracLeft = GasDensThres/GasDens;
            RemovalPos[NNewPar][0] = PID;
            RemovalPos[NNewPar][1] = PGk - Disp_k;
            RemovalPos[NNewPar][2] = PGj - Disp_j;
            RemovalPos[NNewPar][3] = PGi - Disp_i;

            RemovalFlu[NNewPar][0] = GasMFracLeft;
            RemovalFlu[NNewPar][1] = Phi000;
            RemovalFlu[NNewPar][2] = PosX;
            RemovalFlu[NNewPar][3] = PosY;
            RemovalFlu[NNewPar][4] = PosZ;

            NNewPar ++;
         } // # pragma omp critical
      } // pi, pj, pk

      if ( NParAllPatch > 0 )
      {
         for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)      
         {
            delete [] ParAtt_Local[v]; ParAtt_Local[v] = NULL;
         }
      }

   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8) #  pragma omp for schedule( static )

   delete [] Flu_Array_F_In;
   delete [] Mag_Array_F_In;
   delete [] Pot_Array_F;
   } // end of OpenMP parallel region

// Excluding the nearby particles + remove the gas from the cell
// ===========================================================================================================
   MPI_Barrier(MPI_COMM_WORLD);

   int      *GatherNNewPar         = new int [MPI_NRank]; // the number of candidates for each rank
   MPI_Allgather(&NNewPar, 1, MPI_INT, GatherNNewPar, 1, MPI_INT, MPI_COMM_WORLD);

   int       TotalNNewPar          = 0; // get the total number of the candidates
   int      *RecvRemovalFluSize    = new int [MPI_NRank];
   for (int rank=0; rank<MPI_NRank; rank++) 
   {
      TotalNNewPar += GatherNNewPar[rank];
      RecvRemovalFluSize[rank] = 5*GatherNNewPar[rank]; // receive "count" for each MPI rank
   }

   int      *disp                  = new int [MPI_NRank];
   disp[0] = 0;
   for (int rank=1; rank<MPI_NRank; rank++) disp[rank] = disp[rank-1] + RecvRemovalFluSize[rank-1];

   int      SendRemovalFluSize     = NNewPar*5; // send "count" for the current MPI rank
   real    (*GatherRemovalFlu)[5]  = new real [TotalNNewPar][5]; // the array containing all the candidates + their information

   MPI_Allgatherv(RemovalFlu[0], SendRemovalFluSize, MPI_GAMER_REAL, 
                  GatherRemovalFlu[0], RecvRemovalFluSize, disp, MPI_GAMER_REAL, MPI_COMM_WORLD);

   long     *SelNewParPID          = new long [TotalNNewPar]; // PID of the selected paritcles
   long     *SelNewParID           = new long [TotalNNewPar]; // ID (in the current rank) of the selected paritcles
   real DeltaXPar2Par, DeltaYPar2Par, DeltaZPar2Par, DistPar2Par;   // calculate the distance between the two candidated particles
   int SelNNewPar = 0; // the number of selected particles after the following check
   for (int pi=0; pi<NNewPar; pi++)
   {  
      bool CreateHere = true;
      for (int pj=0; pj<TotalNNewPar; pj++)
      {
         DeltaXPar2Par = RemovalFlu[pi][2] - GatherRemovalFlu[pj][2];
         DeltaYPar2Par = RemovalFlu[pi][3] - GatherRemovalFlu[pj][3];
         DeltaZPar2Par = RemovalFlu[pi][4] - GatherRemovalFlu[pj][4];
         DistPar2Par   = SQRT(SQR(DeltaXPar2Par)+SQR(DeltaYPar2Par)+SQR(DeltaZPar2Par));

         if ( DistPar2Par > 2*AccRadius )                       continue;

         // assuming the potential minimum check is fine, the two particles meet the above conditions should have the same potential
         // if (RemovalFlu[pi][1] != GatherRemovalFlu[pj][1])  continue;   // check whether there are other cells with the same potential
         if ((DeltaXPar2Par<0) or (DeltaYPar2Par<0) or (DeltaZPar2Par<0))
         {
            CreateHere = false;
            break;
         }
      } // for (int pj=0; pj<TotalNNewPar; pj++)

      if ( CreateHere )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)
         amr->patch[FluSg][lv][RemovalPos[pi][0]]->fluid[v][RemovalPos[pi][1]][RemovalPos[pi][2]][RemovalPos[pi][3]] *= RemovalFlu[pi][0];

      // add particles to the particle repository
         SelNewParID[SelNNewPar]  = pi;
         SelNewParPID[SelNNewPar] = NewParPID[pi];
         SelNNewPar++;
      }
   } // for (int pi=0; pi<NNewPar; pi++)

   delete [] GatherNNewPar;
   delete [] RecvRemovalFluSize;
   delete [] disp;
   delete [] GatherRemovalFlu;

// Give the selected particles ID
// ===========================================================================================================
   MPI_Barrier(MPI_COMM_WORLD); 

   long    NParAllRank         = amr->Par->NPar_Active_AllRank;; // Total number of active particles summed up over all MPI ranks
   int    *GatherNSelPar       = new int [MPI_NRank]; // the number of selected particles that pass the above check for each rank
   MPI_Allgather(&SelNNewPar, 1, MPI_INT, GatherNSelPar, 1, MPI_INT, MPI_COMM_WORLD);

   int NParPreRank = 0; // the totol number of selected particle before MPI_Rank (this rank)
   for (int r=0; r<MPI_Rank; r++)     NParPreRank += GatherNSelPar[r];

   for (int p=0; p<SelNNewPar; p++)
   {
      int pi = SelNewParID[p];
      NewParAttInt[pi][Idx_ParID] = NParAllRank + NParPreRank + p; // reassign the ID

      // add particles to the particle repository
      NewParRepo[p] = amr->Par->AddOneParticle( NewParAttFlt[pi], NewParAttInt[pi] );
   }

   delete [] SelNewParID;
   delete [] GatherNSelPar;

// Add the selected particles
// ===========================================================================================================
   long   *UniqueParPID  = new long [MaxNewPar]; // Record the non-repeating PID
   int UniqueCount = 0;
   for (int i=0; i<SelNNewPar; i++)
   {
      int j;
      for (j=0; j<i; j++)
      {
         if (SelNewParPID[i] == SelNewParPID[j])
               break;
      }
      if (i==j)
      {
         UniqueParPID[UniqueCount] = SelNewParPID[i];
         UniqueCount ++;
      }
   } // for (int i=0; i<SelNNewPar; i++)

   const long_par *PType = amr->Par->Type;
   int NParInPatch; // particle number in this patch

   for (int i=0; i<UniqueCount; i++)
   {
      const int SPID = UniqueParPID[i];
      long    *ParRepoInPatch      = new long [MaxNewPar]; // ParRepo in the current patch
      NParInPatch = 0;
      for (int p=0; p<SelNNewPar; p++)
      {
         if (SelNewParPID[p] == SPID) 
         {
            ParRepoInPatch[NParInPatch] = NewParRepo[p];
            NParInPatch ++;
         } // if (SelNewParPID[p] == SPID) 
      } // for (int p=0; p<SelNNewPar; p++)

      if ( NParInPatch == 0 )          
      {
         delete [] ParRepoInPatch;
         continue;
      }              

#     ifdef DEBUG_PARTICLE
//    do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//    may change after calling amr->Par->AddOneParticle()
      const real_par *NewParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
      char Comment[100];
      sprintf( Comment, "%s", __FUNCTION__ );
      
      amr->patch[0][lv][SPID]->AddParticle( NParInPatch, ParRepoInPatch, &amr->Par->NPar_Lv[lv],
                                                         PType, NewParPos, amr->Par->NPar_AcPlusInac, Comment );

#    else
      amr->patch[0][lv][SPID]->AddParticle( NParInPatch, ParRepoInPatch, &amr->Par->NPar_Lv[lv], PType );
#    endif

      delete [] ParRepoInPatch;
   } // for (int i=0; i<UniqueCount; i++)
   
   delete [] SelNewParPID;
   delete [] UniqueParPID;

   delete [] RemovalPos;
   delete [] RemovalFlu;
   delete [] NewParAttFlt;
   delete [] NewParAttInt;
   delete [] NewParRepo;
   delete [] NewParPID;

// free memory
   Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_Yes, FaSibBufPatch_No );

// get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : SF_CreateStar_SinkParticle

#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )