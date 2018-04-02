#include "GAMER.h"

static int Table_01( const int SibID, const int Side, const char dim, const int w01, const int w02,
                     const int w10, const int w11, const int w12, const int w20, const int w21 );




//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolateGhostZone
// Description :  Fill up the ghost-zone values by spatial and temporal interpolation
//
// Note        :  1. Work for the function "Prepare_PatchData"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target variables
//                3. Use the input parameter "IntScheme" to control the interpolation scheme
//                4. Invoke the function "Interpolate" for spatial interpolation
//                5. Data preparation order: FLU -> PASSIVE -> DERIVED --> POTE --> GRA_RHO
//                   ** DERIVED must be prepared immediately after FLU and PASSIVE so that both FLU, PASSIVE, and DERIVED
//                      can be prepared at the same time for the non-periodic BC. **
//                6. Use PrepTime to determine the physical time to prepare data
//                   --> Temporal interpolation/extrapolation will be conducted automatically if PrepTime
//                       is NOT equal to the time of data stored previously (e.g., FluSgTime[0/1])
//                7. For simplicity, currently the mode _TEMP returns **pressure/density**, which does NOT include normalization
//                   --> For OPT__FLAG_LOHNER_TEMP only
//                   --> Also note that MinPres is applied to _TEMP when calculating pressure
//
// Parameter   :  lv             : Target "coarse-grid" refinement level
//                PID            : Patch ID at level "lv" used for interpolation
//                IntData        : Array to store the interpolation result
//                SibID          : Sibling index (0~25) used to determine the interpolation region
//                PrepTime       : Target physical time to prepare data
//                GhostSize      : Number of ghost zones
//                IntScheme      : Interpolation scheme
//                                 --> currently supported schemes include
//                                     INT_MINMOD1D : MinMod-1D
//                                     INT_MINMOD3D : MinMod-3D
//                                     INT_VANLEER  : vanLeer
//                                     INT_CQUAD    : conservative quadratic
//                                     INT_QUAD     : quadratic
//                                     INT_CQUAR    : conservative quartic
//                                     INT_QUAR     : quartic
//                NTSib          : Number of target sibling patches along different sibling directions
//                TSib           : Target sibling indices along different sibling directions
//                TVar           : Target variables to be prepared
//                                 --> Supported variables in different models:
//                                     HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELX, _VELY, _VELZ, _PRES, _TEMP,
//                                             [, _POTE]
//                                     MHD   :
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                                 --> _FLUID, _PASSIVE, _TOTAL, and _DERIVED apply to all models
//                NVar_Tot       : Total number of variables to be prepared
//                NVar_Flu       : Number of fluid variables to be prepared
//                                 --> Including passive scalars
//                TFluVarIdxList : List recording the target fluid and passive variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                NVar_Der       : Number of derived variables to be prepared
//                TDerVarList    : List recording the target derived variables
//                IntPhase       : true --> Perform interpolation on rho/phase instead of real/imag parts in ELBDM
//                FluBC          : Fluid boundary condition
//                PotBC          : Gravity boundary condition (not used currently)
//                BC_Face        : Priority of the B.C. along different boundary faces (z>y>x)
//                MinPres        : Minimum allowed pressure
//                DE_Consistency : Ensure the consistency between pressure, total energy density, and the dual-energy variable
//                                 when DUAL_ENERGY is on
//-------------------------------------------------------------------------------------------------------
void InterpolateGhostZone( const int lv, const int PID, real IntData[], const int SibID, const double PrepTime,
                           const int GhostSize, const IntScheme_t IntScheme, const int NTSib[], int *TSib[],
                           const int TVar, const int NVar_Tot, const int NVar_Flu, const int TFluVarIdxList[],
                           const int NVar_Der, const int TDerVarList[], const bool IntPhase,
                           const OptFluBC_t FluBC[], const OptPotBC_t PotBC, const int BC_Face[], const real MinPres,
                           const bool DE_Consistency )
{

// check
#  ifdef GAMER_DEBUG
// nothing to do if GhostSize == 0
   if ( GhostSize == 0 )
   {
      Aux_Message( stderr, "WARNING : GhostSize == 0 !!\n" );
      return;
   }

// temporal interpolation should be unnecessary for the shared time-step integration
   if ( OPT__DT_LEVEL == DT_LEVEL_SHARED  &&  OPT__INT_TIME )
      Aux_Error( ERROR_INFO, "OPT__INT_TIME should be disabled when \"OPT__DT_LEVEL == DT_LEVEL_SHARED\" !!\n" );
#  endif // #ifdef GAMER_DEBUG


#  if   ( MODEL == HYDRO )
   const bool CheckMinPres_No = false;    // we check minimum pressure in the end of Prepare_PatchData()
   const real Gamma_m1        = GAMMA - (real)1.0;
   const bool PrepVx          = ( TVar & _VELX ) ? true : false;
   const bool PrepVy          = ( TVar & _VELY ) ? true : false;
   const bool PrepVz          = ( TVar & _VELZ ) ? true : false;
   const bool PrepPres        = ( TVar & _PRES ) ? true : false;
   const bool PrepTemp        = ( TVar & _TEMP ) ? true : false;

#  elif ( MODEL == MHD   )
#  warning : WAIT MHD !!

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef GRAVITY
   const bool PrepPot      = ( TVar & _POTE     ) ? true : false;
#  endif

#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   real Fluid[NCOMP_FLUID];   // for calculating pressure and temperature only --> don't need NCOMP_TOTAL
#  endif


// set up parameters for the adopted interpolation scheme
   int NSide, CGhost, CSize[3], FSize[3], CSize3D, FSize3D;

   Int_Table( IntScheme, NSide, CGhost );

   const double dh               = amr->dh[lv];
   const int    GhostSize_Padded = GhostSize + (GhostSize&1);
   const int    CGrid            = (GhostSize+1)/2 + 2*CGhost;    // number of coarse grids required for interpolation
   const int    La               = CGrid - CGhost;                // number of coarse grids in the central region

   for (int d=0; d<3; d++)
   {
      CSize[d] = TABLE_01( SibID, 'x'+d, CGrid, PATCH_SIZE+2*CGhost, CGrid );
      FSize[d] = TABLE_01( SibID, 'x'+d, GhostSize_Padded, 2*PATCH_SIZE, GhostSize_Padded );
   }

   CSize3D = CSize[0]*CSize[1]*CSize[2];
   FSize3D = FSize[0]*FSize[1]*FSize[2];


// we assume that we only need ONE coarse-grid patch in each sibling direction
   if ( La > PATCH_SIZE )  Aux_Error( ERROR_INFO, "La (%d) > PATCH_SIZE (%d) !!\n", La, PATCH_SIZE );


// coarse-grid array to store all the data required for interpolation (including the ghost zones in each side)
   real *CData_Ptr = NULL;
   real *CData     = new real [ NVar_Tot*CSize3D ];


// temporal interpolation parameters
   bool FluIntTime;
   int  FluSg, FluSg_IntT;
   real FluWeighting, FluWeighting_IntT;

   if ( NVar_Flu + NVar_Der != 0 ) {
   if      (  Mis_CompareRealValue( PrepTime, amr->FluSgTime[lv][   amr->FluSg[lv] ], NULL, false )  )
   {
      FluIntTime        = false;
      FluSg             = amr->FluSg[lv];
      FluSg_IntT        = NULL_INT;
      FluWeighting      = NULL_REAL;
      FluWeighting_IntT = NULL_REAL;
   }

   else if (  Mis_CompareRealValue( PrepTime, amr->FluSgTime[lv][ 1-amr->FluSg[lv] ], NULL, false )  )
   {
      FluIntTime        = false;
      FluSg             = 1 - amr->FluSg[lv];
      FluSg_IntT        = NULL_INT;
      FluWeighting      = NULL_REAL;
      FluWeighting_IntT = NULL_REAL;
   }

   else
   {
//    check
      if ( OPT__DT_LEVEL == DT_LEVEL_SHARED )
      Aux_Error( ERROR_INFO, "cannot determine FluSg for OPT__DT_LEVEL == DT_LEVEL_SHARED (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                 lv, PrepTime, amr->FluSgTime[lv][0], amr->FluSgTime[lv][1] );

//    print warning messages if temporal extrapolation is required
      const double TimeMin = MIN( amr->FluSgTime[lv][0], amr->FluSgTime[lv][1] );
      const double TimeMax = MAX( amr->FluSgTime[lv][0], amr->FluSgTime[lv][1] );

      if ( TimeMin < 0.0 )
         Aux_Error( ERROR_INFO, "TimeMin (%21.14e) < 0.0 ==> one of the fluid arrays has not been initialized !!\n", TimeMin );

      if ( PrepTime < TimeMin  ||  PrepTime-TimeMax >= 1.0e-12*TimeMax )
         Aux_Message( stderr, "WARNING : temporal extrapolation (lv %d, T_Prep %20.14e, T_Min %20.14e, T_Max %20.14e)\n",
                      lv, PrepTime, TimeMin, TimeMax );

      if ( OPT__INT_TIME )
      {
         FluIntTime        = true;
         FluSg             = 0;
         FluSg_IntT        = 1;
         FluWeighting      =   ( +amr->FluSgTime[lv][FluSg_IntT] - PrepTime )
                             / (  amr->FluSgTime[lv][FluSg_IntT] - amr->FluSgTime[lv][FluSg] );
         FluWeighting_IntT =   ( -amr->FluSgTime[lv][FluSg     ] + PrepTime )
                             / (  amr->FluSgTime[lv][FluSg_IntT] - amr->FluSgTime[lv][FluSg] );
      }

      else
      {
         FluIntTime        = false;
         FluSg             = amr->FluSg[lv]; // set to the current Sg
         FluSg_IntT        = NULL_INT;
         FluWeighting      = NULL_REAL;
         FluWeighting_IntT = NULL_REAL;
      }
   } // Mis_CompareRealValue
   } // if ( NVar_Flu + NVar_Der != 0 )

#  ifdef GRAVITY
   bool PotIntTime;
   int  PotSg, PotSg_IntT;
   real PotWeighting, PotWeighting_IntT;

   if ( PrepPot ) {
   if      (  Mis_CompareRealValue( PrepTime, amr->PotSgTime[lv][   amr->PotSg[lv] ], NULL, false )  )
   {
      PotIntTime        = false;
      PotSg             = amr->PotSg[lv];
      PotSg_IntT        = NULL_INT;
      PotWeighting      = NULL_REAL;
      PotWeighting_IntT = NULL_REAL;
   }

   else if (  Mis_CompareRealValue( PrepTime, amr->PotSgTime[lv][ 1-amr->PotSg[lv] ], NULL, false )  )
   {
      PotIntTime        = false;
      PotSg             = 1 - amr->PotSg[lv];
      PotSg_IntT        = NULL_INT;
      PotWeighting      = NULL_REAL;
      PotWeighting_IntT = NULL_REAL;
   }

   else
   {
//    check
      if ( OPT__DT_LEVEL == DT_LEVEL_SHARED )
      Aux_Error( ERROR_INFO, "cannot determine PotSg for OPT__DT_LEVEL == DT_LEVEL_SHARED (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                 lv, PrepTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );

//    print warning messages if temporal extrapolation is required
      const double TimeMin = MIN( amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );
      const double TimeMax = MAX( amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );

      if ( TimeMin < 0.0 )
         Aux_Error( ERROR_INFO, "TimeMin (%21.14e) < 0.0 ==> one of the potential arrays has not been initialized !!\n", TimeMin );

      if ( PrepTime < TimeMin  ||  PrepTime-TimeMax >= 1.0e-12*TimeMax )
         Aux_Message( stderr, "WARNING : temporal extrapolation (lv %d, T_Prep %20.14e, T_Min %20.14e, T_Max %20.14e)\n",
                      lv, PrepTime, TimeMin, TimeMax );

      if ( OPT__INT_TIME )
      {
         PotIntTime        = true;
         PotSg             = 0;
         PotSg_IntT        = 1;
         PotWeighting      =   ( +amr->PotSgTime[lv][PotSg_IntT] - PrepTime )
                             / (  amr->PotSgTime[lv][PotSg_IntT] - amr->PotSgTime[lv][PotSg] );
         PotWeighting_IntT =   ( -amr->PotSgTime[lv][PotSg     ] + PrepTime )
                             / (  amr->PotSgTime[lv][PotSg_IntT] - amr->PotSgTime[lv][PotSg] );
      }

      else
      {
         PotIntTime        = false;
         PotSg             = amr->PotSg[lv]; // set to the current Sg
         PotSg_IntT        = NULL_INT;
         PotWeighting      = NULL_REAL;
         PotWeighting_IntT = NULL_REAL;
      }
   } // Mis_CompareRealValue
   } // if ( PrepPot )
#  endif // #ifdef GRAVITY



// a. fill up the central region of CData
// ------------------------------------------------------------------------------------------------------------
   int i1, i2, j1, j2, k1, k2, Idx, TFluVarIdx, Disp1[3], Disp2[3], Loop1[3];
   double xyz[3];    // corner coordinates for the user-specified B.C.

   for (int d=0; d<3; d++)
   {
      Loop1[d] = TABLE_01( SibID, 'x'+d, La, PATCH_SIZE, La );
      Disp1[d] = TABLE_01( SibID, 'x'+d, PATCH_SIZE-La, 0, 0 );
      Disp2[d] = TABLE_01( SibID, 'x'+d, 0, CGhost, CGhost );
      xyz  [d] = TABLE_01( SibID, 'x'+d, amr->patch[0][lv][PID]->EdgeL[d] + (0.5+PS1-La)*dh,
                                         amr->patch[0][lv][PID]->EdgeL[d] + (0.5-CGhost)*dh,
                                         amr->patch[0][lv][PID]->EdgeL[d] + (0.5-CGhost)*dh );
   }


// a1. fluid data
   CData_Ptr = CData;

   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[TFluVarIdx][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                          + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][PID]->fluid[TFluVarIdx][k1][j1][i1];

         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }


// a2. derived variables
#  if   ( MODEL == HYDRO )
   if ( PrepVx )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[MOMX][k1][j1][i1] /
                          amr->patch[FluSg][lv][PID]->fluid[DENS][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                          + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMX][k1][j1][i1] /
                                                amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k1][j1][i1] );
         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }

   if ( PrepVy )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[MOMY][k1][j1][i1] /
                          amr->patch[FluSg][lv][PID]->fluid[DENS][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                          + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMY][k1][j1][i1] /
                                                amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k1][j1][i1] );
         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }

   if ( PrepVz )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[MOMZ][k1][j1][i1] /
                          amr->patch[FluSg][lv][PID]->fluid[DENS][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                          + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMZ][k1][j1][i1] /
                                                amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k1][j1][i1] );
         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }

   if ( PrepPres )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

         CData_Ptr[Idx] = CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                           Gamma_m1, CheckMinPres_No, NULL_REAL );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

            CData_Ptr[Idx] = FluWeighting     *CData_Ptr[Idx]
                           + FluWeighting_IntT*CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                                Gamma_m1, CheckMinPres_No, NULL_REAL );
         }

         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }

   if ( PrepTemp )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

         CData_Ptr[Idx] = CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                              Gamma_m1, (MinPres>=0.0), MinPres );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

            CData_Ptr[Idx] = FluWeighting     *CData_Ptr[Idx]
                           + FluWeighting_IntT*CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                                   Gamma_m1, (MinPres>=0.0), MinPres );
         }

         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }

#  elif ( MODEL == MHD   )
#  warning : WAIT MHD !!

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL


#  ifdef GRAVITY
// a3. potential data
   if ( PrepPot )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_Ptr[Idx] = amr->patch[PotSg][lv][PID]->pot[k1][j1][i1];

         if ( PotIntTime ) // temporal interpolation
         CData_Ptr[Idx] =   PotWeighting     *CData_Ptr[Idx]
                          + PotWeighting_IntT*amr->patch[PotSg_IntT][lv][PID]->pot[k1][j1][i1];

         Idx ++;
      }}}

      CData_Ptr += CSize3D;
   }
#  endif // #ifdef GRAVITY




// b. fill up the ghost zone of CData
// ------------------------------------------------------------------------------------------------------------
   int Loop2[3], Disp3[3], Disp4[3], Side, SibPID, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];

   for (int CSib=0; CSib<NTSib[SibID]; CSib++)
   {
      Side   = TSib[SibID][CSib];
      SibPID = amr->patch[0][lv][PID]->sibling[Side];

      for (int d=0; d<3; d++)
      {
         Loop2[d] = Table_01( SibID, Side, 'x'+d, La, CGhost, CGhost, PATCH_SIZE, CGhost, CGhost, La );
         Disp3[d] = Table_01( SibID, Side, 'x'+d, 0, La, 0, CGhost, CGhost+PATCH_SIZE, 0, CGhost );
         Disp4[d] = Table_01( SibID, Side, 'x'+d, PATCH_SIZE-La, 0, PATCH_SIZE-CGhost, 0, 0,
                              PATCH_SIZE-CGhost, 0 );
      }


//    b1. if the target sibling patch exists --> just copy data from the nearby patch at the same level
      if ( SibPID >= 0 )
      {
         CData_Ptr = CData;

//       b1-1. fluid data
         for (int v=0; v<NVar_Flu; v++)
         {
            TFluVarIdx = TFluVarIdxList[v];

            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[TFluVarIdx][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                                + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][SibPID]->fluid[TFluVarIdx][k2][j2][i2];

               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }


//       b1-2. derived variables
#        if   ( MODEL == HYDRO )
         if ( PrepVx )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[MOMX][k2][j2][i2] /
                                amr->patch[FluSg][lv][SibPID]->fluid[DENS][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                                + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMX][k2][j2][i2] /
                                                      amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][k2][j2][i2] );
               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }

         if ( PrepVy )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[MOMY][k2][j2][i2] /
                                amr->patch[FluSg][lv][SibPID]->fluid[DENS][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                                + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMY][k2][j2][i2] /
                                                      amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][k2][j2][i2] );
               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }

         if ( PrepVz )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[MOMZ][k2][j2][i2] /
                                amr->patch[FluSg][lv][SibPID]->fluid[DENS][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_Ptr[Idx] =   FluWeighting     *CData_Ptr[Idx]
                                + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMZ][k2][j2][i2] /
                                                      amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][k2][j2][i2] );
               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }

         if ( PrepPres )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

               CData_Ptr[Idx] = CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                 Gamma_m1, CheckMinPres_No, NULL_REAL );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

                  CData_Ptr[Idx] =  FluWeighting     *CData_Ptr[Idx]
                                  + FluWeighting_IntT*CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                                       Gamma_m1, CheckMinPres_No, NULL_REAL );
               }

               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }

         if ( PrepTemp )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

               CData_Ptr[Idx] = CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                    Gamma_m1, (MinPres>=0.0), MinPres );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

                  CData_Ptr[Idx] =  FluWeighting     *CData_Ptr[Idx]
                                  + FluWeighting_IntT*CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY],
                                                                          Fluid[MOMZ], Fluid[ENGY],
                                                                          Gamma_m1, (MinPres>=0.0), MinPres );
               }

               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }

#        elif ( MODEL == MHD   )
#        warning : WAIT MHD !!

#        elif ( MODEL == ELBDM )
//       no derived variables yet

#        else
#        error : unsupported MODEL !!
#        endif // MODEL


#        ifdef GRAVITY
//       b1-3. potential data
         if ( PrepPot )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_Ptr[Idx] = amr->patch[PotSg][lv][SibPID]->pot[k2][j2][i2];

               if ( PotIntTime ) // temporal interpolation
               CData_Ptr[Idx] =   PotWeighting     *CData_Ptr[Idx]
                                + PotWeighting_IntT*amr->patch[PotSg_IntT][lv][SibPID]->pot[k2][j2][i2];

               Idx ++;
            }}}

            CData_Ptr += CSize3D;
         }
#        endif // #ifdef GRAVITY
      } // if ( SibPID >= 0 )


//    b2. if the target sibling patch does not exist --> something is wrong !!
      else if ( SibPID == -1 )
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (Rank %d, Lv %d, PID %d, Side %d) !!\n",
                    "SibPID", SibPID, MPI_Rank, lv, PID, Side );


//    b3. if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
      else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
      {
         CData_Ptr = CData;

         for (int d=0; d<3; d++)
         {
            BC_Idx_Start[d] = Disp3[d];
            BC_Idx_End  [d] = Loop2[d] + BC_Idx_Start[d] - 1;
         }

         BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

#        ifdef GAMER_DEBUG
         if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
            Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

         if ( FluBC[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "FluBC == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID %d, PID %d, Side %d) !!\n",
                       BC_Sibling, BC_Face[BC_Sibling], SibPID, PID, Side );
#        endif

//       b3-1. fluid B.C.
         if ( NVar_Flu + NVar_Der > 0 )
         {
            switch ( FluBC[ BC_Face[BC_Sibling] ] )
            {
#              if ( MODEL == HYDRO  ||  MODEL == MHD )
               case BC_FLU_OUTFLOW:
                  Hydro_BoundaryCondition_Outflow   ( CData_Ptr, BC_Face[BC_Sibling], NVar_Flu+NVar_Der, CGhost,
                                                      CSize[0], CSize[1], CSize[2], BC_Idx_Start, BC_Idx_End );
               break;

               case BC_FLU_REFLECTING:
                  Hydro_BoundaryCondition_Reflecting( CData_Ptr, BC_Face[BC_Sibling], NVar_Flu,          CGhost,
                                                      CSize[0], CSize[1], CSize[2], BC_Idx_Start, BC_Idx_End,
                                                      TFluVarIdxList, NVar_Der, TDerVarList );
               break;
#              if ( MODEL == MHD )
#              warning : WAIT MHD !!!
#              endif
#              endif

               case BC_FLU_USER:
                  Flu_BoundaryCondition_User        ( CData_Ptr,                      NVar_Flu,
                                                      CSize[0], CSize[1], CSize[2], BC_Idx_Start, BC_Idx_End,
                                                      TFluVarIdxList, PrepTime, dh, xyz, TVar, lv );
               break;

               default:
                  Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
            } // switch ( FluBC[ BC_Face[BC_Sibling] ] )

            CData_Ptr += NVar_Flu*CSize3D;
         } // if ( NVar_Flu > 0 )

//       b3-2. potential B.C.
#        ifdef GRAVITY
         if ( PrepPot )
         {
//          extrapolate potential
            Poi_BoundaryCondition_Extrapolation( CData_Ptr, BC_Face[BC_Sibling], 1, CGhost,
                                                 CSize[0], CSize[1], CSize[2], BC_Idx_Start, BC_Idx_End );

            CData_Ptr += 1*CSize3D;
         }
#        endif // #ifdef GRAVITY

      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )

      else
         Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, Side %d) !!\n", SibPID, PID, Side );

   } // for (int CSib=0; CSib<NTSib[SibID]; CSib++)



// c. interpolation : CData --> IntData
// ------------------------------------------------------------------------------------------------------------
   const bool PhaseUnwrapping_Yes    = true;
   const bool PhaseUnwrapping_No     = false;
   const bool EnsureMonotonicity_Yes = true;
   const bool EnsureMonotonicity_No  = false;
   int CStart[3], CRange[3], FStart[3], NVar_SoFar;

   for (int d=0; d<3; d++)
   {
      CStart[d] = CGhost;
      CRange[d] = CSize[d] - 2*CGhost;
      FStart[d] = 0;
   }


// determine which variables require **monotonic** interpolation
   bool Monotonicity[NVar_Flu];

   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v];

#     if ( MODEL == HYDRO )
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
      /*
      if ( TFluVarIdx == DENS  ||  TFluVarIdx == ENGY  ||  TFluVarIdx >= NCOMP_FLUID )
                                                         Monotonicity[v] = EnsureMonotonicity_Yes;
      else                                               Monotonicity[v] = EnsureMonotonicity_No;
      */
                                                         Monotonicity[v] = EnsureMonotonicity_Yes;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
//    apply monotonic interpolation to density and all passive scalars
      if ( TFluVarIdx != REAL  &&  TFluVarIdx != IMAG )  Monotonicity[v] = EnsureMonotonicity_Yes;
      else                                               Monotonicity[v] = EnsureMonotonicity_No;

#     else
#     error : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION IN THIS NEW MODEL ??
#     endif // MODEL
   }


// interpolation
#  if ( MODEL == ELBDM )
   real *CData_Dens = NULL;
   real *CData_Real = NULL;
   real *CData_Imag = NULL;
   real *FData_Dens = NULL;
   real *FData_Real = NULL;
   real *FData_Imag = NULL;

// c1. interpolation on phase in ELBDM
   if ( IntPhase )
   {
      int DensIdx=-1, RealIdx=-1, ImagIdx=-1;

      for (int v=0; v<NVar_Flu; v++)
      {
         TFluVarIdx = TFluVarIdxList[v];

         if      ( TFluVarIdx == DENS )   DensIdx = v;
         else if ( TFluVarIdx == REAL )   RealIdx = v;
         else if ( TFluVarIdx == IMAG )   ImagIdx = v;
      }

//    check
#     ifdef GAMER_DEBUG
      if ( RealIdx == -1  ||  ImagIdx == -1 )
         Aux_Error( ERROR_INFO, "real and/or imag parts are not found for phase interpolation in ELBDM !!\n" );
#     endif

//    determine the array index to store density
      CData_Dens = CData   + ( (DensIdx==-1) ? ImagIdx : DensIdx )*CSize3D;
      CData_Real = CData   + RealIdx*CSize3D;
      CData_Imag = CData   + ImagIdx*CSize3D;
      FData_Dens = IntData + ( (DensIdx==-1) ? ImagIdx : DensIdx )*FSize3D;
      FData_Real = IntData + RealIdx*FSize3D;
      FData_Imag = IntData + ImagIdx*FSize3D;

//    get the wrapped phase (store in the REAL component) and density (store in the IMAG component)
      real Re, Im;

      for (int t=0; t<CSize3D; t++)
      {
         Re = CData_Real[t];
         Im = CData_Imag[t];

         CData_Real[t] = ATAN2( Im, Re );
         if ( DensIdx == -1 )
         CData_Dens[t] = Re*Re + Im*Im;
      }

//    interpolate density
      Interpolate( CData_Dens, CSize, CStart, CRange, FData_Dens, FSize, FStart, 1, IntScheme,
                   PhaseUnwrapping_No, &EnsureMonotonicity_Yes );

//    interpolate phase
      Interpolate( CData_Real, CSize, CStart, CRange, FData_Real, FSize, FStart, 1, IntScheme,
                   PhaseUnwrapping_Yes, &EnsureMonotonicity_No );
   } // if ( IntPhase )


// c2. interpolation on real/imag parts in ELBDM
   else // if ( IntPhase )
   {
      for (int v=0; v<NVar_Flu; v++)
      Interpolate( CData+CSize3D*v, CSize, CStart, CRange, IntData+FSize3D*v, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, Monotonicity );
   } // if ( IntPhase ) ... else ...

// retrieve real and imaginary parts when phase interpolation is adopted
   if ( IntPhase )
   {
      real Amp, Phase, Rho;

      for (int t=0; t<FSize3D; t++)
      {
         Phase = FData_Real[t];
         Rho   = FData_Dens[t];

//       be careful about the negative density introduced from the round-off errors
//       --> note that we check minimum density in the end of Prepare_PatchData()
         if ( Rho < (real)0.0 )
         {
            FData_Dens[t] = (real)0.0;
            Rho           = (real)0.0;
         }

         Amp           = SQRT( Rho );
         FData_Real[t] = Amp*COS( Phase );
         FData_Imag[t] = Amp*SIN( Phase );
      }
   }

#  else // #if ( MODEL == ELBDM )


// c3. interpolation on original variables for models != ELBDM
   for (int v=0; v<NVar_Flu; v++)
      Interpolate( CData+CSize3D*v, CSize, CStart, CRange, IntData+FSize3D*v, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, Monotonicity );

#  endif // #if ( MODEL == ELBDM ) ... else ...

   NVar_SoFar = NVar_Flu;


// c4. derived variables
#  if   ( MODEL == HYDRO )
// we now apply monotonic interpolation to ALL fluid variables
   if ( PrepVx )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVar_SoFar ++;
   }

   if ( PrepVy )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVar_SoFar ++;
   }

   if ( PrepVz )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVar_SoFar ++;
   }

   if ( PrepPres )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVar_SoFar ++;
   }

   if ( PrepTemp )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVar_SoFar ++;
   }

#  elif ( MODEL == MHD   )
#  warning : WAIT MHD !!

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL


#  ifdef GRAVITY
// c5. interpolation on potential
   if ( PrepPot )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1,
                   IntScheme, PhaseUnwrapping_No, &EnsureMonotonicity_No );
      NVar_SoFar ++;
   }
#  endif

   delete [] CData;


// d. ensure the consistency between pressure, total energy density, and the dual-energy variable
//    when DUAL_ENERGY is on
//    --> we don't have to check the minimum pressure here when DUAL_ENERGY is off
//        --> it's checked in Prepare_PatchData()
// ------------------------------------------------------------------------------------------------------------
#  if (  ( MODEL == HYDRO || MODEL == MHD )  &&  defined DUAL_ENERGY  )
// apply this correction only when preparing all fluid variables
   if (  ( TVar & _TOTAL ) == _TOTAL  &&  DE_Consistency )
   {
      const real _Gamma_m1       = (real)1.0 / Gamma_m1;
      const real UseEnpy2FixEngy = HUGE_NUMBER;

//    assuming that the order of variables stored in IntData is the same as patch->fluid[]
      real *FData_Dens = IntData + DENS*FSize3D;
      real *FData_MomX = IntData + MOMX*FSize3D;
      real *FData_MomY = IntData + MOMY*FSize3D;
      real *FData_MomZ = IntData + MOMZ*FSize3D;
      real *FData_Engy = IntData + ENGY*FSize3D;
      real *FData_Enpy = IntData + ENPY*FSize3D;

      char dummy;    // we do not record the dual-energy status here

      for (int t=0; t<FSize3D; t++)
      {
//       here we ALWAYS use the dual-energy variable to correct the total energy density
//       --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//           the runtime parameter DUAL_ENERGY_SWITCH here
         CPU_DualEnergyFix( FData_Dens[t], FData_MomX[t], FData_MomY[t], FData_MomZ[t], FData_Engy[t], FData_Enpy[t],
                            dummy, Gamma_m1, _Gamma_m1, (MinPres>=0.0), MinPres, UseEnpy2FixEngy );
      }
   } // if (  ( TVar & _TOTAL ) == _TOTAL  &&  DE_Consistency  )
#  endif // if (  ( MODEL == HYDRO || MODEL == MHD )  &&  defined DUAL_ENERGY  )

} // FUNCTION : InterpolateGhostZone



// ============
// |  Tables  |
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  return the loop size and displacement required by the function "InterpolateGhostZone"
//
// Parameter   :  SibID       : Sibling index ONE (0~25)
//                Side        : Sibling index TWO (0~25)
//                dim         : Target spatial direction (x/y/z)
//                w01 ... w21 : Returned values
//-------------------------------------------------------------------------------------------------------
int Table_01( const int SibID, const int Side, const char dim, const int w01, const int w02,
              const int w10, const int w11, const int w12, const int w20, const int w21 )
{

   switch ( dim )
   {
      case 'x':
      {
         switch ( SibID )
         {
            case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
            {
               switch ( Side )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Side", Side );

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w01;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     return w02;
               }
            }

            case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
            {
               switch ( Side )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     return w10;

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w11;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     return w12;
               }
            }

            case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
            {
               switch ( Side )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     return w20;

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w21;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Side", Side );
               }
            }

         } // switch ( SibID )
      } // case 'x':


      case 'y':
      {
         switch ( SibID )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
            {
               switch ( Side )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Side", Side );

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w01;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     return w02;
               }
            }

            case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
            {
               switch ( Side )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     return w10;

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w11;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     return w12;
               }
            }

            case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
            {
               switch ( Side )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     return w20;

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w21;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Side", Side );
               }
            }

         } // switch ( SibID )
      } // case 'y':


      case 'z':
      {
         switch ( SibID )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
            {
               switch ( Side )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Side", Side );

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w01;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     return w02;
               }
            }

            case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
            {
               switch ( Side )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     return w10;

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w11;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     return w12;
               }
            }

            case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
            {
               switch ( Side )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     return w20;

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w21;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Side", Side );
               }
            }

         } // switch ( SibID )
      } // case 'z':


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c !!\n", "dim", dim );
         exit(1);

   } // switch ( dim )

   return NULL_INT;

} // FUNCTION : Table_01
