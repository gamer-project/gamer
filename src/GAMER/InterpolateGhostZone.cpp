#include "GAMER.h"

static int Table_01( const int SibID, const int Side, const char dim, const int w01, const int w02,
                     const int w10, const int w11, const int w12, const int w20, const int w21 );
void SetTempIntPara( const int lv, const int Sg_Current, const double PrepTime, const double Time0, const double Time1,
                     bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );




//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolateGhostZone
// Description :  Fill up the ghost-zone values by spatial and temporal interpolation
//
// Note        :  1. Work for Prepare_PatchData()
//                2. Use the input parameter "TVarCC" and "TVarFC" to control the target variables
//                3. Use the input parameters "IntScheme_CC/FC" to control the interpolation schemes for the
//                   cell-/face-centered variables
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
// Parameter   :  lv                : Target "coarse-grid" refinement level
//                PID               : Patch ID at level "lv" used for interpolation
//                IntData_CC/FC     : Arrays to store the cell-/face-centered interpolation results
//                SibID             : Sibling index (0~25) used to determine the interpolation region
//                PrepTime          : Target physical time to prepare data
//                GhostSize         : Number of ghost zones
//                IntScheme_CC      : Interpolation scheme for the cell-centered variables
//                                    --> Supported schemes include
//                                        INT_MINMOD1D : MinMod-1D
//                                        INT_MINMOD3D : MinMod-3D
//                                        INT_VANLEER  : vanLeer
//                                        INT_CQUAD    : conservative quadratic
//                                        INT_QUAD     : quadratic
//                                        INT_CQUAR    : conservative quartic
//                                        INT_QUAR     : quartic
//                IntScheme_FC      : Interpolation scheme for the face-centered variables
//                                    --> Supported schemes include
//                                        INT_MINMOD1D : MinMod-1D
//                                        INT_VANLEER  : vanLeer
//                                        INT_CQUAD    : conservative quadratic
//                                        INT_CQUAR    : conservative quartic
//                NTSib             : Number of target sibling patches along different sibling directions
//                TSib              : Target sibling indices along different sibling directions
//                TVarCC            : Target cell-centered variables to be prepared
//                                    --> Supported variables in different models:
//                                        HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELX, _VELY, _VELZ, _PRES, _TEMP,
//                                                [, _POTE]
//                                        ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                                    --> _FLUID, _PASSIVE, _TOTAL, and _DERIVED apply to all models
//                NVarCC_Tot        : Total number of cell-centered variables to be prepared
//                NVarCC_Flu        : Number of cell-centered fluid variables to be prepared
//                                    --> Including passive scalars
//                TVarCCIdxList_Flu : List recording the target cell-centered fluid and passive variable indices
//                                    ( = [0 ... NCOMP_TOTAL-1] )
//                NVarCC_Der        : Number of cell-centered derived variables to be prepared
//                TVarCCList_Der    : List recording the target cell-centered derived variables
//                TVarFC            : Target face-centered variables to be prepared
//                                    --> Supported variables in different models:
//                                        HYDRO with MHD : _MAG
//                                        ELBDM          : none
//                                    --> Note that for MHD it does NOT support preparing individual B component
//                                        (e.g., _MAGX, _MAGY, _MAGZ) since the divergenece-free interpolation
//                                        must work on all three components at once
//                NVarFC_Tot        : Total number of face-centered variables to be prepared
//                TVarFCIdxList     : List recording the target face-centered variable indices
//                                    ( = [0 ... NCOMP_MAGNETIC-1] )
//                IntPhase          : true --> Perform interpolation on rho/phase instead of real/imag parts in ELBDM
//                FluBC             : Fluid boundary condition
//                PotBC             : Gravity boundary condition (not used currently)
//                BC_Face           : Priority of the B.C. along different boundary faces (z>y>x)
//                MinPres           : Minimum allowed pressure
//                DE_Consistency    : Ensure the consistency between pressure, total energy density, and the
//                                    dual-energy variable when DUAL_ENERGY is on
//-------------------------------------------------------------------------------------------------------
void InterpolateGhostZone( const int lv, const int PID, real IntData_CC[], real IntData_FC[],
                           const int SibID, const double PrepTime, const int GhostSize,
                           const IntScheme_t IntScheme_CC, const IntScheme_t IntScheme_FC,
                           const int NTSib[], int *TSib[], const int TVarCC, const int NVarCC_Tot,
                           const int NVarCC_Flu, const int TVarCCIdxList_Flu[],
                           const int NVarCC_Der, const int TVarCCList_Der[],
                           const int TVarFC, const int NVarFC_Tot, const int TVarFCIdxList[],
                           const bool IntPhase, const OptFluBC_t FluBC[], const OptPotBC_t PotBC,
                           const int BC_Face[], const real MinPres, const bool DE_Consistency )
{

// check
#  ifdef MHD
   Aux_Error( ERROR_INFO, "%s does not support MHD yet !!\n", __FUNCTION__ );
#  endif
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
   const bool PrepVx          = ( TVarCC & _VELX ) ? true : false;
   const bool PrepVy          = ( TVarCC & _VELY ) ? true : false;
   const bool PrepVz          = ( TVarCC & _VELZ ) ? true : false;
   const bool PrepPres        = ( TVarCC & _PRES ) ? true : false;
   const bool PrepTemp        = ( TVarCC & _TEMP ) ? true : false;
#  ifdef MHD
#  warning : WAIT MHD !!!
#  endif

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef GRAVITY
   const bool PrepPot         = ( TVarCC & _POTE ) ? true : false;
#  endif

#  if ( MODEL == HYDRO )
   real Fluid[NCOMP_FLUID];   // for calculating pressure and temperature only --> don't need NCOMP_TOTAL
#  endif


// set up parameters for the adopted interpolation scheme
   int NSide_CC, CGhost_CC, CSize_CC[3], FSize_CC[3], CSize3D_CC, FSize3D_CC;

   Int_Table( IntScheme_CC, NSide_CC, CGhost_CC );

   const double dh               = amr->dh[lv];
   const int    GhostSize_Padded = GhostSize + (GhostSize&1);        // # of interpolated cells must an even number
   const int    CGrid            = GhostSize_Padded/2 + 2*CGhost_CC; // # of coarse cells required for interpolation
   const int    CGrid_PID        = CGrid - CGhost_CC;                // # of overlapped cells between CGrid and PID

   for (int d=0; d<3; d++)
   {
      CSize_CC[d] = TABLE_01( SibID, 'x'+d, CGrid, PS1+2*CGhost_CC, CGrid );
      FSize_CC[d] = TABLE_01( SibID, 'x'+d, GhostSize_Padded, PS2, GhostSize_Padded );
   }

   CSize3D_CC = CSize_CC[0]*CSize_CC[1]*CSize_CC[2];
   FSize3D_CC = FSize_CC[0]*FSize_CC[1]*FSize_CC[2];


// we assume that we only need ONE coarse-grid patch in each sibling direction
   if ( CGrid_PID > PS1 )
      Aux_Error( ERROR_INFO, "CGrid_PID (%d) > PATCH_SIZE (%d) !!\n", CGrid_PID, PS1 );


// coarse-grid array to store all the data required for interpolation (including the ghost zones in each side)
   real *CData_CC_Ptr = NULL;
   real *CData_CC     = new real [ NVarCC_Tot*CSize3D_CC ];


// temporal interpolation parameters
   bool FluIntTime;
   int  FluSg, FluSg_IntT;
   real FluWeighting, FluWeighting_IntT;

// fluid
   if ( NVarCC_Flu + NVarCC_Der != 0 )
   {
      SetTempIntPara( lv, amr->FluSg[lv], PrepTime, amr->FluSgTime[lv][0], amr->FluSgTime[lv][1],
                      FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );

      if ( FluIntTime  &&  OPT__DT_LEVEL == DT_LEVEL_SHARED )
         Aux_Error( ERROR_INFO, "cannot determine FluSg for OPT__DT_LEVEL == DT_LEVEL_SHARED "
                                "(lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                    lv, PrepTime, amr->FluSgTime[lv][0], amr->FluSgTime[lv][1] );
   }

// magnetic field
#  ifdef MHD
   bool MagIntTime;
   int  MagSg, MagSg_IntT;
   real MagWeighting, MagWeighting_IntT;
#  ifdef MHD
#  warning : REMOVE THESE !!!!!!!!!!!!!!!!!
   int NVar_MagFC=0, NVar_MagCC=0;
#  endif

// check NVarCC_Der as well since we might need B field for calculating, for example, gas pressure
   if ( NVar_MagFC + NVar_MagCC + NVarCC_Der != 0 )
   {
      SetTempIntPara( lv, amr->MagSg[lv], PrepTime, amr->MagSgTime[lv][0], amr->MagSgTime[lv][1],
                      MagIntTime, MagSg, MagSg_IntT, MagWeighting, MagWeighting_IntT );

      if ( MagIntTime  &&  OPT__DT_LEVEL == DT_LEVEL_SHARED )
         Aux_Error( ERROR_INFO, "cannot determine MagSg for OPT__DT_LEVEL == DT_LEVEL_SHARED "
                                "(lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                    lv, PrepTime, amr->MagSgTime[lv][0], amr->MagSgTime[lv][1] );
   }
#  endif // #ifdef MHD

// potential
#  ifdef GRAVITY
   bool PotIntTime;
   int  PotSg, PotSg_IntT;
   real PotWeighting, PotWeighting_IntT;

   if ( PrepPot )
   {
      SetTempIntPara( lv, amr->PotSg[lv], PrepTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1],
                      PotIntTime, PotSg, PotSg_IntT, PotWeighting, PotWeighting_IntT );

      if ( PotIntTime  &&  OPT__DT_LEVEL == DT_LEVEL_SHARED )
         Aux_Error( ERROR_INFO, "cannot determine PotSg for OPT__DT_LEVEL == DT_LEVEL_SHARED "
                                "(lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                    lv, PrepTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );
   }
#  endif // #ifdef GRAVITY



// a. fill up the central region of CData_CC[] and CData_FC[]
// ------------------------------------------------------------------------------------------------------------
   int i1, i2, j1, j2, k1, k2, Idx, TVarCCIdx_Flu, Disp1[3], Disp2[3], Loop1[3];
   double xyz[3];    // corner coordinates for the user-specified B.C.

   for (int d=0; d<3; d++)
   {
      Loop1[d] = TABLE_01( SibID, 'x'+d, CGrid_PID, PS1, CGrid_PID );
      Disp1[d] = TABLE_01( SibID, 'x'+d, PS1-CGrid_PID, 0, 0 );
      Disp2[d] = TABLE_01( SibID, 'x'+d, 0, CGhost_CC, CGhost_CC );
      xyz  [d] = TABLE_01( SibID, 'x'+d, amr->patch[0][lv][PID]->EdgeL[d] + (0.5+PS1-CGrid_PID)*dh,
                                         amr->patch[0][lv][PID]->EdgeL[d] + (0.5-CGhost_CC)*dh,
                                         amr->patch[0][lv][PID]->EdgeL[d] + (0.5-CGhost_CC)*dh );
   }


// a1. fluid data
   CData_CC_Ptr = CData_CC;

   for (int v=0; v<NVarCC_Flu; v++)
   {
      TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[TVarCCIdx_Flu][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                             + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][PID]->fluid[TVarCCIdx_Flu][k1][j1][i1];

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }


// a2. derived variables
#  if   ( MODEL == HYDRO )
   if ( PrepVx )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[MOMX][k1][j1][i1] /
                             amr->patch[FluSg][lv][PID]->fluid[DENS][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                             + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMX][k1][j1][i1] /
                                                   amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k1][j1][i1] );
         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }

   if ( PrepVy )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[MOMY][k1][j1][i1] /
                             amr->patch[FluSg][lv][PID]->fluid[DENS][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                             + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMY][k1][j1][i1] /
                                                   amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k1][j1][i1] );
         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }

   if ( PrepVz )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][PID]->fluid[MOMZ][k1][j1][i1] /
                             amr->patch[FluSg][lv][PID]->fluid[DENS][k1][j1][i1];

         if ( FluIntTime ) // temporal interpolation
         CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                             + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMZ][k1][j1][i1] /
                                                   amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k1][j1][i1] );
         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }

   if ( PrepPres )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

#        ifdef MHD
#        warning : WAIT MHD !!!
         const real EngyB = NULL_REAL;
#        else
         const real EngyB = NULL_REAL;
#        endif
         CData_CC_Ptr[Idx] = CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                              Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

#           ifdef MHD
#           warning : WAIT MHD !!!
            const real EngyB = NULL_REAL;
#           else
            const real EngyB = NULL_REAL;
#           endif
            CData_CC_Ptr[Idx] = FluWeighting     *CData_CC_Ptr[Idx]
                              + FluWeighting_IntT*CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                                   Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB );
         }

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }

   if ( PrepTemp )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

#        ifdef MHD
#        warning : WAIT MHD !!!
         const real EngyB = NULL_REAL;
#        else
         const real EngyB = NULL_REAL;
#        endif
         CData_CC_Ptr[Idx] = CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                 Gamma_m1, (MinPres>=0.0), MinPres, EngyB );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

#           ifdef MHD
#           warning : WAIT MHD !!!
            const real EngyB = NULL_REAL;
#           else
            const real EngyB = NULL_REAL;
#           endif
            CData_CC_Ptr[Idx] = FluWeighting     *CData_CC_Ptr[Idx]
                              + FluWeighting_IntT*CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ],
                                                                      Fluid[ENGY], Gamma_m1, (MinPres>=0.0), MinPres, EngyB );
         }

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }

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
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         CData_CC_Ptr[Idx] = amr->patch[PotSg][lv][PID]->pot[k1][j1][i1];

         if ( PotIntTime ) // temporal interpolation
         CData_CC_Ptr[Idx] =   PotWeighting     *CData_CC_Ptr[Idx]
                             + PotWeighting_IntT*amr->patch[PotSg_IntT][lv][PID]->pot[k1][j1][i1];

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   }
#  endif // #ifdef GRAVITY




// b. fill up the ghost zone of CData_CC[] and CData_FC[]
// ------------------------------------------------------------------------------------------------------------
   int Loop2[3], Disp3[3], Disp4[3], Side, SibPID, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];

   for (int CSib=0; CSib<NTSib[SibID]; CSib++)
   {
      Side   = TSib[SibID][CSib];
      SibPID = amr->patch[0][lv][PID]->sibling[Side];

      for (int d=0; d<3; d++)
      {
         Loop2[d] = Table_01( SibID, Side, 'x'+d, CGrid_PID, CGhost_CC, CGhost_CC, PS1, CGhost_CC, CGhost_CC, CGrid_PID );
         Disp3[d] = Table_01( SibID, Side, 'x'+d, 0, CGrid_PID, 0, CGhost_CC, CGhost_CC+PS1, 0, CGhost_CC );
         Disp4[d] = Table_01( SibID, Side, 'x'+d, PS1-CGrid_PID, 0, PS1-CGhost_CC, 0, 0, PS1-CGhost_CC, 0 );
      }


//    b1. if the target sibling patch exists --> just copy data from the nearby patch at the same level
      if ( SibPID >= 0 )
      {
         CData_CC_Ptr = CData_CC;

//       b1-1. fluid data
         for (int v=0; v<NVarCC_Flu; v++)
         {
            TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[TVarCCIdx_Flu][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                                   + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][SibPID]->fluid[TVarCCIdx_Flu][k2][j2][i2];

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }


//       b1-2. derived variables
#        if   ( MODEL == HYDRO )
         if ( PrepVx )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[MOMX][k2][j2][i2] /
                                   amr->patch[FluSg][lv][SibPID]->fluid[DENS][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                                   + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMX][k2][j2][i2] /
                                                         amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][k2][j2][i2] );
               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }

         if ( PrepVy )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[MOMY][k2][j2][i2] /
                                   amr->patch[FluSg][lv][SibPID]->fluid[DENS][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                                   + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMY][k2][j2][i2] /
                                                         amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][k2][j2][i2] );
               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }

         if ( PrepVz )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_CC_Ptr[Idx] = amr->patch[FluSg][lv][SibPID]->fluid[MOMZ][k2][j2][i2] /
                                   amr->patch[FluSg][lv][SibPID]->fluid[DENS][k2][j2][i2];

               if ( FluIntTime ) // temporal interpolation
               CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                                   + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMZ][k2][j2][i2] /
                                                         amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][k2][j2][i2] );
               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }

         if ( PrepPres )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

#              ifdef MHD
#              warning : WAIT MHD !!!
               const real EngyB = NULL_REAL;
#              else
               const real EngyB = NULL_REAL;
#              endif
               CData_CC_Ptr[Idx] = CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                    Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

#                 ifdef MHD
#                 warning : WAIT MHD !!!
                  const real EngyB = NULL_REAL;
#                 else
                  const real EngyB = NULL_REAL;
#                 endif
                  CData_CC_Ptr[Idx] =  FluWeighting     *CData_CC_Ptr[Idx]
                                     + FluWeighting_IntT*CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ],
                                                                          Fluid[ENGY], Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB );
               }

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }

         if ( PrepTemp )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

#              ifdef MHD
#              warning : WAIT MHD !!!
               const real EngyB = NULL_REAL;
#              else
               const real EngyB = NULL_REAL;
#              endif
               CData_CC_Ptr[Idx] = CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                       Gamma_m1, (MinPres>=0.0), MinPres, EngyB );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NCOMP_FLUID; v++)   Fluid[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

#                 ifdef MHD
#                 warning : WAIT MHD !!!
                  const real EngyB = NULL_REAL;
#                 else
                  const real EngyB = NULL_REAL;
#                 endif
                  CData_CC_Ptr[Idx] =  FluWeighting     *CData_CC_Ptr[Idx]
                                     + FluWeighting_IntT*CPU_GetTemperature( Fluid[DENS], Fluid[MOMX], Fluid[MOMY],
                                                                             Fluid[MOMZ], Fluid[ENGY],
                                                                             Gamma_m1, (MinPres>=0.0), MinPres, EngyB );
               }

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }

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
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               CData_CC_Ptr[Idx] = amr->patch[PotSg][lv][SibPID]->pot[k2][j2][i2];

               if ( PotIntTime ) // temporal interpolation
               CData_CC_Ptr[Idx] =   PotWeighting     *CData_CC_Ptr[Idx]
                                   + PotWeighting_IntT*amr->patch[PotSg_IntT][lv][SibPID]->pot[k2][j2][i2];

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
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
         CData_CC_Ptr = CData_CC;

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
         if ( NVarCC_Flu + NVarCC_Der > 0 )
         {
            switch ( FluBC[ BC_Face[BC_Sibling] ] )
            {
               case BC_FLU_OUTFLOW:
                  Flu_BoundaryCondition_Outflow     ( CData_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu+NVarCC_Der, CGhost_CC,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End );
               break;

#              if ( MODEL == HYDRO )
               case BC_FLU_REFLECTING:
                  Hydro_BoundaryCondition_Reflecting( CData_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu,          CGhost_CC,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End,
                                                      TVarCCIdxList_Flu, NVarCC_Der, TVarCCList_Der );
               break;
#              endif

               case BC_FLU_USER:
                  Flu_BoundaryCondition_User        ( CData_CC_Ptr,                      NVarCC_Flu,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End,
                                                      TVarCCIdxList_Flu, PrepTime, dh, xyz, TVarCC, lv );
               break;

               default:
                  Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
            } // switch ( FluBC[ BC_Face[BC_Sibling] ] )

            CData_CC_Ptr += NVarCC_Flu*CSize3D_CC;
         } // if ( NVarCC_Flu > 0 )

//       b3-2. potential B.C.
#        ifdef GRAVITY
         if ( PrepPot )
         {
//          extrapolate potential
            Poi_BoundaryCondition_Extrapolation( CData_CC_Ptr, BC_Face[BC_Sibling], 1, CGhost_CC,
                                                 CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End );

            CData_CC_Ptr += 1*CSize3D_CC;
         }
#        endif // #ifdef GRAVITY

      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )

      else
         Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, Side %d) !!\n", SibPID, PID, Side );

   } // for (int CSib=0; CSib<NTSib[SibID]; CSib++)



// c. interpolation : CData_CC[] --> IntData_CC[]
// ------------------------------------------------------------------------------------------------------------
   const bool PhaseUnwrapping_Yes    = true;
   const bool PhaseUnwrapping_No     = false;
   const bool EnsureMonotonicity_Yes = true;
   const bool EnsureMonotonicity_No  = false;
   int CStart[3], CRange[3], FStart[3], NVarCC_SoFar;

   for (int d=0; d<3; d++)
   {
      CStart[d] = CGhost_CC;
      CRange[d] = CSize_CC[d] - 2*CGhost_CC;
      FStart[d] = 0;
   }


// determine which variables require **monotonic** interpolation
   bool Monotonicity[NVarCC_Flu];

   for (int v=0; v<NVarCC_Flu; v++)
   {
      TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

#     if   ( MODEL == HYDRO )
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
      /*
      if ( TVarCCIdx_Flu == DENS  ||  TVarCCIdx_Flu == ENGY  ||  TVarCCIdx_Flu >= NCOMP_FLUID )
         Monotonicity[v] = EnsureMonotonicity_Yes;
      else
         Monotonicity[v] = EnsureMonotonicity_No;
      */
         Monotonicity[v] = EnsureMonotonicity_Yes;

#     elif ( MODEL == ELBDM )
//    apply monotonic interpolation to density and all passive scalars
      if ( TVarCCIdx_Flu != REAL  &&  TVarCCIdx_Flu != IMAG )
         Monotonicity[v] = EnsureMonotonicity_Yes;
      else
         Monotonicity[v] = EnsureMonotonicity_No;

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

      for (int v=0; v<NVarCC_Flu; v++)
      {
         TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

         if      ( TVarCCIdx_Flu == DENS )   DensIdx = v;
         else if ( TVarCCIdx_Flu == REAL )   RealIdx = v;
         else if ( TVarCCIdx_Flu == IMAG )   ImagIdx = v;
      }

//    check
#     ifdef GAMER_DEBUG
      if ( RealIdx == -1  ||  ImagIdx == -1 )
         Aux_Error( ERROR_INFO, "real and/or imag parts are not found for phase interpolation in ELBDM !!\n" );
#     endif

//    determine the array index to store density
      CData_Dens = CData_CC   + ( (DensIdx==-1) ? ImagIdx : DensIdx )*CSize3D_CC;
      CData_Real = CData_CC   + RealIdx*CSize3D_CC;
      CData_Imag = CData_CC   + ImagIdx*CSize3D_CC;
      FData_Dens = IntData_CC + ( (DensIdx==-1) ? ImagIdx : DensIdx )*FSize3D_CC;
      FData_Real = IntData_CC + RealIdx*FSize3D_CC;
      FData_Imag = IntData_CC + ImagIdx*FSize3D_CC;

//    get the wrapped phase (store in the REAL component) and density (store in the IMAG component)
      real Re, Im;

      for (int t=0; t<CSize3D_CC; t++)
      {
         Re = CData_Real[t];
         Im = CData_Imag[t];

         CData_Real[t] = ATAN2( Im, Re );
         if ( DensIdx == -1 )
         CData_Dens[t] = Re*Re + Im*Im;
      }

//    interpolate density
      Interpolate( CData_Dens, CSize_CC, CStart, CRange, FData_Dens, FSize_CC, FStart, 1, IntScheme_CC,
                   PhaseUnwrapping_No, &EnsureMonotonicity_Yes );

//    interpolate phase
      Interpolate( CData_Real, CSize_CC, CStart, CRange, FData_Real, FSize_CC, FStart, 1, IntScheme_CC,
                   PhaseUnwrapping_Yes, &EnsureMonotonicity_No );
   } // if ( IntPhase )


// c2. interpolation on real/imag parts in ELBDM
   else // if ( IntPhase )
   {
      for (int v=0; v<NVarCC_Flu; v++)
      Interpolate( CData_CC+CSize3D_CC*v, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*v, FSize_CC, FStart, 1,
                   IntScheme_CC, PhaseUnwrapping_No, Monotonicity );
   } // if ( IntPhase ) ... else ...

// retrieve real and imaginary parts when phase interpolation is adopted
   if ( IntPhase )
   {
      real Amp, Phase, Rho;

      for (int t=0; t<FSize3D_CC; t++)
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
   for (int v=0; v<NVarCC_Flu; v++)
      Interpolate( CData_CC+CSize3D_CC*v, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*v, FSize_CC, FStart, 1,
                   IntScheme_CC, PhaseUnwrapping_No, Monotonicity );

#  endif // #if ( MODEL == ELBDM ) ... else ...

   NVarCC_SoFar = NVarCC_Flu;


// c4. derived variables
#  if   ( MODEL == HYDRO )
// we now apply monotonic interpolation to ALL fluid variables
   if ( PrepVx )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*NVarCC_SoFar,
                   FSize_CC, FStart, 1, IntScheme_CC, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVarCC_SoFar ++;
   }

   if ( PrepVy )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*NVarCC_SoFar,
                   FSize_CC, FStart, 1, IntScheme_CC, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVarCC_SoFar ++;
   }

   if ( PrepVz )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*NVarCC_SoFar,
                   FSize_CC, FStart, 1, IntScheme_CC, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVarCC_SoFar ++;
   }

   if ( PrepPres )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*NVarCC_SoFar,
                   FSize_CC, FStart, 1, IntScheme_CC, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVarCC_SoFar ++;
   }

   if ( PrepTemp )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*NVarCC_SoFar,
                   FSize_CC, FStart, 1, IntScheme_CC, PhaseUnwrapping_No, &EnsureMonotonicity_Yes );
      NVarCC_SoFar ++;
   }

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL


#  ifdef GRAVITY
// c5. interpolation on potential
   if ( PrepPot )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart, CRange, IntData_CC+FSize3D_CC*NVarCC_SoFar,
                   FSize_CC, FStart, 1, IntScheme_CC, PhaseUnwrapping_No, &EnsureMonotonicity_No );
      NVarCC_SoFar ++;
   }
#  endif

   delete [] CData_CC;


// d. ensure the consistency between pressure, total energy density, and the dual-energy variable
//    when DUAL_ENERGY is on
//    --> we don't have to check the minimum pressure here when DUAL_ENERGY is off
//        --> it's checked in Prepare_PatchData()
// ------------------------------------------------------------------------------------------------------------
#  if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
#  ifdef MHD
#  warning : WAIT MHD !!!
#  endif
// apply this correction only when preparing all fluid variables
   if (  ( TVarCC & _TOTAL ) == _TOTAL  &&  DE_Consistency )
   {
      const real _Gamma_m1       = (real)1.0 / Gamma_m1;
      const real UseEnpy2FixEngy = HUGE_NUMBER;

//    assuming that the order of variables stored in IntData_CC[] is the same as patch->fluid[]
      real *FData_Dens = IntData_CC + DENS*FSize3D_CC;
      real *FData_MomX = IntData_CC + MOMX*FSize3D_CC;
      real *FData_MomY = IntData_CC + MOMY*FSize3D_CC;
      real *FData_MomZ = IntData_CC + MOMZ*FSize3D_CC;
      real *FData_Engy = IntData_CC + ENGY*FSize3D_CC;
      real *FData_Enpy = IntData_CC + ENPY*FSize3D_CC;

      char dummy;    // we do not record the dual-energy status here

      for (int t=0; t<FSize3D_CC; t++)
      {
//       here we ALWAYS use the dual-energy variable to correct the total energy density
//       --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//           the runtime parameter DUAL_ENERGY_SWITCH here
         CPU_DualEnergyFix( FData_Dens[t], FData_MomX[t], FData_MomY[t], FData_MomZ[t], FData_Engy[t], FData_Enpy[t],
                            dummy, Gamma_m1, _Gamma_m1, (MinPres>=0.0), MinPres, UseEnpy2FixEngy );
      }
   } // if (  ( TVarCC & _TOTAL ) == _TOTAL  &&  DE_Consistency  )
#  endif // if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )

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
