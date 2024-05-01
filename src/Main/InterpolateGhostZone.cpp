#include "GAMER.h"

static int Table_01( const int FSide, const int CSide, const char dim, const int w01, const int w02,
                     const int w10, const int w11, const int w12, const int w20, const int w21 );
void SetTempIntPara( const int lv, const int Sg0, const double PrepTime, const double Time0, const double Time1,
                     bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );




//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolateGhostZone
// Description :  Fill up the ghost-zone values by spatial and temporal interpolation
//
// Note        :  1. Work for Prepare_PatchData()
//                2. Use the input parameter "TVarCC" and "TVarFC" to control the target variables
//                3. Use the input parameters "IntScheme_CC/FC" to control the interpolation schemes for the
//                   cell-/face-centered variables
//                4. Invoke Interpolate() for spatial interpolation
//                5. Data preparation order: FLU -> PASSIVE -> DERIVED --> POTE --> GRA_RHO
//                   ** DERIVED must be prepared immediately after FLU and PASSIVE so that both FLU, PASSIVE, and DERIVED
//                      can be prepared at the same time for the non-periodic BC. **
//                6. Use PrepTime to determine the physical time to prepare data
//                   --> Temporal interpolation/extrapolation will be conducted automatically if PrepTime
//                       is NOT equal to the time of data stored previously (e.g., FluSgTime[0/1])
//
// Parameter   :  lv                 : Target "coarse-grid" refinement level
//                PID                : Patch ID at level "lv" used for interpolation
//                IntData_CC/FC      : Arrays to store the cell-/face-centered interpolation results
//                IntData_CC_IntTime : Array to store the cell-centered interpolation result when temporal interpolation
//                                     is required
//                                     --> Only used for ELBDM with IntPhase
//                FSide              : Fine-patch sibling index (0~25) for determining the interpolation region
//                PrepTime           : Target physical time to prepare data
//                GhostSize          : Number of ghost zones
//                IntScheme_CC       : Interpolation scheme for the cell-centered variables
//                                     --> Supported schemes include
//                                         INT_MINMOD1D : MinMod-1D
//                                         INT_MINMOD3D : MinMod-3D
//                                         INT_VANLEER  : vanLeer
//                                         INT_CQUAD    : conservative quadratic
//                                         INT_QUAD     : quadratic
//                                         INT_CQUAR    : conservative quartic
//                                         INT_QUAR     : quartic
//                IntScheme_FC       : Interpolation scheme for the face-centered variables
//                                     --> Supported schemes include
//                                         INT_MINMOD1D : MinMod-1D
//                                         INT_VANLEER  : vanLeer
//                                         INT_CQUAD    : conservative quadratic
//                                         INT_CQUAR    : conservative quartic
//                NTSib              : Number of target sibling patches along different sibling directions
//                TSib               : Target sibling indices along different sibling directions
//                TVarCC             : Target cell-centered variables to be prepared
//                                     --> Supported variables in different models:
//                                         HYDRO        : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELX, _VELY, _VELZ, _PRES, _TEMP, _ENTR, _EINT
//                                                        [, _POTE] [, _MAGX_CC, _MAGY_CC, _MAGZ_CC, _MAGE_CC]
//                                         ELBDM_WAVE   : _DENS, _REAL, _IMAG [, _POTE]
//                                         ELBDM_HYBRID : _DENS, _PHAS [, _POTE]
//                                     --> _FLUID, _PASSIVE, _TOTAL, and _DERIVED apply to all models
//                NVarCC_Tot         : Total number of cell-centered variables to be prepared
//                NVarCC_Flu         : Number of cell-centered fluid variables to be prepared
//                                     --> Including passive scalars
//                TVarCCIdxList_Flu  : List recording the target cell-centered fluid and passive variable indices
//                                     ( = [0 ... NCOMP_TOTAL-1] )
//                NVarCC_Der         : Number of cell-centered derived variables to be prepared
//                TVarCCList_Der     : List recording the target cell-centered derived variables
//                TVarFC             : Target face-centered variables to be prepared
//                                     --> Supported variables in different models:
//                                         HYDRO with MHD : _MAG
//                                         ELBDM          : none
//                                     --> Note that for MHD it does NOT support preparing individual B component
//                                         (e.g., _MAGX, _MAGY, _MAGZ) since the divergenece-free interpolation
//                                         must work on all three components at once
//                NVarFC_Tot         : Total number of face-centered variables to be prepared
//                TVarFCIdxList      : List recording the target face-centered variable indices
//                                     ( = [0 ... NCOMP_MAG-1] )
//                IntPhase           : true --> Perform interpolation on rho/phase instead of real/imag parts in ELBDM
//                                     This parameter is useless for the ELBDM hybrid solver, which always interpolates rho/phase
//                FluBC              : Fluid boundary condition
//                PotBC              : Gravity boundary condition (not used currently)
//                BC_Face            : Priority of the B.C. along different boundary faces (z>y>x)
//                MinPres/Temp/Entr  : Minimum allowed pressure/temperature/entropy (<0.0 ==> off)
//                DE_Consistency     : Ensure the consistency between pressure, total energy density, and the
//                                     dual-energy variable when DUAL_ENERGY is on
//                FInterface         : B field on the coarse-fine interfaces for the divergence-preserving interpolation
//-------------------------------------------------------------------------------------------------------
void InterpolateGhostZone( const int lv, const int PID, real IntData_CC[], real IntData_FC[], real IntData_CC_IntTime[],
                           const int FSide, const double PrepTime, const int GhostSize,
                           const IntScheme_t IntScheme_CC, const IntScheme_t IntScheme_FC,
                           const int NTSib[], int *TSib[], const long TVarCC, const int NVarCC_Tot,
                           const int NVarCC_Flu, const int TVarCCIdxList_Flu[],
                           const int NVarCC_Der, const long TVarCCList_Der[],
                           const long TVarFC, const int NVarFC_Tot, const int TVarFCIdxList[],
                           const bool IntPhase, const OptFluBC_t FluBC[], const OptPotBC_t PotBC,
                           const int BC_Face[], const real MinPres, const real MinTemp, const real MinEntr,
                           const bool DE_Consistency, const real *FInterface[6] )
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

   if ( TVarFC != _NONE )
   {
#     ifdef MHD
      if ( TVarFC != _MAG  ||  NVarFC_Tot != NCOMP_MAG )
         Aux_Error( ERROR_INFO, "must work on all three magnetic components at once !!\n" );
#     else
      Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG


// adopt OPT__INT_PRIM and INT_REDUCE_MONO_COEFF when preparing all and only conserved variables (and magnetic field in MHD)
// --> still preserve conservation even with OPT__INT_PRIM because ghost zones do not affect conservation
#  if ( MODEL == HYDRO )
#  ifdef MHD
   const AllCons_t AllCons = ( TVarCC == _TOTAL  &&  TVarFC == _MAG ) ? ALL_CONS_YES : ALL_CONS_NO;
#  else
   const AllCons_t AllCons = ( TVarCC == _TOTAL                     ) ? ALL_CONS_YES : ALL_CONS_NO;
#  endif
   const bool      IntIter = ( AllCons == ALL_CONS_YES );

#  else // if ( MODEL == HYDRO )
   const AllCons_t AllCons = ALL_CONS_NO;
   const bool      IntIter = false;
#  endif // if ( MODEL == HYDRO ) ... else ...


// determine the target fields
#  if   ( MODEL == HYDRO )
   const bool PrepVx      = ( TVarCC & _VELX    ) ? true : false;
   const bool PrepVy      = ( TVarCC & _VELY    ) ? true : false;
   const bool PrepVz      = ( TVarCC & _VELZ    ) ? true : false;
   const bool PrepPres    = ( TVarCC & _PRES    ) ? true : false;
   const bool PrepTemp    = ( TVarCC & _TEMP    ) ? true : false;
   const bool PrepEntr    = ( TVarCC & _ENTR    ) ? true : false;
   const bool PrepEint    = ( TVarCC & _EINT    ) ? true : false;
#  ifdef MHD
   const bool PrepMagX_CC = ( TVarCC & _MAGX_CC ) ? true : false;
   const bool PrepMagY_CC = ( TVarCC & _MAGY_CC ) ? true : false;
   const bool PrepMagZ_CC = ( TVarCC & _MAGZ_CC ) ? true : false;
   const bool PrepMagE_CC = ( TVarCC & _MAGE_CC ) ? true : false;
   const bool PrepMagCC   = ( PrepMagX_CC || PrepMagY_CC || PrepMagZ_CC || PrepMagE_CC );
#  endif

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef GRAVITY
   const bool PrepPot  = ( TVarCC & _POTE ) ? true : false;
#  endif


// fluid variables for the EoS routines
#  if ( MODEL == HYDRO )
#  if ( EOS == EOS_GAMMA  ||  EOS == EOS_ISOTHERMAL )
   const int NFluForEoS = NCOMP_FLUID;    // don't need passsive scalars in EOS_GAMMA/EOS_ISOTHERMAL
#  else
   const int NFluForEoS = NCOMP_TOTAL;
#  endif
   real FluidForEoS[NFluForEoS];
#  endif


// set up parameters for the adopted interpolation scheme
   int NSide_CC_Useless, CGhost_CC, CSize_CC[3], FSize_CC[3], CSize3D_CC, FSize3D_CC;
   int NSide_FC_Useless, CGhost_FC, CSize_FC[3][3], FSize_FC[3][3], CSize3D_FC[3], FSize3D_FC[3];

   if ( NVarCC_Tot > 0 )   Int_Table( IntScheme_CC, NSide_CC_Useless, CGhost_CC );
   else                    CGhost_CC = 0;

   if ( NVarFC_Tot > 0 )   Int_Table( IntScheme_FC, NSide_FC_Useless, CGhost_FC );
   else                    CGhost_FC = 0;

   const double dh                 = amr->dh[lv];
   const int    GhostSize_Padded   = GhostSize + (GhostSize&1);         // # of interpolated cells must an even number
   const int    GhostSize_Padded_2 = GhostSize_Padded / 2;
   const int    CGrid_CC           = GhostSize_Padded_2 + 2*CGhost_CC;  // # of coarse cells required for interpolating
   const int    CGrid_FC_N         = GhostSize_Padded_2 + 1;            // cell-/face-centered variables
   const int    CGrid_FC_T         = GhostSize_Padded_2 + 2*CGhost_FC;  // (N/T = normal/transverse)
   const int    CGrid_CC_PID       = CGrid_CC   - CGhost_CC;            // # of overlapped cells between CGrid_CC/FC and PID
   const int    CGrid_FC_PID       = CGrid_FC_T - CGhost_FC;

   for (int d=0; d<3; d++)
   {
      CSize_CC[d] = TABLE_01( FSide, 'x'+d, CGrid_CC, PS1+2*CGhost_CC, CGrid_CC );
      FSize_CC[d] = TABLE_01( FSide, 'x'+d, GhostSize_Padded, PS2, GhostSize_Padded );

      for (int v=0; v<3; v++)
      {
         if ( v == d )
         {
            CSize_FC[v][d] = TABLE_01( FSide, 'x'+d, CGrid_FC_N, PS1+1, CGrid_FC_N );
            FSize_FC[v][d] = FSize_CC[d] + 1;
         }
         else
         {
            CSize_FC[v][d] = TABLE_01( FSide, 'x'+d, CGrid_FC_T, PS1+2*CGhost_FC, CGrid_FC_T );
            FSize_FC[v][d] = FSize_CC[d];
         }
      }
   }

   CSize3D_CC = CSize_CC[0]*CSize_CC[1]*CSize_CC[2];
   FSize3D_CC = FSize_CC[0]*FSize_CC[1]*FSize_CC[2];
   for (int v=0; v<3; v++)
   {
      CSize3D_FC[v] = CSize_FC[v][0]*CSize_FC[v][1]*CSize_FC[v][2];
      FSize3D_FC[v] = FSize_FC[v][0]*FSize_FC[v][1]*FSize_FC[v][2];
   }

// we assume that we only need ONE coarse-grid patch in each sibling direction
   if ( CGrid_CC_PID > PS1 )
      Aux_Error( ERROR_INFO, "CGrid_CC_PID (%d) > PATCH_SIZE (%d) !!\n", CGrid_CC_PID, PS1 );

   if ( CGrid_FC_PID > PS1 )
      Aux_Error( ERROR_INFO, "CGrid_FC_PID (%d) > PATCH_SIZE (%d) !!\n", CGrid_FC_PID, PS1 );


// coarse-grid data for interpolation (including the ghost zones on each side)
#  if ( MODEL == HYDRO  &&  defined MHD )
// IntIter requires the cell-centered B field
   const int NVarCC_Allocate = ( IntIter ) ? NVarCC_Tot+NCOMP_MAG : NVarCC_Tot;
#  else
   const int NVarCC_Allocate = NVarCC_Tot;
#  endif
   real *CData_CC_Ptr = NULL;
   real *CData_CC     = new real [ NVarCC_Allocate*CSize3D_CC ];
   real **CData_FC    = NULL;

// assuming NVarFC_Tot = either 0 or 3
   CData_FC = new real* [NVarFC_Tot];
   for (int v=0; v<NVarFC_Tot; v++)    CData_FC[v] = new real [ CSize3D_FC[v] ];


// temporal interpolation parameters
   bool FluIntTime;
   int  FluSg, FluSg_IntT;
   real FluWeighting, FluWeighting_IntT;

// fluid
   if ( NVarCC_Flu + NVarCC_Der != 0 )
   {
      const int Sg0 = amr->FluSg[lv];
      SetTempIntPara( lv, Sg0, PrepTime, amr->FluSgTime[lv][Sg0], amr->FluSgTime[lv][1-Sg0],
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

// check PrepPres, PrepTemp, PrepEntr, and PrepEint since they also require B field
   if ( NVarFC_Tot>0 || PrepMagCC || IntIter || PrepPres || PrepTemp || PrepEntr || PrepEint )
   {
      const int Sg0 = amr->MagSg[lv];
      SetTempIntPara( lv, Sg0, PrepTime, amr->MagSgTime[lv][Sg0], amr->MagSgTime[lv][1-Sg0],
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
      const int Sg0 = amr->PotSg[lv];
      SetTempIntPara( lv, Sg0, PrepTime, amr->PotSgTime[lv][Sg0], amr->PotSgTime[lv][1-Sg0],
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
   double xyz_flu[3];   // corner coordinates for the user-specified fluid B.C.

   for (int d=0; d<3; d++)
   {
      Loop1  [d] = TABLE_01( FSide, 'x'+d, CGrid_CC_PID, PS1, CGrid_CC_PID );
      Disp1  [d] = TABLE_01( FSide, 'x'+d, PS1-CGrid_CC_PID, 0, 0 );
      Disp2  [d] = TABLE_01( FSide, 'x'+d, 0, CGhost_CC, CGhost_CC );
      xyz_flu[d] = TABLE_01( FSide, 'x'+d, amr->patch[0][lv][PID]->EdgeL[d] + (0.5+PS1-CGrid_CC_PID)*dh,
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

//       temporal interpolation
//       --> for IntPhase, apply temporal interpolation to density/phase instead of real/imaginary parts for better accuracy
#        if ( MODEL == ELBDM )
#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
//       for fluid patches, we do not require the IntPhase flag and therefore only check whether FluIntTime is set
         if (   ( amr->use_wave_flag[lv] == true  && FluIntTime && !IntPhase )
             || ( amr->use_wave_flag[lv] == false && FluIntTime )  )
#        else
         if ( FluIntTime  &&  !IntPhase )
#        endif
#        else // #if ( MODEL == ELBDM )
         if ( FluIntTime )
#        endif // #if ( MODEL == ELBDM ) ... else ...
         CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                             + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][PID]->fluid[TVarCCIdx_Flu][k1][j1][i1];
         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   } // for (int v=0; v<NVarCC_Flu; v++)


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
   } // if ( PrepVx )

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
   } // if ( PrepVy )

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
   } // if ( PrepVz )

   if ( PrepPres )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

//###REVISE: support dual energy
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg );
#        else
         const real Emag = NULL_REAL;
#        endif
         CData_CC_Ptr[Idx] = Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                             FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                             (MinPres>=(real)0.0), MinPres, Emag,
                                             EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                             EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                             h_EoS_Table, NULL );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg_IntT );
#           else
            const real Emag = NULL_REAL;
#           endif
            CData_CC_Ptr[Idx] =
               FluWeighting     *CData_CC_Ptr[Idx]
             + FluWeighting_IntT*Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                 FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                 (MinPres>=(real)0.0), MinPres, Emag,
                                                 EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                 EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                 h_EoS_Table, NULL );
         }

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   } // if ( PrepPres )

   if ( PrepTemp )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

//###REVISE: support dual energy
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg );
#        else
         const real Emag = NULL_REAL;
#        endif
         CData_CC_Ptr[Idx] = Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                             FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                             (MinTemp>=(real)0.0), MinTemp, Emag,
                                             EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg_IntT );
#           else
            const real Emag = NULL_REAL;
#           endif
            CData_CC_Ptr[Idx] =
               FluWeighting     *CData_CC_Ptr[Idx]
             + FluWeighting_IntT*Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                 FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                 (MinTemp>=(real)0.0), MinTemp, Emag,
                                                 EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         }

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   } // if ( PrepTemp )

#  ifndef SRHD
   if ( PrepEntr )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg );
#        else
         const real Emag = NULL_REAL;
#        endif
         CData_CC_Ptr[Idx] = Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                             FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                             (MinEntr>=(real)0.0), MinEntr, Emag,
                                             EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                             h_EoS_Table );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg_IntT );
#           else
            const real Emag = NULL_REAL;
#           endif
            CData_CC_Ptr[Idx] =
               FluWeighting     *CData_CC_Ptr[Idx]
             + FluWeighting_IntT*Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                 FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                 (MinEntr>=(real)0.0), MinEntr, Emag,
                                                 EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                 h_EoS_Table );
         }

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   } // if ( PrepEntr )
#  endif // #ifndef SRHD

   if ( PrepEint )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k1][j1][i1];

//###REVISE: support dual energy
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg );
#        else
         const real Emag = NULL_REAL;
#        endif
         const bool CheckMinEint_No = false; // floor value is not supported for now
         CData_CC_Ptr[Idx] = Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                             FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                             CheckMinEint_No, NULL_REAL, Emag,
                                             EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

         if ( FluIntTime ) // temporal interpolation
         {
            for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k1][j1][i1];

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i1, j1, k1, MagSg_IntT );
#           else
            const real Emag = NULL_REAL;
#           endif
            CData_CC_Ptr[Idx] =
               FluWeighting     *CData_CC_Ptr[Idx]
             + FluWeighting_IntT*Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                 FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                 CheckMinEint_No, NULL_REAL, Emag,
                                                 EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         }

         Idx ++;
      }}}

      CData_CC_Ptr += CSize3D_CC;
   } // if ( PrepEint )

#  ifdef MHD
   if ( PrepMagCC || IntIter )
   {
      for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
      for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                          Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
      for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

         real B_CC[NCOMP_MAG];
         int  OffsetB = 0;

         MHD_GetCellCenteredBFieldInPatch( B_CC, lv, PID, i1, j1, k1, MagSg );

         if ( PrepMagX_CC || IntIter ) { CData_CC_Ptr[ Idx + OffsetB ] = B_CC[MAGX];  OffsetB += CSize3D_CC; }
         if ( PrepMagY_CC || IntIter ) { CData_CC_Ptr[ Idx + OffsetB ] = B_CC[MAGY];  OffsetB += CSize3D_CC; }
         if ( PrepMagZ_CC || IntIter ) { CData_CC_Ptr[ Idx + OffsetB ] = B_CC[MAGZ];  OffsetB += CSize3D_CC; }
         if ( PrepMagE_CC )            { CData_CC_Ptr[ Idx + OffsetB ] = (real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) ); }

         if ( FluIntTime ) // temporal interpolation
         {
            OffsetB = 0;

            MHD_GetCellCenteredBFieldInPatch( B_CC, lv, PID, i1, j1, k1, MagSg_IntT );

            if ( PrepMagX_CC || IntIter ) {
               CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                               + FluWeighting_IntT*B_CC[MAGX];
               OffsetB += CSize3D_CC;
            }

            if ( PrepMagY_CC || IntIter ) {
               CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                               + FluWeighting_IntT*B_CC[MAGY];
               OffsetB += CSize3D_CC;
            }

            if ( PrepMagZ_CC || IntIter ) {
               CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                               + FluWeighting_IntT*B_CC[MAGZ];
               OffsetB += CSize3D_CC;
            }

            if ( PrepMagE_CC ) {
               CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                               + FluWeighting_IntT*(real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) );
            }
         } // if ( FluIntTime )

         Idx ++;
      }}}

      if ( PrepMagX_CC || IntIter )    CData_CC_Ptr += CSize3D_CC;
      if ( PrepMagY_CC || IntIter )    CData_CC_Ptr += CSize3D_CC;
      if ( PrepMagZ_CC || IntIter )    CData_CC_Ptr += CSize3D_CC;
      if ( PrepMagE_CC )               CData_CC_Ptr += CSize3D_CC;
   } // if ( PrepMagCC )
#  endif // #ifdef MHD


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
   } // if ( PrepPot )
#  endif // #ifdef GRAVITY


// a4. face-centered variables (e.g., magnetic field)
   for (int v=0; v<NVarFC_Tot; v++)
   {
      const int TVarFCIdx = TVarFCIdxList[v];

#     ifdef MHD
//    set array indices
      const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                           ( TVarFCIdx == MAGY ) ? 1 :
                           ( TVarFCIdx == MAGZ ) ? 2 : -1;
#     ifdef GAMER_DEBUG
      if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#     endif

      int ijk_os[3], ijk_oe[3], size_i[3], disp_i[3], idx_i, idx_o, ii, ji, ki;  // s/e=start/end; i/o=in/out

      for (int d=0; d<3; d++)
      {
         if ( d == norm_dir )
         {
            ijk_os[d] = 0;
            ijk_oe[d] = CSize_FC[v][d];
            size_i[d] = PS1 + 1;
            disp_i[d] = TABLE_01( FSide, 'x'+d, PS1-GhostSize_Padded_2, 0, 0 );
         }

         else
         {
            ijk_os[d] = TABLE_01( FSide, 'x'+d, 0, CGhost_FC, CGhost_FC );
            ijk_oe[d] = ijk_os[d] + TABLE_01( FSide, 'x'+d, CGrid_FC_PID, PS1, CGrid_FC_PID );
            size_i[d] = PS1;
            disp_i[d] = TABLE_01( FSide, 'x'+d, PS1-CGrid_FC_PID, 0, 0 ) - ijk_os[d];
         }
      }


//    copy data
      for (int ko=ijk_os[2]; ko<ijk_oe[2]; ko++)   {  ki    = ko + disp_i[2];
      for (int jo=ijk_os[1]; jo<ijk_oe[1]; jo++)   {  ji    = jo + disp_i[1];
                                                      idx_i = IDX321( ijk_os[0]+disp_i[0], ji, ki, size_i[0], size_i[1] );
                                                      idx_o = IDX321( ijk_os[0],           jo, ko, CSize_FC[v][0], CSize_FC[v][1] );
      for (int io=ijk_os[0]; io<ijk_oe[0]; io++)   {

         CData_FC[v][idx_o] = amr->patch[MagSg][lv][PID]->magnetic[TVarFCIdx][idx_i];

         if ( MagIntTime ) // temporal interpolation
         CData_FC[v][idx_o] =   MagWeighting     *CData_FC[v][idx_o]
                              + MagWeighting_IntT*amr->patch[MagSg_IntT][lv][PID]->magnetic[TVarFCIdx][idx_i];

         idx_i ++;
         idx_o ++;
      }}}

#     else
      Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#     endif // #ifdef MHD ... else ...
   } // for (int v=0; v<NVarFC_Tot; v++)



// b. fill up the ghost zone of CData_CC[] and CData_FC[]
// ------------------------------------------------------------------------------------------------------------
   int Loop2[3], Disp3[3], Disp4[3], CSide, SibPID, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];

   for (int s=0; s<NTSib[FSide]; s++)
   {
      CSide  = TSib[FSide][s];
      SibPID = amr->patch[0][lv][PID]->sibling[CSide];

      for (int d=0; d<3; d++)
      {
         Loop2[d] = Table_01( FSide, CSide, 'x'+d, CGrid_CC_PID, CGhost_CC, CGhost_CC, PS1, CGhost_CC, CGhost_CC, CGrid_CC_PID );
         Disp3[d] = Table_01( FSide, CSide, 'x'+d, 0, CGrid_CC_PID, 0, CGhost_CC, CGhost_CC+PS1, 0, CGhost_CC );
         Disp4[d] = Table_01( FSide, CSide, 'x'+d, PS1-CGrid_CC_PID, 0, PS1-CGhost_CC, 0, 0, PS1-CGhost_CC, 0 );
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

//             temporal interpolation
//             --> for IntPhase, apply temporal interpolation to density/phase instead of real/imaginary parts for better accuracy
#              if ( MODEL == ELBDM )
#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               if (   ( amr->use_wave_flag[lv] == true  && FluIntTime && !IntPhase )
                   || ( amr->use_wave_flag[lv] == false && FluIntTime )  )
#              else
               if ( FluIntTime  &&  !IntPhase )
#              endif
#              else // #if ( MODEL == ELBDM )
               if ( FluIntTime )
#              endif // #if ( MODEL == ELBDM ) ... else ...
               CData_CC_Ptr[Idx] =   FluWeighting     *CData_CC_Ptr[Idx]
                                   + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][SibPID]->fluid[TVarCCIdx_Flu][k2][j2][i2];

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         } // for (int v=0; v<NVarCC_Flu; v++)


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
         } // if ( PrepVx )

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
         } // if ( PrepVy )

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
         } // if ( PrepVz )

         if ( PrepPres )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

//###REVISE: support dual energy
#              ifdef MHD
               const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg );
#              else
               const real Emag = NULL_REAL;
#              endif
               CData_CC_Ptr[Idx] = Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                   FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                   (MinPres>=(real)0.0), MinPres, Emag,
                                                   EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                   EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                   h_EoS_Table, NULL );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg_IntT );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  CData_CC_Ptr[Idx] =
                     FluWeighting     *CData_CC_Ptr[Idx]
                   + FluWeighting_IntT*Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                       FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                       (MinPres>=(real)0.0), MinPres, Emag,
                                                       EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                       EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                       h_EoS_Table, NULL );
               }

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         } // if ( PrepPres )

         if ( PrepTemp )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

//###REVISE: support dual energy
#              ifdef MHD
               const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg );
#              else
               const real Emag = NULL_REAL;
#              endif
               CData_CC_Ptr[Idx] = Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                   FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                   (MinTemp>=(real)0.0), MinTemp, Emag,
                                                   EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg_IntT );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  CData_CC_Ptr[Idx] =
                     FluWeighting     *CData_CC_Ptr[Idx]
                   + FluWeighting_IntT*Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                       FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                       (MinTemp>=(real)0.0), MinTemp, Emag,
                                                       EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
               }

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         } // if ( PrepTemp )

#        ifndef SRHD
         if ( PrepEntr )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

#              ifdef MHD
               const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg );
#              else
               const real Emag = NULL_REAL;
#              endif
               CData_CC_Ptr[Idx] = Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                   FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                   (MinEntr>=(real)0.0), MinEntr, Emag,
                                                   EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                   h_EoS_Table );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg_IntT );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  CData_CC_Ptr[Idx] =
                     FluWeighting     *CData_CC_Ptr[Idx]
                   + FluWeighting_IntT*Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                       FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                       (MinEntr>=(real)0.0), MinEntr, Emag,
                                                       EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                       h_EoS_Table );
               }

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         } // if ( PrepEntr )
#        endif // #ifndef SRHD

         if ( PrepEint )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][k2][j2][i2];

//###REVISE: support dual energy
#              ifdef MHD
               const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg );
#              else
               const real Emag = NULL_REAL;
#              endif
               const bool CheckMinEint_No = false; // floor value is not supported for now
               CData_CC_Ptr[Idx] = Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                   FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                   CheckMinEint_No, NULL_REAL, Emag,
                                                   EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

               if ( FluIntTime ) // temporal interpolation
               {
                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][k2][j2][i2];

#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, i2, j2, k2, MagSg_IntT );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  CData_CC_Ptr[Idx] =
                     FluWeighting     *CData_CC_Ptr[Idx]
                   + FluWeighting_IntT*Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                       FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                       CheckMinEint_No, NULL_REAL, Emag,
                                                       EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
               }

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         } // if ( PrepEint )

#        ifdef MHD
         if ( PrepMagCC || IntIter )
         {
            for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
            for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
            for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

               real B_CC[NCOMP_MAG];
               int  OffsetB = 0;

               MHD_GetCellCenteredBFieldInPatch( B_CC, lv, SibPID, i2, j2, k2, MagSg );

               if ( PrepMagX_CC || IntIter ) { CData_CC_Ptr[ Idx + OffsetB ] = B_CC[MAGX];  OffsetB += CSize3D_CC; }
               if ( PrepMagY_CC || IntIter ) { CData_CC_Ptr[ Idx + OffsetB ] = B_CC[MAGY];  OffsetB += CSize3D_CC; }
               if ( PrepMagZ_CC || IntIter ) { CData_CC_Ptr[ Idx + OffsetB ] = B_CC[MAGZ];  OffsetB += CSize3D_CC; }
               if ( PrepMagE_CC )            { CData_CC_Ptr[ Idx + OffsetB ] = (real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) ); }

               if ( FluIntTime ) // temporal interpolation
               {
                  OffsetB = 0;

                  MHD_GetCellCenteredBFieldInPatch( B_CC, lv, SibPID, i2, j2, k2, MagSg_IntT );

                  if ( PrepMagX_CC || IntIter ) {
                     CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                                     + FluWeighting_IntT*B_CC[MAGX];
                     OffsetB += CSize3D_CC;
                  }

                  if ( PrepMagY_CC || IntIter ) {
                     CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                                     + FluWeighting_IntT*B_CC[MAGY];
                     OffsetB += CSize3D_CC;
                  }

                  if ( PrepMagZ_CC || IntIter ) {
                     CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                                     + FluWeighting_IntT*B_CC[MAGZ];
                     OffsetB += CSize3D_CC;
                  }

                  if ( PrepMagE_CC ) {
                     CData_CC_Ptr[ Idx + OffsetB ] =   FluWeighting     *CData_CC_Ptr[ Idx + OffsetB ]
                                                     + FluWeighting_IntT*(real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) );
                  }
               } // if ( FluIntTime )

               Idx ++;
            }}}

            if ( PrepMagX_CC || IntIter )    CData_CC_Ptr += CSize3D_CC;
            if ( PrepMagY_CC || IntIter )    CData_CC_Ptr += CSize3D_CC;
            if ( PrepMagZ_CC || IntIter )    CData_CC_Ptr += CSize3D_CC;
            if ( PrepMagE_CC )               CData_CC_Ptr += CSize3D_CC;
         } // if ( PrepMagCC )
#        endif // #ifdef MHD


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
         } // if ( PrepPot )
#        endif // #ifdef GRAVITY


//       b1-4. face-centered variables (e.g., magnetic field)
         for (int v=0; v<NVarFC_Tot; v++)
         {
            const int TVarFCIdx = TVarFCIdxList[v];

#           ifdef MHD
//          set array indices
            const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                                 ( TVarFCIdx == MAGY ) ? 1 :
                                 ( TVarFCIdx == MAGZ ) ? 2 : -1;
#           ifdef GAMER_DEBUG
            if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#           endif


//          only need ghost zones along the two transverse directions
            if ( CSide >= 6  ||  CSide == norm_dir*2  ||  CSide == norm_dir*2+1 )   continue;


            int ijk_os[3], ijk_oe[3], size_i[3], disp_i[3], idx_i, idx_o, ii, ji, ki;  // s/e=start/end; i/o=in/out

            for (int d=0; d<3; d++)
            {
               if ( d == norm_dir )
               {
                  ijk_os[d] = 0;
                  ijk_oe[d] = CSize_FC[v][d];
                  size_i[d] = PS1 + 1;
                  disp_i[d] = TABLE_01( FSide, 'x'+d, PS1-GhostSize_Padded_2, 0, 0 );
               }

               else
               {
                  ijk_os[d] = Table_01( FSide, CSide, 'x'+d, 0, CGrid_FC_PID, 0, CGhost_FC, CGhost_FC+PS1, 0, CGhost_FC );
                  ijk_oe[d] = ijk_os[d] + Table_01( FSide, CSide, 'x'+d, CGrid_FC_PID, CGhost_FC,
                                                    CGhost_FC, PS1, CGhost_FC, CGhost_FC, CGrid_FC_PID );
                  size_i[d] = PS1;
                  disp_i[d] = Table_01( FSide, CSide, 'x'+d, PS1-CGrid_FC_PID, 0,
                                        PS1-CGhost_FC, 0, 0, PS1-CGhost_FC, 0 ) - ijk_os[d];
               }
            }


//          copy data
            for (int ko=ijk_os[2]; ko<ijk_oe[2]; ko++)   {  ki    = ko + disp_i[2];
            for (int jo=ijk_os[1]; jo<ijk_oe[1]; jo++)   {  ji    = jo + disp_i[1];
                                                            idx_i = IDX321( ijk_os[0]+disp_i[0], ji, ki, size_i[0], size_i[1] );
                                                            idx_o = IDX321( ijk_os[0],           jo, ko, CSize_FC[v][0], CSize_FC[v][1] );
            for (int io=ijk_os[0]; io<ijk_oe[0]; io++)   {

               CData_FC[v][idx_o] = amr->patch[MagSg][lv][SibPID]->magnetic[TVarFCIdx][idx_i];

               if ( MagIntTime ) // temporal interpolation
               CData_FC[v][idx_o] =   MagWeighting     *CData_FC[v][idx_o]
                                    + MagWeighting_IntT*amr->patch[MagSg_IntT][lv][SibPID]->magnetic[TVarFCIdx][idx_i];

               idx_i ++;
               idx_o ++;
            }}}

#           else
            Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#           endif // #ifdef MHD ... else ...
         } // for (int v=0; v<NVarFC_Tot; v++)
      } // if ( SibPID >= 0 )


//    b2. if the target sibling patch does not exist --> something is wrong !!
      else if ( SibPID == -1 )
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (Rank %d, Lv %d, PID %d, CSide %d) !!\n",
                    "SibPID", SibPID, MPI_Rank, lv, PID, CSide );


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
            Aux_Error( ERROR_INFO, "FluBC == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID %d, PID %d, CSide %d) !!\n",
                       BC_Sibling, BC_Face[BC_Sibling], SibPID, PID, CSide );
#        endif

//       b3-1. fluid B.C.
         if ( NVarCC_Flu + NVarCC_Der > 0 )
         {
            switch ( FluBC[ BC_Face[BC_Sibling] ] )
            {
#              if ( MODEL == HYDRO )
               case BC_FLU_OUTFLOW:
                  Hydro_BoundaryCondition_Outflow   ( CData_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu+NVarCC_Der, CGhost_CC,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End );
               break;

               case BC_FLU_REFLECTING:
                  Hydro_BoundaryCondition_Reflecting( CData_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu,            CGhost_CC,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End,
                                                      TVarCCIdxList_Flu, NVarCC_Der, TVarCCList_Der );
               break;

               case BC_FLU_DIODE:
                  Hydro_BoundaryCondition_Diode     ( CData_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu,            CGhost_CC,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End,
                                                      TVarCCIdxList_Flu, NVarCC_Der, TVarCCList_Der );
               break;
#              endif

               case BC_FLU_USER:
                  Flu_BoundaryCondition_User        ( CData_CC_Ptr,                      NVarCC_Flu,            CGhost_CC,
                                                      CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End,
                                                      TVarCCIdxList_Flu, PrepTime, dh, xyz_flu, TVarCC, lv );
               break;

               default:
                  Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
            } // switch ( FluBC[ BC_Face[BC_Sibling] ] )

            CData_CC_Ptr += NVarCC_Flu*CSize3D_CC;
         } // if ( NVarCC_Flu + NVarCC_Der > 0 )


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


//       b3-3. face-centered variables B.C. (e.g., magnetic field)
//             --> work on one component at a time since the array sizes of different components are different
         for (int v=0; v<NVarFC_Tot; v++)
         {
            const int TVarFCIdx = TVarFCIdxList[v];

//          work for MHD only
#           ifdef MHD
//          get the normal direction
            const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                                 ( TVarFCIdx == MAGY ) ? 1 :
                                 ( TVarFCIdx == MAGZ ) ? 2 : -1;
#           ifdef GAMER_DEBUG
            if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#           endif

//          only need ghost zones along the two transverse directions
            if ( CSide >= 6  ||  CSide == norm_dir*2  ||  CSide == norm_dir*2+1 )   continue;

//          set array indices --> correspond to the **cell-centered** array
            int    FC_BC_Idx_Start[3], FC_BC_Idx_End[3], FC_BC_Size[3];
            double xyz_mag[3];   // cell-centered corner coordinates for the user-specified magnetic field B.C.
            for (int d=0; d<3; d++)
            {
               if ( d == norm_dir )
               {
                  FC_BC_Idx_Start[d] = 0;
                  FC_BC_Idx_End  [d] = CSize_FC[v][d] - 2;
                  FC_BC_Size     [d] = CSize_FC[v][d] - 1;
                  xyz_mag        [d] = TABLE_01( FSide, 'x'+d, amr->patch[0][lv][PID]->EdgeL[d] + (0.5+PS1-GhostSize_Padded_2)*dh,
                                                               amr->patch[0][lv][PID]->EdgeL[d] + 0.5*dh,
                                                               amr->patch[0][lv][PID]->EdgeL[d] + 0.5*dh );
               }

               else
               {
                  FC_BC_Idx_Start[d] = Table_01( FSide, CSide, 'x'+d, 0, CGrid_FC_PID, 0, CGhost_FC, CGhost_FC+PS1, 0, CGhost_FC );
                  FC_BC_Idx_End  [d] = Table_01( FSide, CSide, 'x'+d, CGrid_FC_PID, CGhost_FC, CGhost_FC, PS1, CGhost_FC,
                                                 CGhost_FC, CGrid_FC_PID ) + FC_BC_Idx_Start[d] - 1;
                  FC_BC_Size     [d] = CSize_FC[v][d];
                  xyz_mag        [d] = TABLE_01( FSide, 'x'+d, amr->patch[0][lv][PID]->EdgeL[d] + (0.5+PS1-CGrid_FC_PID)*dh,
                                                               amr->patch[0][lv][PID]->EdgeL[d] + (0.5-CGhost_FC)*dh,
                                                               amr->patch[0][lv][PID]->EdgeL[d] + (0.5-CGhost_FC)*dh );
               }
            }

            real *MagDataPtr[NCOMP_MAG] = { NULL, NULL, NULL };
            MagDataPtr[TVarFCIdx] = CData_FC[v];

            switch ( FluBC[ BC_Face[BC_Sibling] ] )
            {
               case BC_FLU_OUTFLOW:
                  MHD_BoundaryCondition_Outflow   ( MagDataPtr, BC_Face[BC_Sibling], 1, CGhost_FC,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &TVarFCIdx );
               break;

               case BC_FLU_REFLECTING:
                  MHD_BoundaryCondition_Reflecting( MagDataPtr, BC_Face[BC_Sibling], 1, CGhost_FC,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &TVarFCIdx );
               break;

               case BC_FLU_DIODE:
                  MHD_BoundaryCondition_Diode     ( MagDataPtr, BC_Face[BC_Sibling], 1, CGhost_FC,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &TVarFCIdx );
               break;

               case BC_FLU_USER:
                  MHD_BoundaryCondition_User      ( MagDataPtr, BC_Face[BC_Sibling], 1,
                                                    FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                    &TVarFCIdx, PrepTime, dh, xyz_mag, lv );
               break;

               default:
                  Aux_Error( ERROR_INFO, "unsupported MHD B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
            } // switch ( FluBC[ BC_Face[BC_Sibling] ] )

#           else
            Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#           endif // #ifdef MHD ... else ...
         } // for (int v=0; v<NVarFC_Tot; v++)

      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )

      else
         Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, CSide %d) !!\n", SibPID, PID, CSide );

   } // for (int s=0; s<NTSib[FSide]; s++)



// c. interpolation : CData_CC[] --> IntData_CC[]
// ------------------------------------------------------------------------------------------------------------
   const bool PhaseUnwrapping_Yes   = true;
   const bool PhaseUnwrapping_No    = false;
   const bool Monotonicity_Yes      = true;
   const bool Monotonicity_No       = false;
   const bool IntOppSign0thOrder_No = false;

   int CStart_CC[3], CRange_CC[3], FStart_CC[3], NVarCC_SoFar;

   for (int d=0; d<3; d++)
   {
      CStart_CC[d] = CGhost_CC;
      CRange_CC[d] = CSize_CC[d] - 2*CGhost_CC;
      FStart_CC[d] = 0;
   }


// determine which variables require **monotonic** interpolation
   bool Monotonicity_CC[NVarCC_Flu];

   for (int v=0; v<NVarCC_Flu; v++)
   {
      TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

#     if   ( MODEL == HYDRO )
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
      /*
      if ( TVarCCIdx_Flu == DENS  ||  TVarCCIdx_Flu == ENGY  ||  TVarCCIdx_Flu >= NCOMP_FLUID )
         Monotonicity_CC[v] = Monotonicity_Yes;
      else
         Monotonicity_CC[v] = Monotonicity_No;
      */
         Monotonicity_CC[v] = Monotonicity_Yes;

#     elif ( MODEL == ELBDM )
//    apply monotonic interpolation to density and all passive scalars
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if (   ( TVarCCIdx_Flu != REAL  &&  TVarCCIdx_Flu != IMAG && amr->use_wave_flag[lv] == true )
          || ( TVarCCIdx_Flu != PHAS  &&  TVarCCIdx_Flu != STUB && amr->use_wave_flag[lv] == false )  )
#     else
      if ( TVarCCIdx_Flu != REAL  &&  TVarCCIdx_Flu != IMAG )
#     endif
         Monotonicity_CC[v] = Monotonicity_Yes;
      else
         Monotonicity_CC[v] = Monotonicity_No;

#     else
#     error : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION IN THIS NEW MODEL ??
#     endif // MODEL
   } // for (int v=0; v<NVarCC_Flu; v++)


// interpolation
// c1. interpolation on face-centered variables
//     --> do it first since we need the cell-centered B field for IntIter
   if ( NVarFC_Tot > 0 )
   {
#     ifdef MHD
//    c1-1. set the array indices
      real *FData_FC[3] = { IntData_FC,
                            IntData_FC + FSize3D_FC[0],
                            IntData_FC + FSize3D_FC[0] + FSize3D_FC[1] };
      int CStart_FC[3][3], CRange_FC[3], FStart_FC[3][3];

      for (int d=0; d<3; d++)
      {
         CRange_FC[d] = CRange_CC[d];

         for (int v=0; v<3; v++)
         {
            CStart_FC[v][d] = ( v == d ) ? 0 : CGhost_FC;
            FStart_FC[v][d] = 0;
         }
      }

//    c1-2. divergence-perserving interpolation
      MHD_InterpolateBField( (const real**)CData_FC, CSize_FC, CStart_FC, CRange_FC,
                             FData_FC, FSize_FC, FStart_FC, FInterface,
                             IntScheme_FC, Monotonicity_Yes );

#     else
      Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#     endif
   } // if ( NVarFC_Tot > 0 )


// c2. interpolation on phase in ELBDM
#  if ( MODEL == ELBDM )

   real *CData_Real = NULL;
   real *CData_Imag = NULL;
   real *CData_Dens = NULL;
   real *CData_Phas = NULL;

   real *FData_Real = NULL;
   real *FData_Imag = NULL;
   real *FData_Dens = NULL;
   real *FData_Phas = NULL;

// parameter IntPhase in hybrid scheme is only relevant where wave scheme is used
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( IntPhase  &&  amr->use_wave_flag[lv] == true )
#  else
   if ( IntPhase )
#  endif
   {
//    determine the array indices
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

//    if we are not preparing the density field:
//    store density in the REAL component and phase in the IMAG component
      if ( DensIdx == -1 )
      {
         CData_Real = CData_CC   + RealIdx*CSize3D_CC;
         CData_Imag = CData_CC   + ImagIdx*CSize3D_CC;
         CData_Dens = CData_Real;
         CData_Phas = CData_Imag;

         FData_Real = IntData_CC + RealIdx*FSize3D_CC;
         FData_Imag = IntData_CC + ImagIdx*FSize3D_CC;
         FData_Dens = FData_Real;
         FData_Phas = FData_Imag;
//    otherwise store density in the DENS component and phase in the REAL component
//    this ensure the two arrays are consecutive in memory which is a necessary requirement for INT_SPECTRAL
      } else {
         CData_Real = CData_CC   + RealIdx*CSize3D_CC;
         CData_Imag = CData_CC   + ImagIdx*CSize3D_CC;
         CData_Dens = CData_CC   + DensIdx*CSize3D_CC;
         CData_Phas = CData_Real;

         FData_Real = IntData_CC + RealIdx*FSize3D_CC;
         FData_Imag = IntData_CC + ImagIdx*FSize3D_CC;
         FData_Dens = IntData_CC + DensIdx*FSize3D_CC;
         FData_Phas = FData_Real;
      } // if ( DensIdx == -1 ) ... else ...

//    get the density and wrapped phase
      for (int t=0; t<CSize3D_CC; t++)
      {
         const real Re = CData_Real[t];
         const real Im = CData_Imag[t];

//###ISSUE: atan2() sometimes returns NaN when both inputs are zero, not sure why ...
//          --> using SATAN2() in Macro.h seems to provide a temporary fix (but needs to be checked further)
//          --> in principle we don't need this fix since atan2(y,x) should be able to handle
//              y=+0/-0, x=+0/-0 (https://en.cppreference.com/w/c/numeric/math/atan2)
         if ( Re == (real)0.0  &&  Im == (real)0.0 )  CData_Phas[t] = (real)0.0;
         else                                         CData_Phas[t] = SATAN2( Im, Re );

         if ( DensIdx == -1 ) // only need to recalculate density if it's not prepared already
         CData_Dens[t] = Re*Re + Im*Im;
      }

      if ( IntScheme_CC == INT_SPECTRAL ) {
//    interpolate density & phase
//    INT_SPECTRAL with PhaseUnwrapping_Yes assumes that the density and phase fields are stored consecutively in memory
      const bool Monotonicity_Spec[2] = { true, false };
      Interpolate( CData_CC, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC, FSize_CC, FStart_CC,
                   2, IntScheme_CC, PhaseUnwrapping_Yes, Monotonicity_Spec, IntOppSign0thOrder_No,
                   ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      } else {
//    interpolate density
      Interpolate( CData_Dens, CSize_CC, CStart_CC, CRange_CC, FData_Dens, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes, IntOppSign0thOrder_No,
                   ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );

//    interpolate phase
      Interpolate( CData_Phas, CSize_CC, CStart_CC, CRange_CC, FData_Phas, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_Yes, &Monotonicity_No, IntOppSign0thOrder_No,
                   ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );

      }

//    temporal interpolation
//    --> apply it to density/phase instead of real/imaginary parts for better accuracy
      if ( FluIntTime )
      {
         const int TVarCCIdxList_Flu_IntTime[2] = { REAL, IMAG };
         const int NVarCC_Flu_IntTime           = 2;

//       a. fill up the central region of CData_CC[] with the data at FluSg_IntT
//       ------------------------------------------------------------------------------------------------------------
         CData_CC_Ptr = CData_CC;

         for (int v=0; v<NVarCC_Flu_IntTime; v++)
         {
            TVarCCIdx_Flu = TVarCCIdxList_Flu_IntTime[v];

            for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
            for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                                Idx = IDX321( Disp2[0], j2, k2, CSize_CC[0], CSize_CC[1] );
            for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

               CData_CC_Ptr[Idx] = amr->patch[FluSg_IntT][lv][PID]->fluid[TVarCCIdx_Flu][k1][j1][i1];

               Idx ++;
            }}}

            CData_CC_Ptr += CSize3D_CC;
         }


//       b. fill up the ghost zone of CData_CC[] with the data at FluSg_IntT
//       ------------------------------------------------------------------------------------------------------------
         for (int s=0; s<NTSib[FSide]; s++)
         {
            CSide  = TSib[FSide][s];
            SibPID = amr->patch[0][lv][PID]->sibling[CSide];

            for (int d=0; d<3; d++)
            {
               Loop2[d] = Table_01( FSide, CSide, 'x'+d, CGrid_CC_PID, CGhost_CC, CGhost_CC, PS1, CGhost_CC, CGhost_CC, CGrid_CC_PID );
               Disp3[d] = Table_01( FSide, CSide, 'x'+d, 0, CGrid_CC_PID, 0, CGhost_CC, CGhost_CC+PS1, 0, CGhost_CC );
               Disp4[d] = Table_01( FSide, CSide, 'x'+d, PS1-CGrid_CC_PID, 0, PS1-CGhost_CC, 0, 0, PS1-CGhost_CC, 0 );
            }

//          b1. if the target sibling patch exists --> just copy data from the nearby patch at the same level
            if ( SibPID >= 0 )
            {
               CData_CC_Ptr = CData_CC;

               for (int v=0; v<NVarCC_Flu_IntTime; v++)
               {
                  TVarCCIdx_Flu = TVarCCIdxList_Flu_IntTime[v];

                  for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
                  for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                      Idx = IDX321( Disp3[0], j1, k1, CSize_CC[0], CSize_CC[1] );
                  for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

                     CData_CC_Ptr[Idx] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[TVarCCIdx_Flu][k2][j2][i2];

                     Idx ++;
                  }}}

                  CData_CC_Ptr += CSize3D_CC;
               }
            } // if ( SibPID >= 0 )

//          b2. if the target sibling patch does not exist --> something is wrong !!
            else if ( SibPID == -1 )
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (Rank %d, Lv %d, PID %d, CSide %d) !!\n",
                          "SibPID", SibPID, MPI_Rank, lv, PID, CSide );

//          b3. if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
            {
               CData_CC_Ptr = CData_CC;

               for (int d=0; d<3; d++)
               {
                  BC_Idx_Start[d] = Disp3[d];
                  BC_Idx_End  [d] = Loop2[d] + BC_Idx_Start[d] - 1;
               }

               BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

               switch ( FluBC[ BC_Face[BC_Sibling] ] )
               {
                  case BC_FLU_USER:
                     Flu_BoundaryCondition_User( CData_CC_Ptr, NVarCC_Flu_IntTime, CGhost_CC,
                                                 CSize_CC[0], CSize_CC[1], CSize_CC[2], BC_Idx_Start, BC_Idx_End,
                                                 TVarCCIdxList_Flu_IntTime, PrepTime, dh, xyz_flu, _REAL|_IMAG, lv );
                  break;

                  default:
                     Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
               } // switch ( FluBC[ BC_Face[BC_Sibling] ] )
            } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )

            else
               Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, CSide %d) !!\n", SibPID, PID, CSide );

         } // for (int s=0; s<NTSib[FSide]; s++)


//       get the density and wrapped phase at FluSg_IntT
#        ifdef GAMER_DEBUG
         if ( IntData_CC_IntTime == NULL )   Aux_Error( ERROR_INFO, "IntData_CC_IntTime == NULL !!\n" );
#        endif

         real *CData_Real_IntTime = CData_CC           + 0*CSize3D_CC;
         real *CData_Imag_IntTime = CData_CC           + 1*CSize3D_CC;
         real *CData_Dens_IntTime = CData_Real_IntTime;
         real *CData_Phas_IntTime = CData_Imag_IntTime;

         real *FData_Dens_IntTime = IntData_CC_IntTime + 0*FSize3D_CC;
         real *FData_Phas_IntTime = IntData_CC_IntTime + 1*FSize3D_CC;

         for (int t=0; t<CSize3D_CC; t++)
         {
            const real Re = CData_Real_IntTime[t];
            const real Im = CData_Imag_IntTime[t];

//###ISSUE: atan2() sometimes returns NaN when both inputs are zero, not sure why ...
//          --> using SATAN2() in Macro.h seems to provide a temporary fix (but needs to be checked further)
//          --> in principle we don't need this fix since atan2(y,x) should be able to handle
//              y=+0/-0, x=+0/-0 (https://en.cppreference.com/w/c/numeric/math/atan2)
            if ( Re == (real)0.0  &&  Im == (real)0.0 )  CData_Phas_IntTime[t] = (real)0.0;
            else                                         CData_Phas_IntTime[t] = SATAN2( Im, Re );

            CData_Dens_IntTime[t] = Re*Re + Im*Im;
         }

         if ( IntScheme_CC == INT_SPECTRAL ) {
//       interpolate density & phase
//       INT_SPECTRAL with PhaseUnwrapping_Yes assumes that the density and phase fields are stored consecutively in memory
         const bool Monotonicity_Spec[2] = { true, false };
         Interpolate( CData_CC, CSize_CC, CStart_CC, CRange_CC,
                      IntData_CC_IntTime, FSize_CC, FStart_CC,
                      2, IntScheme_CC, PhaseUnwrapping_Yes, Monotonicity_Spec, IntOppSign0thOrder_No,
                      ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
         } else {
//       interpolate density
         Interpolate( CData_Dens_IntTime, CSize_CC, CStart_CC, CRange_CC,
                      FData_Dens_IntTime, FSize_CC, FStart_CC,
                      1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes, IntOppSign0thOrder_No,
                      ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );

//       interpolate phase
         Interpolate( CData_Phas_IntTime, CSize_CC, CStart_CC, CRange_CC,
                      FData_Phas_IntTime, FSize_CC, FStart_CC,
                      1, IntScheme_CC, PhaseUnwrapping_Yes, &Monotonicity_No, IntOppSign0thOrder_No,
                      ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );

         }

//       temporal interpolation
         for (int t=0; t<FSize3D_CC; t++)
         {
//          must unwrap phase before interpolating it
            FData_Phas_IntTime[t] = ELBDM_UnwrapPhase( FData_Phas[t], FData_Phas_IntTime[t] );

            FData_Dens[t] = FluWeighting     *FData_Dens        [t]
                          + FluWeighting_IntT*FData_Dens_IntTime[t];
            FData_Phas[t] = FluWeighting     *FData_Phas        [t]
                          + FluWeighting_IntT*FData_Phas_IntTime[t];
         }
      } // FluIntTime


//    density and phase --> real and imaginary parts
      real Dens, Phase, Amp;

      for (int t=0; t<FSize3D_CC; t++)
      {
         Dens  = FData_Dens[t];
         Phase = FData_Phas[t];

//       be careful about the negative density introduced from the round-off errors
//       --> note that we check minimum density in the end of Prepare_PatchData()
         if ( Dens < (real)0.0 )
         {
            FData_Dens[t] = (real)0.0;
            Dens          = (real)0.0;
         }

         Amp           = SQRT( Dens );
         FData_Real[t] = Amp*COS( Phase );
         FData_Imag[t] = Amp*SIN( Phase );
      }
   } // if ( IntPhase )  ||  if ( IntPhase && amr->use_wave_flag[lv] == true ) in hybrid scheme

// c3. interpolation on original variables
   else
#  endif // if ( MODEL == ELBDM )
   {
//    c3-1. prepare the fine-grid, cell-centered B field for IntIter
      real *CMag_CC_IntIter              = NULL;
      real (*FMag_CC_IntIter)[NCOMP_MAG] = NULL;

#     ifdef MHD
      if ( IntIter )
      {
         CMag_CC_IntIter = CData_CC + NCOMP_TOTAL*CSize3D_CC;

         const real *FData_FC[3] = { IntData_FC,
                                     IntData_FC + FSize3D_FC[0],
                                     IntData_FC + FSize3D_FC[0] + FSize3D_FC[1] };

         FMag_CC_IntIter = new real [FSize3D_CC][NCOMP_MAG];

         for (int k=0; k<FSize_CC[2]; k++)
         for (int j=0; j<FSize_CC[1]; j++)
         for (int i=0; i<FSize_CC[0]; i++)
         {
            const int t = IDX321( i, j, k, FSize_CC[0], FSize_CC[1] );

            MHD_GetCellCenteredBField( FMag_CC_IntIter[t], FData_FC[MAGX], FData_FC[MAGY], FData_FC[MAGZ],
                                       FSize_CC[0], FSize_CC[1], FSize_CC[2], i, j, k );
         }
      }
#     endif // MHD

#     if ( MODEL != HYDRO )
      const bool OPT__INT_PRIM = false;
#     endif
      Interpolate( CData_CC, CSize_CC, CStart_CC, CRange_CC, IntData_CC, FSize_CC, FStart_CC, NVarCC_Flu,
                   IntScheme_CC, PhaseUnwrapping_No, Monotonicity_CC, INT_OPP_SIGN_0TH_ORDER,
                   AllCons,
                   (IntIter && OPT__INT_PRIM)?INT_PRIM_YES:INT_PRIM_NO,
                   (IntIter                 )?INT_REDUCE_MONO_COEFF:INT_FIX_MONO_COEFF,
                   CMag_CC_IntIter, FMag_CC_IntIter );

      delete [] FMag_CC_IntIter;
   } // if ( IntPhase )  ||  if ( IntPhase && amr->use_wave_flag[lv] == true ) in hybrid scheme ... else ...

   NVarCC_SoFar = NVarCC_Flu;


// c4. interpolation on derived variables
#  if   ( MODEL == HYDRO )
// we now apply monotonic interpolation to ALL fluid variables
   if ( PrepVx )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }

   if ( PrepVy )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }

   if ( PrepVz )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }

   if ( PrepPres )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }

   if ( PrepTemp )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }

#  ifndef SRHD
   if ( PrepEntr )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }
#  endif

   if ( PrepEint )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }

#  ifdef MHD
   if ( PrepMagCC )
   {
      bool Monotonicity_Mag[4] = { true, true, true, true };
      int NMag = 0;
      if ( PrepMagX_CC )   NMag ++;
      if ( PrepMagY_CC )   NMag ++;
      if ( PrepMagZ_CC )   NMag ++;
      if ( PrepMagE_CC )   NMag ++;

      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   NMag, IntScheme_CC, PhaseUnwrapping_No, Monotonicity_Mag,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar += NMag;
   }
#  endif // #ifdef MHD

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif // MODEL


// c5. interpolation on potential
#  ifdef GRAVITY
   if ( PrepPot )
   {
      Interpolate( CData_CC+CSize3D_CC*NVarCC_SoFar, CSize_CC, CStart_CC, CRange_CC,
                   IntData_CC+FSize3D_CC*NVarCC_SoFar, FSize_CC, FStart_CC,
                   1, IntScheme_CC, PhaseUnwrapping_No, &Monotonicity_No,
                   IntOppSign0thOrder_No, ALL_CONS_NO, INT_PRIM_NO, INT_FIX_MONO_COEFF, NULL, NULL );
      NVarCC_SoFar ++;
   }
#  endif


// free memory
   delete [] CData_CC;
   for (int v=0; v<NVarFC_Tot; v++)    delete [] CData_FC[v];
   delete [] CData_FC;



// d. ensure the consistency between pressure, total energy density, and the dual-energy variable
//    when DUAL_ENERGY is on
//    --> we don't have to check the minimum pressure here when DUAL_ENERGY is off
//        --> it's checked in Prepare_PatchData()
// ------------------------------------------------------------------------------------------------------------
#  if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
// apply this correction only when preparing all fluid variables (and magnetic field)
#  ifdef MHD
   if (  DE_Consistency  &&  ( TVarCC & _TOTAL ) == _TOTAL  &&  TVarFC == _MAG )
#  else
   if (  DE_Consistency  &&  ( TVarCC & _TOTAL ) == _TOTAL )
#  endif
   {
#     ifdef MHD
      const real *FData_FC[3]    = { IntData_FC,
                                     IntData_FC + FSize3D_FC[0],
                                     IntData_FC + FSize3D_FC[0] + FSize3D_FC[1] };
      const int  size_ij         = FSize_CC[0]*FSize_CC[1];
#     endif
      const real UseDual2FixEngy = HUGE_NUMBER;

//    assuming that the order of variables stored in IntData_CC[] is the same as patch->fluid[]
      real *FData_Dens = IntData_CC + DENS*FSize3D_CC;
      real *FData_MomX = IntData_CC + MOMX*FSize3D_CC;
      real *FData_MomY = IntData_CC + MOMY*FSize3D_CC;
      real *FData_MomZ = IntData_CC + MOMZ*FSize3D_CC;
      real *FData_Engy = IntData_CC + ENGY*FSize3D_CC;
      real *FData_Dual = IntData_CC + DUAL*FSize3D_CC;

      char dummy;    // we do not record the dual-energy status here

      for (int t=0; t<FSize3D_CC; t++)
      {
//       compute the cell-centered magnetic energy
#        ifdef MHD
         const int  i     = t % FSize_CC[0];
         const int  j     = t % size_ij / FSize_CC[0];
         const int  k     = t / size_ij;
         const real Emag = MHD_GetCellCenteredBEnergy( FData_FC[MAGX], FData_FC[MAGY], FData_FC[MAGZ],
                                                       FSize_CC[0], FSize_CC[1], FSize_CC[2], i, j, k );
#        else
         const real Emag = NULL_REAL;
#        endif

//       here we ALWAYS use the dual-energy variable to correct the total energy density
//       --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//           the runtime parameter DUAL_ENERGY_SWITCH here
         Hydro_DualEnergyFix( FData_Dens[t], FData_MomX[t], FData_MomY[t], FData_MomZ[t], FData_Engy[t], FData_Dual[t],
                              dummy, EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2], (MinPres>=(real)0.0), MinPres,
                              UseDual2FixEngy, Emag );
      }
   } // if (  DE_Consistency  &&  ( TVarCC & _TOTAL ) == _TOTAL  &&  TVarFC == _MAG )
#  endif // if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )

} // FUNCTION : InterpolateGhostZone



// ============
// |  Tables  |
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  Return the loop size and displacement required by InterpolateGhostZone()
//
// Parameter   :  FSide       : Fine-patch sibling index (0~25)
//                CSide       : Coarse-patch sibling index (0~25)
//                dim         : Target spatial direction (x/y/z)
//                w01 ... w21 : Returned values
//-------------------------------------------------------------------------------------------------------
int Table_01( const int FSide, const int CSide, const char dim, const int w01, const int w02,
              const int w10, const int w11, const int w12, const int w20, const int w21 )
{

   switch ( dim )
   {
      case 'x':
      {
         switch ( FSide )
         {
            case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
            {
               switch ( CSide )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "FSide", FSide, "CSide", CSide );

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w01;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     return w02;
               }
            }

            case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
            {
               switch ( CSide )
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
               switch ( CSide )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     return w20;

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w21;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "FSide", FSide, "CSide", CSide );
               }
            }

         } // switch ( FSide )
      } // case 'x':


      case 'y':
      {
         switch ( FSide )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
            {
               switch ( CSide )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "FSide", FSide, "CSide", CSide );

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w01;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     return w02;
               }
            }

            case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
            {
               switch ( CSide )
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
               switch ( CSide )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     return w20;

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w21;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "FSide", FSide, "CSide", CSide );
               }
            }

         } // switch ( FSide )
      } // case 'y':


      case 'z':
      {
         switch ( FSide )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
            {
               switch ( CSide )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "FSide", FSide, "CSide", CSide );

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w01;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     return w02;
               }
            }

            case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
            {
               switch ( CSide )
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
               switch ( CSide )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     return w20;

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w21;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "FSide", FSide, "CSide", CSide );
               }
            }

         } // switch ( FSide )
      } // case 'z':


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c !!\n", "dim", dim );
         exit(1);

   } // switch ( dim )

   return NULL_INT;

} // FUNCTION : Table_01
