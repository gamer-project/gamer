#include "GAMER.h"

void InterpolateGhostZone( const int lv, const int PID, real IntData_CC[], real IntData_FC[], real IntData_CC_IntTime[],
                           const int FSide, const double PrepTime, const int GhostSize,
                           const IntScheme_t IntScheme_CC, const IntScheme_t IntScheme_FC,
                           const int NTSib[], int *TSib[], const long TVarCC, const int NVarCC_Tot,
                           const int NVarCC_Flu, const int TVarCCIdxList_Flu[],
                           const int NVarCC_Der, const long TVarCCList_Der[],
                           const long TVarFC, const int NVarFC_Tot, const int TVarFCIdxList[],
                           const bool IntPhase, const OptFluBC_t FluBC[], const OptPotBC_t PotBC,
                           const int BC_Face[], const real MinPres, const real MinTemp, const real MinEntr,
                           const bool DE_Consistency, const real *FInterface[6] );
static void SetTargetSibling( int NTSib[], int *TSib[] );
static int  Table_01( const int SibID, const char dim, const int Count, const int GhostSize );
static int  Table_02( const int lv, const int PID, const int Side );
static long Table_03( const int lv, const long GID, const int Side, LB_GlobalPatch* Tree );
void SetTempIntPara( const int lv, const int Sg0, const double PrepTime, const double Time0, const double Time1,
                     bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );
#ifdef MHD
static void MHD_SetFInterface( real *FInt_Data, real *FInt_Ptr[6], const real *Data1PG_FC, const int lv, const int PID0,
                               const int Side, const int GhostSize, const int MagSg, const int MagSg_IntT,
                               const bool MagIntTime, const real MagWeighting, const real MagWeighting_IntT );
#endif

// flag for checking whether Prepare_PatchData_InitParticleDensityArray() has been called
// before preparing either _PAR_DENS or _TOTAL_DENS
#ifdef PARTICLE
bool ParDensArray_Initialized = false;
#endif


// check the divergence-free B field (for debug)
#ifdef MHD
//#  define MHD_CHECK_DIV_B

# ifdef MHD_CHECK_DIV_B
#  ifdef FLOAT8
#  define DIV_B_TOLERANCE  1.0e-12
#  else
#  define DIV_B_TOLERANCE  1.0e-5f
#  endif

static void MHD_CheckDivB( const real *Data1PG_FC, const int GhostSize, const real Tolerance,
                           const int lv, const double PrepTime );
# endif // #ifdef MHD_CHECK_DIV_B
#endif // MHD




//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_PatchData
// Description :  Prepare the uniform data including ghost zones for the target patches or patch groups
//
// Note        :  1. Use the input parameters "TVarCC" and "TVarFC" to control the target cell-centered and
//                   face-centered variables, respectively
//                   --> TVarCC and TVarFC can be any combination of the symbolic constants defined in "Macro.h"
//                       (e.g., "TVarCC = _DENS", "TVarCC = _MOMX|ENGY", "TVarCC = _TOTAL", "TVarFC = _MAGX|_MAGZ")
//                2. If "GhostSize != 0" --> use InterpolateGhostZone() to fill out the
//                   ghost-zone values by spatial interpolation if the corresponding sibling patches do
//                   NOT exist
//                4. Use PrepTime to determine the physical time to prepare data
//                   --> Temporal interpolation/extrapolation will be conducted automatically if PrepTime
//                       is NOT equal to the time of data stored previously (e.g., FluSgTime[0/1])
//                4. Use "patch group" as the preparation unit
//                   --> The data of all patches within the same patch group will be prepared
//                5. It is assumed that both _FLUID and _PASSIVE are stored in the same Sg
//                6. Data are prepared and stored in the order:
//                   _FLUID (where _DENS may be replaced by _TOTAL_DENS) -> _PASSIVE -> _DERIVED --> _POTE --> _PAR_DENS
//                   ** DERIVED must be prepared immediately after FLU and PASSIVE so that both FLU, PASSIVE, and DERIVED
//                      can be prepared at the same time for the non-periodic BC. **
//                7. For _PAR_DENS and _TOTAL_DENS (for PARTICLE only), the rho_ext[] arrays of patches at Lv=lv must be
//                   pre-computed by calling Prepare_PatchData_InitParticleDensityArray() to store the partice mass density
//                   --> amr->patch[0][lv][PID]->rho_ext[]
//                   --> These arrays must be deallocated manually by calling Prepare_PatchData_FreeParticleDensityArray()
//                   --> Before calling this function, one must call
//                       (1) Par_CollectParticle2OneLevel() --> to collect particles from higher levels and from other MPI ranks
//                       (2) Prepare_PatchData_InitParticleDensityArray() --> to initialize all rho_ext[] arrays
//                   --> After calling this function, one must call the following two functions to free memory
//                       (1) Par_CollectParticle2OneLevel_FreeMemory()
//                       (2) Prepare_PatchData_FreeParticleDensityArray()
//                8. Patches stored in PID0_List must be real patches (cannot NOT be buffer patches)
//                9. Option "MHD_CHECK_DIV_B" on the top of this file can be used to check the divergence-free constraint
//                   on the prepared B field
//                   --> Note that when temporal interpolation is required in InterpolateGhostZone(), it will break div(B)=0 for
//                       the interpolated fine-grid B field just outside the central patch group
//                       --> It's because the temporal interpolated B field is in general not equal to the original fine-grid B
//                           field on the coarse-fine interfaces of the central patch group
//                       --> It's OK for the MHD solver since it will still guarantee that the updated B field within the patch group
//                           is divergence free
//
// Parameter   :  lv             : Target refinement level
//                PrepTime       : Target physical time to prepare data
//                OutputCC       : Returned array to store the prepared cell-centered data
//                OutputFC       : Returned array to store the prepared face-centered data
//                GhostSize      : Number of ghost zones to be prepared
//                NPG            : Number of patch groups prepared at a time
//                PID0_List      : List recording the patch indices with LocalID==0 to be prepared
//                TVarCC         : Target cell-centered variables to be prepared
//                                 --> Supported variables in different models:
//                                     HYDRO        : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELX, _VELY, _VELZ, _PRES, _TEMP, _ENTR, _EINT
//                                                    [, _POTE] [, _MAGX_CC, _MAGY_CC, _MAGZ_CC, _MAGE_CC]
//                                     ELBDM_WAVE   : _DENS, _REAL, _IMAG [, _POTE]
//                                     ELBDM_HYBRID : _DENS, _PHAS [, _POTE]
//                                 --> _FLUID, _PASSIVE, _TOTAL, and _DERIVED apply to all models
//                TVarFC         : Target face-centered variables to be prepared
//                                 --> Supported variables in different models:
//                                     HYDRO with MHD : _MAGX, _MAGY, _MAGZ, _MAG
//                                     ELBDM          : none
//                IntScheme_CC   : Interpolation scheme for the cell-centered variables
//                                 --> Supported schemes include
//                                     INT_MINMOD1D : MinMod-1D
//                                     INT_MINMOD3D : MinMod-3D
//                                     INT_VANLEER  : vanLeer
//                                     INT_CQUAD    : conservative quadratic
//                                     INT_QUAD     : quadratic
//                                     INT_CQUAR    : conservative quartic
//                                     INT_QUAR     : quartic
//                                     INT_SPECTRAL : spectral
//                IntScheme_FC   : Interpolation scheme for the face-centered variables
//                                 --> Supported schemes include
//                                     INT_MINMOD1D : MinMod-1D
//                                     INT_VANLEER  : vanLeer
//                                     INT_CQUAD    : conservative quadratic
//                                     INT_CQUAR    : conservative quartic
//                PrepUnit       : Whether or not to separate the prepared data into individual patches
//                                 --> UNIT_PATCH      : prepare data "patch by patch"
//                                     UNIT_PATCHGROUP : prepare data "patch group by patch group"
//                NSide          : Number of sibling directions to prepare data
//                                 --> NSIDE_00 (=  0) : do not prepare any sibling direction (equivalent to GhostSize=0)
//                                     NSIDE_06 (=  6) : prepare only sibling directions 0~5
//                                     NSIDE_26 (= 26) : prepare all sibling directions 0~25
//                IntPhase       : true --> ELBDM_WAVE  : Perform interpolation on rho/phase instead of real/imag parts
//                                                        --> TVarCC must contain _REAL and _IMAG
//                                      --> ELBDM_HYBRID: Perform interpolation on rho/phase instead of real/imag parts on wave levels
//                                                        --> TVarCC must contain _REAL and _IMAG
//                                                        Perform interpolation on rho/phase regardless of IntPhase on fluid levels
//                                                        --> TVarCC must contain _DENS and _PHAS
//                FluBC          : Fluid boundary condition
//                PotBC          : Gravity boundary condition
//                MinDens        : See MinEntr
//                MinPres        : See MinEntr
//                MinTemp        : See MinEntr
//                MinEntr        : Minimum allowed density/pressure/temperature/entropy in the output array (<0.0 ==> off)
//                                 --> MinDens can be applied to both _DENS and _TOTAL_DENS but cannot be applied to _PAR_DENS
//                                 --> Note that when preparing both density and real/imaginary parts for ELBDM, we do NOT
//                                     rescale wave functions after applying MinDens
//                                     --> We can have real^2+imag^2 != density in the prepared data!!
//                                     --> But currently it's not an issue since we never prepare density and wave functions
//                                         at the same time
//                                 --> Currently MinDens is applied in Flu_Prepare(), Flag_Real(), and Poi_Prepare_Rho()
//                                     --> The Guideline is to apply MinDens check only when ghost zones are required
//                                         (because density field is already stored in each patch and we don't want
//                                         Prepare_PatchData() to modify the existing data)
//                                 --> Currently MinPres is applied only in Flag_Real()
//                                     --> The Guideline is to apply MinPres check whenever _PRES or _TEMP is required
//                                         (because pressure field is NOT stored explicitly in each patch and thus existing data
//                                         may still have pressure < MinPres due to round-off errors)
//                DE_Consistency : Ensure the consistency between pressure, total energy density, and the dual-energy variable
//                                 when DUAL_ENERGY is on
//                                 --> Only apply to the ghost-zone interpolation on the assumption that the data stored
//                                     in all patches already satisfy this consistency check
//
// Return      :  OutputCC, OutputFC
//-------------------------------------------------------------------------------------------------------
void Prepare_PatchData( const int lv, const double PrepTime, real *OutputCC, real *OutputFC,
                        const int GhostSize, const int NPG, const int *PID0_List, long TVarCC, long TVarFC,
                        const IntScheme_t IntScheme_CC, const IntScheme_t IntScheme_FC, const PrepUnit_t PrepUnit,
                        const NSide_t NSide, const bool IntPhase, const OptFluBC_t FluBC[], const OptPotBC_t PotBC,
                        const real MinDens, const real MinPres, const real MinTemp, const real MinEntr, const bool DE_Consistency )
{

// nothing to do if there is no target patch group
   if ( NPG == 0 )   return;


// check
// --> do it even when disabling GAMER_DEBUG since this routine is critical and the check is inexpensive
// -----------------------------------------
   long AllVarCC = ( _TOTAL | _DERIVED );
#  ifdef GRAVITY
   AllVarCC |= _POTE;
#  endif
#  ifdef PARTICLE
   AllVarCC |= _PAR_DENS;
   AllVarCC |= _TOTAL_DENS;
#  endif
   if ( TVarCC & ~AllVarCC )   Aux_Error( ERROR_INFO, "unsupported parameter %s = %d !!\n", "TVarCC", TVarCC );

   long AllVarFC = 0;
#  ifdef MHD
   AllVarFC |= _MAG;
#  endif
   if ( TVarFC & ~AllVarFC )   Aux_Error( ERROR_INFO, "unsupported parameter %s = %d !!\n", "TVarFC", TVarFC );

   if ( MinDens >= (real)0.0  &&  MPI_Rank == 0 )
   {
#     if ( MODEL == HYDRO  ||  MODEL == ELBDM )
#     ifdef PARTICLE
      if ( !(TVarCC & _DENS)  &&  !(TVarCC & _TOTAL_DENS) )
         Aux_Message( stderr, "WARNING : MinDens (%13.7e) >= 0.0 but neither _DENS nor _TOTAL_DENS is found !!\n", MinDens );
#     else
      if ( !(TVarCC & _DENS) )
         Aux_Message( stderr, "WARNING : MinDens (%13.7e) >= 0.0 but _DENS is not found !!\n", MinDens );
#     endif
#     else
         Aux_Message( stderr, "WARNING : MinDens (%13.7e) >= 0.0 can only be applied to HYDRO/ELBDM !!\n", MinDens );
#     endif

#     if ( MODEL == ELBDM )
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( amr->use_wave_flag[lv] )
#     endif
      if (  ( TVarCC & _REAL )  ||  ( TVarCC & _IMAG )  )
         Aux_Message( stderr, "WARNING : real and imaginary parts are NOT rescaled after applying the minimum density check !!\n" );
#     endif // #if ( MODEL == ELBDM )
   }

   if ( MinPres >= (real)0.0  &&  MPI_Rank == 0 )
   {
#     if ( MODEL == HYDRO )
      if ( !(TVarCC & _PRES) )
         Aux_Message( stderr, "WARNING : MinPres (%13.7e) >= 0.0 but _PRES is not found !!\n", MinPres );
#     else
         Aux_Message( stderr, "WARNING : MinPres (%13.7e) >= 0.0 can only be applied to HYDRO !!\n", MinPres );
#     endif
   }

   if ( MinTemp >= (real)0.0  &&  MPI_Rank == 0 )
   {
#     if ( MODEL == HYDRO )
      if ( !(TVarCC & _TEMP) )
         Aux_Message( stderr, "WARNING : MinTemp (%13.7e) >= 0.0 but _TEMP is not found !!\n", MinTemp );
#     else
         Aux_Message( stderr, "WARNING : MinTemp (%13.7e) >= 0.0 can only be applied to HYDRO !!\n", MinTemp );
#     endif
   }

   if ( MinEntr >= (real)0.0  &&  MPI_Rank == 0 )
   {
#     if ( MODEL == HYDRO )
      if ( !(TVarCC & _ENTR) )
         Aux_Message( stderr, "WARNING : MinEntr (%13.7e) >= 0.0 but _ENTR is not found !!\n", MinEntr );
#     else
         Aux_Message( stderr, "WARNING : MinEntr (%13.7e) >= 0.0 can only be applied to HYDRO !!\n", MinEntr );
#     endif
   }

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if (   (  ( (TVarCC & _REAL) && !(TVarCC & _IMAG) )  ||  ( !(TVarCC & _REAL) && (TVarCC & _IMAG) )  )
          &&  amr->use_wave_flag[lv]   )
      Aux_Error( ERROR_INFO, "%s() for hybrid scheme currently requires that the real and imaginary parts of the wave function are prepared together !!\n"
                             "This is probably due to operation (b2-3-1) and can probably be easily fixed.\n", __FUNCTION__ );

   if ( ( TVarCC & _DENS )  &&  ( TVarCC & (_REAL | _IMAG) )  &&
        GhostSize > 0  &&  lv > 0  &&  amr->use_wave_flag[lv]  &&  !amr->use_wave_flag[lv-1] )
      Aux_Error( ERROR_INFO, "%s() for hybrid scheme currently does not support preparing the density "
                             "and the real/imaginary part of the wave function simultaneously !!\n"
                             "This is probably due to operation (b2-3-1) and can probably be easily fixed.\n", __FUNCTION__ );
#  endif

   if ( IntPhase )
   {
#     if ( MODEL == ELBDM )
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( amr->use_wave_flag[lv] ) {
#     endif
      if ( !(TVarCC & _REAL)  ||  !(TVarCC & _IMAG) )
      Aux_Error( ERROR_INFO, "real and/or imag parts are not found for phase interpolation in ELBDM !!\n" );
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      } // if ( amr->use_wave_flag[lv] )
#     endif

//    we have assumed in InterpolateGhostZone() that when adopting IntPhase this function will NOT prepare
//    anything other than wave function and, optionally, density
//    --> e.g., one cannot prepare wave function and potential at the same time when enabling IntPhase
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( amr->use_wave_flag[lv] )
#     endif
      if (  TVarCC & ~( _REAL | _IMAG | _DENS )  )
      Aux_Error( ERROR_INFO, "unsupported parameter %s = %d for IntPhase !!\n", "TVarCC", TVarCC );

#     else // #if ( MODEL == ELBDM )
      Aux_Error( ERROR_INFO, "\"interpolation on phase\" is useful only in ELBDM !!\n" );
#     endif // #if ( MODEL == ELBDM ) ... else ...
   }

   if ( FluBC == NULL )    Aux_Error( ERROR_INFO, "FluBC == NULL !!\n" );

   for (int f=0; f<6; f++)
   {
      if ( FluBC[f] != BC_FLU_PERIODIC    &&  FluBC[f] != BC_FLU_OUTFLOW  &&
           FluBC[f] != BC_FLU_REFLECTING  &&  FluBC[f] != BC_FLU_DIODE    &&
           FluBC[f] != BC_FLU_USER )
         Aux_Error( ERROR_INFO, "unsupported parameter FluBC[%d] = %d !!\n", f, FluBC[f] );

#     if ( MODEL != HYDRO )
      if ( FluBC[f] == BC_FLU_OUTFLOW )
         Aux_Error( ERROR_INFO, "outflow boundary condition (OPT__BC_FLU=2) only works with HYDRO !!\n" );

      if ( FluBC[f] == BC_FLU_REFLECTING )
         Aux_Error( ERROR_INFO, "reflecting boundary condition (OPT__BC_FLU=3) only works with HYDRO !!\n" );

      if ( FluBC[f] == BC_FLU_DIODE )
         Aux_Error( ERROR_INFO, "diode boundary condition (OPT__BC_FLU=5) only works with HYDRO !!\n" );
#     endif
   }

   if ( OPT__DT_LEVEL == DT_LEVEL_SHARED  &&  OPT__INT_TIME )
      Aux_Error( ERROR_INFO, "OPT__INT_TIME should be disabled when \"OPT__DT_LEVEL == DT_LEVEL_SHARED\" !!\n" );

   if ( MPI_Rank == 0 )
   if (  ( NSide == NSIDE_00  &&  GhostSize != 0 )  ||  ( NSide != NSIDE_00  &&  GhostSize == 0 )  )
      Aux_Message( stderr, "WARNING : inconsistent NSide (%d) and GhostSize (%d) !!\n", NSide, GhostSize );

#  ifdef PARTICLE
   if ( TVarCC & _PAR_DENS  ||  TVarCC & _TOTAL_DENS )
   {
//    because we only collect particles from nearby 26 sibling patches
      if ( GhostSize > PS1 - amr->Par->GhostSize )
         Aux_Error( ERROR_INFO, "GhostSize (%d) > maximum allowed (%d) when preparing mass density with particles!!\n",
                    GhostSize, PS1 - amr->Par->GhostSize );

      if ( ! ParDensArray_Initialized )
         Aux_Error( ERROR_INFO, "must call Prepare_PatchData_InitParticleDensityArray() in advance !!\n" );
   }

// _DENS, _PAR_DENS, and _TOTAL_DENS do not work together (actually we should be able to support _DENS + _PAR_DENS)
   if (  ( TVarCC & _DENS && TVarCC & _PAR_DENS )  ||  ( TVarCC & _DENS && TVarCC & _TOTAL_DENS )  )
      Aux_Error( ERROR_INFO, "_DENS, _PAR_DENS, and _TOTAL_DENS cannot work together !!\n" );

// for PAR_ONLY, we have _TOTAL_DENS == _PAR_DENS
#  if ( MODEL != PAR_ONLY )
   if ( TVarCC & _TOTAL_DENS  &&  TVarCC & _PAR_DENS )
      Aux_Error( ERROR_INFO, "_DENS, _PAR_DENS, and _TOTAL_DENS cannot work together !!\n" );
#  endif
#  endif // #ifdef PARTICLE

// target patches must be real patches
   for (int TID=0; TID<NPG; TID++)
      if ( PID0_List[TID] < 0  ||  PID0_List[TID] >= amr->NPatchComma[lv][1] )
         Aux_Error( ERROR_INFO, "incorrect target PID %d (NReal = %d) !!\n", PID0_List[TID], amr->NPatchComma[lv][1] );

#  if ( MODEL == HYDRO )
   if (  ( TVarCC & _TEMP )  &&  EoS_DensEint2Temp_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensEint2Temp_CPUPtr == NULL !!\n" );

   if (  ( TVarCC & _ENTR )  &&  EoS_DensEint2Entr_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensEint2Entr_CPUPtr == NULL !!\n" );
#  endif

#  ifdef SRHD
   if ( TVarCC & _ENTR )
      Aux_Error( ERROR_INFO, "SRHD does not support computing entropy !!\n" );
#  endif
// -----------------------------------------
// end of check


   const double dh               = amr->dh[lv];
   const int    PGSize1D_CC      = 2*( PS1 + GhostSize );   // width of a single patch group including ghost zones
   const int    PGSize3D_CC      = CUBE( PGSize1D_CC );
   const int    GhostSize_Padded = GhostSize + (GhostSize&1);
   const int    PGSize1D_FC      = PGSize1D_CC + 1;
   const int    PGSize3D_FC      = PGSize1D_FC*SQR(PGSize1D_CC);

#  if   ( MODEL == HYDRO )
   const bool PrepVx           = ( TVarCC & _VELX    ) ? true : false;
   const bool PrepVy           = ( TVarCC & _VELY    ) ? true : false;
   const bool PrepVz           = ( TVarCC & _VELZ    ) ? true : false;
   const bool PrepPres         = ( TVarCC & _PRES    ) ? true : false;
   const bool PrepTemp         = ( TVarCC & _TEMP    ) ? true : false;
   const bool PrepEntr         = ( TVarCC & _ENTR    ) ? true : false;
   const bool PrepEint         = ( TVarCC & _EINT    ) ? true : false;
#  ifdef MHD
   const bool PrepMagX_CC      = ( TVarCC & _MAGX_CC ) ? true : false;
   const bool PrepMagY_CC      = ( TVarCC & _MAGY_CC ) ? true : false;
   const bool PrepMagZ_CC      = ( TVarCC & _MAGZ_CC ) ? true : false;
   const bool PrepMagE_CC      = ( TVarCC & _MAGE_CC ) ? true : false;
   const bool PrepMagCC        = ( PrepMagX_CC || PrepMagY_CC || PrepMagZ_CC || PrepMagE_CC );
#  endif

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif

#  ifdef GRAVITY
   const bool PrepPot = ( TVarCC & _POTE ) ? true : false;
#  endif

#  ifdef PARTICLE
// note that these two options cannot be turned on at the same time
// --> and we set PrepTotalDens == true ONLY when there are density fields other than particles
// --> for PAR_ONLY mode, we always set PrepTotalDens == false to avoid confusion
#  if ( MODEL == PAR_ONLY )
   const bool PrepParOnlyDens = ( TVarCC & _PAR_DENS  ||  TVarCC & _TOTAL_DENS ) ? true : false;
   const bool PrepTotalDens   = false;
#  else
   const bool PrepParOnlyDens = ( TVarCC & _PAR_DENS   ) ? true : false;
   const bool PrepTotalDens   = ( TVarCC & _TOTAL_DENS ) ? true : false;
#  endif

// turn on _DENS automatically for preparing total density
   if ( PrepTotalDens )    TVarCC |= _DENS;
#  endif // #ifdef PARTICLE


// TVarCCIdxList_Flu: list recording the target cell-centered fluid and passive variable indices (e.g., [0 ... NCOMP_TOTAL-1] )
// TVarCCList_Der   : list recording the target cell-centered derived variable (e.g., _VELX, _PRES)
// TVarFCIdxList    : list recording the target face-centered variable indices (e.g., [0 ... NCOMP_MAG-1])
   const int NVarCC_Der_Max = 20;   // increase it when the maximum number of derived fields exceeds it
   long TVarCCList_Der[NVarCC_Der_Max];
   int  NTSib[26], *TSib[26], NVarCC_Flu, NVarCC_Der, NVarCC_Tot, TVarCCIdxList_Flu[NCOMP_TOTAL];

// set up the target sibling indices for InterpolateGhostZone()
   SetTargetSibling( NTSib, TSib );

// determine the cell-centered fluid components to be prepared
// --> assuming that _VAR_NAME = 1L<<VAR_NAME (e.g., _DENS == 1L<<DENS)
// --> it also determines the order of variables stored in OutputCC (which is the same as patch->fluid[])
   NVarCC_Flu = 0;
   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVarCC & (1L<<v) )    TVarCCIdxList_Flu[ NVarCC_Flu++ ] = v;

   NVarCC_Der = 0;

#  if   ( MODEL == HYDRO )
   if ( PrepVx      )   TVarCCList_Der[ NVarCC_Der ++ ] = _VELX;
   if ( PrepVy      )   TVarCCList_Der[ NVarCC_Der ++ ] = _VELY;
   if ( PrepVz      )   TVarCCList_Der[ NVarCC_Der ++ ] = _VELZ;
   if ( PrepPres    )   TVarCCList_Der[ NVarCC_Der ++ ] = _PRES;
   if ( PrepTemp    )   TVarCCList_Der[ NVarCC_Der ++ ] = _TEMP;
   if ( PrepEntr    )   TVarCCList_Der[ NVarCC_Der ++ ] = _ENTR;
   if ( PrepEint    )   TVarCCList_Der[ NVarCC_Der ++ ] = _EINT;
#  ifdef MHD
   if ( PrepMagX_CC )   TVarCCList_Der[ NVarCC_Der ++ ] = _MAGX_CC;
   if ( PrepMagY_CC )   TVarCCList_Der[ NVarCC_Der ++ ] = _MAGY_CC;
   if ( PrepMagZ_CC )   TVarCCList_Der[ NVarCC_Der ++ ] = _MAGZ_CC;
   if ( PrepMagE_CC )   TVarCCList_Der[ NVarCC_Der ++ ] = _MAGE_CC;
#  endif

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif

   if ( NVarCC_Der > NVarCC_Der_Max )  Aux_Error( ERROR_INFO, "NVarCC_Der (%d) > NVarCC_Der_Max (%d) !!\n", NVarCC_Der, NVarCC_Der_Max );

   NVarCC_Tot = NVarCC_Flu + NVarCC_Der;

#  ifdef GRAVITY
   if ( PrepPot )    NVarCC_Tot ++;
#  endif

// do not increase NVarCC_Tot for PrepTotalDens since _DENS is already turned on automatically for that
#  ifdef PARTICLE
   if ( PrepParOnlyDens )  NVarCC_Tot ++;
#  endif


// determine the face-centered variables to be prepared
// --> currently we do not consider any face-centered derived variable
// --> assuming that _VAR_NAME = 1L<<VAR_NAME (e.g., _MAGX == 1L<<MAGX)
// --> it also determines the order of variables stored in OutputFC (which is the same as patch->magnetic[])
   int NVarFC_Tot = 0;

#  ifdef MHD
   const bool PrepMagFC = ( TVarFC & _MAG ) ? true : false;
   int TVarFCIdxList[NCOMP_MAG];

   for (int v=0; v<NCOMP_MAG; v++)
      if ( TVarFC & (1L<<v) )    TVarFCIdxList[ NVarFC_Tot++ ] = v;

#  else
// currently no other models require any face-centered variable
   int *TVarFCIdxList = NULL;
#  endif


// nothing to do if no target variable is found
   if ( NVarCC_Tot == 0  &&  NVarFC_Tot == 0  &&  MPI_Rank == 0 )
   {
      Aux_Message( stderr, "WARNING : no target variable is found in %s() !!\n", __FUNCTION__ );
      return;
   }


// validate the output pointers
#  ifdef GAMER_DEBUG
   if ( NVarCC_Tot > 0  &&  OutputCC == NULL )
      Aux_Error( ERROR_INFO, "OutputCC == NULL (NVarCC_Tot %d) !!\n", NVarCC_Tot );

   if ( NVarFC_Tot > 0  &&  OutputFC == NULL )
      Aux_Error( ERROR_INFO, "OutputFC == NULL (NVarFC_Tot %d) !!\n", NVarFC_Tot );
#  endif


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
   }

// magnetic field
#  ifdef MHD
   bool MagIntTime;
   int  MagSg, MagSg_IntT;
   real MagWeighting, MagWeighting_IntT;

// check PrepPres, PrepTemp, PrepEntr, and PrepEint since they also require B field
   if ( PrepMagFC || PrepMagCC || PrepPres || PrepTemp || PrepEntr || PrepEint )
   {
      const int Sg0 = amr->MagSg[lv];
      SetTempIntPara( lv, Sg0, PrepTime, amr->MagSgTime[lv][Sg0], amr->MagSgTime[lv][1-Sg0],
                      MagIntTime, MagSg, MagSg_IntT, MagWeighting, MagWeighting_IntT );
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
   }
#  endif // #ifdef GRAVITY


// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   int BC_Face[26], BC_Face_tmp[3];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1  &&  FluBC[BC_Face_tmp[2]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1  &&  FluBC[BC_Face_tmp[1]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1  &&  FluBC[BC_Face_tmp[0]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[0];
      else                                                                             BC_Face[s] = NULL_INT;
   }


// constant settings used by Par_MassAssignment()
#  ifdef MASSIVE_PARTICLES
   const bool InitZero_No       = false;
   const bool Periodic_Check[3] = { FluBC[0]==BC_FLU_PERIODIC, FluBC[2]==BC_FLU_PERIODIC, FluBC[4]==BC_FLU_PERIODIC };
   const int  PeriodicNCell[3]  = { NX0_TOT[0]*(1<<lv),
                                    NX0_TOT[1]*(1<<lv),
                                    NX0_TOT[2]*(1<<lv) };
   const bool UnitDens_No       = false;
   const bool CheckFarAway_Yes  = true;
#  endif


// start to prepare data
#  pragma omp parallel
   {
//    thread-private variables
      int    J, K, I2, J2, K2, Idx1, Idx2, PID0, TVarCCIdx_Flu, TVarFCIdx, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];
      double xyz0[3];   // corner coordinates for the user-specified B.C.

//    fluid variables for the EoS routines
#     if ( MODEL == HYDRO )
#     if ( EOS == EOS_GAMMA  ||  EOS == EOS_ISOTHERMAL )
      const int NFluForEoS = NCOMP_FLUID;    // don't need passsive scalars in EOS_GAMMA/EOS_ISOTHERMAL
#     else
      const int NFluForEoS = NCOMP_TOTAL;
#     endif
      real FluidForEoS[NFluForEoS];
#     endif

//    Data1PG_CC/FC: array to store the prepared cell-centered/face-centered data of one patch group
//                   (including the ghost-zone data)
//    --> for PrepUnit == UNIT_PATCHGROUP, these pointers point to OutputCC/FC directly (which will be set later)
//        for PrepUnit == UNIT_PATCH, these arrays will be copied to different patches in OutputCC/FC later
      real *Data1PG_CC     = ( PrepUnit == UNIT_PATCH ) ? new real [ NVarCC_Tot*PGSize3D_CC ] : NULL;
      real *Data1PG_CC_Ptr = NULL;
      real *Data1PG_FC     = ( PrepUnit == UNIT_PATCH ) ? new real [ NVarFC_Tot*PGSize3D_FC ] : NULL;
      real *Data1PG_FC_Ptr = NULL;


//    IntData_CC/FC: arrays to store the interpolated cell-/face-centered results
//    --> allocate it only once but with the maximum required size to reduce the number of memory allocations
      real *IntData_CC = new real [ NVarCC_Tot*PS2*PS2*(GhostSize_Padded  ) ];
      real *IntData_FC = new real [ NVarFC_Tot*PS2*PS2*(GhostSize_Padded+1) ];


//    B field on the coarse-fine interfaces for the divergence-preserving interpolation
//    --> allocate it only once but with the maximum required size to reduce the number of memory allocations
#     ifdef MHD
      real *FInterface_Data = NULL;

      if ( NVarFC_Tot > 0 )   FInterface_Data = new real [ SQR(PS2) + 4*PS2*GhostSize_Padded ];
#     endif

//    IntData_CC_IntTime: for temporal interpolation on density and phase in ELBDM
#     if ( MODEL == ELBDM )
      real *IntData_CC_IntTime = (  IntPhase  &&  OPT__INT_TIME  &&  lv > 0  &&
                                   !Mis_CompareRealValue( PrepTime, amr->FluSgTime[lv-1][  amr->FluSg[lv-1]], NULL, false )  &&
                                   !Mis_CompareRealValue( PrepTime, amr->FluSgTime[lv-1][1-amr->FluSg[lv-1]], NULL, false )  )
                                 ? new real [ 2*PS2*PS2*GhostSize_Padded ] : NULL;
#     else
      real *IntData_CC_IntTime = NULL;
#     endif


//    prepare eight nearby patches (one patch group) at a time
#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

#        ifdef GAMER_DEBUG
         if ( PID0 < 0  ||  PID0 >= amr->num[lv] )
            Aux_Error( ERROR_INFO, "PID0 (%d) is not within the correct range [%d ... %d] !!\n", PID0, 0, amr->num[lv]-1 );

         if ( PID0%8 != 0 )
            Aux_Error( ERROR_INFO, "PID0 (%d) does not have LocalID == 0 !!\n", PID0 );
#        endif

         for (int d=0; d<3; d++)    xyz0[d] = amr->patch[0][lv][PID0]->EdgeL[d] + (0.5-GhostSize)*dh;

//       Data1PG_CC/FC point to OutputCC/FC directly for PrepUnit == UNIT_PATCHGROUP
         if ( PrepUnit == UNIT_PATCHGROUP )
         {
            Data1PG_CC = OutputCC + TID*NVarCC_Tot*PGSize3D_CC;
            Data1PG_FC = OutputFC + TID*NVarFC_Tot*PGSize3D_FC;
         }


//       a. fill out the central region of Data1PG_CC[]/FC[] (ghost zones will be filled out later)
// ------------------------------------------------------------------------------------------------------------
         for (int LocalID=0; LocalID<8; LocalID++ )
         {
            const int PID    = PID0 + LocalID;
            const int Disp_i = TABLE_02( LocalID, 'x', GhostSize, GhostSize+PS1 );
            const int Disp_j = TABLE_02( LocalID, 'y', GhostSize, GhostSize+PS1 );
            const int Disp_k = TABLE_02( LocalID, 'z', GhostSize, GhostSize+PS1 );

            Data1PG_CC_Ptr = Data1PG_CC;
            Data1PG_FC_Ptr = Data1PG_FC;

//          (a1) fluid data (cell-centered)
            for (int v=0; v<NVarCC_Flu; v++)
            {
               TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][PID]->fluid[TVarCCIdx_Flu][k][j][i];

                  if ( FluIntTime ) // temporal interpolation
                  Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                         + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][PID]->fluid[TVarCCIdx_Flu][k][j][i];
                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            }


//          (a2) derived variables (cell-centered)
#           if   ( MODEL == HYDRO )
            if ( PrepVx )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][PID]->fluid[MOMX][k][j][i] /
                                         amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];

                  if ( FluIntTime ) // temporal interpolation
                  Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                         + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMX][k][j][i] /
                                                               amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k][j][i] );
                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepVx )

            if ( PrepVy )
            {
               for (int k=0; k<PS1; k++)    {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)    {  J    = j + Disp_j;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)    {

                  Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][PID]->fluid[MOMY][k][j][i] /
                                         amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];

                  if ( FluIntTime ) // temporal interpolation
                  Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                         + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMY][k][j][i] /
                                                               amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k][j][i] );
                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepVy )

            if ( PrepVz )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][PID]->fluid[MOMZ][k][j][i] /
                                         amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];

                  if ( FluIntTime ) // temporal interpolation
                  Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                         + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][PID]->fluid[MOMZ][k][j][i] /
                                                               amr->patch[FluSg_IntT][lv][PID]->fluid[DENS][k][j][i] );
                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepVz )

            if ( PrepPres )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//###REVISE: support dual energy
#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  Data1PG_CC_Ptr[Idx1] = Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                         FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                         (MinPres>=(real)0.0), MinPres, Emag,
                                                         EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                         EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                         h_EoS_Table, NULL );

                  if ( FluIntTime ) // temporal interpolation
                  {
                     for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT );
#                    else
                     const real Emag = NULL_REAL;
#                    endif
                     Data1PG_CC_Ptr[Idx1] =
                        FluWeighting     *Data1PG_CC_Ptr[Idx1]
                      + FluWeighting_IntT*Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                          FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                          (MinPres>=(real)0.0), MinPres, Emag,
                                                          EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                          EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                          h_EoS_Table, NULL );
                  }

                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepPres )

            if ( PrepTemp )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//###REVISE: support dual energy
#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  Data1PG_CC_Ptr[Idx1] = Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                         FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                         (MinTemp>=(real)0.0), MinTemp, Emag,
                                                         EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                         EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

                  if ( FluIntTime ) // temporal interpolation
                  {
                     for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT );
#                    else
                     const real Emag = NULL_REAL;
#                    endif
                     Data1PG_CC_Ptr[Idx1] =
                        FluWeighting     *Data1PG_CC_Ptr[Idx1]
                      + FluWeighting_IntT*Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                          FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                          (MinTemp>=(real)0.0), MinTemp, Emag,
                                                          EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                  }

                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepTemp )

#           ifndef SRHD
            if ( PrepEntr )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  Data1PG_CC_Ptr[Idx1] = Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                         FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                         (MinEntr>=(real)0.0), MinEntr, Emag,
                                                         EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                         h_EoS_Table );

                  if ( FluIntTime ) // temporal interpolation
                  {
                     for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT );
#                    else
                     const real Emag = NULL_REAL;
#                    endif
                     Data1PG_CC_Ptr[Idx1] =
                        FluWeighting     *Data1PG_CC_Ptr[Idx1]
                      + FluWeighting_IntT*Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                          FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                          (MinEntr>=(real)0.0), MinEntr, Emag,
                                                          EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                          h_EoS_Table );
                  }

                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepEntr )
#           endif // #ifndef SRHD

            if ( PrepEint )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//###REVISE: support dual energy
#                 ifdef MHD
                  const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
                  const bool CheckMinEint_No = false; // floor value is not supported for now
                  Data1PG_CC_Ptr[Idx1] = Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                         FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                         CheckMinEint_No, NULL_REAL, Emag,
                                                         EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                         EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

                  if ( FluIntTime ) // temporal interpolation
                  {
                     for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][PID]->fluid[v][k][j][i];

#                    ifdef MHD
                     const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT );
#                    else
                     const real Emag = NULL_REAL;
#                    endif
                     Data1PG_CC_Ptr[Idx1] =
                        FluWeighting     *Data1PG_CC_Ptr[Idx1]
                      + FluWeighting_IntT*Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                          FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                          CheckMinEint_No, NULL_REAL, Emag,
                                                          EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                  }

                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepEint )

#           ifdef MHD
            if ( PrepMagCC )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  real B_CC[NCOMP_MAG];
                  int  OffsetB = 0;

                  MHD_GetCellCenteredBFieldInPatch( B_CC, lv, PID, i, j, k, MagSg );

                  if ( PrepMagX_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = B_CC[MAGX];  OffsetB += PGSize3D_CC; }
                  if ( PrepMagY_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = B_CC[MAGY];  OffsetB += PGSize3D_CC; }
                  if ( PrepMagZ_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = B_CC[MAGZ];  OffsetB += PGSize3D_CC; }
                  if ( PrepMagE_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = (real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) ); }

                  if ( FluIntTime ) // temporal interpolation
                  {
                     OffsetB = 0;

                     MHD_GetCellCenteredBFieldInPatch( B_CC, lv, PID, i, j, k, MagSg_IntT );

                     if ( PrepMagX_CC ) {
                        Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                           + FluWeighting_IntT*B_CC[MAGX];
                        OffsetB += PGSize3D_CC;
                     }

                     if ( PrepMagY_CC ) {
                        Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                           + FluWeighting_IntT*B_CC[MAGY];
                        OffsetB += PGSize3D_CC;
                     }

                     if ( PrepMagZ_CC ) {
                        Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                           + FluWeighting_IntT*B_CC[MAGZ];
                        OffsetB += PGSize3D_CC;
                     }

                     if ( PrepMagE_CC ) {
                        Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                           + FluWeighting_IntT*(real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) );
                     }
                  } // if ( FluIntTime )

                  Idx1 ++;
               }}}

               if ( PrepMagX_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
               if ( PrepMagY_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
               if ( PrepMagZ_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
               if ( PrepMagE_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepMagCC )
#           endif // #ifdef MHD


#           elif ( MODEL == ELBDM )
//          no derived variables yet


#           else
#           error : unsupported MODEL !!
#           endif // MODEL


#           ifdef GRAVITY
//          (a3) potential data (cell-centered)
            if ( PrepPot )
            {
               for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
               for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                             Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=0; i<PS1; i++)  {

                  Data1PG_CC_Ptr[Idx1] = amr->patch[PotSg][lv][PID]->pot[k][j][i];

                  if ( PotIntTime ) // temporal interpolation
                  Data1PG_CC_Ptr[Idx1] =   PotWeighting     *Data1PG_CC_Ptr[Idx1]
                                         + PotWeighting_IntT*amr->patch[PotSg_IntT][lv][PID]->pot[k][j][i];
                  Idx1 ++;
               }}}

               Data1PG_CC_Ptr += PGSize3D_CC;
            } // if ( PrepPot )
#           endif // #ifdef GRAVITY


//          (a4) face-centered variables (e.g., magnetic field)
            for (int v=0; v<NVarFC_Tot; v++)
            {
               TVarFCIdx = TVarFCIdxList[v];

#              ifdef MHD

//             set array indices
               const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                                    ( TVarFCIdx == MAGY ) ? 1 :
                                    ( TVarFCIdx == MAGZ ) ? 2 : -1;
#              ifdef GAMER_DEBUG
               if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#              endif

               int ijk_s[3], ijk_e[3], size_i[3], size_o[3], idx_i, idx_o;    // s/e=start/end; i/o=in/out

               for (int d=0; d<3; d++)
               {
                  ijk_s [d] = 0;
                  ijk_e [d] = PS1;
                  size_i[d] = PS1;
                  size_o[d] = PGSize1D_CC;
               }

//             avoid redundant assignment of the patch interface values by copying the left interface
//             of only the left patch
               ijk_s [norm_dir] += TABLE_02( LocalID, 'x'+norm_dir, 0, 1 );
               ijk_e [norm_dir] ++;
               size_i[norm_dir] ++;
               size_o[norm_dir] ++;


//             copy data
               for (int k=ijk_s[2]; k<ijk_e[2]; k++)  {  K     = k + Disp_k;
               for (int j=ijk_s[1]; j<ijk_e[1]; j++)  {  J     = j + Disp_j;
                                                         idx_o = IDX321( ijk_s[0]+Disp_i, J, K, size_o[0], size_o[1] );
                                                         idx_i = IDX321( ijk_s[0],        j, k, size_i[0], size_i[1] );
               for (int i=ijk_s[0]; i<ijk_e[0]; i++)  {

                  Data1PG_FC_Ptr[idx_o] = amr->patch[MagSg][lv][PID]->magnetic[TVarFCIdx][idx_i];

                  if ( MagIntTime ) // temporal interpolation
                  Data1PG_FC_Ptr[idx_o] =   MagWeighting     *Data1PG_FC_Ptr[idx_o]
                                          + MagWeighting_IntT*amr->patch[MagSg_IntT][lv][PID]->magnetic[TVarFCIdx][idx_i];
                  idx_i ++;
                  idx_o ++;
               }}}

#              else
               Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#              endif // #ifdef MHD ... else ...

               Data1PG_FC_Ptr += PGSize3D_FC;
            } // for (int v=0; v<NVarFC_Tot; v++)

         } // for (int LocalID=0; LocalID<8; LocalID++ )


//       b. fill out the ghost zones of Data1PG_CC[]/FC[]
// ------------------------------------------------------------------------------------------------------------
//       direct memory copy
         for (int Side=0; Side<NSide; Side++)
         {
//          nothing to do if no ghost zone is required
            if ( GhostSize == 0 )   break;


            const int SibPID0 = Table_02( lv, PID0, Side );    // the 0th patch of the sibling patch group

//          (b1) if the target sibling patch exists --> just copy data from the nearby patches at the same level
            if ( SibPID0 >= 0 )
            {
               int loop[3], disp2[3];
               for (int d=0; d<3; d++)
               {
                  loop [d] = TABLE_01( Side, 'x'+d, GhostSize, PS1, GhostSize );
                  disp2[d] = TABLE_01( Side, 'x'+d, PS1-GhostSize, 0, 0 );
               }

//###OPTIMIZATION: simplify TABLE_03 and TABLE_04
               for (int Count=0; Count<TABLE_04(Side); Count++)
               {
                  const int LocalID = TABLE_03( Side, Count );
                  const int SibPID  = SibPID0 + LocalID;

                  int disp[3];
                  for (int d=0; d<3; d++)    disp[d] = Table_01( Side, 'x'+d, Count, GhostSize );

                  Data1PG_CC_Ptr = Data1PG_CC;
                  Data1PG_FC_Ptr = Data1PG_FC;

//                (b1-1) fluid data (cell-centered)
                  for (int v=0; v<NVarCC_Flu; v++)
                  {
                     TVarCCIdx_Flu = TVarCCIdxList_Flu[v];

                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][SibPID]->fluid[TVarCCIdx_Flu][K2][J2][I2];

                        if ( FluIntTime ) // temporal interpolation
                        Data1PG_CC_Ptr[Idx1] =
                           FluWeighting     *Data1PG_CC_Ptr[Idx1]
                         + FluWeighting_IntT*amr->patch[FluSg_IntT][lv][SibPID]->fluid[TVarCCIdx_Flu][K2][J2][I2];

                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // for (int v=0; v<NVarCC_Flu; v++)


//                (b1-2) derived variables (cell-centered)
#                 if   ( MODEL == HYDRO )
                  if ( PrepVx )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][SibPID]->fluid[MOMX][K2][J2][I2] /
                                               amr->patch[FluSg][lv][SibPID]->fluid[DENS][K2][J2][I2];

                        if ( FluIntTime ) // temporal interpolation
                        Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                               + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMX][K2][J2][I2] /
                                                                     amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][K2][J2][I2] );
                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepVx )

                  if ( PrepVy )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][SibPID]->fluid[MOMY][K2][J2][I2] /
                                               amr->patch[FluSg][lv][SibPID]->fluid[DENS][K2][J2][I2];

                        if ( FluIntTime ) // temporal interpolation
                        Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                               + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMY][K2][J2][I2] /
                                                                     amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][K2][J2][I2] );
                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepVy )

                  if ( PrepVz )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        Data1PG_CC_Ptr[Idx1] = amr->patch[FluSg][lv][SibPID]->fluid[MOMZ][K2][J2][I2] /
                                               amr->patch[FluSg][lv][SibPID]->fluid[DENS][K2][J2][I2];

                        if ( FluIntTime ) // temporal interpolation
                        Data1PG_CC_Ptr[Idx1] =   FluWeighting     *Data1PG_CC_Ptr[Idx1]
                                               + FluWeighting_IntT*( amr->patch[FluSg_IntT][lv][SibPID]->fluid[MOMZ][K2][J2][I2] /
                                                                     amr->patch[FluSg_IntT][lv][SibPID]->fluid[DENS][K2][J2][I2] );
                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepVz )

                  if ( PrepPres )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][K2][J2][I2];

//###REVISE: support dual energy
#                       ifdef MHD
                        const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg );
#                       else
                        const real Emag = NULL_REAL;
#                       endif
                        Data1PG_CC_Ptr[Idx1] = Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                               FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                               (MinPres>=(real)0.0), MinPres, Emag,
                                                               EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                               EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                               h_EoS_Table, NULL );

                        if ( FluIntTime ) // temporal interpolation
                        {
                           for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][K2][J2][I2];

#                          ifdef MHD
                           const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg_IntT );
#                          else
                           const real Emag = NULL_REAL;
#                          endif
                           Data1PG_CC_Ptr[Idx1] =
                              FluWeighting     *Data1PG_CC_Ptr[Idx1]
                            + FluWeighting_IntT*Hydro_Con2Pres( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                                FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                                (MinPres>=(real)0.0), MinPres, Emag,
                                                                EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                                EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                                h_EoS_Table, NULL );
                        }

                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepPres )

                  if ( PrepTemp )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][K2][J2][I2];

//###REVISE: support dual energy
#                       ifdef MHD
                        const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg );
#                       else
                        const real Emag = NULL_REAL;
#                       endif
                        Data1PG_CC_Ptr[Idx1] = Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                               FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                               (MinTemp>=(real)0.0), MinTemp, Emag,
                                                               EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                               EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

                        if ( FluIntTime ) // temporal interpolation
                        {
                           for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][K2][J2][I2];

#                          ifdef MHD
                           const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg_IntT );
#                          else
                           const real Emag = NULL_REAL;
#                          endif
                           Data1PG_CC_Ptr[Idx1] =
                              FluWeighting     *Data1PG_CC_Ptr[Idx1]
                            + FluWeighting_IntT*Hydro_Con2Temp( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                                FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                                (MinTemp>=(real)0.0), MinTemp, Emag,
                                                                EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                        }

                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepTemp )

#                 ifndef SRHD
                  if ( PrepEntr )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][K2][J2][I2];

#                       ifdef MHD
                        const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg );
#                       else
                        const real Emag = NULL_REAL;
#                       endif
                        Data1PG_CC_Ptr[Idx1] = Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                               FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                               (MinEntr>=(real)0.0), MinEntr, Emag,
                                                               EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                               h_EoS_Table );

                        if ( FluIntTime ) // temporal interpolation
                        {
                           for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][K2][J2][I2];

#                          ifdef MHD
                           const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg_IntT );
#                          else
                           const real Emag = NULL_REAL;
#                          endif
                           Data1PG_CC_Ptr[Idx1] =
                              FluWeighting     *Data1PG_CC_Ptr[Idx1]
                            + FluWeighting_IntT*Hydro_Con2Entr( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                                FluidForEoS[MOMZ], FluidForEoS[ENGY], FluidForEoS+NCOMP_FLUID,
                                                                (MinEntr>=(real)0.0), MinEntr, Emag,
                                                                EoS_DensEint2Entr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                                h_EoS_Table );
                        }

                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepEntr )
#                 endif // #ifndef SRHD

                  if ( PrepEint )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg][lv][SibPID]->fluid[v][K2][J2][I2];

//###REVISE: support dual energy
#                       ifdef MHD
                        const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg );
#                       else
                        const real Emag = NULL_REAL;
#                       endif
                        const bool CheckMinEint_No = false; // floor value is not supported for now
                        Data1PG_CC_Ptr[Idx1] = Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                               FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                               CheckMinEint_No, NULL_REAL, Emag,
                                                               EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                               EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

                        if ( FluIntTime ) // temporal interpolation
                        {
                           for (int v=0; v<NFluForEoS; v++)    FluidForEoS[v] = amr->patch[FluSg_IntT][lv][SibPID]->fluid[v][K2][J2][I2];

#                          ifdef MHD
                           const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, SibPID, I2, J2, K2, MagSg_IntT );
#                          else
                           const real Emag = NULL_REAL;
#                          endif
                           Data1PG_CC_Ptr[Idx1] =
                              FluWeighting     *Data1PG_CC_Ptr[Idx1]
                            + FluWeighting_IntT*Hydro_Con2Eint( FluidForEoS[DENS], FluidForEoS[MOMX], FluidForEoS[MOMY],
                                                                FluidForEoS[MOMZ], FluidForEoS[ENGY],
                                                                CheckMinEint_No, NULL_REAL, Emag,
                                                                EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                        }

                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepEint )

#                 ifdef MHD
                  if ( PrepMagCC )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        real B_CC[NCOMP_MAG];
                        int  OffsetB = 0;

                        MHD_GetCellCenteredBFieldInPatch( B_CC, lv, SibPID, I2, J2, K2, MagSg );

                        if ( PrepMagX_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = B_CC[MAGX];  OffsetB += PGSize3D_CC; }
                        if ( PrepMagY_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = B_CC[MAGY];  OffsetB += PGSize3D_CC; }
                        if ( PrepMagZ_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = B_CC[MAGZ];  OffsetB += PGSize3D_CC; }
                        if ( PrepMagE_CC ) { Data1PG_CC_Ptr[ Idx1 + OffsetB ] = (real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) ); }

                        if ( FluIntTime ) // temporal interpolation
                        {
                           OffsetB = 0;

                           MHD_GetCellCenteredBFieldInPatch( B_CC, lv, SibPID, I2, J2, K2, MagSg_IntT );

                           if ( PrepMagX_CC ) {
                              Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                                 + FluWeighting_IntT*B_CC[MAGX];
                              OffsetB += PGSize3D_CC;
                           }

                           if ( PrepMagY_CC ) {
                              Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                                 + FluWeighting_IntT*B_CC[MAGY];
                              OffsetB += PGSize3D_CC;
                           }

                           if ( PrepMagZ_CC ) {
                              Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                                 + FluWeighting_IntT*B_CC[MAGZ];
                              OffsetB += PGSize3D_CC;
                           }

                           if ( PrepMagE_CC ) {
                              Data1PG_CC_Ptr[ Idx1 + OffsetB ] =   FluWeighting     *Data1PG_CC_Ptr[ Idx1 + OffsetB ]
                                                                 + FluWeighting_IntT*(real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) );
                           }
                        } // if ( FluIntTime )

                        Idx1 ++;
                     }}}

                     if ( PrepMagX_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
                     if ( PrepMagY_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
                     if ( PrepMagZ_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
                     if ( PrepMagE_CC )   Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepMagCC )
#                 endif // #ifdef MHD


#                 elif ( MODEL == ELBDM )
//                no derived variables yet

#                 else
#                 error : unsupported MODEL !!
#                 endif // MODEL


#                 ifdef GRAVITY
//                (b1-3) potential data (cell-centered)
                  if ( PrepPot )
                  {
                     for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                     for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                      Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                     for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {

                        Data1PG_CC_Ptr[Idx1] = amr->patch[PotSg][lv][SibPID]->pot[K2][J2][I2];

                        if ( PotIntTime ) // temporal interpolation
                        Data1PG_CC_Ptr[Idx1] =   PotWeighting     *Data1PG_CC_Ptr[Idx1]
                                               + PotWeighting_IntT*amr->patch[PotSg_IntT][lv][SibPID]->pot[K2][J2][I2];
                        Idx1 ++;
                     }}}

                     Data1PG_CC_Ptr += PGSize3D_CC;
                  } // if ( PrepPot )
#                 endif // #ifdef GRAVITY


//                (b1-4) face-centered variables
                  for (int v=0; v<NVarFC_Tot; v++)
                  {
                     TVarFCIdx = TVarFCIdxList[v];

#                    ifdef MHD
//                   set array indices
                     const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                                          ( TVarFCIdx == MAGY ) ? 1 :
                                          ( TVarFCIdx == MAGZ ) ? 2 : -1;
#                    ifdef GAMER_DEBUG
                     if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#                    endif

                     int ijk_s[3], ijk_e[3], size_i[3], size_o[3], disp_o[3], idx_i, idx_o;  // s/e=start/end; i/o=in/out

                     for (int d=0; d<3; d++)
                     {
                        ijk_s [d] = disp2[d];
                        ijk_e [d] = disp2[d] + loop[d];
                        size_i[d] = PS1;
                        size_o[d] = PGSize1D_CC;
                        disp_o[d] = disp[d] - disp2[d];
                     }

//                   avoid redundant assignment of the patch interface values **within the same patch group**
//                   by copying the left interface of only the left patch
//                   --> but always assign the interface values of patches **outside** the target patch group (PID0)
//                       --> to provide B field on the C-F interfaces for the divergence-free interpolation
//                       --> but it will result in redundant data copy on the F-F interfaces
//                           --> the checks of (Side<6) below only avoid the redundant copies on the F-F interfaces
//                               between the central 8 patches and their sibling patches but do not remove the
//                               redundant copies on the F-F interfaces between those sibling patches
                     ijk_s [norm_dir] += TABLE_01( Side, 'x'+norm_dir, 0, TABLE_02(LocalID,'x'+norm_dir,0,1), (Side<6)?1:0 );
                     ijk_e [norm_dir] += TABLE_01( Side, 'x'+norm_dir, (Side<6)?0:1, 1, 1 );
                     size_i[norm_dir] ++;
                     size_o[norm_dir] ++;


//                   copy data
                     for (int k=ijk_s[2]; k<ijk_e[2]; k++)  {  K     = k + disp_o[2];
                     for (int j=ijk_s[1]; j<ijk_e[1]; j++)  {  J     = j + disp_o[1];
                                                               idx_o = IDX321( ijk_s[0]+disp_o[0], J, K, size_o[0], size_o[1] );
                                                               idx_i = IDX321( ijk_s[0],           j, k, size_i[0], size_i[1] );
                     for (int i=ijk_s[0]; i<ijk_e[0]; i++)  {

                        Data1PG_FC_Ptr[idx_o] = amr->patch[MagSg][lv][SibPID]->magnetic[TVarFCIdx][idx_i];

                        if ( MagIntTime ) // temporal interpolation
                        Data1PG_FC_Ptr[idx_o] =
                           MagWeighting     *Data1PG_FC_Ptr[idx_o]
                         + MagWeighting_IntT*amr->patch[MagSg_IntT][lv][SibPID]->magnetic[TVarFCIdx][idx_i];

                        idx_i ++;
                        idx_o ++;
                     }}}

#                    else
                     Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#                    endif // #ifdef MHD ... else ...

                     Data1PG_FC_Ptr += PGSize3D_FC;
                  } // for (int v=0; v<NVarFC_Tot; v++)

               } // for (int Count=0; Count<TABLE_04(Side); Count++)
            } // if ( SibPID0 >= 0 )
         } // for (int Side=0; Side<NSide; Side++)


//       interpolation or boundary condition
//       --> apply interpolation AFTER copying data from all existing sibling patches to Data1PG_CC[]
//           since the divergence-free interpolation on the magnetic field requires these data
         for (int Side=0; Side<NSide; Side++)
         {
//          nothing to do if no ghost zone is required
            if ( GhostSize == 0 )   break;


            const int SibPID0 = Table_02( lv, PID0, Side );    // the 0th patch of the sibling patch group

//          (b2) if the target sibling patch does not exist --> interpolate from patches at level lv-1
            if ( SibPID0 == -1 )
            {
//             interpolation should never be applied to the base level
#              ifdef GAMER_DEBUG
               if ( lv == 0 )    Aux_Error( ERROR_INFO, "performing interpolation at the base level !!\n" );
#              endif


//             get the array size to store the interpolation result
               int FSize[3];
               for (int d=0; d<3; d++)  FSize[d] = TABLE_01( Side, 'x'+d, GhostSize_Padded, PS2, GhostSize_Padded );

               real *IntData_CC_Ptr = NULL;
               real *IntData_FC_Ptr = NULL;


//             (b2-1) determine the target PID at lv-1
               const int FaPID    = amr->patch[0][lv][PID0]->father;
               const int FaSibPID = amr->patch[0][lv-1][FaPID]->sibling[Side];

#              ifdef GAMER_DEBUG
               if ( FaSibPID < 0 )  Aux_Error( ERROR_INFO, "FaSibPID = %d < 0 (lv %d, PID0 %d, FaPID %d, sib %d) !!\n",
                                               FaSibPID, lv, PID0, FaPID, Side );
#              endif


//             (b2-2) collect the fine-grid magnetic field on the coarse-fine interpolation boundaries
//             for the divergence-preserving interpolation
               real *FInterface_Ptr[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
#              ifdef MHD
               if ( NVarFC_Tot > 0 )
               MHD_SetFInterface( FInterface_Data, FInterface_Ptr, Data1PG_FC, lv, PID0, Side, GhostSize,
                                  MagSg, MagSg_IntT, MagIntTime, MagWeighting, MagWeighting_IntT );
#              endif

               long TVarCCBuffer = TVarCC;
               int NVarCC_FluBuffer = NVarCC_Flu;
               int TVarCCIdxList_FluBuffer[NCOMP_TOTAL];
               memcpy( TVarCCIdxList_FluBuffer, TVarCCIdxList_Flu, sizeof(TVarCCIdxList_FluBuffer) );


#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
//             if fluid scheme is used on level lv-1 and wave scheme is used on level lv,
//             set target variable for lv-1 to dens and phase and later convert to real and imaginary parts
//             for preparing patch with real and imaginary part on level lv
               const bool ConvertWaveToFluid = (  amr->use_wave_flag[lv]  &&  !amr->use_wave_flag[lv-1]  &&  ( TVarCC & (_REAL | _IMAG) )  );
               if ( ConvertWaveToFluid ) {
//                preserve other flags while replacing _REAL and _IMAG with _DENS and _PHAS
                  TVarCCBuffer     = ( TVarCC & ~(_REAL | _IMAG) ) | (_DENS | _PHAS);
                  NVarCC_FluBuffer = 0;
                  for (int v=0; v<NCOMP_TOTAL; v++)
                     if ( TVarCCBuffer & (1L<<v) )
                        TVarCCIdxList_FluBuffer[ NVarCC_FluBuffer++ ] = v;
               }
#              endif

//             (b2-3) perform interpolation and store the results in IntData_CC[] and IntData_FC[]
               InterpolateGhostZone( lv-1, FaSibPID, IntData_CC, IntData_FC, IntData_CC_IntTime, Side, PrepTime, GhostSize,
                                     IntScheme_CC, IntScheme_FC, NTSib, TSib, TVarCCBuffer, NVarCC_Tot, NVarCC_FluBuffer,
                                     TVarCCIdxList_FluBuffer, NVarCC_Der, TVarCCList_Der, TVarFC, NVarFC_Tot, TVarFCIdxList,
                                     IntPhase, FluBC, PotBC, BC_Face, MinPres, MinTemp, MinEntr, DE_Consistency,
                                     (const real **)FInterface_Ptr );


//             (b2-3-1) convert density and phase in IntData_CC[] in hybrid scheme if we interpolate from fluid to wave scheme
//             --> do not convert NUseless-cell-wide useless data returned by InterpolateGhostZone()
               const int NUseless = GhostSize & 1;
               int loop[3], disp1[3], disp2[3];

               for (int d=0; d<3; d++)
               {
                  loop [d] = TABLE_01( Side, 'x'+d, GhostSize, PS2, GhostSize );
                  disp1[d] = TABLE_01( Side, 'x'+d, 0, GhostSize, GhostSize+PS2 );
                  disp2[d] = TABLE_01( Side, 'x'+d, NUseless, 0, 0 );
               }

               Data1PG_CC_Ptr = Data1PG_CC;
               IntData_CC_Ptr = IntData_CC;

#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
//             set target variables correctly to convert density and phase to real and imaginary parts
               if ( ConvertWaveToFluid ) {
//                density and phase --> real and imaginary parts
                  const int FSize3D_CC = FSize[0]*FSize[1]*FSize[2];
                  real Dens, Phase, Amp;

                  for (int k=0; k<loop[2]; k++) {  K2 = k + disp2[2];
                  for (int j=0; j<loop[1]; j++) {  J2 = j + disp2[1];
                                                   Idx2 = IDX321( disp2[0], J2, K2, FSize[0], FSize[1] );
                  for (int i=0; i<loop[0]; i++) {
                     Dens  = IntData_CC_Ptr[ Idx2 ];
                     Phase = IntData_CC_Ptr[ Idx2 + FSize3D_CC ];
                     Amp   = SQRT( Dens );

                     IntData_CC_Ptr[ Idx2               ] = Amp*COS( Phase );
                     IntData_CC_Ptr[ Idx2  + FSize3D_CC ] = Amp*SIN( Phase );
                     Idx2 ++;
                  }}}
               }
#              endif

//             (b2-4) copy cell-centered data from IntData_CC[] to Data1PG_CC[]
//             --> must get rid of NUseless-cell-wide useless data returned by InterpolateGhostZone()
               for (int v=0; v<NVarCC_Tot; v++)
               {
                  for (int k=0; k<loop[2]; k++) {  K = k + disp1[2];  K2 = k + disp2[2];
                  for (int j=0; j<loop[1]; j++) {  J = j + disp1[1];  J2 = j + disp2[1];
                                                   Idx1 = IDX321( disp1[0], J,  K,  PGSize1D_CC, PGSize1D_CC );
                                                   Idx2 = IDX321( disp2[0], J2, K2, FSize[0], FSize[1] );
                  for (int i=0; i<loop[0]; i++) {

                     Data1PG_CC_Ptr[ Idx1 ++ ] = IntData_CC_Ptr[ Idx2 ++ ];

                  }}}

                  Data1PG_CC_Ptr += PGSize3D_CC;
                  IntData_CC_Ptr += FSize[0]*FSize[1]*FSize[2];
               }


//             (b2-5) copy face-centered data from IntData_FC[] to Data1PG_FC[]
//             --> must get rid of NUseless-cell-wide useless data returned by InterpolateGhostZone()
               Data1PG_FC_Ptr = Data1PG_FC;
               IntData_FC_Ptr = IntData_FC;

               for (int v=0; v<NVarFC_Tot; v++)
               {
                  TVarFCIdx = TVarFCIdxList[v];

#                 ifdef MHD

//                set array indices
                  const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                                       ( TVarFCIdx == MAGY ) ? 1 :
                                       ( TVarFCIdx == MAGZ ) ? 2 : -1;
#                 ifdef GAMER_DEBUG
                  if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#                 endif

                  int ijk_s[3], ijk_e[3], size_i[3], size_o[3], disp_o[3], idx_i, idx_o;  // s/e=start/end; i/o=in/out

                  for (int d=0; d<3; d++)
                  {
                     ijk_s [d] = disp2[d];
                     ijk_e [d] = disp2[d] + loop[d];
                     size_i[d] = FSize[d];
                     size_o[d] = PGSize1D_CC;
                     disp_o[d] = disp1[d] - disp2[d];
                  }

//                avoid redundant assignment of the patch interface values
                  ijk_s [norm_dir] += TABLE_01( Side, 'x'+norm_dir, 0, 0, 1 );
                  ijk_e [norm_dir] += TABLE_01( Side, 'x'+norm_dir, 0, 1, 1 );
                  size_i[norm_dir] ++;
                  size_o[norm_dir] ++;


//                copy data
                  for (int k=ijk_s[2]; k<ijk_e[2]; k++)  {  K     = k + disp_o[2];
                  for (int j=ijk_s[1]; j<ijk_e[1]; j++)  {  J     = j + disp_o[1];
                                                            idx_o = IDX321( ijk_s[0]+disp_o[0], J, K, size_o[0], size_o[1] );
                                                            idx_i = IDX321( ijk_s[0],           j, k, size_i[0], size_i[1] );
                  for (int i=ijk_s[0]; i<ijk_e[0]; i++)  {

                     Data1PG_FC_Ptr[idx_o] = IntData_FC_Ptr[idx_i];

                     idx_i ++;
                     idx_o ++;
                  }}}

                  IntData_FC_Ptr += size_i[0]*size_i[1]*size_i[2];
#                 else
                  Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#                 endif // #ifdef MHD ... else ...

                  Data1PG_FC_Ptr += PGSize3D_FC;
               } // for (int v=0; v<NVarFC_Tot; v++)

            } // else if ( SibPID0 == -1 )


//          (b3) if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( SibPID0 <= SIB_OFFSET_NONPERIODIC )
            {
               Data1PG_CC_Ptr = Data1PG_CC;
               Data1PG_FC_Ptr = Data1PG_FC;

               for (int d=0; d<3; d++)
               {
                  BC_Idx_Start[d] = TABLE_01( Side, 'x'+d, 0, GhostSize, GhostSize+PS2 );
                  BC_Idx_End  [d] = TABLE_01( Side, 'x'+d, GhostSize, PS2, GhostSize ) + BC_Idx_Start[d] - 1;
               }

               BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID0;

#              ifdef GAMER_DEBUG
               if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
                  Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

               if ( FluBC[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
                  Aux_Error( ERROR_INFO, "FluBC == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID0 %d, PID0 %d, Side %d) !!\n",
                             BC_Sibling, BC_Face[BC_Sibling], SibPID0, PID0, Side );
#              endif

//             (b3-1) fluid B.C.
               if ( TVarCC & (_TOTAL|_DERIVED) )
               {
                  switch ( FluBC[ BC_Face[BC_Sibling] ] )
                  {
#                    if ( MODEL == HYDRO )
                     case BC_FLU_OUTFLOW:
                        Hydro_BoundaryCondition_Outflow   ( Data1PG_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu+NVarCC_Der, GhostSize,
                                                            PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End );
                     break;

                     case BC_FLU_REFLECTING:
                        Hydro_BoundaryCondition_Reflecting( Data1PG_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu,            GhostSize,
                                                            PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                            TVarCCIdxList_Flu, NVarCC_Der, TVarCCList_Der );
                     break;

                     case BC_FLU_DIODE:
                        Hydro_BoundaryCondition_Diode     ( Data1PG_CC_Ptr, BC_Face[BC_Sibling], NVarCC_Flu,            GhostSize,
                                                            PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                            TVarCCIdxList_Flu, NVarCC_Der, TVarCCList_Der );
                     break;
#                    endif

                     case BC_FLU_USER:
                        Flu_BoundaryCondition_User        ( Data1PG_CC_Ptr,                      NVarCC_Flu,            GhostSize,
                                                            PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                            TVarCCIdxList_Flu, PrepTime, dh, xyz0, TVarCC, lv );
                     break;

                     default:
                        Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
                  } // switch ( FluBC[ BC_Face[BC_Sibling] ] )

                  Data1PG_CC_Ptr += NVarCC_Flu*PGSize3D_CC;
               } // if ( TVarCC & (_TOTAL|_DERIVED) )


//             (b3-2) potential B.C.
#              ifdef GRAVITY
               if ( PrepPot )
               {
//                extrapolate potential
                  Poi_BoundaryCondition_Extrapolation( Data1PG_CC_Ptr, BC_Face[BC_Sibling], 1, GhostSize,
                                                       PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End );

                  Data1PG_CC_Ptr += 1*PGSize3D_CC;
               } // if ( PrepPot )
#              endif // #ifdef GRAVITY


//             (b3-3) face-centered variables B.C. (i.e., magnetic field)
#              ifdef MHD
               if ( NVarFC_Tot > 0 )
               {
                  real *MagDataPtr[NCOMP_MAG] = { NULL, NULL, NULL };
                  for (int v=0; v<NVarFC_Tot; v++)
                  {
                     MagDataPtr[ TVarFCIdxList[v] ] = Data1PG_FC_Ptr;
                     Data1PG_FC_Ptr                += PGSize3D_FC;
                  }

                  switch ( FluBC[ BC_Face[BC_Sibling] ] )
                  {
//                   input PGSize1D_CC instead PGSize1D_FC for the MHD boundary condition routines
                     case BC_FLU_OUTFLOW:
                        MHD_BoundaryCondition_Outflow   ( MagDataPtr, BC_Face[BC_Sibling], NVarFC_Tot, GhostSize,
                                                          PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                          TVarFCIdxList );
                     break;

                     case BC_FLU_REFLECTING:
                        MHD_BoundaryCondition_Reflecting( MagDataPtr, BC_Face[BC_Sibling], NVarFC_Tot, GhostSize,
                                                          PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                          TVarFCIdxList );
                     break;

                     case BC_FLU_DIODE:
                        MHD_BoundaryCondition_Diode     ( MagDataPtr, BC_Face[BC_Sibling], NVarFC_Tot, GhostSize,
                                                          PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                          TVarFCIdxList );
                     break;

                     case BC_FLU_USER:
                        MHD_BoundaryCondition_User      ( MagDataPtr, BC_Face[BC_Sibling], NVarFC_Tot,
                                                          PGSize1D_CC, PGSize1D_CC, PGSize1D_CC, BC_Idx_Start, BC_Idx_End,
                                                          TVarFCIdxList, PrepTime, dh, xyz0, lv );
                     break;

                     default:
                        Aux_Error( ERROR_INFO, "unsupported MHD B.C. (%d) !!\n", FluBC[ BC_Face[BC_Sibling] ] );
                  } // switch ( FluBC[ BC_Face[BC_Sibling] ] )
               } // if ( NVarFC_Tot > 0 )
#              endif // #ifdef MHD

            } // else if ( SibPID0 <= SIB_OFFSET_NONPERIODIC )


            else if ( SibPID0 < 0 )
               Aux_Error( ERROR_INFO, "SibPID0 == %d (lv %d, PID0 %d, Side %d) !!\n", SibPID0, lv, PID0, Side );

         } // for (int Side=0; Side<NSide; Side++)


//       c. deposit particle mass onto grids
// ------------------------------------------------------------------------------------------------------------
#        ifdef MASSIVE_PARTICLES
         if ( PrepParOnlyDens || PrepTotalDens )
         {
//          (c1) determine the array index for DENS
            real *ArrayDens = NULL;
            int   DensIdx   = -1;

            if ( PrepTotalDens )
            {
               for (int v=0; v<NVarCC_Flu; v++)
               {
                  if ( TVarCCIdxList_Flu[v] == DENS )
                  {
                     DensIdx = v;
                     break;
                  }
               }
            }

            else
               DensIdx = NVarCC_Tot - 1;  // store particle-only density in the last element

#           ifdef DEBUG_PARTICLE
            if ( DensIdx == -1 )    Aux_Error( ERROR_INFO, "DensIdx == -1 !!\n" );
#           endif

            ArrayDens = Data1PG_CC + DensIdx*PGSize3D_CC;

//          initialize density array as zero if there are no other density fields
            if ( PrepParOnlyDens )  for (int t=0; t<PGSize3D_CC; t++)   ArrayDens[t] = (real)0.0;


//          (c2) deposit particle mass in the central eight patches
            for (int LocalID=0; LocalID<8; LocalID++ )
            {
               const int PID = PID0 + LocalID;

//             skip patches without particles
               if ( amr->patch[0][lv][PID]->rho_ext == NULL  ||
                    amr->patch[0][lv][PID]->rho_ext[0][0][0] == RHO_EXT_NEED_INIT )    continue;

//             calculate the offset between rho_ext[] and ArrayDens[]
               const int Disp_i = TABLE_02( LocalID, 'x', GhostSize-RHOEXT_GHOST_SIZE, GhostSize+PS1-RHOEXT_GHOST_SIZE );
               const int Disp_j = TABLE_02( LocalID, 'y', GhostSize-RHOEXT_GHOST_SIZE, GhostSize+PS1-RHOEXT_GHOST_SIZE );
               const int Disp_k = TABLE_02( LocalID, 'z', GhostSize-RHOEXT_GHOST_SIZE, GhostSize+PS1-RHOEXT_GHOST_SIZE );

//             take care of the case with GhostSize < RHOEXT_GHOST_SIZE
               const int is = ( Disp_i >= 0 ) ? 0 : -Disp_i;
               const int js = ( Disp_j >= 0 ) ? 0 : -Disp_j;
               const int ks = ( Disp_k >= 0 ) ? 0 : -Disp_k;
               const int ie = ( RHOEXT_NXT+Disp_i <= PGSize1D_CC ) ? RHOEXT_NXT-1 : PGSize1D_CC-Disp_i-1;
               const int je = ( RHOEXT_NXT+Disp_j <= PGSize1D_CC ) ? RHOEXT_NXT-1 : PGSize1D_CC-Disp_j-1;
               const int ke = ( RHOEXT_NXT+Disp_k <= PGSize1D_CC ) ? RHOEXT_NXT-1 : PGSize1D_CC-Disp_k-1;

//             add particle density to the total density array
               for (int k=ks; k<=ke; k++)  {  K    = k + Disp_k;
               for (int j=js; j<=je; j++)  {  Idx1 = IDX321( is+Disp_i, j+Disp_j, K, PGSize1D_CC, PGSize1D_CC );
               for (int i=is; i<=ie; i++)  {

                  ArrayDens[ Idx1 ++ ] += amr->patch[0][lv][PID]->rho_ext[k][j][i];

               }}}
            } // for (int LocalID=0; LocalID<8; LocalID++ )


//          (c3) deposit particle mass from the sibling patches
            if ( amr->Par->GhostSize > 0  ||  GhostSize > 0  ||  amr->Par->PredictPos )
            for (int Side=0; Side<26; Side++)
            {
               const int SibPID0 = Table_02( lv, PID0, Side );    // the 0th patch of the sibling patch group

//             (c3-1) if the target sibling patch exists --> loop over nearby patches at the same level
               if ( SibPID0 >= 0 )
               {
//                LG = LargeGhost : for the case with GhostSize >= RHOEXT_NXT
                  const int is_LG = TABLE_01( Side, 'x', RHOEXT_NXT-GhostSize-RHOEXT_GHOST_SIZE, 0, 0 );
                  const int js_LG = TABLE_01( Side, 'y', RHOEXT_NXT-GhostSize-RHOEXT_GHOST_SIZE, 0, 0 );
                  const int ks_LG = TABLE_01( Side, 'z', RHOEXT_NXT-GhostSize-RHOEXT_GHOST_SIZE, 0, 0 );
                  const int ie_LG = TABLE_01( Side, 'x', RHOEXT_NXT-1, RHOEXT_NXT-1, GhostSize+RHOEXT_GHOST_SIZE-1 );
                  const int je_LG = TABLE_01( Side, 'y', RHOEXT_NXT-1, RHOEXT_NXT-1, GhostSize+RHOEXT_GHOST_SIZE-1 );
                  const int ke_LG = TABLE_01( Side, 'z', RHOEXT_NXT-1, RHOEXT_NXT-1, GhostSize+RHOEXT_GHOST_SIZE-1 );

//###OPTIMIZATION: simplify TABLE_03 and TABLE_04
                  for (int Count=0; Count<TABLE_04(Side); Count++)
                  {
                     const int SibPID = TABLE_03( Side, Count ) + SibPID0;

//                   skip patches without particles
                     if ( amr->patch[0][lv][SibPID]->rho_ext == NULL  ||
                          amr->patch[0][lv][SibPID]->rho_ext[0][0][0] == RHO_EXT_NEED_INIT )    continue;

//                   calculate the offset between rho_nxt and ArrayDens
                     const int Disp_i = Table_01( Side, 'x', Count, GhostSize-RHOEXT_GHOST_SIZE ) - is_LG;
                     const int Disp_j = Table_01( Side, 'y', Count, GhostSize-RHOEXT_GHOST_SIZE ) - js_LG;
                     const int Disp_k = Table_01( Side, 'z', Count, GhostSize-RHOEXT_GHOST_SIZE ) - ks_LG;

//                   take care of the case with GhostSize < RHOEXT_GHOST_SIZE
                     const int is = ( is_LG+Disp_i >= 0 ) ? is_LG : -Disp_i;
                     const int js = ( js_LG+Disp_j >= 0 ) ? js_LG : -Disp_j;
                     const int ks = ( ks_LG+Disp_k >= 0 ) ? ks_LG : -Disp_k;
                     const int ie = ( ie_LG+Disp_i < PGSize1D_CC ) ? ie_LG : PGSize1D_CC-Disp_i-1;
                     const int je = ( je_LG+Disp_j < PGSize1D_CC ) ? je_LG : PGSize1D_CC-Disp_j-1;
                     const int ke = ( ke_LG+Disp_k < PGSize1D_CC ) ? ke_LG : PGSize1D_CC-Disp_k-1;

//                   add particle density to the total density array
                     for (int k=ks; k<=ke; k++)  {  K    = k + Disp_k;
                     for (int j=js; j<=je; j++)  {  Idx1 = IDX321( is+Disp_i, j+Disp_j, K, PGSize1D_CC, PGSize1D_CC );
                     for (int i=is; i<=ie; i++)  {

                        ArrayDens[ Idx1 ++ ] += amr->patch[0][lv][SibPID]->rho_ext[k][j][i];

                     }}}
                  } // for (int Count=0; Count<TABLE_04(Side); Count++)
               } // if ( SibPID0 >= 0 )


//             (c3-2) if the target sibling patch does not exist --> find the father sibling patch
               else if ( SibPID0 == -1 )
               {
//                root-level patches should never have sibling ID == -1
#                 ifdef DEBUG_PARTICLE
                  if ( lv == 0 )    Aux_Error( ERROR_INFO, "SibPID0 == -1 at the root level !!\n" );
#                 endif

                  const int    FaPID    = amr->patch[0][lv][PID0]->father;
                  const int    FaSibPID = amr->patch[0][lv-1][FaPID]->sibling[Side];
                  const double EdgeL[3] = { amr->patch[0][lv][PID0]->EdgeL[0] - GhostSize*dh,
                                            amr->patch[0][lv][PID0]->EdgeL[1] - GhostSize*dh,
                                            amr->patch[0][lv][PID0]->EdgeL[2] - GhostSize*dh };
                  long      *ParList = NULL;
                  int        NPar;
                  bool       UseInputMassPos;
                  real_par **InputMassPos = NULL;
                  long_par **InputType    = NULL;

#                 ifdef DEBUG_PARTICLE
                  if ( FaSibPID < 0 )  Aux_Error( ERROR_INFO, "FaSibPID = %d < 0 (lv %d, PID0 %d, FaPID %d, sib %d) !!\n",
                                                  FaSibPID, lv, PID0, FaPID, Side );

                  if ( amr->patch[0][lv-1][FaSibPID]->son != -1 )
                     Aux_Error( ERROR_INFO, "FaSibPID->son = %d != -1 (lv %d, PID0 %d, FaPID %d, FaSibPID %d, sib %d) !!\n",
                                amr->patch[0][lv-1][FaSibPID]->son, lv, PID0, FaPID, FaSibPID, Side );
#                 endif

//                (c3-2-1) determine the number of particles and the particle list (father-sibling patch must NOT have son)
                  if ( FaSibPID < amr->NPatchComma[lv-1][1] )
                  {
                     NPar            = amr->patch[0][lv-1][FaSibPID]->NPar;
                     ParList         = amr->patch[0][lv-1][FaSibPID]->ParList;
                     UseInputMassPos = false;
                     InputMassPos    = NULL;
                     InputType       = NULL;
                  }

                  else
                  {
#                    ifdef LOAD_BALANCE
                     NPar            = amr->patch[0][lv-1][FaSibPID]->NPar_Copy;
                     ParList         = NULL;
                     UseInputMassPos = true;
                     InputMassPos    = amr->patch[0][lv-1][FaSibPID]->ParAttFlt_Copy;
                     InputType       = amr->patch[0][lv-1][FaSibPID]->ParAttInt_Copy;
#                    else
                     Aux_Error( ERROR_INFO, "FaSibPID (%d) is not a real patch (NReal %d) !!\n",
                                FaSibPID, amr->NPatchComma[lv-1][1] );
#                    endif
                  }

#                 ifdef DEBUG_PARTICLE
                  if ( NPar < 0 )
                     Aux_Error( ERROR_INFO, "NPar (%d) has not been calculated (lv %d, FaSibPID %d) !!\n",
                                NPar, lv-1, FaSibPID );

                  if ( NPar > 0 )
                  {
                     if ( UseInputMassPos )
                     {
                        if ( InputMassPos[PAR_MASS] == NULL  ||  InputMassPos[PAR_POSX] == NULL  ||
                             InputMassPos[PAR_POSY] == NULL  ||  InputMassPos[PAR_POSZ] == NULL )
                           Aux_Error( ERROR_INFO, "InputMassPos[0/1/2/3] == NULL for NPar (%d) > 0 (lv %d, FaSibPID %d) !!\n",
                                      NPar, lv-1, FaSibPID );
                        if ( InputType[PAR_TYPE] == NULL )
                           Aux_Error( ERROR_INFO, "InputType[0] == NULL for NPar (%d) > 0 (lv %d, FaSibPID %d) !!\n",
                                      NPar, lv-1, FaSibPID );
                     }

                     else if ( ParList == NULL )
                        Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, FaSibPID %d) !!\n",
                                   NPar, lv-1, FaSibPID );
                  }
#                 endif // #ifdef DEBUG_PARTICLE

//                (c3-2-2) deposit particle mass onto grids (**from particles in the father-sibling patches**)
//                         --> need to take care of the periodicity here since particles may have position far
//                             away from the target patch boundaries (i.e., patch->EdgeL/R)
//                             --> Periodic_Check, CheckFarAway_Yes
                  if ( NPar > 0 )
                  Par_MassAssignment( ParList, NPar, amr->Par->Interp, ArrayDens, PGSize1D_CC, EdgeL, dh,
                                      (amr->Par->PredictPos && !UseInputMassPos), PrepTime, InitZero_No,
                                      Periodic_Check, PeriodicNCell, UnitDens_No, CheckFarAway_Yes,
                                      UseInputMassPos, InputMassPos, InputType );
               } // else if ( SibPID0 == -1 )
            } // for (int Side=0; Side<26; Side++) if ( amr->Par->GhostSize > 0  ||  GhostSize > 0 )
         } // if ( PrepParOnlyDens || PrepTotalDens )
#        endif // #ifdef MASSIVE_PARTICLES


//       d. checks
// ------------------------------------------------------------------------------------------------------------
//       (d1) minimum density
//       --> note that it's unnecessary to check negative passive scalars thanks to the monotonic interpolation
#        if ( MODEL == HYDRO  ||  MODEL == ELBDM )
         if ( MinDens >= (real)0.0 )
         {
//          note that _DENS is turned on automatically for _TOTAL_DENS (and total density is stored in DENS)
            if ( TVarCC & _DENS )
            {
//             assuming that the order of variables stored in OutputCC is the same as patch->fluid[]
               const int DensIdx = DENS;
               real *ArrayDens = Data1PG_CC + DensIdx*PGSize3D_CC;

//             apply minimum density
//             --> note that for ELBDM it will result in dens != real^2 + imag^2
               for (int t=0; t<PGSize3D_CC; t++)   ArrayDens[t] = FMAX( ArrayDens[t], MinDens );
            }
         } // if ( MinDens >= (real)0.0 )
#        endif // #if ( MODEL == HYDRO  ||  MODEL == ELBDM )


//       (d2) divergence-free magnetic field
#        if ( defined MHD  &&  defined MHD_CHECK_DIV_B )
         if ( PrepMagFC )  MHD_CheckDivB( Data1PG_FC, GhostSize, DIV_B_TOLERANCE, lv, PrepTime );
#        endif


//       e. copy data from Data1PG_CC[] to OutputCC[]
// ------------------------------------------------------------------------------------------------------------
         if ( PrepUnit == UNIT_PATCH ) // separate the prepared patch group data into individual patches
         {
            const int PSize1D_CC = PS1 + 2*GhostSize;    // width of a single patch including ghost zones
            const int PSize3D_CC = CUBE(PSize1D_CC);
            const int PSize1D_FC = PSize1D_CC + 1;
            const int PSize3D_FC = PSize1D_FC*SQR(PSize1D_CC);

            real *OutputCC_Ptr = NULL;
            real *OutputFC_Ptr = NULL;


            for (int LocalID=0; LocalID<8; LocalID++)
            {
               const int N      = 8*TID + LocalID;
               const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1 );
               const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1 );
               const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1 );

//             cell-centered variables
               Data1PG_CC_Ptr = Data1PG_CC;
               OutputCC_Ptr   = OutputCC + N*NVarCC_Tot*PSize3D_CC;
               Idx2           = 0;

               for (int v=0; v<NVarCC_Tot; v++)
               {
                  for (int k=Disp_k; k<Disp_k+PSize1D_CC; k++)
                  for (int j=Disp_j; j<Disp_j+PSize1D_CC; j++)
                  {
                     Idx1 = IDX321( Disp_i, j, k, PGSize1D_CC, PGSize1D_CC );

                     for (int i=0; i<PSize1D_CC; i++)    OutputCC_Ptr[ Idx2 ++ ] = Data1PG_CC_Ptr[ Idx1 ++ ];
                  }

                  Data1PG_CC_Ptr += PGSize3D_CC;
               }


//             face-centered variables
               Data1PG_FC_Ptr = Data1PG_FC;
               OutputFC_Ptr   = OutputFC + N*NVarFC_Tot*PSize3D_FC;
               Idx2           = 0;

               for (int v=0; v<NVarFC_Tot; v++)
               {
                  TVarFCIdx = TVarFCIdxList[v];

#                 ifdef MHD

//                set array indices
                  int size_p[3], size_pg[3];    // p=patch, pg=patch_group


                  const int norm_dir = ( TVarFCIdx == MAGX ) ? 0 :
                                       ( TVarFCIdx == MAGY ) ? 1 :
                                       ( TVarFCIdx == MAGZ ) ? 2 : -1;
#                 ifdef GAMER_DEBUG
                  if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#                 endif

                  for (int d=0; d<3; d++)
                  {
                     if ( d == norm_dir )
                     {
                        size_p [d] = PSize1D_FC;
                        size_pg[d] = PGSize1D_FC;
                     }

                     else
                     {
                        size_p [d] = PSize1D_CC;
                        size_pg[d] = PGSize1D_CC;
                     }
                  }


//                copy data
                  for (int k=Disp_k; k<Disp_k+size_p[2]; k++)
                  for (int j=Disp_j; j<Disp_j+size_p[1]; j++)
                  {
                     Idx1 = IDX321( Disp_i, j, k, size_pg[0], size_pg[1] );

                     for (int i=0; i<size_p[0]; i++)  OutputFC_Ptr[ Idx2 ++ ] = Data1PG_FC_Ptr[ Idx1 ++ ];
                  }

#                 else
                  Aux_Error( ERROR_INFO, "currently only MHD supports face-centered variables !!" );
#                 endif // #ifdef MHD ... else ...

                  Data1PG_FC_Ptr += PGSize3D_FC;
               } // for (int v=0; v<NVarFC_Tot; v++)

            } // for (int LocalID=0; LocalID<8; LocalID++)
         } // if ( PrepUnit == UNIT_PATCH )

      } // for (int TID=0; TID<NPG; TID++)

      if ( PrepUnit == UNIT_PATCH )
      {
         delete [] Data1PG_CC;
         delete [] Data1PG_FC;
      }
      delete [] IntData_CC;
      delete [] IntData_FC;
      delete [] IntData_CC_IntTime;

#     ifdef MHD
      delete [] FInterface_Data;
#     endif

   } // end of OpenMP parallel region


// free memory
   for (int s=0; s<26; s++)   delete [] TSib[s];

} // FUNCTION : Prepare_PatchData



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTempIntPara
// Description :  Set the temporal interpolation parameters
//
// Note        :  1. Invoked by Prepare_PatchData() and InterpolateGhostZone()
//                2. Apply to fluid, potential, and magnetic field
//                3. Use call-by-reference to set the returned parameters
//
// Parameter   :  lv            : Target refinement level
//                Sg0           : Current Sg
//                PrepTime      : Target physical time to prepare data
//                Time0         : Physical time of Sg0
//                Time1         : Physical time of 1-Sg0
//                IntTime       : Whether or not the temporal interpolation is required
//                Sg            : Sg if temporal interpolation is not required
//                Sg_Int        : Sg if temporal interpolation is required
//                Weighting     : Weighting for the data stored in Sg
//                Weighting_Int : Weighting for the data stored in Sg_Int
//
// Return      :  IntTime, Sg, Sg_Int, Weighting, Weighting_Int
//-------------------------------------------------------------------------------------------------------
void SetTempIntPara( const int lv, const int Sg0, const double PrepTime, const double Time0, const double Time1,
                     bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT )
{

   const int Sg1 = 1 - Sg0;

   if      (  Mis_CompareRealValue( PrepTime, Time0, NULL, false )  )
   {
      IntTime        = false;
      Sg             = Sg0;
      Sg_IntT        = NULL_INT;
      Weighting      = NULL_REAL;
      Weighting_IntT = NULL_REAL;
   }

   else if (  Mis_CompareRealValue( PrepTime, Time1, NULL, false )  )
   {
      IntTime        = false;
      Sg             = Sg1;
      Sg_IntT        = NULL_INT;
      Weighting      = NULL_REAL;
      Weighting_IntT = NULL_REAL;
   }

   else
   {
//    print warning messages if temporal extrapolation is required
      const double TimeMin = MIN( Time0, Time1 );
      const double TimeMax = MAX( Time0, Time1 );

      if ( TimeMin < 0.0 )
         Aux_Error( ERROR_INFO, "TimeMin (%21.14e) < 0.0 ==> one of the data arrays has not been initialized !!\n", TimeMin );

      if (  ( PrepTime < TimeMin  ||  PrepTime-TimeMax >= 1.0e-12*TimeMax )  &&  MPI_Rank == 0  )
         Aux_Message( stderr, "WARNING : performing temporal extrapolation (lv %d, T_Prep %20.14e, T_Min %20.14e, T_Max %20.14e)\n",
                      lv, PrepTime, TimeMin, TimeMax );

      if ( OPT__INT_TIME )
      {
         IntTime        = true;
         Sg             = Sg0;
         Sg_IntT        = Sg1;
         Weighting      = ( +Time1 - PrepTime ) / ( Time1 - Time0 );
         Weighting_IntT = ( -Time0 + PrepTime ) / ( Time1 - Time0 );
      }

      else
      {
         IntTime        = false;
         Sg             = Sg0;
         Sg_IntT        = NULL_INT;
         Weighting      = NULL_REAL;
         Weighting_IntT = NULL_REAL;
      }
   } // Mis_CompareRealValue

} // FUNCTION : SetTempIntPara



// ============
// |  Tables  |
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  Return the displacement for Prepare_PatchData()
//
// Parameter   :  SibID     : Sibling index (0~25)
//                dim       : Target spatial direction (x/y/z)
//                Count     : Patch counter (0~3)
//                GhostSize : Number of ghost zones
//-------------------------------------------------------------------------------------------------------
int Table_01( const int SibID, const char dim, const int Count, const int GhostSize )
{

   switch ( dim )
   {
      case 'x':
      {
         switch ( SibID )
         {
            case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
               return 0;

            case 2: case 3:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 4: case 5:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 10: case 11: case 12: case 13:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
               return GhostSize + 2*PS1;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'x':


      case 'y':
      {
         switch ( SibID )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
               return 0;

            case 0: case 1:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 4: case 5:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 14: case 15: case 16: case 17:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
               return GhostSize + 2*PS1;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'y':


      case 'z':
      {
         switch ( SibID )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
               return 0;

            case 0: case 1:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 2: case 3:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 6: case 7: case 8: case 9:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PS1;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n",
                                "SibID", SibID, "Count", Count );
               }
            }

            case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
               return GhostSize + 2*PS1;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'z':


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c !!\n", "dim", dim );
         exit(1);
   } // switch ( dim )

   return NULL_INT;

} // FUNCTION : Table_01



//-------------------------------------------------------------------------------------------------------
// Function    :  Table_02
// Description :  Return the patch ID of the 0th patch (local ID = 0) of the sibling patch group
//
// Note        :  Work for Prepare_PatchData()
//
// Parameter   :  lv   : Target refinement level
//                PID  : Target patch ID to find its sibling patches
//                Side : Sibling index (0~25)
//-------------------------------------------------------------------------------------------------------
int Table_02( const int lv, const int PID, const int Side )
{

   int Sib;

   switch ( Side )
   {
      case 0:
         Sib = amr->patch[0][lv][PID  ]->sibling[0];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 1:
         Sib = amr->patch[0][lv][PID+1]->sibling[1];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 2:
         Sib = amr->patch[0][lv][PID  ]->sibling[2];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 3:
         Sib = amr->patch[0][lv][PID+2]->sibling[3];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 4:
         Sib = amr->patch[0][lv][PID  ]->sibling[4];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 5:
         Sib = amr->patch[0][lv][PID+3]->sibling[5];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 6:
         Sib = amr->patch[0][lv][PID  ]->sibling[6];
         if ( Sib >= 0 )  return Sib-4;
         else             return Sib;

      case 7:
         Sib = amr->patch[0][lv][PID+1]->sibling[7];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 8:
         Sib = amr->patch[0][lv][PID+2]->sibling[8];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 9:
         Sib = amr->patch[0][lv][PID+4]->sibling[9];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 10:
         Sib = amr->patch[0][lv][PID  ]->sibling[10];
         if ( Sib >= 0 )  return Sib-5;
         else             return Sib;

      case 11:
         Sib = amr->patch[0][lv][PID+2]->sibling[11];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 12:
         Sib = amr->patch[0][lv][PID+3]->sibling[12];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 13:
         Sib = amr->patch[0][lv][PID+5]->sibling[13];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 14:
         Sib = amr->patch[0][lv][PID  ]->sibling[14];
         if ( Sib >= 0 )  return Sib-6;
         else             return Sib;

      case 15:
         Sib = amr->patch[0][lv][PID+3]->sibling[15];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 16:
         Sib = amr->patch[0][lv][PID+1]->sibling[16];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 17:
         Sib = amr->patch[0][lv][PID+6]->sibling[17];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 18:
         Sib = amr->patch[0][lv][PID  ]->sibling[18];
         if ( Sib >= 0 )  return Sib-7;
         else             return Sib;

      case 19:
         Sib = amr->patch[0][lv][PID+1]->sibling[19];
         if ( Sib >= 0 )  return Sib-5;
         else             return Sib;

      case 20:
         Sib = amr->patch[0][lv][PID+2]->sibling[20];
         if ( Sib >= 0 )  return Sib-6;
         else             return Sib;

      case 21:
         Sib = amr->patch[0][lv][PID+4]->sibling[21];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 22:
         Sib = amr->patch[0][lv][PID+3]->sibling[22];
         if ( Sib >= 0 )  return Sib-4;
         else             return Sib;

      case 23:
         Sib = amr->patch[0][lv][PID+6]->sibling[23];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 24:
         Sib = amr->patch[0][lv][PID+5]->sibling[24];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 25:
         Sib = amr->patch[0][lv][PID+7]->sibling[25];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Side", Side );
         exit(1);

   } // switch ( Side )

   return NULL_INT;

} // FUNCTION : Table_02



//-------------------------------------------------------------------------------------------------------
// Function    :  Table_03
// Description :  Return the patch GID of the 0th patch (local ID = 0) of the sibling patch group
//
// Note        :  Work for Prepare_PatchData_HasWaveCounterpart()
//
// Parameter   :  lv   : Target refinement level
//                GID  : Target patch GID to find its sibling patches
//                Side : Sibling index (0~25)
//                Tree : Array of LB_GlobalPatches with GID information
//
// Return      :  Sib
//-------------------------------------------------------------------------------------------------------
long Table_03( const int lv, const long GID, const int Side, LB_GlobalTree& Tree )
{

   int Sib;

   switch ( Side )
   {
      case 0:
         Sib = Tree[GID  ].sibling[0];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 1:
         Sib = Tree[GID+1].sibling[1];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 2:
         Sib = Tree[GID  ].sibling[2];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 3:
         Sib = Tree[GID+2].sibling[3];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 4:
         Sib = Tree[GID  ].sibling[4];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 5:
         Sib = Tree[GID+3].sibling[5];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 6:
         Sib = Tree[GID  ].sibling[6];
         if ( Sib >= 0 )  return Sib-4;
         else             return Sib;

      case 7:
         Sib = Tree[GID+1].sibling[7];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 8:
         Sib = Tree[GID+2].sibling[8];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 9:
         Sib = Tree[GID+4].sibling[9];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 10:
         Sib = Tree[GID  ].sibling[10];
         if ( Sib >= 0 )  return Sib-5;
         else             return Sib;

      case 11:
         Sib = Tree[GID+2].sibling[11];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 12:
         Sib = Tree[GID+3].sibling[12];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 13:
         Sib = Tree[GID+5].sibling[13];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 14:
         Sib = Tree[GID  ].sibling[14];
         if ( Sib >= 0 )  return Sib-6;
         else             return Sib;

      case 15:
         Sib = Tree[GID+3].sibling[15];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 16:
         Sib = Tree[GID+1].sibling[16];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 17:
         Sib = Tree[GID+6].sibling[17];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 18:
         Sib = Tree[GID  ].sibling[18];
         if ( Sib >= 0 )  return Sib-7;
         else             return Sib;

      case 19:
         Sib = Tree[GID+1].sibling[19];
         if ( Sib >= 0 )  return Sib-5;
         else             return Sib;

      case 20:
         Sib = Tree[GID+2].sibling[20];
         if ( Sib >= 0 )  return Sib-6;
         else             return Sib;

      case 21:
         Sib = Tree[GID+4].sibling[21];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 22:
         Sib = Tree[GID+3].sibling[22];
         if ( Sib >= 0 )  return Sib-4;
         else             return Sib;

      case 23:
         Sib = Tree[GID+6].sibling[23];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 24:
         Sib = Tree[GID+5].sibling[24];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 25:
         Sib = Tree[GID+7].sibling[25];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Side", Side );
         exit(1);

   } // switch ( Side )

   return NULL_INT;

} // FUNCTION : Table_03



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTargetSibling
// Description :  Set the target sibling directions for preparing the ghost-zone data at the coarse-grid level
//
// Note        :  1. Work for Prepare_PatchData()
//                2. TSib needs to be deallocated manually
//                3. Sibling directions recorded in TSib must be in ascending numerical order for filling the
//                   non-periodic ghost-zone data in InterpolateGhostZone()
//                   --> Therefore, this function CANNOT be applied in LB_RecordExchangeDataPatchID(), in which
//                       case "SetTargetSibling" and "SetReceiveSibling" must be declared consistently
//
// Parameter   :  NTSib : Number of target sibling patches along different sibling directions
//                TSib  : Target sibling indices along different sibling directions
//-------------------------------------------------------------------------------------------------------
void SetTargetSibling( int NTSib[], int *TSib[] )
{

   for (int t= 0; t< 6; t++)  NTSib[t] = 17;
   for (int t= 6; t<18; t++)  NTSib[t] = 11;
   for (int t=18; t<26; t++)  NTSib[t] =  7;

   for (int s=0; s<26; s++)   TSib[s] = new int [ NTSib[s] ];

   TSib[ 0][ 0] =  1;
   TSib[ 0][ 1] =  2;
   TSib[ 0][ 2] =  3;
   TSib[ 0][ 3] =  4;
   TSib[ 0][ 4] =  5;
   TSib[ 0][ 5] =  7;
   TSib[ 0][ 6] =  9;
   TSib[ 0][ 7] = 10;
   TSib[ 0][ 8] = 11;
   TSib[ 0][ 9] = 12;
   TSib[ 0][10] = 13;
   TSib[ 0][11] = 16;
   TSib[ 0][12] = 17;
   TSib[ 0][13] = 19;
   TSib[ 0][14] = 21;
   TSib[ 0][15] = 23;
   TSib[ 0][16] = 25;

   TSib[ 1][ 0] =  0;
   TSib[ 1][ 1] =  2;
   TSib[ 1][ 2] =  3;
   TSib[ 1][ 3] =  4;
   TSib[ 1][ 4] =  5;
   TSib[ 1][ 5] =  6;
   TSib[ 1][ 6] =  8;
   TSib[ 1][ 7] = 10;
   TSib[ 1][ 8] = 11;
   TSib[ 1][ 9] = 12;
   TSib[ 1][10] = 13;
   TSib[ 1][11] = 14;
   TSib[ 1][12] = 15;
   TSib[ 1][13] = 18;
   TSib[ 1][14] = 20;
   TSib[ 1][15] = 22;
   TSib[ 1][16] = 24;

   TSib[ 2][ 0] =  0;
   TSib[ 2][ 1] =  1;
   TSib[ 2][ 2] =  3;
   TSib[ 2][ 3] =  4;
   TSib[ 2][ 4] =  5;
   TSib[ 2][ 5] =  8;
   TSib[ 2][ 6] =  9;
   TSib[ 2][ 7] = 11;
   TSib[ 2][ 8] = 13;
   TSib[ 2][ 9] = 14;
   TSib[ 2][10] = 15;
   TSib[ 2][11] = 16;
   TSib[ 2][12] = 17;
   TSib[ 2][13] = 20;
   TSib[ 2][14] = 21;
   TSib[ 2][15] = 24;
   TSib[ 2][16] = 25;

   TSib[ 3][ 0] =  0;
   TSib[ 3][ 1] =  1;
   TSib[ 3][ 2] =  2;
   TSib[ 3][ 3] =  4;
   TSib[ 3][ 4] =  5;
   TSib[ 3][ 5] =  6;
   TSib[ 3][ 6] =  7;
   TSib[ 3][ 7] = 10;
   TSib[ 3][ 8] = 12;
   TSib[ 3][ 9] = 14;
   TSib[ 3][10] = 15;
   TSib[ 3][11] = 16;
   TSib[ 3][12] = 17;
   TSib[ 3][13] = 18;
   TSib[ 3][14] = 19;
   TSib[ 3][15] = 22;
   TSib[ 3][16] = 23;

   TSib[ 4][ 0] =  0;
   TSib[ 4][ 1] =  1;
   TSib[ 4][ 2] =  2;
   TSib[ 4][ 3] =  3;
   TSib[ 4][ 4] =  5;
   TSib[ 4][ 5] =  6;
   TSib[ 4][ 6] =  7;
   TSib[ 4][ 7] =  8;
   TSib[ 4][ 8] =  9;
   TSib[ 4][ 9] = 12;
   TSib[ 4][10] = 13;
   TSib[ 4][11] = 15;
   TSib[ 4][12] = 17;
   TSib[ 4][13] = 22;
   TSib[ 4][14] = 23;
   TSib[ 4][15] = 24;
   TSib[ 4][16] = 25;

   TSib[ 5][ 0] =  0;
   TSib[ 5][ 1] =  1;
   TSib[ 5][ 2] =  2;
   TSib[ 5][ 3] =  3;
   TSib[ 5][ 4] =  4;
   TSib[ 5][ 5] =  6;
   TSib[ 5][ 6] =  7;
   TSib[ 5][ 7] =  8;
   TSib[ 5][ 8] =  9;
   TSib[ 5][ 9] = 10;
   TSib[ 5][10] = 11;
   TSib[ 5][11] = 14;
   TSib[ 5][12] = 16;
   TSib[ 5][13] = 18;
   TSib[ 5][14] = 19;
   TSib[ 5][15] = 20;
   TSib[ 5][16] = 21;

   TSib[ 6][ 0] =  1;
   TSib[ 6][ 1] =  3;
   TSib[ 6][ 2] =  4;
   TSib[ 6][ 3] =  5;
   TSib[ 6][ 4] =  9;
   TSib[ 6][ 5] = 11;
   TSib[ 6][ 6] = 13;
   TSib[ 6][ 7] = 16;
   TSib[ 6][ 8] = 17;
   TSib[ 6][ 9] = 21;
   TSib[ 6][10] = 25;

   TSib[ 7][ 0] =  0;
   TSib[ 7][ 1] =  3;
   TSib[ 7][ 2] =  4;
   TSib[ 7][ 3] =  5;
   TSib[ 7][ 4] =  8;
   TSib[ 7][ 5] = 11;
   TSib[ 7][ 6] = 13;
   TSib[ 7][ 7] = 14;
   TSib[ 7][ 8] = 15;
   TSib[ 7][ 9] = 20;
   TSib[ 7][10] = 24;

   TSib[ 8][ 0] =  1;
   TSib[ 8][ 1] =  2;
   TSib[ 8][ 2] =  4;
   TSib[ 8][ 3] =  5;
   TSib[ 8][ 4] =  7;
   TSib[ 8][ 5] = 10;
   TSib[ 8][ 6] = 12;
   TSib[ 8][ 7] = 16;
   TSib[ 8][ 8] = 17;
   TSib[ 8][ 9] = 19;
   TSib[ 8][10] = 23;

   TSib[ 9][ 0] =  0;
   TSib[ 9][ 1] =  2;
   TSib[ 9][ 2] =  4;
   TSib[ 9][ 3] =  5;
   TSib[ 9][ 4] =  6;
   TSib[ 9][ 5] = 10;
   TSib[ 9][ 6] = 12;
   TSib[ 9][ 7] = 14;
   TSib[ 9][ 8] = 15;
   TSib[ 9][ 9] = 18;
   TSib[ 9][10] = 22;

   TSib[10][ 0] =  0;
   TSib[10][ 1] =  1;
   TSib[10][ 2] =  3;
   TSib[10][ 3] =  5;
   TSib[10][ 4] =  8;
   TSib[10][ 5] =  9;
   TSib[10][ 6] = 13;
   TSib[10][ 7] = 15;
   TSib[10][ 8] = 17;
   TSib[10][ 9] = 24;
   TSib[10][10] = 25;

   TSib[11][ 0] =  0;
   TSib[11][ 1] =  1;
   TSib[11][ 2] =  2;
   TSib[11][ 3] =  5;
   TSib[11][ 4] =  6;
   TSib[11][ 5] =  7;
   TSib[11][ 6] = 12;
   TSib[11][ 7] = 15;
   TSib[11][ 8] = 17;
   TSib[11][ 9] = 22;
   TSib[11][10] = 23;

   TSib[12][ 0] =  0;
   TSib[12][ 1] =  1;
   TSib[12][ 2] =  3;
   TSib[12][ 3] =  4;
   TSib[12][ 4] =  8;
   TSib[12][ 5] =  9;
   TSib[12][ 6] = 11;
   TSib[12][ 7] = 14;
   TSib[12][ 8] = 16;
   TSib[12][ 9] = 20;
   TSib[12][10] = 21;

   TSib[13][ 0] =  0;
   TSib[13][ 1] =  1;
   TSib[13][ 2] =  2;
   TSib[13][ 3] =  4;
   TSib[13][ 4] =  6;
   TSib[13][ 5] =  7;
   TSib[13][ 6] = 10;
   TSib[13][ 7] = 14;
   TSib[13][ 8] = 16;
   TSib[13][ 9] = 18;
   TSib[13][10] = 19;

   TSib[14][ 0] =  1;
   TSib[14][ 1] =  2;
   TSib[14][ 2] =  3;
   TSib[14][ 3] =  5;
   TSib[14][ 4] =  7;
   TSib[14][ 5] =  9;
   TSib[14][ 6] = 12;
   TSib[14][ 7] = 13;
   TSib[14][ 8] = 17;
   TSib[14][ 9] = 23;
   TSib[14][10] = 25;

   TSib[15][ 0] =  1;
   TSib[15][ 1] =  2;
   TSib[15][ 2] =  3;
   TSib[15][ 3] =  4;
   TSib[15][ 4] =  7;
   TSib[15][ 5] =  9;
   TSib[15][ 6] = 10;
   TSib[15][ 7] = 11;
   TSib[15][ 8] = 16;
   TSib[15][ 9] = 19;
   TSib[15][10] = 21;

   TSib[16][ 0] =  0;
   TSib[16][ 1] =  2;
   TSib[16][ 2] =  3;
   TSib[16][ 3] =  5;
   TSib[16][ 4] =  6;
   TSib[16][ 5] =  8;
   TSib[16][ 6] = 12;
   TSib[16][ 7] = 13;
   TSib[16][ 8] = 15;
   TSib[16][ 9] = 22;
   TSib[16][10] = 24;

   TSib[17][ 0] =  0;
   TSib[17][ 1] =  2;
   TSib[17][ 2] =  3;
   TSib[17][ 3] =  4;
   TSib[17][ 4] =  6;
   TSib[17][ 5] =  8;
   TSib[17][ 6] = 10;
   TSib[17][ 7] = 11;
   TSib[17][ 8] = 14;
   TSib[17][ 9] = 18;
   TSib[17][10] = 20;

   TSib[18][ 0] =  1;
   TSib[18][ 1] =  3;
   TSib[18][ 2] =  5;
   TSib[18][ 3] =  9;
   TSib[18][ 4] = 13;
   TSib[18][ 5] = 17;
   TSib[18][ 6] = 25;

   TSib[19][ 0] =  0;
   TSib[19][ 1] =  3;
   TSib[19][ 2] =  5;
   TSib[19][ 3] =  8;
   TSib[19][ 4] = 13;
   TSib[19][ 5] = 15;
   TSib[19][ 6] = 24;

   TSib[20][ 0] =  1;
   TSib[20][ 1] =  2;
   TSib[20][ 2] =  5;
   TSib[20][ 3] =  7;
   TSib[20][ 4] = 12;
   TSib[20][ 5] = 17;
   TSib[20][ 6] = 23;

   TSib[21][ 0] =  0;
   TSib[21][ 1] =  2;
   TSib[21][ 2] =  5;
   TSib[21][ 3] =  6;
   TSib[21][ 4] = 12;
   TSib[21][ 5] = 15;
   TSib[21][ 6] = 22;

   TSib[22][ 0] =  1;
   TSib[22][ 1] =  3;
   TSib[22][ 2] =  4;
   TSib[22][ 3] =  9;
   TSib[22][ 4] = 11;
   TSib[22][ 5] = 16;
   TSib[22][ 6] = 21;

   TSib[23][ 0] =  0;
   TSib[23][ 1] =  3;
   TSib[23][ 2] =  4;
   TSib[23][ 3] =  8;
   TSib[23][ 4] = 11;
   TSib[23][ 5] = 14;
   TSib[23][ 6] = 20;

   TSib[24][ 0] =  1;
   TSib[24][ 1] =  2;
   TSib[24][ 2] =  4;
   TSib[24][ 3] =  7;
   TSib[24][ 4] = 10;
   TSib[24][ 5] = 16;
   TSib[24][ 6] = 19;

   TSib[25][ 0] =  0;
   TSib[25][ 1] =  2;
   TSib[25][ 2] =  4;
   TSib[25][ 3] =  6;
   TSib[25][ 4] = 10;
   TSib[25][ 5] = 14;
   TSib[25][ 6] = 18;

} // FUNCTION : SetTargetSibling



#ifdef MASSIVE_PARTICLES

// flag for checking whether Par_CollectParticle2OneLevel() has been called
// --> declared in Par_CollectParticle2OneLevel.cpp
extern bool Particle_Collected;

//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_PatchData_InitParticleDensityArray
// Description :  Initialize rho_ext[] used by Prepare_PatchData() for preparing particle density on grids
//
// Note        :  1. This function is currently called by Gra_AdvanceDt(), Output_DumpData_Total(),
//                   Output_DumpData_Total_HDF5(), Output_BasePowerSpectrum(), Aux_FindExtrema(),
//                   Aux_FindWeightedAverageCenter(), and Aux_ComputeProfile()
//                2. Apply to both real and buffer patches
//                3. For patches without any particle, this routine ensures either rho_ext == NULL or
//                   rho_ext[0][0][0] == RHO_EXT_NEED_INIT
//                4. rho_ext[] is always stored in Sg==0
//                5. Must call Par_CollectParticle2OneLevel() in advance to collect particles from higher levels
//                   and from other MPI ranks
//                6. Call Prepare_PatchData_FreeParticleDensityArray() to free resources after Prepare_PatchData()
//                7. rho_ext[] does NOT necessary store the correct particle mass density,
//                   especially for cells adjacent to the C-C and C-F boundaries
//                   --> It's because this routine doesn't consider particles in nearby patches, which are handled
//                       by Prepare_PatchData()
//                   --> So rho_ext[] should only be used within Prepare_PatchData()
//                8. Note that computing the particle density of a given patch requires computing rho_ext[] in
//                   nearby patches too. As a result, computing rho_ext[] in Prepare_PatchData(), as in the previous
//                   approach, would prevent from applying OpenMP parallelization *outside* Prepare_PatchData().
//                   This is the major reason why we switch to pre-compute rho_ext[] of all patches here instead of in
//                   Prepare_PatchData().
//                9. Only the master thread in OpenMP is allowed to call this routine
//                   --> Also note that adding "#pragma omp master" when calling this routine may significantly impact
//                       its performance since nested parallelism is disabled by default in Init_OpenMP()
//                   --> Therefore, in general, we do not recommend calling this routine within an OpenMP parallel region
//
// Parameter   :  lv       : Target refinement level
//                PrepTime : Target physical time for predicting particle position
//                           --> Must pass the same value to Prepare_PatchData()
//-------------------------------------------------------------------------------------------------------
void Prepare_PatchData_InitParticleDensityArray( const int lv, const double PrepTime )
{

// check
   if ( ! Particle_Collected )
      Aux_Error( ERROR_INFO, "must call Par_CollectParticle2OneLevel() in advance !!\n" );

#  ifdef OPENMP
   const int TID = omp_get_thread_num();
   if ( TID != 0 )   Aux_Error( ERROR_INFO, "only the master thread is allowed to call %s() (thread ID %d) !!\n", __FUNCTION__, TID );
#  endif


// constant settings used by Par_MassAssignment()
   const double dh              = amr->dh[lv];
   const bool   InitZero_Yes    = true;
   const bool   Periodic_No[3]  = { false, false, false };
   const bool   UnitDens_No     = false;
   const bool   CheckFarAway_No = false;

#  pragma omp parallel
   {
//    thread-private variables
      long      *ParList = NULL;
      int        NPar;
      double     EdgeL[3];
      bool       UseInputMassPos;
      real_par **InputMassPos = NULL;
      long_par **InputType    = NULL;

//    loop over all patches including buffer patches
#     pragma omp for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
      {
//       determine the number of particles and the particle list
         if ( amr->patch[0][lv][PID]->son == -1  &&  PID < amr->NPatchComma[lv][1] )
         {
            NPar            = amr->patch[0][lv][PID]->NPar;
            ParList         = amr->patch[0][lv][PID]->ParList;
            UseInputMassPos = false;
            InputMassPos    = NULL;
            InputType       = NULL;

#           ifdef DEBUG_PARTICLE
            if ( amr->patch[0][lv][PID]->NPar_Copy != -1 )
               Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                          lv, PID, amr->patch[0][lv][PID]->NPar_Copy );
#           endif
         }

         else
         {
//          note that amr->patch[0][lv][PID]->NPar>0 is still possible
            NPar            = amr->patch[0][lv][PID]->NPar_Copy;
#           ifdef LOAD_BALANCE
            ParList         = NULL;
            UseInputMassPos = true;
            InputMassPos    = amr->patch[0][lv][PID]->ParAttFlt_Copy;
            InputType       = amr->patch[0][lv][PID]->ParAttInt_Copy;
#           else
            ParList         = amr->patch[0][lv][PID]->ParList_Copy;
            UseInputMassPos = false;
            InputMassPos    = NULL;
            InputType       = NULL;
#           endif
         }

#        ifdef DEBUG_PARTICLE
         if ( NPar < 0 )
         {
            bool Pass = true;

//          exclude buffer patches not adjacent to a real patch
#           ifdef LOAD_BALANCE
            for (int p=0; p<amr->Par->R2B_Buff_NPatchTotal[lv][0]; p++) {
               if ( PID == amr->Par->R2B_Buff_PIDList[lv][0][p] ) {
                  Pass = false;
                  break;
               }
            }
#           else
            Pass = false;
#           endif

            if ( !Pass )   Aux_Error( ERROR_INFO, "NPar (%d) has not been calculated (lv %d, PID %d) !!\n",
                                      NPar, lv, PID );
         } // if ( NPar < 0 )
#        endif // #ifdef DEBUG_PARTICLE

         if ( NPar > 0 )
         {
#           ifdef DEBUG_PARTICLE
            if ( UseInputMassPos )
            {
               if ( InputMassPos[PAR_MASS] == NULL  ||  InputMassPos[PAR_POSX] == NULL  ||
                    InputMassPos[PAR_POSY] == NULL  ||  InputMassPos[PAR_POSZ] == NULL )
                  Aux_Error( ERROR_INFO, "InputMassPos[0/1/2/3] == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                             NPar, lv, PID );
               if ( InputType[PAR_TYPE] == NULL )
                  Aux_Error( ERROR_INFO, "InputType[0] == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                             NPar, lv, PID );
            }

            else if ( ParList == NULL )
               Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                          NPar, lv, PID );
#           endif

//          set the left edge of rho_ext[]
            const double RhoExtGhostPhySize = RHOEXT_GHOST_SIZE*dh;
            for (int d=0; d<3; d++)    EdgeL[d] = amr->patch[0][lv][PID]->EdgeL[d] - RhoExtGhostPhySize;

//          allocate rho_ext[]
            if ( amr->patch[0][lv][PID]->rho_ext == NULL )    amr->patch[0][lv][PID]->dnew();

//          deposit particle mass onto grids (**from particles in their home patch**)
//          --> don't have to worry about the periodicity (even for external buffer patches) here since
//              (1) all input particles should be close to the target patches even with position prediction
//              (2) amr->patch[0][lv][PID]->EdgeL/R already assumes periodicity for external buffer patches
//              --> Periodic_No, CheckFarAway_No
//          --> must initialize rho_ext[] as zero by InitZero_Yes
            Par_MassAssignment( ParList, NPar, amr->Par->Interp, amr->patch[0][lv][PID]->rho_ext[0][0], RHOEXT_NXT,
                                EdgeL, dh, (amr->Par->PredictPos && !UseInputMassPos), PrepTime, InitZero_Yes,
                                Periodic_No, NULL, UnitDens_No, CheckFarAway_No, UseInputMassPos, InputMassPos, InputType );
         } // if ( NPar > 0 )

         else
         {
//          set rho_ext[0][0][0] = RHO_EXT_NEED_INIT to indicate that it hasn't been set yet
            if ( amr->patch[0][lv][PID]->rho_ext != NULL )
               amr->patch[0][lv][PID]->rho_ext[0][0][0] = RHO_EXT_NEED_INIT;
         } // if ( NPar > 0 ) ... else ...
      } // for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   } // end of OpenMP parallel region

// set flag to true to indicate that this function has been called
   ParDensArray_Initialized = true;

} // FUNCTION : Prepare_PatchData_InitParticleDensityArray



//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_PatchData_FreeParticleDensityArray
// Description :  Free rho_ext[] allocated by Prepare_PatchData() temporarily for storing the partice mass density
//
// Note        :  1. This function is currently called by Gra_AdvanceDt(), Output_DumpData_Total(),
//                   Output_DumpData_Total_HDF5(), Output_BasePowerSpectrum(), Aux_FindExtrema(),
//                   Aux_FindWeightedAverageCenter(), and Aux_ComputeProfile()
//                2. Apply to buffer patches as well
//                3. Do not free memory if OPT__REUSE_MEMORY is on
//                   --> rho_ext[] will only be free'd together with other data arrays (e.g., fluid, pot)
//                4. Only the master thread in OpenMP is allowed to call this routine
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Prepare_PatchData_FreeParticleDensityArray( const int lv )
{

#  ifdef OPENMP
   const int TID = omp_get_thread_num();
   if ( TID != 0 )   Aux_Error( ERROR_INFO, "only the master thread is allowed to call %s() (thread ID %d) !!\n", __FUNCTION__, TID );
#  endif

// free memory for all patches (both real and buffer) if OPT__REUSE_MEMORY is off
   if ( ! OPT__REUSE_MEMORY ) {
      for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
         amr->patch[0][lv][PID]->ddelete();
   }

// set flag to false to indicate that Prepare_PatchData_InitParticleDensityArray() has not been called
   ParDensArray_Initialized = false;

} // FUNCTION : Prepare_PatchData_FreeParticleDensityArray
#endif // #ifdef MASSIVE_PARTICLES



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_SetFInterface
// Description :  Collect the fine-grid magnetic field on the coarse-fine interpolation boundaries
//                for the divergence-preserving interpolation
//
// Note        :  1. Collected B field (i.e., the FInt_Ptr[] array) will be passed to
//                   MHD_InterpolateBField() when invoking InterpolateGhostZone()
//                2. Currently the temporal interpolation, although supported, is not actually used
//                   --> MagIntTime is always false
//                3. Since the interpolated ghost zones must be an even number (i.e., GhostSize_Padded),
//                   one cannot use Data1PG_FC[] to get all the required fine-grid magnetic field when
//                   GhostSize is an odd number
//                   --> We only copy data from Data1PG_FC[] on the interfaces between the central
//                       patch group and it's sibling patches
//                   --> For the B field on the interfaces outside the central patch group, we recollect
//                       data from nearby fine patches
//                4. FInt_Data[] is preallocated to avoid frequent memory allocation/deallocation
//
// Parameter   :  FInt_Data         : Array to store the fine-grid magnetic field to be returned
//                FInt_Ptr          : Pointer arrays pointing to FInt_Data[]
//                Data1PG_FC        : Array storing the already prepared fine-grid magnetic field
//                lv                : Target refinement level
//                PID0              : 0th PID of the central patch group on lv
//                Side              : Target sibling direction relative to PID0
//                GhostSize         : Number of ghost zones to be prepared
//                MagSg             : Sandglass of the magnetic field
//                MagSg_IntT        : Sandglass for conducting temporal interpolation on the magnetic field
//                MagIntTime        : Whether or not to perform temporal interpolation on the magnetic field
//                MagWeighting      : Weighting of data stored in MagSg      when MagIntTime is on
//                MagWeighting_IntT : Weighting of data stored in MagSg_IntT when MagIntTime is on
//
// Return      :  FInt_Data, FInt_Ptr
//-------------------------------------------------------------------------------------------------------
void MHD_SetFInterface( real *FInt_Data, real *FInt_Ptr[6], const real *Data1PG_FC, const int lv, const int PID0,
                        const int Side, const int GhostSize, const int MagSg, const int MagSg_IntT,
                        const bool MagIntTime, const real MagWeighting, const real MagWeighting_IntT )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv == 0 )
      Aux_Error( ERROR_INFO, "%s should NOT be applied to the base level !!\n", __FUNCTION__ );
#  endif


   const int FaPID            = amr->patch[0][lv][PID0]->father;
   const int FaSibPID         = amr->patch[0][lv-1][FaPID]->sibling[Side];
   const int GhostSize_Padded = GhostSize + (GhostSize&1);
   const int PGSize1D_CC      = 2*( PS1 + GhostSize );
   const int PGSize1D_FC      = PGSize1D_CC + 1;
   const int PGSize3D_FC      = PGSize1D_FC*SQR(PGSize1D_CC);

   const real *Data1PG_FC_Ptr = NULL;
   int FaSibSibPID, norm_dir, sign, FInt_Side, Offset=0;
   int LCR[3], loop[3], disp_i[3], disp_o[3], size_i[3], size_o[3], ijk_i[3], ijk_o[3], idx_i, idx_o;    // i/o=in/out


// check
#  ifdef GAMER_DEBUG
   if ( FaSibPID < 0 )
      Aux_Error( ERROR_INFO, "FaSibPID = %d < 0 (lv %d) !!\n", FaSibPID, lv );

   if ( amr->patch[0][lv-1][FaSibPID]->son != -1 )
      Aux_Error( ERROR_INFO, "son = %d != -1 (lv %d, FaSibPID %d) !!\n",
                 amr->patch[0][lv-1][FaSibPID]->son, lv, FaSibPID );
#  endif


// iterate over the six faces of the target ghost-zone region
   for (int f=0; f<6; f++)
   {
      norm_dir    = f/2;   // [0,0,1,1,2,2]
      sign        = f&1;   // [0,1,0,1,0,1]
      FInt_Ptr[f] = NULL;  // initialize as NULL --> coarse-coarse interface

//    nothing to do on the left/right faces of left/right patches
      if (  TABLE_01( Side, 'x'+norm_dir, 0, NULL_INT, 1 ) == sign  )   continue;


      FaSibSibPID = amr->patch[0][lv-1][FaSibPID]->sibling[f];

#     ifdef GAMER_DEBUG
//    FaSibSibPID < -1 is possible due to non-periodic boundary conditions, but it
//    cannot be -1 due to the proper-nesting constraint
      if ( FaSibSibPID == -1 )   Aux_Error( ERROR_INFO, "FaSibSibPID == -1 !!\n" );
#     endif

//    check if the target face is a coarse-fine interface
//    --> note that [FaSibSibPID]->son can be < -1 since the son may live abroad
      if ( FaSibSibPID >= 0  &&  amr->patch[0][lv-1][FaSibSibPID]->son != -1 )
      {
//       1. get the sibling direction index relative to the central patch group
//          --> i.e., between FaPID and FaSibSibPID
         FInt_Side = -2;

//       target sibling->sibling patch group == central patch group
         if ( FaSibSibPID == FaPID )
            FInt_Side = -1;

         else
         {
            for (int s=0; s<26; s++)
            {
               if ( amr->patch[0][lv-1][FaPID]->sibling[s] == FaSibSibPID )
               {
                  const int LR[3] = {  TABLE_01( s, 'x', -1, 123, +1 ),
                                       TABLE_01( s, 'y', -1, 123, +1 ),
                                       TABLE_01( s, 'z', -1, 123, +1 )  };

//                this check is necessary when there are only two patches along any periodic direction
                  if (  TABLE_01( Side, 'x', +1, 456, -1 ) == LR[0]  ||
                        TABLE_01( Side, 'y', +1, 456, -1 ) == LR[1]  ||
                        TABLE_01( Side, 'z', +1, 456, -1 ) == LR[2]  ||
                        TABLE_01(    f, 'x', +1, 456, -1 ) == LR[0]  ||
                        TABLE_01(    f, 'y', +1, 456, -1 ) == LR[1]  ||
                        TABLE_01(    f, 'z', +1, 456, -1 ) == LR[2]    )
                     continue;

                  else
                  {
                     FInt_Side = s;
                     break;
                  }
               }
            }
         }

         if ( FInt_Side == -2 )  Aux_Error( ERROR_INFO, "cannot determine the sibling direction index !!\n" );


//       2. copy data to FInt_Data[] --> similar to step (b1) in Prepare_PatchData()
         FInt_Ptr[f] = FInt_Data + Offset;

//       2-1. copy data from the central patches
//            --> note that these data have already been stored in Data1PG_FC[]
         if ( FInt_Side == -1 )
         {
//          set array indices
            for (int d=0; d<3; d++)
            {
               if ( d == norm_dir )
               {
                  size_i[d] = PGSize1D_FC;
                  size_o[d] = 1;
                  disp_i[d] = GhostSize + (1-sign)*PS2;
               }

               else
               {
                  size_i[d] = PGSize1D_CC;
                  size_o[d] = PS2;
                  disp_i[d] = GhostSize;
               } // if ( d == norm_dir ) ... else ...
            }  // for (int d=0; d<3; d++)

//          copy data
            Data1PG_FC_Ptr = Data1PG_FC + norm_dir*PGSize3D_FC;
            idx_o = 0;

            for (int k=0; k<size_o[2]; k++) { ijk_i[2] = k + disp_i[2];
            for (int j=0; j<size_o[1]; j++) { ijk_i[1] = j + disp_i[1];
                                              idx_i = IDX321( disp_i[0], ijk_i[1], ijk_i[2], size_i[0], size_i[1] );
            for (int i=0; i<size_o[0]; i++) {

//             no temporal interpolation since it has already been applied to Data1PG_FC_Ptr[] if necessary
               FInt_Ptr[f][idx_o] = Data1PG_FC_Ptr[idx_i];

               idx_i ++;
               idx_o ++;
            }}} // i,j,k
         } // if ( FInt_Side == -1 )


//       2-2. copy data from the sibling patches
         else // FInt_Side = 0~25
         {
//          set array indices
            for (int d=0; d<3; d++)
            {
               LCR[d] = TABLE_01( FInt_Side, 'x'+d, -1, 0, 1 );

               if ( d == norm_dir )
               {
                  size_i[d] = PS1 + 1;
                  size_o[d] = 1;
                  loop  [d] = 1;
                  disp_i[d] = ( sign == 0 ) ? PS1 : 0;
               }

               else
               {
                  size_i[d] = PS1;

                  switch ( LCR[d] )
                  {
                     case -1:
                        size_o[d] = GhostSize_Padded;
                        loop  [d] = GhostSize_Padded;
                        disp_i[d] = PS1 - GhostSize_Padded;
                        break;

                     case 0:
                        size_o[d] = PS2;
                        loop  [d] = PS1;
                        disp_i[d] = 0;
                        break;

                     case 1:
                        size_o[d] = GhostSize_Padded;
                        loop  [d] = GhostSize_Padded;
                        disp_i[d] = 0;
                        break;

                     default:
                        Aux_Error( ERROR_INFO, "incorrect LCR[%d] = %d !!\n", d, LCR[d] );
                        exit( -1 );
                  }
               } // if ( d == norm_dir ) ... else ...
            }  // for (int d=0; d<3; d++)

            for (int Count=0; Count<TABLE_04(FInt_Side); Count++)
            {
//             note that we should not get SibPID0 by amr->patch[0][lv-1][FaSibSibPID]->son
//             since the latter can be < -1 for sons living abroad
               const int LocalID = TABLE_03( FInt_Side, Count );
               const int SibPID0 = Table_02( lv, PID0, FInt_Side );
               const int SibPID  = SibPID0 + LocalID;
#              ifdef GAMER_DEBUG
               if ( SibPID0 <= -1 )    Aux_Error( ERROR_INFO, "SibPID0 = %d <= -1 !!\n", SibPID0 );
#              endif

//             skip patches not adjacent to the target coarse-fine interface
               if (  TABLE_02( LocalID, 'x'+norm_dir, 0, 1 ) == sign  )    continue;

//             set array indices
               for (int d=0; d<3; d++)
               {
                  if (  d == norm_dir  ||  LCR[d] != 0 )    disp_o[d] = 0;
                  else                                      disp_o[d] = TABLE_02( LocalID, 'x'+d, 0, PS1 );
               }

//             copy data
               for (int k=0; k<loop[2]; k++) { ijk_i[2] = k + disp_i[2];   ijk_o[2] = k + disp_o[2];
               for (int j=0; j<loop[1]; j++) { ijk_i[1] = j + disp_i[1];   ijk_o[1] = j + disp_o[1];
                                               idx_i = IDX321( disp_i[0], ijk_i[1], ijk_i[2], size_i[0], size_i[1] );
                                               idx_o = IDX321( disp_o[0], ijk_o[1], ijk_o[2], size_o[0], size_o[1] );
               for (int i=0; i<loop[0]; i++) {

                  FInt_Ptr[f][idx_o] = amr->patch[MagSg][lv][SibPID]->magnetic[norm_dir][idx_i];

                  if ( MagIntTime ) // temporal interpolation
                  FInt_Ptr[f][idx_o] =
                     MagWeighting     *FInt_Ptr[f][idx_o]
                   + MagWeighting_IntT*amr->patch[MagSg_IntT][lv][SibPID]->magnetic[norm_dir][idx_i];

                  idx_i ++;
                  idx_o ++;
               }}} // i,j,k
            } // for (int Count=0; Count<TABLE_04(FInt_Side); Count++)
         } // if ( FInt_Side == -1 ) ... else ...

         Offset += size_o[0]*size_o[1]*size_o[2];
      } // check C-F interface
   } // for (int f=0; f<6; f++)

} // FUNCTION : MHD_SetFInterface



#ifdef MHD_CHECK_DIV_B
//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_CheckDivB
// Description :  Check if the divergence of the prepared magnetic field exceeds a given threshold
//
// Note        :  1. To enable this check, define MHD_CHECK_DIV_B and set the tolerance value
//                   in DIV_B_TOLERANCE manually on the top of this file
//                   --> Currently this check is disabled even when GAMER_DEBUG is on
//                2. This check may fail when both spatial and temporal interpolations are required
//                   to prepare the ghost-zone B field
//                   --> Because the area-averaged fine-grid B field on a coarse-fine interface !=
//                       **temporally interpolated** coarse-grid B field on the same interface
//                   --> Coarse cells adjacent to a C-F boundary is NOT divergence-free
//                       --> Since we will use the original fine-grid data on this C-F interface
//                   --> Moreover, the adopted interpolation scheme for B field is divergence-preserving
//                       instead of divergence-free
//
// Parameter   :  Data1PG_FC : Array storing the prepared B field to be checked
//                GhostSize  : Number of ghost zones
//                Tolerance  : Tolerance relative error in div(B)
//                             --> Relative error is defined as |div(B)|*dh/<B>
//                lv         : Target refinement level
//                PrepTime   : Target physical time
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void MHD_CheckDivB( const real *Data1PG_FC, const int GhostSize, const real Tolerance,
                    const int lv, const double PrepTime )
{

   const int  PGSize1D_CC = 2*( PS1 + GhostSize );
   const int  PGSize1D_FC = PGSize1D_CC + 1;
   const int  PGSize3D_FC = PGSize1D_FC*SQR(PGSize1D_CC);
   const int  didx_Bx     = 1;
   const int  didx_By     = PGSize1D_CC;
   const int  didx_Bz     = SQR(PGSize1D_CC);
   const real one_six     = 1.0/6.0;

   const real *Bx = Data1PG_FC + MAGX*PGSize3D_FC;
   const real *By = Data1PG_FC + MAGY*PGSize3D_FC;
   const real *Bz = Data1PG_FC + MAGZ*PGSize3D_FC;

   real BxL, BxR, ByL, ByR, BzL, BzR, AveB, DivB, DivB_max=-1.0;
   int  idx_Bx, idx_By, idx_Bz, i_max, j_max, k_max;


// find the cell with the maximum div(B)
   for (int k=0; k<PGSize1D_CC; k++ )
   for (int j=0; j<PGSize1D_CC; j++ )
   for (int i=0; i<PGSize1D_CC; i++ )
   {
      idx_Bx = IDX321( i, j, k, PGSize1D_FC, PGSize1D_CC );
      idx_By = IDX321( i, j, k, PGSize1D_CC, PGSize1D_FC );
      idx_Bz = IDX321( i, j, k, PGSize1D_CC, PGSize1D_CC );

      BxL = Bx[ idx_Bx           ];
      ByL = By[ idx_By           ];
      BzL = Bz[ idx_Bz           ];
      BxR = Bx[ idx_Bx + didx_Bx ];
      ByR = By[ idx_By + didx_By ];
      BzR = Bz[ idx_Bz + didx_Bz ];

      AveB = ( BxR + ByR + BzR + BxL + ByL + BzL ) * one_six;
      DivB = ( BxR + ByR + BzR - BxL - ByL - BzL );
      DivB = FABS( DivB );

      if ( DivB > DivB_max )
      {
         DivB_max = DivB;
         i_max    = i;
         j_max    = j;
         k_max    = k;
      }
   } // i,j,k


// warning if the maximum div(B) exceeds the tolerance value
   if ( DivB_max > Tolerance )
      Aux_Message( stderr, "WARNING : max div(B) = %20.14e at [%2d,%2d,%2d] (lv %d, PrepTime %20.14e, GhostSize %d) !!\n",
                   DivB_max, i_max-GhostSize, j_max-GhostSize, k_max-GhostSize, lv, PrepTime, GhostSize );

} // FUNCTION : MHD_CheckDivB
#endif // #ifdef MHD_CHECK_DIV_B

#endif // #ifdef MHD



#if ( ELBDM_SCHEME == ELBDM_HYBRID )
//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_PatchData_HasWaveCounterpart
// Description :  Prepare a uniform array indicating which cells have a refined wave counterpart including
//                ghost zones for the target patches or patch groups
//
// Note        :  1. Use "patch group" as the preparation unit
//                   --> The data of all patches within the same patch group will be prepared
//                2. Patches stored in PID0_List must be real patches (must NOT be buffer patches)
//
// Parameter   :  lv                   : Target refinement level
//                h_HasWaveCounterpart : Array to store the prepared booleans indicating which cells have
//                                       wave counterparts on refined levels
//                GhostSize            : Number of ghost zones to be prepared
//                NPG                  : Number of patch groups prepared at a time
//                PID0_List            : List recording the patch indices with LocalID==0 to be prepared
//                NSide                : Number of sibling directions to prepare data
//                                       --> NSIDE_00 (=  0) : do not prepare any sibling direction (equivalent to GhostSize=0)
//                                           NSIDE_06 (=  6) : prepare only sibling directions 0~5
//                                           NSIDE_26 (= 26) : prepare all sibling directions 0~25
//                GlobalTree           : LB_GlobalTree object for indexing patches with their GID
//-------------------------------------------------------------------------------------------------------
void Prepare_PatchData_HasWaveCounterpart( const int lv, bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                           const int GhostSize, const int NPG, const int *PID0_List,
                                           const NSide_t NSide, LB_GlobalTree* GlobalTree )
{

// nothing to do if there is no target patch group
   if ( NPG == 0 )   return;

   const double dh               = amr->dh[lv];
   const int    PGSize1D_CC      = 2*( PS1 + GhostSize );   // width of a single patch group including ghost zones

#  pragma omp parallel
   {
//    thread-private variables
      int    J, K, I2, J2, K2, Idx1, Idx2, PID0;
      long   GID0;

//    prepare eight nearby patches (one patch group) at a time
#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];
         GID0 = GlobalTree->PID2GID( PID0, lv );

//       1. reset h_HasWaveCounterpart
         for (int i=0; i<CUBE(HYB_NXT); i++)  {
               h_HasWaveCounterpart[TID][i] = false;
         }


//       2. fill out the central region of h_HasWaveCounterpart
         for (int LocalID=0; LocalID<8; LocalID++ )
         {
            const long GID    = GID0 + LocalID;
            const int  Disp_i = TABLE_02( LocalID, 'x', GhostSize, GhostSize+PS1 );
            const int  Disp_j = TABLE_02( LocalID, 'y', GhostSize, GhostSize+PS1 );
            const int  Disp_k = TABLE_02( LocalID, 'z', GhostSize, GhostSize+PS1 );

            for (int k=0; k<PS1; k++)  {  K    = k + Disp_k;
            for (int j=0; j<PS1; j++)  {  J    = j + Disp_j;
                                          Idx1 = IDX321( Disp_i, J, K, PGSize1D_CC, PGSize1D_CC );
            for (int i=0; i<PS1; i++)  {
               h_HasWaveCounterpart[TID][Idx1] = ELBDM_HasWaveCounterpart( i, j, k, GID, GID, *GlobalTree );
               Idx1 ++;
            }}}
         } // for (int LocalID=0; LocalID<8; LocalID++ )


//       3. fill out the ghost zones of h_HasWaveCounterpart by direct memory copy
         for (int Side=0; Side<NSide; Side++)
         {
//          nothing to do if no ghost zone is required
            if ( GhostSize == 0 )   break;

            const long SibGID0 = Table_03( lv, GID0, Side, *GlobalTree );  // the 0th patch of the sibling patch group

//          if the target sibling patch exists, just copy data from the nearby patches at the same level
//          --> if it does not exist, the respective cells do not have any children or lie outside the simulation domain
//          --> in both cases, we leave h_HasWaveCounterpart as false
            if ( SibGID0 >= 0 )
            {
               int loop[3], disp2[3];
               for (int d=0; d<3; d++)
               {
                  loop [d] = TABLE_01( Side, 'x'+d, GhostSize, PS1, GhostSize );
                  disp2[d] = TABLE_01( Side, 'x'+d, PS1-GhostSize, 0, 0 );
               }

               for (int Count=0; Count<TABLE_04(Side); Count++)
               {
                  const int  LocalID = TABLE_03( Side, Count );
                  const long SibGID  = SibGID0 + LocalID;

                  int disp[3];
                  for (int d=0; d<3; d++)    disp[d] = Table_01( Side, 'x'+d, Count, GhostSize );

                  for (int k=0; k<loop[2]; k++)  { K = k + disp[2];   K2 = k + disp2[2];
                  for (int j=0; j<loop[1]; j++)  { J = j + disp[1];   J2 = j + disp2[1];
                                                   Idx1 = IDX321( disp[0], J, K, PGSize1D_CC, PGSize1D_CC );
                  for (I2=disp2[0]; I2<disp2[0]+loop[0]; I2++) {
                     h_HasWaveCounterpart[TID][Idx1] = ELBDM_HasWaveCounterpart( I2, J2, K2, SibGID, SibGID, *GlobalTree );
                     Idx1 ++;
                  }}}

               } // for (int Count=0; Count<TABLE_04(Side); Count++)
            } // if ( SibPID0 >= 0 )
         } // for (int Side=0; Side<NSide; Side++)
      } // for (int TID=0; TID<NPG; TID++)
   } // end of OpenMP parallel region
} // FUNCTION : Prepare_PatchData_HasWaveCounterpart
#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
