#include "GAMER.h"

void Flag_Grandson( const int lv, const int PID, const int LocalID );
void Prepare_for_Lohner( const OptLohnerForm_t Form, const real *Var1D, real *Ave1D, real *Slope1D, const int NVar );

#if ( MODEL == ELBDM )
void Prepare_for_Spectral_Criterion( const real *Var1D, real& Cond1D );
#endif

static bool Flag_IterateCells( const int Mode, const int lv, const int PID, const real dv,
                               const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real MagCC[][PS1][PS1][PS1],
                               const real Vel[][PS1][PS1][PS1], const real Pres[][PS1][PS1], const real Lrtz[][PS1][PS1],
                               const real LCool[][PS1][PS1], const real (*Cs2)[PS1][PS1],
                               const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                               const real ParCount[][PS1][PS1], const real ParDens[][PS1][PS1],
                               const real *Interf_Var, const real Spectral_Cond,
                               const int SibID_Array[][3][3], const int FlagBuf, const real JeansCoeff_Factor );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Real
// Description :  Flag the real patches at level "lv" according to the given refinement criteria
//
// Note        :  1. Buffer patches are flagged by Flag_Buffer()
//                   --> But they can still be flagged by this function due to the non-zero
//                   (FLAG_BUFFER_SIZE, FLAG_BUFFER_SIZE_MAXM1_LV, FLAG_BUFFER_SIZE_MAXM2_LV) and the grandson check
//                3. To add new refinement criteria, please edit Flag_Check()
//                4. Prepare_for_Lohner() is defined in Flag_Lohner.cpp
//
// Parameter   :  lv        : Target refinement level to be flagged
//                UseLBFunc : Use the load-balance alternative functions for the grandson check and exchanging
//                            the buffer flags (useless if LOAD_BALANCE is off)
//                            --> USELB_YES : use the load-balance alternative functions
//                                USELB_NO  : do not use the load-balance alternative functions
//                                            --> useful for LOAD_BALANCE during the initialization
//-------------------------------------------------------------------------------------------------------
void Flag_Real( const int lv, const UseLBFunc_t UseLBFunc )
{

// check
   if ( lv == NLEVEL-1 )
      Aux_Error( ERROR_INFO, "function <%s> should NOT be applied to the finest level\" !!\n", __FUNCTION__ );

// user-specified operations before flagging
   if ( Flag_UserWorkBeforeFlag_Ptr != NULL )   Flag_UserWorkBeforeFlag_Ptr( Time[lv], lv );


// initialize all flags as false
#  pragma omp parallel for schedule( static )
   for (int PID=0; PID<amr->num[lv]; PID++)  amr->patch[0][lv][PID]->flag = false;


   const int SibID_Array[3][3][3]       = {  { {18, 10, 19}, {14,   4, 16}, {20, 11, 21} },
                                             { { 6,  2,  7}, { 0, 999,  1}, { 8,  3,  9} },
                                             { {22, 12, 23}, {15,   5, 17}, {24, 13, 25} }  };    // sibling indices
   const int  FlagBuf                   = ( lv == MAX_LEVEL-1 ) ? FLAG_BUFFER_SIZE_MAXM1_LV :
                                          ( lv == MAX_LEVEL-2 ) ? FLAG_BUFFER_SIZE_MAXM2_LV :
                                                                  FLAG_BUFFER_SIZE;
   const real dv                        = CUBE( amr->dh[lv] );
   const bool IntPhase_No               = false;                     // for invoking Prepare_PatchData()
   const bool DE_Consistency_No         = false;                     // for invoking Prepare_PatchData()
   const int  NPG                       = 1;                         // for invoking Prepare_PatchData()

// Lohner criterion
   const int  Lohner_NGhost             = 2;                         // number of ghost cells for the Lohner error estimator
   const int  Lohner_NCell              = PS1 + 2*Lohner_NGhost;     // size of the variable array for Lohner
   const int  Lohner_NAve               = Lohner_NCell - 2;          // size of the average array for Lohner
   const int  Lohner_NSlope             = Lohner_NAve;               // size of the slope array for Lohner
   const IntScheme_t Lohner_IntScheme   = INT_MINMOD1D;              // interpolation scheme for Lohner

#  if ( MODEL == ELBDM )
// interference criterion
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   const int  Interf_NGhost             = 1;                         // number of ghost cells for the interference criterion
   const int  Interf_NCell              = PS1 + 2*Interf_NGhost;     // size of the input array
   const IntScheme_t Interf_IntScheme   = INT_CQUAD;                 // interpolation scheme
#  endif

// spectral refinement criterion
   const int  Spectral_NGhost           = 1;                         // number of ghost cells for the spectral refinement criterion
   const int  Spectral_NCell            = PS2 + 2*Spectral_NGhost;   // size of the input array
   const IntScheme_t Spectral_IntScheme = INT_CQUAD;                 // interpolation scheme
#  endif // # if ( MODEL == ELBDM )

#  if ( MODEL == HYDRO  &&  defined GRAVITY )
#  ifdef COMOVING
   const real JeansCoeff_Factor       = M_PI/( SQR(FlagTable_Jeans[lv])*NEWTON_G*Time[lv] ); // flag if dh^2 > JeansCoeff_Factor*Gamma*Pres/Dens^2
#  else
   const real JeansCoeff_Factor       = M_PI/( SQR(FlagTable_Jeans[lv])*NEWTON_G          );
#  endif
#  else
   const real JeansCoeff_Factor       = NULL_REAL;
#  endif // #if ( MODEL == HYDRO  &&  defined GRAVITY ) ... else ...
#  ifndef GRAVITY
   const OptPotBC_t OPT__BC_POT       = BC_POT_NONE;
#  endif
#  ifdef PARTICLE
   const bool PredictPos_No           = false;                 // used by Par_MassAssignment()
   const bool InitZero_Yes            = true;
   const bool Periodic_No[3]          = { false, false, false };
   const bool UnitDens_Yes            = true;
   const bool UnitDens_No             = false;
   const bool CheckFarAway_No         = false;
   const bool SibBufPatch_No          = false;
   const bool FaSibBufPatch_No        = false;
   const bool JustCountNPar_Yes       = true;
   const bool JustCountNPar_No        = false;
   const bool TimingSendPar_No        = false;
#  endif

// flag-free region used by OPT__NO_FLAG_NEAR_BOUNDARY
// --> must be set precisely on the target level for OPT__UM_IC_DOWNGRADE
   const int  NoRefineBoundaryRegion  = ( OPT__NO_FLAG_NEAR_BOUNDARY ) ? PS1*( 1<<(NLEVEL-lv) )*( (1<<lv)-1 ) : NULL_INT;


// set the variables for the Lohner's error estimator and interference criterion
   int  Lohner_NVar=0, Lohner_Stride=0;
   long Lohner_TVar=0;
   int  Interf_NVar=0, Interf_Stride=0;
   int  Spectral_NVar=0;
   real MinDens=-1.0, MinPres=-1.0, MinTemp=-1.0, MinEntr=-1.0;  // default is to disable all floors

#  if   ( MODEL == HYDRO )
   if ( OPT__FLAG_LOHNER_DENS )  {  Lohner_NVar++;   Lohner_TVar |= _DENS;   MinDens = MIN_DENS;  }
   if ( OPT__FLAG_LOHNER_ENGY )  {  Lohner_NVar++;   Lohner_TVar |= _ENGY;                        }
   if ( OPT__FLAG_LOHNER_PRES )  {  Lohner_NVar++;   Lohner_TVar |= _PRES;   MinPres = MIN_PRES;  }
   if ( OPT__FLAG_LOHNER_TEMP )  {  Lohner_NVar++;   Lohner_TVar |= _TEMP;   MinTemp = MIN_TEMP;  }
   if ( OPT__FLAG_LOHNER_ENTR )  {  Lohner_NVar++;   Lohner_TVar |= _ENTR;   MinEntr = MIN_ENTR;  }
#  ifdef COSMIC_RAY
   if ( OPT__FLAG_LOHNER_CRAY )  {  Lohner_NVar++;   Lohner_TVar |= _CRAY;                        }
#  endif

#  elif ( MODEL == ELBDM )
   if ( OPT__FLAG_LOHNER_DENS )
   {
//    use Lohner criterion on wave levels
      if ( amr->use_wave_flag[lv] ) {
         Lohner_NVar = 2;
         Lohner_TVar = _REAL | _IMAG;
//    do not use Lohner criterion on fluid levels
      } else {
         Lohner_NVar = 0;
      }
   }

   if ( OPT__FLAG_SPECTRAL )
   {
//    use spectral criterion on wave levels
      if ( amr->use_wave_flag[lv] ) {
         Spectral_NVar = 2;
//    do not use spectral criterion on fluid levels
      } else {
         Spectral_NVar = 0;
      }
   }

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( OPT__FLAG_INTERFERENCE )
   {
//    use interference criterion on fluid levels
      if ( !amr->use_wave_flag[lv] ) {
         Interf_NVar = 2;
//    do not use interference criterion on wave levels
      } else {
         Interf_NVar = 0;
      }
      Interf_Stride = Interf_NVar*CUBE(Interf_NCell); // stride of array for one interference criterion patch
   }
#  endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

   Lohner_Stride = Lohner_NVar*Lohner_NCell*Lohner_NCell*Lohner_NCell;  // stride of array for one Lohner patch


// 0. collect particles to **real** patches at lv
#  ifdef PARTICLE
   long ColParFltAtt=_NONE, ColParIntAtt=_NONE;

   if ( OPT__FLAG_NPAR_CELL  ||  OPT__FLAG_PAR_MASS_CELL ) {
      ColParFltAtt |= _PAR_MASS | _PAR_POSX | _PAR_POSY | _PAR_POSZ;
      ColParIntAtt |= _PAR_TYPE;
   }

   if ( OPT__FLAG_PAR_TARGET != FLAG_PAR_NONE ) {
      ColParIntAtt |= _PAR_FLAG;
   }

   if ( ColParFltAtt != _NONE  ||  ColParIntAtt != _NONE )
      Par_CollectParticle2OneLevel( lv, ColParFltAtt, ColParIntAtt, PredictPos_No,
                                    NULL_REAL, SibBufPatch_No, FaSibBufPatch_No, JustCountNPar_No,
                                    TimingSendPar_No );

// Par_CollectParticle2OneLevel() with JustCountNPar_No also sets NPar_Copy for each patch
// --> call Par_CollectParticle2OneLevel() with JustCountNPar_Yes only when no particle attributes need to be collected
   else if ( OPT__FLAG_NPAR_PATCH != 0 )
      Par_CollectParticle2OneLevel( lv, _NONE, _NONE, PredictPos_No,
                                    NULL_REAL, SibBufPatch_No, FaSibBufPatch_No, JustCountNPar_Yes,
                                    TimingSendPar_No );
#  endif


//###ISSUE: use atomic ??
#  pragma omp parallel
   {
      const real (*Fluid)[PS1][PS1][PS1] = NULL;
      real (*Pot )[PS1][PS1]             = NULL;
      real (*MagCC)[PS1][PS1][PS1]       = NULL;
      real (*Vel)[PS1][PS1][PS1]         = NULL;
      real (*Pres)[PS1][PS1]             = NULL;
      real (*Cs2)[PS1][PS1]              = NULL;
      real (*Lrtz)[PS1][PS1]             = NULL;
      real (*LCool)[PS1][PS1]            = NULL;
      real (*ParCount)[PS1][PS1]         = NULL;   // declare as **real** to be consistent with Par_MassAssignment()
      real (*ParDens )[PS1][PS1]         = NULL;
      real *Lohner_Var                   = NULL;   // array storing the variables for Lohner
      real *Lohner_Ave                   = NULL;   // array storing the averages of Lohner_Var for Lohner
      real *Lohner_Slope                 = NULL;   // array storing the slopes of Lohner_Var for Lohner
      real *Interf_Var                   = NULL;   // array storing the density and phase for the interference criterion
      real *Spectral_Var                 = NULL;   // array storing a patch group of real and imaginary parts for the spectral criterion
      real  Spectral_Cond                = 0.0;    // variable storing the magnitude of the largest coefficient for the spectral criterion
      real *Grackle_TCool                = NULL;   // array storing a patch group of grackle cooling time

#     if ( MODEL == HYDRO )
      bool NeedPres = false;
      bool NeedCs2  = false;
      if ( OPT__FLAG_PRES_GRADIENT )   NeedPres = true;
#     ifdef GRAVITY
      if ( OPT__FLAG_JEANS )           NeedPres = true;
      if ( OPT__FLAG_JEANS )           NeedCs2  = true;
#     endif
#     ifdef SUPPORT_GRACKLE
      if ( OPT__FLAG_COOLING_LEN )     NeedCs2  = true;
#     endif
      if ( NeedCs2 )                   NeedPres = true;

#     ifdef MHD
      if ( OPT__FLAG_CURRENT || NeedPres )   MagCC    = new real [3][PS1][PS1][PS1];
#     endif
#     ifdef SRHD
      if ( OPT__FLAG_LRTZ_GRADIENT )         Lrtz     = new real    [PS1][PS1][PS1];
#     endif
#     ifdef SUPPORT_GRACKLE
      if ( OPT__FLAG_COOLING_LEN )           LCool    = new real    [PS1][PS1][PS1];
#     endif
      if ( OPT__FLAG_VORTICITY )             Vel      = new real [3][PS1][PS1][PS1];
      if ( NeedPres )                        Pres     = new real    [PS1][PS1][PS1];
      if ( NeedCs2 )                         Cs2      = new real    [PS1][PS1][PS1];
#     endif // HYDRO

#     ifdef PARTICLE
      if ( OPT__FLAG_NPAR_CELL )             ParCount = new real    [PS1][PS1][PS1];
      if ( OPT__FLAG_PAR_MASS_CELL )         ParDens  = new real    [PS1][PS1][PS1];
#     endif

#     if ( MODEL == ELBDM )
      if ( Spectral_NVar > 0 )
         Spectral_Var = new real [ Spectral_NVar*CUBE(Spectral_NCell) ];   // prepare one patch group
#     endif

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( Interf_NVar > 0 )
         Interf_Var   = new real [ 8*Interf_NVar*CUBE(Interf_NCell) ];  // 8: number of local patches
#     endif

      if ( Lohner_NVar > 0 )
      {
         Lohner_Var   = new real [ 8*Lohner_NVar*CUBE(Lohner_NCell)  ]; // 8: number of local patches
         Lohner_Ave   = new real [ 3*Lohner_NVar*CUBE(Lohner_NAve)   ]; // 3: X/Y/Z of 1 patch
         Lohner_Slope = new real [ 3*Lohner_NVar*CUBE(Lohner_NSlope) ]; // 3: X/Y/Z of 1 patch
      }

#     ifdef SUPPORT_GRACKLE
      if ( OPT__FLAG_COOLING_LEN )
         Grackle_TCool = new real [ CUBE(PS2) ];   // prepare one patch group
#     endif


//    loop over all REAL patches (the buffer patches will be flagged only due to the FlagBuf
//    extension or the grandson check)
#     pragma omp for schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
//       1. precompute various per-patch-group quantities required by the selected refinement flag checks
//       1-1. prepare the ghost-zone data for Lohner
         if ( Lohner_NVar > 0 )
            Prepare_PatchData( lv, Time[lv], Lohner_Var, NULL, Lohner_NGhost, NPG, &PID0, Lohner_TVar, _NONE,
                               Lohner_IntScheme, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No, OPT__BC_FLU, OPT__BC_POT,
                               MinDens, MinPres, MinTemp, MinEntr, DE_Consistency_No );


//       1-2. prepare the ghost-zone data for interference criterion
#        if ( MODEL == ELBDM )
         if ( Spectral_NVar > 0 )
         {
            Prepare_PatchData( lv, Time[lv], Spectral_Var, NULL, Spectral_NGhost, NPG, &PID0, _REAL|_IMAG, _NONE,
                               Spectral_IntScheme, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No, OPT__BC_FLU, OPT__BC_POT,
                               MinDens, MinPres, MinTemp, MinEntr, DE_Consistency_No );

//          evaluate the spectral refinement criterion
            Prepare_for_Spectral_Criterion( Spectral_Var, Spectral_Cond );
         }
#        endif

#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         if ( Interf_NVar > 0 )
            Prepare_PatchData( lv, Time[lv], Interf_Var, NULL, Interf_NGhost, NPG, &PID0, _DENS|_PHAS, _NONE,
                               Interf_IntScheme, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No, OPT__BC_FLU, OPT__BC_POT,
                               MinDens, MinPres, MinTemp, MinEntr, DE_Consistency_No );
#        endif


//       1-3. cooling length
#        ifdef SUPPORT_GRACKLE
         if ( OPT__FLAG_COOLING_LEN )
         {
//          Grackle_Calculate() involves shared variables and is not thread-safe
//          --> use critical directive to avoid thread racing
#           pragma omp critical
            Grackle_Calculate( Grackle_TCool, _GRACKLE_TCOOL, lv, NPG, &PID0 );
         }
#        endif


//       loop over all local patches within the same patch group
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int PID = PID0 + LocalID;

//          2. skip this patch if any refinement pre-check fails
            if (  ! Flag_Precheck( lv, PID, NoRefineBoundaryRegion )  )    continue;


            Fluid = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
#           ifdef GRAVITY
            Pot   = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;
#           endif


//          3. precompute various per-patch quantities required by the selected refinement flag checks
#           if ( MODEL == HYDRO )
#           ifdef MHD
//          3-1. evaluate cell-centered B field
            if ( OPT__FLAG_CURRENT || NeedPres )
            {
               real MagCC_1Cell[NCOMP_MAG];

               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
                  MHD_GetCellCenteredBFieldInPatch( MagCC_1Cell, lv, PID, i, j, k, amr->MagSg[lv] );

                  for (int v=0; v<NCOMP_MAG; v++)  MagCC[v][k][j][i] = MagCC_1Cell[v];
               }
            } // if ( OPT__FLAG_CURRENT || NeedPres )
#           endif // #ifdef MHD


//          3-2. evaluate velocity
            if ( OPT__FLAG_VORTICITY )
            {
               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
                  const real _Dens = (real)1.0 / Fluid[DENS][k][j][i];

                  Vel[0][k][j][i] = Fluid[MOMX][k][j][i]*_Dens;
                  Vel[1][k][j][i] = Fluid[MOMY][k][j][i]*_Dens;
                  Vel[2][k][j][i] = Fluid[MOMZ][k][j][i]*_Dens;
               }
            } // if ( OPT__FLAG_VORTICITY )


//          3-3. evaluate pressure
            if ( NeedPres )
            {
               const bool CheckMinPres_Yes = true;

               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
//                if applicable, compute pressure from the dual-energy variable to reduce the round-off errors
#                 ifdef DUAL_ENERGY

#                 if   ( DUAL_ENERGY == DE_ENPY )
                  Pres[k][j][i] = Hydro_DensDual2Pres( Fluid[DENS][k][j][i], Fluid[DUAL][k][j][i],
                                                       EoS_AuxArray_Flt[1], CheckMinPres_Yes, MIN_PRES );
#                 elif ( DUAL_ENERGY == DE_EINT )
#                 error : DE_EINT is NOT supported yet !!
#                 endif

#                 else // #ifdef DUAL_ENERGY

#                 ifdef MHD
                  const real Emag = (real)0.5*(  SQR( MagCC[MAGX][k][j][i] )
                                               + SQR( MagCC[MAGY][k][j][i] )
                                               + SQR( MagCC[MAGZ][k][j][i] )  );
#                 else
                  const real Emag = NULL_REAL;
#                 endif
#                 if ( EOS != EOS_GAMMA  &&  EOS != EOS_ISOTHERMAL  &&  NCOMP_PASSIVE > 0 )
                  real Passive[NCOMP_PASSIVE];
                  for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Fluid[ NCOMP_FLUID + v ][k][j][i];
#                 else
                  const real *Passive = NULL;
#                 endif

                  Pres[k][j][i] = Hydro_Con2Pres( Fluid[DENS][k][j][i], Fluid[MOMX][k][j][i], Fluid[MOMY][k][j][i],
                                                  Fluid[MOMZ][k][j][i], Fluid[ENGY][k][j][i], Passive,
                                                  CheckMinPres_Yes, MIN_PRES, PassiveFloorMask, Emag,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                                  NULL );
#                 endif // #ifdef DUAL_ENERGY ... else ...
               } // k,j,i
            } // if ( NeedPres )


//          3-4. evaluate sound speed squared
            if ( NeedCs2 )
            {
               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
#                 if ( EOS != EOS_GAMMA  &&  EOS != EOS_ISOTHERMAL  &&  NCOMP_PASSIVE > 0 )
                  real Passive[NCOMP_PASSIVE];
                  for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Fluid[ NCOMP_FLUID + v ][k][j][i];
#                 else
                  const real *Passive = NULL;
#                 endif

                  Cs2[k][j][i] = EoS_DensPres2CSqr_CPUPtr( Fluid[DENS][k][j][i], Pres[k][j][i], Passive,
                                                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
               } // k,j,i
            } // if ( NeedCs2 )


#           ifdef SRHD
//          3-5. evaluate Lorentz factor
            if ( OPT__FLAG_LRTZ_GRADIENT )
            {
               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
                  real HTilde, Factor, U1, U2, U3;
                  real Cons[NCOMP_FLUID] = { Fluid[DENS][k][j][i], Fluid[MOMX][k][j][i], Fluid[MOMY][k][j][i],
                                             Fluid[MOMZ][k][j][i], Fluid[ENGY][k][j][i] };

#                 ifdef CHECK_UNPHYSICAL_IN_FLUID
                  Hydro_IsUnphysical( UNPHY_MODE_CONS, Cons, NULL_REAL,
                                      EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                      EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                      PassiveFloorMask, ERROR_INFO, UNPHY_VERBOSE );
#                 endif

                  HTilde = Hydro_Con2HTilde( Cons, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                  Factor = Cons[0]*((real)1.0 + HTilde);
                  U1     = Cons[1]/Factor;
                  U2     = Cons[2]/Factor;
                  U3     = Cons[3]/Factor;

                  Lrtz[k][j][i] = SQRT( (real)1.0 + SQR(U1) + SQR(U2) + SQR(U3) );
               } // i,j,k
            } // if ( OPT__FLAG_LRTZ_GRADIENT )
#           endif // #ifdef SRHD


#           ifdef SUPPORT_GRACKLE
//          3-6. evaluate cooling length
            if ( OPT__FLAG_COOLING_LEN )
            {
               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
                  const real TCool = Grackle_TCool[ (LocalID*CUBE(PS1) + k*SQR(PS1) + j*PS1 + i) ];
//                positive cooling time returned by Grackle corresponds to heating, and vice versa
//                --> cooling length is defined as (|cooling time| x sound speed) and is only applied to cooling
                  LCool[k][j][i] = ( TCool < (real)0.0 ) ? ( FABS(TCool) * SQRT(Cs2[k][j][i]) ) : HUGE_NUMBER;
               } // k,j,i
            } // if ( OPT__FLAG_COOLING_LEN )
#           endif // #ifdef SUPPORT_GRACKLE
#           endif // #if ( MODEL == HYDRO )


//          3-7. evaluate the averages and slopes along x/y/z for Lohner
            if ( Lohner_NVar > 0 )
               Prepare_for_Lohner( OPT__FLAG_LOHNER_FORM, Lohner_Var+LocalID*Lohner_Stride, Lohner_Ave, Lohner_Slope,
                                   Lohner_NVar );


//          3-8. count the number of particles and/or particle mass density on each cell
#           ifdef PARTICLE
            if ( OPT__FLAG_NPAR_CELL  ||  OPT__FLAG_PAR_MASS_CELL )
            {
               long      *ParList = NULL;
               int        NParThisPatch;
               bool       UseInputMassPos;
               real_par **InputMassPos = NULL;
               long_par **InputType    = NULL;

//             determine the number of particles and the particle list
               if ( amr->patch[0][lv][PID]->son == -1 )
               {
                  NParThisPatch   = amr->patch[0][lv][PID]->NPar;
                  ParList         = amr->patch[0][lv][PID]->ParList;
                  UseInputMassPos = false;
                  InputMassPos    = NULL;
                  InputType       = NULL;

#                 ifdef DEBUG_PARTICLE
                  if ( amr->patch[0][lv][PID]->NPar_Copy != -1 )
                     Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                                lv, PID, amr->patch[0][lv][PID]->NPar_Copy );
#                 endif
               }

               else
               {
                  NParThisPatch   = amr->patch[0][lv][PID]->NPar_Copy;
#                 ifdef LOAD_BALANCE
                  ParList         = NULL;
                  UseInputMassPos = true;
                  InputMassPos    = amr->patch[0][lv][PID]->ParAttFlt_Copy;
                  InputType       = amr->patch[0][lv][PID]->ParAttInt_Copy;
#                 else
                  ParList         = amr->patch[0][lv][PID]->ParList_Copy;
                  UseInputMassPos = false;
                  InputMassPos    = NULL;
                  InputType       = NULL;
#                 endif

#                 ifdef DEBUG_PARTICLE
                  if ( amr->patch[0][lv][PID]->NPar != 0 )
                     Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d != 0 !!\n",
                                lv, PID, amr->patch[0][lv][PID]->NPar );
#                 endif
               }

#              ifdef DEBUG_PARTICLE
               if ( NParThisPatch < 0 )
                  Aux_Error( ERROR_INFO, "NPar (%d) has not been calculated (lv %d, PID %d) !!\n",
                             NParThisPatch, lv, PID );

               if ( NParThisPatch > 0 )
               {
                  if ( UseInputMassPos )
                  {
                     if ( InputMassPos[PAR_MASS] == NULL  ||  InputMassPos[PAR_POSX] == NULL  ||
                          InputMassPos[PAR_POSY] == NULL  ||  InputMassPos[PAR_POSZ] == NULL )
                        Aux_Error( ERROR_INFO, "InputMassPos[0/1/2/3] == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                                   NParThisPatch, lv, PID );
                     if ( InputType[PAR_TYPE] == NULL )
                        Aux_Error( ERROR_INFO, "InputType[0] == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                                   NParThisPatch, lv, PID );
                  }

                  else if ( ParList == NULL )
                  Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                             NParThisPatch, lv, PID );
               }
#              endif

//             deposit particle mass onto grids
//             --> for OPT__FLAG_NPAR_CELL, set UnitDens_Yes
//             --> for OPT__FLAG_PAR_MASS_CELL, set UnitDens_No and note that Par_MassAssignment() returns
//                 **density** instead of mass
//                 --> must multiply with the cell volume before checking the particle **mass** refinement criterion
               if ( OPT__FLAG_NPAR_CELL )
               Par_MassAssignment( ParList, NParThisPatch, PAR_INTERP_NGP, ParCount[0][0], PS1,
                                   amr->patch[0][lv][PID]->EdgeL, amr->dh[lv], PredictPos_No, NULL_REAL,
                                   InitZero_Yes, Periodic_No, NULL, UnitDens_Yes, CheckFarAway_No,
                                   UseInputMassPos, InputMassPos, InputType );

               if ( OPT__FLAG_PAR_MASS_CELL )
               Par_MassAssignment( ParList, NParThisPatch, PAR_INTERP_NGP, ParDens [0][0], PS1,
                                   amr->patch[0][lv][PID]->EdgeL, amr->dh[lv], PredictPos_No, NULL_REAL,
                                   InitZero_Yes, Periodic_No, NULL, UnitDens_No,  CheckFarAway_No,
                                   UseInputMassPos, InputMassPos, InputType );
            } // if ( OPT__FLAG_NPAR_CELL  ||  OPT__FLAG_PAR_MASS_CELL )
#           endif // #ifdef PARTICLE


//          4. flag patches for refinement
//          --> check patch-based criteria before cell-based criteria for better performance (except for OPT__FLAG_INTERFERENCE)

            bool NextPatch = false;

//          4-1. ELBDM interference check
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            if ( OPT__FLAG_INTERFERENCE  &&  !amr->use_wave_flag[lv]  &&  lv < MAX_LEVEL  &&  !NextPatch )
               NextPatch = Flag_IterateCells( 1, lv, PID, dv, Fluid, Pot, MagCC, Vel, Pres, Lrtz, LCool, Cs2,
                                              Lohner_Var+LocalID*Lohner_Stride, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                                              ParCount, ParDens, Interf_Var+LocalID*Interf_Stride, Spectral_Cond,
                                              SibID_Array, FlagBuf, JeansCoeff_Factor );
#           endif


#           ifdef PARTICLE
//          4-2. flag based on the number of particles per patch (which doesn't need to go through all cells one-by-one)
            if ( OPT__FLAG_NPAR_PATCH != 0  &&  lv < MAX_LEVEL  &&  !NextPatch )
            {
               const int NParFlag = FlagTable_NParPatch[lv];
               int NParThisPatch;

               if ( amr->patch[0][lv][PID]->son == -1 )  NParThisPatch = amr->patch[0][lv][PID]->NPar;
               else                                      NParThisPatch = amr->patch[0][lv][PID]->NPar_Copy;

#              ifdef DEBUG_PARTICLE
               if ( NParThisPatch < 0 )
                  Aux_Error( ERROR_INFO, "NPar (%d) has not been calculated (lv %d, PID %d) !!\n",
                             NParThisPatch, lv, PID );
#              endif

               if ( NParThisPatch > NParFlag )
               {
//                flag itself
                  amr->patch[0][lv][PID]->flag = true;

//                flag all siblings for OPT__FLAG_NPAR_PATCH == 2
                  if ( OPT__FLAG_NPAR_PATCH == 2 )
                  {
                     for (int s=0; s<26; s++)
                     {
                        const int SibPID = amr->patch[0][lv][PID]->sibling[s];

#                       ifdef DEBUG_PARTICLE
                        if ( SibPID == -1 )
                           Aux_Error( ERROR_INFO, "SibPID == -1 --> proper-nesting check failed !!\n" );

                        if ( SibPID <= SIB_OFFSET_NONPERIODIC  &&  OPT__NO_FLAG_NEAR_BOUNDARY )
                           Aux_Error( ERROR_INFO, "SibPID (%d) <= %d when OPT__NO_FLAG_NEAR_BOUNDARY is on !!\n",
                                      SibPID, SIB_OFFSET_NONPERIODIC );
#                       endif

//                      note that we can have SibPID <= SIB_OFFSET_NONPERIODIC when OPT__NO_FLAG_NEAR_BOUNDARY == false
                        if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
                     }

//                   skip the remaining refinement flag checks for better performance
//                   --> this is safe because all sibling patches have already been flagged,
//                       and the FLAG_BUFFER_SIZE* checks are therefore no longer needed
                     NextPatch = true;
                  } // if ( OPT__FLAG_NPAR_PATCH == 2 )
               } // if ( NParThisPatch > NParFlag )
            } // if ( OPT__FLAG_NPAR_PATCH != 0 )


//          4-3. check if this patch contains any particles flagged for refinement
            if (  ( OPT__FLAG_PAR_TARGET == FLAG_PAR_MUST || OPT__FLAG_PAR_TARGET == FLAG_PAR_BOTH )  &&
                  ( ! amr->patch[0][lv][PID]->flag || OPT__FLAG_PAR_TARGET_SIB )  &&
                  lv < MAX_LEVEL  &&  !NextPatch )
            {
               if (  Par_Flag_TargetParticle( lv, PID, FLAG_PAR_MUST )  )
               {
                  amr->patch[0][lv][PID]->flag = true;

//                flag all siblings for OPT__FLAG_PAR_TARGET_SIB
                  if ( OPT__FLAG_PAR_TARGET_SIB )
                  {
                     for (int s=0; s<26; s++)
                     {
                        const int SibPID = amr->patch[0][lv][PID]->sibling[s];

#                       ifdef DEBUG_PARTICLE
                        if ( SibPID == -1 )
                           Aux_Error( ERROR_INFO, "SibPID == -1 --> proper-nesting check failed !!\n" );

                        if ( SibPID <= SIB_OFFSET_NONPERIODIC  &&  OPT__NO_FLAG_NEAR_BOUNDARY )
                           Aux_Error( ERROR_INFO, "SibPID (%d) <= %d when OPT__NO_FLAG_NEAR_BOUNDARY is on !!\n",
                                      SibPID, SIB_OFFSET_NONPERIODIC );
#                       endif

//                      note that we can have SibPID <= SIB_OFFSET_NONPERIODIC when OPT__NO_FLAG_NEAR_BOUNDARY == false
                        if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
                     }

//                   skip the remaining refinement flag checks for better performance
//                   --> this is safe because all sibling patches have already been flagged,
//                       and the FLAG_BUFFER_SIZE* checks are therefore no longer needed
                     NextPatch = true;
                  } // if ( OPT__FLAG_PAR_TARGET_SIB )
               } // if (  Par_Flag_TargetParticle( lv, PID, FLAG_PAR_MUST )  )
            } // if (  ( ! amr->patch[0][lv][PID]->flag || OPT__FLAG_PAR_TARGET_SIB )  && ... )
#           endif // #ifdef PARTICLE


//          4-4. cell-based refinement criteria
            if ( lv < MAX_LEVEL  &&  !NextPatch )
               NextPatch = Flag_IterateCells( 2, lv, PID, dv, Fluid, Pot, MagCC, Vel, Pres, Lrtz, LCool, Cs2,
                                              Lohner_Var+LocalID*Lohner_Stride, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                                              ParCount, ParDens, Interf_Var+LocalID*Interf_Stride, Spectral_Cond,
                                              SibID_Array, FlagBuf, JeansCoeff_Factor );


//          5. check the derefinement criterion of Lohner if required
//          --> do it separately from all other refinement criteria since it
//              (a) does not flag sibling patches and
//              (b) only applies to patches with sons but have been marked for derefinement
//          --> it can suppress derefinement (by having derefinement thresholds lower than refinement thresholds)
            if ( Lohner_NVar > 0  &&  FlagTable_Lohner[lv][1] < FlagTable_Lohner[lv][0]  &&
                 !amr->patch[0][lv][PID]->flag  &&  amr->patch[0][lv][PID]->son != -1 )
            {
               bool Skip = false;

               for (int k=0; k<PS1; k++)  {  if ( Skip )  break;
               for (int j=0; j<PS1; j++)  {  if ( Skip )  break;
               for (int i=0; i<PS1; i++)  {  if ( Skip )  break;

//                check Lohner only if density is greater than the minimum threshold
#                 ifdef DENS
                  if ( Fluid[DENS][k][j][i] >= FlagTable_Lohner[lv][4] )
#                 endif
                  if (  Flag_Lohner( i, j, k, OPT__FLAG_LOHNER_FORM,
                                     Lohner_Var+LocalID*Lohner_Stride, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                                     FlagTable_Lohner[lv][1], FlagTable_Lohner[lv][2], FlagTable_Lohner[lv][3] )  )
                  {
//                   flag itself
                     amr->patch[0][lv][PID]->flag = true;

//                   skip all remaining cells
                     Skip = true;
                  }
               }}} // i,j,k
            } // if ( Lohner_NVar > 0 ... )
         } // for (int LocalID=0; LocalID<8; LocalID++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


      delete [] MagCC;
      delete [] Vel;
      delete [] Pres;
      delete [] Cs2;
      delete [] Lrtz;
      delete [] LCool;
      delete [] ParCount;
      delete [] ParDens;
      delete [] Lohner_Var;
      delete [] Lohner_Ave;
      delete [] Lohner_Slope;
      delete [] Interf_Var;
      delete [] Spectral_Var;
      delete [] Grackle_TCool;

   } // OpenMP parallel region


// 6. free memory allocated by Par_CollectParticle2OneLevel()
#  ifdef PARTICLE
   Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_No, FaSibBufPatch_No );
#  endif


// 7. apply the proper-nesting constraint again (should also apply to the buffer patches)
// --> necessary because of the flag buffers
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->num[lv]; PID++)
   {
      for (int sib=0; sib<26; sib++)
      {
//       do not check if sibling[]<-1 to allow for refinement around boundaries first
         if ( amr->patch[0][lv][PID]->sibling[sib] == -1 )
         {
            amr->patch[0][lv][PID]->flag = false;

//          enforce proper-nesting constraint for use_wave_flag
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            amr->patch[0][lv][PID]->switch_to_wave_flag = false;
#           endif

            break;
         }
      }

//    check further if refinement around boundaries is forbidden
      if ( OPT__NO_FLAG_NEAR_BOUNDARY  &&  amr->patch[0][lv][PID]->flag )
      {
        for (int d=0; d<3; d++)
        {
           int CornerL = amr->patch[0][lv][PID]->corner[d];
           int CornerR = CornerL + Mis_Cell2Scale( PS1, lv );

           if ( CornerL <= 0                + NoRefineBoundaryRegion  ||
                CornerR >= amr->BoxScale[d] - NoRefineBoundaryRegion    )
            {
               amr->patch[0][lv][PID]->flag = false;

#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               amr->patch[0][lv][PID]->switch_to_wave_flag = false;
#              endif

               break;
            }
         }
      }
   } // for (int PID=0; PID<amr->num[lv]; PID++)


// 8. invoke the load-balance functions
#  ifdef LOAD_BALANCE
   if ( UseLBFunc == USELB_YES )
   {
//    grandson check
      if ( lv < NLEVEL-2 )    LB_GrandsonCheck( lv );

//    exchange the flagged buffer patches
      LB_ExchangeFlaggedBuffer( lv );

      return;
   }
#  endif


// 9. grandson check
   int SonPID;

   if ( lv < NLEVEL-2 )
   {
//###ISSUE: use atomic ??
//#     pragma omp parallel for private( SonPID )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         SonPID = amr->patch[0][lv][PID]->son;

         if ( SonPID != -1 )
         {
            for (int LocalID=0; LocalID<8; LocalID++)
            {
               if ( amr->patch[0][lv+1][SonPID+LocalID]->son != -1 )    // if grandson exists
               {
//                flag the corresponding siblings of patch PID
                  Flag_Grandson( lv, PID, LocalID );

//                flag the patch PID
                  amr->patch[0][lv][PID]->flag = true;
               }
            }
         }
      }
   }  // if ( lv < NLEVEL-2 )


// 10. set up the BounFlag_NList and BounFlag_PosList
   Buf_RecordBoundaryFlag( lv );

} // FUNCTION : Flag_Real



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Grandson
// Description :  Properly flag the siblings of patch "PID" at level "lv" to ensure that the proper-nesting
//                condition is satisfied at level "lv+2"
//
// Note        :  1. Properly take care of the non-periodic BC which does not allocate buffer patches
//                2. No need to validate the proper-nesting constraint when flagging the sibling patches on level "lv"
//                   (in other words, these sibling patches on lv must be allowed to be flagged) since the existence of
//                   grandson patches should already enforce this constraint
//
// Parameter   :  lv      : Target level to be flagged
//                PID     : Target patch ID at level "lv"
//                LocalID : Index of son (0~7) which has its own son (grandson of patch "PID" at level "lv")
//-------------------------------------------------------------------------------------------------------
void Flag_Grandson( const int lv, const int PID, const int LocalID )
{

   int SibPID;

   switch ( LocalID )
   {
      case 0:
//       must ensure SibPID >= 0 due to the non-periodic B.C., for which we can have SibPID = SIB_OFFSET_NONPERIODIC-sib < -1
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[18];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[18];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 1:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[19];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[19];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 2:
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[20];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[20];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 3:
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[22];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[22];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 4:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[21];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[21];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 5:
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[24];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[24];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 6:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[23];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[23];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      case 7:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[25];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[25];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LocalID", LocalID );

   } // switch ( LocalID )

} // FUNCTION : Flag_Grandson



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_IterateCells
// Description :  Check whether any cell within the target patch satisfies the cell-based refinement criteria
//
// Note        :  1. Invoke either ELBDM_Flag_Interference() or Flag_Check(), depending on the "Mode" parameter
//                   --> ELBDM_Flag_Interference() is separated from Flag_Check() for the following reasons:
//                       (1) The ELBDM interference check must be performed before all other refinement checks
//                           in order to set switch_to_wave_flag correctly
//                       (2) Patch-based flag checks are performed before the remaining cell-based checks in Flag_Check()
//                           to improve performance
//                           --> Specifically, flag checks are performed in the following order:
//                           a. Cell-based check of ELBDM_Flag_Interference()
//                           b. Patch-based flag checks (e.g., OPT__FLAG_NPAR_PATCH and OPT__FLAG_PAR_TARGET)
//                           c. Cell-based flag checks in Flag_Check()
//                2. Invoked by Flag_Real()
//
// Parameter   :  Mode              : 1 --> Call ELBDM_Flag_Interference()
//                                    2 --> Call Flag_Check()
//                lv                : Target refinement level
//                PID               : Target patch ID
//                dv                : Cell volume at the target level
//                Fluid             : Input fluid array (with NCOMP_TOTAL components)
//                Pot               : Input potential array
//                MagCC             : Input cell-centered B field array
//                Vel               : Input velocity array
//                Pres              : Input pressure array
//                Lrtz              : Input Lorentz factor array
//                LCool             : Input cooling length array
//                Cs2               : Input sound speed squared array
//                Lohner_Var        : Input array storing the variables for the Lohner error estimator
//                Lohner_Ave        : Input array storing the averages for the Lohner error estimator
//                Lohner_Slope      : Input array storing the slopes for the Lohner error estimator
//                Lohner_NVar       : Number of variables stored in Lohner_Ave and Lohner_Slope
//                ParCount          : Input array storing the number of particles on each cell
//                                    (note that it has the **real** type)
//                ParDens           : Input array storing the particle mass density on each cell
//                Interf_Var        : Input array storing the density and phase for the interference condition
//                Spectral_Cond     : Input variable storing the spectral refinement condition
//                SibID_Array       : sibling indices for flag buffers
//                FlagBuf           : Flag buffer size
//                JeansCoeff_Factor : Factor for OPT__FLAG_JEANS
//
// Return      :  true  : The flag buffers of the target patch have flagged all 26 sibling patches
//                false : Otherwise
//-------------------------------------------------------------------------------------------------------
bool Flag_IterateCells( const int Mode, const int lv, const int PID, const real dv,
                        const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real MagCC[][PS1][PS1][PS1],
                        const real Vel[][PS1][PS1][PS1], const real Pres[][PS1][PS1], const real Lrtz[][PS1][PS1],
                        const real LCool[][PS1][PS1], const real (*Cs2)[PS1][PS1],
                        const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                        const real ParCount[][PS1][PS1], const real ParDens[][PS1][PS1],
                        const real *Interf_Var, const real Spectral_Cond,
                        const int SibID_Array[][3][3], const int FlagBuf, const real JeansCoeff_Factor )
{

   if ( Mode != 1  &&  Mode != 2 )  Aux_Error( ERROR_INFO, "unsupported Mode = %d !!\n", Mode );


// do not flag patches at or above MAX_LEVEL
   if ( lv >= MAX_LEVEL )  return false;


// loop over all cells within the target patch
   for (int k=0; k<PS1; k++)  {  const int k_start = ( k - FlagBuf < 0    ) ? 0 : 1;
                                 const int k_end   = ( k + FlagBuf >= PS1 ) ? 2 : 1;
   for (int j=0; j<PS1; j++)  {  const int j_start = ( j - FlagBuf < 0    ) ? 0 : 1;
                                 const int j_end   = ( j + FlagBuf >= PS1 ) ? 2 : 1;
   for (int i=0; i<PS1; i++)  {  const int i_start = ( i - FlagBuf < 0    ) ? 0 : 1;
                                 const int i_end   = ( i + FlagBuf >= PS1 ) ? 2 : 1;

//    retrieve the adiabatic index for Jeans length refinement criterion
#     if ( MODEL == HYDRO  &&  defined GRAVITY )
      const real JeansCoeff = ( OPT__FLAG_JEANS )
                            ? JeansCoeff_Factor * Cs2[k][j][i] * Fluid[DENS][k][j][i] / Pres[k][j][i]
                            : NULL_REAL;
#     else
      const real JeansCoeff = NULL_REAL;
#     endif


//    check if the target cell satisfies any cell-based refinement criterion (useless pointers are always == NULL)
      bool Flag = false;

//    ELBDM interference check
      if ( Mode == 1 )
      {
#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         Flag = ELBDM_Flag_Interference( i, j, k, Interf_Var,
                                         FlagTable_Interference[lv][0], FlagTable_Interference[lv][1],
                                         FlagTable_Interference[lv][2], FlagTable_Interference[lv][3]>0.5 );

//       switch to wave solver when refining to ELBDM_FIRST_WAVE_LEVEL
         if ( Flag  &&  lv+1 >= ELBDM_FIRST_WAVE_LEVEL )    amr->patch[0][lv][PID]->switch_to_wave_flag = true;
#        endif
      }

//    other criteria
      else
      {
         Flag = Flag_Check( lv, PID, i, j, k, dv, Fluid, Pot, MagCC, Vel, Pres, Lrtz, LCool,
                            Lohner_Var, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                            ParCount, ParDens, JeansCoeff, Spectral_Cond );
      }


//    flag patches
      if ( Flag )
      {
//       flag itself
         amr->patch[0][lv][PID]->flag = true;

//       flag sibling patches according to the size of FlagBuf
         for (int kk=k_start; kk<=k_end; kk++)
         for (int jj=j_start; jj<=j_end; jj++)
         for (int ii=i_start; ii<=i_end; ii++)
         {
            const int SibID = SibID_Array[kk][jj][ii];

            if ( SibID != 999 )
            {
               const int SibPID = amr->patch[0][lv][PID]->sibling[SibID];

#              ifdef GAMER_DEBUG
               if ( SibPID == -1 )
                  Aux_Error( ERROR_INFO, "SibPID == -1 --> proper-nesting check failed !!\n" );

               if ( SibPID <= SIB_OFFSET_NONPERIODIC  &&  OPT__NO_FLAG_NEAR_BOUNDARY )
                  Aux_Error( ERROR_INFO, "SibPID (%d) <= %d when OPT__NO_FLAG_NEAR_BOUNDARY is on !!\n",
                             SibPID, SIB_OFFSET_NONPERIODIC );
#              endif

//             note that we can have SibPID <= SIB_OFFSET_NONPERIODIC when OPT__NO_FLAG_NEAR_BOUNDARY == false
               if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

//             switch_to_wave_flag should be consistent with the flag buffer
#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               if ( Mode == 1  &&  amr->patch[0][lv][PID]->switch_to_wave_flag  &&  SibPID >= 0 )
                  amr->patch[0][lv][SibPID]->switch_to_wave_flag = true;
#              endif
            } // if ( SibID != 999 )
         } // ii,jj,kk

//       to improve performance, return true immediately if all 26 siblings have been flagged by the target patch's flag buffers
         if ( k_end - k_start == 2  &&
              j_end - j_start == 2  &&
              i_end - i_start == 2   )    return true;
      } // if ( Flag )
   }}} // k, j, i

   return false;

} // FUNCTION : Flag_IterateCells
