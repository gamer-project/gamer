#include "GAMER.h"
#include "CUFLU.h"


static IntSchemeFunc_t Int_SelectScheme( const IntScheme_t IntScheme );

#if ( MODEL == HYDRO )
static void Interpolate_Iterate( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                                 real FData[], const int FSize[3], const int FStart[3],
                                 const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase,
                                 const bool Monotonic[], const bool OppSign0thOrder,
                                 const IntPrim_t IntPrim, const ReduceOrFixMonoCoeff_t ReduceMonoCoeff,
                                 const real CMag[], const real FMag[][NCOMP_MAG] );
#endif

void Int_MinMod1D  ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_MinMod3D  ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_vanLeer   ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_CQuadratic( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_Quadratic ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_CQuartic  ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_Quartic   ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
#ifdef SUPPORT_SPECTRAL_INT
void Int_Spectral  ( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                     real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Interpolate
// Description :  Perform spatial interpolation
//
// Note        :  1. Use IntScheme to determine the interpolation scheme
//                2. CData[] may be overwritten
//                3. Switch to Interpolate_Iterate() when enabling ReduceMonoCoeff
//                   --> Only applicable for HYDRO with AllCons==true
//
// Parameter   :  CData           : Input coarse-grid array (which may be overwritten)
//                CSize           : Size of CData[]
//                CStart          : (x,y,z) starting indices to perform interpolation on CData[]
//                CRange          : Number of coarse cells along each direction to perform interpolation
//                FData           : Output fine-grid array
//                FSize           : Size of FData[]
//                FStart          : (x,y,z) starting indices to store the interpolation results
//                NComp           : Number of components in the CData and FData array
//                IntScheme       : Interpolation scheme
//                                  --> currently supported schemes include
//                                      INT_MINMOD3D : MinMod-3D
//                                      INT_MINMOD1D : MinMod-1D
//                                      INT_VANLEER  : vanLeer
//                                      INT_CQUAD    : conservative quadratic
//                                      INT_QUAD     : quadratic
//                                      INT_CQUAR    : conservative quartic
//                                      INT_QUAR     : quartic
//                UnwrapPhase     : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//                Monotonic       : Ensure that all interpolation results are monotonic
//                                  --> Useful for interpolating positive-definite variables, such as density, energy, ...
//                OppSign0thOrder : Apply 0th-order interpolation if the values to be interpolated change
//                                  signs in adjacent cells
//                                  --> See Int_MinMod1D() for details
//                AllCons         : Input fields include all conserved hydro variables (i.e., _TOTAL)
//                                  --> For determining whether ReduceMonoCoeff and IntPrim are applicable
//                                  --> HYDRO only and must have NComp==NCOMP_TOTAL
//                IntPrim         : Whether or not switch from conserved to primitive variables when interpolation fails
//                                  --> Must enable AllCons
//                ReduceMonoCoeff : (true/false) --> (reduce/fix) the monotonic coefficient when interpolation fails
//                                  --> Must enable AllCons
//                C/FMag_IntIter  : Coarse/Fine-grid, cell-centered B field for Interpolate_Iterate()
//                                  --> Note that the two arrays have different format
//
// Return      :  FData[]
//-------------------------------------------------------------------------------------------------------
void Interpolate( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData[], const int FSize[3], const int FStart[3],
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase,
                  const bool Monotonic[], const bool OppSign0thOrder, const bool AllCons,
                  const IntPrim_t IntPrim, const ReduceOrFixMonoCoeff_t ReduceMonoCoeff,
                  const real CMag_IntIter[], const real FMag_IntIter[][NCOMP_MAG] )
{

// check
#  ifdef GAMER_DEBUG
   int NGhost, NSide;
   Int_Table( IntScheme, NSide, NGhost );

   for (int d=0; d<3; d++)
   {
      if ( CSize[d] < 0 )  Aux_Error( ERROR_INFO, "CSize[%d] = %d < 0 !!\n", d, CSize[d] );
      if ( FSize[d] < 0 )  Aux_Error( ERROR_INFO, "FSize[%d] = %d < 0 !!\n", d, FSize[d] );
      if ( CStart[d] < NGhost  ||  CStart[d] >= CSize[d]-NGhost )
         Aux_Error( ERROR_INFO, "incorrect CStart[%d] = %d (Min = %d, Max = %d) !!\n",
                    d, CStart[d], NGhost, CSize[d]-NGhost-1 );
      if ( FStart[d] < 0  ||  FStart[d] >= FSize[d]-1 )
         Aux_Error( ERROR_INFO, "incorrect FStart[%d] = %d (Min = %d, Max = %d) !!\n",
                    d, FStart[d], 0, FSize[d]-2 );
      if ( CStart[d]+CRange[d] >= CSize[d]-NGhost+1 )
         Aux_Error( ERROR_INFO, "incorrect CStart[%d] (%d) + CRange[%d] (%d) = %d (Max = %d) !!\n",
                    d, CStart[d], d, CRange[d], CStart[d]+CRange[d], CSize[d]-NGhost );
   }

   if ( UnwrapPhase )
   {
#     if ( MODEL == ELBDM )
      if ( IntScheme == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "unsupported phase interpolation scheme (%d) !!\n", IntScheme );
#     else
      Aux_Error( ERROR_INFO, "phase unwrapping is useful in ELBDM model only !!\n" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG

#  if ( MODEL == HYDRO )
   if ( ReduceMonoCoeff || IntPrim )
   {
      if ( !AllCons )
         Aux_Error( ERROR_INFO, "AllCons == false for ReduceMonoCoeff/IntPrim !!\n" );

      if ( NComp != NCOMP_TOTAL )
         Aux_Error( ERROR_INFO, "NComp (%d) != NCOMP_TOTAL (%d) for ReduceMonoCoeff/IntPrim !!\n", NComp, NCOMP_TOTAL );

#     ifdef MHD
      if ( CMag_IntIter == NULL  &&  IntPrim )
         Aux_Error( ERROR_INFO, "CMag_IntIter == NULL for IntPrim !!\n" );

      if ( FMag_IntIter == NULL )
         Aux_Error( ERROR_INFO, "FMag_IntIter == NULL for ReduceMonoCoeff/IntPrim !!\n" );
#     endif
   }
#  endif // HYDRO


// determine whether or not to switch to Interpolate_Iterate()
#  if ( MODEL == HYDRO )
   if ( ReduceMonoCoeff || IntPrim )
      Interpolate_Iterate( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, IntScheme,
                           UnwrapPhase, Monotonic, OppSign0thOrder, IntPrim, ReduceMonoCoeff,
                           CMag_IntIter, FMag_IntIter );

   else
#  endif
   {
//    select an interpolation scheme
      IntSchemeFunc_t IntSchemeFunc = Int_SelectScheme( IntScheme );

#     ifdef GAMER_DEBUG
      if ( IntSchemeFunc == NULL )  Aux_Error( ERROR_INFO, "IntSchemeFunc == NULL!!\n" );
#     endif

//    perform interpolation
      IntSchemeFunc( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp,
                     UnwrapPhase, Monotonic, INT_MONO_COEFF, OppSign0thOrder );
   }

} // FUNCTION : Interpolate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  Interpolate_Iterate
// Description :  Perform spatial interpolation with iterations when the interpolation fails
//
// Note        :  1. Use IntScheme to determine the interpolation scheme
//                2. Locally reduce the monotonic coefficient when detecting unphysical results
//                   --> Enabled by ReduceMonoCoeff
//                   --> Not applicable for MinMod-3D, MinMod-1D, and vanLeer since they do not use monotonic coefficient
//                   --> Must have NComp==NCOMP_TOTAL
//                3. Switch from conserved to primitive variables when detecting unphysical results
//                   --> Enabled by IntPrim
//                   --> Must have NComp==NCOMP_TOTAL
//                   --> It will break conservation when used in grid refinement
//                       --> Refine() or LB_Refine_AllocateNewPatch()
//                       But it will preserve conservation when used in ghost-zone interpolation since ghost zones
//                       do not affect conservation
//                       --> InterpolateGhostZone()
//                4. Procedure
//                   a. Interpolate conserved variables with the original monotonic coefficient
//                   b. [IntPrim] If interpolation fails, interpolate primitive variables with the
//                      original monotonic coefficient
//                   c. [ReduceMonoCoeff] If interpolation fails again, interpolate conserved variables
//                      (or primitive variables when enabling IntPrim) with a reduced monotonic coefficient
//                      until either interpolation succeeds or the monotonic coefficient becomes zero
//                5. CData[] may be overwritten
//                6. Only applicable for HYDRO
//                7. When enabling INTERP_MASK (in Macro.h), only iterate on cells with unphysical results
//
// Parameter   :  See Interpolate()
//
// Return      :  FData[]
//-------------------------------------------------------------------------------------------------------
void Interpolate_Iterate( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                          real FData[], const int FSize[3], const int FStart[3],
                          const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase,
                          const bool Monotonic[], const bool OppSign0thOrder,
                          const IntPrim_t IntPrim, const ReduceOrFixMonoCoeff_t ReduceMonoCoeff,
                          const real CMag[], const real FMag[][NCOMP_MAG] )
{

//###REVISE: support FStart[*] != 0
// check
   for (int d=0; d<3; d++)
      if ( FStart[d] != 0 )   Aux_Error( ERROR_INFO, "FStart[%d] = %d != 0 !!\n", d, FStart[d] );


   const int  CSize3D         = CSize[0]*CSize[1]*CSize[2];
   const int  FSize3D         = FSize[0]*FSize[1]*FSize[2];
   const int  MonoMaxIter     = ( ReduceMonoCoeff ) ? MONO_MAX_ITER : 0;
   const int  MaxIter         = ( IntPrim ) ? MonoMaxIter+1 : MonoMaxIter;
   const bool JeansMinPres_No = false;

   int  Iteration;
   real IntMonoCoeff;
   bool Fail_AnyCell, FData_is_Prim, ContinueIteration;
   real Cons[NCOMP_TOTAL_PLUS_MAG], Temp[NCOMP_TOTAL_PLUS_MAG];   // must include B field


// select an interpolation scheme
   IntSchemeFunc_t IntSchemeFunc = Int_SelectScheme( IntScheme );

#  ifdef GAMER_DEBUG
   if ( IntSchemeFunc == NULL )  Aux_Error( ERROR_INFO, "IntSchemeFunc == NULL!!\n" );
#  endif

   real *FData_tmp = new real [NCOMP_TOTAL*FSize3D];

#  ifdef INTERP_MASK
   bool *Mask      = new bool [FSize3D];
   for (int i=0; i<FSize3D; i++)    Mask[i] = UNMASKED;
#  endif


// start iterations
   Iteration = 0;

   do {
//    1. apply the original monotonic coefficient to conserved variables
      if ( Iteration == 0 )
      {
         FData_is_Prim = false;
         IntMonoCoeff  = (real)INT_MONO_COEFF;
      }


//    2. interpolate primitive variables with the original monotonic coefficient
      else if ( Iteration == 1  &&  IntPrim )
      {
//       conserved --> primitive
         for (int i=0; i<CSize3D; i++)
         {
//          assuming **all** elements of CData[] and CMag[] are filled in (i.e., no unused cells)
            for (int v=0; v<NCOMP_TOTAL; v++)   Cons[v] = CData[ CSize3D*v + i ];
#           ifdef MHD
            for (int v=0; v<NCOMP_MAG; v++)
            {
               const real B = CMag[ CSize3D*v + i ];
               Cons[ MAG_OFFSET + v ] = B;

//             abort if the coarse-grid B field is unphysical
               if ( ! Aux_IsFinite(B) )   Aux_Error( ERROR_INFO, "unphysical coarse-grid B field (B%d = %14.7e) !!\n", v, B );
            }
#           endif

            Hydro_Con2Pri( Cons, Temp, MIN_PRES, PassiveFloorMask,
                           OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                           JeansMinPres_No, NULL_REAL,
                           EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL, NULL );

//          no need to copy the magnetic field here
            for (int v=0; v<NCOMP_TOTAL; v++)   CData[ CSize3D*v + i ] = Temp[v];
         } // for (int i=0; i<CSize3D; i++)

         FData_is_Prim = true;
      } // else if ( Iteration == 1  &&  IntPrim )


//    3. reduce the original monotonic coefficient
      else
      {
//       vanLeer, MinMod-3D, and MinMod-1D do not use monotonic coefficient --> terminate
         if ( IntScheme == INT_VANLEER  ||  IntScheme == INT_MINMOD3D  ||  IntScheme == INT_MINMOD1D )
            Aux_Error( ERROR_INFO, "INT_VANLEER/INT_MINMOD3D/INT_MINMOD1D do not support MONO_MAX_ITER != 0 !!\n" );

//       no need to worry about MonoMaxIter==0 since it is guaranteed to be positive here
         IntMonoCoeff -= (real)INT_MONO_COEFF / (real)MonoMaxIter;

//       ensure IntMonoCoeff is non-negative
         IntMonoCoeff = FMAX( IntMonoCoeff, (real)0.0 );
      }


//    4. perform interpolation
      IntSchemeFunc( CData, CSize, CStart, CRange, FData_tmp, FSize, FStart, NComp,
                     UnwrapPhase, Monotonic, IntMonoCoeff, OppSign0thOrder );


      Fail_AnyCell = false;

      for (int i=0; i<FSize3D; i++)
      {
//       skip masked cells
#        ifdef INTERP_MASK
         if ( Mask[i] == MASKED )   continue;
#        endif


//       Temp[] can store either conserved or primitive variables
         for (int v=0; v<NCOMP_TOTAL; v++)   Temp[v] = FData_tmp[ FSize3D*v + i ];
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG;   v++)   Temp[ MAG_OFFSET + v ] = FMag[i][v];
#        endif


//       5. check unphysical results
//       5-1. abort if the fine-grid B field is unphysical
#        ifdef MHD
         const real Emag = (real)0.5*( SQR(FMag[i][MAGX]) + SQR(FMag[i][MAGY]) + SQR(FMag[i][MAGZ]) );
         if ( ! Aux_IsFinite(Emag) )   Aux_Error( ERROR_INFO, "unphysical fine-grid B energy (%14.7e) !!\n", Emag );
#        else
         const real Emag = NULL_REAL;
#        endif


//       use dual-energy fix before the general check
#        ifdef DUAL_ENERGY
         const bool CheckMinPres_No = false;
         const real UseDual2FixEngy = HUGE_NUMBER;
         char dummy;    // we do not record the dual-energy status here

         if ( !FData_is_Prim )
            Hydro_DualEnergyFix( Temp[DENS], Temp[MOMX], Temp[MOMY], Temp[MOMZ], Temp[ENGY], Temp[DUAL],
                                 dummy, EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2],
                                 CheckMinPres_No, NULL_REAL, PassiveFloorMask, UseDual2FixEngy, Emag );
#        endif


//       5-2. general check
         bool Fail_ThisCell
            = Hydro_IsUnphysical( (FData_is_Prim)?UNPHY_MODE_PRIM:UNPHY_MODE_CONS, Temp, Emag,
                                  EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                  PassiveFloorMask, ERROR_INFO, UNPHY_SILENCE );


//       5-3. additional check
         real Eint=NULL_REAL;
//       check the Eint --> Pres conversion for general EoS
#        if ( EOS != EOS_GAMMA  &&  EOS != EOS_COSMIC_RAY  &&  !defined BAROTROPIC_EOS )
#           define CHECK_E2P
#        endif
#        ifdef CHECK_E2P
         real Pres=NULL_REAL;
#        endif

         if ( !Fail_ThisCell )
         {
            if ( FData_is_Prim )
            {
//             check internal energy
               if ( EoS_DensPres2Eint_CPUPtr != NULL ) {
//                convert passive scalars from mass fraction back to mass density
#                 if ( NCOMP_PASSIVE > 0 )
                  real Passive[NCOMP_PASSIVE];

                  for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Temp[ NCOMP_FLUID + v ];

                  if ( OPT__INT_FRAC_PASSIVE_LR )
                     for (int v=0; v<PassiveIntFrac_NVar; v++)    Passive[ PassiveIntFrac_VarIdx[v] ] *= Temp[DENS];
#                 else
                  const real *Passive = NULL;
#                 endif

                  Eint = EoS_DensPres2Eint_CPUPtr( Temp[DENS], Temp[ENGY], Passive,
                                                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#                 ifdef CHECK_E2P
                  Pres = EoS_DensEint2Pres_CPUPtr( Temp[DENS], Eint,       Passive,
                                                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#                 endif

//                internal energy cannot be negative (even within machine precision) since a pressure floor has been applied
//                when calling Hydro_Con2Pri()
                  if (  Hydro_IsUnphysical_Single( Eint, "interpolated internal energy", (real)0.0, HUGE_NUMBER,
                                                   ERROR_INFO, UNPHY_SILENCE )  )
                     Fail_ThisCell = true;

#                 ifdef CHECK_E2P
                  if (  Hydro_IsUnphysical_Single( Pres, "interpolated pressure",        (real)0.0, HUGE_NUMBER,
                                                   ERROR_INFO, UNPHY_SILENCE )  )
                     Fail_ThisCell = true;
#                 endif
               } // if ( EoS_DensPres2Eint_CPUPtr != NULL )
            } // if ( FData_is_Prim )

            else
            {
//             one can add additional checks for conserved variables here
            } // if ( FData_is_Prim ) ... else ...
         } // if ( !Fail_ThisCell )

         Fail_AnyCell |= Fail_ThisCell;


//       6. store results
//       6-1. skip failed cells
         if ( Fail_ThisCell )
         {
            if ( Iteration == MaxIter )
            {
               Aux_Message( stderr, "ERROR : %s() failed !!\n", __FUNCTION__ );
               Aux_Message( stderr, "NComp=%d, IntScheme=%d, UnwrapPhase=%d, Monotonic=%d, OppSign0thOrder=%d\n",
                            NComp, IntScheme, UnwrapPhase, Monotonic[0], OppSign0thOrder );
               Aux_Message( stderr, "FData_is_Prim=%d, Iter=%d, IntMonoCoeff=%13.7e\n", FData_is_Prim, Iteration, IntMonoCoeff );

               Aux_Message( stderr, "Fluid: " );
               for (int v=0; v<NCOMP_TOTAL; v++)   Aux_Message( stderr, " [%d]=%14.7e", v, Temp[v] );
               Aux_Message( stderr, "\n" );

#              ifdef MHD
               Aux_Message( stderr, "B field: " );
               for (int v=0; v<NCOMP_MAG; v++)     Aux_Message( stderr, " [%d]=%14.7e", v, FMag[i][v] );
               Aux_Message( stderr, " Emag=%14.7e\n", Emag );
#              endif

//             output additional information
               if ( FData_is_Prim ) {
//                output Eint only if it has been recalculated
                  if ( Eint != NULL_REAL )   Aux_Message( stderr, "Eint=%14.7e\n", Eint );
               }

               else {
                  const real CheckMinPres_No = false;
                  const real Pres = Hydro_Con2Pres( Temp[DENS], Temp[MOMX], Temp[MOMY], Temp[MOMZ], Temp[ENGY], Temp+NCOMP_FLUID,
                                                    CheckMinPres_No, NULL_REAL, PassiveFloorMask, Emag,
                                                    EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, &Eint );
                  Aux_Message( stderr, "Eint=%14.7e, Pres=%14.7e\n", Eint, Pres );
               }

               MPI_Exit();    // abort the simulation if interpolation fails
            } // if ( Iteration == MaxIter )

#           ifndef INTERP_MASK
            break;   // no need to check remaining cells if not using mask
#           endif
         } // if ( Fail_ThisCell )


//       6-2. store the correct results
         else
         {
//          primitive --> conserved
            if ( FData_is_Prim ) {
               Hydro_Pri2Con( Temp, Cons, OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar,
                              PassiveIntFrac_VarIdx, EoS_DensPres2Eint_CPUPtr,
                              EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                              EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

#              ifdef GAMER_DEBUG
               if (  Hydro_IsUnphysical( UNPHY_MODE_CONS, Cons, Emag,
                                         EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                         EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                         PassiveFloorMask, ERROR_INFO, UNPHY_VERBOSE )  )
                  Aux_Error( ERROR_INFO, "unphysical interpolated energy in %s() !!\n", __FUNCTION__ );
#              endif
            }

            else {
               for (int v=0; v<NCOMP_TOTAL; v++)   Cons[v] = Temp[v];
            }

//          no need to copy the magnetic field here
            for (int v=0; v<NCOMP_TOTAL; v++)   FData[ FSize3D*v + i ] = Cons[v];

#           ifdef INTERP_MASK
            Mask[i] = MASKED;
#           endif
         } // if ( Fail_ThisCell ) ... else ...
      } // for (int i=0; i<FSize3D; i++)


//    7. decide whether to abort the iteration
      if ( Fail_AnyCell  &&  Iteration < MaxIter ) {

//       if any fine cell remains failed, unmask all eight fine cells with the same parent cell
//       --> ensure conservation when disabling IntPrim
#        ifdef INTERP_MASK
         typedef bool (*vla)[ FSize[1] ][ FSize[0] ];
         vla Mask3D = ( vla )Mask;

         for (int k=0; k<FSize[2]; k+=2)  {  const int kp = k+1;
         for (int j=0; j<FSize[1]; j+=2)  {  const int jp = j+1;
         for (int i=0; i<FSize[0]; i+=2)  {  const int ip = i+1;

            if ( Mask3D[k ][j ][i ] == UNMASKED  ||
                 Mask3D[k ][j ][ip] == UNMASKED  ||
                 Mask3D[k ][jp][i ] == UNMASKED  ||
                 Mask3D[kp][j ][i ] == UNMASKED  ||
                 Mask3D[k ][jp][ip] == UNMASKED  ||
                 Mask3D[kp][jp][i ] == UNMASKED  ||
                 Mask3D[kp][j ][ip] == UNMASKED  ||
                 Mask3D[kp][jp][ip] == UNMASKED   )
            {
               Mask3D[k ][j ][i ] = UNMASKED;
               Mask3D[k ][j ][ip] = UNMASKED;
               Mask3D[k ][jp][i ] = UNMASKED;
               Mask3D[kp][j ][i ] = UNMASKED;
               Mask3D[k ][jp][ip] = UNMASKED;
               Mask3D[kp][jp][i ] = UNMASKED;
               Mask3D[kp][j ][ip] = UNMASKED;
               Mask3D[kp][jp][ip] = UNMASKED;
            }
         }}}
#        endif // #ifdef INTERP_MASK

         ContinueIteration = true;
      } // if ( Fail_AnyCell  &&  Iteration < MaxIter )

      else {
         ContinueIteration = false;
      } // if ( Fail_AnyCell  &&  Iteration < MaxIter ) ... else ...


//    8. counter increment
      Iteration ++;

   } while ( ContinueIteration );


// check if there is any missing cell
#  if ( defined GAMER_DEBUG  &&  defined INTERP_MASK )
   for (int i=0; i<FSize3D; i++)
      if ( Mask[i] == UNMASKED )    Aux_Error( ERROR_INFO, "Mask[%d] == UNMASKED !!\n", i );
#  endif


// 9. free resource
   delete [] FData_tmp;
#  ifdef INTERP_MASK
   delete [] Mask;
#  endif

} // FUNCTION : Interpolate_Iterate
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Int_SelectScheme
// Description :  Select a spatial interpolation scheme
//
// Note        :  Use the input parameter "IntScheme" to determine the adopted interpolation scheme
//
// Parameter   :  IntScheme : Interpolation scheme
//                            --> Currently supported schemes include
//                                INT_MINMOD3D : MinMod-3D
//                                INT_MINMOD1D : MinMod-1D
//                                INT_VANLEER  : vanLeer
//                                INT_CQUAD    : conservative quadratic
//                                INT_QUAD     : quadratic
//                                INT_CQUAR    : conservative quartic
//                                INT_QUAR     : quartic
//
// Return      :  IntSchemeFunc_t
//-------------------------------------------------------------------------------------------------------
static IntSchemeFunc_t Int_SelectScheme( const IntScheme_t IntScheme )
{

   switch ( IntScheme )
   {
      case INT_MINMOD3D :  return Int_MinMod3D;    break;
      case INT_MINMOD1D :  return Int_MinMod1D;    break;
      case INT_VANLEER  :  return Int_vanLeer;     break;
      case INT_CQUAD    :  return Int_CQuadratic;  break;
      case INT_QUAD     :  return Int_Quadratic;   break;
      case INT_CQUAR    :  return Int_CQuartic;    break;
      case INT_QUAR     :  return Int_Quartic;     break;
      case INT_SPECTRAL :
#                          ifdef SUPPORT_SPECTRAL_INT
                           return Int_Spectral;    break;
#                          else
                           Aux_Error( ERROR_INFO, "must enable \"SUPPORT_SPECTRAL_INT\" to use spectral interpolation (%d) !!\n",
                                      INT_SPECTRAL );
                           return NULL;            break;
#                          endif // # ifdef SUPPORT_SPECTRAL_INT
      default           :  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
   }

   return NULL;

} // FUNCTION : Int_SelectScheme
