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

   const int CSize3D     = CSize[0]*CSize[1]*CSize[2];
   const int FSize3D     = FSize[0]*FSize[1]*FSize[2];
   const int MonoMaxIter = ( ReduceMonoCoeff ) ? MONO_MAX_ITER : 0;
   const int MaxIter     = ( IntPrim ) ? MonoMaxIter+1 : MonoMaxIter;

   int  Iteration    = 0;
   real IntMonoCoeff = NULL_REAL;
   bool GotFailCell  = false;

   const bool JeansMinPres_No = false;
   bool FData_is_Prim = false;
   real Cons[NCOMP_TOTAL_PLUS_MAG], Prim[NCOMP_TOTAL_PLUS_MAG];   // must include B field


// select an interpolation scheme
   IntSchemeFunc_t IntSchemeFunc = Int_SelectScheme( IntScheme );

#  ifdef GAMER_DEBUG
   if ( IntSchemeFunc == NULL )  Aux_Error( ERROR_INFO, "IntSchemeFunc == NULL!!\n" );
#  endif


// start iterations
   do {
//    1. adopt the original monotonic coefficient first
      if ( Iteration == 0 )   IntMonoCoeff = (real)INT_MONO_COEFF;


//    2. interpolate primitive variables with the original monotonic coefficient
      else if ( Iteration == 1  &&  IntPrim )
      {
         for (int i=0; i<CSize3D; i++)
         {
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

            Hydro_Con2Pri( Cons, Prim, MIN_PRES,
                           OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                           JeansMinPres_No, NULL_REAL,
                           EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

//          no need to copy the magnetic field here
            for (int v=0; v<NCOMP_TOTAL; v++)   CData[ CSize3D*v + i ] = Prim[v];
         }

         FData_is_Prim = true;
      }


//    3. reduce the original monotonic coefficient
      else
      {
//       as vanLeer, MinMod-3D, and MinMod-1D do not use monotonic coefficient, we break the loop immediately
         if ( IntScheme == INT_VANLEER  ||  IntScheme == INT_MINMOD3D  ||  IntScheme == INT_MINMOD1D )    break;

//       no need to worry about MonoMaxIter==0 since it is guaranteed to be positive here
         IntMonoCoeff -= (real)INT_MONO_COEFF / (real)MonoMaxIter;

//       ensure IntMonoCoeff is non-negative
         IntMonoCoeff = FMAX( IntMonoCoeff, (real)0.0 );
      }


//    4. perform interpolation
      IntSchemeFunc( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp,
                     UnwrapPhase, Monotonic, IntMonoCoeff, OppSign0thOrder );


//    5. check unphysical results
      GotFailCell = false;

      for (int i=0; i<FSize3D; i++)
      {
         real Temp[NCOMP_TOTAL];
         for (int v=0; v<NCOMP_TOTAL; v++)   Temp[v] = FData[ FSize3D*v + i ];

//       5-1. check the interpolation results without EoS conversion
         GotFailCell = Hydro_CheckUnphysical( (FData_is_Prim)?UNPHY_MODE_PRIM:UNPHY_MODE_CONS, Temp, NULL,
                                              ERROR_INFO, UNPHY_SILENCE );


//       5-2. check the interpolation results with EoS conversion
//            --> only check either pressure of internal energy for now
         const bool FailBeforeEoS = GotFailCell;
         real Eint=NULL_REAL, Pres=NULL_REAL;

         if ( !GotFailCell )
         {
            if ( FData_is_Prim )
            {
//             convert passive scalars from mass fraction back to mass density
#              if ( NCOMP_PASSIVE > 0 )
               real Passive[NCOMP_PASSIVE];

               for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Temp[ NCOMP_FLUID + v ];

               if ( OPT__INT_FRAC_PASSIVE_LR )
                  for (int v=0; v<PassiveIntFrac_NVar; v++)    Passive[ PassiveIntFrac_VarIdx[v] ] *= Temp[DENS];
#              else
               const real *Passive = NULL;
#              endif

               Eint = EoS_DensPres2Eint_CPUPtr( Temp[DENS], Temp[ENGY], Passive,
                                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

               if (  Hydro_CheckUnphysical( UNPHY_MODE_SING, &Eint, "interpolated internal energy", ERROR_INFO, UNPHY_SILENCE )  )
                  GotFailCell = true;
            }

            else
            {
               const bool CheckMinPres_No = false;
#              ifdef MHD
               const real Emag            = (real)0.5*( SQR(FMag[i][MAGX]) + SQR(FMag[i][MAGY]) + SQR(FMag[i][MAGZ]) );

//             abort if the fine-grid B field is unphysical
               if ( ! Aux_IsFinite(Emag) )   Aux_Error( ERROR_INFO, "unphysical fine-grid B energy (%14.7e) !!\n", Emag );

#              else
               const real Emag            = NULL_REAL;
#              endif

               Pres = Hydro_Con2Pres( Temp[DENS], Temp[MOMX], Temp[MOMY], Temp[MOMZ], Temp[ENGY], Temp+NCOMP_FLUID,
                                      CheckMinPres_No, NULL_REAL, Emag,
                                      EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

               if (  Hydro_CheckUnphysical( UNPHY_MODE_SING, &Pres, "interpolated pressure", ERROR_INFO, UNPHY_SILENCE )  )
                  GotFailCell = true;
            }
         } // if ( !GotFailCell )

         if ( GotFailCell )
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
               Aux_Message( stderr, "\n" );
#              endif

               if ( !FailBeforeEoS )
               {
                  if ( FData_is_Prim )    Aux_Message( stderr, "Eint=%14.7e\n", Eint );
                  else                    Aux_Message( stderr, "Pres=%14.7e\n", Pres );
               }

               MPI_Exit();    // abort the simulation if interpolation fails
            }

            break;
         } // if ( GotFailCell )
      } // for (int i=0; i<FSize3D; i++)


//    6. counter increment
      Iteration++;

   } while ( GotFailCell  &&  Iteration <= MaxIter );


// transform FData[] storing primitive variables back to conserved variables
   if ( FData_is_Prim )
   {
      for (int i=0; i<FSize3D; i++)
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   Prim[              v ] = FData[ FSize3D*v + i ];
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG;   v++)   Prim[ MAG_OFFSET + v ] = FMag[i][v];
#        endif

         Hydro_Pri2Con( Prim, Cons, OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                        EoS_DensPres2Eint_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

//       no need to copy the magnetic field here
         for (int v=0; v<NCOMP_TOTAL; v++)   FData[ FSize3D*v + i ] = Cons[v];
      }
   }

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
      default           :  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
   }

   return NULL;

} // FUNCTION : Int_SelectScheme
