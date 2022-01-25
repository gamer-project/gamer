#include "GAMER.h"


static Int_Scheme_t Int_SelectScheme( const IntScheme_t IntScheme );

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
// Description :  Perform spatial interpolation using different schemes with a locally reducing min-mod coefficient
//                when encountering a unphysical cell.
//
// Note        :  1. Use the input parameter "IntScheme" to determine the adopted interpolation scheme
//                2. Locally reduce the min-mod coefficient when unphysical results are found in the interpolated ghost zones
//                   or newly allocate patches
//                   --> we do not take this remedy for MinMod-3D, MinMod-1D, and vanLeer,
//                       since they do not involve min-mod coefficient.
//
//                3. Remedial strategy:
//                   Case1: Allocate new patches.
//                          i.e. LB_Refine_AllocateNewPatch() or Refine() invokes Interpolate()
//
//                          1a. Interpolate conserved variables with original min-mod coefficient
//
//                              1b. If failed cell was found
//                              --> Interpolate conserved variables with reducing min-mod coefficient
//                                  until the min-mod coefficient is reduced to zero
//
//                   Case2: Interpolate ghost zone
//                          i.e. InterpolateGhostZone() invokes Interpolate()
//
//                          2a. Interpolate conserved variables with original min-mod coefficient
//
//                              2b. If failed cell was found
//                              --> Interpolate primitive variables with original min-mod coefficient
//
//                              2c. If failed cell was found again
//                              --> Interpolate primitive variables with reducing min-mod coefficient
//                                  until the min-mod coefficient is reduced to zero
//
//                4.  Interpolating primitive variables still preserves conservation because ghost zones do not
//                    affect conservation.
//
//                5.  Reducing min-mod coefficient is enabled by
//                    (i)  ReduceMinModCoeff == INT_REDUCE_MINMOD_COEFF/true
//                         and
//                    (ii) NComp == NCOMP_TOTAL
//
// Parameter   :  CData             : Input coarse-grid array
//                CSize             : Size of CData[]
//                CStart            : (x,y,z) starting indices to perform interpolation on CData[]
//                CRange            : Number of coarse cells along each direction to perform interpolation
//                FData             : Output fine-grid array
//                FSize             : Size of FData[]
//                FStart            : (x,y,z) starting indices to store the interpolation results
//                NComp             : Number of components in the CData and FData array
//                IntScheme         : Interpolation scheme
//                                    --> currently supported schemes include
//                                        INT_MINMOD3D : MinMod-3D
//                                        INT_MINMOD1D : MinMod-1D
//                                        INT_VANLEER  : vanLeer
//                                        INT_CQUAD    : conservative quadratic
//                                        INT_QUAD     : quadratic
//                                        INT_CQUAR    : conservative quartic
//                                        INT_QUAR     : quartic
//                UnwrapPhase       : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//                Monotonic         : Ensure that all interpolation results are monotonic
//                                    --> Useful for interpolating positive-definite variables, such as density, energy, ...
//                OppSign0thOrder   : Apply 0th-order interpolation if the values to be interpolated change
//                                    signs in adjacent cells
//                                    --> See Int_MinMod1D() for details
//                IntPrim           : (true/false)  --> interpolate (primitive/non-primitive) variables
//                ReduceMinModCoeff : (true/false) --> (reduce/fix) min-mod coefficient
//-------------------------------------------------------------------------------------------------------
void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData [], const int FSize[3], const int FStart[3],
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase,
                  const bool Monotonic[], const bool OppSign0thOrder,
                  const ReduceOrFixMinModCoeff_t IntPrim, const ReduceOrFixMinModCoeff_t ReduceMinModCoeff )
{

//   check
#    ifdef GAMER_DEBUG
     int NGhost, NSide;
     Int_Table( IntScheme, NSide, NGhost );

     for (int d=0; d<3; d++)
     {
        if ( CSize[d] < 0 )     Aux_Error( ERROR_INFO, "CSize[%d] = %d < 0 !!\n", d, CSize[d] );
        if ( FSize[d] < 0 )     Aux_Error( ERROR_INFO, "FSize[%d] = %d < 0 !!\n", d, FSize[d] );
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
#       if ( MODEL == ELBDM )
        if ( IntScheme == INT_MINMOD1D )
        Aux_Error( ERROR_INFO, "unsupported phase interpolation scheme (%d) !!\n", IntScheme );
#       else
        Aux_Error( ERROR_INFO, "phase unwrapping is useful in ELBDM model only !!\n" );
#       endif
     }
#    endif // #ifdef GAMER_DEBUG

     if ( ReduceMinModCoeff  &&  NComp != NCOMP_TOTAL )
        Aux_Error( ERROR_INFO, "NComp (%d) != NCOMP_TOTAL (%d) for ReduceMinModCoeff == INT_REDUCE_MINMOD_COEFF !!\n", NComp, NCOMP_TOTAL );

     int Iteration = 0;
     real IntMonoCoeff = NULL_REAL;
     Int_Scheme_t Int_Scheme_FunPtr = NULL;
     bool GotFailCell = false;
     const int CSize3D = CSize[0]*CSize[1]*CSize[2];
     const int FSize3D = FSize[0]*FSize[1]*FSize[2];
     real Cons[NCOMP_TOTAL], Prim[NCOMP_TOTAL], Array[NCOMP_TOTAL];
     bool FData_is_Prim = false;
     const bool JeansMinPres_No = false;



//   select an interpolation scheme and assign it to Int_Scheme_FunPtr()
     Int_Scheme_FunPtr = Int_SelectScheme( IntScheme );


#    ifdef GAMER_DEBUG
     if ( Int_Scheme_FunPtr == NULL ) Aux_Error( ERROR_INFO, "Int_Scheme_FunPtr == NULL!!\n" );
#    endif


     if ( ReduceMinModCoeff )
     {

        do {


//         0. initialize GotFailCell as false
           GotFailCell = false;


//         1. interpolate ghost-zones or newly allocated patches with the original min-mod coefficient
           if ( Iteration == 0 )
           {
              IntMonoCoeff = (real)INT_MONO_COEFF;
           }


//         2. interpolate primitive variables with the original min-mod coefficient
//            --> this step is only for interpolating ghost zone.
           else if ( Iteration == 1 && IntPrim )
           {
//            As vanLeer, MinMod-3D, and MinMod-1D do not involve min-mod coefficient, we break the loop immediately
//            and do nothing when encountering unphysical results
              if ( IntScheme == INT_VANLEER || IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) break;

              for (int i=0; i<CSize3D; i++)
              {
                for (int v = 0 ; v < NCOMP_TOTAL ;v++) Cons[v] = CData[CSize3D*v+i];
                Hydro_Con2Pri( Cons, Prim, MIN_PRES,
                               OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                               JeansMinPres_No, NULL_REAL,
                               EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                               EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

                for (int v = 0 ; v < NCOMP_TOTAL ;v++) CData[CSize3D*v+i] = Prim[v];
              }

              FData_is_Prim = true;
           }


//         3. reduce the original min-mod coefficient
           else
           {
              // process always skip this block when MINMOD_MAX_ITER == 0
              // --> no division by zero occurs
              IntMonoCoeff -= (real)INT_MONO_COEFF / (real)MINMOD_MAX_ITER;

              // ensure IntMonoCoeff is non-negative
              IntMonoCoeff = FMAX( IntMonoCoeff, (real)0.0 );
           }


//         4. perform interpolation
           for (int v=0; v<NComp; v++)
              Int_Scheme_FunPtr( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                                 FSize, FStart, 1, UnwrapPhase, Monotonic, IntMonoCoeff, OppSign0thOrder );


//         5. check failed cell
           for ( int i = 0 ;i < FSize3D; i++ )
           {

              for (int v = 0 ; v < NCOMP_TOTAL ;v++) Array[v] = FData[FSize3D*v+i];

              GotFailCell = Hydro_CheckUnphysical( (FData_is_Prim)?UNPHY_MODE_PRIM:UNPHY_MODE_CONS, Array, NULL,
                                                   __FILE__, __FUNCTION__, __LINE__, UNPHY_SILENCE );
              if ( GotFailCell )    break;

           }


//         6. counter increment
           Iteration++;



        } while ( GotFailCell  &&  Iteration <= MINMOD_MAX_ITER );

     }
     else
     {
         for (int v=0; v<NComp; v++)
            Int_Scheme_FunPtr( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                               FSize, FStart, 1, UnwrapPhase, Monotonic, INT_MONO_COEFF, OppSign0thOrder );
     }


//   transform FData[] storing primitive variables back to conserved variables
     if ( FData_is_Prim )
     {
         for (int i=0; i<FSize3D; i++)
         {
            for (int v = 0 ; v < NCOMP_TOTAL ;v++) Prim[v] = FData[FSize3D*v+i];

            Hydro_Pri2Con( Prim, Cons,
                           OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                           EoS_DensPres2Eint_CPUPtr,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

            for (int v = 0 ; v < NCOMP_TOTAL ;v++) FData[FSize3D*v+i] = Cons[v];
         }
     }




//   check unphysical results
#    ifdef CHECK_UNPHYSICAL_IN_FLUID
     for ( int i = 0 ;i < FSize3D; i++ )
     {
        for (int v = 0 ; v < NCOMP_TOTAL ;v++) Cons[v] = FData[FSize3D*v+i];

        if ( Hydro_CheckUnphysical( UNPHY_MODE_CONS, Cons, NULL,  __FILE__, __FUNCTION__, __LINE__, UNPHY_VERBOSE ) )
          printf("NComp=%d, IntScheme=%d, UnwrapPhase=%d, Monotonic=%d, OppSign0thOrder=%d, IntPrim=%d, ReduceMinModCoeff=%d",
                  NComp,    IntScheme,    UnwrapPhase,    Monotonic,    OppSign0thOrder,    IntPrim,    ReduceMinModCoeff );
     }
#    endif
}










//-------------------------------------------------------------------------------------------------------
// Function    :  Int_SelectScheme
// Description :  Select a spatial interpolation scheme
//
// Note        :  Use the input parameter "IntScheme" to determine the adopted interpolation scheme
//
// Parameter   :  IntScheme : Interpolation scheme
//                            --> currently supported schemes include
//                                INT_MINMOD3D : MinMod-3D
//                                INT_MINMOD1D : MinMod-1D
//                                INT_VANLEER  : vanLeer
//                                INT_CQUAD    : conservative quadratic
//                                INT_QUAD     : quadratic
//                                INT_CQUAR    : conservative quartic
//                                INT_QUAR     : quartic
//
// Return      :  Int_Scheme_t Int_Scheme_FunPtr()
//-------------------------------------------------------------------------------------------------------
static Int_Scheme_t Int_SelectScheme( const IntScheme_t IntScheme )
{


   switch ( IntScheme )
   {
      case INT_MINMOD3D :
         return Int_MinMod3D;
         break;

      case INT_MINMOD1D :
         return Int_MinMod1D;
         break;

      case INT_VANLEER :
         return Int_vanLeer;
         break;

      case INT_CQUAD :
         return Int_CQuadratic;
         break;

      case INT_QUAD :
         return Int_Quadratic;
         break;

      case INT_CQUAR :
         return Int_CQuartic;
         break;

      case INT_QUAR :
         return Int_Quartic;
         break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
   } // switch ( IntScheme )


   return NULL;

} // FUNCTION : Interpolate
