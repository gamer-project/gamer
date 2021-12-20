#include "GAMER.h"


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
// Function    :  Int_SelectScheme
// Description :  Select spatial interpolation schemes
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
Int_Scheme_t Int_SelectScheme( const IntScheme_t IntScheme )
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


} // FUNCTION : Interpolate
