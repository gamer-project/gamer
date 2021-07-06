#include "GAMER.h"


//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolateWithAdaptiveMinmod
// Description :  Perform spatial interpolation using different schemes with reducing min-mod coefficient
//                when encountering unphysical cell in patch.
//
// Note        :  1. Use the input parameter "IntScheme" to determine the adopted interpolation scheme
//                2. Reduce the min-mod coefficient when unphysical results are found in interpolted patches
//                   --> we do not take this remedy for MinMod-3D, MinMod-1D, and vanLeer,
//                       since they do not involve min-mod coefficient.
//                3. Remedial strategy: interpolate conserved varibales with original min-mod coefficient
//
//                                      if failed cell found       
//                                      --> interpolate primitive varibales with original min-mod coefficient
//
//                                      if failed cell found again
//                                      --> interpolate primitive varibales with reducing min-mod coefficient
//                                          until the min-mod coefficient is reduced to zero
//                                      
//                4.  Interpolate primive variables still preserves conservation because ghost zones do not
//                    affect conservation.
//
// Parameter   :  CData           : Input coarse-grid array
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
//-------------------------------------------------------------------------------------------------------
void InterpolateWithAdaptiveMinmod( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                                    real FData [], const int FSize[3], const int FStart[3],
                                    const int NComp, const int TVar, const IntScheme_t IntScheme, const bool UnwrapPhase,
                                    const bool Monotonic[], const bool OppSign0thOrder )
{
     int itr = -1;
     real IntMonoCoeff = NULL_REAL;
     bool GotFailCell = false;
     const int Max = 3;
     const int CSize3D = CSize[0]*CSize[1]*CSize[2];
     const int FSize3D = FSize[0]*FSize[1]*FSize[2];
     real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID];

#    if ( MODEL == HYDRO  &&  defined GRAVITY )
     const real JeansMinPres_Coeff = ( JEANS_MIN_PRES ) ?
                                     NEWTON_G*SQR(JEANS_MIN_PRES_NCELL*amr->dh[JEANS_MIN_PRES_LEVEL])/(GAMMA*M_PI) : NULL_REAL;
#    else
     const real JEANS_MIN_PRES     = false;
     const real JeansMinPres_Coeff = NULL_REAL;
#    endif


     if ( TVar == _TOTAL )
     {

        do {
              if ( itr == -1 )
              {
                 IntMonoCoeff = INT_MONO_COEFF;
              }
              else if ( itr == 0 )
              {

                 // vanLeer, MinMod-3D, and MinMod-1D do not use min-mod coefficient, so we break the loop immediately
                 // and do nothing when encountering unphysical results
                 if ( IntScheme == INT_VANLEER || IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) break;

                 IntMonoCoeff = INT_MONO_COEFF;

                 for (int i=0; i<CSize3D; i++)
                 {
                   for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = CData[CSize3D*v+i];

                   Hydro_Con2Pri( Cons, Prim, MIN_PRES,
                                  OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                                  JEANS_MIN_PRES, JeansMinPres_Coeff,
                                  EoS.DensEint2Pres_FuncPtr, EoS.DensPres2Eint_FuncPtr,
                                  EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table, NULL );

                   for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];
                 }
              }
              else
              {
                 real IntMonoCoeff_Min = (real)0.0;

//               reduce min-mod coefficient
                 IntMonoCoeff -= itr * ( (real)INT_MONO_COEFF - IntMonoCoeff_Min ) / (real) Max ;
              }

//            interpolation
              for (int v=0; v<NComp; v++)
              Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                           FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, IntMonoCoeff, OppSign0thOrder );


//            check
              if ( itr == -1 )
              {
                   for ( int i = 0 ;i < FSize3D; i++ )
                   {
                      for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

                      if ( Hydro_CheckUnphysical( Cons, NULL, NULL, NULL, NULL,  __FILE__, __FUNCTION__, __LINE__, false ) )
                      {
                        GotFailCell = true;
                        break;
                      }
                   }
              }
              else if ( itr == 0 )
              {
                  for ( int i = 0 ;i < FSize3D; i++ )
                  {
                     for (int v = 0 ; v < NCOMP_FLUID ;v++) Prim[v] = FData[FSize3D*v+i];

                      if ( Hydro_CheckUnphysical( NULL, Prim, NULL, NULL, NULL,  __FILE__, __FUNCTION__, __LINE__, false ) )
                      {
                        GotFailCell = true;
                        break;
                      }
                  }
              }


              itr++;

           } while (GotFailCell && itr <= Max );
     }
     else
     {
         for (int v=0; v<NComp; v++)
         Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                      FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, INT_MONO_COEFF, OppSign0thOrder );
     }

     if ( itr > 0 )
     {
         for (int i=0; i<FSize3D; i++)
         {
            for (int v = 0 ; v < NCOMP_FLUID ;v++) Prim[v] = FData[FSize3D*v+i];

            Hydro_Pri2Con( Prim, Cons,
                           OPT__INT_FRAC_PASSIVE_LR, PassiveIntFrac_NVar, PassiveIntFrac_VarIdx,
                           EoS.DensPres2Eint_FuncPtr,
                           EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table, NULL );

            for (int v = 0 ; v < NCOMP_FLUID ;v++) FData[FSize3D*v+i] = Cons[v];
         }
     }

//   check unphysical results
#    ifdef CHECK_UNPHYSICAL_IN_FLUID
     if ( TVar == _TOTAL )
     {
       for ( int i = 0 ;i < FSize3D; i++ )
       {
          for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

          if ( Hydro_CheckUnphysical( Cons, NULL, NULL, NULL, NULL,  __FILE__, __FUNCTION__, __LINE__, true ) )
          {
            printf("itr=%d\n", itr);
            exit(EXIT_FAILURE);
          }
       }
     }
#    endif
}
