#include "GAMER.h"

//void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
//                  real FData [], const int FSize[3], const int FStart[3],
//                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic[],
//                  const real IntMonoCoeffconst, bool OppSign0thOrder );


void InterpolateWithAdaptiveMinmod( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                                    real FData [], const int FSize[3], const int FStart[3],
                                    const int NComp, const int TVar, const IntScheme_t IntScheme, const bool UnwrapPhase,
                                    const bool Monotonic[], const bool OppSign0thOrder )
{
     int itr = -1;
     real IntMonoCoeff;
     bool state = false;
     const int Max = 3;
     const int CSize3D = CSize[0]*CSize[1]*CSize[2];
     const int FSize3D = FSize[0]*FSize[1]*FSize[2];
     real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID];

     if ( TVar == _TOTAL )
     {

        do {
              if ( itr == -1 )
              {
                 IntMonoCoeff = INT_MONO_COEFF;
              }
              else if ( itr == 0 )
              {
                 IntMonoCoeff = INT_MONO_COEFF;

                 for (int i=0; i<CSize3D; i++)
                 {
                   for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = CData[CSize3D*v+i];

                   Hydro_Con2Pri( Cons, Prim, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                                  (real)NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                  EoS_AuxArray_Int, h_EoS_Table, NULL, NULL );

                   for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];
                 }
              }
              else
              {
                 real Mono_Min = (real)0.0;

//               adaptive IntMonoCoeff
                 IntMonoCoeff -= itr * ( INT_MONO_COEFF - Mono_Min ) / (real) Max ;
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

                      if ( SRHD_CheckUnphysical( Cons, NULL, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, false  ) )
                      {
                        state = true;
                        break;
                      }
                   }
              }
              else if ( itr == 0 )
              {
                  for ( int i = 0 ;i < FSize3D; i++ )
                  {
                     for (int v = 0 ; v < NCOMP_FLUID ;v++) Prim[v] = FData[FSize3D*v+i];

                      if ( SRHD_CheckUnphysical( NULL, Prim, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, false  ) )
                      {
                        state = true;
                        break;
                      }
                  }
              }

              if ( ( IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) && itr == 0 ) break;

              itr++;

        } while (state && itr <= Max );
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

            Hydro_Pri2Con( Prim, Cons, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                           EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                           EoS_AuxArray_Int, h_EoS_Table, NULL );

            for (int v = 0 ; v < NCOMP_FLUID ;v++) FData[FSize3D*v+i] = Cons[v];
         }
     }

//   check minimum energy
#    ifdef CHECK_FAILED_CELL_IN_FLUID
     if ( TVar == _TOTAL )
     {
       for ( int i = 0 ;i < FSize3D; i++ )
       {
          for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

          if ( SRHD_CheckUnphysical( Cons, NULL, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  ) )
          {
            printf("itr=%d\n", itr);
            exit(EXIT_FAILURE);
          }
       }
     }
#    endif
}
