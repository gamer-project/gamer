#include "GAMER.h"


#if ( MODEL == HYDRO )
bool Unphysical( const real Fluid[], const real Gamma_m1, const int CheckMinEngyOrPres );

void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
				    const bool NormPassive, const int NNorm, const int NormIdx[],
					const bool JeansMinPres, const real JeansMinPres_Coeff );

void Hydro_Pri2Con( const real In[], real Out[], const real _Gamma_m1,
                    const bool NormPassive, const int NNorm, const int NormIdx[] );
#endif

void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData [], const int FSize[3], const int FStart[3],
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic[], const real IntMonoCoeff );


void AdaptiveInterpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                          real FData [], const int FSize[3], const int FStart[3],
                          const int NComp, const int TVar, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic[] )

{
#    if ( MODEL == HYDRO )
     const bool NormPassive_No  = false;
	 const bool JeansMinPres_No = false;                       
	 const int  CheckMinEngy = 0;
	 const int  CheckMinPres = 1;
	 const real Gamma_m1 = (real)GAMMA - (real)1.0;
	 const real _Gamma_m1 = (real)1.0/Gamma_m1;
#    elif ( MODEL == SR_HYDRO )
#    endif


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

#                  if ( MODEL == HYDRO )
				   Hydro_Con2Pri( Cons, Prim, Gamma_m1, MIN_PRES, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL);
#                  elif ( MODEL == SR_HYDRO )				   
                   SRHydro_Con2Pri(Cons, Prim, GAMMA, MIN_TEMP);
#                  endif
                                                                                                       
                   for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];      
                 }
              }
              else
              { 
                 real IntMonoCoeff;
                 real Mono_Min = (real)0.0;

//               adaptive IntMonoCoeff
                 IntMonoCoeff -= itr * ( INT_MONO_COEFF - Mono_Min ) / (real) Max;
              }

//            interpolation
              for (int v=0; v<NComp; v++)
              Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                           FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, IntMonoCoeff );
 
   
//            check
              if ( itr == -1 )
              {
                   for ( int i = 0 ;i < FSize3D; i++ )
                   {
                      for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

#                     if ( MODEL == HYDRO )
					  if (Unphysical(Cons, Gamma_m1, CheckMinPres))
#                     elif ( MODEL == SR_HYDRO )
                      if (SRHydro_CheckUnphysical(Cons, NULL, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, false))
#                     endif
                       {
                          state = true;
                          break; 
					   } else state = false;
                   }
              }
			   
//              no need to check primitive variables
//
//              else if ( itr == 0 )
//              {
//                  for ( int i = 0 ;i < FSize3D; i++ )
//                  {
//                     for (int v = 0 ; v < NCOMP_FLUID ;v++) Prim[v] = FData[FSize3D*v+i];
//
//                     if (SRHydro_CheckUnphysical(NULL, Prim, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, false))
//                      {
//                         state = true;
//                         break; 
//                      } else state = false;
//                  }
//              }


              if ( ( IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) && itr == 0 ) break;

              itr++;
  
        } while (state && itr <= Max );
     }
     else
     {
         for (int v=0; v<NComp; v++)
         Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                      FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, INT_MONO_COEFF );
     }

     if ( itr > 0 )
     {
         for (int i=0; i<FSize3D; i++)
         {
            for (int v = 0 ; v < NCOMP_FLUID ;v++) Prim[v] = FData[FSize3D*v+i];
#           if ( MODEL == HYDRO )
			Hydro_Pri2Con(Prim, Cons, _Gamma_m1, NormPassive_No, NULL_INT, NULL);
#           elif ( MODEL == SR_HYDRO )
            SRHydro_Pri2Con(Prim, Cons, GAMMA);
#           endif
            for (int v = 0 ; v < NCOMP_FLUID ;v++) FData[FSize3D*v+i] = Cons[v];
         }
     }

//   check minimum energy
#    ifdef CHECK_NEGATIVE_IN_FLUID
     if ( TVar == _TOTAL )
     {
        for ( int i = 0 ;i < FSize3D; i++ )
        {
           for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

           if (SRHydro_CheckUnphysical(Cons, NULL, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, false))
            {
               state = true;
               break; 
            } else state = false;
        }
     }
#    endif
}
