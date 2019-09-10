#include "GAMER.h"

void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData [], const int FSize[3], const int FStart[3],
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic[], const real IntMonoCoeff );


void AdaptiveInterpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                          real FData [], const int FSize[3], const int FStart[3],
                          const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic[] )

{
     int itr = -1;
     real IntMonoCoeff;
     bool state = false;
     const int Max = 3;
     const int CSize3D = CSize[0]*CSize[1]*CSize[2];
     const int FSize3D = FSize[0]*FSize[1]*FSize[2];
     real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID];

     do {
           switch ( itr )
           {
              case -1:
                IntMonoCoeff = INT_MONO_COEFF;
                break;
  
              case  0:
                IntMonoCoeff = INT_MONO_COEFF;

                for (int i=0; i<CSize3D; i++) 
                { 
                  for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = CData[CSize3D*v+i];
                                                                                                      
                  SRHydro_Con2Pri(Cons, Prim, GAMMA, MIN_TEMP);                                       
                                                                                                      
                  for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];      
                }
                break;
  
              default:
                int Max = 3;
                real IntMonoCoeff;
                real Mono_Min = (real)0.0;

//              adaptive IntMonoCoeff
                IntMonoCoeff -= itr * ( INT_MONO_COEFF - Mono_Min ) / (real) Max ;
           }
//         interpolation
           for (int v=0; v<NComp; v++)
           Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                        FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, IntMonoCoeff );
 
   
//           check
//           for (int k=0; k<FSize[2]; k++)
//           for (int j=0; j<FSize[1]; j++)
//           for (int i=0; i<FSize[0]; i++)
//           {
//              for (int v = 0 ; v < NCOMP_FLUID;v++) 
//               Prim[v] = *(FData + i*NCOMP_FLUID*FSize[2]*FSize[1]+j*FSize[2]*NCOMP_TOTAL+k*NCOMP_TOTAL+v);
//        
//              if(SRHydro_CheckUnphysical(NULL, Prim, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, true))
//              {
//               i = j = k = FSize[0]; // break nested loop
//               state = true;
//               break;
//              }else state = false;
//           }
        
           if ( ( IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) && itr == 0 ) break;

//           itr++;
  
     } while (state && itr <= Max );

     switch ( itr )
     {
        case -1:
        break;
        case 0:
          for (int k=0; k<FSize[2]; k++)
          for (int j=0; j<FSize[1]; j++)
          for (int i=0; i<FSize[0]; i++)
          {                 
            for (int v = 0 ; v < NCOMP_FLUID;v++) 
             Prim[v] = *(FData + i*NCOMP_FLUID*FSize[2]*FSize[1]+j*FSize[2]*NCOMP_TOTAL+k*NCOMP_TOTAL+v);

            SRHydro_Pri2Con(Prim, Cons, GAMMA);
                            
            for (int v = 0 ; v < NCOMP_FLUID;v++) 
             *(FData + i*NCOMP_FLUID*FSize[2]*FSize[1]+j*FSize[2]*NCOMP_TOTAL+k*NCOMP_TOTAL+v) = Cons[v];
          }               
          break;
     }





//   check minimum energy
#    ifdef CHECK_NEGATIVE_IN_FLUID
     for (int k=0; k<FSize[2]; k++)
     for (int j=0; j<FSize[1]; j++)
     for (int i=0; i<FSize[0]; i++)
     {
     for (int v = 0 ; v < NCOMP_FLUID;v++) Con[v] = *(FData + v*NCOMP_FLUID*FSize[2]*FSize[1]+k*FSize[2]*FSize[1]+j*FSize[1]+i);
     if( SRHydro_CheckUnphysical(Con, NULL, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, true)) exit(EXIT_FAILURE);
     }

#    endif
}
