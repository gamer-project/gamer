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

//-------------------------------------------------------------------------------------------------------
// Function    :  AdaptiveInterpolate
// Description :
// Note        :
// Parameter   : [ 1] CData       : Input coarse-grid array
//               [ 2] CSize       : Size of the CData array
//               [ 3] CStart      : (x,y,z) starting indices to perform interpolation on the CData array
//               [ 4] CRange      : Number of grids in each direction to perform interpolation
//               [ 5] FData       : Output fine-grid array
//               [ 6] FSize       : Size of the FData array
//               [ 7] FStart      : (x,y,z) starting indcies to store the interpolation results
//               [ 8] TVar        : target variables to be interpolated
//               [ 9] NComp       : Number of components in the CData and FData array
//               [10] IntScheme   : Interpolation scheme
//                                  --> currently supported schemes include
//                                      INT_MINMOD3D : MinMod-3D
//                                      INT_MINMOD1D : MinMod-1D
//                                      INT_VANLEER  : vanLeer
//                                      INT_CQUAD    : conservative quadratic
//                                      INT_QUAD     : quadratic
//                                      INT_CQUAR    : conservative quartic
//                                      INT_QUAR     : quartic
//               [11] UnwrapPhase : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//               [12] Monotonic   : Ensure that all interpolation results are monotonic
//                                  --> Useful for interpolating positive-definite variables, such as density, energy, ...
//               [13] Strategy    : the strategy for interpolation (1/2)
//                                  --> 1  : step1: interpolate conserved variables with default mono coefficient
//                                           step2: interpolate primitive variables with default mono coefficient
//
//                                  --> 2  : step1: interpolate conserved variables with default mono coefficient
//                                           step2: interpolate conserved variables with reduced mono coefficient
//                                           step3: interpolate primitive variables with default mono coefficient
//
//-------------------------------------------------------------------------------------------------------

void AdaptiveInterpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                          real FData [], const int FSize[3], const int FStart[3],
                          const int NComp, const int TVar, const IntScheme_t IntScheme, const bool UnwrapPhase, 
						  const bool Monotonic[] , int Strategy )
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


     int itr = 0;

     real IntMonoCoeff;
     const real MinMonoCoeff = (real)0.0;
     const int MaxItr = 8;

     bool state = false;
     const int CSize3D = CSize[0]*CSize[1]*CSize[2];
     const int FSize3D = FSize[0]*FSize[1]*FSize[2];

     real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID];

     if ( TVar == _TOTAL )
     {
        if ( Strategy == 1 )
		{
            IntMonoCoeff = INT_MONO_COEFF;


            for (int v=0; v<NComp; v++)
            Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                         FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, IntMonoCoeff );
	
		    // check unphysical cell in fined patches
            for ( int i = 0 ;i < FSize3D; i++ )
            {
               for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

#              if ( MODEL == HYDRO )
               if (Unphysical(Cons, Gamma_m1, CheckMinPres))
#              elif ( MODEL ==SR_HYDRO )
               if (SRHydro_CheckUnphysical(Cons, NULL, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, false))
#              else
#              error: ERROR!
#              endif
               {
                  state = true;
                  break; 
			   } 
			   else state = false;
            }

			// if unphysical cell was found, interpolate primitive variables instantly
            if ( state )
			{
			   // conserved --> primitive
               for (int i=0; i<CSize3D; i++) 
               { 
                  for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = CData[CSize3D*v+i];

#                 if ( MODEL == HYDRO )
                  Hydro_Con2Pri( Cons, Prim, Gamma_m1, MIN_PRES, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL);
#                 elif ( MODEL == SR_HYDRO )				   
                  SRHydro_Con2Pri(Cons, Prim, GAMMA, MIN_TEMP);
#                 else
#                 error: ERROR!
#                 endif
                                                                                                      
                  for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];      
               }
               
			   	
               for (int v=0; v<NComp; v++)
               Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                            FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, IntMonoCoeff );

               // primitive --> conserved
               for (int i=0; i<FSize3D; i++)
               {
                  for (int v = 0 ; v < NCOMP_FLUID ;v++) Prim[v] = FData[FSize3D*v+i];
#                 if ( MODEL == HYDRO )
                  Hydro_Pri2Con( Prim, Cons, _Gamma_m1, NormPassive_No, NULL_INT, NULL );
#                 elif ( MODEL == SR_HYDRO )
                  SRHydro_Pri2Con(Prim, Cons, GAMMA);
#                 else
#                 error: ERROR!
#                 endif
                  for (int v = 0 ; v < NCOMP_FLUID ;v++) FData[FSize3D*v+i] = Cons[v];
               }
			}

		
		}
		else // Strategy == 2
		{

           IntMonoCoeff = INT_MONO_COEFF;


           do {
                 if ( ( IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) && itr == 1 ) break;

                 // interpolation
                 for (int v=0; v<NComp; v++)
                 Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                              FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, IntMonoCoeff );


                 // check
                 if ( itr <= MaxItr )
                 {

                      for ( int i = 0 ;i < FSize3D; i++ )
                      {
                         for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

#                        if ( MODEL == HYDRO )
                         if (Unphysical(Cons, Gamma_m1, CheckMinPres))
#                        elif ( MODEL ==SR_HYDRO )
                         if (SRHydro_CheckUnphysical(Cons, NULL, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, false))
#                        else
#                        error: ERROR!
#                        endif
                         {
                            state = true;
//                          reduce mono coefficient for interpolation
                            IntMonoCoeff -= ( (real)INT_MONO_COEFF - MinMonoCoeff ) / (real) MaxItr;
                            break; 
	                     } 
	                     else state = false;
                      }
                 }
				 else if ( itr == MaxItr && state )
				 {
                      IntMonoCoeff = INT_MONO_COEFF;

                      for (int i=0; i<CSize3D; i++) 
                      { // conserved --> primitive
                        for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = CData[CSize3D*v+i];

#                       if ( MODEL == HYDRO )
                        Hydro_Con2Pri( Cons, Prim, Gamma_m1, MIN_PRES, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL);
#                       elif ( MODEL == SR_HYDRO )				   
                        SRHydro_Con2Pri(Cons, Prim, GAMMA, MIN_TEMP);
#                       else
#                       error: ERROR!
#                       endif

                        for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];      
                      }
                 }

                 itr++;
  
           } while (state && itr <= MaxItr+1 );

           if( itr == MaxItr+2 || ( ( IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) && itr == 1 ) )
		   {
               for (int i=0; i<CSize3D; i++) 
               { // primitive --> conserved
                 for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = CData[CSize3D*v+i];

#                if ( MODEL == HYDRO )
                 Hydro_Pri2Con( Prim, Cons, _Gamma_m1, NormPassive_No, NULL_INT, NULL );
#                elif ( MODEL == SR_HYDRO )
                 SRHydro_Pri2Con(Prim, Cons, GAMMA);
#                else
#                error: ERROR!
#                endif

                 for (int v = 0 ; v < NCOMP_FLUID ;v++) CData[CSize3D*v+i] = Prim[v];      
               }
		   } // if( itr == MaxItr+2 || ( ( IntScheme == INT_MINMOD3D || IntScheme == INT_MINMOD1D ) && itr == 1 ) )

        } // else Strategy == 2

	 } // if ( TVar == _TOTAL )
	 else
	 {
        for (int v=0; v<NComp; v++)
        Interpolate( CData+v*CSize3D, CSize, CStart, CRange, FData+v*FSize3D,
                     FSize, FStart, 1, IntScheme, UnwrapPhase, Monotonic, INT_MONO_COEFF );
	 }




//   check minimum energy
#    ifdef CHECK_NEGATIVE_IN_FLUID
     if ( TVar == _TOTAL )
     {
        for ( int i = 0 ;i < FSize3D; i++ )
        {
           for (int v = 0 ; v < NCOMP_FLUID ;v++) Cons[v] = FData[FSize3D*v+i];

#           if ( MODEL == HYDRO )
            if (Unphysical(Cons, Gamma_m1, CheckMinPres))
#           elif ( MODEL == SR_HYDRO )
            if (SRHydro_CheckUnphysical(Cons, NULL, GAMMA, MIN_TEMP, __FUNCTION__, __LINE__, true))
#           else
#           error: ERROR!
#           endif
            {
                exit(0);
            } //  if (Unphysical(Cons, Gamma_m1, CheckMinPres))
        }
     } // if ( TVar == _TOTAL )
#    endif
}
