#ifndef SRHYDRO_PROTOOTYPES_H
#define SRHYDRO_PROTOOTYPES_H


 void SRHydro_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL],
                              const int NIn, const int NGhost,
			      const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
			      const real EP_Coeff, const real dt, const real dh, const real MinDens,
                              const real MinTemp, int iteration );

 void SRHydro_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma );

 void SRHydro_Con2Pri( const real In[], real Out[], const real Gamma);

 void SRHydro_Pri2Con( const real In[], real Out[], const real Gamma);

 void SRHydro_ComputeFlux( const real FC_Var[][6][NCOMP_TOTAL], real FC_Flux[][3][NCOMP_TOTAL], const int NFlux, const int Gap,
                             const real Gamma, const bool CorrHalfVel, const real Pot_USG[], const double Corner[],
                             const real dt, const real dh, const double Time, const OptGravityType_t GravityType,
                             const double ExtAcc_AuxArray[], const real MinTemp );

 bool SRHydro_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                                const real Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                                const real Gamma, const real MinDens, const real MinTemp, const real DualEnergySwitch,
                                const bool NormPassive, const int NNorm, const int NormIdx[] );

 void SRHydro_StoreFlux( real Flux_Array[][NCOMP_TOTAL][ PS2*PS2 ], const real FC_Flux[][3][NCOMP_TOTAL] );
#if ( RSOLVER == HLLE )
 void SRHydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinTemp );
#elif ( RSOLVER == HLLC )
 void SRHydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinTemp );
#endif

void SRHydro_Rotate3D( real InOut[], const int XYZ, const bool Forward );


real SRHydro_CheckMinTemp (const real InTemp, const real MinTemp);

void SRHydro_4Velto3Vel( const real In[], real Out[] );
void SRHydro_3Velto4Vel( const real In[], real Out[] );

real SRHydro_CheckMinTempInEngy (const real Con[], const real MinTemp, const real Gamma);

bool SRHydro_CheckUnphysical( const real Con[], const real Pri[], const real Gamma, const char s[], const int line, bool show);

#endif
