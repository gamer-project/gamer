#ifndef CPU_PROTOOTYPES_H
#define CPU_PROTOOTYPES_H


 void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                                    const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                    const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres );

 void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres );

 void CPU_Con2Pri( const real In[], real Out[], const real Gamma);

 void CPU_Pri2Con( const real In[], real Out[], const real Gamma);

 void CPU_ComputeFlux( const real FC_Var[][6][NCOMP_TOTAL], real FC_Flux[][3][NCOMP_TOTAL], const int NFlux, const int Gap,
                             const real Gamma, const bool CorrHalfVel, const real Pot_USG[], const double Corner[],
                             const real dt, const real dh, const double Time, const OptGravityType_t GravityType,
                             const double ExtAcc_AuxArray[], const real MinPres );

 void CPU_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                                const real Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                                const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                                const bool NormPassive, const int NNorm, const int NormIdx[] );

 void CPU_StoreFlux( real Flux_Array[][NCOMP_TOTAL][ PS2*PS2 ], const real FC_Flux[][3][NCOMP_TOTAL] );
#if   ( RSOLVER == EXACT )
 void CPU_RiemannSolver_Exact( const int XYZ, real eival_out[], real L_star_out[], real R_star_out[],
                                     real Flux_Out[], const real L_In[], const real R_In[], const real Gamma );
#elif ( RSOLVER == ROE )
 void CPU_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                   const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLE )
 void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLC )
 void CPU_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );
#endif
real CPU_CheckMinPres( const real InPres, const real MinPres );

void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward );

#endif

void CPU_4Velto3Vel( const real In[], real Out[] );
void CPU_3Velto4Vel( const real In[], real Out[] );

real CPU_CheckMinPresInEngy (const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                             const real Gamma_m1, const real _Gamma_m1, const real MinPres);
