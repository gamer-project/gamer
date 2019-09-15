#ifndef SRHYDRO_PROTOOTYPES_H
#define SRHYDRO_PROTOOTYPES_H


void SRHydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                                       real g_PriVar   [][ CUBE(FLU_NXT) ],
                                       real g_FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                       real g_Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                                 const int NIn, const int NGhost, const real Gamma,
                                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                 const real dt, const real dh, const real MinDens, const real MinTemp );

void SRHydro_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma, const real MinTemp );

real SRHydro_Con2Pri( const real In[], real Out[], const real Gamma, const real MinTemp );

void SRHydro_Pri2Con( const real In[], real Out[], const real Gamma);

void SRHydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                          const int Gap, const real Gamma, const real MinTemp, 
                          const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] );

void SRHydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                             const real g_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ], const real dt, const real dh,
                             const real Gamma, const real MinDens, const real MinTemp, int *state );

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

bool SRHydro_CheckUnphysical( const real Con[], const real Pri[], const real Gamma, const real MinTemp, const char s[], const int line, bool show);
real VectorDotProduct( real V1, real V2, real V3);


#endif
