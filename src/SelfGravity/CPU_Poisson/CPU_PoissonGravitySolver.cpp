#include "GAMER.h"

#if ( defined GRAVITY  &&  !defined GPU )



#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];


// Poisson solver prototypes
#if   ( POT_SCHEME == SOR )
void CPU_PoissonSolver_SOR( const real Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                            const real Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                  real Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                            const int NPatchGroup, const real dh, const int Min_Iter, const int Max_Iter,
                            const real Omega, const real Poi_Coeff, const IntScheme_t IntScheme );

#elif ( POT_SCHEME == MG  )
void CPU_PoissonSolver_MG( const real Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                           const real Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                 real Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                           const int NPatchGroup, const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                           const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                           const IntScheme_t IntScheme );
#endif // POT_SCHEME


// Gravity solver prototypes
#if   ( MODEL == HYDRO )
void CPU_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const int NPatchGroup,
   const real dt, const real dh, const bool P5_Gradient,
   const OptGravityType_t GravityType, const double c_ExtAcc_AuxArray[],
   const double TimeNew, const double TimeOld, const real MinEint );

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
void CPU_ELBDMGravitySolver(       real Flu_Array[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE],
                             const real Pot_Array[][GRA_NXT][GRA_NXT][GRA_NXT],
                             const double Corner_Array[][3],
                             const int NPatchGroup, const real EtaDt, const real dh, const real Lambda,
                             const bool ExtPot, const double Time, const double ExtPot_AuxArray[] );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonGravitySolver
// Description :  Invoke the CPU_PoissonSolver and/or CPU_GravitySolver to evaluate the potential and/or
//                advance the fluid variables by the gravitational acceleration for a group of patches
//
// Parameter   :  h_Rho_Array          : Host array storing the input density
//                h_Pot_Array_In       : Host array storing the input "coarse-grid" potential for interpolation
//                h_Pot_Array_Out      : Host array to store the output potential
//                h_Flu_Array          : Host array to store the fluid variables for the Gravity solver
//                h_Corner_Array       : Host array storing the physical corner coordinates of each patch
//                h_Pot_Array_USG      : Host array storing the prepared potential for UNSPLIT_GRAVITY
//                h_Flu_Array_USG      : Host array storing the prepared density + momentum for UNSPLIT_GRAVITY
//                h_DE_Array           : Host array storing the dual-energy status (for both input and output)
//                NPatchGroup          : Number of patch groups evaluated simultaneously by GPU
//                dt                   : Time interval to advance solution
//                dh                   : Cell size
//                SOR_Min_Iter         : Minimum number of iterations for SOR
//                SOR_Max_Iter         : Maximum number of iterations for SOR
//                SOR_Omega            : Over-relaxation parameter
//                MG_Max_Iter          : Maximum number of iterations for multigrid
//                MG_NPre_Smooth       : Number of pre-smoothing steps for multigrid
//                MG_NPos_tSmooth      : Number of post-smoothing steps for multigrid
//                MG_Tolerated_Error   : Maximum tolerated error for multigrid
//                Poi_Coeff            : Coefficient in front of the RHS in the Poisson eq.
//                IntScheme            : Interpolation scheme for potential
//                                       --> currently supported schemes include
//                                           INT_CQUAD : conservative quadratic interpolation
//                                           INT_QUAD  : quadratic interpolation
//                P5_Gradient          : Use 5-points stencil to evaluate the potential gradient
//                ELBDM_Eta            : Particle mass / Planck constant in ELBDM
//                ELBDM_Lambda         : Quartic self-interaction coefficient in ELBDM
//                Poisson              : true --> invoke the Poisson solver
//                GraAcc               : true --> invoke the Gravity solver
//                GravityType          : Types of gravity --> self-gravity, external gravity, both
//                TimeNew              : Physical time at the current  step (for the external gravity solver)
//                TimeOld              : Physical time at the previous step (for the external gravity solver in UNSPLIT_GRAVITY)
//                ExtPot               : Add the external potential
//                MinEint              : Minimum allowed internal energy (== MIN_PRES / (GAMMA-1))
//
// Useless parameters in HYDRO : ELBDM_Eta, ELBDM_Lambda
// Useless parameters in ELBDM : P5_Gradient
//
// Return      :  h_Pot_Array_Out, h_Flu_Array
//-------------------------------------------------------------------------------------------------------
void CPU_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                               const real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                     real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                     real h_Flu_Array    [][GRA_NIN][PS1][PS1][PS1],
                               const double h_Corner_Array[][3],
                               const real h_Pot_Array_USG[][USG_NXT_G][USG_NXT_G][USG_NXT_G],
                               const real h_Flu_Array_USG[][GRA_NIN-1][PS1][PS1][PS1],
                                     char h_DE_Array     [][PS1][PS1][PS1],
                               const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter,
                               const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                               const int MG_NPre_Smooth, const int MG_NPost_Smooth, const real MG_Tolerated_Error,
                               const real Poi_Coeff, const IntScheme_t IntScheme, const bool P5_Gradient,
                               const real ELBDM_Eta, const real ELBDM_Lambda, const bool Poisson, const bool GraAcc,
                               const OptGravityType_t GravityType, const double TimeNew, const double TimeOld,
                               const bool ExtPot, const real MinEint )
{

// check
#  ifdef GAMER_DEBUG
   if ( Poisson  &&  ( IntScheme != INT_CQUAD  &&  IntScheme != INT_QUAD )  )
   Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );

   if ( GraAcc )
   {
      if (  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH || ExtPot )  &&  h_Corner_Array == NULL  )
      Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );

#     ifdef UNSPLIT_GRAVITY
      if (  ( GravityType == GRAVITY_SELF || GravityType == GRAVITY_BOTH )  &&  h_Pot_Array_USG == NULL  )
      Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );

      if ( h_Flu_Array_USG == NULL )
      Aux_Error( ERROR_INFO, "h_Flu_Array_USG == NULL !!\n" );
#     endif

#     ifdef DUAL_ENERGY
      if ( h_DE_Array == NULL )
      Aux_Error( ERROR_INFO, "h_DE_Array == NULL !!\n" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG


// Poisson solver
   if ( Poisson )
   {
#     if   ( POT_SCHEME == SOR )

      CPU_PoissonSolver_SOR( h_Rho_Array, h_Pot_Array_In, h_Pot_Array_Out, NPatchGroup, dh,
                             SOR_Min_Iter, SOR_Max_Iter, SOR_Omega,
                             Poi_Coeff, IntScheme );

#     elif ( POT_SCHEME == MG  )

      CPU_PoissonSolver_MG ( h_Rho_Array, h_Pot_Array_In, h_Pot_Array_Out, NPatchGroup, dh,
                             MG_Max_Iter, MG_NPre_Smooth, MG_NPost_Smooth, MG_Tolerated_Error,
                             Poi_Coeff, IntScheme );

#     else

#     error : ERROR : unsupported CPU Poisson solver !!

#     endif
   } // if ( Poisson )


// Gravity solver
   if ( GraAcc )
   {
#     if   ( MODEL == HYDRO )
      CPU_HydroGravitySolver( (real(*)[GRA_NIN][ CUBE(PS1) ])   h_Flu_Array,
                              (real(*)[ CUBE(GRA_NXT) ])        h_Pot_Array_Out,
                                                                h_Corner_Array,
                              (real(*)[ CUBE(USG_NXT_G) ])      h_Pot_Array_USG,
                              (real(*)[GRA_NIN-1][ CUBE(PS1) ]) h_Flu_Array_USG,
                              (char(*)[ CUBE(PS1) ])            h_DE_Array,
                              NPatchGroup, dt, dh, P5_Gradient, GravityType,
                              ExtAcc_AuxArray, TimeNew, TimeOld, MinEint );

#     elif ( MODEL == MHD )
#     error : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      CPU_ELBDMGravitySolver( h_Flu_Array, h_Pot_Array_Out, h_Corner_Array, NPatchGroup, ELBDM_Eta*dt, dh, ELBDM_Lambda,
                              ExtPot, TimeNew, ExtPot_AuxArray );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   } // if ( GraAcc )

} // FUNCTION : CPU_PoissonGravitySolver



#endif // #if ( defined GRAVITY  &&  !defined GPU )
