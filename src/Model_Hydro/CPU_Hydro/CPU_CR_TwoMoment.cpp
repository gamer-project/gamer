#include "CUFLU.h"

#ifdef CR_TWOMOMENT




//-------------------------------------------------------------------------------------------------------
// Function    : RotateVec
// Description : Rotate a vector from lab frame to B-aligned frame
//
// Note        : 1. B-aligned frame: x' along B, y' and z' perpendicular to B
//               2. The rotation is defined by two angles: theta (B to z-axis) and phi (Bxy to x-axis)
//
// Parameter   : sint : sin(theta) = |Bxy|/|B|
//               cost : cos(theta) = Bz/|B|
//               sinp : sin(phi) = By/|Bxy|
//               cosp : cos(phi) = Bx/|Bxy|
//               v1, v2, v3 : input/output vector components
//
// Reference   : Athena++ src/utils/rotate_vectors.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void RotateVec( const real sint, const real cost,
                       const real sinp, const real cosp,
                       real &v1, real &v2, real &v3 )
{
// The two rotation matrices:
// R_1 (around z by phi):
// [cos_p  sin_p 0]
// [-sin_p cos_p 0]
// [0       0    1]
//
// R_2 (around y by theta-pi/2):
// [sin_t  0 cos_t]
// [0      1    0]
// [-cos_t 0 sin_t]

// First apply R1, then apply R2
   real newv1 =  cosp * v1 + sinp * v2;
   v2         = -sinp * v1 + cosp * v2;

// now apply R2
   v1         =  sint * newv1 + cost * v3;
   real newv3 = -cost * newv1 + sint * v3;
   v3 = newv3;

} // FUNCTION : RotateVec



//-------------------------------------------------------------------------------------------------------
// Function    : InvRotateVec
// Description : Rotate a vector from B-aligned frame back to lab frame (inverse of RotateVec)
//
// Note        : 1. See RotateVec()
//
// Parameter   : sint, cost, sinp, cosp : same as RotateVec
//               v1, v2, v3 : input/output vector components
//
// Reference   : Athena++ src/utils/rotate_vectors.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void InvRotateVec( const real sint, const real cost,
                          const real sinp, const real cosp,
                          real &v1, real &v2, real &v3 )
{
// R_1^-1 (around z by -phi):
// [cos_p  -sin_p 0]
// [sin_p  cos_p  0]
// [0       0     1]
//
// R_2^-1 (around y by pi/2-theta):
// [sin_t  0 -cos_t]
// [0      1    0  ]
// [cos_t  0 sin_t ]

// First apply R2^-1, then apply R1^-1
   real newv1 = sint * v1 - cost * v3;
   v3         = cost * v1 + sint * v3;

// now apply R1^-1
   v1         = cosp * newv1 - sinp * v2;
   real newv2 = sinp * newv1 + cosp * v2;
   v2 = newv2;

} // FUNCTION : InvRotateVec



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeBFieldAngles
// Description : Compute the rotation angles from B field direction
//
// Note        : 1. Same conventions as RotateVec and InvRotateVec
//
// Parameter   : Bx, By, Bz : magnetic field components
//               sint, cost, sinp, cosp : output angles (see RotateVec)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_ComputeBFieldAngles( const real Bx, const real By, const real Bz,
                                    real &sint, real &cost, real &sinp, real &cosp )
{
   const real bxby = SQRT( SQR(Bx) + SQR(By) );
   const real btot = SQRT( SQR(Bx) + SQR(By) + SQR(Bz) );

   if ( btot > TINY_NUMBER ) {
      sint = bxby / btot;
      cost = Bz / btot;
   } else {
      sint = (real)1.0;
      cost = (real)0.0;
   }

   if ( bxby > TINY_NUMBER ) {
      sinp = By / bxby;
      cosp = Bx / bxby;
   } else {
      sinp = (real)0.0;
      cosp = (real)1.0;
   }

} // FUNCTION : CR_ComputeBFieldAngles



//-------------------------------------------------------------------------------------------------------
// Function    : CR_UpdateStreaming_OneCell
// Description : Update streaming opacity (sigma_adv) and velocity (v_adv) for a single cell
//
// Note        : 1. sigma_adv is the opacity parallel to B
//               2. Invoked by CR_UpdateStreaming() and CR_UpdateOpacity()
//
// Parameter   : Ec          : CR energy density
//               rho         : gas density
//               Bx, By, Bz  : magnetic field components
//               grad_pc[3]  : gradient of CR pressure
//               vmax        : effective speed of light
//               sigma_adv   : output streaming opacity
//               v_adv[3]    : output streaming velocity components
//               MicroPhy    : Microphysics object
//
// Reference   : Athena++ src/cr/cr.cpp DefaultStreaming()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_UpdateStreaming_OneCell( const real Ec, const real rho,
                              const real Bx, const real By, const real Bz,
                              const real grad_pc[3], const real vmax,
                              real &sigma_adv, real v_adv[3], const MicroPhy_t *MicroPhy )
{
   const real invlim = (real)1.0 / vmax;
   const real bsq  = SQR(Bx) + SQR(By) + SQR(Bz);
   const real btot = SQRT( bsq );
   const real inv_sqrt_rho = (real)1.0 / SQRT( rho );
   const real va = btot * inv_sqrt_rho;  // Alfven velocity
   const real va1 = Bx * inv_sqrt_rho;
   const real va2 = By * inv_sqrt_rho;
   const real va3 = Bz * inv_sqrt_rho;
   const real b_grad_pc = Bx * grad_pc[0] + By * grad_pc[1] + Bz * grad_pc[2];

// determine sign of B dot grad(Pc)
   real dpc_sign = (real)0.0;
   if ( b_grad_pc > TINY_NUMBER )
      dpc_sign = (real)1.0;
   else if ( -b_grad_pc > TINY_NUMBER )
      dpc_sign = (real)-1.0;

// compute streaming velocity: v_adv = -sign(B dot grad Pc) * v_Alfven * b_hat
   v_adv[0] = -va1 * dpc_sign;
   v_adv[1] = -va2 * dpc_sign;
   v_adv[2] = -va3 * dpc_sign;

// compute streaming opacity (parallel to B):
// sigma_adv = |B dot grad Pc| / (|B| * v_A * (4/3) * (1/vmax) * Ec)
   if ( va > TINY_NUMBER && Ec > TINY_NUMBER ) {
      sigma_adv = FABS(b_grad_pc) / ( btot * va * ((real)4.0/(real)3.0) * invlim * Ec );
   } else {
      sigma_adv = MicroPhy->CR_max_opacity;
   }

} // FUNCTION : CR_UpdateStreaming_OneCell



//-------------------------------------------------------------------------------------------------------
// Function    : CR_UpdateStreaming
// Description : Update streaming opacity and velocity using flux divergence
//
// Note        : 1. Invoke CR_UpdateStreaming_OneCell()
//               2. Invoked by CPU_FluidSolver_MHM() for both half and full step after the flux term
//               3. grad_pc is computed from flux divergence (g_Flux)
//               4. Callers pass the same array to g_Output and g_CellVar, and the same values to NVar_*
//                  and *_offset; the separate parameters are kept for a future dedicated ADV_* buffer
//
// Parameter   : g_Output     : Array to store the updated streaming fields (ADV_SIGMA, ADV_V*)
//               g_CellVar    : Array storing input cell-centered DENS and CR_E
//               g_CC_B       : Array storing input cell-centered B field
//               g_Flux       : Array storing the input fluxes [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ]
//               NFlux        : Stride for accessing g_Flux[] (N_HF_FLUX or N_FL_FLUX)
//               NVar_Out     : Stride for accessing g_Output[] (FLU_NXT or N_HF_VAR)
//               NVar_In      : Stride for accessing g_CellVar[] (FLU_NXT or N_HF_VAR)
//               NVar_B       : Stride for accessing g_CC_B[] (FLU_NXT or N_HF_VAR)
//               out_offset   : Offset between flux cell index and output index
//               in_offset    : Offset between flux cell index and input index
//               dh           : Cell size
//               MicroPhy     : Microphysics object
//
// Return      : g_Output[ADV_*]
//
// Reference   : Athena++ src/cr/cr.cpp DefaultStreaming()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_UpdateStreaming( real g_Output[][ CUBE(FLU_NXT) ],
                       const real g_CellVar[][ CUBE(FLU_NXT) ],
                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                       const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                       const int NFlux, const int NVar_Out, const int NVar_In, const int NVar_B,
                       const int out_offset, const int in_offset,
                       const real dh, const MicroPhy_t *MicroPhy )
{
   const real vmax = MicroPhy->CR_vmax;
   const real _dh  = (real)1.0 / dh;
   const int  didx_flux[3] = { 1, NFlux, SQR(NFlux) };
   const int  CRF_v[3] = { CR_F1, CR_F2, CR_F3 };

// loop bounds
   const int cell_offset  = 1;  // skip boundary cells
   const int cell_size_i  = NFlux - 2*cell_offset;
   const int cell_size_j  = NFlux - 2*cell_offset;
   const int cell_size_k  = NFlux - 2*cell_offset;
   const int cell_size_ij = cell_size_i * cell_size_j;

   CGPU_LOOP( idx_cell, cell_size_i*cell_size_j*cell_size_k )
   {
//    flux cell indices
      const int i_cell = idx_cell % cell_size_i + cell_offset;
      const int j_cell = idx_cell % cell_size_ij / cell_size_i + cell_offset;
      const int k_cell = idx_cell / cell_size_ij + cell_offset;
      const int idx_flux = IDX321( i_cell, j_cell, k_cell, NFlux, NFlux );

//    output indices with out_offset
      const int i_out = i_cell + out_offset;
      const int j_out = j_cell + out_offset;
      const int k_out = k_cell + out_offset;
      const int idx_out = IDX321( i_out, j_out, k_out, NVar_Out, NVar_Out );

//    input indices with in_offset
      const int i_in = i_cell + in_offset;
      const int j_in = j_cell + in_offset;
      const int k_in = k_cell + in_offset;
      const int idx_in = IDX321( i_in, j_in, k_in, NVar_In, NVar_In );

//    B field index
      const int idx_B = IDX321( i_in, j_in, k_in, NVar_B, NVar_B );

//    compute grad_pc from flux divergence
      real grad_pc[3];
      for (int n=0; n<3; n++) {
         grad_pc[n] = (real)0.0;
         for (int d=0; d<3; d++)
            grad_pc[n] += g_Flux[d][ CRF_v[n] ][idx_flux] - g_Flux[d][ CRF_v[n] ][ idx_flux - didx_flux[d] ];
         grad_pc[n] *= _dh / vmax;
      }

      const real Ec  = g_CellVar[CR_E ][idx_in];
      const real rho = g_CellVar[DENS][idx_in];
      const real Bx  = g_CC_B[0][idx_B];
      const real By  = g_CC_B[1][idx_B];
      const real Bz  = g_CC_B[2][idx_B];

//    compute new sigma_adv and v_adv
      real sigma_adv_new;
      real v_adv_new[3];
      CR_UpdateStreaming_OneCell( Ec, rho, Bx, By, Bz, grad_pc, vmax, sigma_adv_new, v_adv_new, MicroPhy );

      g_Output[ADV_SIGMA][idx_out] = sigma_adv_new;
      g_Output[ADV_VX   ][idx_out] = v_adv_new[0];
      g_Output[ADV_VY   ][idx_out] = v_adv_new[1];
      g_Output[ADV_VZ   ][idx_out] = v_adv_new[2];

   } // CGPU_LOOP

} // FUNCTION : CR_UpdateStreaming



//-------------------------------------------------------------------------------------------------------
// Function    : CR_UpdateOpacity
// Description : Update streaming opacity and velocity using the central-difference CR pressure gradient
//
// Note        : 1. Invoke CR_UpdateStreaming_OneCell()
//               2. Invoked by CPU_FluidSolver_MHM() before the half-step flux and by Hydro_RiemannPredict()
//                  after the half-step source term (for the full-step)
//                  --> also invoked by CPU_FluidSolver_MHM() to restore ADV_* before a retry
//               3. grad_pc is computed with central differences of CR pressure
//               4. g_Output[] may be conserved (g_Flu_Array_In) or primitive (g_PriVar_Half): safe
//                  because only DENS, CR_E, and ADV_* are accessed
//
// Parameter   : g_Output     : Array to store the updated streaming fields (ADV_SIGMA, ADV_V*)
//                              --> also storing input cell-centered DENS and CR_E
//               g_CC_B       : Array storing input cell-centered B field
//               NVar_Out     : Stride for accessing g_Output[] (FLU_NXT or N_HF_VAR)
//               NVar_In      : Stride for accessing g_Output[] (FLU_NXT or N_HF_VAR)
//               NVar_B       : Stride for accessing g_CC_B[] (FLU_NXT or N_HF_VAR)
//               out_offset   : Offset between loop index and output index
//               in_offset    : Offset between loop index and input index
//               NSize        : Number of cells to update along each direction
//               dh           : Cell size
//               MicroPhy     : Microphysics object
//
// Return      : g_Output[ADV_*]
//
// Reference   : Athena++ src/cr/cr.cpp DefaultOpacity()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_UpdateOpacity( real g_Output[][ CUBE(FLU_NXT) ],
                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                       const int NVar_Out, const int NVar_In, const int NVar_B,
                       const int out_offset, const int in_offset, const int NSize,
                       const real dh, const MicroPhy_t *MicroPhy )
{
   const real vmax   = MicroPhy->CR_vmax;
   const real _2dh   = (real)0.5 / dh;

   const int size_ij = NSize * NSize;

   CGPU_LOOP( idx_cell, NSize*NSize*NSize )
   {
//    cell indices
      const int i_cell = idx_cell % NSize;
      const int j_cell = idx_cell % size_ij / NSize;
      const int k_cell = idx_cell / size_ij;

//    output indices
      const int i_out = i_cell + out_offset;
      const int j_out = j_cell + out_offset;
      const int k_out = k_cell + out_offset;
      const int idx_out = IDX321( i_out, j_out, k_out, NVar_Out, NVar_Out );

//    input indices for central difference
      const int i_in = i_cell + in_offset;
      const int j_in = j_cell + in_offset;
      const int k_in = k_cell + in_offset;
      const int idx_in = IDX321( i_in, j_in, k_in, NVar_In, NVar_In );

//    neighbor indices for central difference
      const int idx_ip1 = IDX321( i_in+1, j_in,   k_in,   NVar_In, NVar_In );
      const int idx_im1 = IDX321( i_in-1, j_in,   k_in,   NVar_In, NVar_In );
      const int idx_jp1 = IDX321( i_in,   j_in+1, k_in,   NVar_In, NVar_In );
      const int idx_jm1 = IDX321( i_in,   j_in-1, k_in,   NVar_In, NVar_In );
      const int idx_kp1 = IDX321( i_in,   j_in,   k_in+1, NVar_In, NVar_In );
      const int idx_km1 = IDX321( i_in,   j_in,   k_in-1, NVar_In, NVar_In );

//    B field index
      const int idx_B = IDX321( i_in, j_in, k_in, NVar_B, NVar_B );

//    compute grad(Pc) using central differences
//    Pc = Ec / 3, so grad(Pc) = (1/3) * grad(Ec)
      real grad_pc[3];
      grad_pc[0] = ( g_Output[CR_E][idx_ip1] - g_Output[CR_E][idx_im1] ) / (real)3.0 * _2dh;
      grad_pc[1] = ( g_Output[CR_E][idx_jp1] - g_Output[CR_E][idx_jm1] ) / (real)3.0 * _2dh;
      grad_pc[2] = ( g_Output[CR_E][idx_kp1] - g_Output[CR_E][idx_km1] ) / (real)3.0 * _2dh;

//    get cell-centered values
      const real Ec  = g_Output[CR_E ][idx_in];
      const real rho = g_Output[DENS][idx_in];
      const real Bx  = g_CC_B[0][idx_B];
      const real By  = g_CC_B[1][idx_B];
      const real Bz  = g_CC_B[2][idx_B];

//    compute new sigma_adv and v_adv
      real sigma_adv;
      real v_adv[3];
      CR_UpdateStreaming_OneCell( Ec, rho, Bx, By, Bz, grad_pc, vmax, sigma_adv, v_adv, MicroPhy );

      g_Output[ADV_SIGMA][idx_out] = sigma_adv;
      g_Output[ADV_VX   ][idx_out] = v_adv[0];
      g_Output[ADV_VY   ][idx_out] = v_adv[1];
      g_Output[ADV_VZ   ][idx_out] = v_adv[2];

   } // CGPU_LOOP

} // FUNCTION : CR_UpdateOpacity



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeVdiff
// Description : Compute the cell-centered diffusion velocity (v_diff) along a given direction
//
// Note        : 1. Invoked by CR_TwoMomentFlux_HalfStep() and CR_TwoMomentFlux_FullStep()
//
// Parameter   : sigma_adv   : streaming opacity (parallel to B)
//               Bx, By, Bz  : cell-centered magnetic field components
//               Ec          : CR energy density
//               rho         : gas density
//               vmax        : effective speed of light
//               dh          : cell size
//               fdir        : flux direction (0/1/2 = x/y/z)
//               MicroPhy    : Microphysics object
//
// Return      : v_diff along the specified direction (fdir)
//
// Reference   : Athena++ src/cr/integrators/cr_transport.cpp CRIntegrator::CalculateFluxes()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static real CR_ComputeVdiff( const real sigma_adv,
                             const real Bx, const real By, const real Bz,
                             const real Ec, const real rho,
                             const real vmax, const real dh, const int fdir,
                             const MicroPhy_t *MicroPhy )
{
   const real edd = (real)1.0 / (real)3.0;
   const real sigma_adv_perp = MicroPhy->CR_max_opacity;
   const real sigma_diff = MicroPhy->CR_sigma;
   const real sigma_diff_perp = MicroPhy->CR_sigma_perp;

// compute B-field angles for rotation
   real sint, cost, sinp, cosp;
   CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

// In B-aligned frame: x parallel to B
   real sigma_x = sigma_diff;
   real sigma_y = sigma_diff_perp;
   real sigma_z = sigma_diff_perp;
// if streaming is enabled: total sigma = combine sigma_diff with sigma_adv
   if ( MicroPhy->CR_stream ) {
      sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv );
      sigma_y = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
      sigma_z = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
   }

// compute tau and diffv for each B-aligned direction
   real tau_x = MicroPhy->CR_taufact * sigma_x * dh;
   tau_x = tau_x * tau_x / ( (real)2.0 * edd );
   real diffv_x;
   if ( tau_x < MicroPhy->CR_tau_asym_lim )
      diffv_x = SQRT( (real)1.0 - (real)0.5 * tau_x );
   else
      diffv_x = SQRT( ( (real)1.0 - EXP(-tau_x) ) / tau_x );

   real tau_y = MicroPhy->CR_taufact * sigma_y * dh;
   tau_y = tau_y * tau_y / ( (real)2.0 * edd );
   real diffv_y;
   if ( tau_y < MicroPhy->CR_tau_asym_lim )
      diffv_y = SQRT( (real)1.0 - (real)0.5 * tau_y );
   else
      diffv_y = SQRT( ( (real)1.0 - EXP(-tau_y) ) / tau_y );

   real tau_z = MicroPhy->CR_taufact * sigma_z * dh;
   tau_z = tau_z * tau_z / ( (real)2.0 * edd );
   real diffv_z;
   if ( tau_z < MicroPhy->CR_tau_asym_lim )
      diffv_z = SQRT( (real)1.0 - (real)0.5 * tau_z );
   else
      diffv_z = SQRT( ( (real)1.0 - EXP(-tau_z) ) / tau_z );

   real vdiff_1 = vmax * SQRT(edd) * diffv_x;
   real vdiff_2 = vmax * SQRT(edd) * diffv_y;
   real vdiff_3 = vmax * SQRT(edd) * diffv_z;

// rotate vdiff_1/2/3 from B-aligned frame to lab frame
   InvRotateVec( sint, cost, sinp, cosp, vdiff_1, vdiff_2, vdiff_3 );

   vdiff_1 = FABS( vdiff_1 );
   vdiff_2 = FABS( vdiff_2 );
   vdiff_3 = FABS( vdiff_3 );

// add CR sound speed if CR_VEL_FLX_FLAG is on
   const real cr_sound = MicroPhy->CR_vel_flx_flag * SQRT( ((real)4.0/(real)9.0) * Ec / rho );
   vdiff_1 += cr_sound;
   vdiff_2 += cr_sound;
   vdiff_3 += cr_sound;

// return component along flux direction
   if ( fdir == 0 )      return vdiff_1;
   else if ( fdir == 1 ) return vdiff_2;
   else                  return vdiff_3;

} // FUNCTION : CR_ComputeVdiff



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeHLLEFlux
// Description : Compute HLLE flux for CR variables at a cell interface
//
// Parameter   : Ec_L/R      : CR energy density (left/right)
//               Fc_L/R[3]   : reduced CR flux Fc/vmax (left/right)
//               v_L/R       : gas velocity along flux direction (left/right)
//               vdiff_L/R   : diffusion velocity (left/right)
//               vmax        : effective speed of light
//               fdir        : flux direction (0/1/2 = x/y/z)
//               flux_E      : output flux for CR energy
//               flux_F[3]   : output flux for CR flux components
//
// Reference   : Athena++ src/cr/integrators/cr_flux.cpp CRIntegrator::CRFlux()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_ComputeHLLEFlux( const real Ec_L, const real Ec_R,
                                const real Fc_L[3], const real Fc_R[3],
                                const real v_L, const real v_R,
                                const real vdiff_L, const real vdiff_R,
                                const real vmax,
                                const int fdir,
                                real &flux_E, real flux_F[3] )
{
// Eddington factors (isotropic)
   const real edd   = (real)1.0 / (real)3.0;   // diagonal component
   const real edd12 = (real)0.0;               // off-diagonal components
   const real edd13 = (real)0.0;
   const real edd23 = (real)0.0;

// compute wave speeds
   const real meanadv  = (real)0.5 * ( v_L + v_R );
   const real meandiffv = (real)0.5 * ( vdiff_L + vdiff_R );

   real al = FMIN( meanadv - meandiffv, v_L - vdiff_L );
   real ar = FMAX( meanadv + meandiffv, v_R + vdiff_R );

// clamp to the max CR signal speed vmax*sqrt(edd)
   const real sqrt_edd = SQRT( edd );
   ar = FMIN( ar,  vmax * sqrt_edd );
   al = FMAX( al, -vmax * sqrt_edd );

   const real bp = FMAX( ar, (real)0.0 );   // positive wave speed
   const real bm = FMIN( al, (real)0.0 );   // negative wave speed

// compute L/R fluxes minus wave transport: F_L - bm*U_L, F_R - bp*U_R
   real fl_e, fr_e;
   real fl_f1, fr_f1, fl_f2, fr_f2, fl_f3, fr_f3;

// energy equation flux: vmax * F_c - b * E_c
   fl_e = vmax * Fc_L[fdir] - bm * Ec_L;
   fr_e = vmax * Fc_R[fdir] - bp * Ec_R;

// flux equation fluxes depend on direction
   if ( fdir == 0 ) {       // x-direction (CRF1)
      fl_f1 = vmax * edd   * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd   * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd12 * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd12 * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd13 * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd13 * Ec_R - bp * Fc_R[2];
   }
   else if ( fdir == 1 ) {  // y-direction (CRF2)
      fl_f1 = vmax * edd12 * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd12 * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd   * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd   * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd23 * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd23 * Ec_R - bp * Fc_R[2];
   }
   else {                   // z-direction (CRF3)
      fl_f1 = vmax * edd13 * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd13 * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd23 * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd23 * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd   * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd   * Ec_R - bp * Fc_R[2];
   }

// HLLE flux formula
   real tmp = (real)0.0;
   if ( FABS(bm - bp) > TINY_NUMBER )
      tmp = (real)0.5 * (bp + bm) / (bp - bm);

   flux_E    = (real)0.5 * (fl_e  + fr_e ) + (fl_e  - fr_e ) * tmp;
   flux_F[0] = (real)0.5 * (fl_f1 + fr_f1) + (fl_f1 - fr_f1) * tmp;
   flux_F[1] = (real)0.5 * (fl_f2 + fr_f2) + (fl_f2 - fr_f2) * tmp;
   flux_F[2] = (real)0.5 * (fl_f3 + fr_f3) + (fl_f3 - fr_f3) * tmp;

} // FUNCTION : CR_ComputeHLLEFlux



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_HalfStep
//
// Description : Compute the half-step CR two-moment HLLE fluxes
//
// Note        : 1. Donor-cell reconstruction, since Athena++ forces it at VL2 stage 1
//               2. vdiff is evaluated at cell centers with the stored ADV_SIGMA
//               3. Overwrites g_Flux_Half[CR_*] computed by Hydro_RiemannPredict_Flux()
//                  and zeroes g_Flux_Half[ADV_*] (recomputed each step, never advected)
//
// Parameter   : g_ConVar    : Array storing the input cell-centered conserved variables
//               g_Flux_Half : Array storing the face-centered fluxes (CR_* and ADV_* slots are overwritten)
//               g_CC_B      : Array storing the input cell-centered B field
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[CR_*/ADV_*]
//
// Reference   : Athena++ src/cr/integrators/cr_transport.cpp, cr_flux.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentFlux_HalfStep( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                      real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_CC_B[][ CUBE(FLU_NXT) ],
                                const real dh, const MicroPhy_t *MicroPhy )
{
   const int  didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int  flux_offset  = 1;  // skip the extra transverse fluxes used only by the CT electric field

   for (int d=0; d<3; d++)
   {
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i   = N_HF_FLUX-1;              size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = 0;                        j_offset = flux_offset;              k_offset    = flux_offset;
                  break;

         case 1 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-1;              size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = flux_offset;              j_offset = 0;                        k_offset    = flux_offset;
                  break;

         case 2 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-1;
                  i_offset = flux_offset;              j_offset = flux_offset;              k_offset    = 0;
                  break;
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux   = idx % size_i           + i_offset;
         const int j_flux   = idx % size_ij / size_i + j_offset;
         const int k_flux   = idx / size_ij          + k_offset;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

//       conserved variable and cell-centered magnetic field index
         const int i_cvar   = i_flux;
         const int j_cvar   = j_flux;
         const int k_cvar   = k_flux;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       get left and right cell indices
         const int idx_L = idx_cvar;
         const int idx_R = idx_cvar + didx_cvar[d];

//       1. get CR variables (donor-cell data reconstruction)
         const real Ec_L   = g_ConVar[CR_E ][idx_L];
         const real Ec_R   = g_ConVar[CR_E ][idx_R];
         const real Fc_L[3] = { g_ConVar[CR_F1][idx_L], g_ConVar[CR_F2][idx_L], g_ConVar[CR_F3][idx_L] };
         const real Fc_R[3] = { g_ConVar[CR_F1][idx_R], g_ConVar[CR_F2][idx_R], g_ConVar[CR_F3][idx_R] };

//       2. get density and velocity
         const real rho_L = g_ConVar[DENS][idx_L];
         const real rho_R = g_ConVar[DENS][idx_R];
         const real v_L   = g_ConVar[MOMX+d][idx_L] / rho_L;
         const real v_R   = g_ConVar[MOMX+d][idx_R] / rho_R;

//       3. compute vdiff
//       3-1. get cell-centered magnetic field for left and right cells
         const real Bx_L = g_CC_B[0][idx_L];
         const real By_L = g_CC_B[1][idx_L];
         const real Bz_L = g_CC_B[2][idx_L];

         const real Bx_R = g_CC_B[0][idx_R];
         const real By_R = g_CC_B[1][idx_R];
         const real Bz_R = g_CC_B[2][idx_R];

//       3-2. get sigma_adv (set by CR_UpdateOpacity())
         const real sigma_adv_L = g_ConVar[ADV_SIGMA][idx_L];
         const real sigma_adv_R = g_ConVar[ADV_SIGMA][idx_R];

//       3-3. compute vdiff for the left/right cells
         const real vmax = MicroPhy->CR_vmax;  // Should be optimized
         const real vdiff_L = CR_ComputeVdiff( sigma_adv_L, Bx_L, By_L, Bz_L, Ec_L, rho_L, vmax, dh, d, MicroPhy );
         const real vdiff_R = CR_ComputeVdiff( sigma_adv_R, Bx_R, By_R, Bz_R, Ec_R, rho_R, vmax, dh, d, MicroPhy );

//       4. compute HLLE flux
         real flux_E, flux_F[3];
         CR_ComputeHLLEFlux( Ec_L, Ec_R, Fc_L, Fc_R, v_L, v_R, vdiff_L, vdiff_R, vmax, d,
                             flux_E, flux_F );

//       5. overwrite the CR fluxes from the hydro solver
         g_Flux_Half[d][CR_E ][idx_flux] = flux_E;
         g_Flux_Half[d][CR_F1][idx_flux] = flux_F[0];
         g_Flux_Half[d][CR_F2][idx_flux] = flux_F[1];
         g_Flux_Half[d][CR_F3][idx_flux] = flux_F[2];

//       6. zero the ADV_* fluxes so ADV_* are not advected
         g_Flux_Half[d][ADV_SIGMA][idx_flux] = (real)0.0;
         g_Flux_Half[d][ADV_VX   ][idx_flux] = (real)0.0;
         g_Flux_Half[d][ADV_VY   ][idx_flux] = (real)0.0;
         g_Flux_Half[d][ADV_VZ   ][idx_flux] = (real)0.0;

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)


} // FUNCTION : CR_TwoMomentFlux_HalfStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_FullStep
//
// Description : Compute the full-step CR two-moment HLLE fluxes
//
// Note        : 1. Fluxes are computed from reconstructed left/right states
//               2. vdiff is computed from cell-centered g_PriVar_Half
//               3. Overwrites g_FC_Flux[CR_*] computed by Hydro_ComputeFlux()
//                  and zeroes g_FC_Flux[ADV_*] (recomputed each step, never advected)
//
// Parameter   : g_FC_Var      : Array storing the input reconstructed face-centered conserved variables
//               g_PriVar_Half : Array storing the input cell-centered half-step variables
//               g_FC_Flux     : Array storing the full-step fluxes (CR_* and ADV_* slots are overwritten)
//               NFlux         : Stride for accessing g_FC_Flux[]
//               NSkip_N       : Number of fluxes to skip at the beginning and end in the normal direction
//               NSkip_T       : Number of fluxes to skip at the beginning and end in the transverse direction
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_FC_Flux[CR_*/ADV_*]
//
// Reference   : Athena++ src/cr/integrators/cr_transport.cpp, cr_flux.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentFlux_FullStep( const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                 const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                       real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 const int NFlux, const int NSkip_N, const int NSkip_T,
                                 const real dh, const MicroPhy_t *MicroPhy )
{
   const int  didx_fc[3]   = { 1, N_FC_VAR, SQR(N_FC_VAR) };
   const int  didx_pvar[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };

// offset from g_FC_Var index to g_PriVar_Half index
   const int fc2pvar_offset = ( N_HF_VAR - N_FC_VAR ) / 2;

   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;        // left face index
      const int faceR = faceL + 1;  // right face index
      int idx_fc_s[3], idx_flux_e[3];

      switch ( d )
      {
         case 0 : idx_fc_s  [0] = NSkip_N;              idx_fc_s  [1] = NSkip_T;              idx_fc_s  [2] = NSkip_T;
                  idx_flux_e[0] = N_FC_VAR-1-2*NSkip_N; idx_flux_e[1] = N_FC_VAR-2*NSkip_T;   idx_flux_e[2] = N_FC_VAR-2*NSkip_T;
                  break;

         case 1 : idx_fc_s  [0] = NSkip_T;              idx_fc_s  [1] = NSkip_N;              idx_fc_s  [2] = NSkip_T;
                  idx_flux_e[0] = N_FC_VAR-2*NSkip_T;   idx_flux_e[1] = N_FC_VAR-1-2*NSkip_N; idx_flux_e[2] = N_FC_VAR-2*NSkip_T;
                  break;

         case 2 : idx_fc_s  [0] = NSkip_T;              idx_fc_s  [1] = NSkip_T;              idx_fc_s  [2] = NSkip_N;
                  idx_flux_e[0] = N_FC_VAR-2*NSkip_T;   idx_flux_e[1] = N_FC_VAR-2*NSkip_T;   idx_flux_e[2] = N_FC_VAR-1-2*NSkip_N;
                  break;
      } // switch ( d )

      const int size_ij = idx_flux_e[0]*idx_flux_e[1];

      CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
      {
//       flux index
         const int i_flux   = idx % idx_flux_e[0];
         const int j_flux   = idx % size_ij / idx_flux_e[0];
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

//       face-centered variable index for accessing g_FC_Var
         const int i_fc     = i_flux + idx_fc_s[0];
         const int j_fc     = j_flux + idx_fc_s[1];
         const int k_fc     = k_flux + idx_fc_s[2];
         const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR );

//       get left and right cell indices for g_PriVar_Half to calculate vdiff
         const int i_pvar   = i_fc + fc2pvar_offset;
         const int j_pvar   = j_fc + fc2pvar_offset;
         const int k_pvar   = k_fc + fc2pvar_offset;
         const int idx_pvar = IDX321( i_pvar, j_pvar, k_pvar, N_HF_VAR, N_HF_VAR );
         const int idx_pvar_L = idx_pvar;
         const int idx_pvar_R = idx_pvar + didx_pvar[d];

//       1. get the reconstructed CR variables
//          --> L/R state = faceR of the left cell / faceL of the right cell
         const real Ec_L   = g_FC_Var[faceR][CR_E ][idx_fc];
         const real Ec_R   = g_FC_Var[faceL][CR_E ][idx_fc + didx_fc[d]];
         const real Fc_L[3] = { g_FC_Var[faceR][CR_F1][idx_fc], g_FC_Var[faceR][CR_F2][idx_fc], g_FC_Var[faceR][CR_F3][idx_fc] };
         const real Fc_R[3] = { g_FC_Var[faceL][CR_F1][idx_fc + didx_fc[d]], g_FC_Var[faceL][CR_F2][idx_fc + didx_fc[d]], g_FC_Var[faceL][CR_F3][idx_fc + didx_fc[d]] };

//       2. get the reconstructed density and velocity
         const real rho_L = g_FC_Var[faceR][DENS][idx_fc];
         const real rho_R = g_FC_Var[faceL][DENS][idx_fc + didx_fc[d]];
         const real v_L   = g_FC_Var[faceR][MOMX+d][idx_fc] / rho_L;
         const real v_R   = g_FC_Var[faceL][MOMX+d][idx_fc + didx_fc[d]] / rho_R;

//       3. compute vdiff with cell-centered variables (donor-cell)
//       3-1. get the variables from Hydro_RiemannPredict()
         const real Bx_L = g_PriVar_Half[MAG_OFFSET+0][idx_pvar_L];
         const real By_L = g_PriVar_Half[MAG_OFFSET+1][idx_pvar_L];
         const real Bz_L = g_PriVar_Half[MAG_OFFSET+2][idx_pvar_L];

         const real Bx_R = g_PriVar_Half[MAG_OFFSET+0][idx_pvar_R];
         const real By_R = g_PriVar_Half[MAG_OFFSET+1][idx_pvar_R];
         const real Bz_R = g_PriVar_Half[MAG_OFFSET+2][idx_pvar_R];

         const real Ec_cc_L  = g_PriVar_Half[CR_E][idx_pvar_L];
         const real Ec_cc_R  = g_PriVar_Half[CR_E][idx_pvar_R];
         const real rho_cc_L = g_PriVar_Half[DENS][idx_pvar_L];
         const real rho_cc_R = g_PriVar_Half[DENS][idx_pvar_R];

//       3-2. get sigma_adv (set by CR_UpdateOpacity() after the half-step source)
         const real sigma_adv_L = g_PriVar_Half[ADV_SIGMA][idx_pvar_L];
         const real sigma_adv_R = g_PriVar_Half[ADV_SIGMA][idx_pvar_R];

//       3-3. compute vdiff for the left/right cells
         const real vmax = MicroPhy->CR_vmax;  // Should be optimized
         const real vdiff_L = CR_ComputeVdiff( sigma_adv_L, Bx_L, By_L, Bz_L, Ec_cc_L, rho_cc_L, vmax, dh, d, MicroPhy );
         const real vdiff_R = CR_ComputeVdiff( sigma_adv_R, Bx_R, By_R, Bz_R, Ec_cc_R, rho_cc_R, vmax, dh, d, MicroPhy );

//       4. compute HLLE flux
         real flux_E, flux_F[3];
         CR_ComputeHLLEFlux( Ec_L, Ec_R, Fc_L, Fc_R, v_L, v_R, vdiff_L, vdiff_R, vmax, d,
                             flux_E, flux_F );

//       5. overwrite the CR fluxes from the hydro solver
         g_FC_Flux[d][CR_E ][idx_flux] = flux_E;
         g_FC_Flux[d][CR_F1][idx_flux] = flux_F[0];
         g_FC_Flux[d][CR_F2][idx_flux] = flux_F[1];
         g_FC_Flux[d][CR_F3][idx_flux] = flux_F[2];

//       6. zero the ADV_* fluxes so ADV_* are not advected
         g_FC_Flux[d][ADV_SIGMA][idx_flux] = (real)0.0;
         g_FC_Flux[d][ADV_VX   ][idx_flux] = (real)0.0;
         g_FC_Flux[d][ADV_VY   ][idx_flux] = (real)0.0;
         g_FC_Flux[d][ADV_VZ   ][idx_flux] = (real)0.0;

      } // CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
   } // for (int d=0; d<3; d++)


} // FUNCTION : CR_TwoMomentFlux_FullStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_HalfStep
//
// Description : Apply implicit source term update for CR two-moment equations at half-step
//
// Note        : 1. Invoked by Hydro_RiemannPredict() after the flux-divergence update
//               2. Gas rho/momentum are taken from the post-flux OneCell[]
//               3. Applies the gas back-reaction when CR_source is on
//               4. The ec_source term uses grad(Pc) from g_Flux_Half, velocity from t^n g_ConVar_In
//                  (= Athena++ stage-1 pre-stage w) just to match Athena++
//               5. Rotation uses the t^n B_n[], not the half-step OneCell[MAG_OFFSET+*]
//                  (= Athena++ stage-1 b_angle) just to match Athena++
//
// Parameter   : OneCell     : Single-cell fluid array to update
//               g_ConVar_In : Array storing the input cell-centered conserved variables
//               g_Flux_Half : Array storing the input face-centered fluxes
//               B_n         : t^n cell-centered B field [3] (see Note 5)
//               idx_in      : Index for accessing g_ConVar_In[]
//               idx_flux    : Index for accessing g_Flux_Half[]
//               didx_flux   : Index increment of g_Flux_Half[]
//               dt          : Full time step
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : OneCell[] (modified CR and gas fields)
//
// Reference   : Athena++ src/cr/integrators/cr_source.cpp, cr_transport.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_HalfStep( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                            const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                            const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real B_n[],
                            const int idx_in, const int idx_flux, const int didx_flux[3],
                            const real dt, const real dh, const MicroPhy_t *MicroPhy )
{
// CR parameters
   const real vmax   = MicroPhy->CR_vmax;
   const real invlim = (real)1.0 / vmax;
   const real dt_source = (real)0.5 * dt;

   const bool CR_source    = MicroPhy->CR_source;      // flag to enable back-reaction to gas
   const bool CR_stream    = MicroPhy->CR_stream;      // flag to enable streaming
   const bool CR_Ec_source = MicroPhy->CR_Ec_source;   // flag to include the CR energy source term

// 1. get the CR state
//    --> floor CR_E here since the passive-scalar floor in Hydro_RiemannPredict() is applied only after this function
   real ec  = FMAX( OneCell[CR_E], TINY_NUMBER );
   real fc1 = OneCell[CR_F1];
   real fc2 = OneCell[CR_F2];
   real fc3 = OneCell[CR_F3];

// 2. get gas density and velocity
//    --> floor rho here since the density floor is applied after this function
   const real rho = FMAX( OneCell[DENS], TINY_NUMBER );
   real v1 = OneCell[MOMX] / rho;
   real v2 = OneCell[MOMY] / rho;
   real v3 = OneCell[MOMZ] / rho;

// 3. get the t^n B field (see Note 5)
   const real Bx = B_n[MAGX];
   const real By = B_n[MAGY];
   const real Bz = B_n[MAGZ];

// 4. compute B-field angles for rotation
   real sint, cost, sinp, cosp;
   CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

// 5. compute the total velocity
//    --> ADV_* in g_ConVar_In[] are set by CR_UpdateStreaming()
   const real v_adv_x = g_ConVar_In[ADV_VX][idx_in];
   const real v_adv_y = g_ConVar_In[ADV_VY][idx_in];
   const real v_adv_z = g_ConVar_In[ADV_VZ][idx_in];

   real vtot1 = v1;
   real vtot2 = v2;
   real vtot3 = v3;
   if ( CR_stream ) {
      vtot1 += v_adv_x;
      vtot2 += v_adv_y;
      vtot3 += v_adv_z;
   }

// 6. save original CR energy and flux for the floor check and gas back-reaction
   const real ec_old  = ec;
   const real fc1_old = fc1;
   const real fc2_old = fc2;
   const real fc3_old = fc3;

// 7. rotate all vectors to B-aligned frame
   real fr1 = fc1, fr2 = fc2, fr3 = fc3;

#  ifdef MHD
   RotateVec( sint, cost, sinp, cosp, v1, v2, v3 );
   RotateVec( sint, cost, sinp, cosp, fr1, fr2, fr3 );
   RotateVec( sint, cost, sinp, cosp, vtot1, vtot2, vtot3 );

// perpendicular work is added via ec_source (step 11)
// the vtot2/vtot3 terms below are kept for generality
   vtot2 = (real)0.0;
   vtot3 = (real)0.0;
#  endif

// 8. compute effective opacity
   const real sigma_diff = MicroPhy->CR_sigma;
   const real sigma_diff_perp = MicroPhy->CR_sigma_perp;
   const real sigma_adv_para = g_ConVar_In[ADV_SIGMA][idx_in];
   const real sigma_adv_perp = MicroPhy->CR_max_opacity;

   real sigma_x = sigma_diff;
   real sigma_y = sigma_diff_perp;
   real sigma_z = sigma_diff_perp;
   if ( CR_stream ) {
      sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
      sigma_y = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
      sigma_z = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
   }

// 9. build implicit matrix and solve
//     Source terms:
//     dEc/dt = -vtot * sigma * (Fc - v*Ec*(4/3)/vmax)
//     dFc/dt = -vmax * sigma * (Fc - v*Ec*(4/3)/vmax)

   const real rhs1 = ec;
   const real rhs2 = fr1;
   const real rhs3 = fr2;
   const real rhs4 = fr3;

// coefficients for Ec equation (row 1)
   const real coef_11 = (real)1.0 - dt_source * sigma_x * vtot1 * v1 * invlim * ((real)4.0/(real)3.0)
                                  - dt_source * sigma_y * vtot2 * v2 * invlim * ((real)4.0/(real)3.0)
                                  - dt_source * sigma_z * vtot3 * v3 * invlim * ((real)4.0/(real)3.0);
   const real coef_12 = dt_source * sigma_x * vtot1;
   const real coef_13 = dt_source * sigma_y * vtot2;
   const real coef_14 = dt_source * sigma_z * vtot3;

// coefficients for Fr1 equation (row 2)
   const real coef_21 = -dt_source * v1 * sigma_x * ((real)4.0/(real)3.0);
   const real coef_22 = (real)1.0 + dt_source * vmax * sigma_x;

// coefficients for Fr2 equation (row 3)
   const real coef_31 = -dt_source * v2 * sigma_y * ((real)4.0/(real)3.0);
   const real coef_33 = (real)1.0 + dt_source * vmax * sigma_y;

// coefficients for Fr3 equation (row 4)
   const real coef_41 = -dt_source * v3 * sigma_z * ((real)4.0/(real)3.0);
   const real coef_44 = (real)1.0 + dt_source * vmax * sigma_z;

// solve by substitution
   const real e_coef = coef_11 - coef_12 * coef_21 / coef_22
                               - coef_13 * coef_31 / coef_33
                               - coef_14 * coef_41 / coef_44;

   real new_ec = rhs1 - coef_12 * rhs2 / coef_22
                      - coef_13 * rhs3 / coef_33
                      - coef_14 * rhs4 / coef_44;
   new_ec /= e_coef;

// CR_Ec_source=0 drops the Ec source term --> dEc/dt = -div(Fc) (Jiang & Oh 2018 Eq. 19)
   if ( !CR_Ec_source )   new_ec = rhs1;

// back-substitute to get new flux
   real newfr1 = ( rhs2 - coef_21 * new_ec ) / coef_22;
   real newfr2 = ( rhs3 - coef_31 * new_ec ) / coef_33;
   real newfr3 = ( rhs4 - coef_41 * new_ec ) / coef_44;

// 10. rotate back to lab frame
#  ifdef MHD
   InvRotateVec( sint, cost, sinp, cosp, newfr1, newfr2, newfr3 );
#  endif

// 11. compute perpendicular work ec_source = v_perp*grad(Pc)
// --> grad(Pc)[n] = div(flux of CR_F(n+1))/vmax, as in CR_UpdateStreaming()
#  ifdef MHD
   const real _dh = (real)1.0 / dh;
   const int  CRF_v[3] = { CR_F1, CR_F2, CR_F3 };
   real grad_pc[3];
   for (int n=0; n<3; n++) {
      grad_pc[n] = (real)0.0;
      for (int d=0; d<3; d++)
         grad_pc[n] += g_Flux_Half[d][ CRF_v[n] ][idx_flux] - g_Flux_Half[d][ CRF_v[n] ][ idx_flux - didx_flux[d] ];
      grad_pc[n] *= _dh * invlim;
   }

   real dpcdx_B = grad_pc[0], dpcdy_B = grad_pc[1], dpcdz_B = grad_pc[2];
   RotateVec( sint, cost, sinp, cosp, dpcdx_B, dpcdy_B, dpcdz_B );

// compute perpendicular gas velocity (t^n, see Note 4)
   const real rho_n = g_ConVar_In[DENS][idx_in];
   real v1_B = g_ConVar_In[MOMX][idx_in] / rho_n;
   real v2_B = g_ConVar_In[MOMY][idx_in] / rho_n;
   real v3_B = g_ConVar_In[MOMZ][idx_in] / rho_n;
   RotateVec( sint, cost, sinp, cosp, v1_B, v2_B, v3_B );

   const real ec_source = v2_B * dpcdy_B + v3_B * dpcdz_B;
   if ( CR_Ec_source )   new_ec += dt_source * ec_source;  // TODO: Wrap this step in if ( CR_Ec_source )
#  endif

// 12. revert CR energy if negative
   if ( new_ec < TINY_NUMBER )
      new_ec = ec_old;

// 13. back-reaction to gas momentum and energy
   if ( CR_source ) {
//    momentum change = -(new_Fc - old_Fc) / vmax
      OneCell[MOMX] += -( newfr1 - fc1_old ) * invlim;
      OneCell[MOMY] += -( newfr2 - fc2_old ) * invlim;
      OneCell[MOMZ] += -( newfr3 - fc3_old ) * invlim;

//    energy change = -(new_Ec - old_Ec)
//    --> CR energy lost goes to gas thermal energy
      real new_eg = OneCell[ENGY] - ( new_ec - ec_old );
      if ( new_eg < (real)0.0 )
         new_eg = OneCell[ENGY];
      OneCell[ENGY] = new_eg;
   }

// 14. update CR fields
   OneCell[CR_E ] = new_ec;
   OneCell[CR_F1] = newfr1;
   OneCell[CR_F2] = newfr2;
   OneCell[CR_F3] = newfr3;

} // FUNCTION : CR_TwoMomentSource_HalfStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_FullStep
//
// Description : Apply implicit source term update for CR two-moment equations at full-step
//
// Note        : 1. Invoked by CPU_FluidSolver_MHM() after Hydro_FullStepUpdate()
//               2. Gas rho/momentum are taken from the post-flux g_Flu_Array_Out[]
//               3. Applies the gas back-reaction when CR_source is on
//               4. The ec_source term uses grad(Pc) from g_Flux, velocity from half-step g_PriVar_Half
//                  (= Athena++ stage-2 pre-stage w) just to match Athena++
//               5. Rotation uses the half-step B in g_PriVar_Half[MAG_OFFSET+*]
//                  (= Athena++ stage-2 b_angle)
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered half-step primitive variables
//               g_Output      : Array storing the post-flux fluid data to update
//               g_Flux        : Array storing the input face-centered fluxes
//               dt            : Full time step
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_Output[] (modified CR and gas fields)
//
// Reference   : Athena++ src/cr/integrators/cr_source.cpp, cr_transport.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real dt, const real dh, const MicroPhy_t *MicroPhy )
{
   const int  didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };

// offset between output array and primitive array
   const int  pvar_offset  = ( N_HF_VAR - PS2 ) / 2;

// CR parameters
   const real vmax   = MicroPhy->CR_vmax;
   const real invlim = (real)1.0 / vmax;

   const bool CR_source    = MicroPhy->CR_source;      // flag to enable back-reaction to gas
   const bool CR_stream    = MicroPhy->CR_stream;      // flag to enable streaming
   const bool CR_Ec_source = MicroPhy->CR_Ec_source;   // flag to include the CR energy source term

   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
//    output array indices
      const int i_out   = idx_out % PS2;
      const int j_out   = idx_out % SQR(PS2) / PS2;
      const int k_out   = idx_out / SQR(PS2);

//    primitive variable indices
      const int i_pvar  = i_out + pvar_offset;
      const int j_pvar  = j_out + pvar_offset;
      const int k_pvar  = k_out + pvar_offset;
      const int idx_pvar = IDX321( i_pvar, j_pvar, k_pvar, N_HF_VAR, N_HF_VAR );

//    1. get the CR state
//       --> CR_E is already floored by Hydro_FullStepUpdate(); FMAX is only a safeguard
      real ec  = FMAX( g_Output[CR_E ][idx_out], TINY_NUMBER );  // Should we use MIN_DENS instead?
      real fc1 = g_Output[CR_F1][idx_out];
      real fc2 = g_Output[CR_F2][idx_out];
      real fc3 = g_Output[CR_F3][idx_out];

//    2. get gas density and velocity
//       --> floor rho here since Hydro_FullStepUpdate() does not floor it (deferred to Flu_Close())
      const real rho = FMAX( g_Output[DENS][idx_out], TINY_NUMBER );
      real v1 = g_Output[MOMX][idx_out] / rho;
      real v2 = g_Output[MOMY][idx_out] / rho;
      real v3 = g_Output[MOMZ][idx_out] / rho;

//    3. get the half-step B field (see Note 5)
      const real Bx = g_PriVar_Half[MAG_OFFSET+MAGX][idx_pvar];
      const real By = g_PriVar_Half[MAG_OFFSET+MAGY][idx_pvar];
      const real Bz = g_PriVar_Half[MAG_OFFSET+MAGZ][idx_pvar];

//    4. compute B-field angles for rotation
      real sint, cost, sinp, cosp;
      CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

//    5. compute the total velocity
//       --> ADV_* in g_PriVar_Half[] are set by CR_UpdateStreaming()
      const real v_adv_x = g_PriVar_Half[ADV_VX][idx_pvar];
      const real v_adv_y = g_PriVar_Half[ADV_VY][idx_pvar];
      const real v_adv_z = g_PriVar_Half[ADV_VZ][idx_pvar];

      real vtot1 = v1;
      real vtot2 = v2;
      real vtot3 = v3;
      if ( CR_stream ) {
         vtot1 += v_adv_x;
         vtot2 += v_adv_y;
         vtot3 += v_adv_z;
      }

//    6. save original CR energy and flux for the floor check and gas back-reaction
      const real ec_old  = ec;
      const real fc1_old = fc1;
      const real fc2_old = fc2;
      const real fc3_old = fc3;

//    7. rotate all vectors to B-aligned frame
      real fr1 = fc1, fr2 = fc2, fr3 = fc3;

#     ifdef MHD
      RotateVec( sint, cost, sinp, cosp, v1, v2, v3 );
      RotateVec( sint, cost, sinp, cosp, fr1, fr2, fr3 );
      RotateVec( sint, cost, sinp, cosp, vtot1, vtot2, vtot3 );

//    perpendicular work is added via ec_source (step 11)
//    the vtot2/vtot3 terms below are kept for generality
      vtot2 = (real)0.0;
      vtot3 = (real)0.0;
#     endif

//    8. compute effective opacity
      const real sigma_diff = MicroPhy->CR_sigma;
      const real sigma_diff_perp = MicroPhy->CR_sigma_perp;
      const real sigma_adv_para = g_PriVar_Half[ADV_SIGMA][idx_pvar];
      const real sigma_adv_perp = MicroPhy->CR_max_opacity;

      real sigma_x = sigma_diff;
      real sigma_y = sigma_diff_perp;
      real sigma_z = sigma_diff_perp;
      if ( CR_stream ) {
         sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
         sigma_y = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
         sigma_z = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
      }

//    9. build implicit matrix and solve
//        Source terms:
//        dEc/dt = -vtot * sigma * (Fc - v*Ec*(4/3)/vmax)
//        dFc/dt = -vmax * sigma * (Fc - v*Ec*(4/3)/vmax)

      const real rhs1 = ec;
      const real rhs2 = fr1;
      const real rhs3 = fr2;
      const real rhs4 = fr3;

//    coefficients for Ec equation (row 1)
      const real coef_11 = (real)1.0 - dt * sigma_x * vtot1 * v1 * invlim * ((real)4.0/(real)3.0)
                                     - dt * sigma_y * vtot2 * v2 * invlim * ((real)4.0/(real)3.0)
                                     - dt * sigma_z * vtot3 * v3 * invlim * ((real)4.0/(real)3.0);
      const real coef_12 = dt * sigma_x * vtot1;
      const real coef_13 = dt * sigma_y * vtot2;
      const real coef_14 = dt * sigma_z * vtot3;

//    coefficients for Fr1 equation (row 2)
      const real coef_21 = -dt * v1 * sigma_x * ((real)4.0/(real)3.0);
      const real coef_22 = (real)1.0 + dt * vmax * sigma_x;

//    coefficients for Fr2 equation (row 3)
      const real coef_31 = -dt * v2 * sigma_y * ((real)4.0/(real)3.0);
      const real coef_33 = (real)1.0 + dt * vmax * sigma_y;

//    coefficients for Fr3 equation (row 4)
      const real coef_41 = -dt * v3 * sigma_z * ((real)4.0/(real)3.0);
      const real coef_44 = (real)1.0 + dt * vmax * sigma_z;

//    solve by substitution
      const real e_coef = coef_11 - coef_12 * coef_21 / coef_22
                                  - coef_13 * coef_31 / coef_33
                                  - coef_14 * coef_41 / coef_44;

      real new_ec = rhs1 - coef_12 * rhs2 / coef_22
                         - coef_13 * rhs3 / coef_33
                         - coef_14 * rhs4 / coef_44;
      new_ec /= e_coef;

//    CR_Ec_source=0 drops the Ec source term --> dEc/dt = -div(Fc) (Jiang & Oh 2018 Eq. 19)
      if ( !CR_Ec_source )   new_ec = rhs1;

//    back-substitute to get new flux
      real newfr1 = ( rhs2 - coef_21 * new_ec ) / coef_22;
      real newfr2 = ( rhs3 - coef_31 * new_ec ) / coef_33;
      real newfr3 = ( rhs4 - coef_41 * new_ec ) / coef_44;

//    10. rotate back to lab frame
#     ifdef MHD
      InvRotateVec( sint, cost, sinp, cosp, newfr1, newfr2, newfr3 );
#     endif

//    11. compute perpendicular work ec_source = v_perp*grad(Pc)
//    --> grad(Pc)[n] = div(flux of CR_F(n+1))/vmax, as in CR_UpdateStreaming()
#     ifdef MHD
      const real _dh = (real)1.0 / dh;
//    +1: skip the extra transverse fluxes used only by the CT electric field (as in Hydro_FullStepUpdate())
      const int  idx_flux = IDX321( i_out+1, j_out+1, k_out+1, N_FL_FLUX, N_FL_FLUX );
      const int  CRF_v[3] = { CR_F1, CR_F2, CR_F3 };
      real grad_pc[3];
      for (int n=0; n<3; n++) {
         grad_pc[n] = (real)0.0;
         for (int d=0; d<3; d++)
            grad_pc[n] += g_Flux[d][ CRF_v[n] ][idx_flux] - g_Flux[d][ CRF_v[n] ][ idx_flux - didx_flux[d] ];
         grad_pc[n] *= _dh * invlim;
      }

      real dpcdx_B = grad_pc[0], dpcdy_B = grad_pc[1], dpcdz_B = grad_pc[2];
      RotateVec( sint, cost, sinp, cosp, dpcdx_B, dpcdy_B, dpcdz_B );

//    compute perpendicular gas velocity (half-step primitive, see Note 4)
      real v1_B = g_PriVar_Half[1][idx_pvar];
      real v2_B = g_PriVar_Half[2][idx_pvar];
      real v3_B = g_PriVar_Half[3][idx_pvar];
      RotateVec( sint, cost, sinp, cosp, v1_B, v2_B, v3_B );

      const real ec_source = v2_B * dpcdy_B + v3_B * dpcdz_B;
      if ( CR_Ec_source )   new_ec += dt * ec_source;  // TODO: Wrap this step in if ( CR_Ec_source )
#     endif

//    12. revert CR energy if negative
      if ( new_ec < TINY_NUMBER )
         new_ec = ec_old;

//    13. back-reaction to gas momentum and energy
      if ( CR_source ) {
//       momentum change = -(new_Fc - old_Fc) / vmax
         g_Output[MOMX][idx_out] += -( newfr1 - fc1_old ) * invlim;
         g_Output[MOMY][idx_out] += -( newfr2 - fc2_old ) * invlim;
         g_Output[MOMZ][idx_out] += -( newfr3 - fc3_old ) * invlim;

//       energy change = -(new_Ec - old_Ec)
//       --> CR energy lost goes to gas thermal energy
         real new_eg = g_Output[ENGY][idx_out] - ( new_ec - ec_old );
         if ( new_eg < (real)0.0 )
            new_eg = g_Output[ENGY][idx_out];
         g_Output[ENGY][idx_out] = new_eg;
      }

//    14. update CR fields
      g_Output[CR_E ][idx_out] = new_ec;
      g_Output[CR_F1][idx_out] = newfr1;
      g_Output[CR_F2][idx_out] = newfr2;
      g_Output[CR_F3][idx_out] = newfr3;

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

} // FUNCTION : CR_TwoMomentSource_FullStep

#endif // #ifdef CR_TWOMOMENT
