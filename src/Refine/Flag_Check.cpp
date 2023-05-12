#include "GAMER.h"

static bool Check_Gradient( const int i, const int j, const int k, const real Input[], const double Threshold );
static bool Check_Curl( const int i, const int j, const int k,
                        const real vx[][PS1][PS1], const real vy[][PS1][PS1], const real vz[][PS1][PS1],
                        const double Threshold );
static bool ELBDM_Flag_VolumeFracQP( const real Cond[], const double Threshold_QP, const double Threshold_VolumeFraction );
extern bool (*Flag_Region_Ptr)( const int i, const int j, const int k, const int lv, const int PID );
extern bool (*Flag_User_Ptr)( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );


//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Check
// Description :  Check if the target cell (i,j,k) satisfies the refinement criteria
//
// Note        :  1. Useless input arrays are set to NULL (e.g, Pot[] if GRAVITY is off)
//                2. For OPT__FLAG_USER, the function pointer "Flag_User_Ptr" must be set by a
//                   test problem initializer
//
// Parameter   :  lv           : Target refinement level
//                PID          : Target patch ID
//                i,j,k        : Indices of the target cell
//                dv           : Cell volume at the target level
//                Fluid        : Input fluid array (with NCOMP_TOTAL components)
//                Pot          : Input potential array
//                MagCC        : Input cell-centered B field array
//                Vel          : Input velocity array
//                Pres         : Input pressure array
//                Lohner_Ave   : Input array storing the averages for the Lohner error estimator
//                Lohner_Slope : Input array storing the slopes for the Lohner error estimator
//                Lohner_NVar  : Number of variables stored in Lohner_Ave and Lohner_Slope
//                ParCount     : Input array storing the number of particles on each cell
//                               (note that it has the **real** type)
//                ParDens      : Input array storing the particle mass density on each cell
//                JeansCoeff   : Pi*GAMMA/(SafetyFactor^2*G), where SafetyFactor = FlagTable_Jeans[lv]
//                               --> Flag if dh^2 > JeansCoeff*Pres/Dens^2
//                Phase_Slope  : Input array storing the slope of the phase field to determine whether the dB wavelength is resolved
//                Interf_Cond  : Input array storing the dimensionless quantum pressures for the interference condition
//
// Return      :  "true"  if any  of the refinement criteria is satisfied
//                "false" if none of the refinement criteria is satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_Check( const int lv, const int PID, const int i, const int j, const int k, const real dv,
                 const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real MagCC[][PS1][PS1][PS1],
                 const real Vel[][PS1][PS1][PS1], const real Pres[][PS1][PS1],
                 const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                 const real ParCount[][PS1][PS1], const real ParDens[][PS1][PS1], const real JeansCoeff,
                 const real *Interf_Cond)
{

   bool Flag = false;

// check whether the input cell is within the regions allowed to be refined
// ===========================================================================================
   if ( OPT__FLAG_REGION )
   {
      if ( Flag_Region_Ptr == NULL )   Aux_Error( ERROR_INFO, "Flag_Region_Ptr == NULL for OPT__FLAG_REGION !!\n" );

      if (  !Flag_Region_Ptr( i, j, k, lv, PID )  )    return false;
   }


// check ELBDM interference
// if FlagTable_Interference[lv][2] > 0 we set use wave flag if interference > FlagTable_Interference[lv][0],[1]
// if FlagTable_Interference[lv][2] < 0 we just use it as refinement criterion for fluid patches
// ===========================================================================================
#  if ( MODEL == ELBDM )
   if ( OPT__FLAG_INTERFERENCE )
   {
      real (*Var)  [PS1 ][PS1 ][PS1 ] = ( real(*) [PS1 ][PS1 ][PS1 ] )  Interf_Cond;
      bool FlagIntQP = false, FlagIntPhaseDiscont = false, dBWavelengthNotResolved = false;

      if ( FlagTable_Interference[lv][1] > 0 )
         FlagIntQP  =  ELBDM_Flag_VolumeFracQP( Interf_Cond, FlagTable_Interference[lv][0], FlagTable_Interference[lv][1]);
      else
         FlagIntQP  =  ELBDM_Flag_Interference( i, j, k, Interf_Cond,                       FlagTable_Interference[lv][0]);

      Flag |= FlagIntQP;


//    For wave solver: Check whether dB wavelength is resolved using FlagTable_Interference[lv][2] as maximum phase difference between neighbouring cells

#     if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )
      if ( amr->use_wave_flag[lv] ) {
#     endif // # if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )

      dBWavelengthNotResolved = ELBDM_Flag_Interference( i, j, k, Interf_Cond + 2 * CUBE(PS1), FlagTable_Interference[lv][2] );

      Flag |= dBWavelengthNotResolved;

#     if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )
//    For fluid solver: Check phase curvature using FlagTable_Interference[lv][2] as maximum allowed phase curvature
      } else {
         FlagIntPhaseDiscont  =  ELBDM_Flag_Interference( i, j, k, Interf_Cond + CUBE(PS1), FlagTable_Interference[lv][2]);

         Flag |= FlagIntPhaseDiscont;
      }
#     endif // # if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )

#     if ( ELBDM_SCHEME == HYBRID )
      if ( Flag &&  FlagTable_Interference[lv][3] >= 0.0 ) {

#        ifdef GAMER_DEBUG
//       we distinguish two cases:
//       if the slope of the phase field is lower than PI, we know that the dB wavelength will be resolved after refinement and we can safely switch to the wave scheme
//       if the slope is bigger, this may be due to a real phase discontinuity because of a 2 pi winding of the phase at a point of zero density
//       we check for the latter with FlagInterferenceTwo
//       if the curvature of the phase is high at a point where the phase jump is bigger than pi, we assume that it is a real jump and still switch to the wave scheme
         bool dBResolvedAfterRefine = ! ELBDM_Flag_Interference( i, j, k, Interf_Cond + 2 * CUBE(PS1), M_PI ) || FlagIntQP;
         if ( !dBResolvedAfterRefine ) {
            const int    Idx               = 2 * PS1*PS1*PS1 + k*PS1*PS1 + j*PS1 + i;
            const real   PhaseDifference   = Interf_Cond[Idx];
//          convert coordinates of cell to global integer coordinate system
            int coordinates[3] = {i, j, k};

            for ( int l = 0; l < 3; ++l ) {
               coordinates[l] *= amr->scale[lv];
               coordinates[l] += amr->patch[0][lv][PID]->corner[l];
            }
            Aux_Message( stdout, "WARNING: Switching to wave solver, but dB wavelength not resolved at lv %d i %d j %d k %d for phase jump = %4.2f on rank %d\n", lv, coordinates[0], coordinates[1], coordinates[2], PhaseDifference, MPI_Rank);

         }
#        endif

         amr->patch[0][lv][PID]->use_wave_flag =  true;

      }
#     endif // # if ( ELBDM_SCHEME == HYBRID )

      if ( Flag )
            return Flag;
   }
#  endif // # if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )


#  ifdef PARTICLE
// check the number of particles on each cell
// ===========================================================================================
   if ( OPT__FLAG_NPAR_CELL )
   {
      Flag |= ( ParCount[k][j][i] > FlagTable_NParCell[lv] );
      if ( Flag )    return Flag;
   }


// check the particle mass on each cell
// ===========================================================================================
   if ( OPT__FLAG_PAR_MASS_CELL )
   {
      Flag |= ( ParDens[k][j][i]*dv > FlagTable_ParMassCell[lv] );
      if ( Flag )    return Flag;
   }
#  endif


#  ifdef DENS
// check density magnitude
// ===========================================================================================
if ( OPT__FLAG_RHO )
   {
      Flag |= ( Fluid[DENS][k][j][i] > FlagTable_Rho[lv] );
      if ( Flag )    return Flag;
   }


// check density gradient
// ===========================================================================================
   if ( OPT__FLAG_RHO_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Fluid[DENS][0][0][0], FlagTable_RhoGradient[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check pressure gradient
// ===========================================================================================
#  if ( MODEL == HYDRO )
   if ( OPT__FLAG_PRES_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Pres[0][0][0], FlagTable_PresGradient[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check vorticity
// ===========================================================================================
#  if ( MODEL == HYDRO )
   if ( OPT__FLAG_VORTICITY )
   {
      Flag |= Check_Curl( i, j, k, Vel[0], Vel[1], Vel[2], FlagTable_Vorticity[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check current density in MHD
// ===========================================================================================
#  ifdef MHD
   if ( OPT__FLAG_CURRENT )
   {
      Flag |= Check_Curl( i, j, k, MagCC[0], MagCC[1], MagCC[2], FlagTable_Current[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check Jeans length
// ===========================================================================================
#  if ( MODEL == HYDRO  &&  defined GRAVITY )
   if ( OPT__FLAG_JEANS )
   {
#     ifdef GAMER_DEBUG
      if ( Pres == NULL )  Aux_Error( ERROR_INFO, "Pres == NULL !!\n" );
#     endif

      const real Dens_1Cell = Fluid[DENS][k][j][i];
      const real Pres_1Cell = Pres       [k][j][i];

      Flag |= (  SQR(amr->dh[lv]) > JeansCoeff*Pres_1Cell/SQR( Dens_1Cell )  );
      if ( Flag )    return Flag;
   }
#  endif


// check ELBDM energy density
// ===========================================================================================
#  if ( MODEL == ELBDM )
#  if ( ELBDM_SCHEME == HYBRID )
   if ( amr->use_wave_flag[lv] )
#  endif
   if ( OPT__FLAG_ENGY_DENSITY )
   {
      Flag |= ELBDM_Flag_EngyDensity( i, j, k, &Fluid[REAL][0][0][0], &Fluid[IMAG][0][0][0],
                                      FlagTable_EngyDensity[lv][0], FlagTable_EngyDensity[lv][1] );
      if ( Flag )    return Flag;
   }
#  endif



// check Lohner's error estimator
// ===========================================================================================

#  if ( ELBDM_SCHEME == HYBRID )
   if ( amr->use_wave_flag[lv] )
#  endif
   if ( Lohner_NVar > 0 )
   {
//    check Lohner only if density is greater than the minimum threshold
#     ifdef DENS
      if ( Fluid[DENS][k][j][i] >= FlagTable_Lohner[lv][4] )
#     endif
      Flag |= Flag_Lohner( i, j, k, OPT__FLAG_LOHNER_FORM, Lohner_Var, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                           FlagTable_Lohner[lv][0], FlagTable_Lohner[lv][2], FlagTable_Lohner[lv][3] );
      if ( Flag )    return Flag;
   }


// check user-defined criteria
// ===========================================================================================
   if ( OPT__FLAG_USER )
   {
      if ( Flag_User_Ptr != NULL )
      {
         Flag |= Flag_User_Ptr( i, j, k, lv, PID, FlagTable_User[lv] );
         if ( Flag )    return Flag;
      }

      else
         Aux_Error( ERROR_INFO, "Flag_User_Ptr == NULL for OPT__FLAG_USER !!\n" );
   }


   return Flag;

} // FUNCTION : Flag_Check



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Gradient
// Description :  Check if the gradient of the input data at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Size of the array "Input" should be PATCH_SIZE^3
//                2. For cells adjacent to the patch boundaries, only first-order approximation is adopted
//                   to estimate gradient. Otherwise, second-order approximation is adopted.
//                   --> Advantage: NO need to prepare the ghost-zone data for the target patch
//
// Parameter   :  i,j,k     : Indices of the target cell in the array "Input"
//                Input     : Input array
//                Threshold : Threshold for the flag operation
//
// Return      :  "true"  if the gradient is larger           than the given threshold
//                "false" if the gradient is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool Check_Gradient( const int i, const int j, const int k, const real Input[], const double Threshold )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0 ||  j >= PS1  ||  k < 0  ||  k >= PS1   )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif


   const int ijk[3]  = { i, j, k };
   const int Idx     = k*PS1*PS1 + j*PS1 + i;
   const int dIdx[3] = { 1, PS1, PS1*PS1 };

   int  Idx_p, Idx_m;
   real _dh, Self, Gradient;
   bool Flag = false;

   for (int d=0; d<3; d++)
   {
      switch ( ijk[d] )
      {
         case 0     : Idx_m = Idx;           Idx_p = Idx+dIdx[d];    _dh = (real)1.0;  break;
         case PS1-1 : Idx_m = Idx-dIdx[d];   Idx_p = Idx;            _dh = (real)1.0;  break;
         default    : Idx_m = Idx-dIdx[d];   Idx_p = Idx+dIdx[d];    _dh = (real)0.5;  break;
      }

      Self     = Input[Idx];
      Gradient = _dh*( Input[Idx_p] - Input[Idx_m] );
      Flag    |= (  FABS( Gradient/Self ) > Threshold  );

      if ( Flag )    return Flag;
   } // for (int d=0; d<3; d++)

   return Flag;

} // FUNCTION : Check_Gradient



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Curl
// Description :  Check if the curl of the input vector at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Flag if |curl(v)|*dh/|v| > threshold
//                2. For cells adjacent to the patch boundaries, only first-order approximation is adopted
//                   to estimate derivatives. Otherwise, second-order approximation is adopted.
//                   --> Advantage: NO need to prepare the ghost-zone data for the target patch
//                3. Size of the input arrays "vx/y/z" should be PATCH_SIZE^3
//                   --> They should store **cell-centered** values
//
// Parameter   :  i,j,k     : Target array indices
//                vx/y/z    : Input vectors
//                Threshold : Refinement threshold
//
// Return      :  "true"  if |curl(v)|*dh/|v| >  threshold
//                "false" if |curl(v)|*dh/|v| <= threshold
//-------------------------------------------------------------------------------------------------------
bool Check_Curl( const int i, const int j, const int k,
                 const real vx[][PS1][PS1], const real vy[][PS1][PS1], const real vz[][PS1][PS1],
                 const double Threshold )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0 ||  j >= PS1  ||  k < 0  ||  k >= PS1   )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   int  im, ip, jm, jp, km, kp;
   real _dx, _dy, _dz, v2, wx, wy, wz, w2;
   bool Flag = false;


// check if the target cell indices are adjacent to the patch boundaries
   if      ( i == 0     ) { im = 0;       ip = 1;       _dx = (real)1.0; }
   else if ( i == PS1-1 ) { im = PS1-2;   ip = PS1-1;   _dx = (real)1.0; }
   else                   { im = i-1;     ip = i+1;     _dx = (real)0.5; }

   if      ( j == 0     ) { jm = 0;       jp = 1;       _dy = (real)1.0; }
   else if ( j == PS1-1 ) { jm = PS1-2;   jp = PS1-1;   _dy = (real)1.0; }
   else                   { jm = j-1;     jp = j+1;     _dy = (real)0.5; }

   if      ( k == 0     ) { km = 0;       kp = 1;       _dz = (real)1.0; }
   else if ( k == PS1-1 ) { km = PS1-2;   kp = PS1-1;   _dz = (real)1.0; }
   else                   { km = k-1;     kp = k+1;     _dz = (real)0.5; }


// calculate magnitude
   v2 = SQR( vx[k][j][i] ) + SQR( vy[k][j][i] ) + SQR( vz[k][j][i] );


// calculate w=curl(v)*dh
   wx = _dy*( vz[k ][jp][i ] - vz[k ][jm][i ] ) - _dz*( vy[kp][j ][i ] - vy[km][j ][i ] );
   wy = _dz*( vx[kp][j ][i ] - vx[km][j ][i ] ) - _dx*( vz[k ][j ][ip] - vz[k ][j ][im] );
   wz = _dx*( vy[k ][j ][ip] - vy[k ][j ][im] ) - _dy*( vx[k ][jp][i ] - vx[k ][jm][i ] );
   w2 = SQR(wx) + SQR(wy) + SQR(wz);


// flag if |curl(v)|*dh/|v| > threshold
   Flag = ( w2/v2 > SQR(Threshold) );

   return Flag;

} // FUNCTION : Check_Curl




//-------------------------------------------------------------------------------------------------------
// Function    :  Check_VolumeFracQP
// Description :  Check if the quantum pressure of more than Threshold_VolumeFraction of the cells in the input patch exceed Threshold_QP
//
// Note        :  1, Size of the input array "Cond" should be PATCH_SIZE^3
//
// Parameter   :  Cond      : Input vector
//                Threshold_QP : QP threshold
//                Threshold_VolumeFraction : Volume Fraction Threshold ( 0.0 - 1.0)
//
// Return      :  "true"  if sum | cell > Threshold_QP | / NCell > Theshold_VolumeFraction
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool ELBDM_Flag_VolumeFracQP( const real Cond[], const double Threshold_QP, const double Threshold_VolumeFraction )
{
   const float NTotal = PS1 * PS1 * PS1;
   float NExceedThreshold = 0;

   for (int Idx = 0; Idx < NTotal; ++Idx)
         if ( Cond[Idx] > Threshold_QP ) ++NExceedThreshold;

   float ratio =  NExceedThreshold / NTotal;

   return ratio > Threshold_VolumeFraction;

} // FUNCTION : Check_VolumeFracQP
