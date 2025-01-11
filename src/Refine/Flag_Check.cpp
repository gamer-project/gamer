#include "GAMER.h"

static bool Check_Gradient( const int i, const int j, const int k, const real Input[], const double Threshold );
static bool Check_Curl( const int i, const int j, const int k,
                        const real vx[][PS1][PS1], const real vy[][PS1][PS1], const real vz[][PS1][PS1],
                        const double Threshold );
static bool Check_Angular_Max( const int i, const int j, const int k, const int lv, const int PID,
                               const double CenX, const double CenY, const double CenZ,
                               const double AngRes_Max, const double AngRes_Max_R );
static bool Check_Angular_Min( const int i, const int j, const int k, const int lv, const int PID,
                               const double CenX, const double CenY, const double CenZ,
                               const double AngRes_Min );
static bool Check_Radial( const int i, const int j, const int k, const int lv, const int PID,
                          const double CenX, const double CenY, const double CenZ, const double Refine_Rad );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Check
// Description :  Check if the target cell (i,j,k) satisfies the refinement criteria
//
// Note        :  1. Useless input arrays are set to NULL (e.g, Pot[] if GRAVITY is off)
//                2. For OPT__FLAG_USER, the function pointer "Flag_User_Ptr" must be set by a
//                   test problem initializer
//
// Parameter   :  lv            : Target refinement level
//                PID           : Target patch ID
//                i,j,k         : Indices of the target cell
//                dv            : Cell volume at the target level
//                Fluid         : Input fluid array (with NCOMP_TOTAL components)
//                Pot           : Input potential array
//                MagCC         : Input cell-centered B field array
//                Vel           : Input velocity array
//                Pres          : Input pressure array
//                Lrtz          : Input Lorentz factor array
//                Lohner_Ave    : Input array storing the averages for the Lohner error estimator
//                Lohner_Slope  : Input array storing the slopes for the Lohner error estimator
//                Lohner_NVar   : Number of variables stored in Lohner_Ave and Lohner_Slope
//                ParCount      : Input array storing the number of particles on each cell
//                                (note that it has the **real** type)
//                ParDens       : Input array storing the particle mass density on each cell
//                JeansCoeff    : Pi*GAMMA/(SafetyFactor^2*G), where SafetyFactor = FlagTable_Jeans[lv]
//                                --> Flag if dh^2 > JeansCoeff*Pres/Dens^2
//                Interf_Var    : Input array storing the density and phase for the interference condition
//                Spectral_Cond : Input variable storing the spectral refinement condition
//
// Return      :  "true"  if any  of the refinement criteria is satisfied
//                "false" if none of the refinement criteria is satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_Check( const int lv, const int PID, const int i, const int j, const int k, const real dv,
                 const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real MagCC[][PS1][PS1][PS1],
                 const real Vel[][PS1][PS1][PS1], const real Pres[][PS1][PS1], const real Lrtz[][PS1][PS1],
                 const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                 const real ParCount[][PS1][PS1], const real ParDens[][PS1][PS1], const real JeansCoeff,
                 const real *Interf_Var, const real Spectral_Cond )
{

   bool Flag = false;


// *******************************************************************************************
// refinement flags must be checked in the following order
// 1. no-refinement criteria --> exclude patches outside the regions allowed for refinement
// 2. OPT__FLAG_INTERFERENCE --> ensure amr->patch[0][lv][PID]->switch_to_wave_flag is set correctly
// 3. refinement criteria
// *******************************************************************************************


// *****************************
// 1. no-refinement criteria
// *****************************

// check whether the input cell is within the regions allowed to be refined
// ===========================================================================================
   if ( OPT__FLAG_REGION )
   {
      if ( Flag_Region_Ptr == NULL )   Aux_Error( ERROR_INFO, "Flag_Region_Ptr == NULL for OPT__FLAG_REGION !!\n" );

      if (  !Flag_Region_Ptr( i, j, k, lv, PID )  )    return false;
   }


// check maximum angular resolution
// ===========================================================================================
   if ( OPT__FLAG_ANGULAR )
   {
      bool Within = Check_Angular_Max( i, j, k, lv, PID, FLAG_ANGULAR_CEN_X, FLAG_ANGULAR_CEN_Y,
                                       FLAG_ANGULAR_CEN_Z, FlagTable_Angular[lv][0],
                                       FlagTable_Angular[lv][2] );
      if ( ! Within )   return false;
   }



// *****************************
// 2. OPT__FLAG_INTERFERENCE
// *****************************

// ELBDM interference check must be performed before any other refinement checks in order to set switch_to_wave_flag correctly
// ===========================================================================================
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( OPT__FLAG_INTERFERENCE  &&  !amr->use_wave_flag[lv] )
   {
      Flag |= ELBDM_Flag_Interference( i, j, k, Interf_Var, FlagTable_Interference[lv][0], FlagTable_Interference[lv][1],
                                       FlagTable_Interference[lv][2], FlagTable_Interference[lv][3]>0.5 );

//    switch to wave solver when refining to ELBDM_FIRST_WAVE_LEVEL
      if ( Flag  &&  lv+1 >= ELBDM_FIRST_WAVE_LEVEL )    amr->patch[0][lv][PID]->switch_to_wave_flag = true;

      if ( Flag )    return Flag;
   }
#  endif



// *****************************
// 3. refinement criteria
// *****************************

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


// check Lorentz factor gradient in SRHD
// ===========================================================================================
#  if ( MODEL == HYDRO  &&  defined SRHD )
   if ( OPT__FLAG_LRTZ_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Lrtz[0][0][0], FlagTable_LrtzGradient[lv] );
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


// check cosmic ray
// ===========================================================================================
#  ifdef COSMIC_RAY
   if ( OPT__FLAG_CRAY )
   {
      Flag |= ( Fluid[CRAY][k][j][i] > FlagTable_CRay[lv] );
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
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   if ( OPT__FLAG_ENGY_DENSITY )
   {
      Flag |= ELBDM_Flag_EngyDensity( i, j, k, &Fluid[REAL][0][0][0], &Fluid[IMAG][0][0][0],
                                      FlagTable_EngyDensity[lv][0], FlagTable_EngyDensity[lv][1] );
      if ( Flag )    return Flag;
   }
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } // if ( amr->use_wave_flag[lv] )
#  endif
#  endif // ELBDM


// check ELBDM spectral criterion
// ===========================================================================================
#  if ( MODEL == ELBDM )
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   if ( OPT__FLAG_SPECTRAL )
   {
      Flag |= ( Spectral_Cond > FlagTable_Spectral[lv][0] );
      if ( Flag )    return Flag;
   }
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } // if ( amr->use_wave_flag[lv] )
#  endif
#  endif // ELBDM


// check Lohner's error estimator
// ===========================================================================================
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
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
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } // if ( amr->use_wave_flag[lv] )
#  endif


// check minimum angular resolution
// ===========================================================================================
   if ( OPT__FLAG_ANGULAR )
   {
      Flag |= Check_Angular_Min( i, j, k, lv, PID, FLAG_ANGULAR_CEN_X, FLAG_ANGULAR_CEN_Y,
                                 FLAG_ANGULAR_CEN_Z, FlagTable_Angular[lv][1] );
      if ( Flag )    return Flag;
   }


// check radial resolution
// ===========================================================================================
   if ( OPT__FLAG_RADIAL )
   {
      Flag |= Check_Radial( i, j, k, lv, PID, FLAG_RADIAL_CEN_X, FLAG_RADIAL_CEN_Y,
                            FLAG_RADIAL_CEN_Z, FlagTable_Radial[lv] );
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
// Function    :  Check_Angular_Max
// Description :  Check if dh/R at cell (i,j,k) is smaller than the maximum angular resolution
//
// Note        :  1. Enabled by the runtime option "OPT__FLAG_ANGULAR"
//
// Parameter   :  i,j,k        : Target cell indices in the patch amr->patch[0][lv][PID]
//                lv           : Refinement level of the target patch
//                PID          : ID of the target patch
//                CenX/Y/Z     : x/y/z-coordinate of the center for calculating angular resolution
//                AngRes_Max   : Maximum allowed angular resolution (in radians)
//                AngRes_Max_R : Minimum radius to apply AngRes_Max
//
// Return      :  "true/false"  if the input cell "is/is not" within the region allowed for refinement
//-------------------------------------------------------------------------------------------------------
bool Check_Angular_Max( const int i, const int j, const int k, const int lv, const int PID,
                        const double CenX, const double CenY, const double CenZ,
                        const double AngRes_Max, const double AngRes_Max_R )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   const double dh     = amr->dh[lv];                                         // cell size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,     // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double dR [3] = { Pos[0]-CenX, Pos[1]-CenY, Pos[2]-CenZ };
   const double R      = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   return ( AngRes_Max < 0.0  ||  2.0 * R * AngRes_Max <= dh  ||  R <= AngRes_Max_R );

} // FUNCTION : Check_Angular_Max



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Angular_Min
// Description :  Check if dh/R at cell (i,j,k) is larger than the minimum angular resolution
//
// Note        :  1. Enabled by the runtime option "OPT__FLAG_ANGULAR"
//
// Parameter   :  i,j,k        : Target cell indices in the patch amr->patch[0][lv][PID]
//                lv           : Refinement level of the target patch
//                PID          : ID of the target patch
//                CenX/Y/Z     : x/y/z-coordinate of the center for calculating angular resolution
//                AngRes_Min   : Minimum allowed angular resolution (in radians)
//
// Return      :  "true"  if the minimum angular resolution is not reached
//                "false" if the minimum angular resolution is     reached
//-------------------------------------------------------------------------------------------------------
bool Check_Angular_Min( const int i, const int j, const int k, const int lv, const int PID,
                        const double CenX, const double CenY, const double CenZ,
                        const double AngRes_Min )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   const double dh     = amr->dh[lv];                                         // cell size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,     // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double dR [3] = { Pos[0]-CenX, Pos[1]-CenY, Pos[2]-CenZ };
   const double R      = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   return ( AngRes_Min >= 0.0  &&  R * AngRes_Min < dh );

} // FUNCTION : Check_Angular_Min



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Radial
// Description :  Check if the cell is within the given radius
//
// Note        :  1. Enabled by the runtime option "OPT__FLAG_RADIAL"
//
// Parameter   :  i,j,k       : Target cell indices in the patch amr->patch[0][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                CenX/Y/Z    : x/y/z-coordinate of the center for calculating radial resolution
//                Refine_Rad  : Radius at level lv within which grids are refined
//
// Return      :  "true"  if r <  Refine_Rad
//                "false" if r >= Refine_Rad or Refine_Rad is not set
//-------------------------------------------------------------------------------------------------------
bool Check_Radial( const int i, const int j, const int k, const int lv, const int PID,
                   const double CenX, const double CenY, const double CenZ, const double Refine_Rad )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   const double dh     = amr->dh[lv];                                         // cell size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,     // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double dR [3] = { Pos[0]-CenX, Pos[1]-CenY, Pos[2]-CenZ };
   const double R      = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

// do not refine if the target radius is not set
   if ( Refine_Rad < 0.0 )   return false;

// refine the region within r < Refine_Rad and the innermost cells
   if ( R < dh  ||  R < Refine_Rad )   return true;

   return false;

} // FUNCTION : Check_Radial
