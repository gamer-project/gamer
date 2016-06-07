#include "Copyright.h"
#include "GAMER.h"

static bool Check_Gradient( const int i, const int j, const int k, const real Input[], const double Threshold );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Check
// Description :  Check if the targeted cell (i,j,k) satisfies the refinement criteria
//
// Note        :  Useless input arrays are set to NULL 
//                (e.g, Pot if GRAVITY is off, Pres if OPT__FLAG_PRES_GRADIENT is off)
//
// Parameter   :  lv             : Targeted refinement level 
//                PID            : Targeted patch ID
//                i,j,k          : Indices of the targeted cell
//                Fluid          : Input fluid array (with NCOMP components)
//                Pot            : Input potential array
//                Pres           : Input pressure array
//                Lohner_Ave     : Input array storing the averages for the Lohner error estimator
//                Lohner_Slope   : Input array storing the slopes for the Lohner error estimator
//                Lohner_NVar    : Number of variables stored in Lohner_Ave and Lohner_Slope
//                ParCount       : Input array storing the number of particles on each cell
//                                 (note that it has the **real** type)
//
// Return      :  "true"  if any  of the refinement criteria is satisfied
//                "false" if none of the refinement criteria is satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_Check( const int lv, const int PID, const int i, const int j, const int k, 
                 const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real Pres[][PS1][PS1],
                 const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                 const real ParCount[][PS1][PS1] )
{
   
   bool Flag = false;

// check whether the input cell is within the regions allowed to be refined
// ===========================================================================================
   if ( OPT__FLAG_REGION )
   {
      if (  !Flag_Region( i, j, k, lv, PID )  )    return false;
   }


#  ifdef PARTICLE
// check the number of particles on each cell
// ===========================================================================================
   if ( OPT__FLAG_NPAR_CELL )
   {
      Flag |= ( ParCount[k][j][i] > FlagTable_NParCell[lv] );
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
#  if   ( MODEL == HYDRO )
   if ( OPT__FLAG_PRES_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Pres[0][0][0], FlagTable_PresGradient[lv] );
      if ( Flag )    return Flag;
   }
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL


// check ELBDM energy density
// ===========================================================================================
#  if ( MODEL == ELBDM )
   if ( OPT__FLAG_ENGY_DENSITY )
   {
      Flag |= ELBDM_Flag_EngyDensity( i, j, k, &Fluid[REAL][0][0][0], &Fluid[IMAG][0][0][0], 
                                      FlagTable_EngyDensity[lv][0], FlagTable_EngyDensity[lv][1] );
      if ( Flag )    return Flag;
   }
#  endif


// check Lohner's error estimator
// ===========================================================================================
   if ( Lohner_NVar > 0 )
   {
      Flag |= Flag_Lohner( i, j, k, OPT__FLAG_LOHNER_FORM, Lohner_Var, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                           FlagTable_Lohner[lv][0], FlagTable_Lohner[lv][1], FlagTable_Lohner[lv][2] );
      if ( Flag )    return Flag;
   }


// check user-defined criteria
// ===========================================================================================
   if ( OPT__FLAG_USER )
   {
      Flag |= Flag_UserCriteria( i, j, k, lv, PID, FlagTable_User[lv] );
      if ( Flag )    return Flag;
   }


   return Flag;

} // FUNCTION : Flag_Check



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Gradient
// Description :  Check if the gradient of the input data at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Size of the array "Input" should be PATCH_SIZE^3 
//                2. For cells adjacent to the patch boundary, only first-order approximation is adopted 
//                   to estimate gradient. Otherwise, second-order approximation is adopted.
//                   --> Do NOT need to prepare the ghost-zone data for the targeted patch
//
// Parameter   :  i,j,k       : Indices of the targeted cell in the array "Input"
//                Input       : Input array
//                Threshold   : Threshold for the flag operation
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
      Flag    |= (  fabs( Gradient/Self ) > Threshold  );

      if ( Flag )    return Flag;
   } // for (int d=0; d<3; d++)

   return Flag;

} // FUNCTION : Check_Gradient
