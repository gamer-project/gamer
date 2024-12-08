#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Lohner
// Description :  Check if the numerical error at the cell (i,j,k) estimated by Lohner's prescription
//                exceeds the given threshold
//
// Note        :  1. Invoked by the function "Flag_Check"
//                2. Adopt the modified version in FLASH4 and MPI_AMRVAC
//                3. The Ave1D and Slope1D arrays must be prepared in advance by invoking the function "Prepare_for_Lohner"
//                4. The sizes of the arrays Ave1D and  Slope1D must be NVar*3*(PS1+2)^3
//
// Parameter   :  i,j,k       : Target cell indices in the target patch
//                Form        : Form of the Lohner's error estimator (FLASH4/form-invariant + Stencil = 1/2)
//                Var1D       : Array storing the input variables for the Lohner error estimator (used only for stencil == 1)
//                Ave1D       : Array storing the input averages for the Lohner error estimator
//                Slope1D     : Array storing the input slopes for the Lohner error estimator
//                NVar        : Number of variables stored in Ave and Slope arrays (must be 2 in ELBDM)
//                Threshold   : Refinement threshold
//                Filter      : Filter parameter for preventing refinement of small ripples
//                Soften      : Minimum number in the denominator --> error = sqrt( N/max(D,Soften) ), where
//                              N and D are numerator and denominator in the Lohner's formula, respectively
//
// Return      :  "true"  if the estimated error is larger           than the given threshold
//                "false" if the estimated error is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool Flag_Lohner( const int i, const int j, const int k, const OptLohnerForm_t Form, const real *Var1D, const real *Ave1D,
                  const real *Slope1D, const int NVar, const double Threshold, const double Filter, const double Soften )
{

// return if there is no target variable (e.g., for fluid levels in ELBDM hybrid scheme)
   if ( NVar == 0 )  return false;


// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0  ||  j >= PS1  ||  k < 0  ||  k >= PS1  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );

#  if ( ELBDM_SCHEME == ELBDM_WAVE )
   if ( NVar != 2 )  Aux_Error( ERROR_INFO, "NVar (%d) != 2 when ELBDM_SCHEME == ELBDM_WAVE !!\n", NVar );
#  endif
#  endif // #ifdef GAMER_DEBUG


   const int NCell   = PS1 + 4;                                                     // size of the array Var
   const int NAve    = PS1 + 2;                                                     // size of the array Ave
   const int NSlope  = NAve;                                                        // size of the array Slope
   const int Stencil = ( Form==LOHNER_FLASH2 || Form==LOHNER_FORM_INV2 ) ? 2 : 1;   // stencil size

// ii,jj,kk: for the arrays "Ave and Slope"
   const int ii  =  i + 1;   const int jj  =  j + 1;   const int kk  =  k + 1;
   const int iip = ii + 1;   const int jjp = jj + 1;   const int kkp = kk + 1;
   const int iim = ii - 1;   const int jjm = jj - 1;   const int kkm = kk - 1;

// i2,j2,k2: for the array "Var"
   const int i2  =  i + 2;   const int j2  =  j + 2;   const int k2  =  k + 2;   // for stencil = 1 only
   const int i2p = i2 + 1;   const int j2p = j2 + 1;   const int k2p = k2 + 1;   // for stencil = 1 only
   const int i2m = i2 - 1;   const int j2m = j2 - 1;   const int k2m = k2 - 1;   // for stencil = 1 only

   const real Filter_4  = (real)0.2500*Filter;                                   // for stencil = 1 only
   const real Soften_16 = (real)0.0625*Soften;                                   // for stencil = 1 only


   real Der2_xx, Der2_yy, Der2_zz;           // tensor Der2_ij = d2(Var)/didj (i/j=x,y,z)
   real Der2_xy, Der2_yz, Der2_zx;           // Der2_xy = Der2_yx, ...

// off-diagonal terms are useful for LOHNER_FLASH1/2 only
   real Der1_xx, Der1_yy, Der1_zz;
   real Der1_xy, Der1_xz;                    // tensor Der1_ij = d(Var)/di averaged over j (i/j=x,y,z)
   real Der1_yx, Der1_yz;                    // Der1_yx != Der1_xy
   real Der1_zx, Der1_zy;

   real Filter_xx, Filter_yy, Filter_zz;     // filter tensor along different directions
   real Filter_xy, Filter_yz, Filter_zx;     // Filter_xy = Filter_yx, ...

   real Nume=0.0, Deno=0.0, Error;
   bool Flag;


// convert the 1D arrays
   real (*Var)     [NCell ][NCell ][NCell ] = ( real(*)   [NCell ][NCell ][NCell ] )  Var1D; // for stencil = 1 only
   real (*Ave  )[3][NAve  ][NAve  ][NAve  ] = ( real(*)[3][NAve  ][NAve  ][NAve  ] )  Ave1D;
   real (*Slope)[3][NSlope][NSlope][NSlope] = ( real(*)[3][NSlope][NSlope][NSlope] )Slope1D;


   for (int v=0; v<NVar; v++)
   {
//    1. numerator
      if ( Stencil == 2 )
      {
         Der2_xx = Slope[v][0][kk ][jj ][iip] - Slope[v][0][kk ][jj ][iim];
         Der2_yy = Slope[v][1][kk ][jjp][ii ] - Slope[v][1][kk ][jjm][ii ];
         Der2_zz = Slope[v][2][kkp][jj ][ii ] - Slope[v][2][kkm][jj ][ii ];
         Der2_xy = Slope[v][1][kk ][jj ][iip] - Slope[v][1][kk ][jj ][iim];
         Der2_yz = Slope[v][2][kk ][jjp][ii ] - Slope[v][2][kk ][jjm][ii ];
         Der2_zx = Slope[v][0][kkp][jj ][ii ] - Slope[v][0][kkm][jj ][ii ];
      }

      else // Stencil == 1
      {
         Der2_xx = Var[v][k2 ][j2 ][i2p] - (real)2.0*Var[v][k2][j2][i2] + Var[v][k2 ][j2 ][i2m];
         Der2_yy = Var[v][k2 ][j2p][i2 ] - (real)2.0*Var[v][k2][j2][i2] + Var[v][k2 ][j2m][i2 ];
         Der2_zz = Var[v][k2p][j2 ][i2 ] - (real)2.0*Var[v][k2][j2][i2] + Var[v][k2m][j2 ][i2 ];

         Der2_xy = (real)0.25*( + Var[v][k2 ][j2p][i2p] + Var[v][k2 ][j2m][i2m]
                                - Var[v][k2 ][j2p][i2m] - Var[v][k2 ][j2m][i2p] );
         Der2_yz = (real)0.25*( + Var[v][k2p][j2p][i2 ] + Var[v][k2m][j2m][i2 ]
                                - Var[v][k2m][j2p][i2 ] - Var[v][k2p][j2m][i2 ] );
         Der2_zx = (real)0.25*( + Var[v][k2p][j2 ][i2p] + Var[v][k2m][j2 ][i2m]
                                - Var[v][k2p][j2 ][i2m] - Var[v][k2m][j2 ][i2p] );
      } // if ( Stencil == 2 ) ... else ...


      Nume += SQR(Der2_xx) + SQR(Der2_yy) + SQR(Der2_zz) + (real)2.0*( SQR(Der2_xy) + SQR(Der2_yz) + SQR(Der2_zx) );



//    2a. denominator: gradient
      if ( Stencil == 2 )
      {
         Der1_xx = FABS( Slope[v][0][kk ][jj ][iip] ) + FABS( Slope[v][0][kk ][jj ][iim] );
         Der1_yy = FABS( Slope[v][1][kk ][jjp][ii ] ) + FABS( Slope[v][1][kk ][jjm][ii ] );
         Der1_zz = FABS( Slope[v][2][kkp][jj ][ii ] ) + FABS( Slope[v][2][kkm][jj ][ii ] );

         if ( Form == LOHNER_FLASH2 ) {
         Der1_xy = FABS( Slope[v][0][kk ][jjp][ii ] ) + FABS( Slope[v][0][kk ][jjm][ii ] );
         Der1_xz = FABS( Slope[v][0][kkp][jj ][ii ] ) + FABS( Slope[v][0][kkm][jj ][ii ] );

         Der1_yx = FABS( Slope[v][1][kk ][jj ][iip] ) + FABS( Slope[v][1][kk ][jj ][iim] );
         Der1_yz = FABS( Slope[v][1][kkp][jj ][ii ] ) + FABS( Slope[v][1][kkm][jj ][ii ] );

         Der1_zx = FABS( Slope[v][2][kk ][jj ][iip] ) + FABS( Slope[v][2][kk ][jj ][iim] );
         Der1_zy = FABS( Slope[v][2][kk ][jjp][ii ] ) + FABS( Slope[v][2][kk ][jjm][ii ] ); }
      }

      else // Stencil == 1
      {
         Der1_xx = MAX(  Slope[v][0][kk ][jj ][iip]   ,       Slope[v][0][kk ][jj ][iim]    );
         Der1_yy = MAX(  Slope[v][1][kk ][jjp][ii ]   ,       Slope[v][1][kk ][jjm][ii ]    );
         Der1_zz = MAX(  Slope[v][2][kkp][jj ][ii ]   ,       Slope[v][2][kkm][jj ][ii ]    );

         if ( Form == LOHNER_FLASH1 ) {
         Der1_xy = MAX(  Slope[v][0][kk ][jjp][ii ]   ,       Slope[v][0][kk ][jjm][ii ]    );
         Der1_xz = MAX(  Slope[v][0][kkp][jj ][ii ]   ,       Slope[v][0][kkm][jj ][ii ]    );

         Der1_yx = MAX(  Slope[v][1][kk ][jj ][iip]   ,       Slope[v][1][kk ][jj ][iim]    );
         Der1_yz = MAX(  Slope[v][1][kkp][jj ][ii ]   ,       Slope[v][1][kkm][jj ][ii ]    );

         Der1_zx = MAX(  Slope[v][2][kk ][jj ][iip]   ,       Slope[v][2][kk ][jj ][iim]    );
         Der1_zy = MAX(  Slope[v][2][kk ][jjp][ii ]   ,       Slope[v][2][kk ][jjm][ii ]    ); }
      } // if ( Stencil == 2 ) ... else ...


//    2b. denominator: filter
      if ( Stencil == 2 )
      {
         Filter_xx = Filter*( Ave[v][0][kk ][jj ][iip] + Ave[v][0][kk ][jj ][iim] );
         Filter_yy = Filter*( Ave[v][1][kk ][jjp][ii ] + Ave[v][1][kk ][jjm][ii ] );
         Filter_zz = Filter*( Ave[v][2][kkp][jj ][ii ] + Ave[v][2][kkm][jj ][ii ] );
         Filter_xy = Filter*( Ave[v][0][kk ][jjp][ii ] + Ave[v][0][kk ][jjm][ii ] );
         Filter_yz = Filter*( Ave[v][1][kkp][jj ][ii ] + Ave[v][1][kkm][jj ][ii ] );
         Filter_zx = Filter*( Ave[v][2][kk ][jj ][iip] + Ave[v][2][kk ][jj ][iim] );
      }

      else // Stencil == 1
      {
//       Ave[][1][][][] & Ave[][2][][][] are useless for Stencil == 1
         Filter_xx = Filter_4*( Ave[v][0][kk ][jj ][iip] + (real)2.0*Ave[v][0][kk][jj][ii] + Ave[v][0][kk ][jj ][iim] );
         Filter_yy = Filter_4*( Ave[v][0][kk ][jjp][ii ] + (real)2.0*Ave[v][0][kk][jj][ii] + Ave[v][0][kk ][jjm][ii ] );
         Filter_zz = Filter_4*( Ave[v][0][kkp][jj ][ii ] + (real)2.0*Ave[v][0][kk][jj][ii] + Ave[v][0][kkm][jj ][ii ] );

         Filter_xy = Filter_4*( + Ave[v][0][kk ][jj ][iip] + Ave[v][0][kk ][jjp][ii ]
                                + Ave[v][0][kk ][jj ][iim] + Ave[v][0][kk ][jjm][ii ] );
         Filter_yz = Filter_4*( + Ave[v][0][kk ][jjp][ii ] + Ave[v][0][kkp][jj ][ii ]
                                + Ave[v][0][kk ][jjm][ii ] + Ave[v][0][kkm][jj ][ii ] );
         Filter_zx = Filter_4*( + Ave[v][0][kkp][jj ][ii ] + Ave[v][0][kk ][jj ][iip]
                                + Ave[v][0][kkm][jj ][ii ] + Ave[v][0][kk ][jj ][iim] );
      } // if ( Stencil == 2 ) ... else ...


//    2c. denominator: sum
      if ( Form == LOHNER_FLASH1  || Form == LOHNER_FLASH2 )
         Deno += SQR( Der1_xx + Filter_xx ) + SQR( Der1_xy + Filter_xy ) + SQR( Der1_xz + Filter_zx ) +
                 SQR( Der1_yx + Filter_xy ) + SQR( Der1_yy + Filter_yy ) + SQR( Der1_yz + Filter_yz ) +
                 SQR( Der1_zx + Filter_zx ) + SQR( Der1_zy + Filter_yz ) + SQR( Der1_zz + Filter_zz );

      else // form-invariant
         Deno += (real)3.0*( SQR(Der1_xx) + SQR(Der1_yy) + SQR(Der1_zz) ) +
                 SQR(Filter_xx) + SQR(Filter_xy) + SQR(Filter_zx) +
                 SQR(Filter_xy) + SQR(Filter_yy) + SQR(Filter_yz) +
                 SQR(Filter_zx) + SQR(Filter_yz) + SQR(Filter_zz);


//    3. check the flag
#     if ( MODEL == ELBDM )
      if ( v == NVar - 1 )    // we check "SQRT(REAL^2 + IMAG^2)" in ELBDM
#     endif
      {
         if ( Stencil == 2 )
         Error = SQRT( Nume / MAX(Deno, Soften   ) );
         else
         Error = SQRT( Nume / MAX(Deno, Soften_16) );

         Flag  = Error > Threshold;

         if ( Flag )    return true;
         else
         {
            Nume = 0.0;
            Deno = 0.0;
         }
      }
   } // for (int v=0; v<NVar; v++)


   return false;

} // FUNCTION : Flag_Lohner



//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_for_Lohner
// Description :  Evaluate slopes and averages along x/y/z for the Lohner error estimator
//
// Note        :  1. This function is called in "Flag_Real" before looping over all cells in the patch in order to
//                   achieve higher performance
//                2. Evaluate slope by the discrete central difference:
//                      stencil == 2: slope_x(i,j,k) = var(i+1,j,k) - var(i-1,j,k)
//                      stencil == 1: slope_x(i,j,k) = MAX(  ABS( var(i+1,j,k) - var(i  ,j,k) ),
//                                                           ABS( var(i,  j,k) - var(i-1,j,k) ) )
//                                     ~ half of the stencil==2 case
//                3. The averages along different directions are evaluated as:
//                      stencil == 2: ave_x(i,j,k) = FABS( var(i+1,j,k) ) + FABS( var(i-1,j,k) )
//                      stencil == 1: ave_x(i,j,k) = FABS( var(i  ,j,k) )
//                                     ~ half of the stencil==2 case
//                4. Do not take into account the physical size of each cell since the Lohner error estimator
//                   is dimensionless
//                5. The sizes of the arrays (Var1D, Ave1D, Slope1D) must be NVar*( (PS1+4)^3, 3*(PS1+2)^3, 3*(PS1+2)^3 )
//
// Parameter   :  Form    : Form of the Lohner's error estimator (FLASH/form-invariant + Stencil = 1/2)
//                Var1D   : Array storing the input variables for the Lohner error estimator
//                Ave1D   : Array to store the output average of Var for the Lohner error estimator
//                Slope1D : Array to store the output slopes of Var for the Lohner error estimator
//                NVar    : Number of variables stored in "Var, Ave, and Slope"
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Prepare_for_Lohner( const OptLohnerForm_t Form, const real *Var1D, real *Ave1D, real *Slope1D, const int NVar )
{

   const int NCell  = PS1 + 4;   // size of the array Var
   const int NAve   = PS1 + 2;   // size of the array Ave
   const int NSlope = NAve;      // size of the array Slope

   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;
   real Slope_L[3], Slope_R[3];

// convert the 1D arrays
   real (*Var)     [NCell ][NCell ][NCell ] = ( real(*)   [NCell ][NCell ][NCell ] )  Var1D;
   real (*Ave  )[3][NAve  ][NAve  ][NAve  ] = ( real(*)[3][NAve  ][NAve  ][NAve  ] )  Ave1D;
   real (*Slope)[3][NSlope][NSlope][NSlope] = ( real(*)[3][NSlope][NSlope][NSlope] )Slope1D;


// stencil = 2
   if ( Form == LOHNER_FLASH2  ||  Form == LOHNER_FORM_INV2 )
   {
      for (int v=0; v<NVar; v++)    {
      for (int k=0; k<NAve; k++)    {  kk = k + 1;   kkp = kk + 1;   kkm = kk - 1;
      for (int j=0; j<NAve; j++)    {  jj = j + 1;   jjp = jj + 1;   jjm = jj - 1;
      for (int i=0; i<NAve; i++)    {  ii = i + 1;   iip = ii + 1;   iim = ii - 1;

         Ave  [v][0][k][j][i] = FABS( Var[v][kk ][jj ][iip] ) + FABS( Var[v][kk ][jj ][iim] );
         Ave  [v][1][k][j][i] = FABS( Var[v][kk ][jjp][ii ] ) + FABS( Var[v][kk ][jjm][ii ] );
         Ave  [v][2][k][j][i] = FABS( Var[v][kkp][jj ][ii ] ) + FABS( Var[v][kkm][jj ][ii ] );

         Slope[v][0][k][j][i] =       Var[v][kk ][jj ][iip]   -       Var[v][kk ][jj ][iim];
         Slope[v][1][k][j][i] =       Var[v][kk ][jjp][ii ]   -       Var[v][kk ][jjm][ii ];
         Slope[v][2][k][j][i] =       Var[v][kkp][jj ][ii ]   -       Var[v][kkm][jj ][ii ];

      }}}} // v,k,j,i
   }


// stencil = 1
   else
   {
      for (int v=0; v<NVar; v++)    {
      for (int k=0; k<NAve; k++)    {  kk = k + 1;   kkp = kk + 1;   kkm = kk - 1;
      for (int j=0; j<NAve; j++)    {  jj = j + 1;   jjp = jj + 1;   jjm = jj - 1;
      for (int i=0; i<NAve; i++)    {  ii = i + 1;   iip = ii + 1;   iim = ii - 1;

//       Ave_y/z are useless when stencil = 1
         Ave  [v][0][k][j][i] = FABS( Var[v][kk ][jj ][ii ] );
         /*
         Ave  [v][1][k][j][i] = FABS( Var[v][kk ][jj ][ii ] );
         Ave  [v][2][k][j][i] = FABS( Var[v][kk ][jj ][ii ] );
         */

//       Slope is positive-definite when stencil = 1 (ok since we don't use Slope to estimate Der2)
         Slope_L[0]           = FABS( Var[v][kk ][jj ][ii ]   -       Var[v][kk ][jj ][iim] );
         Slope_L[1]           = FABS( Var[v][kk ][jj ][ii ]   -       Var[v][kk ][jjm][ii ] );
         Slope_L[2]           = FABS( Var[v][kk ][jj ][ii ]   -       Var[v][kkm][jj ][ii ] );
         Slope_R[0]           = FABS( Var[v][kk ][jj ][iip]   -       Var[v][kk ][jj ][ii ] );
         Slope_R[1]           = FABS( Var[v][kk ][jjp][ii ]   -       Var[v][kk ][jj ][ii ] );
         Slope_R[2]           = FABS( Var[v][kkp][jj ][ii ]   -       Var[v][kk ][jj ][ii ] );

         for (int d=0; d<3; d++)    Slope[v][d][k][j][i] = MAX( Slope_R[d], Slope_L[d] );

      }}}} // v,k,j,i
   } // if ( Form == LOHNER_FLASH2  ||  Form == LOHNER_FORM_INV2 ) ... else ...


} // FUNCTION : Prepare_for_Lohner

