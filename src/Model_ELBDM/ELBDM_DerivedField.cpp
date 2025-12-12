#include "GAMER.h"

#if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_differentiation
// Description :  Calculate gradients and Laplacians for ELBDM fields using Richardson extrapolation
//
// Note        :  1. Uses Richardson extrapolation for derivatives
//                2. Standard second-order central difference code is commented out but preserved (less accurate)
//                3. Calculates derivatives for density, real part, imaginary part, and sqrt(density)
//
// Parameter   :  FieldIn     : Input field array
//                k,j,i       : Grid indices of the target cell
//                GradD       : Output array for density gradients [3]
//                GradR       : Output array for real part gradients [3]
//                GradI       : Output array for imaginary part gradients [3]
//                LapD        : Output Laplacian of density
//                LapR        : Output Laplacian of real part
//                LapI        : Output Laplacian of imaginary part
//                Lapf        : Output Laplacian of sqrt(density)
//                dh          : Grid spacing
//                NGhost      : Number of ghost cells
//
// Return      :  GradD, GradR, GradI, LapD, LapR, LapI, Lapf
//-------------------------------------------------------------------------------------------------------
void ELBDM_differentiation( const real FieldIn[], const int k, const int j, const int i,
                            real GradD[3], real GradR[3], real GradI[3],
                            real &LapD, real &LapR, real &LapI, real &Lapf, const real dh, const int NGhost )
{
// finite difference coefficients
   const real _2dh   = 0.5/dh;
   const real _dh2   = 1.0/SQR(dh);

// neighboring grid indices
   const int  im     = i - 1;
   const int  ip     = i + 1;
   const int  jm     = j - 1;
   const int  jp     = j + 1;
   const int  km     = k - 1;
   const int  kp     = k + 1;
   const int  imm    = i - 2;
   const int  ipp    = i + 2;
   const int  jmm    = j - 2;
   const int  jpp    = j + 2;
   const int  kmm    = k - 2;
   const int  kpp    = k + 2;

// cast input array to 3D field structure
   typedef real (*Field3d)[PS1+2*NGhost ][PS1+2*NGhost ][PS1+2*NGhost ];
   Field3d Field = ( Field3d )FieldIn;

// field values at the target cell
   real Dens    = Field[DENS    ][k][j][i];
   real Real    = Field[REAL    ][k][j][i];
   real Imag    = Field[IMAG    ][k][j][i];


// Richardson extrapolation 2nd order
   GradD[0] = (    4*_2dh*( Field[DENS][k  ][j  ][ip ] - Field[DENS][k  ][j  ][im ] )
               - 0.5*_2dh*( Field[DENS][k  ][j  ][ipp] - Field[DENS][k  ][j  ][imm] ))/3;
   GradD[1] = (    4*_2dh*( Field[DENS][k  ][jp ][i  ] - Field[DENS][k  ][jm ][i  ] )
               - 0.5*_2dh*( Field[DENS][k  ][jpp][i  ] - Field[DENS][k  ][jmm][i  ] ))/3;
   GradD[2] = (    4*_2dh*( Field[DENS][kp ][j  ][i  ] - Field[DENS][km ][j  ][i  ] )
               - 0.5*_2dh*( Field[DENS][kpp][j  ][i  ] - Field[DENS][kmm][j  ][i  ] ))/3;

   GradR[0] = (    4*_2dh*( Field[REAL][k  ][j  ][ip ] - Field[REAL][k  ][j  ][im ] )
               - 0.5*_2dh*( Field[REAL][k  ][j  ][ipp] - Field[REAL][k  ][j  ][imm] ))/3;
   GradR[1] = (    4*_2dh*( Field[REAL][k  ][jp ][i  ] - Field[REAL][k  ][jm ][i  ] )
               - 0.5*_2dh*( Field[REAL][k  ][jpp][i  ] - Field[REAL][k  ][jmm][i  ] ))/3;
   GradR[2] = (    4*_2dh*( Field[REAL][kp ][j  ][i  ] - Field[REAL][km ][j  ][i  ] )
               - 0.5*_2dh*( Field[REAL][kpp][j  ][i  ] - Field[REAL][kmm][j  ][i  ] ))/3;

   GradI[0] = (    4*_2dh*( Field[IMAG][k  ][j  ][ip ] - Field[IMAG][k  ][j  ][im ] )
               - 0.5*_2dh*( Field[IMAG][k  ][j  ][ipp] - Field[IMAG][k  ][j  ][imm] ))/3;
   GradI[1] = (    4*_2dh*( Field[IMAG][k  ][jp ][i  ] - Field[IMAG][k  ][jm ][i  ] )
               - 0.5*_2dh*( Field[IMAG][k  ][jpp][i  ] - Field[IMAG][k  ][jmm][i  ] ))/3;
   GradI[2] = (    4*_2dh*( Field[IMAG][kp ][j  ][i  ] - Field[IMAG][km ][j  ][i  ] )
               - 0.5*_2dh*( Field[IMAG][kpp][j  ][i  ] - Field[IMAG][kmm][j  ][i  ] ))/3;

   LapD    =  (    4*_dh2*( Field[DENS][k  ][j  ][ip ] + Field[DENS][k  ][jp ][i  ] + Field[DENS][kp ][j  ][i  ] +
                              Field[DENS][k  ][j  ][im ] + Field[DENS][k  ][jm ][i  ] + Field[DENS][km ][j  ][i  ] - 6.0*Dens )
               - 0.25*_dh2*( Field[DENS][k  ][j  ][ipp] + Field[DENS][k  ][jpp][i  ] + Field[DENS][kpp][j  ][i  ] +
                              Field[DENS][k  ][j  ][imm] + Field[DENS][k  ][jmm][i  ] + Field[DENS][kmm][j  ][i  ] - 6.0*Dens ) )/3;

   LapR     = (    4*_dh2*( Field[REAL][k  ][j  ][ip ] + Field[REAL][k  ][jp ][i  ] + Field[REAL][kp ][j  ][i  ] +
                              Field[REAL][k  ][j  ][im ] + Field[REAL][k  ][jm ][i  ] + Field[REAL][km ][j  ][i  ] - 6.0*Real )
               - 0.25*_dh2*( Field[REAL][k  ][j  ][ipp] + Field[REAL][k  ][jpp][i  ] + Field[REAL][kpp][j  ][i  ] +
                              Field[REAL][k  ][j  ][imm] + Field[REAL][k  ][jmm][i  ] + Field[REAL][kmm][j  ][i  ] - 6.0*Real ) )/3;

   LapI     = (    4*_dh2*( Field[IMAG][k  ][j  ][ip ] + Field[IMAG][k  ][jp ][i  ] + Field[IMAG][kp ][j  ][i  ] +
                              Field[IMAG][k  ][j  ][im ] + Field[IMAG][k  ][jm ][i  ] + Field[IMAG][km ][j  ][i  ] - 6.0*Imag )
               - 0.25*_dh2*( Field[IMAG][k  ][j  ][ipp] + Field[IMAG][k  ][jpp][i  ] + Field[IMAG][kpp][j  ][i  ] +
                              Field[IMAG][k  ][j  ][imm] + Field[IMAG][k  ][jmm][i  ] + Field[IMAG][kmm][j  ][i  ] - 6.0*Imag ) )/3;

   Lapf     = (    4*_dh2*( SQRT(Field[DENS][k  ][j  ][ip ]) + SQRT(Field[DENS][k  ][jp ][i  ]) + SQRT(Field[DENS][kp ][j  ][i  ]) +
                              SQRT(Field[DENS][k  ][j  ][im ]) + SQRT(Field[DENS][k  ][jm ][i  ]) + SQRT(Field[DENS][km ][j  ][i  ]) - 6.0*SQRT(Dens) )
               - 0.25*_dh2*( SQRT(Field[DENS][k  ][j  ][ipp]) + SQRT(Field[DENS][k  ][jpp][i  ]) + SQRT(Field[DENS][kpp][j  ][i  ]) +
                              SQRT(Field[DENS][k  ][j  ][imm]) + SQRT(Field[DENS][k  ][jmm][i  ]) + SQRT(Field[DENS][kmm][j  ][i  ]) - 6.0*SQRT(Dens) ) )/3;

/*
   // calculate gradients using second-order central differences (commented out but preserved for reference)
   GradD[0] = _2dh*( Field[DENS][k ][j ][ip] - Field[DENS][k ][j ][im] );
   GradD[1] = _2dh*( Field[DENS][k ][jp][i ] - Field[DENS][k ][jm][i ] );
   GradD[2] = _2dh*( Field[DENS][kp][j ][i ] - Field[DENS][km][j ][i ] );

   GradR[0] = _2dh*( Field[REAL][k ][j ][ip] - Field[REAL][k ][j ][im] );
   GradR[1] = _2dh*( Field[REAL][k ][jp][i ] - Field[REAL][k ][jm][i ] );
   GradR[2] = _2dh*( Field[REAL][kp][j ][i ] - Field[REAL][km][j ][i ] );

   GradI[0] = _2dh*( Field[IMAG][k ][j ][ip] - Field[IMAG][k ][j ][im] );
   GradI[1] = _2dh*( Field[IMAG][k ][jp][i ] - Field[IMAG][k ][jm][i ] );
   GradI[2] = _2dh*( Field[IMAG][kp][j ][i ] - Field[IMAG][km][j ][i ] );

   // calculate Laplacians using second-order central differences (commented out but preserved for reference)
   LapD    = ( Field[DENS][k ][j ][ip] + Field[DENS][k ][jp][i ] + Field[DENS][kp][j ][i ] +
               Field[DENS][k ][j ][im] + Field[DENS][k ][jm][i ] + Field[DENS][km][j ][i ] -
               6.0*Dens )*_dh2;
   LapR     = ( Field[REAL][k ][j ][ip] + Field[REAL][k ][jp][i ] + Field[REAL][kp][j ][i ] +
               Field[REAL][k ][j ][im] + Field[REAL][k ][jm][i ] + Field[REAL][km][j ][i ] -
               6.0*Real )*_dh2;
   LapI     = ( Field[IMAG][k ][j ][ip] + Field[IMAG][k ][jp][i ] + Field[IMAG][kp][j ][i ] +
               Field[IMAG][k ][j ][im] + Field[IMAG][k ][jm][i ] + Field[IMAG][km][j ][i ] -
               6.0*Imag )*_dh2;
   Lapf     = ( SQRT(Field[DENS][k ][j ][ip]) + SQRT(Field[DENS][k ][jp][i ]) + SQRT(Field[DENS][kp][j ][i ]) +
               SQRT(Field[DENS][k ][j ][im]) + SQRT(Field[DENS][k ][jm][i ]) + SQRT(Field[DENS][km][j ][i ]) -
               6.0*SQRT(Dens) )*_dh2;
*/

} // FUNCTION : ELBDM_differentiation



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_DerivedField
// Description :  Calculate derived fields for ELBDM model
//
// Note        :  1. FieldID=0: bulk velocity v = (R*dI/dx - I*dR/dx)/rho*hbar/m
//                2. FieldID=1: thermal velocity w = (R*dR/dx + I*dI/dx)/rho*hbar/m
//                3. FieldID=2: quantum pressure Q = -1/2*(Laplacian(f)/f)*hbar^2/m^2
//                4. FieldID=3: stress tensor Sigma_(i,j) components
//
// Parameter   :  ELBDMOut    : Output array for the derived field
//                ELBDMIn     : Input ELBDM field array
//                FieldID     : Target derived field ID (0-3)
//                direction   : Spatial direction or tensor component index
//                NGhost      : Number of ghost cells
//                dh          : Grid spacing
//
// Return      :  ELBDMOut
//-------------------------------------------------------------------------------------------------------
void ELBDM_DerivedField( real ELBDMOut[], const real ELBDMIn[], int  FieldID,
                     int direction, const int NGhost, const real dh )
{

// check
   if ( FieldID < 0 || FieldID > 3 )
      Aux_Error( ERROR_INFO, "incorrect FieldID (%d) !!\n", FieldID );

   if ( direction < 0 || direction > 5 )
      Aux_Error( ERROR_INFO, "incorrect direction (%d) !!\n", direction );

// cast input/output arrays to 3D structures
   typedef real (*ELBDM_in)[PS1+2*NGhost ][PS1+2*NGhost ][PS1+2*NGhost ];
   typedef real (*ELBDM_out)[PS1 ][PS1 ][PS1 ];
   ELBDM_in ELBDMIn3d = ( ELBDM_in )ELBDMIn;
   ELBDM_out ELBDMOut3d = ( ELBDM_out )ELBDMOut;

// conversion factor: hbar/m
   const real _Eta   = 1.0/ELBDM_ETA;

// loop over all interior cells
   for (int k=NGhost; k<PS1+NGhost; k++)
   {
      for (int j=NGhost; j<PS1+NGhost; j++)
      {
         for (int i=NGhost; i<PS1+NGhost; i++)
         {
            real GradD[3], GradR[3], GradI[3];
            real LapD, LapR, LapI, Lapf;

//          field values at current cell
            real Dens    = ELBDMIn3d[DENS    ][k][j][i];
            real Real    = ELBDMIn3d[REAL    ][k][j][i];
            real Imag    = ELBDMIn3d[IMAG    ][k][j][i];
            real _Dens = 1/Dens;

//          calculate derivatives
            ELBDM_differentiation( ELBDMIn, k, j, i, GradD, GradR, GradI, LapD, LapR, LapI, Lapf, dh, NGhost );

//          calculate the requested derived field
            if ( FieldID == 0)
            {
//             bulk velocity: v = (R*dI/dx - I*dR/dx)/rho*hbar/m
               ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = _Eta*_Dens*( Real*GradI[direction] - Imag*GradR[direction] );
            }
            else if ( FieldID == 1)
            {
//             thermal velocity: w = (R*dR/dx + I*dI/dx)/rho*hbar/m
               ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = _Eta*_Dens*( Real*GradR[direction] + Imag*GradI[direction] );
            }
            else if ( FieldID == 2)
            {
//             quantum pressure: Q = -1/2*(Laplacian(f)/f)*hbar^2/m^2
               ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = FABS(-0.5*Lapf*_Eta*_Eta*SQRT(_Dens));
            }
            else if ( FieldID == 3)
            {
//             stress tensor: Sigma_(i,j) = 0.25*hbar^2/m^2 (1/rho*drho/dx_i*drho/dx_j - delta_ij*laplacian(rho))
               if (direction <3)        ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = FABS(0.25*_Eta*_Eta*( _Dens*GradD[direction]*GradD[direction] - LapD ));   // diagonal components
               else if (direction == 3) ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = FABS(0.25*_Eta*_Eta*( _Dens*GradD[0]*GradD[1] ));   // xy
               else if (direction == 4) ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = FABS(0.25*_Eta*_Eta*( _Dens*GradD[1]*GradD[2] ));   // yz
               else if (direction == 5) ELBDMOut3d[0][k-NGhost][j-NGhost][i-NGhost] = FABS(0.25*_Eta*_Eta*( _Dens*GradD[0]*GradD[2] ));   // xz
            }

         } // for i
      } // for j
   } // for k

} // FUNCTION : ELBDM_DerivedField

#endif // #if ( MODEL == ELBDM )

