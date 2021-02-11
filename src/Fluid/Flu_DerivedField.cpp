#include "GAMER.h"




#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_DerivedField_DivVel
// Description :  Define the built-in derived field "divergence(velocity)"
//
// Note        :  1. FluIn[] must include all NCOMP_TOTAL fluid fields
//             :  2. Array size:
//                      Out  [] : NFieldOut*(NCellInX-2*NGhost)*(NCellInY-2*NGhost)*(NCellInZ-2*NGhost)
//                      FluIn[] : NCOMP_TOTAL*NCellInX*NCellInY*NCellInZ
//                      MagIn[] : NCOMP_MAG  *NCellInX*NCellInY*NCellInZ
//                3. NFieldOut must NOT be larger than DER_NOUT_MAX defined in Macro.h (default = 10)
//
// Parameter   :  Out          : Output array
//                FluIn        : Input fluid array
//                MagIn        : Input magnetic field array (cell-centered)
//                NFieldOut    : Number of output fields
//                NCellInX/Y/Z : Size of FluIn[] and MagIn[] along each direction
//                               --> Including ghost zones
//                NGhost       : Number of ghost cells
//                dh           : Cell width
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
void Flu_DerivedField_DivVel( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                              const int NCellInX, const int NCellInY, const int NCellInZ,
                              const int NGhost, const double dh )
{

// determine the output array size
   const int NCellOutX = NCellInX - 2*NGhost;
   const int NCellOutY = NCellInY - 2*NGhost;
   const int NCellOutZ = NCellInZ - 2*NGhost;


// check
#  ifdef GAMER_DEBUG
// general
   if ( Out   == NULL )             Aux_Error( ERROR_INFO, "Out == NULL !!\n" );
   if ( FluIn == NULL )             Aux_Error( ERROR_INFO, "FluIn == NULL !!\n" );
#  ifdef MHD
// if ( MagIn == NULL )             Aux_Error( ERROR_INFO, "MagIn == NULL !!\n" );
#  endif
   if ( NCellOutX <= 0 )            Aux_Error( ERROR_INFO, "NCellOutX (%d) <= 0 !!\n", NCellOutX );
   if ( NCellOutY <= 0 )            Aux_Error( ERROR_INFO, "NCellOutY (%d) <= 0 !!\n", NCellOutY );
   if ( NCellOutZ <= 0 )            Aux_Error( ERROR_INFO, "NCellOutZ (%d) <= 0 !!\n", NCellOutZ );
   if ( NFieldOut > DER_NOUT_MAX )  Aux_Error( ERROR_INFO, "NFieldOut (%d) > DER_NOUT_MAX (%d) !!\n", NFieldOut, DER_NOUT_MAX );

// specific to this function
   if ( NGhost < 1 )                Aux_Error( ERROR_INFO, "NGhost (%d) < 1 !!\n", NGhost );
   if ( NFieldOut != 1 )            Aux_Error( ERROR_INFO, "NFieldOut (%d) != 1 !!\n", NFieldOut );
#  endif // #ifdef GAMER_DEBUG


// 1D arrays -> 3D arrays
   real (*FluIn3D)[NCellInZ ][NCellInY ][NCellInX ] = ( real (*)[NCellInZ ][NCellInY ][NCellInX ] )FluIn;
   real (*Out3D  )[NCellOutZ][NCellOutY][NCellOutX] = ( real (*)[NCellOutZ][NCellOutY][NCellOutX] )Out;


// fill in the output array
   const real _2dh = (real)0.5 / (real)dh;
   real vx_p, vx_m, vy_p, vy_m, vz_p, vz_m;

   for (int ko=0; ko<NCellOutZ; ko++)  {  const int ki = ko + NGhost;
   for (int jo=0; jo<NCellOutY; jo++)  {  const int ji = jo + NGhost;
   for (int io=0; io<NCellOutX; io++)  {  const int ii = io + NGhost;

      vx_p = FluIn3D[MOMX][ki  ][ji  ][ii+1] / FluIn3D[DENS][ki  ][ji  ][ii+1];
      vx_m = FluIn3D[MOMX][ki  ][ji  ][ii-1] / FluIn3D[DENS][ki  ][ji  ][ii-1];
      vy_p = FluIn3D[MOMY][ki  ][ji+1][ii  ] / FluIn3D[DENS][ki  ][ji+1][ii  ];
      vy_m = FluIn3D[MOMY][ki  ][ji-1][ii  ] / FluIn3D[DENS][ki  ][ji-1][ii  ];
      vz_p = FluIn3D[MOMZ][ki+1][ji  ][ii  ] / FluIn3D[DENS][ki+1][ji  ][ii  ];
      vz_m = FluIn3D[MOMZ][ki-1][ji  ][ii  ] / FluIn3D[DENS][ki-1][ji  ][ii  ];

      Out3D[0][ko][jo][io] = _2dh*( (vx_p-vx_m) + (vy_p-vy_m) + (vz_p-vz_m) );

   }}} // k,j,i

} // FUNCTION : Flu_DerivedField_DivVel
#endif // #if ( MODEL == HYDRO )
