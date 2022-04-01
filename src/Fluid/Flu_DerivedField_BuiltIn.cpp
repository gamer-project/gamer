#include "GAMER.h"




#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_DerivedField_DivVel
// Description :  Define the built-in derived field "divergence(velocity)"
//
// Note        :  1. Array size:
//                      Out  [] : NFieldOut*(NCellInX-2*NGhost)*(NCellInY-2*NGhost)*(NCellInZ-2*NGhost)
//                      FluIn[] : NCOMP_TOTAL*NCellInX*NCellInY*NCellInZ
//                      MagIn[] : NCOMP_MAG  *NCellInX*NCellInY*NCellInZ
//                2. FluIn[] stores all the conserved variables and passive scalars
//                   --> But this function only uses density and momentum
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
   typedef real (*vla_in)[NCellInZ ][NCellInY ][NCellInX ];
   typedef real (*vla_out)[NCellOutZ][NCellOutY][NCellOutX];
   vla_in FluIn3D = ( vla_in )FluIn;
   vla_out Out3D = ( vla_out )Out;


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



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_DerivedField_Mach
// Description :  Define the built-in derived field "Mach number"
//
// Note        :  1. Array size:
//                      Out  [] : NFieldOut*(NCellInX-2*NGhost)*(NCellInY-2*NGhost)*(NCellInZ-2*NGhost)
//                      FluIn[] : NCOMP_TOTAL*NCellInX*NCellInY*NCellInZ
//                      MagIn[] : NCOMP_MAG  *NCellInX*NCellInY*NCellInZ
//                2. FluIn[] stores all the conserved variables and passive scalars
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
void Flu_DerivedField_Mach( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
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
   if ( MagIn == NULL )             Aux_Error( ERROR_INFO, "MagIn == NULL !!\n" );
#  endif
   if ( NCellOutX <= 0 )            Aux_Error( ERROR_INFO, "NCellOutX (%d) <= 0 !!\n", NCellOutX );
   if ( NCellOutY <= 0 )            Aux_Error( ERROR_INFO, "NCellOutY (%d) <= 0 !!\n", NCellOutY );
   if ( NCellOutZ <= 0 )            Aux_Error( ERROR_INFO, "NCellOutZ (%d) <= 0 !!\n", NCellOutZ );
   if ( NFieldOut > DER_NOUT_MAX )  Aux_Error( ERROR_INFO, "NFieldOut (%d) > DER_NOUT_MAX (%d) !!\n", NFieldOut, DER_NOUT_MAX );

// specific to this function
   if ( NGhost < 0 )                Aux_Error( ERROR_INFO, "NGhost (%d) < 0 !!\n", NGhost );
   if ( NFieldOut != 1 )            Aux_Error( ERROR_INFO, "NFieldOut (%d) != 1 !!\n", NFieldOut );
#  endif // #ifdef GAMER_DEBUG


// 1D arrays -> 3D arrays
   typedef real (*vla_in)[NCellInZ ][NCellInY ][NCellInX ];
   vla_in FluIn3D = ( vla_in )FluIn;
#  ifdef MHD
   vla_in MagIn3D = ( vla_in )MagIn;
#  endif
   typedef real (*vla_out)[NCellOutZ][NCellOutY][NCellOutX];
   vla_out Out3D = ( vla_out )Out;


// fill in the output array
   const bool CheckMinPres_Yes = true;

   real fluid[NCOMP_TOTAL], _Dens, Vx, Vy, Vz, V2, Emag, Pres, Cs2, Mach;
#  ifdef MHD
   real B[NCOMP_MAG];
#  endif

   for (int ko=0; ko<NCellOutZ; ko++)  {  const int ki = ko + NGhost;
   for (int jo=0; jo<NCellOutY; jo++)  {  const int ji = jo + NGhost;
   for (int io=0; io<NCellOutX; io++)  {  const int ii = io + NGhost;

      for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = FluIn3D[v][ki][ji][ii];

      _Dens = (real)1.0 / fluid[DENS];
      Vx    = fluid[MOMX] * _Dens;
      Vy    = fluid[MOMY] * _Dens;
      Vz    = fluid[MOMZ] * _Dens;
      V2    = SQR( Vx ) + SQR( Vy ) + SQR( Vz );

#     ifdef MHD
      for (int v=0; v<NCOMP_MAG; v++)  B[v] = MagIn3D[v][ki][ji][ii];
      Emag  = (real)0.5*(  SQR( B[MAGX] ) + SQR( B[MAGY] ) + SQR( B[MAGZ] )  );
#     else
      Emag  = NULL_REAL;
#     endif

      Pres  = Hydro_Con2Pres( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], fluid+NCOMP_FLUID,
                              CheckMinPres_Yes, MIN_PRES, Emag,
                              EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
      Cs2   = EoS_DensPres2CSqr_CPUPtr( fluid[DENS], Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                        h_EoS_Table );
      Mach  = SQRT( V2 / Cs2 );

      Out3D[0][ko][jo][io] = Mach;

   }}} // k,j,i

} // FUNCTION : Flu_DerivedField_Mach
#endif // #if ( MODEL == HYDRO )
