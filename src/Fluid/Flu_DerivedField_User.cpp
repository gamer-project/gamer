#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Flu_DerivedField_User_Template( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                                            const int NCellInX, const int NCellInY, const int NCellInZ,
                                            const int NGhost, const double dh );

// this function pointer must be set by Init_DerivedField_User_Ptr()
void (*Flu_DerivedField_User_Ptr)( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                                   const int NCellInX, const int NCellInY, const int NCellInZ,
                                   const int NGhost, const double dh ) = NULL;

// this function pointer must be set by a test problem initializer
void (*Init_DerivedField_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_DerivedField_User_Template
// Description :  Template of user-defined derived fields
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
void Flu_DerivedField_User_Template( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                                     const int NCellInX, const int NCellInY, const int NCellInZ,
                                     const int NGhost, const double dh )
{

// determine the output array size
   const int NCellOutX = NCellInX - 2*NGhost;
   const int NCellOutY = NCellInY - 2*NGhost;
   const int NCellOutZ = NCellInZ - 2*NGhost;


// check
#  ifdef GAMER_DEBUG
   if ( Out   == NULL )             Aux_Error( ERROR_INFO, "Out == NULL !!\n" );
   if ( FluIn == NULL )             Aux_Error( ERROR_INFO, "FluIn == NULL !!\n" );
#  ifdef MHD
   if ( MagIn == NULL )             Aux_Error( ERROR_INFO, "MagIn == NULL !!\n" );
#  endif
   if ( NCellOutX <= 0 )            Aux_Error( ERROR_INFO, "NCellOutX (%d) <= 0 !!\n", NCellOutX );
   if ( NCellOutY <= 0 )            Aux_Error( ERROR_INFO, "NCellOutY (%d) <= 0 !!\n", NCellOutY );
   if ( NCellOutZ <= 0 )            Aux_Error( ERROR_INFO, "NCellOutZ (%d) <= 0 !!\n", NCellOutZ );
   if ( NFieldOut > DER_NOUT_MAX )  Aux_Error( ERROR_INFO, "NFieldOut (%d) > DER_NOUT_MAX (%d) !!\n", NFieldOut, DER_NOUT_MAX );
   if ( NGhost < 0 )                Aux_Error( ERROR_INFO, "NGhost (%d) < 0 !!\n", NGhost );
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
   /*
   real fluid[NCOMP_TOTAL];
#  ifdef MHD
   real B[NCOMP_MAG];
#  endif

   for (int ko=0; ko<NCellOutZ; ko++)  {  const int ki = ko + NGhost;
   for (int jo=0; jo<NCellOutY; jo++)  {  const int ji = jo + NGhost;
   for (int io=0; io<NCellOutX; io++)  {  const int ii = io + NGhost;

      for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = FluIn3D[v][ki][ji][ii];
#     ifdef MHD
      for (int v=0; v<NCOMP_MAG; v++)     B    [v] = MagIn3D[v][ki][ji][ii];
#     endif

      Out3D[0][ko][jo][io] = ...;

   }}} // k,j,i
   */

} // FUNCTION : Flu_DerivedField_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_DerivedField_User_Template
// Description :  Initialize the user-defined derived fields
//
// Note        :  1. Set the following variables
//                   a. UserDerField_Num
//                   b. UserDerField_Label
//                   c. UserDerField_Unit
//                   d. Flu_DerivedField_User_Ptr
//                2. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  See the Note above
//-------------------------------------------------------------------------------------------------------
void Init_DerivedField_User_Template()
{

   /*
// 1. set the number of derived fields
//    --> must not exceed DER_NOUT_MAX defined in Macro.h (default = 10)
   UserDerField_Num = 2;


// 2. set the field labels
   UserDerField_Label = new char [UserDerField_Num][MAX_STRING];

   sprintf( UserDerField_Label[0], "%s", "UserDerField_Label0" );
   sprintf( UserDerField_Label[1], "%s", "UserDerField_Label1" );


// 3. set the field units
   UserDerField_Unit = new char [UserDerField_Num][MAX_STRING];

   sprintf( UserDerField_Unit[0], "%s", "UserDerField_Unit0" );
   sprintf( UserDerField_Unit[1], "%s", "UserDerField_Unit1" );


// 4. set the major function pointer
   Flu_DerivedField_User_Ptr = Flu_DerivedField_User_Template;
   */

} // FUNCTION : Init_DerivedField_User_Template
