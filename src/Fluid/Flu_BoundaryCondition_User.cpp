#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void BC_User( real fluid[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*BC_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] ) = BC_User;

#ifdef MHD
extern void (*BC_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                   const int lv, double AuxArray[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by Flu_BoundaryCondition_User() using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Always return NCOMP_TOTAL variables
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//                   --> Just like the fluid initialization routine
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_User( real fluid[], const double x, const double y, const double z, const double Time,
              const int lv, double AuxArray[] )
{

// put your B.C. here
// ##########################################################################################################
// Example 1 : set to time-independent values for HYDRO
   /*
   const double C[3] = { 0.5*amr->BoxSize[0],
                         0.5*amr->BoxSize[1],
                         0.5*amr->BoxSize[2] };
   const real Height = 100.0;
   const real Width  =  64.0;
   const real Gamma2 = real( 1.0/GAMMA/(GAMMA-1.0) );
   const real Cs     = 1.0;
   const real Rho0   = 1.0;

   fluid[DENS] = Rho0 + Height*EXP(  -( SQR(x-C[0]) + SQR(y-C[1]) + SQR(z-C[2]) ) / SQR(Width)  );
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = Cs*Cs*fluid[DENS]*Gamma2 + (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

// remember to set passive scalars as well
   fluid[EINT] = XXX;
   */


// Example 2 : set to time-dependent values for HYDRO

// ##########################################################################################################

} // FUNCTION : BC_User



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_BoundaryCondition_User
// Description :  Fill up the ghost-zone values by the user-specified boundary condition
//
// Note        :  1. Work for Prepare_PatchData(), InterpolateGhostZone(), Refine(), and LB_Refine_GetNewRealPatchList()
//                2. The function pointer "BC_User_Ptr" points to BC_User() by default but may be overwritten
//                   by various test problem initializers
//                3. User-defined boundary conditions for the magnetic field are set in MHD_BoundaryCondition_User()
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables are NOT included)
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                Time           : Current physical time
//                dh             : Cell size
//                Corner         : Physcial coordinates at the center of the cell (0,0,0) --> Array[0]
//                TVar           : Target variables to be prepared --> only used for preparing the derived variables
//                lv             : Refinement level
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void Flu_BoundaryCondition_User( real *Array, const int NVar_Flu, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                                 const int TFluVarIdxList[], const double Time, const double dh, const double *Corner,
                                 const int TVar, const int lv )
{

// check
   if ( BC_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "BC_User_Ptr == NULL !!\n" );


   const double x0 = Corner[0] + (double)Idx_Start[0]*dh;   // starting x,y,z coordinates
   const double y0 = Corner[1] + (double)Idx_Start[1]*dh;
   const double z0 = Corner[2] + (double)Idx_Start[2]*dh;

#  if   ( MODEL == HYDRO )
#  ifdef MHD
   const double dh_2             = 0.5*dh;
#  endif
   const bool   CheckMinPres_Yes = true;
   const real   Gamma_m1         = GAMMA - (real)1.0;
   const bool   PrepVx           = ( TVar & _VELX ) ? true : false;
   const bool   PrepVy           = ( TVar & _VELY ) ? true : false;
   const bool   PrepVz           = ( TVar & _VELZ ) ? true : false;
   const bool   PrepPres         = ( TVar & _PRES ) ? true : false;
   const bool   PrepTemp         = ( TVar & _TEMP ) ? true : false;

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif


// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   int    i, j, k, v2;
   real   BVal[NCOMP_TOTAL];
   double x, y, z;

   for (k=Idx_Start[2], z=z0; k<=Idx_End[2]; k++, z+=dh)
   for (j=Idx_Start[1], y=y0; j<=Idx_End[1]; j++, y+=dh)
   for (i=Idx_Start[0], x=x0; i<=Idx_End[0]; i++, x+=dh)
   {
//    1. primary variables
//    get the boundary values of all NCOMP_TOTAL fields
      BC_User_Ptr( BVal, x, y, z, Time, lv, NULL );

//    add the magnetic energy for MHD
#     if ( MODEL == HYDRO )
#     ifdef MHD
      real EngyB, BxL, BxR, Bx, ByL, ByR, By, BzL, BzR, Bz, B3v[NCOMP_MAG];

      BC_BField_User_Ptr( B3v, x-dh_2, y,      z,      Time, lv, NULL );   BxL = B3v[MAGX];
      BC_BField_User_Ptr( B3v, x+dh_2, y,      z,      Time, lv, NULL );   BxR = B3v[MAGX];
      BC_BField_User_Ptr( B3v, x,      y-dh_2, z,      Time, lv, NULL );   ByL = B3v[MAGY];
      BC_BField_User_Ptr( B3v, x,      y+dh_2, z,      Time, lv, NULL );   ByR = B3v[MAGY];
      BC_BField_User_Ptr( B3v, x,      y,      z-dh_2, Time, lv, NULL );   BzL = B3v[MAGZ];
      BC_BField_User_Ptr( B3v, x,      y,      z+dh_2, Time, lv, NULL );   BzR = B3v[MAGZ];

      Bx    = (real)0.5*( BxL + BxR );
      By    = (real)0.5*( ByL + ByR );
      Bz    = (real)0.5*( BzL + BzR );
      EngyB = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );

      BVal[ENGY] += EngyB;

#     else
      const real EngyB = NULL_REAL;
#     endif
#     endif // #ifdef ( MODEL == HYDRO )

//    store results to the output array
      for (int v=0; v<NVar_Flu; v++)   Array3D[v][k][j][i] = BVal[ TFluVarIdxList[v] ];


//    2. derived variables
      v2 = NVar_Flu;

#     if   ( MODEL == HYDRO )
      if ( PrepVx   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMX] / BVal[DENS];
      if ( PrepVy   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMY] / BVal[DENS];
      if ( PrepVz   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMZ] / BVal[DENS];
      if ( PrepPres )   Array3D[ v2 ++ ][k][j][i] = Hydro_GetPressure( BVal[DENS], BVal[MOMX], BVal[MOMY], BVal[MOMZ], BVal[ENGY],
                                                                       Gamma_m1, CheckMinPres_Yes, MIN_PRES, EngyB );
      if ( PrepTemp )   Array3D[ v2 ++ ][k][j][i] = Hydro_GetTemperature( BVal[DENS], BVal[MOMX], BVal[MOMY], BVal[MOMZ], BVal[ENGY],
                                                                          Gamma_m1, CheckMinPres_Yes, MIN_PRES, EngyB );

#     elif ( MODEL == ELBDM )
//    no derived variables yet

#     else
#     error : unsupported MODEL !!
#     endif
   } // k,j,i

} // FUNCTION : Flu_BoundaryCondition_User
