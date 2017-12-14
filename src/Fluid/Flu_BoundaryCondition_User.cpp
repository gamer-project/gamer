#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void BC_User( real fluid[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*BC_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] ) = BC_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Always return NCOMP_TOTAL variables
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
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

// please put your B.C. here
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
// Note        :  1. Work for the functions "Prepare_PatchData, InterpolateGhostZone, Refine, LB_Refine_AllocateNewPatch"
//                2. The function pointer "BC_User_Ptr" points to "BC_User()" by default but may be overwritten
//                   by various test problem initializers
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables are NOT included)
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                Time           : Current physical time
//                dh             : Grid size
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
   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = GAMMA - (real)1.0;
   const bool PrepVx           = ( TVar & _VELX ) ? true : false;
   const bool PrepVy           = ( TVar & _VELY ) ? true : false;
   const bool PrepVz           = ( TVar & _VELZ ) ? true : false;
   const bool PrepPres         = ( TVar & _PRES ) ? true : false;
   const bool PrepTemp         = ( TVar & _TEMP ) ? true : false;

#  elif ( MODEL == MHD   )
#  warning : WAIT MHD !!

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
      BC_User_Ptr( BVal, x, y, z, Time, lv, NULL );

      for (int v=0; v<NVar_Flu; v++)   Array3D[v][k][j][i] = BVal[ TFluVarIdxList[v] ];


//    derived variables
      v2 = NVar_Flu;

#     if   ( MODEL == HYDRO )
      if ( PrepVx   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMX] / BVal[DENS];
      if ( PrepVy   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMY] / BVal[DENS];
      if ( PrepVz   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMZ] / BVal[DENS];
      if ( PrepPres )   Array3D[ v2 ++ ][k][j][i] = CPU_GetPressure( BVal[DENS], BVal[MOMX], BVal[MOMY], BVal[MOMZ], BVal[ENGY],
                                                                     Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      if ( PrepTemp )   Array3D[ v2 ++ ][k][j][i] = CPU_GetTemperature( BVal[DENS], BVal[MOMX], BVal[MOMY], BVal[MOMZ], BVal[ENGY],
                                                                        Gamma_m1, CheckMinPres_Yes, MIN_PRES );

#     elif ( MODEL == MHD   )
#     warning : WAIT MHD !!

#     elif ( MODEL == ELBDM )
//    no derived variables yet

#     else
#     error : unsupported MODEL !!
#     endif
   } // k,j,i

} // FUNCTION : Flu_BoundaryCondition_User
