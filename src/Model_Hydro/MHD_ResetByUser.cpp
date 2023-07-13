#include "GAMER.h"

#ifdef MHD


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static double MHD_ResetByUser_VecPot_Template( const double x, const double y, const double z, const double Time,
                                               const double dt, const int lv, const char Component, double AuxArray[] );
static double MHD_ResetByUser_BField_Template( const double x, const double y, const double z, const double Time,
                                               const double dt, const int lv, const char Component, double AuxArray[], const double B_in,
                                               const bool UseVecPot, const real *Ax, const real *Ay, const real *Az,
                                               const int i, const int j, const int k );

// this function pointer must be set by a test problem initializer
double (*MHD_ResetByUser_VecPot_Ptr)( const double x, const double y, const double z, const double Time,
                                      const double dt, const int lv, const char Component, double AuxArray[] ) = NULL;
double (*MHD_ResetByUser_BField_Ptr)( const double x, const double y, const double z, const double Time,
                                      const double dt, const int lv, const char Component, double AuxArray[], const double B_in,
                                      const bool UseVecPot, const real *Ax, const real *Ay, const real *Az,
                                      const int i, const int j, const int k ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ResetByUser_VecPot_Template
// Description :  Function template to reset the magnetic field by a vector potential
//
// Note        :  1. Invoked by Flu_ResetByUser_API_Default() and Model_Init_ByFunction_AssignData() using the
//                   function pointer "MHD_ResetByUser_VecPot_Ptr", which must be set by a test problem initializer
//                   --> Model_Init_ByFunction_AssignData(): constructing an initial condition
//                       Flu_ResetByUser_API_Default()     : during evolution
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                dt        : Time interval to advance solution
//                lv        : Target refinement level
//                Component : Component of the output magnetic vector potential
//                            --> Supported components: 'x', 'y', 'z'
//                AuxArray  : Auxiliary array
//                            --> Useless since it is currently fixed to NULL
//
// Return      :  "Component"-component of the magnetic vector potential at (x, y, z, Time)
//-------------------------------------------------------------------------------------------------------
double MHD_ResetByUser_VecPot_Template( const double x, const double y, const double z, const double Time,
                                        const double dt, const int lv, const char Component, double AuxArray[] )
{

   const double amp  = 1.0e2;   // amplitude
   const double r0   = 1.0e-2;  // scale radius
   const double tmin = 1.0e-3;  // starting time
   const double tmax = 2.0e-3;  // ending time

   const double dx   = x - amr->BoxCenter[0];
   const double dy   = y - amr->BoxCenter[1];
   const double dz   = z - amr->BoxCenter[2];
   const double r    = sqrt( SQR(dx) + SQR(dy) + SQR(dz) );
   const double Az   = amp*exp( -r/r0 );

   double A;

   switch ( Component )
   {
      case 'x' : A = 0.0;  break;
      case 'y' : A = 0.0;  break;
      case 'z' : A = ( Time >= tmin  &&  Time <= tmax ) ? Az : 0.0;  break;

      default  : Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
   }

   return A;

} // FUNCTION : MHD_ResetByUser_VecPot_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ResetByUser_BField_Template
// Description :  Function template to reset the magnetic field
//
// Note        :  1. Invoked by Flu_ResetByUser_API_Default() and Model_Init_ByFunction_AssignData() using the
//                   function pointer "MHD_ResetByUser_BField_Ptr", which must be set by a test problem initializer
//                   --> Model_Init_ByFunction_AssignData(): constructing an initial condition
//                       Flu_ResetByUser_API_Default()     : during evolution
//                2. Enabled by the runtime options "OPT__RESET_FLUID" or "OPT__RESET_FLUID_INIT"
//                3. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                4. Total energy density will be updated automatically
//                5. To ensure divergence-free to the machine precision, one should either use vector potential
//                   by setting MHD_ResetByUser_VecPot_Ptr or compute the *area-average* B field in this routine
//                   --> For example, Bx should be the value averaged over
//                       (x,y,z) = (x, [y-0.5*amr->dh[lv] ... y+0.5*amr->dh[lv]], [z-0.5*amr->dh[lv] ... z+0.5*amr->dh[lv]])
//                   --> When using vector potential, the divergence-free condition can still break on the coarse-fine interfaces
//                       --> Will be improved in the future
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                dt        : Time interval to advance solution
//                lv        : Target refinement level
//                Component : Component of the output magnetic field
//                            --> Supported components: 'x', 'y', 'z'
//                AuxArray  : Auxiliary array
//                            --> Useless since it is currently fixed to NULL
//                B_in      : Original magnetic field before resetting
//                UseVecPot : Use vector potential to compute the B field
//                Ax/y/z    : Vector potential arrays used by UseVecPot
//                i/j/k     : Array indices used by UseVecPot
//
// Return      :  "Component"-component of the magnetic field at (x, y, z, Time)
//-------------------------------------------------------------------------------------------------------
double MHD_ResetByUser_BField_Template( const double x, const double y, const double z, const double Time,
                                        const double dt, const int lv, const char Component, double AuxArray[], const double B_in,
                                        const bool UseVecPot, const real *Ax, const real *Ay, const real *Az,
                                        const int i, const int j, const int k )
{

   const double amp  = 1.0e2;   // amplitude
   const double r0   = 1.0e-2;  // scale radius
   const double tmin = 1.0e-3;  // starting time
   const double tmax = 2.0e-3;  // ending time

   const double dh   = amr->dh[lv];
   const double dx   = x - amr->BoxCenter[0];
   const double dy   = y - amr->BoxCenter[1];
   const double dz   = z - amr->BoxCenter[2];
   const double r    = sqrt( SQR(dx) + SQR(dy) + SQR(dz) );
   double B_out;

   if ( UseVecPot )
   {
      switch ( Component )
      {
         case 'x' :  B_out = B_in + dt*MHD_ResetByUser_A2Bx( Ax, Ay, Az, i, j, k, dh );    break;
         case 'y' :  B_out = B_in + dt*MHD_ResetByUser_A2By( Ax, Ay, Az, i, j, k, dh );    break;
         case 'z' :  B_out = B_in + dt*MHD_ResetByUser_A2Bz( Ax, Ay, Az, i, j, k, dh );    break;
         default  :  Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
      }
   } // if ( UseVecPot )

   else
   {
      if ( Time < tmin  ||  Time > tmax )
         B_out = B_in;

      else
         switch ( Component )
         {
            case 'x' :  B_out = B_in - dt*amp*dy/(r0*r)*exp( -r/r0 );   break;
            case 'y' :  B_out = B_in + dt*amp*dx/(r0*r)*exp( -r/r0 );   break;
            case 'z' :  B_out = B_in;                                   break;
            default  :  Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
         }
   } // if ( UseVecPot ) ... else ...

   return B_out;

} // FUNCTION : MHD_ResetByUser_BField_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ResetByUser_A2Bx/y/z
// Description :  Convert vector potential to B field
//
// Note        :  1. Called by MHD_ResetByUser_BField_Ptr
//
// Parameter   :  Ax/y/z : Vector potential arrays
//                i/j/k  : Array indices
//                dh     : Cell width
//
// Return      :  Bx/By/Bz
//-------------------------------------------------------------------------------------------------------
real MHD_ResetByUser_A2Bx( const real *Ax, const real *Ay, const real *Az,
                           const int i, const int j, const int k, const double dh )
{
   const int  idx  = IDX321( i, j,   k,   PS1+1, PS1+1 );
   const int  idxj = IDX321( i, j+1, k,   PS1+1, PS1+1 );
   const int  idxk = IDX321( i, j,   k+1, PS1+1, PS1+1 );
   const real Bx   = ( Az[idxj] - Az[idx] - Ay[idxk] + Ay[idx] ) / dh;

   return Bx;
} // FUNCTION : MHD_ResetByUser_A2Bx



real MHD_ResetByUser_A2By( const real *Ax, const real *Ay, const real *Az,
                           const int i, const int j, const int k, const double dh )
{
   const int  idx  = IDX321( i,   j, k,   PS1+1, PS1+1 );
   const int  idxi = IDX321( i+1, j, k,   PS1+1, PS1+1 );
   const int  idxk = IDX321( i,   j, k+1, PS1+1, PS1+1 );
   const real By   = ( Ax[idxk] - Ax[idx] - Az[idxi] + Az[idx] ) / dh;

   return By;
} // FUNCTION : MHD_ResetByUser_A2By



real MHD_ResetByUser_A2Bz( const real *Ax, const real *Ay, const real *Az,
                           const int i, const int j, const int k, const double dh )
{
   const int  idx  = IDX321( i,   j,   k, PS1+1, PS1+1 );
   const int  idxi = IDX321( i+1, j,   k, PS1+1, PS1+1 );
   const int  idxj = IDX321( i,   j+1, k, PS1+1, PS1+1 );
   const real Bz   = ( Ay[idxi] - Ay[idx] - Ax[idxj] + Ax[idx] ) / dh;

   return Bz;
} // FUNCTION : MHD_ResetByUser_A2Bz



#endif // #ifdef MHD
