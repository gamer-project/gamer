#include "GAMER.h"

#ifdef MHD

extern double Blast_ResetB_amp;
extern double Blast_ResetB_r0;
extern double Blast_ResetB_tmin;
extern double Blast_ResetB_tmax;




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ResetByUser_VecPot_BlastWave
// Description :  Function to reset the magnetic field by a vector potential
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
double MHD_ResetByUser_VecPot_BlastWave( const double x, const double y, const double z, const double Time,
                                         const double dt, const int lv, const char Component, double AuxArray[] )
{

   const double amp  = Blast_ResetB_amp;
   const double r0   = Blast_ResetB_r0;
   const double tmin = Blast_ResetB_tmin;
   const double tmax = Blast_ResetB_tmax;
   const double dx   = x - amr->BoxCenter[0];
   const double dy   = y - amr->BoxCenter[1];
   const double dz   = z - amr->BoxCenter[2];
   const double r    = sqrt( SQR(dx) + SQR(dy) + SQR(dz) );
   const double Az   = amp*exp( -SQR(r/r0) );

   double A;

   switch ( Component )
   {
      case 'x' : A = 0.0;  break;
      case 'y' : A = 0.0;  break;
      case 'z' : A = ( Time >= tmin  &&  Time <= tmax ) ? Az : 0.0;  break;

      default  : Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
   }

   return A;

} // FUNCTION : MHD_ResetByUser_VecPot_BlastWave



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ResetByUser_BField_BlastWave
// Description :  Function to reset the magnetic field
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
double MHD_ResetByUser_BField_BlastWave( const double x, const double y, const double z, const double Time,
                                         const double dt, const int lv, const char Component, double AuxArray[], const double B_in,
                                         const bool UseVecPot, const real *Ax, const real *Ay, const real *Az,
                                         const int i, const int j, const int k )
{

   const double amp  = Blast_ResetB_amp;
   const double r0   = Blast_ResetB_r0;
   const double tmin = Blast_ResetB_tmin;
   const double tmax = Blast_ResetB_tmax;
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
            case 'x' :  B_out = B_in - dt*2.0*amp*dy/SQR(r0)*exp( -SQR(r/r0) );  break;
            case 'y' :  B_out = B_in + dt*2.0*amp*dx/SQR(r0)*exp( -SQR(r/r0) );  break;
            case 'z' :  B_out = B_in;                                            break;
            default  :  Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
         }
   } // if ( UseVecPot ) ... else ...

   return B_out;

} // FUNCTION : MHD_ResetByUser_BField_BlastWave



#endif // #ifdef MHD
