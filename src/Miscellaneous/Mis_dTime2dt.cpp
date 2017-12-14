#include "GAMER.h"



// conversion method
// ======================================================================
#define RK4
//#define EXPAND_THIRD


// inline function evaluating the first derivative of dt over dTime
// ======================================================================
#if ( defined RK4  &&  defined COMOVING )
inline double RK4_FirstDerivative( const double Time_In )
{
   return pow(  OMEGA_M0*pow( Time_In, 3.0 ) + (1.0-OMEGA_M0)*pow( Time_In, 6.0 ),  -0.5  );
}
#endif // #if ( defined RK4  &&  defined COMOVING )




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_dTime2dt
// Description :  Convert the input physical time interval (dTime) to the evolution time-step (dt)
//
// Note        :  1. Physical coordinates : dTime == dt
//                   Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
//                2. Two methods are supported here:
//                   (1) Fourth-order Runge-Kutta method (default)
//                   (2) Expand dt to the third power
//
// Parameter   :  Time_In  : Current physical time
//                dTime_In : Current physical time interval
//
// Return      :  Evolution time-step (dt_Out)
//-------------------------------------------------------------------------------------------------------
double Mis_dTime2dt( const double Time_In, const double dTime_In )
{

   double dt_Out;    // variable to be returned


#  ifdef COMOVING

// Method 1 : RK4
// =============================================================================================================
#  ifdef RK4
   const double k1 = dTime_In*RK4_FirstDerivative( Time_In + 0.0*dTime_In );
   const double k2 = dTime_In*RK4_FirstDerivative( Time_In + 0.5*dTime_In );
   const double k3 = k2;
   const double k4 = dTime_In*RK4_FirstDerivative( Time_In + 1.0*dTime_In );

   dt_Out = ( k1 + 2.0*k2 + 2.0*k3 + k4 ) / 6.0;


// Method 2 : expand dt to the third power
// =============================================================================================================
#  elif ( defined EXPAND_THIRD )
   double dt_dTime1, dt_dTime2, dt_dTime3;   // first/second/third derivatives

   dt_dTime1 =       pow(      OMEGA_M0*pow(Time_In, 3.0) +      (1.0-OMEGA_M0)*pow(Time_In, 6.0),  -0.5  );
   dt_dTime2 = -0.50*pow(      OMEGA_M0*pow(Time_In, 3.0) +      (1.0-OMEGA_M0)*pow(Time_In, 6.0),  -1.5  ) *
                        (  3.0*OMEGA_M0*pow(Time_In, 2.0) +  6.0*(1.0-OMEGA_M0)*pow(Time_In, 5.0)         );
   dt_dTime3 =  0.75*pow(      OMEGA_M0*pow(Time_In, 3.0) +      (1.0-OMEGA_M0)*pow(Time_In, 6.0),  -2.5  ) *
                     pow(  3.0*OMEGA_M0*pow(Time_In, 2.0) +  6.0*(1.0-OMEGA_M0)*pow(Time_In, 5.0),   2.0  ) +
               -0.50*pow(      OMEGA_M0*pow(Time_In, 3.0) +      (1.0-OMEGA_M0)*pow(Time_In, 6.0),  -1.5  ) *
                        (  6.0*OMEGA_M0*    Time_In       + 30.0*(1.0-OMEGA_M0)*pow(Time_In, 4.0)         );

// expand dt to the 3rd power
   dt_Out = dt_dTime1*dTime_In + 0.5*dt_dTime2*dTime_In*dTime_In + 1.0/6.0*dt_dTime3*dTime_In*dTime_In*dTime_In;

#  else
#  error : ERROR : NO METHOD OF CONVERTING DTIME TO DT IS FOUND !!
#  endif // conversion method


#  else // #ifndef COMOVING

// nothing to convert in the physical frame
// =============================================================================================================
   dt_Out = dTime_In;

#  endif // #ifdef COMOVING ... else ...


   return dt_Out;

} // FUNCTION : Mis_dTime2dt
