#include "GAMER.h"
#include <typeinfo>




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_CompareRealValue
// Description :  Check if two input floating-point values are within the tolerable round-off errors
//
// Note        :  Work for both single and double precisions.
//                Explicit template instantiation is put in the end of this file
//
// Parameter   :  Input1   : The first  input time
//                Input2   : The second input time
//                comment  : You can put the location where this function is invoked in this string
//                Verbose  : Output the warning message if the check fails
//
// Return      :  true (false) --> Input 1 and 2 are (not) within the tolerable round-off errors
//-------------------------------------------------------------------------------------------------------
template <typename T>
bool Mis_CompareRealValue( const T Input1, const T Input2, const char *comment, const bool Verbose )
{

   const double TolErr = ( typeid(T) == typeid(double) ) ? 1.0e-12 : 1.0e-5;

   double RelErr;

   if      ( Input1 == Input2 )                 RelErr = 0.0;
   else if ( Input1 == 0.0 || Input2 == 0.0 )   RelErr = fabs(  Input1 - Input2 );
   else                                         RelErr = fabs( (Input1 - Input2)/Input1 );

   if ( RelErr > TolErr )
   {
      if ( Verbose )
      {
         Aux_Message( stderr, "WARNING : \"%s\" : <%s> FAILED at rank %d !!\n",
                      comment, __FUNCTION__, MPI_Rank );
         Aux_Message( stderr, "          Input1 = %20.14e vs. Input2 = %20.14e --> RelErr = %20.14e !!\n",
                      Input1, Input2, RelErr );
      }

      return false;
   }

   else
      return true;

} // FUNCTION : Mis_CompareRealValue



// explicit template instantiation
template bool Mis_CompareRealValue <float > ( const float  Input1, const float  Input2, const char *comment, const bool Verbose );
template bool Mis_CompareRealValue <double> ( const double Input1, const double Input2, const char *comment, const bool Verbose );
