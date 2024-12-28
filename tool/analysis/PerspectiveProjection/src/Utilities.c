#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  checkNAN
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
bool checkNAN( real input, const char FunctionName[], const int line )
{

   if ( input != input )
   {
      printf( "nan is at %s:%d\n", FunctionName, line );
      return true;
   }

   if ( input > __DBL_MAX__  ||  input < __DBL_MIN__ )
   {
      printf( "input is -inf/inf %s:%d\n", FunctionName, line );
      return true;
   }

   return false;
} // FUNCTION : checkNAN



//-------------------------------------------------------------------------------------------------------
// Function    :  checkMinus
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
bool checkMinus( real input, const char FunctionName[], const int line )
{
   if ( input >= 0.0 )   return false;

   printf( "minus is at %s:%d\n", FunctionName, line );
   return true;
} // FUNCTION : checkMinus



//-------------------------------------------------------------------------------------------------------
// Function    :  checkInt32Overflow
// Description :
// Note        :
// Parameter   :
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void checkInt32Overflow( int32_t a, int32_t b, int operation, int line )
{
   bool overflow = false;

   if ( operation == '+' )
   {
      if ( ( b > 0 ) && ( a > INT32_MAX-b ) )   overflow = true;
      if ( ( b < 0 ) && ( a < INT32_MIN-b ) )   overflow = true;
   }
   else if ( operation == '*' )
   {
      if (  b != 0 && a > INT32_MAX / b  )      overflow = true;
      if (  b != 0 && a < INT32_MIN / b  )      overflow = true;
   }
   else
   {
      ERROR_EXIT( 0, "ERROR : something wrong !!\n" );
   }

   if ( overflow )   ERROR_EXIT( 0, "ERROR : Integer overflow !! a=%d, b=%d, line=%d\n", a, b, line );
} // FUNCTION : checkInt32Overflow



//-------------------------------------------------------------------------------------------------------
// Function    :  checkmemoryContiguous
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
bool checkmemoryContiguous( real *Ptr, int sizeOf, int Length )
{
   long d = (Ptr+1)-(Ptr);

   for (int i=2; i<Length; i++)
      if ( (Ptr+i) - (Ptr+i-1) != d )
         return false;
         // printf("%ld\n", (Ptr+i) - (Ptr+i-1));

   return true;
} // FUNCTION : checkmemoryContiguous



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputBinary
// Description :
// Note        :
// Parameter   :
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void OutputBinary( void *Ptr, int size, int count, char Name [] )
{
   FILE *pFile;

   pFile = fopen( Name, "wb" );

   fwrite( Ptr, size, count, pFile );

   fclose( pFile );
} // FUNCTION : OutputBinary
