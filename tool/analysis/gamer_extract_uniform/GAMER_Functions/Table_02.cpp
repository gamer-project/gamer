#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_02
// Description :  Return w0, or w1 according to the inputs "LocalID" and "dim"
//
// Note        :  This table is particularly useful when the required values depended on the direction
//                of the local ID of patch within a patch group
//-------------------------------------------------------------------------------------------------------
int TABLE_02( const int LocalID, const char dim, const int w0, const int w1 )
{

   switch ( dim )
   {
      case 'x':
         switch ( LocalID )
         {
            case 0: case 2: case 3: case 5:
               return w0;

            case 1: case 4: case 6: case 7:
               return w1;

            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %c, %s = %d\" !!\n",
                        "dim", dim, "LocalID", LocalID );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }

      case 'y':
         switch ( LocalID )
         {
            case 0: case 1: case 3: case 6:
               return w0;

            case 2: case 4: case 5: case 7:
               return w1;

            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %c, %s = %d\" !!\n",
                        "dim", dim, "LocalID", LocalID );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }

      case 'z':
         switch ( LocalID )
         {
            case 0: case 1: case 2: case 4:
               return w0;

            case 3: case 5: case 6: case 7:
               return w1;

            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %c, %s = %d\" !!\n",
                        "dim", dim, "LocalID", LocalID );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
         }

      default:
         fprintf( stderr, "ERROR : \"incorrect parameter %s = %c\" !!\n", "dim", dim );
         fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         MPI_Exit();
         exit( -1 );
   }

}





