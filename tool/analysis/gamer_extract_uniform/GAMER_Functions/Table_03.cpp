#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_03
// Description :  Return the information of patch ID for the fucntions "Prepare_InputGPUArray",
//                "SiblingSearch_Base", "Aux_Check_FluxAllocate", and "Flu_AllocateFluxArray_Buffer"
//
// Parameter   :  SibID : sibling index (0~25)
//                Count : patch counter (0~3)
//-------------------------------------------------------------------------------------------------------
int TABLE_03( const int SibID, const int Count )
{

   switch ( SibID )
   {
      case 0:
      {
         switch( Count )
         {
            case 0:  return 1;
            case 1:  return 6;
            case 2:  return 4;
            case 3:  return 7;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 1:
      {
         switch ( Count )
         {
            case 0:  return 0;
            case 1:  return 3;
            case 2:  return 2;
            case 3:  return 5;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 2:
      {
         switch ( Count )
         {
            case 0:  return 2;
            case 1:  return 4;
            case 2:  return 5;
            case 3:  return 7;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 3:
      {
         switch ( Count )
         {
            case 0:  return 0;
            case 1:  return 1;
            case 2:  return 3;
            case 3:  return 6;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 4:
      {
         switch ( Count )
         {
            case 0:  return 3;
            case 1:  return 5;
            case 2:  return 6;
            case 3:  return 7;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 5:
      {
         switch ( Count )
         {
            case 0:  return 0;
            case 1:  return 2;
            case 2:  return 1;
            case 3:  return 4;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 6:
      {
         switch ( Count )
         {
            case 0:  return 4;
            case 1:  return 7;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 7:
      {
         switch ( Count )
         {
            case 0:  return 2;
            case 1:  return 5;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 8:
      {
         switch ( Count )
         {
            case 0:  return 1;
            case 1:  return 6;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 9:
      {
         switch ( Count )
         {
            case 0:  return 0;
            case 1:  return 3;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 10:
      {
         switch ( Count )
         {
            case 0:  return 5;
            case 1:  return 7;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 11:
      {
         switch ( Count )
         {
            case 0:  return 3;
            case 1:  return 6;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 12:
      {
         switch ( Count )
         {
            case 0:  return 2;
            case 1:  return 4;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 13:
      {
         switch ( Count )
         {
            case 0:  return 0;
            case 1:  return 1;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 14:
      {
         switch ( Count )
         {
            case 0:  return 6;
            case 1:  return 7;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 15:
      {
         switch ( Count )
         {
            case 0:  return 1;
            case 1:  return 4;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 16:
      {
         switch ( Count )
         {
            case 0:  return 3;
            case 1:  return 5;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 17:
      {
         switch ( Count )
         {
            case 0:  return 0;
            case 1:  return 2;
            default:
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d, %s = %d\" !!\n",
                        "SibID", SibID, "Count", Count );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }
      }

      case 18:    return 7;
      case 19:    return 5;
      case 20:    return 6;
      case 21:    return 3;
      case 22:    return 4;
      case 23:    return 2;
      case 24:    return 1;
      case 25:    return 0;

      default:
      {
         fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "SibID", SibID );
         fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         MPI_Exit();
         exit( -1 );
      }

   } // switch ( SibID )

}


