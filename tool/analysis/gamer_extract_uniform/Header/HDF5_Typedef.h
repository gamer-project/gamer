#ifndef __HDF5_TYPEDEF_H__
#define __HDF5_TYPEDEF_H__


/*===========================================================================
Data structures defined here are mainly used for creating the compound
datatypes in the HDF5 format
===========================================================================*/

#include "hdf5.h"

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif




//-------------------------------------------------------------------------------------------------------
// Structure   :  Info_t
// Description :  Data structure for storing the simulation information
//
// Note        :  1. Used by SetHDF5Info()
//-------------------------------------------------------------------------------------------------------
struct Info_t
{

   int    DumpID;
   int    GridDimension[3];
   int    WithUnit;
   int    OutputMode;

   double Time;
   double CellWidth;
   double SubdomainSize[3];
   double SubdomainLeftEdge[3];
   double Unit_L;
   double Unit_M;
   double Unit_T;
   double Unit_V;
   double Unit_D;
   double Unit_E;
   double Unit_P;

#  if ( MODEL == ELBDM )
   double ELBDM_Mass;
#  endif

}; // struct Info_t



//-------------------------------------------------------------------------------------------------------
// Function    :  SyncHDF5File
// Description :  Force NFS to synchronize (dump data to the disk before return)
//
// Note        :  1. Synchronization is accomplished by first opening the target file with the "appending" mode
//                   and then close it immediately
//                2. It's an inline function
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
inline void SyncHDF5File( const char *FileName )
{

   FILE *File      = fopen( FileName, "ab" );
   int CloseStatus = fclose( File );

   if ( CloseStatus != 0 )
      Aux_Message( stderr, "WARNING : failed to close the file \"%s\" for synchronization !!\n", FileName );

} // FUNCTION : SyncHDF5File



#endif // #ifndef __HDF5_TYPEDEF_H__
