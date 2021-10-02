#ifndef __HDF5_TYPEDEF_H__
#define __HDF5_TYPEDEF_H__


/*===========================================================================
Data structures defined here are mainly used for creating the compound
datatypes in the HDF5 format
===========================================================================*/

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif




//-------------------------------------------------------------------------------------------------------
// Structure   :  PatchInfo_t
// Description :  Data structure for outputting the patch information
//
// Note        :  1. Father, Son and Sibling store the "global identification (GID)" instead of the "local
//                   patch identification (PID)"
//                   --> Patches at the same level at all MPI ranks will have different GID but can have
//                       the same PID
//                2. LB_Idx will be stored in all parallelization modes (SERIAL, NO_LOAD_BALANCE, LOAD_BALANCE)
//-------------------------------------------------------------------------------------------------------
struct PatchInfo_t
{

   int    Lv;
   int    GID;
   int    Father;
   int    Son;
   int    Sibling[26];
   int    Corner[3];

   ulong  PaddedCr1D;
   long   LBIdx;

   double EdgeL[3];
   double EdgeR[3];

}; // struct PatchInfo_t



//-------------------------------------------------------------------------------------------------------
// Structure   :  PatchData_t
// Description :  Data structure for outputting the patch data
//
// Note        :  1. "Pote" array MUST be the LAST member
//-------------------------------------------------------------------------------------------------------
struct PatchData_t
{

#  if   ( MODEL == HYDRO )
   real Dens[PS1][PS1][PS1];
   real MomX[PS1][PS1][PS1];
   real MomY[PS1][PS1][PS1];
   real MomZ[PS1][PS1][PS1];
   real Engy[PS1][PS1][PS1];

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real Dens[PS1][PS1][PS1];
   real Real[PS1][PS1][PS1];
   real Imag[PS1][PS1][PS1];

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

// Pote MUST be the LAST member of the PatchData_t structure (since it's output only if OPT__OUTPUT_POT == true)
   real Pote[PS1][PS1][PS1];

}; // struct PatchData_t



#endif // #ifndef __HDF5_TYPEDEF_H__
