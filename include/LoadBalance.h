#ifndef __LOADBALANCE_H__
#define __LOADBALANCE_H__



#include "Macro.h"




//-------------------------------------------------------------------------------------------------------
// Structure   :  LB_t
// Description :  Data structure for the load-balance parallelization
//
// Data_Member :  MPI_NRank               : Number of MPI ranks ( == global variable "MPI_NRank" )
//                OverlapMPI_FluSyncN     : Number of patches with LocalID==0 which will NOT be overlapped with
//                                          the MPI communication (fluid solver)
//                OverlapMPI_FluSyncPID0  : Patch indices with LocalID==0 which will NOT be overlapped with
//                                          the MPI communication (fluid solver)
//                OverlapMPI_FluAsyncN    : Number of patches with LocalID==0 which will be overlapped with
//                                          the MPI communication (fluid solver)
//                OverlapMPI_FluAsyncPID0 : Patch indices with LocalID==0 which will be overlapped with
//                                          the MPI communication (fluid solver)
//                OverlapMPI_PotSyncN     : Number of patches with LocalID==0 which will NOT be overlapped with
//                                          the MPI communication (gravity solver)
//                OverlapMPI_PotSyncPID0  : Patch indices with LocalID==0 which will NOT be overlapped with
//                                          the MPI communication (gravity solver)
//                OverlapMPI_PotAsyncN    : Number of patches with LocalID==0 which will be overlapped with
//                                          the MPI communication (gravity solver)
//                OverlapMPI_PotAsyncPID0 : Patch indices with LocalID==0 which will be overlapped with
//                                          the MPI communication (gravity solver)
//                WLI                     : Weighted load-imbalance factor of patches at all levels
//                WLI_Max                 : WLI threshold for redistributing patches at all levels
//                Par_Weight              : Load-balance weighting of one particle over one cell
//                                          --> Weighting of each patch is estimated as "PATCH_SIZE^3 + NParThisPatch*Par_Weight"
//                CutPoint                : Cut points in the space filling curve
//                IdxList_Real            : Sorted LB_Idx list of all real patches
//                IdxList_Real_IdxTable   : Index table for LB_IdxList_Real
//                PaddedCr1DList          : Sorted PaddedCr1D list of all patches (real + buffer)
//                PaddedCr1DList_IdxTable : Index table for LB_PaddedC1DrList
//
//                SendH_NList             : Number of patches for sending hydrodynamic data
//                SendH_IDList            : Patch indices for sending hydrodynamic data
//                SendH_SibList           : Sibling indices for sending hydrodynamic data
//                SendH_SibDiffList       : Sibling indices for sending hydrodynamic data after grid refinement
//                SendH_LBIdxList         : Load-balance indices for sending hydrodynamic data
//                RecvH_NList             : Number of patches for receiving hydrodynamic data
//                RecvH_IDList            : Patch indices for receiving hydrodynamic data
//                RecvH_IDList_IdxTable   : Index table for LB_RecvH_IDList
//                RecvH_SibList           : Sibling indices for receiving hydrodynamic data
//                RecvH_SibDiffList       : Sibling indices for receiving hydrodynamic data after grid refinement
//                RecvH_LBIdxList         : Load-balance indices for receiving hydrodynamic data
//                RecvH_PCr1D             : Padded 1D corner for receiving hydrodynamic data
//                RecvH_PCr1D_IdxTable    : Index table for RecvH_PCr1D
//
//                SendX_NList             : Number of patches for sending hydrodynamic data after fix-up
//                SendX_NResList          : Number of patches for sending hydrodynamic data after restriction
//                SendX_IDList            : Patch indices for sending hydrodynamic data after fix-up
//                SendX_SibList           : Sibling indices for sending hydrodynamic data after fix-up
//                RecvX_NList             : Number of patches for receiving hydrodynamic data after fix-up
//                RecvX_NResList          : Number of patches for receiving hydrodynamic data after restriction
//                RecvX_IDList            : Patch indices for receiving hydrodynamic data after fix-up
//                RecvX_SibList           : Sibling indices for receiving hydrodynamic data after fix-up
//
//                SendR_NList             : Number of patches for sending restricted hydro data
//                SendR_IDList            : Patch indices for sending restricted hydro data
//                SendR_IDList_IdxTable   : Index table for LB_SendR_IDList
//                RecvR_NList             : Number of patches for receiving restricted hydro data
//                RecvR_IDList            : Patch indices for receiving restricted hydro data
//
//                SendF_NList             : Number of flux arrays for sending
//                SendF_IDList            : Patch indices for sending fluxes
//                SendF_SibList           : Sibling directions for sending fluxes
//                RecvF_NList             : Number of flux arrays for receiving
//                RecvF_IDList            : Patch indices for receiving fluxes
//                RecvF_IDList_IdxTable   : Index table for LB_RecvF_IDList_IdxTable
//                RecvF_SibList           : Sibling directions for receiving fluxes
//
//                SendG_NList             : Number of patches for sending potential data
//                SendG_IDList            : Patch indices for sending potential data
//                SendG_SibList           : Sibling indices for sending potential data
//                SendG_SibDiffList       : Sibling indices for sending potential data after grid refinement
//                SendG_LBIdxList         : Load-balance indices for sending potential data
//                RecvG_NList             : Number of patches for receiving potential data
//                RecvG_IDList            : Patch indices for receiving potential data
//                RecvG_IDList_IdxTable   : Index table for LB_RecvH_IDList
//                RecvG_SibList           : Sibling indices for receiving potential data
//                RecvG_SibDiffList       : Sibling indices for receiving potential data after grid refinement
//                RecvG_LBIdxList         : Load-balance indices for receiving potential data
//                RecvG_PCr1D             : Padded 1D corner for receiving potential data
//                RecvG_PCr1D_IdxTable    : Index table for RecvG_PCr1D
//
// Method      :  LB_t  : Constructor
//               ~LB_t  : Destructor
//                reset : Reset all pointers and counters for re-constructing load-balance parallelization
//-------------------------------------------------------------------------------------------------------
struct LB_t
{

// data members
// ===================================================================================
   int    MPI_NRank;
   int    OverlapMPI_FluSyncN    [NLEVEL];
   int   *OverlapMPI_FluSyncPID0 [NLEVEL];
   int    OverlapMPI_FluAsyncN   [NLEVEL];
   int   *OverlapMPI_FluAsyncPID0[NLEVEL];
#  ifdef GRAVITY
   int    OverlapMPI_PotSyncN    [NLEVEL];
   int   *OverlapMPI_PotSyncPID0 [NLEVEL];
   int    OverlapMPI_PotAsyncN   [NLEVEL];
   int   *OverlapMPI_PotAsyncPID0[NLEVEL];
#  endif
   double WLI;
   double WLI_Max;
#  ifdef PARTICLE
   double Par_Weight;
#  endif
   long  *CutPoint               [NLEVEL];
   long  *IdxList_Real           [NLEVEL];
   int   *IdxList_Real_IdxTable  [NLEVEL];
   ulong *PaddedCr1DList         [NLEVEL];
   int   *PaddedCr1DList_IdxTable[NLEVEL];

   int   *SendH_NList            [NLEVEL];
   int  **SendH_IDList           [NLEVEL];
   int  **SendH_SibList          [NLEVEL];
   int  **SendH_SibDiffList      [NLEVEL];
   long **SendH_LBIdxList        [NLEVEL];
   int   *RecvH_NList            [NLEVEL];
   int  **RecvH_IDList           [NLEVEL];
   int  **RecvH_IDList_IdxTable  [NLEVEL];
   int  **RecvH_SibList          [NLEVEL];
   int  **RecvH_SibDiffList      [NLEVEL];
   long **RecvH_LBIdxList        [NLEVEL];
  ulong **RecvH_PCr1D            [NLEVEL];
   int  **RecvH_PCr1D_IdxTable   [NLEVEL];

   int  *SendX_NList             [NLEVEL];
   int  *SendX_NResList          [NLEVEL];
   int **SendX_IDList            [NLEVEL];
   int **SendX_SibList           [NLEVEL];
   int  *RecvX_NList             [NLEVEL];
   int  *RecvX_NResList          [NLEVEL];
   int **RecvX_IDList            [NLEVEL];
   int **RecvX_SibList           [NLEVEL];

   int  *SendR_NList             [NLEVEL];
   int **SendR_IDList            [NLEVEL];
   int **SendR_IDList_IdxTable   [NLEVEL];
   int  *RecvR_NList             [NLEVEL];
   int **RecvR_IDList            [NLEVEL];

   int  *SendF_NList             [NLEVEL];
   int **SendF_IDList            [NLEVEL];
   int **SendF_SibList           [NLEVEL];
   int  *RecvF_NList             [NLEVEL];
   int **RecvF_IDList            [NLEVEL];
   int **RecvF_IDList_IdxTable   [NLEVEL];
   int **RecvF_SibList           [NLEVEL];

#  ifdef GRAVITY
   int   *SendG_NList            [NLEVEL];
   int  **SendG_IDList           [NLEVEL];
   int  **SendG_SibList          [NLEVEL];
   int  **SendG_SibDiffList      [NLEVEL];
   long **SendG_LBIdxList        [NLEVEL];
   int   *RecvG_NList            [NLEVEL];
   int  **RecvG_IDList           [NLEVEL];
   int  **RecvG_IDList_IdxTable  [NLEVEL];
   int  **RecvG_SibList          [NLEVEL];
   int  **RecvG_SibDiffList      [NLEVEL];
   long **RecvG_LBIdxList        [NLEVEL];
  ulong **RecvG_PCr1D            [NLEVEL];
   int  **RecvG_PCr1D_IdxTable   [NLEVEL];
#  endif



   //===================================================================================
   // Constructor :  LB_t
   // Description :  Constructor of the structure "LB_t"
   //
   // Note        :  1. Allocate memory for pointers whose sizes depend on the number of MPI ranks
   //                2. Initialize pointers as NULL and counters as zero.
   //                3. "IdxList_Real, IdxList_Real_IdxTable, PaddedCr1DList, and
   //                   PaddedCr1DList_IdxTable", whose sizes can not be determined during
   //                   initialization, are NOT allocated with memory
   //
   // Parameter   :  NRank             : Number of MPI ranks
   //                Input__WLI_Max    : WLI_Max loaded from the input parameter file
   //                Input__Par_Weight : Par_Weight loaded from the input parameter file
   //===================================================================================
   LB_t( const int NRank, const double Input__WLI_Max, const double Input__Par_Weight )
   {

      MPI_NRank  = NRank;
      WLI        = NULL_REAL;
      WLI_Max    = Input__WLI_Max;
#     ifdef PARTICLE
      Par_Weight = Input__Par_Weight;
#     endif

      for (int lv=0; lv<NLEVEL; lv++)
      {
         OverlapMPI_FluSyncN    [lv] = 0;
         OverlapMPI_FluSyncPID0 [lv] = NULL;
         OverlapMPI_FluAsyncN   [lv] = 0;
         OverlapMPI_FluAsyncPID0[lv] = NULL;
#        ifdef GRAVITY
         OverlapMPI_PotSyncN    [lv] = 0;
         OverlapMPI_PotSyncPID0 [lv] = NULL;
         OverlapMPI_PotAsyncN   [lv] = 0;
         OverlapMPI_PotAsyncPID0[lv] = NULL;
#        endif
         CutPoint               [lv] = new long [MPI_NRank+1];
         IdxList_Real           [lv] = NULL;
         IdxList_Real_IdxTable  [lv] = NULL;
         PaddedCr1DList         [lv] = NULL;
         PaddedCr1DList_IdxTable[lv] = NULL;

         SendH_NList            [lv] = new int   [MPI_NRank];
         SendH_IDList           [lv] = new int*  [MPI_NRank];
         SendH_SibList          [lv] = new int*  [MPI_NRank];
         SendH_SibDiffList      [lv] = new int*  [MPI_NRank];
         SendH_LBIdxList        [lv] = new long* [MPI_NRank];
         RecvH_NList            [lv] = new int   [MPI_NRank];
         RecvH_IDList           [lv] = new int*  [MPI_NRank];
         RecvH_IDList_IdxTable  [lv] = new int*  [MPI_NRank];
         RecvH_SibList          [lv] = new int*  [MPI_NRank];
         RecvH_SibDiffList      [lv] = new int*  [MPI_NRank];
         RecvH_LBIdxList        [lv] = new long* [MPI_NRank];
         RecvH_PCr1D            [lv] = new ulong*[MPI_NRank];
         RecvH_PCr1D_IdxTable   [lv] = new int*  [MPI_NRank];

         SendX_NList            [lv] = new int   [MPI_NRank];
         SendX_NResList         [lv] = new int   [MPI_NRank];
         SendX_IDList           [lv] = new int*  [MPI_NRank];
         SendX_SibList          [lv] = new int*  [MPI_NRank];
         RecvX_NList            [lv] = new int   [MPI_NRank];
         RecvX_NResList         [lv] = new int   [MPI_NRank];
         RecvX_IDList           [lv] = new int*  [MPI_NRank];
         RecvX_SibList          [lv] = new int*  [MPI_NRank];

         SendR_NList            [lv] = new int   [MPI_NRank];
         SendR_IDList           [lv] = new int*  [MPI_NRank];
         SendR_IDList_IdxTable  [lv] = new int*  [MPI_NRank];
         RecvR_NList            [lv] = new int   [MPI_NRank];
         RecvR_IDList           [lv] = new int*  [MPI_NRank];

         SendF_NList            [lv] = new int   [MPI_NRank];
         SendF_IDList           [lv] = new int*  [MPI_NRank];
         SendF_SibList          [lv] = new int*  [MPI_NRank];
         RecvF_NList            [lv] = new int   [MPI_NRank];
         RecvF_IDList           [lv] = new int*  [MPI_NRank];
         RecvF_IDList_IdxTable  [lv] = new int*  [MPI_NRank];
         RecvF_SibList          [lv] = new int*  [MPI_NRank];

#        ifdef GRAVITY
         SendG_NList            [lv] = new int   [MPI_NRank];
         SendG_IDList           [lv] = new int*  [MPI_NRank];
         SendG_SibList          [lv] = new int*  [MPI_NRank];
         SendG_SibDiffList      [lv] = new int*  [MPI_NRank];
         SendG_LBIdxList        [lv] = new long* [MPI_NRank];
         RecvG_NList            [lv] = new int   [MPI_NRank];
         RecvG_IDList           [lv] = new int*  [MPI_NRank];
         RecvG_IDList_IdxTable  [lv] = new int*  [MPI_NRank];
         RecvG_SibList          [lv] = new int*  [MPI_NRank];
         RecvG_SibDiffList      [lv] = new int*  [MPI_NRank];
         RecvG_LBIdxList        [lv] = new long* [MPI_NRank];
         RecvG_PCr1D            [lv] = new ulong*[MPI_NRank];
         RecvG_PCr1D_IdxTable   [lv] = new int*  [MPI_NRank];
#        endif

         for (int r=0; r<MPI_NRank; r++)
         {
            CutPoint             [lv][r] = -1;

            SendH_NList          [lv][r] = 0;
            SendH_IDList         [lv][r] = NULL;
            SendH_SibList        [lv][r] = NULL;
            SendH_SibDiffList    [lv][r] = NULL;
            SendH_LBIdxList      [lv][r] = NULL;
            RecvH_NList          [lv][r] = 0;
            RecvH_IDList         [lv][r] = NULL;
            RecvH_IDList_IdxTable[lv][r] = NULL;
            RecvH_SibList        [lv][r] = NULL;
            RecvH_SibDiffList    [lv][r] = NULL;
            RecvH_LBIdxList      [lv][r] = NULL;
            RecvH_PCr1D          [lv][r] = NULL;
            RecvH_PCr1D_IdxTable [lv][r] = NULL;

            SendX_NList          [lv][r] = 0;
            SendX_NResList       [lv][r] = 0;
            SendX_IDList         [lv][r] = NULL;
            SendX_SibList        [lv][r] = NULL;
            RecvX_NList          [lv][r] = 0;
            RecvX_NResList       [lv][r] = 0;
            RecvX_IDList         [lv][r] = NULL;
            RecvX_SibList        [lv][r] = NULL;

            SendR_NList          [lv][r] = 0;
            SendR_IDList         [lv][r] = NULL;
            RecvR_NList          [lv][r] = 0;
            SendR_IDList_IdxTable[lv][r] = NULL;
            RecvR_IDList         [lv][r] = NULL;

            SendF_NList          [lv][r] = 0;
            SendF_IDList         [lv][r] = NULL;
            SendF_SibList        [lv][r] = NULL;
            RecvF_NList          [lv][r] = 0;
            RecvF_IDList         [lv][r] = NULL;
            RecvF_IDList_IdxTable[lv][r] = NULL;
            RecvF_SibList        [lv][r] = NULL;

#           ifdef GRAVITY
            SendG_NList          [lv][r] = 0;
            SendG_IDList         [lv][r] = NULL;
            SendG_SibList        [lv][r] = NULL;
            SendG_SibDiffList    [lv][r] = NULL;
            SendG_LBIdxList      [lv][r] = NULL;
            RecvG_NList          [lv][r] = 0;
            RecvG_IDList         [lv][r] = NULL;
            RecvG_IDList_IdxTable[lv][r] = NULL;
            RecvG_SibList        [lv][r] = NULL;
            RecvG_SibDiffList    [lv][r] = NULL;
            RecvG_LBIdxList      [lv][r] = NULL;
            RecvG_PCr1D          [lv][r] = NULL;
            RecvG_PCr1D_IdxTable [lv][r] = NULL;
#           endif
         } // for (int r=0; r<MPI_NRank; r++)
      } // for (int lv=0; lv<NLEVEL; lv++)

   } // METHOD : LB_t



   //===================================================================================
   // Constructor :  ~LB_t
   // Description :  Destructor of the structure "LB_t"
   //
   // Note        :  1. Deallocate memory previously allocated
   //                2. Memory allocated by "malloc/realloc" are deallocated by "free"
   //===================================================================================
   ~LB_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
//       release memory whose size depends on the number of patches at each rank
         reset( lv );


//       release memory whose size is independent of the number of patches at each rank
//       miscellaneous
         if ( CutPoint [lv] != NULL )  delete [] CutPoint[lv];
         CutPoint[lv] = NULL;

//       NList
         if ( SendH_NList   [lv] != NULL )   delete [] SendH_NList   [lv];
         if ( RecvH_NList   [lv] != NULL )   delete [] RecvH_NList   [lv];
         if ( SendX_NList   [lv] != NULL )   delete [] SendX_NList   [lv];
         if ( SendX_NResList[lv] != NULL )   delete [] SendX_NResList[lv];
         if ( RecvX_NList   [lv] != NULL )   delete [] RecvX_NList   [lv];
         if ( RecvX_NResList[lv] != NULL )   delete [] RecvX_NResList[lv];
         if ( SendR_NList   [lv] != NULL )   delete [] SendR_NList   [lv];
         if ( RecvR_NList   [lv] != NULL )   delete [] RecvR_NList   [lv];
         if ( SendF_NList   [lv] != NULL )   delete [] SendF_NList   [lv];
         if ( RecvF_NList   [lv] != NULL )   delete [] RecvF_NList   [lv];
#        ifdef GRAVITY
         if ( SendG_NList   [lv] != NULL )   delete [] SendG_NList   [lv];
         if ( RecvG_NList   [lv] != NULL )   delete [] RecvG_NList   [lv];
#        endif

         SendH_NList   [lv] = NULL;
         RecvH_NList   [lv] = NULL;
         SendX_NList   [lv] = NULL;
         SendX_NResList[lv] = NULL;
         RecvX_NList   [lv] = NULL;
         RecvX_NResList[lv] = NULL;
         SendR_NList   [lv] = NULL;
         RecvR_NList   [lv] = NULL;
         SendF_NList   [lv] = NULL;
         RecvF_NList   [lv] = NULL;
#        ifdef GRAVITY
         SendG_NList   [lv] = NULL;
         RecvG_NList   [lv] = NULL;
#        endif

//       IDList
         if ( SendH_IDList         [lv] != NULL )  delete [] SendH_IDList         [lv];
         if ( SendH_SibList        [lv] != NULL )  delete [] SendH_SibList        [lv];
         if ( SendH_SibDiffList    [lv] != NULL )  delete [] SendH_SibDiffList    [lv];
         if ( SendH_LBIdxList      [lv] != NULL )  delete [] SendH_LBIdxList      [lv];
         if ( RecvH_IDList         [lv] != NULL )  delete [] RecvH_IDList         [lv];
         if ( RecvH_IDList_IdxTable[lv] != NULL )  delete [] RecvH_IDList_IdxTable[lv];
         if ( RecvH_SibList        [lv] != NULL )  delete [] RecvH_SibList        [lv];
         if ( RecvH_SibDiffList    [lv] != NULL )  delete [] RecvH_SibDiffList    [lv];
         if ( RecvH_LBIdxList      [lv] != NULL )  delete [] RecvH_LBIdxList      [lv];
         if ( RecvH_PCr1D          [lv] != NULL )  delete [] RecvH_PCr1D          [lv];
         if ( RecvH_PCr1D_IdxTable [lv] != NULL )  delete [] RecvH_PCr1D_IdxTable [lv];
         if ( SendX_IDList         [lv] != NULL )  delete [] SendX_IDList         [lv];
         if ( SendX_SibList        [lv] != NULL )  delete [] SendX_SibList        [lv];
         if ( RecvX_IDList         [lv] != NULL )  delete [] RecvX_IDList         [lv];
         if ( RecvX_SibList        [lv] != NULL )  delete [] RecvX_SibList        [lv];
         if ( SendR_IDList         [lv] != NULL )  delete [] SendR_IDList         [lv];
         if ( SendR_IDList_IdxTable[lv] != NULL )  delete [] SendR_IDList_IdxTable[lv];
         if ( RecvR_IDList         [lv] != NULL )  delete [] RecvR_IDList         [lv];
         if ( SendF_IDList         [lv] != NULL )  delete [] SendF_IDList         [lv];
         if ( SendF_SibList        [lv] != NULL )  delete [] SendF_SibList        [lv];
         if ( RecvF_IDList         [lv] != NULL )  delete [] RecvF_IDList         [lv];
         if ( RecvF_IDList_IdxTable[lv] != NULL )  delete [] RecvF_IDList_IdxTable[lv];
         if ( RecvF_SibList        [lv] != NULL )  delete [] RecvF_SibList        [lv];
#        ifdef GRAVITY
         if ( SendG_IDList         [lv] != NULL )  delete [] SendG_IDList         [lv];
         if ( SendG_SibList        [lv] != NULL )  delete [] SendG_SibList        [lv];
         if ( SendG_SibDiffList    [lv] != NULL )  delete [] SendG_SibDiffList    [lv];
         if ( SendG_LBIdxList      [lv] != NULL )  delete [] SendG_LBIdxList      [lv];
         if ( RecvG_IDList         [lv] != NULL )  delete [] RecvG_IDList         [lv];
         if ( RecvG_IDList_IdxTable[lv] != NULL )  delete [] RecvG_IDList_IdxTable[lv];
         if ( RecvG_SibList        [lv] != NULL )  delete [] RecvG_SibList        [lv];
         if ( RecvG_SibDiffList    [lv] != NULL )  delete [] RecvG_SibDiffList    [lv];
         if ( RecvG_LBIdxList      [lv] != NULL )  delete [] RecvG_LBIdxList      [lv];
         if ( RecvG_PCr1D          [lv] != NULL )  delete [] RecvG_PCr1D          [lv];
         if ( RecvG_PCr1D_IdxTable [lv] != NULL )  delete [] RecvG_PCr1D_IdxTable [lv];
#        endif

         SendH_IDList         [lv] = NULL;
         SendH_SibList        [lv] = NULL;
         SendH_SibDiffList    [lv] = NULL;
         SendH_LBIdxList      [lv] = NULL;
         RecvH_IDList         [lv] = NULL;
         RecvH_IDList_IdxTable[lv] = NULL;
         RecvH_SibList        [lv] = NULL;
         RecvH_SibDiffList    [lv] = NULL;
         RecvH_LBIdxList      [lv] = NULL;
         RecvH_PCr1D          [lv] = NULL;
         RecvH_PCr1D_IdxTable [lv] = NULL;
         SendX_IDList         [lv] = NULL;
         SendX_SibList        [lv] = NULL;
         RecvX_IDList         [lv] = NULL;
         RecvX_SibList        [lv] = NULL;
         SendR_IDList         [lv] = NULL;
         SendR_IDList_IdxTable[lv] = NULL;
         RecvR_IDList         [lv] = NULL;
         SendF_IDList         [lv] = NULL;
         SendF_SibList        [lv] = NULL;
         RecvF_IDList         [lv] = NULL;
         RecvF_IDList_IdxTable[lv] = NULL;
         RecvF_SibList        [lv] = NULL;
#        ifdef GRAVITY
         SendG_IDList         [lv] = NULL;
         SendG_SibList        [lv] = NULL;
         SendG_SibDiffList    [lv] = NULL;
         SendG_LBIdxList      [lv] = NULL;
         RecvG_IDList         [lv] = NULL;
         RecvG_IDList_IdxTable[lv] = NULL;
         RecvG_SibList        [lv] = NULL;
         RecvG_SibDiffList    [lv] = NULL;
         RecvG_LBIdxList      [lv] = NULL;
         RecvG_PCr1D          [lv] = NULL;
         RecvG_PCr1D_IdxTable [lv] = NULL;
#        endif
      } // for (int lv=0; lv<NLEVEL; lv++)

   } // METHOD : ~LB_t



   //===================================================================================
   // Constructor :  reset
   // Description :  Reset memory for re-distributing all patches at all levels
   //
   // Note        :  1. Pointers whose sizes depend on the number of MPI ranks are NOT deallocated
   //                2. This function must NOT reset CutPoint[] since in LB_Init_LoadBalance() we need to
   //                   call LB->reset() AFTER setting CutPoint[]
   //                3. WLI will be reset to NULL_REAL even though this function just resets
   //                   load-balance variables on a specific level
   //
   // Parameter   :  lv : Target refinement level
   //===================================================================================
   void reset( const int lv )
   {
      WLI = NULL_REAL;

//    miscellaneous
      if ( OverlapMPI_FluSyncPID0 [lv] != NULL )   delete [] OverlapMPI_FluSyncPID0 [lv];
      if ( OverlapMPI_FluAsyncPID0[lv] != NULL )   delete [] OverlapMPI_FluAsyncPID0[lv];
#     ifdef GRAVITY
      if ( OverlapMPI_PotSyncPID0 [lv] != NULL )   delete [] OverlapMPI_PotSyncPID0 [lv];
      if ( OverlapMPI_PotAsyncPID0[lv] != NULL )   delete [] OverlapMPI_PotAsyncPID0[lv];
#     endif
      if ( IdxList_Real           [lv] != NULL )   delete [] IdxList_Real           [lv];
      if ( IdxList_Real_IdxTable  [lv] != NULL )   delete [] IdxList_Real_IdxTable  [lv];
      if ( PaddedCr1DList         [lv] != NULL )   free(     PaddedCr1DList         [lv] );
      if ( PaddedCr1DList_IdxTable[lv] != NULL )   free(     PaddedCr1DList_IdxTable[lv] );

      OverlapMPI_FluSyncPID0 [lv] = NULL;
      OverlapMPI_FluAsyncPID0[lv] = NULL;
#     ifdef GRAVITY
      OverlapMPI_PotSyncPID0 [lv] = NULL;
      OverlapMPI_PotAsyncPID0[lv] = NULL;
#     endif
      IdxList_Real           [lv] = NULL;
      IdxList_Real_IdxTable  [lv] = NULL;
      PaddedCr1DList         [lv] = NULL;
      PaddedCr1DList_IdxTable[lv] = NULL;

      for (int r=0; r<MPI_NRank; r++)
      {
//       NList
         SendH_NList   [lv][r] = 0;
         RecvH_NList   [lv][r] = 0;
         SendX_NList   [lv][r] = 0;
         SendX_NResList[lv][r] = 0;
         RecvX_NList   [lv][r] = 0;
         RecvX_NResList[lv][r] = 0;
         SendR_NList   [lv][r] = 0;
         RecvR_NList   [lv][r] = 0;
         SendF_NList   [lv][r] = 0;
         RecvF_NList   [lv][r] = 0;
#        ifdef GRAVITY
         SendG_NList   [lv][r] = 0;
         RecvG_NList   [lv][r] = 0;
#        endif

//       IDList
         if ( r == 0 )
         {
         if ( SendH_IDList         [lv][r] != NULL )  delete [] SendH_IDList         [lv][r];
         if ( SendH_SibList        [lv][r] != NULL )  delete [] SendH_SibList        [lv][r];
         if ( SendH_SibDiffList    [lv][r] != NULL )  delete [] SendH_SibDiffList    [lv][r];
         if ( SendH_LBIdxList      [lv][r] != NULL )  delete [] SendH_LBIdxList      [lv][r];
         if ( RecvH_IDList_IdxTable[lv][r] != NULL )  delete [] RecvH_IDList_IdxTable[lv][r];
         if ( RecvH_SibList        [lv][r] != NULL )  delete [] RecvH_SibList        [lv][r];
         if ( RecvH_SibDiffList    [lv][r] != NULL )  delete [] RecvH_SibDiffList    [lv][r];
         if ( RecvH_LBIdxList      [lv][r] != NULL )  delete [] RecvH_LBIdxList      [lv][r];
         if ( RecvH_PCr1D          [lv][r] != NULL )  delete [] RecvH_PCr1D          [lv][r];
         if ( RecvH_PCr1D_IdxTable [lv][r] != NULL )  delete [] RecvH_PCr1D_IdxTable [lv][r];
         }
         if ( RecvH_IDList         [lv][r] != NULL )  free(     RecvH_IDList         [lv][r] );
         if ( SendX_IDList         [lv][r] != NULL )  delete [] SendX_IDList         [lv][r];
         if ( SendX_SibList        [lv][r] != NULL )  delete [] SendX_SibList        [lv][r];
         if ( RecvX_IDList         [lv][r] != NULL )  delete [] RecvX_IDList         [lv][r];
         if ( RecvX_SibList        [lv][r] != NULL )  delete [] RecvX_SibList        [lv][r];
         if ( SendR_IDList         [lv][r] != NULL )  free(     SendR_IDList         [lv][r] );
         if ( SendR_IDList_IdxTable[lv][r] != NULL )  delete [] SendR_IDList_IdxTable[lv][r];
         if ( RecvR_IDList         [lv][r] != NULL )  delete [] RecvR_IDList         [lv][r];
         if ( SendF_IDList         [lv][r] != NULL )  delete [] SendF_IDList         [lv][r];
         if ( SendF_SibList        [lv][r] != NULL )  delete [] SendF_SibList        [lv][r];
         if ( RecvF_IDList         [lv][r] != NULL )  free(     RecvF_IDList         [lv][r] );
         if ( RecvF_IDList_IdxTable[lv][r] != NULL )  delete [] RecvF_IDList_IdxTable[lv][r];
         if ( RecvF_SibList        [lv][r] != NULL )  delete [] RecvF_SibList        [lv][r];
#        ifdef GRAVITY
         if ( r == 0 )
         {
         if ( SendG_IDList         [lv][r] != NULL )  delete [] SendG_IDList         [lv][r];
         if ( SendG_SibList        [lv][r] != NULL )  delete [] SendG_SibList        [lv][r];
         if ( SendG_SibDiffList    [lv][r] != NULL )  delete [] SendG_SibDiffList    [lv][r];
         if ( SendG_LBIdxList      [lv][r] != NULL )  delete [] SendG_LBIdxList      [lv][r];
         if ( RecvG_IDList_IdxTable[lv][r] != NULL )  delete [] RecvG_IDList_IdxTable[lv][r];
         if ( RecvG_SibList        [lv][r] != NULL )  delete [] RecvG_SibList        [lv][r];
         if ( RecvG_SibDiffList    [lv][r] != NULL )  delete [] RecvG_SibDiffList    [lv][r];
         if ( RecvG_LBIdxList      [lv][r] != NULL )  delete [] RecvG_LBIdxList      [lv][r];
         if ( RecvG_PCr1D          [lv][r] != NULL )  delete [] RecvG_PCr1D          [lv][r];
         if ( RecvG_PCr1D_IdxTable [lv][r] != NULL )  delete [] RecvG_PCr1D_IdxTable [lv][r];
         }
         if ( RecvG_IDList         [lv][r] != NULL )  free(     RecvG_IDList         [lv][r] );
#        endif

         SendH_IDList         [lv][r] = NULL;
         SendH_SibList        [lv][r] = NULL;
         SendH_SibDiffList    [lv][r] = NULL;
         SendH_LBIdxList      [lv][r] = NULL;
         RecvH_IDList         [lv][r] = NULL;
         RecvH_IDList_IdxTable[lv][r] = NULL;
         RecvH_SibList        [lv][r] = NULL;
         RecvH_SibDiffList    [lv][r] = NULL;
         RecvH_LBIdxList      [lv][r] = NULL;
         RecvH_PCr1D          [lv][r] = NULL;
         RecvH_PCr1D_IdxTable [lv][r] = NULL;
         SendX_IDList         [lv][r] = NULL;
         SendX_SibList        [lv][r] = NULL;
         RecvX_IDList         [lv][r] = NULL;
         RecvX_SibList        [lv][r] = NULL;
         SendR_IDList         [lv][r] = NULL;
         SendR_IDList_IdxTable[lv][r] = NULL;
         RecvR_IDList         [lv][r] = NULL;
         SendF_IDList         [lv][r] = NULL;
         SendF_SibList        [lv][r] = NULL;
         RecvF_IDList         [lv][r] = NULL;
         RecvF_IDList_IdxTable[lv][r] = NULL;
         RecvF_SibList        [lv][r] = NULL;
#        ifdef GRAVITY
         SendG_IDList         [lv][r] = NULL;
         SendG_SibList        [lv][r] = NULL;
         SendG_SibDiffList    [lv][r] = NULL;
         SendG_LBIdxList      [lv][r] = NULL;
         RecvG_IDList         [lv][r] = NULL;
         RecvG_IDList_IdxTable[lv][r] = NULL;
         RecvG_SibList        [lv][r] = NULL;
         RecvG_SibDiffList    [lv][r] = NULL;
         RecvG_LBIdxList      [lv][r] = NULL;
         RecvG_PCr1D          [lv][r] = NULL;
         RecvG_PCr1D_IdxTable [lv][r] = NULL;
#        endif
      } // for (int r=0; r<MPI_NRank; r++)

   } // METHOD : reset


}; // struct LB_t



#endif // #ifndef __LOADBALANCE_H__
