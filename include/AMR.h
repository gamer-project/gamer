#ifndef __AMR_H__
#define __AMR_H__



#include "Macro.h"
#include "Patch.h"

#ifdef PARTICLE
#  include "Particle.h"
#endif

#ifndef SERIAL
#  include "ParaVar.h"
#  ifdef LOAD_BALANCE
#  include "LoadBalance.h"
#  endif
#else
#  ifdef LOAD_BALANCE
#  error ERROR : options LOAD_BALANCE and SERIAL should NOT be turned on at the same time
#  endif
#endif

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  AMR_t
// Description :  Data structure of the AMR implementation
//
// Data Member :  patch       : Pointers of all patches
//                num         : Number of patches (real patch + buffer patch) at each level
//                scale       : Grid scale at each level (grid size normalized to that at the finest level)
//                FluSg       : Sandglass of the current fluid     data [0/1]
//                PotSg       : Sandglass of the current potential data [0/1]
//                FluSgTime   : Physical time of FluSg
//                PotSgTime   : Physical time of PotSg
//                NPatchComma : (1) SERIAL: [1] = [2] = ... = [26] = num[lv] = total number of patches
//                              (2) Parallel, but no LOAD_BALANCE:
//                                  [ 0, start of buffer patches [s=0], start of buffer patches [s=1]
//                                   ... start of buffer patches [s=25], total # of patches (=num[lv]) ]
//                              (3) Parallel, with LOAD_BALANCE:
//                                  [ 0, start of sibling-buffer patches, start of father-buffer patches,
//                                    total # of patches (=num[lv]), same as [3], ...]
//                               --> In all cases, [ 1] gives the total number of "real" patches
//                                                 [26] gives the total number of "real+buffer" patches
//                BoxEdgeL    : Simulation box left  edge in the adopted coordinate system
//                BoxEdgeR    : Simulation box right edge in the adopted coordinate system
//                BoxCenter   : Simulation box center     in the adopted coordinate system
//                dh          : Grid size at each level
//                BoxSize     : Simulation box physical size
//                BoxScale    : Simulation box scale
//                WithFlux    : Whether of not to allocate the flux arrays at all coarse-fine boundaries
//                Par         : Particle data
//                ParaVar     : Variables for parallelization
//                LB          : Variables for load-balance
//                ResPower2   : ceil(  log2( effective resolution at lv )  ) --> mainly used by LOAD_BALANCE
//                NUpdateLv   : Number of updates at each level in one global time-step
//                              --> Do not take into account the number of patches and particles at each level
//                              --> Mainly used for estimating the weighted load-imbalance factor to determine
//                                  when to redistribute all patches (when LOAD_BALANCE is on)
//
// Method      :  AMR_t    : Constructor
//               ~AMR_t    : Destructor
//                pnew     : Allocate one patch
//                pdelete  : Deallocate one patch
//                Lvdelete : Deallocate all patches in the given level
//-------------------------------------------------------------------------------------------------------
struct AMR_t
{

// data members
// ===================================================================================
   patch_t    *patch[2][NLEVEL][MAX_PATCH];

#  ifdef PARTICLE
   Particle_t *Par;
#  endif

#  ifndef SERIAL
   ParaVar_t  *ParaVar;
#  ifdef LOAD_BALANCE
   LB_t       *LB;
#  endif
#  endif

   int    num         [NLEVEL];
   int    scale       [NLEVEL];
   int    FluSg       [NLEVEL];
   double FluSgTime   [NLEVEL][2];
#  ifdef GRAVITY
   int    PotSg       [NLEVEL];
   double PotSgTime   [NLEVEL][2];
#  endif
   int    NPatchComma [NLEVEL][28];
   double dh          [NLEVEL];
   int    ResPower2   [NLEVEL];
   double BoxEdgeL    [3];
   double BoxEdgeR    [3];
   double BoxCenter   [3];
   double BoxSize     [3];
   int    BoxScale    [3];
   bool   WithFlux;
   long   NUpdateLv   [NLEVEL];



   //===================================================================================
   // Constructor :  AMR_t
   // Description :  Constructor of the structure "AMR_t"
   //
   // Note        :  Initialize the data members
   //===================================================================================
   AMR_t()
   {

      for (int lv=0; lv<NLEVEL; lv++)
      {
         num  [lv] = 0;
         scale[lv] = 1<<(NLEVEL-1-lv);
         FluSg[lv] = 0;
#        ifdef GRAVITY
         PotSg[lv] = FluSg[lv];
#        endif

//       initialized as arbitrary "negative" number to indicate that they have not been set yet
//       --> these will be reset by Init_SetDefaultParameter and Init_ByRestart_*
         FluSgTime[lv][   FluSg[lv] ] = -__FLT_MAX__;
         FluSgTime[lv][ 1-FluSg[lv] ] = -__FLT_MAX__;
#        ifdef GRAVITY
         PotSgTime[lv][   PotSg[lv] ] = -__FLT_MAX__;
         PotSgTime[lv][ 1-PotSg[lv] ] = -__FLT_MAX__;
#        endif
      }

      for (int Sg=0; Sg<2; Sg++)
      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<MAX_PATCH; PID++)
         patch[Sg][lv][PID] = NULL;

      for (int lv=0; lv<NLEVEL; lv++)
      for (int m=0; m<28; m++)
         NPatchComma[lv][m] = 0;

#     ifdef PARTICLE
      Par = NULL;
#     endif

#     ifndef SERIAL
      ParaVar = new ParaVar_t;
#     endif

#     ifdef LOAD_BALANCE
      LB = NULL;
#     endif

      WithFlux = false;

   } // METHOD : AMR_t



   //===================================================================================
   // Constructor :  ~AMR_t
   // Description :  Destructor of the structure "AMR_t"
   //
   // Note        :  Deallocate memory and reset parameters
   //===================================================================================
   ~AMR_t()
   {

      const bool ReusePatchMemory_No = false;
      for (int lv=0; lv<NLEVEL; lv++)  Lvdelete( lv, ReusePatchMemory_No );

#     ifdef PARTICLE
      if ( Par != NULL )
      {
         delete Par;
         Par = NULL;
      }
#     endif

#     ifndef SERIAL
      if ( ParaVar != NULL )
      {
         delete ParaVar;
         ParaVar = NULL;
      }
#     endif

#     ifdef LOAD_BALANCE
      if ( LB != NULL )
      {
         delete LB;
         LB = NULL;
      }
#     endif

   } // METHOD : ~AMR_t



   //===================================================================================
   // Method      :  pnew
   // Description :  allocate a single patch
   //
   // Note        :  1. Each patch contains two patch pointers --> SANDGLASS (Sg) = 0 / 1
   //                2. Sg = 0 : Store both data and relation (father,son.sibling,corner,flag,flux)
   //                   Sg = 1 : Store only data
   //
   // Parameter   :  lv       : Target refinement level
   //                x,y,z    : Physical coordinates of the patch corner
   //                FaPID    : Patch ID of the parent patch at level "lv-1"
   //                FluData  : true --> Allocate hydrodynamic array "fluid"
   //                PotData  : true --> Allocate potential array "pot"
   //===================================================================================
   void pnew( const int lv, const int x, const int y, const int z, const int FaPID, const bool FluData,
              const bool PotData )
   {

      const int NewPID = num[lv];

      if ( NewPID > MAX_PATCH-1 )
         Aux_Error( ERROR_INFO, "exceed MAX_PATCH (%d) => please reset it in the Makefile !!\n", MAX_PATCH );

//    allocate new patches if there are no inactive patches
      if ( patch[0][lv][NewPID] == NULL )
      {
#        ifdef GAMER_DEBUG
//       assuming Sg=0/1 are always allocated and deallocated together
         if ( patch[1][lv][NewPID] != NULL )
            Aux_Error( ERROR_INFO, "conflicting patch allocation (Lv %d, PID %d, FaPID %d) !!\n", lv, NewPID, FaPID );
#        endif

         patch[0][lv][NewPID] = new patch_t( x, y, z, FaPID, FluData, PotData, FluData, lv, BoxScale, dh[TOP_LEVEL] );
         patch[1][lv][NewPID] = new patch_t( 0, 0, 0,    -1, FluData, PotData,   false, lv, BoxScale, dh[TOP_LEVEL] );
      }

//    reactivate inactive patches
      else
      {
#        ifdef GAMER_DEBUG
//       assuming Sg=0/1 are always allocated and deallocated together
         if ( patch[1][lv][NewPID] == NULL )
            Aux_Error( ERROR_INFO, "conflicting patch allocation (Lv %d, PID %d, FaPID %d) !!\n", lv, NewPID, FaPID );

         for (int Sg=0; Sg<2; Sg++)
            if ( patch[Sg][lv][NewPID]->Active )
               Aux_Error( ERROR_INFO, "this patch has been activated already (Lv %d, PID %d, FaPID %d, Sg %d) !!\n",
                          lv, NewPID, FaPID, Sg );
#        endif

//       do NOT initialize field pointers as NULL since they may be allocated already
         const bool InitPtrAsNull_No = false;

         patch[0][lv][NewPID]->Activate( x, y, z, FaPID, FluData, PotData, FluData, lv, BoxScale, dh[TOP_LEVEL], InitPtrAsNull_No );
         patch[1][lv][NewPID]->Activate( 0, 0, 0,    -1, FluData, PotData,   false, lv, BoxScale, dh[TOP_LEVEL], InitPtrAsNull_No );
      } // if ( patch[0][lv][NewPID] == NULL ) ... else ...

      num[lv] ++;

   } // METHOD : pnew



   //===================================================================================
   // Method      :  pdelete
   // Description :  Deallocate a single patch
   //
   // Note        :  1. This function should NOT be applied to the base-level patches (unless
   //                   the option "LOAD_BALANCE" is turned on, in which the base-level patches need
   //                   to be redistributed)
   //                2. This function will also deallocate the flux arrays of the target patch
   //                3. Delete a patch with son is forbidden
   //                4. Delete a patch with home particles is forbidden
   //
   // Parameter   :  lv          : Target refinement level
   //                PID         : Patch ID to be removed
   //                ReuseMemory : true  --> mark patch as inactive, but do not deallocate memory (for OPT__REUSE_MEMORY)
   //                              false --> deallocate patch
   //===================================================================================
   void pdelete( const int lv, const int PID, const bool ReuseMemory )
   {

#     ifdef GAMER_DEBUG
      if ( patch[0][lv][PID] == NULL  ||  patch[1][lv][PID] == NULL )
         Aux_Error( ERROR_INFO, "delete a non-existing patch (Lv %d, PID %d) !!\n", lv, PID );

      if ( patch[0][lv][PID]->son != -1 )
         Aux_Error( ERROR_INFO, "delete a patch with son (Lv %d, PID %d, SonPID %d) !!\n",
                    lv, PID, patch[0][lv][PID]->son );
#     ifdef PARTICLE
      if ( patch[0][lv][PID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "delete a patch with home particles (Lv %d, PID %d, NPar %d) !!\n",
                    lv, PID, patch[0][lv][PID]->NPar );
#     endif
#     endif // #ifdef GAMER_DEBUG

      if ( ReuseMemory )
      {
#        ifdef GAMER_DEBUG
         for (int Sg=0; Sg<2; Sg++)
            if ( !patch[Sg][lv][PID]->Active )
               Aux_Error( ERROR_INFO, "this patch has been deactivated already (Lv %d, PID %d, Sg %d) !!\n", lv, PID, Sg );
#        endif

         patch[0][lv][PID]->Active = false;
         patch[1][lv][PID]->Active = false;

//       always deallocate flux arrays because
//       (1) we use them to determine which patches require the flux fix-up operation (see Flu_FixUp.cpp)
//       (2) flux arrays do not consume much memory (at most 6/32, where 6 = 6 faces and 32 = patch group size*two sg)
//       (3) different patches may require flux arrays along different directions, and thus allocating memory pool for flux arrays
//           can be inefficient and less useful
         patch[0][lv][PID]->fdelete();
      }

      else
      {
         delete patch[0][lv][PID];
         delete patch[1][lv][PID];

         patch[0][lv][PID] = NULL;
         patch[1][lv][PID] = NULL;
      } // if ( ReuseMemory ) ... else ...

      num[lv] --;

#     ifdef GAMER_DEBUG
      if ( num[lv] < 0 )
         Aux_Error( ERROR_INFO, "num[%d] = %d < 0 !!\n", lv, num[lv] );
#     endif

   } // METHOD : pdelete



   //===================================================================================
   // Method      :  Lvdelete
   // Description :  Deallocate all patches in the target level and initialize all
   //                parameters as the default values
   //
   // Note        :  1. This function will delete a patch even if it has sons (and particles)
   //                2. This function will scan over amr->num[lv] patches
   //                3. The variables "scale, FluSg, PotSg, and dh" will NOT be modified
   //
   // Parameter   :  lv               : Target refinement level
   //                ReusePatchMemory : true  --> mark patch as inactive, but do not deallocate memory (for OPT__REUSE_MEMORY)
   //                                   false --> deallocate patch
   //===================================================================================
   void Lvdelete( const int lv, const bool ReusePatchMemory )
   {

#     ifdef GAMER_DEBUG
      if ( lv < 0  ||  lv >= NLEVEL )
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "lv", lv );
#     endif

//    store the total number of patches to be deleted since pdelete will modify num
      const int NPatch2Delete = num[lv];
      for (int PID=0; PID<NPatch2Delete; PID++)
      {
//       reset son=-1 to skip the check in pdelete
         patch[0][lv][PID]->son = -1;

         pdelete( lv, PID, ReusePatchMemory );
      }

      for (int m=0; m<28; m++)   NPatchComma[lv][m] = 0;

#     ifndef SERIAL
      if ( ParaVar != NULL )     ParaVar->Lvdelete( lv );
#     endif

   } // METHOD : Lvdelete


}; // struct AMR_t



#endif // #ifndef __AMR_H__
