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
//                dh          : Grid size at each level
//                BoxSize     : Simulation box physical size
//                BoxScale    : Simulation box scale
//                WithFlux    : Whether of not to allocate the flux arrays at all coarse-fine boundaries
//                Par         : Particle data
//                ParaVar     : Variables for parallelization
//                LB          : Variables for load-balance
//                ResPower2   : ceil(  log2( effective resolution at lv )  ) --> mainly used by LOAD_BALANCE
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
   patch_t   *patch[2][NLEVEL][MAX_PATCH];

#  ifdef PARTICLE
   Particle_t  *Par;
#  endif

#  ifndef SERIAL
   ParaVar_t *ParaVar;
#  ifdef LOAD_BALANCE
   LB_t      *LB;
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
   double BoxSize     [3];
   int    BoxScale    [3];
   bool   WithFlux;
   


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

         FluSgTime[lv][   FluSg[lv] ] = 0.0;          // must be initialized as zero
         FluSgTime[lv][ 1-FluSg[lv] ] = -__FLT_MAX__; // initialized as negative for debug
#        ifdef GRAVITY
         PotSgTime[lv][   PotSg[lv] ] = 0.0;          // must be initialized as zero
         PotSgTime[lv][ 1-PotSg[lv] ] = -__FLT_MAX__; // initialized as negative for debug
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

      for (int lv=0; lv<NLEVEL; lv++)  Lvdelete( lv );

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
   // Note        :  a. Each patch contains two patch pointers --> SANDGLASS (Sg) = 0 / 1 
   //                b. Sg = 0 : Store both data and relation (father,son.sibling,corner,flag,flux)
   //                   Sg = 1 : Store only data 
   //
   // Parameter   :  lv       : Targeted refinement level
   //                x,y,z    : Physical coordinates of the patch corner
   //                FaPID    : Patch ID of the parent patch at level "lv-1"
   //                FluData  : true --> Allocate hydrodynamic array "fluid"
   //                PotData  : true --> Allocate potential array "pot"
   //===================================================================================
   void pnew( const int lv, const int x, const int y, const int z, const int FaPID, const bool FluData,
              const bool PotData )
   {

      if ( num[lv]+1 > MAX_PATCH )
         Aux_Error( ERROR_INFO, "exceed MAX_PATCH (%d) => please reset it in the Makefile !!\n", MAX_PATCH );

#     ifdef GAMER_DEBUG
      if ( patch[0][lv][num[lv]] != NULL  ||  patch[1][lv][num[lv]] != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing patch (Lv %d, PID %d, FaPID %d) !!\n",
                    lv, num[lv], FaPID );
#     endif

      patch[0][lv][ num[lv] ] = new patch_t( x, y, z, FaPID, FluData, PotData, lv, BoxScale, dh[TOP_LEVEL] );
      patch[1][lv][ num[lv] ] = new patch_t( 0, 0, 0,    -1, FluData, PotData, lv, BoxScale, dh[TOP_LEVEL] );

      num[lv] ++;

   } // METHOD : pnew



   //===================================================================================
   // Method      :  pdelete
   // Description :  Deallocate a single patch 
   //
   // Note        :  a. This function should NOT be applied to the base-level patches (unless
   //                   the option "LOAD_BALANCE" is turned on, in which the base-level patches need
   //                   to be redistributed)
   //                b. This function will also deallocate the flux arrays of the targeted patch
   //                c. Delete a patch with son is forbidden
   //                d. Delete a patch with home particles is forbidden
   //
   // Parameter   :  lv  : The targeted refinement level
   //                PID : The patch ID to be removed
   //===================================================================================
   void pdelete( const int lv, const int PID )
   {

#     ifdef GAMER_DEBUG
#     ifndef LOAD_BALANCE
      if ( lv == 0 )
         Aux_Error( ERROR_INFO, "delete a base-level patch !!\n" );
#     endif

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

      delete patch[0][lv][PID];
      delete patch[1][lv][PID];

      patch[0][lv][PID] = NULL;
      patch[1][lv][PID] = NULL;

      num[lv] --;

#     ifdef GAMER_DEBUG
      if ( num[lv] < 0 )
         Aux_Error( ERROR_INFO, "num[%d] = %d < 0 !!\n", lv, num[lv] );
#     endif

   } // METHOD : pdelete



   //===================================================================================
   // Method      :  Lvdelete
   // Description :  Deallocate all patches in the targeted level and initialize all
   //                parameters as the default values
   //
   // Note        :  a. This function will delete a patch even if it has sons
   //                b. This function will scan over amr->num[lv] patches
   //                c. The variables "scale, FluSg, PotSg, and dh" will NOT be modified
   //
   // Parameter   :  lv : Targeted refinement level
   //===================================================================================
   void Lvdelete( const int lv )
   {

#     ifdef GAMER_DEBUG
      if ( lv < 0  ||  lv >= NLEVEL )
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "lv", lv );
#     endif

      for (int Sg=0; Sg<2; Sg++)
      for (int PID=0; PID<num[lv]; PID++)
      {
#        ifdef GAMER_DEBUG
         if ( patch[Sg][lv][PID] == NULL )  
            Aux_Error( ERROR_INFO, "amr->patch[%d][%d][%d] does not exist (==NULL) !!\n", Sg, lv, PID );
#        endif

         delete patch[Sg][lv][PID];
         patch[Sg][lv][PID] = NULL;
      }

      num[lv] = 0;

      for (int m=0; m<28; m++)   NPatchComma[lv][m] = 0;

#     ifndef SERIAL
      if ( ParaVar != NULL )     ParaVar->Lvdelete( lv );
#     endif

   } // METHOD : Lvdelete


}; // struct AMR_t



#endif // #ifndef __AMR_H__
