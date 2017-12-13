#ifndef __PATCH_H__
#define __PATCH_H__



#include "Macro.h"

#ifdef PARTICLE
#  include <cmath>
#  include "Particle.h"
#endif


void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
#ifdef LOAD_BALANCE
ulong Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] );
long  LB_Corner2Index( const int lv, const int Corner[], const Check_t Check );
#endif




//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t 
// Description :  Data structure of a single patch 
//
// Data Member :  fluid           : Fluid variables (mass density, momentum density x, y ,z, energy density) 
//                passive         : Passively advected variables (e.g., metal density)
//                pot             : Potential
//                pot_ext         : Potential with GRA_GHOST_SIZE ghost cells on each side
//                                  --> Allocated only if STORE_POT_GHOST is on (used for Par->ImproveAcc only)
//                                  --> Ghost-zone potential are obtained from the Poisson solver directly
//                                      (not from exchanging potential between sibling patches)
//                                  --> Currently the base-level patches do not store pot_ext
//                flux[6]         : Fluid flux (for the flux-correction operation)
//                flux_passive[6] : Passive variable flux (for the flux-correction operation)
//                flux_debug[6]   : Fluid flux for the debug mode (ensuring that the round-off errors are 
//                                  exactly the same in different parallelization parameters/strategies)
//                corner[3]       : Grid indices of the cell at patch corner
//                                  --> Note that for an external patch its recorded "corner" will lie outside
//                                      the simulation domain. In other words, periodicity is NOT used to
//                                      convert the corner to that of the corresponding real patch
//                                  --> For example, external patches can have corner < 0
//                sibling[26]     : Patch IDs of the 26 sibling patches (-1->no sibling; -1XX->external)
//
//                                  NOTE FOR NON-PERIODIC BOUNDARY CONDITIONS:
//                                  If a target sibling patch is an external patch (which lies outside the simulation domain),
//                                  the corresponding sibling index is set to "SIB_OFFSET_NONPERIODIC-Sibling", where Sibling 
//                                  represents the sibling direction of the boundary region.
//                                  (e.g., if amr->sibling[XX] lies outside the -y boundary (boundary sibling = 2), we set 
//                                   amr->sibling[XX] = SIB_OFFSET_NONPERIODIC-2 = -102 for SIB_OFFSET_NONPERIODIC==-100)
//                                  --> All patches without sibling patches along the XX direction will have 
//                                      "amr->sibling[XX] < 0"
//
//                father          : Patch ID of the father patch
//                son             : Patch ID of the child patch (-1->no son)
//
//                                  LOAD_BALANCE NOTE: 
//                                  For patches with sons living abroad (not at the same rank as iteself), their son indices
//                                  are set to "SON_OFFSET_LB-SonRank", where SonRank represents the MPI rank where their 
//                                  sons live. (e.g., if the real son at rank 123, the son index is set to 
//                                  SON_OFFSET_LB-123=-1123 for SON_OFFSET_LB==-1000)
//                                  --> All patches with sons will have SonPID != -1
//
//                flag            : Refinement flag (true/false)
//                EdgeL/R         : Left and right edge of the patch
//                PaddedCr1D      : 1D corner coordiniate padded with two base-level patches on each side 
//                                  in each direction, normalized to the finest-level patch scale (PATCH_SIZE)
//                                  --> each PaddedCr1D defines a unique 3D position
//                                  --> patches at different levels with the same PaddedCr1D have the same 
//                                      3D corner coordinates
//                LB_Idx          : Space-filling-curve index for load balance
//                NPar            : Number of particles belonging to this leaf patch
//                ParListSize     : Size of the array ParList (ParListSize can be >= NPar) 
//                ParList         : List recording the IDs of all particles belonging to this leaf patch
//                NPar_Desc       : Number of particles belonging to the descendants (sons, grandsons, ...) of this patch
//                ParList_Desc    : List recording the IDs of all particles belonging to the descendants of this patch
//                NPar_Escp       : Number of particles escaping from this patch
//                ParList_Escp    : List recording the IDs of all particles escaping from this patch
//
// Method      :  patch_t         : Constructor 
//               ~patch_t         : Destructor
//                fnew            : Allocate one flux array 
//                fdelete         : Deallocate one flux array
//                hnew            : Allocate hydrodynamic array
//                hdelete         : Deallocate hydrodynamic array
//                gnew            : Allocate potential array
//                gdelete         : Deallocate potential array
//                AddParticle     : Add particles to the particle list
//                RemoveParticle  : Remove particles from the particle list
//-------------------------------------------------------------------------------------------------------
struct patch_t
{

// data members
// ===================================================================================
   real (*fluid)  [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
#  if ( NPASSIVE > 0 )
   real (*passive)[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
#  endif

#  ifdef GRAVITY
   real (*pot)[PATCH_SIZE][PATCH_SIZE];

#  ifdef STORE_POT_GHOST
   real (*pot_ext)[GRA_NXT][GRA_NXT];
#  endif
#  endif

   real (*flux        [6])[PATCH_SIZE][PATCH_SIZE];
#  if ( NPASSIVE > 0 )
   real (*flux_passive[6])[PATCH_SIZE][PATCH_SIZE];
#  endif
#  ifdef GAMER_DEBUG
   real (*flux_debug  [6])[PATCH_SIZE][PATCH_SIZE];
#  endif

   int  corner[3];
   int  sibling[26];
   int  father;
   int  son;
   bool flag;
   real EdgeL[3];
   real EdgeR[3];

#  ifdef LOAD_BALANCE
   ulong PaddedCr1D;
   long  LB_Idx;
#  endif

#  ifdef PARTICLE
   int   NPar;
   int   ParListSize;
   long *ParList;

   int   NPar_Desc;
   long *ParList_Desc;

   int   NPar_Escp[26];
   long *ParList_Escp[26];
#  endif



   //===================================================================================
   // Constructor :  patch_t 
   // Description :  Constructor of the structure "patch_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  x,y,z    : Scale indices of the patch corner
   //                FaPID    : Patch ID of the father patch    
   //                FluData  : true --> Allocate hydrodynamic array(s) "fluid" (and "passive")
   //                PotData  : true --> Allocate potential array "pot" (has no effect if "GRAVITY" is turned off)
   //                                    "pot_ext" will be allocated as well if STORE_POT_GHOST is on
   //                lv       : Refinement level of the newly created patch
   //                BoxScale : Simulation box scale
   //                dh_min   : Cell size at the maximum level
   //===================================================================================
   patch_t( const int x, const int y, const int z, const int FaPID, const bool FluData, const bool PotData,
            const int lv, const int BoxScale[], const real dh_min )
   {

      corner[0] = x; 
      corner[1] = y;
      corner[2] = z;
      father    = FaPID;
      son       = -1;

#     ifdef LOAD_BALANCE
      const int Padded              = 1<<NLEVEL;
      const int BoxNScale_Padded[3] = { BoxScale[0]/PATCH_SIZE + 2*Padded,
                                        BoxScale[1]/PATCH_SIZE + 2*Padded,
                                        BoxScale[2]/PATCH_SIZE + 2*Padded }; // normalized and padded box scale
      int Cr_Padded[3];
      for (int d=0; d<3; d++)    Cr_Padded[d] = corner[d]/PATCH_SIZE + Padded;

#     ifdef GAMER_DEBUG
      for (int d=0; d<3; d++)
      {
         if ( Cr_Padded[d] < 0  ||  Cr_Padded[d] >= BoxNScale_Padded[d] )
            Aux_Error( ERROR_INFO, "incorrect Cr_Padded[%d]=(%d) --> out of range [%d ... %d) !!\n",
                       d, Cr_Padded[d], 0, BoxNScale_Padded[d] );
      }
#     endif

      PaddedCr1D = Mis_Idx3D2Idx1D( BoxNScale_Padded, Cr_Padded );
      LB_Idx     = LB_Corner2Index( lv, corner, CHECK_OFF );   // this number can be wrong for external patches
#     endif // #ifdef LOAD_BALANCE
      
      for (int s=0; s<26; s++ )  sibling[s] = -1;     // -1 <--> NO sibling

      flag    = false;
      fluid   = NULL;
#     if ( NPASSIVE > 0 )
      passive = NULL;
#     endif
#     ifdef GRAVITY
      pot     = NULL;
#     ifdef STORE_POT_GHOST
      pot_ext = NULL;
#     endif
#     endif

//    set the patch edge
      const int PScale = PS1*( 1<<(TOP_LEVEL-lv) );
      for (int d=0; d<3; d++)
      {
         EdgeL[d] = corner[d]*dh_min;
         EdgeR[d] = ( corner[d] + PScale )*dh_min;
      }

      for (int s=0; s<6; s++)    
      {
         flux        [s] = NULL;
#        if ( NPASSIVE > 0 )
         flux_passive[s] = NULL;
#        endif
#        ifdef GAMER_DEBUG
         flux_debug  [s] = NULL;
#        endif
      }

      if ( FluData )    hnew();
#     ifdef GRAVITY
      if ( PotData )    gnew();
#     endif

#     ifdef PARTICLE
      NPar         = 0;          // must be initialized as 0
      ParListSize  = 0;          // must be initialized as 0
      ParList      = NULL;

      NPar_Desc    = -1;         // -1 : indicating that NPar_Desc is not calculated yet
      ParList_Desc = NULL;

      for (int s=0; s<26; s++)
      {
         NPar_Escp   [s] = -1;   // -1 : indicating that NPar_Escp is not calculated yet
         ParList_Escp[s] = NULL;
      }
#     endif

   } // METHOD : patch_t



   //===================================================================================
   // Destructor  :  ~patch_t 
   // Description :  Destructor of the structure "patch_t"
   //
   // Note        :  Deallocate flux and data arrays
   //===================================================================================
   ~patch_t()
   {

      fdelete();
      hdelete();
#     ifdef GRAVITY
      gdelete();
#     endif

#     ifdef PARTICLE
      if ( ParList != NULL )
      {
         free( ParList );
         ParList = NULL;
      }

      if ( ParList_Desc != NULL )
      {
         delete [] ParList_Desc;
         ParList_Desc = NULL;
      }

      for (int s=0; s<26; s++)
      if ( ParList_Escp[s] != NULL )
      {
         free( ParList_Escp[s] );
         ParList_Escp[s] = NULL;
      }
#     endif

   } // METHOD : ~patch_t



   //===================================================================================
   // Method      :  fnew 
   // Description :  Allocate flux array in the given direction
   //
   // Note        :  Flux array is initialized as zero
   //
   // Parameter   :  SibID : Targeted ID of the flux array (0,1,2,3,4,5) <--> (-x,+x,-y,+y,-z,+z) 
   //===================================================================================
   void fnew( const int SibID )
   {

#     ifdef GAMER_DEBUG
      if ( SibID < 0  ||  SibID > 5 )
         Aux_Error( ERROR_INFO, "incorrect input in the member function \"fnew\" !!\n" );

      if ( flux[SibID] != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing flux array (sibling = %d) !!\n", SibID );

      if ( flux_debug[SibID] != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing flux_debug array (sibling = %d) !!\n", SibID );
#     endif

      flux[SibID] = new real [NFLUX+NPASSIVE][PATCH_SIZE][PATCH_SIZE];
            
      for(int v=0; v<NFLUX+NPASSIVE; v++)
      for(int m=0; m<PATCH_SIZE; m++)
      for(int n=0; n<PATCH_SIZE; n++)
         flux[SibID][v][m][n] = 0.0;

#     if ( NPASSIVE > 0 )
      flux_passive[SibID] = flux[SibID] + NFLUX;
#     endif

#     ifdef GAMER_DEBUG
      flux_debug[SibID] = new real [NFLUX+NPASSIVE][PATCH_SIZE][PATCH_SIZE];
            
      for(int v=0; v<NFLUX+NPASSIVE; v++)
      for(int m=0; m<PATCH_SIZE; m++)
      for(int n=0; n<PATCH_SIZE; n++)
         flux_debug[SibID][v][m][n] = 0.0;
#     endif

   } // METHOD : fnew



   //===================================================================================
   // Method      :  fdelete
   // Description :  Deallocate all flux arrays allocated previously 
   //===================================================================================
   void fdelete()
   {

      for (int s=0; s<6; s++)
      {
         if ( flux[s] != NULL )
         {
            delete [] flux[s];
            flux[s] = NULL;

#           if ( NPASSIVE > 0 )
            flux_passive[s] = NULL;
#           endif

#           ifdef GAMER_DEBUG
            delete [] flux_debug[s];
            flux_debug[s] = NULL;
#           endif
         }
      }

   } // METHOD : fdelete



   //===================================================================================
   // Method      :  hnew 
   // Description :  Allocate hydrodynamic array
   //
   // Note        :  An error message will be displayed if array has already been allocated
   //===================================================================================
   void hnew()
   {

#     ifdef GAMER_DEBUG
      if ( fluid != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing fluid array !!\n" );

#     if ( NPASSIVE > 0 )
      if ( passive != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing passive array !!\n" );
#     endif
#     endif // GAMER_DEBUG

//    fluid = new real [NCOMP+NPASSIVE][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
      fluid = new real [NCOMP+NPASSIVE+1][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
      fluid[0][0][0][0] = -1; 

#     if ( NPASSIVE > 0 )
      passive = fluid + NCOMP;
#     endif

   } // METHOD : hnew



   //===================================================================================
   // Method      :  hdelete
   // Description :  Deallocate hydrodynamic array
   //===================================================================================
   void hdelete()
   {

      if ( fluid != NULL )    
      {
         delete [] fluid;
         fluid = NULL;
      }

#     if ( NPASSIVE > 0 )
      passive = NULL;
#     endif

   } // METHOD : hdelete



#  ifdef GRAVITY
   //===================================================================================
   // Method      :  gnew 
   // Description :  Allocate potential array
   //
   // Note        :  An error message will be displayed if array has already been allocated
   //===================================================================================
   void gnew()
   {

#     ifdef GAMER_DEBUG
      if ( pot != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing pot array !!\n" );

#     ifdef STORE_POT_GHOST
      if ( pot_ext != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing pot_ext array !!\n" );
#     endif
#     endif

      pot     = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
#     ifdef STORE_POT_GHOST
      pot_ext = new real [GRA_NXT][GRA_NXT][GRA_NXT];
#     endif

   } // METHOD : gnew



   //===================================================================================
   // Method      :  gdelete
   // Description :  Deallocate potential array
   //===================================================================================
   void gdelete()
   {

      if ( pot != NULL )    
      {
         delete [] pot;
         pot = NULL;
      }

#     ifdef STORE_POT_GHOST
      if ( pot_ext != NULL )    
      {
         delete [] pot_ext;
         pot_ext = NULL;
      }
#     endif

   } // METHOD : gdelete
#  endif // #ifdef GRAVITY



#  ifdef PARTICLE
   //===================================================================================
   // Method      :  AddParticle
   // Description :  Add new particles into the particle list
   //
   // Note        :  1. Particles in the list must already lie inside the target patch
   //                   --> It will be checked in the debug mode
   //                2. If ParList already exists, new particles will be attached to the existing particle list 
   //                3. Also update the total number of particles at the target level
   //
   // Parameter   :  NNew     : Number of new particles to be added 
   //                NewList  : List storing the indices of new particles
   //                NPar_Lv  : Pointer to amr->Par->NPar_Lv[TargetLv]
   //                ParPos   : Particle position list             (debug only)
   //                NParTot  : Total number of existing particles (debug only)
   //                Comment  : Message for the debug mode         (debug only)
   //===================================================================================
#  ifdef DEBUG_PARTICLE
   void AddParticle( const int NNew, const long *NewList, long *NPar_Lv,
                     const real **ParPos, const long NParTot, const char *Comment )
#  else
   void AddParticle( const int NNew, const long *NewList, long *NPar_Lv )
#  endif
   {

//    check
#     ifdef DEBUG_PARTICLE
//    check 1: no strange settings
      if ( NewList == NULL  &&  NNew != 0 )  Aux_Error( ERROR_INFO, "\"%s\": NewList == NULL !!\n", Comment );
      if ( NNew < 0 )                        Aux_Error( ERROR_INFO, "\"%s\": NNew (%d) < 0 !!\n",   Comment, NNew );
      if ( NPar < 0 )                        Aux_Error( ERROR_INFO, "\"%s\": NPar (%d) < 0 !!\n",   Comment, NPar );

//    check 2: particle indices in NewList are NEW 
      for (int q=0; q<NPar; q++)
      for (int p=0; p<NNew; p++)
      {
         if ( ParList[q] == NewList[p] )  
            Aux_Error( ERROR_INFO, "\"%s\": repeated particle index (NewList[%d] = %ld) !!\n", Comment, p, NewList[p] );
      }

//    check 3: no duplicate particle index in NewList
      for (int q=0; q<NNew; q++)
      for (int p=0; p<NNew; p++)
      {
         if ( p != q  &&  NewList[q] == NewList[p] )  
            Aux_Error( ERROR_INFO, "\"%s\": duplicate particle index (NewList[%d] = NewList[%d] = %ld) !!\n", 
                       Comment, p, q, NewList[p] );
      }

//    check 4: whether the new particles indeed lie inside the target patch
      for (int p=0; p<NNew; p++)
      {
         const long ParID = NewList[p];
         
         for (int d=0; d<3; d++)
         {
            if ( ParPos[d][ParID] < EdgeL[d]  ||  ParPos[d][ParID] >= EdgeR[d] )
               Aux_Error( ERROR_INFO, "\"%s\": wrong home patch (L/R edge = %13.6e/%13.6e, pos[%d][%d] = %13.6e) !!\n",
                          Comment, EdgeL[d], EdgeR[d], d, ParID, ParPos[d][ParID] );
         }
      }

//    check 5: new particle IDs do not exceed the total number of particles
      for (int p=0; p<NNew; p++)
      {
         const long ParID = NewList[p];

         if ( ParID >= NParTot )    
            Aux_Error( ERROR_INFO, "\"%s\": Particle ID (%ld) >= Total number of active + inactive particles (%ld) !!\n", 
                       Comment, ParID, NParTot );
      }

      if ( NPar_Lv == NULL)   Aux_Error( ERROR_INFO, "NPar_Lv == NULL !!\n" );
#     endif // #ifdef DEBUG_PARTICLE


//    nothing to do if NNew == 0
      if ( NNew == 0 )  return;


//    allocate enough memory for the particle list array
      const int NPar_New = NPar + NNew;
      if ( NPar_New > ParListSize )   
      {
         ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*NPar_New );
         ParList     = (long*)realloc( ParList, ParListSize*sizeof(long) );
      }

//    record the new particle indices
      for (int p=0; p<NNew; p++)    ParList[ NPar + p ] = NewList[p];

//    update the particle number
      NPar = NPar_New;

//    update the particle number at the target level
      *NPar_Lv += NNew;

   } // METHOD : AddParticle



   //===================================================================================
   // Method      :  RemoveParticle
   // Description :  Remove particles from the particle list
   //
   // Note        :  1. The numbers stored in RemoveList are "array indices of ParList" 
   //                   instead of "particle indices stored in ParList"
   //                2. The numbers stored in RemoveList must be in ascending numerical order
   //                3. Also update the total number of particles at the target level
   //
   // Parameter   :  NRemove     : Number of particles to be removed 
   //                RemoveList  : List storing the array indices of "ParList" to be removed
   //                              **array indices, NOT particle indices**
   //                NPar_Lv     : Pointer to amr->Par->NPar_Lv[TargetLv]
   //                RemoveAll   : true --> remove all particle in this patch
   //===================================================================================
   void RemoveParticle( const int NRemove, const int *RemoveList, long *NPar_Lv, const bool RemoveAll )
   {

//    removing all particles is easy
      if ( RemoveAll )
      {
//       update the particle number at the target level
#        ifdef DEBUG_PARTICLE
         if ( NPar_Lv == NULL)   Aux_Error( ERROR_INFO, "NPar_Lv == NULL !!\n" );
#        endif

         *NPar_Lv -= NPar;

#        ifdef DEBUG_PARTICLE
         if ( *NPar_Lv < 0 )  Aux_Error( ERROR_INFO, "NPar_Lv = %ld < 0 !!\n", *NPar_Lv );
#        endif


//       remove all particles
         NPar        = 0;
         ParListSize = 0;

         if ( ParList != NULL )
         {
            free( ParList );
            ParList = NULL;
         }

         return;
      } // if ( RemoveAll )


//    check
#     ifdef DEBUG_PARTICLE
//    check 1: no strange settings
      if ( RemoveList == NULL  &&  NRemove != 0 )  Aux_Error( ERROR_INFO, "RemoveList == NULL !!\n" );
      if ( NRemove < 0 )                           Aux_Error( ERROR_INFO, "NRemove (%d) < 0 !!\n", NRemove );
      if ( NPar < 0 )                              Aux_Error( ERROR_INFO, "NPar (%d) < 0 !!\n", NPar );
      if ( NRemove > NPar )                        Aux_Error( ERROR_INFO, "NRemove (%d) > NPar (%d) !!\n", NRemove, NPar );

//    check 2: indices in RemoveList do not exceed the number of existing particles
      for (int p=0; p<NRemove; p++)
      {
         if ( RemoveList[p] >= NPar  ||  RemoveList[p] < 0 )
            Aux_Error( ERROR_INFO, "incorrect RemoveList[%d]=%d, (NPar=%d) !!\n", p, RemoveList[p], NPar );
      }

//    check 3: no duplicate index in RemoveList
      for (int q=0; q<NRemove; q++)
      for (int p=0; p<NRemove; p++)
      {
         if ( p != q  &&  RemoveList[q] == RemoveList[p] )  
            Aux_Error( ERROR_INFO, "duplicate index (RemoveList[%d] = RemoveList[%d] = %d) !!\n", 
                       p, q, RemoveList[p] );
      }

//    check 4: RemoveList is in ascending numerical order
      for (int p=1; p<NRemove; p++)
      {
         if ( RemoveList[p] <= RemoveList[p-1] )
            Aux_Error( ERROR_INFO, "RemoveList is NOT in ascending numerical order ([%d]=%d <= [%d]=%d) !!\n", 
                       p, RemoveList[p], p-1, RemoveList[p-1] );
      }

      if ( NPar_Lv == NULL)   Aux_Error( ERROR_INFO, "NPar_Lv == NULL !!\n" );
#     endif // #ifdef DEBUG_PARTICLE


//    nothing to do if NRemove == 0
      if ( NRemove == 0 )  return;


//    remove particles
      int LastP;
      for (int p=NRemove-1; p>=0; p--)
      {
         LastP = NPar - 1;

         if ( RemoveList[p] != LastP )   ParList[ RemoveList[p] ] = ParList[ LastP ];

         NPar --;
      }


//    release some memory if we waste too much
      const int ParListSize_Min = (int)floor( PARLIST_REDUCE_FACTOR*ParListSize );

      if ( NPar <= ParListSize_Min )
      {
         ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*NPar );
         ParList     = (long*)realloc( ParList, ParListSize*sizeof(long) );
      }


//    update the particle number at the target level
      *NPar_Lv -= NRemove;

#     ifdef DEBUG_PARTICLE
      if ( *NPar_Lv < 0 )  Aux_Error( ERROR_INFO, "NPar_Lv = %ld < 0 !!\n", *NPar_Lv );
#     endif

   } // METHOD : RemoveParticle
#  endif // #ifdef PARTICLE

   
}; // struct patch_t



#endif // #ifndef __PATCH_H__
