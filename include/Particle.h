#ifndef __PARTICLE_H__
#define __PARTICLE_H__



#ifndef PARTICLE
#  error : ERROR : PARTICLE is not defined !!
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_PARTICLE
#endif

// Factors used in "Add(One)Particle" and "Remove(One)Particle" to resize the allocated arrays
// PARLIST_GROWTH_FACTOR must >= 1.0; PARLIST_REDUCE_FACTOR must <= 1.0
#  define PARLIST_GROWTH_FACTOR     1.1
#  define PARLIST_REDUCE_FACTOR     0.8

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  Particle_t
// Description :  Data structure of particles
//
// Data Member :  ParListSize             : Size of the particle data arrays (must be >= NPar_AcPlusInac)
//                InactiveParListSize     : Size of the inactive particle list (InactiveParList)
//                NPar_Active_AllRank     : Total number of active particles summed up over all MPI ranks
//                NPar_AcPlusInac         : Total number of particles (active + inactive) in this MPI rank
//                NPar_Active             : Total number of active particles in this MPI rank
//                                          --> Inactive particles: particles removed from simulations because of
//                                           (1) lying outside the active region (mass == PAR_INACTIVE_OUTSIDE)
//                                           (2) being sent to other MPI ranks   (mass == PAR_INACTIVE_MPI)
//                NPar_Inactive           : Total number of inactive particles in this MPI rank
//                NPar_Lv                 : Total number of active particles at each level in this MPI rank
//                Init                    : Initialization methods (1/2/3 --> call function/restart/load from file)
//                ParICFormat             : Data format of the particle initialization file (1=[att][id], 2=[id][att])
//                ParICMass               : Assign this mass to all particles for Init=3
//                ParICType               : Assign this type to all particles for Init=3
//                Interp                  : Mass/acceleration/velocity interpolation scheme (NGP,CIC,TSC)
//                InterpTracer            : Mass/acceleration/velocity interpolation scheme for tracers (NGP,CIC,TSC)
//                Integ                   : Integration scheme (PAR_INTEG_EULER, PAR_INTEG_KDK)
//                IntegTracer             : Integration scheme for tracers (PAR_INTEG_EULER, PAR_INTEG_RK2)
//                ImproveAcc              : Improve force accuracy around the patch boundaries
//                                          (by using potential in the patch ghost zone instead of nearby patch
//                                          or interpolation)
//                PredictPos              : Predict particle position during mass assignment
//                TracerVelCorr           : Apply velocity correction term for tracer particles in regions where
//                                          the velocity gradient is large
//                RemoveCell              : remove particles RemoveCell-base-level-cells away from the boundary
//                                          (for non-periodic BC only)
//                GhostSize               : Number of ghost zones required for the interpolation scheme of massive particles
//                GhostSizeTracer         : Number of ghost zones required for the interpolation scheme of tracer  particles
//                AttributeFlt            : Pointer arrays to different particle floating-point attributes (Mass, Pos, Vel, ...)
//                AttributeInt            : Pointer arrays to different particle integer        attributes (Type)
//                InactiveParList         : List of inactive particle IDs
//                Mesh_Attr               : Pointer arrays to different mesh quantities mapped onto tracer particles
//                Mesh_Attr_Num           : Number of mesh quantities mapped onto tracer particles
//                Mesh_Attr_Idx           : Field indices of mesh quantities mapped onto tracer particles
//                Mesh_Attr_Label         : Field labels of mesh quantities mapped onto tracer particles
//                R2B_Real_NPatchTotal    : see R2B_Buff_NPatchTotal
//                R2B_Real_NPatchEachRank : see R2B_Buff_NPatchEachRank
//                R2B_Real_PIDList        : see R2B_Buff_PIDList
//                R2B_Buff_NPatchTotal    : Number of buffer patches to receive particles from the corresponding real patches
//                                          --> R2B : Real to Buffer
//                                          --> Mainly for the Poisson solver
//                                          --> For LOAD_BALANCE only
//                R2B_Buff_NPatchEachRank : Number of buffer patches to receive data from different MPI ranks
//                R2B_Buff_PIDList        : Buffer patch IDs list to receive particles from other ranks
//                                          --> Mainly for the Poisson solver
//                                          --> For LOAD_BALANCE only
//                                          --> It has the dimension [NLEVEL][2]
//                                              [lv][0/1] is for receiving particles at lv/lv-1 around real patches at lv
//                                              --> Mainly for the Poisson solver at lv (i.e., for calculating the
//                                                  total density field at lv)
//                                              --> More specific,
//                                                  [lv][0] is for receiving the particles of sibling-buffer patches at lv
//                                                  adjacent to real patches at lv
//                                                  [lv][1] is for receiving the particles of father-sibling-buffer patches at lv-1
//                                                  adjacent to real patches at lv
//                                              --> Both are constructed by calling Par_LB_RecordExchangeParticlePatchID( lv )
//                B2R_Real_NPatchTotal    : Similar to R2B_Real_NPatchTotal,    but for sending particles from buffer to real patches
//                B2R_Real_NPatchEachRank : Similar to R2B_Real_NPatchEachRank, ...
//                B2R_Real_PIDList        : Similar to R2B_Real_PIDList,        ...
//                B2R_Buff_NPatchTotal    : Similar to R2B_Buff_NPatchTotal,    ...
//                B2R_Buff_NPatchEachRank : Similar to R2B_Buff_NPatchEachRank, ...
//                B2R_Buff_PIDList        : Similar to R2B_Buff_PIDList,        ...
//                F2S_Send_NPatchTotal    : Similar to R2B_Real_NPatchTotal,    but for sending particles from father to son patches
//                F2S_Send_NPatchEachRank : Similar to R2B_Real_NPatchEachRank, ...
//                F2S_Send_PIDList        : Similar to R2B_Real_PIDList,        ...
//                F2S_Recv_NPatchTotal    : Similar to R2B_Buff_NPatchTotal,    ...
//                F2S_Recv_NPatchEachRank : Similar to R2B_Buff_NPatchEachRank, ...
//                F2S_Recv_PIDList        : Similar to R2B_Buff_PIDList,        ...
//                Mass                    : Particle mass
//                                          Mass < 0.0 --> this particle has been removed from simulations
//                                                     --> PAR_INACTIVE_OUTSIDE: fly outside the simulation box
//                                                         PAR_INACTIVE_MPI:     sent to other MPI ranks
//                Pos                     : Particle position
//                Vel                     : Particle velocity
//                Time                    : Particle physical time
//                Type                    : Particle type (e.g., tracer, generic, dark matter, star)
//                Acc                     : Particle acceleration (only when STORE_PAR_ACC is on)
//
// Method      :  Particle_t        : Constructor
//               ~Particle_t        : Destructor
//                InitRepo          : Initialize particle repository
//                AddOneParticle    : Add one new particle into the particle list
//                RemoveOneParticle : Remove one particle from the particle list
//-------------------------------------------------------------------------------------------------------
struct Particle_t
{

// data members
// ===================================================================================
   long          ParListSize;
   long          InactiveParListSize;
   long          NPar_Active_AllRank;
   long          NPar_AcPlusInac;
   long          NPar_Active;
   long          NPar_Inactive;
   long          NPar_Lv[NLEVEL];
   ParInit_t     Init;
   ParICFormat_t ParICFormat;
   double        ParICMass;
   int           ParICType;
   ParInteg_t    Integ;
   ParInterp_t   Interp;
   TracerInteg_t IntegTracer;
   ParInterp_t   InterpTracer;
   bool          ImproveAcc;
   bool          PredictPos;
   bool          TracerVelCorr;
   double        RemoveCell;
   int           GhostSize;
   int           GhostSizeTracer;
   real_par     *AttributeFlt[PAR_NATT_FLT_TOTAL];
   long_par     *AttributeInt[PAR_NATT_INT_TOTAL];
   long         *InactiveParList;
   real_par    **Mesh_Attr;
   int           Mesh_Attr_Num;
   long         *Mesh_Attr_Idx;
   char        **Mesh_Attr_Label;

#  ifdef LOAD_BALANCE
   int           R2B_Real_NPatchTotal   [NLEVEL][2];
   int          *R2B_Real_NPatchEachRank[NLEVEL][2];
   int          *R2B_Real_PIDList       [NLEVEL][2];
   int           R2B_Buff_NPatchTotal   [NLEVEL][2];
   int          *R2B_Buff_NPatchEachRank[NLEVEL][2];
   int          *R2B_Buff_PIDList       [NLEVEL][2];

   int           B2R_Real_NPatchTotal   [NLEVEL][2];
   int          *B2R_Real_NPatchEachRank[NLEVEL][2];
   int          *B2R_Real_PIDList       [NLEVEL][2];
   int           B2R_Buff_NPatchTotal   [NLEVEL][2];
   int          *B2R_Buff_NPatchEachRank[NLEVEL][2];
   int          *B2R_Buff_PIDList       [NLEVEL][2];

   int           F2S_Send_NPatchTotal   [NLEVEL];
   int          *F2S_Send_NPatchEachRank[NLEVEL];
   int          *F2S_Send_PIDList       [NLEVEL];
   int           F2S_Recv_NPatchTotal   [NLEVEL];
   int          *F2S_Recv_NPatchEachRank[NLEVEL];
   int          *F2S_Recv_PIDList       [NLEVEL];
#  endif // #ifdef LOAD_BALANCE

   real_par     *Mass;
   real_par     *PosX;
   real_par     *PosY;
   real_par     *PosZ;
   real_par     *VelX;
   real_par     *VelY;
   real_par     *VelZ;
   real_par     *Time;
#  ifdef STORE_PAR_ACC
   real_par     *AccX;
   real_par     *AccY;
   real_par     *AccZ;
#  endif
   long_par     *Type;


   //===================================================================================
   // Constructor :  Particle_t
   // Description :  Constructor of the structure "Particle_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Particle_t()
   {

      NPar_Active_AllRank = -1;
      NPar_AcPlusInac     = -1;
      Init                = PAR_INIT_NONE;
      ParICFormat         = PAR_IC_FORMAT_NONE;
      ParICMass           = -1.0;
      ParICType           = -1;
      Interp              = PAR_INTERP_NONE;
      InterpTracer        = PAR_INTERP_NONE;
      Integ               = PAR_INTEG_NONE;
      IntegTracer         = TRACER_INTEG_NONE;
      ImproveAcc          = true;
      PredictPos          = true;
      TracerVelCorr       = false;
      RemoveCell          = -999.9;
      GhostSize           = -1;
      GhostSizeTracer     = -1;

      for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv[lv] = 0;

      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   AttributeFlt[v] = NULL;

      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   AttributeInt[v] = NULL;

      InactiveParList = NULL;
      Mesh_Attr       = NULL;
      Mesh_Attr_Num   = 0;
      Mesh_Attr_Idx   = NULL;
      Mesh_Attr_Label = NULL;

#     ifdef LOAD_BALANCE
      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int t=0; t<2; t++) {
         R2B_Real_NPatchTotal   [lv][t] = 0;
         R2B_Real_NPatchEachRank[lv][t] = NULL;
         R2B_Real_PIDList       [lv][t] = NULL;
         R2B_Buff_NPatchTotal   [lv][t] = 0;
         R2B_Buff_NPatchEachRank[lv][t] = NULL;
         R2B_Buff_PIDList       [lv][t] = NULL;

         B2R_Real_NPatchTotal   [lv][t] = 0;
         B2R_Real_NPatchEachRank[lv][t] = NULL;
         B2R_Real_PIDList       [lv][t] = NULL;
         B2R_Buff_NPatchTotal   [lv][t] = 0;
         B2R_Buff_NPatchEachRank[lv][t] = NULL;
         B2R_Buff_PIDList       [lv][t] = NULL;
         }

         F2S_Send_NPatchTotal   [lv]    = 0;
         F2S_Send_NPatchEachRank[lv]    = NULL;
         F2S_Send_PIDList       [lv]    = NULL;
         F2S_Recv_NPatchTotal   [lv]    = 0;
         F2S_Recv_NPatchEachRank[lv]    = NULL;
         F2S_Recv_PIDList       [lv]    = NULL;
      } // for (int lv=0; lv<NLEVEL; lv++)
#     endif // #ifdef LOAD_BALANCE

      Mass = NULL;
      PosX = NULL;
      PosY = NULL;
      PosZ = NULL;
      VelX = NULL;
      VelY = NULL;
      VelZ = NULL;
      Time = NULL;
      Type = NULL;
#     ifdef STORE_PAR_ACC
      AccX = NULL;
      AccY = NULL;
      AccZ = NULL;
#     endif

   } // METHOD : Particle_t



   //===================================================================================
   // Destructor  :  ~Particle_t
   // Description :  Destructor of the structure "Particle_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Particle_t()
   {

      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
         if ( AttributeFlt[v] != NULL )   free( AttributeFlt[v] );

      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
         if ( AttributeInt[v] != NULL )   free( AttributeInt[v] );

      if ( InactiveParList != NULL )   free( InactiveParList );

      for (int v=0; v<Mesh_Attr_Num; v++)
      {
         if ( Mesh_Attr[v]       != NULL )   free( Mesh_Attr[v] );
         if ( Mesh_Attr_Label[v] != NULL )   free( Mesh_Attr_Label[v] );
      }

      if ( Mesh_Attr       != NULL )   free( Mesh_Attr );
      if ( Mesh_Attr_Idx   != NULL )   free( Mesh_Attr_Idx );
      if ( Mesh_Attr_Label != NULL )   free( Mesh_Attr_Label );

#     ifdef LOAD_BALANCE
      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int t=0; t<2; t++) {
         if ( R2B_Real_NPatchEachRank[lv][t] != NULL )   delete [] R2B_Real_NPatchEachRank[lv][t];
         if ( R2B_Buff_NPatchEachRank[lv][t] != NULL )   delete [] R2B_Buff_NPatchEachRank[lv][t];
         if ( R2B_Real_PIDList       [lv][t] != NULL )   delete [] R2B_Real_PIDList       [lv][t];
         if ( R2B_Buff_PIDList       [lv][t] != NULL )   free(     R2B_Buff_PIDList       [lv][t] );

         if ( B2R_Real_NPatchEachRank[lv][t] != NULL )   delete [] B2R_Real_NPatchEachRank[lv][t];
         if ( B2R_Buff_NPatchEachRank[lv][t] != NULL )   delete [] B2R_Buff_NPatchEachRank[lv][t];
         if ( B2R_Real_PIDList       [lv][t] != NULL )   delete [] B2R_Real_PIDList       [lv][t];
         if ( B2R_Buff_PIDList       [lv][t] != NULL )   free(     B2R_Buff_PIDList       [lv][t] );
         } // for (int t=0; t<2; t++)

         if ( F2S_Send_NPatchEachRank[lv] != NULL )      delete [] F2S_Send_NPatchEachRank[lv];
         if ( F2S_Recv_NPatchEachRank[lv] != NULL )      delete [] F2S_Recv_NPatchEachRank[lv];
         if ( F2S_Send_PIDList       [lv] != NULL )      free(     F2S_Send_PIDList       [lv] );
         if ( F2S_Recv_PIDList       [lv] != NULL )      delete [] F2S_Recv_PIDList       [lv];
      } // for (int lv=0; lv<NLEVE; lv++)
#     endif // LOAD_BALANCE

   } // METHOD : ~Particle_t



   //===================================================================================
   // Constructor :  InitRepo
   // Description :  Initialize particle repository
   //                --> All variables related to the number of particles
   //
   // Note        :  1. Initialize both "NPar_AcPlusInac, NPar_Active and ParListSize" as NPar_Input
   //                   --> Assuming no inactive particles (i.e., NPar_Inactive = 0)
   //                2. For LOAD_BALANCE, some lists recording the information for exchanging
   //                   particles between different ranks are also allocated here
   //
   // Parameter   :  NPar_Input : Total number of active particles
   //                NRank      : Total number of MPI ranks
   //===================================================================================
   void InitRepo( const long NPar_Input, const int NRank )
   {

//    check
      if ( NPar_Input < 0 )   Aux_Error( ERROR_INFO, "NPar_Input (%ld) < 0 !!\n", NPar_Input );

//    initialize variables related to the number of particles
      NPar_AcPlusInac     = NPar_Input;
      NPar_Active         = NPar_Input;                  // assuming all particles are active initially
      NPar_Inactive       = 0;
      ParListSize         = NPar_Input;                  // set ParListSize = NPar_AcPlusInac at the beginning
      InactiveParListSize = MAX( 1, ParListSize/100 );   // set arbitrarily (but must > 0)

//    allocate arrays (use malloc so that realloc can be used later to resize the array)
//    --> free memory first since other functions (e.g., LB_Init_LoadBalance()) will call InitRepo() again
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
      {
         if ( AttributeFlt[v] != NULL )   free( AttributeFlt[v] );
         AttributeFlt[v] = (real_par*)malloc( ParListSize*sizeof(real_par) );
      }
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
      {
         if ( AttributeInt[v] != NULL )   free( AttributeInt[v] );
         AttributeInt[v] = (long_par*)malloc( ParListSize*sizeof(long_par) );
      }

      if ( InactiveParList != NULL )   free( InactiveParList );
      InactiveParList = (long*)malloc( InactiveParListSize*sizeof(long) );

#     ifdef LOAD_BALANCE
      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int t=0; t<2; t++) {
            if ( R2B_Real_NPatchEachRank[lv][t] == NULL )   R2B_Real_NPatchEachRank[lv][t] = new int [NRank];
            if ( R2B_Buff_NPatchEachRank[lv][t] == NULL )   R2B_Buff_NPatchEachRank[lv][t] = new int [NRank];

            for (int r=0; r<NRank; r++)
            {
               R2B_Real_NPatchEachRank[lv][t][r] = 0;
               R2B_Buff_NPatchEachRank[lv][t][r] = 0;
            }


            if ( B2R_Real_NPatchEachRank[lv][t] == NULL )   B2R_Real_NPatchEachRank[lv][t] = new int [NRank];
            if ( B2R_Buff_NPatchEachRank[lv][t] == NULL )   B2R_Buff_NPatchEachRank[lv][t] = new int [NRank];

            for (int r=0; r<NRank; r++)
            {
               B2R_Real_NPatchEachRank[lv][t][r] = 0;
               B2R_Buff_NPatchEachRank[lv][t][r] = 0;
            }
         } // for (int t=0; t<2; t++)


         if ( F2S_Send_NPatchEachRank[lv] == NULL )   F2S_Send_NPatchEachRank[lv] = new int [NRank];
         if ( F2S_Recv_NPatchEachRank[lv] == NULL )   F2S_Recv_NPatchEachRank[lv] = new int [NRank];

         for (int r=0; r<NRank; r++)
         {
            F2S_Send_NPatchEachRank[lv][r] = 0;
            F2S_Recv_NPatchEachRank[lv][r] = 0;
         }
      } // for (int lv=0; lv<NLEVEL; lv++)
#     endif // #ifdef LOAD_BALANCE

//    set pointers
      Mass = AttributeFlt[PAR_MASS];
      PosX = AttributeFlt[PAR_POSX];
      PosY = AttributeFlt[PAR_POSY];
      PosZ = AttributeFlt[PAR_POSZ];
      VelX = AttributeFlt[PAR_VELX];
      VelY = AttributeFlt[PAR_VELY];
      VelZ = AttributeFlt[PAR_VELZ];
      Time = AttributeFlt[PAR_TIME];
#     ifdef STORE_PAR_ACC
      AccX = AttributeFlt[PAR_ACCX];
      AccY = AttributeFlt[PAR_ACCY];
      AccZ = AttributeFlt[PAR_ACCZ];
#     endif
      Type = AttributeInt[PAR_TYPE];

   } // METHOD : InitRepo



   //===================================================================================
   // Method      :  AddOneParticle
   // Description :  Add ONE new particle into the particle list
   //
   // Note        :  1. This function will modify several global variables
   //                   --> One must be careful when implementing OpenMP to avoid data race
   //                   --> For example, use the **critical** construct
   //                2. Please provide all necessary information of the new particle
   //                3. If there are inactive particles, their IDs will be reassigned
   //                   to the newly added particles
   //                4. Note that the global variable "AveDensity_Init" will NOT be recalculated
   //                   automatically here
   //
   // Parameter   :  NewAttFlt : Array storing the floating-point attributes of new particles
   //                NewAttInt : Array storing the integer        attributes of new particles
   //
   // Return      :  Index of the new particle (ParID)
   //===================================================================================
   long AddOneParticle( const real_par *NewAttFlt, const long_par *NewAttInt )
   {

//    check
#     ifdef DEBUG_PARTICLE
      if ( NPar_AcPlusInac < 0 ) Aux_Error( ERROR_INFO, "NPar_AcPlusInac (%ld) < 0 !!\n", NPar_AcPlusInac );

      if ( NewAttFlt == NULL )   Aux_Error( ERROR_INFO, "NewAttFlt == NULL !!\n" );

      if ( NewAttInt == NULL )   Aux_Error( ERROR_INFO, "NewAttInt == NULL !!\n" );

      if ( NewAttFlt[PAR_MASS] < (real_par)0.0 )
         Aux_Error( ERROR_INFO, "Adding an inactive particle (mass = %21.14e) !!\n", NewAttFlt[PAR_MASS] );

      if ( NewAttFlt[PAR_POSX] != NewAttFlt[PAR_POSX] ||
           NewAttFlt[PAR_POSY] != NewAttFlt[PAR_POSY] ||
           NewAttFlt[PAR_POSZ] != NewAttFlt[PAR_POSZ]   )
         Aux_Error( ERROR_INFO, "Adding a particle with strange position (%21.14e, %21.14e, %21.14e) !!\n",
                    NewAttFlt[PAR_POSX], NewAttFlt[PAR_POSY], NewAttFlt[PAR_POSZ] );

      if ( NewAttInt[PAR_TYPE] < (long_par)0  ||  NewAttInt[PAR_TYPE] >= (long_par)PAR_NTYPE )
         Aux_Error( ERROR_INFO, "Incorrect particle type (%ld) !!\n", (long)NewAttInt[PAR_TYPE] );
#     endif


//    1. determine the target particle ID
      long ParID;

//    1-1. reuse an inactive particle ID
      if ( NPar_Inactive > 0 )
      {
         ParID = InactiveParList[ NPar_Inactive-1 ];
         NPar_Inactive --;

#        ifdef DEBUG_PARTICLE
         if ( ParID < 0  ||  ParID >= NPar_AcPlusInac )
            Aux_Error( ERROR_INFO, "Incorrect ParID (%ld), NPar_AcPlusInac = %ld !!\n", ParID, NPar_AcPlusInac );
#        endif
      }

//    1-2. add a new particle ID
      else
      {
//       allocate enough memory for the particle variable array
         if ( NPar_AcPlusInac >= ParListSize )
         {
            ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*(ParListSize+1) );

            for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   AttributeFlt[v] = (real_par*)realloc( AttributeFlt[v], ParListSize*sizeof(real_par) );
            for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   AttributeInt[v] = (long_par*)realloc( AttributeInt[v], ParListSize*sizeof(long_par) );

            Mass = AttributeFlt[PAR_MASS];
            PosX = AttributeFlt[PAR_POSX];
            PosY = AttributeFlt[PAR_POSY];
            PosZ = AttributeFlt[PAR_POSZ];
            VelX = AttributeFlt[PAR_VELX];
            VelY = AttributeFlt[PAR_VELY];
            VelZ = AttributeFlt[PAR_VELZ];
            Time = AttributeFlt[PAR_TIME];
#           ifdef STORE_PAR_ACC
            AccX = AttributeFlt[PAR_ACCX];
            AccY = AttributeFlt[PAR_ACCY];
            AccZ = AttributeFlt[PAR_ACCZ];
#           endif
            Type = AttributeInt[PAR_TYPE];
         }

         ParID = NPar_AcPlusInac;
         NPar_AcPlusInac ++;
      } // if ( NPar_Inactive > 0 ) ... else ...


//    2. record the data of new particles
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   AttributeFlt[v][ParID] = NewAttFlt[v];
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   AttributeInt[v][ParID] = NewAttInt[v];


//    3. update the total number of active particles (assuming all new particles are active)
      NPar_Active ++;


//    4. return the new particle index
      return ParID;

   } // METHOD : AddOneParticle



   //===================================================================================
   // Method      :  RemoveOneParticle
   // Description :  Remove ONE particle from the particle list
   //
   // Note        :  1. This function will modify NPar_Active and NPar_Inactive
   //                   --> Since these are global variables, one must be careful for the OpenMP
   //                       implementation to avoid data race
   //                   --> For example, use the **critical** construct
   //                2. Particles being removed will NOT be actually deleted in the memory.
   //                   Instead, their masses are assigned to Marker, which must be negative
   //                   in order to be distinguishable from active particles.
   //                   --> IDs of all removed particles will be recorded in InactiveParList[]
   //                   --> These IDs will be reassigned to new particles added later on when
   //                       calling AddOneParticle()
   //                3. Note that the global variable "AveDensity_Init" will NOT be recalculated
   //                   automatically here
   //
   // Parameter   :  ParID  : Particle ID to be removed
   //                Marker : Value assigned to the mass of the particle being removed
   //                         (PAR_INACTIVE_OUTSIDE or PAR_INACTIVE_MPI)
   //
   // Return      :  None
   //===================================================================================
   void RemoveOneParticle( const long ParID, const real_par Marker )
   {

//    check
#     ifdef DEBUG_PARTICLE
      if ( ParID < 0  ||  ParID >= NPar_AcPlusInac )
         Aux_Error( ERROR_INFO, "Wrong ParID (%ld) !!\n", ParID );

      if ( Marker != PAR_INACTIVE_OUTSIDE  &&  Marker != PAR_INACTIVE_MPI )
         Aux_Error( ERROR_INFO, "Unsupported Marker (%14.7e) !!\n", Marker );
#     endif


//    1. allocate enough memory for InactiveParList
      if ( NPar_Inactive >= InactiveParListSize )
      {
         InactiveParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*(InactiveParListSize+1) );

         InactiveParList = (long*)realloc( InactiveParList, InactiveParListSize*sizeof(long) );
      }


//    2. record the particle ID to be removed
      InactiveParList[NPar_Inactive] = ParID;


//    3. remove the target particle
      Mass[ParID] = Marker;

      NPar_Active   --;
      NPar_Inactive ++;

#     ifdef DEBUG_PARTICLE
      if ( NPar_Active + NPar_Inactive != NPar_AcPlusInac )
         Aux_Error( ERROR_INFO, "NPar_Active (%ld) + NPar_Inactive (%ld) != NPar_AcPlusInac (%ld) !!\n",
                    NPar_Active, NPar_Inactive, NPar_AcPlusInac );
#     endif

   } // METHOD : RemoveOneParticle


}; // struct Particle_t



#endif // #ifndef __PARTICLE_H__
