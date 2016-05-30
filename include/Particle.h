#ifndef __PARTICLE_H__
#define __PARTICLE_H__



#ifndef PARTICLE
#  error : ERROR : PARTICLE is not defined !!
#endif

//#ifdef GAMER_DEBUG
#  define DEBUG_PARTICLE
//#endif

// Factors used in "AddParticle" and "RemoveParticle" to determine the ratio between ParListSize and NPar
// PARLIST_GROWTH_FACTOR must >= 1.0; PARLIST_REDUCE_FACTOR must <= 1.0
#  define PARLIST_GROWTH_FACTOR     1.1
#  define PARLIST_REDUCE_FACTOR     0.8

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  Particle_t 
// Description :  Data structure of particles
//
// Data Member :  ParListSize : Size of the particle data arrays (must be >= NPar) 
//                NPar        : Total number of particles (active + inactive)
//                NPar_Active : Total number of active particles 
//                              (inactive: particles removed from simulation likely because they lie outside the active region)
//                NPar_Lv     : Total number of active particles at each level
//                Init        : Initialization methods (1/2/3 --> call function/restart/load from file)
//                Interp      : Mass/acceleration interpolation scheme (NGP,CIC,TSC)
//                Integ       : Integration scheme (PAR_INTEG_EULER, PAR_INTEG_KDK)
//                ImproveAcc  : Improve force accuracy around the patch boundaries
//                              (by using potential in the patch ghost zone instead of nearby patch or interpolation)
//                PredictPos  : Predict particle position during mass assignment
//                RemoveCell  : remove particles RemoveCell-base-level-cells away from the boundary (for non-periodic BC only
//                GhostSize   : Number of ghost zones required for interpolation scheme
//                ParVar      : Pointer arrays to different particle variables (Mass, Pos, Vel, ...)
//                Passive     : Pointer arrays to different passive variables (e.g., metalicity)
//                Mass        : Particle mass
//                              Mass < 0.0 --> this particle has been removed from the simulation box
//                                         --> -1.0: fly outside the simulation box
//                Pos         : Particle position
//                Vel         : Particle velocity
//                Time        : Particle physical time
//
// Method      :  Particle_t     : Constructor 
//               ~Particle_t     : Destructor
//                InitVar        : Initialize arrays and some variables
//                AddOneParticle : Add new particles into the particle list
//-------------------------------------------------------------------------------------------------------
struct Particle_t
{

// data members
// ===================================================================================
   long        ParListSize;
   long        NPar;
   long        NPar_Active;
   long        NPar_Lv[NLEVEL];
   ParInit_t   Init;
   ParInterp_t Interp;
   ParInteg_t  Integ;
   bool        ImproveAcc;
   bool        PredictPos;
   double      RemoveCell;
   int         GhostSize;
   real       *ParVar [NPAR_VAR    ];
   real       *Passive[NPAR_PASSIVE];

   real       *Mass;
   real       *PosX;
   real       *PosY;
   real       *PosZ;
   real       *VelX;
   real       *VelY;
   real       *VelZ;
   real       *Time;


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

      NPar       = -1;
      Init       = PAR_INIT_NONE;
      Interp     = PAR_INTERP_NONE;
      Integ      = PAR_INTEG_NONE;
      ImproveAcc = true;
      PredictPos = true;
      RemoveCell = -999.9;

      for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv[lv] = 0;

      for (int v=0; v<NPAR_VAR;     v++)  ParVar [v] = NULL;
      for (int v=0; v<NPAR_PASSIVE; v++)  Passive[v] = NULL;

      Mass = NULL;
      PosX = NULL;
      PosY = NULL;
      PosZ = NULL;
      VelX = NULL;
      VelY = NULL;
      VelZ = NULL;
      Time = NULL;

   } // METHOD : Particle_t



   //===================================================================================
   // Destructor  :  ~Particle_t 
   // Description :  Destructor of the structure "Particle_t"
   //
   // Note        :  Free memory 
   //===================================================================================
   ~Particle_t()
   {

      for (int v=0; v<NPAR_VAR; v++)
      {
         if ( ParVar[v] != NULL )   
         {
            free( ParVar[v] );
            ParVar[v] = NULL;
         }
      }

      for (int v=0; v<NPAR_PASSIVE; v++)
      {
         if ( Passive[v] != NULL )   
         {
            free( Passive[v] );
            Passive[v] = NULL;
         }
      }

   } // METHOD : ~Particle_t



   //===================================================================================
   // Constructor :  InitVar
   // Description :  Initialize the particle attribute arrays and some other variables
   //
   // Note        :  1. NPar must be set properly (>=0) in advance 
   //                2. Initialize both "NPar_Active" and "NPar_Active" as NPar
   //                3. Set "GhostSize" ("Interp" must be set in advance)
   //
   // Parameter   :  None
   //===================================================================================
   void InitVar()
   {

//    check
      if ( NPar < 0 )                     Aux_Error( ERROR_INFO, "NPar (%d) < 0 !!\n", NPar );
      if ( Interp == PAR_INTERP_NONE )    Aux_Error( ERROR_INFO, "Interp == NONE !!\n" );

//    initialize NPar_Active and ParListSize
      NPar_Active = NPar;  // assuming all particles are active initially
      ParListSize = NPar;  // set ParListSize = NPar at the beginning

//    set the number of ghost zones for the interpolation scheme
      switch ( Interp )
      {
         case ( PAR_INTERP_NGP ): GhostSize = 0;   break;
         case ( PAR_INTERP_CIC ): GhostSize = 1;   break;
         case ( PAR_INTERP_TSC ): GhostSize = 1;   break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }

//    allocate arrays (use malloc so that realloc can be used later to resize the array)
      for (int v=0; v<NPAR_VAR;     v++)  ParVar [v] = (real*)malloc( ParListSize*sizeof(real) );
      for (int v=0; v<NPAR_PASSIVE; v++)  Passive[v] = (real*)malloc( ParListSize*sizeof(real) );

      Mass = ParVar[PAR_MASS];
      PosX = ParVar[PAR_POSX];
      PosY = ParVar[PAR_POSY];
      PosZ = ParVar[PAR_POSZ];
      VelX = ParVar[PAR_VELX];
      VelY = ParVar[PAR_VELY];
      VelZ = ParVar[PAR_VELZ];
      Time = ParVar[PAR_TIME];

   } // METHOD : InitVar



   //===================================================================================
   // Method      :  AddOneParticle
   // Description :  Add ONE new particle into the particle list
   //
   // Note        :  Please provide all necessary information of the new particle
   //
   // Parameter   :  NewVar      : Array storing the         variables of new particles
   //                NewPassive  : Array storing the passive variables of new particles
   //===================================================================================
   void AddOneParticle( const real *NewVar, const real *NewPassive )
   {

      const int NNew = 1;

//    check
#     ifdef DEBUG_PARTICLE
      if ( NPar < 0 )            Aux_Error( ERROR_INFO, "NPar (%d) < 0 !!\n", NPar );
      if ( NewVar     == NULL )  Aux_Error( ERROR_INFO, "NewVar == NULL !!\n" );
#     if ( NPAR_PASSIVE > 0 )
      if ( NewPassive == NULL )  Aux_Error( ERROR_INFO, "NewVar == NULL !!\n" );
#     endif
#     endif


//    1. nothing to do if NNew == 0
//    if ( NNew == 0 )  return;


//    2. allocate enough memory for the particle variable array
      const int NPar_New = NPar + NNew;

      if ( NPar_New > ParListSize )   
      {
         ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*NPar_New );

         for (int v=0; v<NPAR_VAR;     v++)  ParVar [v] = (real*)realloc( ParVar [v], ParListSize*sizeof(real) );
         for (int v=0; v<NPAR_PASSIVE; v++)  Passive[v] = (real*)realloc( Passive[v], ParListSize*sizeof(real) );

         Mass = ParVar[PAR_MASS];
         PosX = ParVar[PAR_POSX];
         PosY = ParVar[PAR_POSY];
         PosZ = ParVar[PAR_POSZ];
         VelX = ParVar[PAR_VELX];
         VelY = ParVar[PAR_VELY];
         VelZ = ParVar[PAR_VELZ];
         Time = ParVar[PAR_TIME];
      }


//    3. record the data of new particles
      long ParID;

      for (int p=0; p<NNew; p++)    
      {
         ParID = NPar + p;

         /*
         for (int v=0; v<NPAR_VAR;     v++)  ParVar [v][ParID] = NewVar    [p][v];
         for (int v=0; v<NPAR_PASSIVE; v++)  Passive[v][ParID] = NewPassive[p][v];
         */
         for (int v=0; v<NPAR_VAR;     v++)  ParVar [v][ParID] = NewVar    [v];
         for (int v=0; v<NPAR_PASSIVE; v++)  Passive[v][ParID] = NewPassive[v];
      }


//    4. update the particle number (assuming all new particles are active)
      NPar        += NNew;
      NPar_Active += NNew;

   } // METHOD : AddOneParticle

   
}; // struct Particle_t 



#endif // #ifndef __PARTICLE_H__
