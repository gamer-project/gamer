#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

static real GetMaxVelocity( const int lv, bool excludeWave);
static void GetMaxLevel(patch_t *patch, int lv, int& maxlv);

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Velocity
// Description :  Estimate the evolution time-step by via the CFL condition from the Hamilton-Jacobi equation
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. dt = 0.5 * dx * m/hbar / sum_i |v_i|
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Velocity( const int lv )
{

   real   MaxdS_dx, NoExclusion;
   double dt, dt2;


// get the velocity dS/dx / dx as first derivative of phase
   MaxdS_dx     = GetMaxVelocity( lv, true);
   //NoExclusion  = GetMaxVelocity( lv, false);


   MaxdS_dx += (MaxdS_dx == 0) ? 1e-8 : 0;
   //NoExclusion += (NoExclusion == 0) ? 1e-8 : 0;

// get the time-step
   dt  = 1 / MaxdS_dx * 0.5 * ELBDM_ETA * DT__VELOCITY;
   //dt2 = 1 / NoExclusion * 0.5 * ELBDM_ETA * DT__VELOCITY;

   //if (MPI_Rank == 0)
   //printf("V with Excl: %f No Excl: %f dt w. Excl: %f dt n. Excl. %f\n", MaxdS_dx, NoExclusion, dt, dt2);

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Velocity




bool hasWaveCounterpart(int I, int J, int K, long GID0, LB_GlobalPatch* tree) {
   int lv = tree[GID0].level; 

   //printf("Checking I %d J %d K %d on lv %d with GID %ld\n", I, J, K, lv, GID0);
// convert coordinates of cell to global integer coordinate system 
   int coordinates[3] = {I, J, K};
   for ( int l = 0; l < 3; ++l ) 
      coordinates[l] = tree[GID0].corner[l] + coordinates[l] * amr->scale[lv];
   
   //printf("parent corner x %d y %d z %d with side length %d\n", 
   //tree[GID0].corner[0], tree[GID0].corner[1] , tree[GID0].corner[2], 
   //PS1 * amr->scale[lv]);

// during first iteration, we can cover the eight patches in the patch group
// imagine that patches in patch group are children of a common father even at level 0
// we set currentLv = 0 in loop and iterate over patch group
// once we find the patch that the coordinate is in, we iterate over that patch's children

   int  currentLv  = lv - 1;
   bool isInside   = true; 

   long currentGID = GID0;
   long sonGID     = GID0; 

// traverse the tree until we get to leave nodes
   while ( sonGID != -1 ) {
      currentLv += 1;

//    loop over patches in patch block or sons and check which patch the cell {K, J, I} belongs to
      for (int GID=sonGID; GID<sonGID+8; GID++)
      {
         isInside = true; 
//       check whether cell {I, J, K} is within g,iven patch
         for ( int l = 0; l < 3; ++l ) {
            if (coordinates[l] < tree[GID].corner[l] || coordinates[l] >=  tree[GID].corner[l] + PS1 * amr->scale[currentLv])
               isInside = false;
         }

//       if it is within the given patch, go to the patch's son
         if ( isInside ) { 
            currentGID = GID;
            sonGID     = tree[GID].son;
            break;
         } else if ( !isInside && (GID == (sonGID + 7)) ) {
            Aux_Error(ERROR_INFO, "Coordinates in hasWaveCounter part not in correct patch block!\n");
            return false;
         }
      }
   }

   return amr->use_wave_flag[currentLv]; 
}


//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxVelocity
// Description :  Evaluate the maximum first spatial derivative of phase for the time-step estimation
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Velocity()
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MaxdS_dx
//-------------------------------------------------------------------------------------------------------
real GetMaxVelocity( const int lv, bool excludeWave )
{
   // Maximum velocity for calculation of time step zero for wave patches
   // since we are only concerned with a velocity dependent time step criterion for the fluid solver
   if ( amr->use_wave_flag[lv] )
      return 0;


   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;
   const int  NGhost            = 1;               // number of ghost zones to calculate dS_dx
   const int  Size_Flu          = PS2 + 2*NGhost;  // size of the array Flu_Array
   const int  NPG               = 1;               // number of patch groups (must be ONE here)
   const int  NCOMP1            = 1;               // we only need to retrieve phase for velocity criterion
   real (*Flu_Array)[NCOMP1][Size_Flu][Size_Flu][Size_Flu] = NULL;
   real (*Wave_Patch_Array)[NCOMP1][Size_Flu][Size_Flu][Size_Flu] = NULL;

   real dS_dx, MaxdS_dx, _dh, _dh2, GradS[3];
   int im, ip, jm, jp, km, kp, I, J, K;

   MaxdS_dx = 0.0;
   _dh      = (real)1.0/amr->dh[lv];
   _dh2     = (real)0.5*_dh;

   LB_PatchCount pc;

// construct global tree structure 
   LB_GlobalPatch* global_tree = LB_AllgatherTree(pc);

/*
   for(int i = 0; i < MPI_NRank; i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (i == MPI_Rank) {
         printf("Rank %d \n\n", MPI_Rank);
         printf("Adress of pointer %p\n", global_tree);
      
         printf("\nGlobal Tree\n");
      
         printf("NPatchAllLv: %ld\n", NPatchAllLv);
      
         for (int i = 0; i < NLEVEL; ++i) {
            printf("lv: %d GID_Offset: %ld\n", i, GID_Offset[i]);
         }   
      
      
         printf("\nThe GIDs\n");
      
         for (long i = 0; i < NPatchAllLv; ++i ) {
            printf("lv %d GID %d father %d son %d\n", global_tree[i].level, i, global_tree[i].father, global_tree[i].son);
            }
         }
   }
*/
#  pragma omp parallel private( Flu_Array, dS_dx, GradS, \
                                im, ip, jm, jp, km, kp, I, J, K)
   {
      Flu_Array = new real [NPG][NCOMP1][Size_Flu][Size_Flu][Size_Flu];

//    loop over all patches
#     pragma omp for reduction( max:MaxdS_dx ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)
      {
         
//       prepare phase with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][0][0][0][0], NULL, NGhost, NPG, &PID0, _PHAS, _NONE,
                            INT_MINMOD1D, INT_NONE, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


         long GID0 = PID0 + pc.GID_Offset[lv]; 

//       evaluate dS_dt
         for (int k=NGhost; k<Size_Flu-NGhost; k++)    {  km = k - 1;    kp = k + 1;   K = k - NGhost;
         for (int j=NGhost; j<Size_Flu-NGhost; j++)    {  jm = j - 1;    jp = j + 1;   J = j - NGhost;
         for (int i=NGhost; i<Size_Flu-NGhost; i++)    {  im = i - 1;    ip = i + 1;   I = i - NGhost;

//       skip velocities of cells that have a wave counterpart on the refined levels

         bool calculateVelocity = true;
         if (excludeWave) calculateVelocity = !hasWaveCounterpart(I, J, K, GID0, global_tree );

//       check whether leave node corresponding to cell uses phase scheme
         if ( calculateVelocity ) {

            GradS[0]   = _dh2 * ( Flu_Array[0][0][k ][j ][ip] - Flu_Array[0][0][k ][j ][im] );
            GradS[1]   = _dh2 * ( Flu_Array[0][0][k ][jp][i ] - Flu_Array[0][0][k ][jm][i ] );
            GradS[2]   = _dh2 * ( Flu_Array[0][0][kp][j ][i ] - Flu_Array[0][0][km][j ][i ] );

            dS_dx      = _dh * (FABS(GradS[0]) + FABS(GradS[1]) + FABS(GradS[2]));
            MaxdS_dx   = MAX( MaxdS_dx, dS_dx );

         }

         }}} // i,j,k
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)

      delete [] Flu_Array;
   } // OpenMP parallel region

// clean up global tree
   delete [] global_tree;


// get the maximum potential in all ranks
   real MaxdS_dx_AllRank;
#  ifdef FLOAT8
   MPI_Allreduce( &MaxdS_dx, &MaxdS_dx_AllRank, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#  else
   MPI_Allreduce( &MaxdS_dx, &MaxdS_dx_AllRank, 1, MPI_FLOAT,  MPI_MAX, MPI_COMM_WORLD );
#  endif


// check
   if ( MaxdS_dx_AllRank == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MaxdS_dx == 0.0 at lv %d !!\n", lv );


   return MaxdS_dx_AllRank;

} // FUNCTION : GetMaxVelocity

#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

/*
# else 
#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

static real GetMaxVelocity( const int lv );
static void GetMaxLevel(patch_t *patch, int lv, int& maxlv);

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Velocity
// Description :  Estimate the evolution time-step by via the CFL condition from the Hamilton-Jacobi equation
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. dt = 0.5 * dx * m/hbar / sum_i |v_i|
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Velocity( const int lv )
{

   real   MaxdS_dx;
   double dt;


// get the velocity dS/dx / dx as first derivative of phase
   MaxdS_dx  = GetMaxVelocity( lv );
   if (MPI_Rank == 0)
   printf("Max velocity: %f\n", MaxdS_dx);
   MaxdS_dx += (MaxdS_dx == 0) ? 1e-8 : 0;

// get the time-step
   dt = 1 / MaxdS_dx * 0.5 * ELBDM_ETA * DT__VELOCITY;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Velocity



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxVelocity
// Description :  Evaluate the maximum first spatial derivative of phase for the time-step estimation
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Velocity()
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MaxdS_dx
//-------------------------------------------------------------------------------------------------------
real GetMaxVelocity( const int lv )
{
   // Maximum velocity for calculation of time step zero for wave patches
   // since we are only concerned with a velocity dependent time step criterion for the fluid solver
   if ( amr->use_wave_flag[lv] )
      return 0;


   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;
   const int  NGhost            = 1;               // number of ghost zones to calculate dS_dx
   const int  Size_Flu          = PS2 + 2*NGhost;  // size of the array Flu_Array
   const int  NPG               = 1;               // number of patch groups (must be ONE here)
   const int  NCOMP1            = 1;               // we only need to retrieve phase for velocity criterion
   real (*Flu_Array)[NCOMP1][Size_Flu][Size_Flu][Size_Flu] = NULL;
   real (*Wave_Patch_Array)[NCOMP1][Size_Flu][Size_Flu][Size_Flu] = NULL;

   real dS_dx, MaxdS_dx, _dh, _dh2, GradS[3];
   int im, ip, jm, jp, km, kp, I, J, K;
   int maxlv = 0;
   int maxFluidLevel = 0;
   for (int i = 1; i < NLEVEL; ++i) {
      if ( !amr->use_wave_flag[i] )
         ++maxFluidLevel;
      else 
         break;
   }

   bool skipPatchgroup;


   MaxdS_dx = 0.0;
   _dh      = (real)1.0/amr->dh[lv];
   _dh2     = (real)0.5*_dh;



#  pragma omp parallel private( Flu_Array, dS_dx, GradS, \
                                im, ip, jm, jp, km, kp, I, J, K, maxlv, skipPatchgroup )
   {
      Flu_Array = new real [NPG][NCOMP1][Size_Flu][Size_Flu][Size_Flu];

//    loop over all patches
#     pragma omp for reduction( max:MaxdS_dx ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)
      {
         
//       prepare phase with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][0][0][0][0], NULL, NGhost, NPG, &PID0, _PHAS, _NONE,
                            INT_MINMOD1D, INT_NONE, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );



//       evaluate dS_dt
         for (int k=NGhost; k<Size_Flu-NGhost; k++)    {  km = k - 1;    kp = k + 1;   K = k - NGhost;
         for (int j=NGhost; j<Size_Flu-NGhost; j++)    {  jm = j - 1;    jp = j + 1;   J = j - NGhost;
         for (int i=NGhost; i<Size_Flu-NGhost; i++)    {  im = i - 1;    ip = i + 1;   I = i - NGhost;

//       skip velocities of cells that have a wave counterpart on the refined levels

//       convert coordinates of cell to global integer coordinate system 
         int coordinates[3] = {I, J, K};

         for ( int l = 0; l < 3; ++l ) {
            coordinates[l] *= amr->scale[lv];
            coordinates[l] += amr->patch[0][lv][PID0]->corner[l];
         }

         int curLv  = lv - 1, curPID = PID0, sonPID = PID0; 
         bool inside = true; 

//       traverse the tree until we get to leave nodes
         while ( sonPID != -1 ) {
            curPID = sonPID; 
            curLv += 1;

//          loop over patches in patch block or sons and check which patch the cell {K, J, I} belongs to
            for (int PID=curPID; PID<curPID+8; PID++)
            {
               inside = true; 

//             check whether cell {I, J, K} is within given patch
               for ( int l = 0; l < 3; ++l ) {
                  if (coordinates[l] < amr->patch[0][curLv][PID]->corner[l] || coordinates[l] >=  amr->patch[0][curLv][PID]->corner[l] + PS1 * amr->scale[curLv])
                     inside = false;
               }

//             if it is within the given patch, go to the patch's son
               if ( inside ) { 
                  curPID = PID;
                  sonPID = amr->patch[0][curLv][curPID]->son;
                  break;
               }
            }
         }

//       check whether leave node corresponding to cell uses phase scheme
         if ( !amr->use_wave_flag[curLv] ) {

            GradS[0]   = _dh2 * ( Flu_Array[0][0][k ][j ][ip] - Flu_Array[0][0][k ][j ][im] );
            GradS[1]   = _dh2 * ( Flu_Array[0][0][k ][jp][i ] - Flu_Array[0][0][k ][jm][i ] );
            GradS[2]   = _dh2 * ( Flu_Array[0][0][kp][j ][i ] - Flu_Array[0][0][km][j ][i ] );

            dS_dx      = _dh * (FABS(GradS[0]) + FABS(GradS[1]) + FABS(GradS[2]));
            MaxdS_dx   = MAX( MaxdS_dx, dS_dx );

         }

         }}} // i,j,k
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)

      delete [] Flu_Array;
   } // OpenMP parallel region


// get the maximum potential in all ranks
   real MaxdS_dx_AllRank;
#  ifdef FLOAT8
   MPI_Allreduce( &MaxdS_dx, &MaxdS_dx_AllRank, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#  else
   MPI_Allreduce( &MaxdS_dx, &MaxdS_dx_AllRank, 1, MPI_FLOAT,  MPI_MAX, MPI_COMM_WORLD );
#  endif


// check
   if ( MaxdS_dx_AllRank == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MaxdS_dx == 0.0 at lv %d !!\n", lv );


   return MaxdS_dx_AllRank;

} // FUNCTION : GetMaxVelocity


void GetMaxLevel(patch_t *patch, int lv, int& maxlv) {
   int SonPID = patch->son;

   if ( SonPID != -1 )
   {
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         GetMaxLevel ( amr->patch[0][lv+1][SonPID+LocalID], lv + 1, maxlv);
      }
   } 
   
   maxlv = MAX( maxlv, lv );      
}

#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
#endif */