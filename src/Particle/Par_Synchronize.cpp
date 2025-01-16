#include "GAMER.h"

#ifdef PARTICLE


// static global variables
static double     CurrentSyncTime   = -1.0;
static long       Backup_NPar       = -1;
static long      *Backup_ParID      = NULL;
static real_par (*Backup_ParAtt)[7] = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Synchronize
// Description :  Synchronize all particles to a given time
//
// Note        :  1. For SyncOption = PAR_SYNC_TEMP, particle position and velocity are only temporarily
//                   synchronized and can be restored by calling Par_Synchronize_Restore()
//                2. For SyncOption = PAR_SYNC_FORCE, particle position and velocity are permanently
//                   synchronized and cannot be restored
//                3. One can use the static global variable "CurrentSyncTime" to check whether particles
//                   have already been synchronized (avoid duplicate synchronization when, for example,
//                   outputting data and calculating total energy at the same step)
//                4. Must work with STORE_PAR_ACC
//                5. Particles may cross patch boundaries after synchronization
//                   --> One may need to call Par_PassParticle2Sibling() and Par_PassParticle2Son_MultiPatch() to
//                       properly transfer particles between patches, especially for SyncOption = PAR_SYNC_FORCE
//                6. Currently it's only invoked by Flu_CorrAfterAllSync()
//
// Parameter   :  SyncTime   : Target synchronization time
//                SyncOption : PAR_SYNC_NONE / PAR_SYNC_TEMP / PAR_SYNC_FORCE
//
// Return      :  0/1/2 --> success/PAR_SYNC_NONE/routine has already been called previously
//-------------------------------------------------------------------------------------------------------
int Par_Synchronize( const double SyncTime, const ParSync_t SyncOption )
{

// check
   if ( SyncOption == PAR_SYNC_NONE )  return 1;

   if ( SyncOption != PAR_SYNC_TEMP  &&  SyncOption != PAR_SYNC_FORCE )
      Aux_Error( ERROR_INFO, "unsupported SyncOption = %d !!\n", SyncOption );

#  if ( defined MASSIVE_PARTICLES  &&  !defined STORE_PAR_ACC )
   if ( SyncOption == PAR_SYNC_TEMP  ||  SyncOption == PAR_SYNC_FORCE )
      Aux_Error( ERROR_INFO, "please turn on STORE_PAR_ACC in the makefile for particle synchronization !!\n" );
#  endif


// nothing to do if this routine has been called previously with the same synchronization time
   if ( SyncTime == CurrentSyncTime )  return 2;


// allocate the backup array
   long MemUnit, MemSize;

   if ( SyncOption == PAR_SYNC_TEMP )
   {
      if ( Backup_NPar >= 0 )    Aux_Error( ERROR_INFO, "backup arrays have already been allocated !!\n" );

      MemUnit       = MAX( 1, amr->Par->NPar_Active/100 );  // set arbitrarily (but must > 0)
      MemSize       = MemUnit;
      Backup_NPar   = 0;
      Backup_ParID  = ( long *          )malloc(   MemSize*sizeof(long)     );
      Backup_ParAtt = ( real_par (*)[7] )malloc( 7*MemSize*sizeof(real_par) );  // 7 = pos*3, vel*3, time
   }


// synchronize all active particles
         real_par *ParTime   =   amr->Par->Time;
         long_par *ParType   =   amr->Par->Type;
         real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         real_par *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
#  ifdef STORE_PAR_ACC
   const real_par *ParAcc[3] = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };
#  else
   const real_par *ParAcc[3] = { NULL, NULL, NULL };
#  endif

// convert SyncTime from double to real_par in advance to be consistent with particle time
   const real_par SyncTime_Real = (real_par)SyncTime;

   real_par dt;
#  ifdef COMOVING
   bool     dTime2dt, Initialized=false;
   real_par ParTime_Prev=NULL_REAL;
   real_par dt_Prev     =NULL_REAL;
#  endif

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
//    skip inactive particles
      if ( amr->Par->Mass[p] < 0.0 )   continue;

//    skip massive particles when enabling OPT__FREEZE_PAR
      if ( OPT__FREEZE_PAR  &&  ParType[p] != PTYPE_TRACER )   continue;

      if (  ! Mis_CompareRealValue( SyncTime_Real, ParTime[p], NULL, false )  )
      {
//       backup data before synchronization
         if ( SyncOption == PAR_SYNC_TEMP )
         {
//          allocate enough memory for the backup array
            if ( Backup_NPar >= MemSize )
            {
               MemSize      += MemUnit;
               Backup_ParID  = ( long *          )realloc( Backup_ParID,    MemSize*sizeof(long)     );
               Backup_ParAtt = ( real_par (*)[7] )realloc( Backup_ParAtt, 7*MemSize*sizeof(real_par) );
            }

            Backup_ParID [ Backup_NPar ]    = p;
            Backup_ParAtt[ Backup_NPar ][0] = ParPos[0][p];
            Backup_ParAtt[ Backup_NPar ][1] = ParPos[1][p];
            Backup_ParAtt[ Backup_NPar ][2] = ParPos[2][p];
            Backup_ParAtt[ Backup_NPar ][3] = ParVel[0][p];
            Backup_ParAtt[ Backup_NPar ][4] = ParVel[1][p];
            Backup_ParAtt[ Backup_NPar ][5] = ParVel[2][p];
            Backup_ParAtt[ Backup_NPar ][6] = ParTime  [p];

            Backup_NPar ++;
         }


//       synchronize particles
         dt = SyncTime_Real - ParTime[p];

//       convert time-step for comoving
#        ifdef COMOVING
         if ( Initialized )
            dTime2dt    = ( ParTime[p] != ParTime_Prev );

         else
         {
            dTime2dt    = true;
            Initialized = true;
         }

//       avoid redundant calculations
         if ( dTime2dt )
         {
            dt           = (real_par)Mis_dTime2dt( (double)ParTime[p], (double)dt );
            dt_Prev      = dt;
            ParTime_Prev = ParTime[p];
         }

         else
            dt = dt_Prev;
#        endif // #ifdef COMOVING

         for (int d=0; d<3; d++)
         {
            ParPos[d][p] += ParVel[d][p]*dt;
//          only accelerate massive particles
            if ( ParType[p] != PTYPE_TRACER )
               ParVel[d][p] += ParAcc[d][p]*dt;
         }

         ParTime[p] = SyncTime_Real;
      } // if (  ! Mis_CompareRealValue( SyncTime_Real, ParTime[p], NULL, false )  )
   } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)


   CurrentSyncTime = SyncTime;

   return 0;

} // FUNCTION : Par_Synchronize



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Synchronize_Restore
// Description :  Restore particle position, velocity, and time to the values before synchronization
//
// Note        :  1. Only work if Par_Synchronize has been invoked with SyncOption == PAR_SYNC_TEMP
//                2. Use the input SyncTime and the global static variable CurrentSyncTime to determine
//                   whether particles have been synchronized previously
//
// Parameter   :  SyncTime : Target synchronization time
//-------------------------------------------------------------------------------------------------------
void Par_Synchronize_Restore( const double SyncTime )
{

// check whether particles have been synchronized previously
   if ( SyncTime != CurrentSyncTime )
      Aux_Error( ERROR_INFO, "SyncTime (%20.14e) != CurrentSyncTime (%20.14e) for %s !!\n",
                 SyncTime, CurrentSyncTime, __FUNCTION__ );

   if ( Backup_NPar < 0 )  Aux_Error( ERROR_INFO, "backup arrays have NOT been allocated !!\n" );


// restore particle attributes (position, velocity, and time)
   real_par *ParTime   =   amr->Par->Time;
   real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real_par *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   long ParID;

   for (long p=0; p<Backup_NPar; p++)
   {
      ParID            = Backup_ParID [p];
      ParPos[0][ParID] = Backup_ParAtt[p][0];
      ParPos[1][ParID] = Backup_ParAtt[p][1];
      ParPos[2][ParID] = Backup_ParAtt[p][2];
      ParVel[0][ParID] = Backup_ParAtt[p][3];
      ParVel[1][ParID] = Backup_ParAtt[p][4];
      ParVel[2][ParID] = Backup_ParAtt[p][5];
      ParTime  [ParID] = Backup_ParAtt[p][6];
   }


// free memory and set global static variables to indicate that particles are not synchronized anymore
   free( Backup_ParID  );
   free( Backup_ParAtt );

   CurrentSyncTime = -1.0;
   Backup_NPar     = -1;
   Backup_ParID    = NULL;
   Backup_ParAtt   = NULL;

} // FUNCTION : Par_Synchronize_Restore



#endif // #ifdef PARTICLE
