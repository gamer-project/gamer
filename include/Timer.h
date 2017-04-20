#ifndef __TIMER_H__
#define __TIMER_H__



#include "sys/time.h"

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );


//-------------------------------------------------------------------------------------------------------
// Structure   :  Timer_t
// Description :  Data structure for measuring the elapsed time
//
// Data Member :  Status    : The status of each timer : (false / true) <--> (stop / ticking)
//                Time      : The variable recording the elapsed time (in microseconds)
//                WorkingID : The currently working ID of the array "Time"
//                NTimer    : The number of timers
//
// Method      :  Timer_t   : Constructor
//               ~Timer_t   : Destructor
//                Start     : Start timing
//                Stop      : Stop timing
//                GetValue  : Get the elapsed time recorded in timer (in seconds)
//                Reset     : Reset timer
//-------------------------------------------------------------------------------------------------------
struct Timer_t
{

// data members
// ===================================================================================
   bool  *Status;
   ulong *Time;
   uint  WorkingID;
   uint  NTimer;



   //===================================================================================
   // Constructor :  Timer_t
   // Description :  Constructor of the structure "Timer_t"
   //
   // Note        :  Allocate and initialize all data member
   //                a. The WorkingID is initialized as "0"
   //                b. The Status is initialized as "false"
   //                c. The timer is initialized as "0"
   //
   // Parameter   :  N : The number of timers to be allocated
   //===================================================================================
   Timer_t( const uint N )
   {
      NTimer    = N;
      WorkingID = 0;
      Time      = new ulong [NTimer];
      Status    = new bool  [NTimer];

      for (uint t=0; t<NTimer; t++)
      {
         Time  [t] = 0;
         Status[t] = false;
      }
   }



   //===================================================================================
   // Destructor  :  ~Timer_t
   // Description :  Destructor of the structure "Timer_t"
   //
   // Note        :  Release memory
   //===================================================================================
   ~Timer_t()
   {
      delete [] Time;
      delete [] Status;
   }



   //===================================================================================
   // Method      :  Start
   // Description :  Start timing and set status as "true"
   //===================================================================================
   void Start()
   {
#     ifdef GAMER_DEBUG
      if ( WorkingID >= NTimer )
         Aux_Error( ERROR_INFO, "timer is exhausted !!\n" );

      if ( Status[WorkingID] )
         Aux_Message( stderr, "WARNING : timer has already been started (WorkingID = %u) !!\n", WorkingID );
#     endif

      timeval tv;
      gettimeofday( &tv, NULL );

      Time[WorkingID] = tv.tv_sec*1000000 + tv.tv_usec - Time[WorkingID];

      Status[WorkingID] = true;
   }



   //===================================================================================
   // Method      :  Stop
   // Description :  Stop timing and set status as "false"
   //
   // Note        :  Increase the WorkingID if the input parameter "Next == true"
   //
   // Parameter   :  Next : (true / false) --> (increase / not increase) the WorkingID
   //===================================================================================
   void Stop( const bool Next )
   {
#     ifdef GAMER_DEBUG
      if ( WorkingID >= NTimer )
         Aux_Error( ERROR_INFO, "timer is exhausted !!\n" );

      if ( !Status[WorkingID] )
         Aux_Message( stderr, "WARNING : timer has NOT been started (WorkingID = %u) !!\n", WorkingID );
#     endif

      timeval tv;
      gettimeofday( &tv, NULL );

      Time[WorkingID] = tv.tv_sec*1000000 + tv.tv_usec - Time[WorkingID];

      Status[WorkingID] = false;

      if ( Next )    WorkingID++;
   }



   //===================================================================================
   // Method      :  GetValue
   // Description :  Get the elapsed time (in seconds) recorded in the timer "TargetID"
   //
   // Parameter   :  TargetID : The ID of the target timer
   //===================================================================================
   double GetValue( const uint TargetID )
   {
#     ifdef GAMER_DEBUG
      if ( TargetID >= NTimer )
         Aux_Error( ERROR_INFO, "timer does NOT exist (TargetID = %u, NTimer = %u) !!\n",
                    TargetID, NTimer );

      if ( Status[TargetID] )
         Aux_Message( stderr, "WARNING : timer is still ticking (TargetID = %u) !!\n", TargetID );
#     endif

      return Time[TargetID]*1.e-6f;
   }



   //===================================================================================
   // Method      :  Reset
   // Description :  Reset all timers and set WorkingID as "0"
   //
   // Note        :  All timing information previously recorded will be lost after
   //                invoking this function
   //===================================================================================
   void Reset()
   {
      for (uint t=0; t<NTimer; t++)
      {
#        ifdef GAMER_DEBUG
         if ( Status[t] )
            Aux_Message( stderr, "WARNING : resetting a ticking timer (WorkingID = %u) !!\n", t );
#        endif

         Time[t]   = 0;
         WorkingID = 0;
      }
   }


}; // struct Timer_t



// macro for timing functions
#ifdef TIMING

#  define TIMING_FUNC( call, timer, next )                        \
   {                                                              \
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );   \
      timer->Start();                                             \
      call;                                                       \
      timer->Stop( next );                                        \
   }

#else

#  define TIMING_FUNC( call, timer, next )   call

#endif


// macro for timing solvers
#if ( defined TIMING_SOLVER  &&  defined TIMING )

#  ifdef GPU
#     define GPU_SYNC()  CUAPI_Synchronize()
#  else
#     define GPU_SYNC()
#  endif

#  define TIMING_SYNC( call, timer )            \
   {                                            \
      timer->Start();                           \
      call;                                     \
      GPU_SYNC();                               \
      timer->Stop( false );                     \
   }

#else

#  define TIMING_SYNC( call, timer )   call

#endif



#endif // #ifndef __TIMER_H__
