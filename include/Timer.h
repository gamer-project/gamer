#ifndef __TIMER_H__
#define __TIMER_H__



#include "sys/time.h"

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  Timer_t
// Description :  Data structure for measuring the elapsed time
//
// Data Member :  Status : (false / true) <--> (stop / ticking)
//                Time   : Variable recording the elapsed time (in microseconds)
//
// Method      :  Timer_t  : Constructor
//               ~Timer_t  : Destructor
//                Start    : Start timing
//                Stop     : Stop timing
//                GetValue : Get the elapsed time recorded in timer (in seconds)
//                Reset    : Reset timer
//-------------------------------------------------------------------------------------------------------
struct Timer_t
{

// data members
// ===================================================================================
   bool  Status;
   ulong Time;



   //===================================================================================
   // Constructor :  Timer_t
   // Description :  Constructor of the structure "Timer_t"
   //
   // Note        :  Initialize all data members
   //===================================================================================
   Timer_t()
   {
      Time   = 0;
      Status = false;
   }



   //===================================================================================
   // Destructor  :  ~Timer_t
   // Description :  Destructor of the structure "Timer_t"
   //
   // Note        :  Release memory
   //===================================================================================
   ~Timer_t()
   {
   }



   //===================================================================================
   // Method      :  Start
   // Description :  Start timing and set status as "true"
   //
   // Note        :  1. Timer must not already be running
   //                2. The same timer can be started and stopped multiple times and the
   //                   results will be accumulated
   //                   --> Invoke Reset() to reset a timer
   //===================================================================================
   void Start()
   {
#     ifdef GAMER_DEBUG
      if ( Status )  Aux_Message( stderr, "WARNING : timer has already been started !!\n" );
#     endif

      timeval tv;
      gettimeofday( &tv, NULL );

      Time   = tv.tv_sec*1000000 + tv.tv_usec - Time;
      Status = true;
   }



   //===================================================================================
   // Method      :  Stop
   // Description :  Stop timing and set status as "false"
   //
   // Note        :  Timer must already be running
   //===================================================================================
   void Stop()
   {
#     ifdef GAMER_DEBUG
      if ( !Status )    Aux_Message( stderr, "WARNING : timer has NOT been started !!\n" );
#     endif

      timeval tv;
      gettimeofday( &tv, NULL );

      Time   = tv.tv_sec*1000000 + tv.tv_usec - Time;
      Status = false;
   }



   //===================================================================================
   // Method      :  GetValue
   // Description :  Get the elapsed time (in seconds) recorded in the timer
   //
   // Note        :  Timer must not be running
   //===================================================================================
   double GetValue()
   {
#     ifdef GAMER_DEBUG
      if ( Status )  Aux_Message( stderr, "WARNING : timer is still ticking !!\n" );
#     endif

      return Time*1.0e-6;
   }



   //===================================================================================
   // Method      :  Reset
   // Description :  Reset the timer
   //
   // Note        :  1. The timing information previously recorded will be lost after
   //                   invoking this function
   //                2. Timer must not be running
   //===================================================================================
   void Reset()
   {
#     ifdef GAMER_DEBUG
      if ( Status )  Aux_Message( stderr, "WARNING : resetting a ticking timer !!\n" );
#     endif

      Time = 0;
   }


}; // struct Timer_t



// macro for timing functions
#ifdef TIMING

#  define TIMING_FUNC( call, timer )                              \
   {                                                              \
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );   \
      timer->Start();                                             \
      call;                                                       \
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );   \
      timer->Stop();                                              \
   }

#else

#  define TIMING_FUNC( call, timer )   call

#endif


// macro for timing solvers
#if ( defined TIMING_SOLVER  &&  defined TIMING )

#  ifdef GPU
#     define GPU_SYNC()  CUAPI_Synchronize()
#  else
#     define GPU_SYNC()
#  endif

#  define TIMING_SYNC( call, timer )                              \
   {                                                              \
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );   \
      timer->Start();                                             \
      call;                                                       \
      GPU_SYNC();                                                 \
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );   \
      timer->Stop();                                              \
   }

#else

#  define TIMING_SYNC( call, timer )   call

#endif



#endif // #ifndef __TIMER_H__
