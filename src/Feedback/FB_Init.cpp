#include "GAMER.h"

#ifdef FEEDBACK



// prototypes of built-in feedbacks
void FB_Init_SNe();
void FB_Init_Resolved_SNeII();


// user-specified feedback to be set by a test problem initializer
void (*FB_Init_User_Ptr)() = NULL;


// random number generators
RandomNumber_t *FB_RNG = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init
// Description :  Initialize feedback
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set "FB_Init_User_Ptr" in a test problem initializer for a user-specified feedback
//                3. Allocate NT thread-safe random number generators, where NT is the number of OpenMP threads
//                   --> Freed by FB_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if at least one feedback is activated
   FB_Any = ( FB_SNE  ||  FB_RESOLVED_SNEII  ||  FB_USER );

   if ( ! FB_Any  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : all feedback options (e.g., FB_SNE) are inactivated for FEEDBACK !!\n" );


// call the initialization routines of different feedbacks
   if ( FB_SNE )             FB_Init_SNe();

   if ( FB_RESOLVED_SNEII )  FB_Init_Resolved_SNeII();

   if ( FB_USER )
   {
      if ( FB_Init_User_Ptr != NULL )  FB_Init_User_Ptr();
      else                             Aux_Error( ERROR_INFO, "FB_Init_User_Ptr == NULL for FB_USER !!\n" );

//    check if FB_User_Ptr is set properly by FB_Init_User_Ptr()
      if ( FB_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "FB_User_Ptr == NULL for FB_USER !!\n" );
   }


// allocate the random number generators for all OpenMP threads
   int NT;
#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {  NT = omp_get_num_threads();  }
#  else
   {  NT = 1;                      }
#  endif

   FB_RNG = new RandomNumber_t( NT );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : FB_Init



#endif // #ifdef FEEDBACK
