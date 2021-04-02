#include "GAMER.h"

#ifdef FEEDBACK



// prototypes of built-in feedbacks
void FB_Init_SNE();


// user-specified feedback to be set by a test problem initializer
void (*FB_Init_User_Ptr)() = NULL;

extern void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                            real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init
// Description :  Initialize feedback
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set "FB_Init_User_Ptr" in a test problem initializer for a user-specified feedback
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   if ( FB_USER )
   {
      if ( FB_Init_User_Ptr != NULL )  FB_Init_User_Ptr();
      else                             Aux_Error( ERROR_INFO, "FB_Init_User_Ptr == NULL for FB_USER !!\n" );

//    check if FB_User_Ptr is set properly by FB_Init_User_Ptr()
      if ( FB_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "FB_User_Ptr == NULL for FB_USER !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : FB_Init



#endif // #ifdef FEEDBACK
