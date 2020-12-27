#include "GAMER.h"



// prototypes of built-in source terms
void Src_Init_Deleptonization();

// this function pointer can be set by a test problem initializer for a non-built-in source term
void (*Src_Init_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Init
// Description :  Initialize the source terms
//
// Note        :  1. Invoked by Init_GAMER()
//                2. For a non-built-in source term, "Src_Init_User_Ptr" must be set by a test problem initializer
//                   in advance
//                3. Set SRC_TERMS
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_Init()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if at least one source term is activated
   if ( SRC_TERMS.Deleptonization  ||
        SRC_TERMS.User
      )
      SRC_TERMS.Any = true;
   else
      SRC_TERMS.Any = false;


// initialize all function pointers as NULL
   SRC_TERMS.User_FuncPtr = NULL;


// set auxiliary parameters
   for (int d=0; d<3; d++)    SRC_TERMS.BoxCenter[d] = amr->BoxCenter[d];

   SRC_TERMS.Unit_L = UNIT_L;
   SRC_TERMS.Unit_M = UNIT_M;
   SRC_TERMS.Unit_T = UNIT_T;
   SRC_TERMS.Unit_V = UNIT_V;
   SRC_TERMS.Unit_D = UNIT_D;
   SRC_TERMS.Unit_E = UNIT_E;
   SRC_TERMS.Unit_P = UNIT_P;
#  ifdef MHD
   SRC_TERMS.Unit_B = UNIT_B;
#  endif


// initialize source terms
   if ( SRC_TERMS.Deleptonization )
      Src_Init_Deleptonization();

   if ( SRC_TERMS.User )
   {
      if ( Src_Init_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "Src_Init_User_Ptr == NULL !!\n" );

      Src_Init_User_Ptr();

//    check if the source-term function is set properly
      if ( SRC_TERMS.User_FuncPtr == NULL )  Aux_Error( ERROR_INFO, "SrcTerms.User_FuncPtr == NULL !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Src_Init
