#include "GAMER.h"



// prototypes of built-in source terms
#if ( MODEL == HYDRO )
void Src_Init_Deleptonization();
void Src_Init_ExactCooling();
#endif

// this function pointer can be set by a test problem initializer for a user-specified source term
void (*Src_Init_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Init
// Description :  Initialize the source terms
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set SrcTerms
//                3. Set "Src_Init_User_Ptr" in a test problem initializer for a user-specified source term
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_Init()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if at least one source term is activated
   if (
#       if ( MODEL == HYDRO )
        SrcTerms.Deleptonization  ||
        SrcTerms.ExactCooling     ||
#       endif
        SrcTerms.User
      )
      SrcTerms.Any = true;
   else
      SrcTerms.Any = false;


// set auxiliary parameters
   for (int d=0; d<3; d++)    SrcTerms.BoxCenter[d] = amr->BoxCenter[d];

   SrcTerms.Unit_L = UNIT_L;
   SrcTerms.Unit_M = UNIT_M;
   SrcTerms.Unit_T = UNIT_T;
   SrcTerms.Unit_V = UNIT_V;
   SrcTerms.Unit_D = UNIT_D;
   SrcTerms.Unit_E = UNIT_E;
   SrcTerms.Unit_P = UNIT_P;
#  ifdef MHD
   SrcTerms.Unit_B = UNIT_B;
#  endif


// initialize all function pointers as NULL
#  if ( MODEL == HYDRO )
   SrcTerms.Dlep_FuncPtr              = NULL;
   SrcTerms.Dlep_CPUPtr               = NULL;
#  ifdef GPU
   SrcTerms.Dlep_GPUPtr               = NULL;
#  endif
   SrcTerms.Dlep_AuxArrayDevPtr_Flt   = NULL;
   SrcTerms.Dlep_AuxArrayDevPtr_Int   = NULL;
   SrcTerms.Dlep_Profile_DataDevPtr   = NULL;
   SrcTerms.Dlep_Profile_RadiusDevPtr = NULL;
#  endif

#  if ( MODEL == HYDRO )
   SrcTerms.EC_FuncPtr            = NULL;
   SrcTerms.EC_CPUPtr             = NULL;
#  ifdef GPUEC
   SrcTerms.EC_GPUPtr             = NULL;
#  endif
   SrcTerms.EC_AuxArrayDevPtr_Flt = NULL;
   SrcTerms.EC_AuxArrayDevPtr_Int = NULL;
   SrcTerms.EC_TEF_lambda_DevPtr  = NULL;
   SrcTerms.EC_TEF_alpha_DevPtr   = NULL;
   SrcTerms.EC_TEFc_DevPtr        = NULL;
#  endif

   SrcTerms.User_FuncPtr              = NULL;
   SrcTerms.User_CPUPtr               = NULL;
#  ifdef GPU
   SrcTerms.User_GPUPtr               = NULL;
#  endif
   SrcTerms.User_AuxArrayDevPtr_Flt   = NULL;
   SrcTerms.User_AuxArrayDevPtr_Int   = NULL;


// initialize all source terms
// (1) deleptonization
#  if ( MODEL == HYDRO )
   if ( SrcTerms.Deleptonization )
   {
      Src_Init_Deleptonization();

//    check if the source-term function is set properly
      if ( SrcTerms.Dlep_FuncPtr == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.Dlep_FuncPtr == NULL !!\n" );
      if ( SrcTerms.Dlep_CPUPtr  == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.Dlep_CPUPtr  == NULL !!\n" );
#     ifdef GPU
      if ( SrcTerms.Dlep_GPUPtr  == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.Dlep_GPUPtr  == NULL !!\n" );
#     endif
   }
#  endif

// (2) exact cooling
#  if ( MODEL == HYDRO )
   if ( SrcTerms.ExactCooling )
   { 
      Src_Init_ExactCooling();

//    check if the source-term function is set properly
      if ( SrcTerms.EC_FuncPtr == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.EC_FuncPtr == NULL !!\n" );
      if ( SrcTerms.EC_CPUPtr  == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.EC_CPUPtr  == NULL !!\n" );
#     ifdef GPU
      if ( SrcTerms.EC_GPUPtr  == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.EC_GPUPtr  == NULL !!\n" );
#     endif         
   } 
#  endif    

// (3) user-specified source term
   if ( SrcTerms.User )
   {
      if ( Src_Init_User_Ptr == NULL )       Aux_Error( ERROR_INFO, "Src_Init_User_Ptr == NULL !!\n" );

      Src_Init_User_Ptr();

//    check if the source-term function is set properly
      if ( SrcTerms.User_FuncPtr == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.User_FuncPtr == NULL !!\n" );
      if ( SrcTerms.User_CPUPtr  == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.User_CPUPtr  == NULL !!\n" );
#     ifdef GPU
      if ( SrcTerms.User_GPUPtr  == NULL )   Aux_Error( ERROR_INFO, "SrcTerms.User_GPUPtr  == NULL !!\n" );
#     endif
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Src_Init
