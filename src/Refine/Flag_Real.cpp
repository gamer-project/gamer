#include "Copyright.h"
#include "GAMER.h"

void Flag_Grandson( const int lv, const int PID, const int LocalID );
void Prepare_for_Lohner( const OptLohnerForm_t Form, const real *Var1D, real *Ave1D, real *Slope1D, const int NVar );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Real
// Description :  Flag the real patches at level "lv" according to the given refinement criteria
//
// Note        :  1. Flag operation of the buffer patches is performed by the function "Flag_Buffer"
//                2. In this function, the buffer patches may still be flagged due to the FLAG_BUFFER_SIZE
//                   extension and the grandson check
//                3. To add new refinement criteria, please edit the function "Flag_Check"
//                4. Definition of the function "Prepare_for_Lohner" is put in the file "Flag_Lohner"
//
// Parameter   :  lv          : Targeted refinement level to be flagged
//                UseLBFunc   : Use the load-balance alternative functions for the grandson check and exchanging
//                              the buffer flags (useless if LOAD_BALANCE is off)
//                              --> USELB_YES : use the load-balance alternative functions
//                                  USELB_NO  : do not use the load-balance alternative functions
//-------------------------------------------------------------------------------------------------------
void Flag_Real( const int lv, const UseLBFunc_t UseLBFunc )
{

// check
   if ( lv == NLEVEL-1 )
      Aux_Error( ERROR_INFO, "function <%s> should NOT be applied to the finest level\" !!\n", __FUNCTION__ );


// initialize all flags as false
#  pragma omp parallel for
   for (int PID=0; PID<amr->num[lv]; PID++)  amr->patch[0][lv][PID]->flag = false;


// set sibling indices
   const int SibID_Array[3][3][3]     = {  { {18, 10, 19}, {14,   4, 16}, {20, 11, 21} }, 
                                           { { 6,  2,  7}, { 0, 999,  1}, { 8,  3,  9} }, 
                                           { {22, 12, 23}, {15,   5, 17}, {24, 13, 25} }  };
   const bool IntPhase_No             = false;                 // for invoking "Prepare_PatchData"
   const bool GetTotDens_No           = false;                 // for invoking "Prepare_PatchData"
   const int  NPG                     = 1;                     // for invoking "Prepare_PatchData"
   const int  Lohner_NGhost           = 2;                     // number of ghost cells for the Lohner error estimator
   const int  Lohner_NCell            = PS1 + 2*Lohner_NGhost; // size of the variable array for Lohner
   const int  Lohner_NAve             = Lohner_NCell - 2;      // size of the average array for Lohner
   const int  Lohner_NSlope           = Lohner_NAve;           // size of the slope array for Lohner
   const IntScheme_t Lohner_IntScheme = INT_MINMOD1D;          // interpolation scheme for Lohner
#  ifndef GRAVITY
   const OptPotBC_t OPT__BC_POT       = BC_POT_NONE;
#  endif
#  ifdef PARTICLE
   const bool PredictPos_No           = false;                 // used by Par_MassAssignment
   const bool InitZero_Yes            = true;
   const bool Periodic_No             = false;
   const bool UnitDens_Yes            = true;
   const bool CheckFarAway_No         = false;
#  endif

//###NOTE: no refinement is allowed near the simulation boundary if the isolated BC for self-gravity is selected
#  ifdef GRAVITY
   const bool NoRefineNearBoundary    =  (  ( OPT__GRAVITY_TYPE == GRAVITY_SELF || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  && 
                                            OPT__BC_POT == BC_POT_ISOLATED  ) ? true : false;
#  else
   const bool NoRefineNearBoundary    = false;
#  endif
   const int  NoRefineBoundaryRegion  = ( NoRefineNearBoundary ) ? PS1*( 1<<(NLEVEL-lv) )*( (1<<lv)-1 ) : NULL_INT;


// set the variables for the Lohner's error estimator
   int Lohner_NVar=0, Lohner_TVar=0, Lohner_Stride;

#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   if ( OPT__FLAG_LOHNER_DENS )  {  Lohner_NVar++;    Lohner_TVar |= _DENS;   }
   if ( OPT__FLAG_LOHNER_ENGY )  {  Lohner_NVar++;    Lohner_TVar |= _ENGY;   }
   if ( OPT__FLAG_LOHNER_PRES )  {  Lohner_NVar++;    Lohner_TVar |= _PRES;   }

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( OPT__FLAG_LOHNER_DENS )
   {
      Lohner_NVar = 2;
      Lohner_TVar = _REAL | _IMAG;
   }

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

   Lohner_Stride = Lohner_NVar*Lohner_NCell*Lohner_NCell*Lohner_NCell;  // stride of array for one patch


// collect particles from all descendant patches
#  ifdef PARTICLE
   if ( OPT__FLAG_NPAR_CELL )    Par_CollectParticleFromDescendant( lv );
#  endif
      

//###ISSUE: use atomic ??      
#  pragma omp parallel
   {
      const real (*Fluid)[PS1][PS1][PS1] = NULL;
      real (*Pres)[PS1][PS1]             = NULL;
      real (*Pot )[PS1][PS1]             = NULL;
      real (*Lohner_Var)                 = NULL;   // array storing the variables for Lohner
      real (*Lohner_Ave)                 = NULL;   // array storing the averages of Lohner_Var for Lohner
      real (*Lohner_Slope)               = NULL;   // array storing the slopes of Lohner_Var for Lohner
      real (*ParCount)[PS1][PS1]         = NULL;   // declare as **real** to be consistent with Par_MassAssignment

      int  i_start, i_end, j_start, j_end, k_start, k_end, SibID, SibPID, PID;
      bool ProperNesting, NextPatch;

#     if   ( MODEL == HYDRO )
      if ( OPT__FLAG_PRES_GRADIENT )   Pres = new real [PS1][PS1][PS1];
#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     endif // MODEL

#     ifdef PARTICLE
      if ( OPT__FLAG_NPAR_CELL )    ParCount = new real [PS1][PS1][PS1];
#     endif

      if ( Lohner_NVar > 0 )    
      { 
         Lohner_Var   = new real [ 8*Lohner_NVar*Lohner_NCell *Lohner_NCell *Lohner_NCell  ]; // 8: number of local patches
         Lohner_Ave   = new real [ 3*Lohner_NVar*Lohner_NAve  *Lohner_NAve  *Lohner_NAve   ]; // 3: X/Y/Z of 1 patch
         Lohner_Slope = new real [ 3*Lohner_NVar*Lohner_NSlope*Lohner_NSlope*Lohner_NSlope ]; // 3: X/Y/Z of 1 patch
      }


//    loop over all REAL patches (the buffer patches will be flagged only due to the FLAG_BUFFER_SIZE
//    extension or the grandson check )
#     pragma omp for schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
//       prepare the ghost-zone data for Lohner
         if ( Lohner_NVar > 0 )    
            Prepare_PatchData( lv, Time[lv], Lohner_Var, Lohner_NGhost, NPG, &PID0, Lohner_TVar, 
                               Lohner_IntScheme, UNIT_PATCH, NSIDE_26, IntPhase_No, OPT__BC_FLU, OPT__BC_POT, GetTotDens_No );


//       loop over all local patches within the same patch group
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            PID = PID0 + LocalID;

//          check the proper-nesting condition
            ProperNesting = true;

            for (int sib=0; sib<26; sib++)
            {
//             for the non-periodic BC we have sibling index = SIB_OFFSET_NONPERIODIC - sibling = -1000 - sibling < 0
               if (  (  NoRefineNearBoundary && amr->patch[0][lv][PID]->sibling[sib] < 0   )  ||
                     ( !NoRefineNearBoundary && amr->patch[0][lv][PID]->sibling[sib] == -1 )     )
               {
                  ProperNesting = false;
                  break;
               }
            }

//          check further --> necessary for OPT__UM_START_DOWNGRADE to avoid refinement near the boundary
            if ( ProperNesting  &&  NoRefineNearBoundary )
            {
               for (int d=0; d<3; d++)
               {
                  int CornerL = amr->patch[0][lv][PID]->corner[d];
                  int CornerR = CornerL + Mis_Cell2Scale( PS1, lv );

                  if ( CornerL <= 0                + NoRefineBoundaryRegion  ||
                       CornerR >= amr->BoxScale[d] - NoRefineBoundaryRegion     )
                  {
                     ProperNesting = false;
                     break;
                  }
               }
            }


//          do flag check only if 26 siblings all exist (proper-nesting constraint)
            if ( ProperNesting )
            {
               NextPatch = false;
               Fluid     = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
#              ifdef GRAVITY
               Pot       = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;
#              endif


//             evaluate pressure
#              if   ( MODEL == HYDRO )
               if ( OPT__FLAG_PRES_GRADIENT )
               {
                  const real Gamma_m1 = GAMMA - (real)1.0;
                  real (*FluData)[PS1][PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
                  real Ek; 

                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     Ek = (real)0.5*( FluData[MOMX][k][j][i]*FluData[MOMX][k][j][i] + 
                                      FluData[MOMY][k][j][i]*FluData[MOMY][k][j][i] +
                                      FluData[MOMZ][k][j][i]*FluData[MOMZ][k][j][i] ) / FluData[DENS][k][j][i];

                     Pres[k][j][i] = Gamma_m1 * ( FluData[ENGY][k][j][i] - Ek );
                  }
               }

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!
#              endif // MODEL


//             evaluate the averages and slopes along x/y/z for Lohner
               if ( Lohner_NVar > 0 )    
                  Prepare_for_Lohner( OPT__FLAG_LOHNER_FORM, Lohner_Var+LocalID*Lohner_Stride, Lohner_Ave, Lohner_Slope,
                                      Lohner_NVar );


//             count the number of particles on each cell
#              ifdef PARTICLE
               if ( OPT__FLAG_NPAR_CELL )
               {
                  long *ParList = NULL;
                  int   NPar;

//                determine the number of particles and the particle list
                  if ( amr->patch[0][lv][PID]->son == -1 )
                  {
                     NPar    = amr->patch[0][lv][PID]->NPar;
                     ParList = amr->patch[0][lv][PID]->ParList;
                  }

                  else
                  {
                     NPar    = amr->patch[0][lv][PID]->NPar_Desc;
                     ParList = amr->patch[0][lv][PID]->ParList_Desc;
                  }

#                 ifdef DEBUG_PARTICLE
                  if ( NPar < 0 )
                     Aux_Error( ERROR_INFO, "NPar (%d) has not been calculated (lv %d, PID %d) !!\n",
                                NPar, lv, PID );

                  if ( NPar > 0  &&  ParList == NULL )
                     Aux_Error( ERROR_INFO, "ParList == NULL for NPar_Desc (%d) > 0 (lv %d, PID %d) !!\n",
                                NPar, lv, PID );
#                 endif

//                deposit particles mass on grids (assuming unit density)
                  Par_MassAssignment( ParList, NPar, PAR_INTERP_NGP, ParCount[0][0], PS1,
                                      amr->patch[0][lv][PID]->EdgeL, amr->dh[lv], PredictPos_No, NULL_REAL,
                                      InitZero_Yes, Periodic_No, NULL, UnitDens_Yes, CheckFarAway_No );
               } // if ( OPT__FLAG_NPAR_CELL )
#              endif // #ifdef PARTICLE


//             loop over all cells within the target patch
               for (int k=0; k<PS1; k++)  {  if ( NextPatch )  break;
                                             k_start = ( k - FLAG_BUFFER_SIZE < 0    ) ? 0 : 1;
                                             k_end   = ( k + FLAG_BUFFER_SIZE >= PS1 ) ? 2 : 1;

               for (int j=0; j<PS1; j++)  {  if ( NextPatch )  break; 
                                             j_start = ( j - FLAG_BUFFER_SIZE < 0    ) ? 0 : 1;
                                             j_end   = ( j + FLAG_BUFFER_SIZE >= PS1 ) ? 2 : 1;

               for (int i=0; i<PS1; i++)  {  if ( NextPatch )  break;
                                             i_start = ( i - FLAG_BUFFER_SIZE < 0    ) ? 0 : 1;
                                             i_end   = ( i + FLAG_BUFFER_SIZE >= PS1 ) ? 2 : 1;

//                check if the targeted cell satisfies the refinement criteria (useless pointers are always == NULL)
                  if (  lv < MAX_LEVEL  &&  Flag_Check( lv, PID, i, j, k, Fluid, Pot, Pres, Lohner_Var+LocalID*Lohner_Stride,
                                                        Lohner_Ave, Lohner_Slope, Lohner_NVar, ParCount )  )
                  {
//                   flag itself
                     amr->patch[0][lv][PID]->flag = true;

//                   flag sibling patches according to the size of FLAG_BUFFER_SIZE
                     for (int kk=k_start; kk<=k_end; kk++)
                     for (int jj=j_start; jj<=j_end; jj++)
                     for (int ii=i_start; ii<=i_end; ii++)
                     {
                        SibID = SibID_Array[kk][jj][ii];

                        if ( SibID != 999 )  
                        {
                           SibPID = amr->patch[0][lv][PID]->sibling[SibID];

#                          ifdef GAMER_DEBUG
                           if ( SibPID == -1 )  
                              Aux_Error( ERROR_INFO, "SibPID == -1 --> proper-nesting check failed !!\n" );

                           if ( SibPID <= SIB_OFFSET_NONPERIODIC  &&  NoRefineNearBoundary )
                              Aux_Error( ERROR_INFO, "SibPID == %d when NoRefineNearBoundary is on !!\n", SibPID );
#                          endif

                           if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
                        }
                     }

//                   for FLAG_BUFFER_SIZE == PATCH_SIZE, once a cell is flagged, all 26 siblings will be flagged
                     if ( FLAG_BUFFER_SIZE == PS1 )   NextPatch = true;

                  } // check flag
               }}} // k, j, i


//             flag based on the number particles per patch (which doesn't need to go through all cells one-by-one)
#              ifdef PARTICLE
               if ( lv < MAX_LEVEL  &&  OPT__FLAG_NPAR_PATCH != 0 )
               {
                  const int NParFlag = FlagTable_NParPatch[lv];
                  bool Flag;

                  if ( amr->patch[0][lv][PID]->son == -1 )  Flag = amr->patch[0][lv][PID]->NPar             > NParFlag;
                  else                                      Flag = Par_CountParticleInDescendant( lv, PID ) > NParFlag;

                  if ( Flag )
                  {
//                   flag itself
                     amr->patch[0][lv][PID]->flag = true;

//                   flag all siblings for OPT__FLAG_NPAR_PATCH == 2
                     if ( OPT__FLAG_NPAR_PATCH == 2 )
                     {
                        for (int s=0; s<26; s++)
                        {
                           SibPID = amr->patch[0][lv][PID]->sibling[s];

#                          ifdef GAMER_DEBUG
                           if ( SibPID == -1 )  
                              Aux_Error( ERROR_INFO, "SibPID == -1 --> proper-nesting check failed !!\n" );

                           if ( SibPID <= SIB_OFFSET_NONPERIODIC  &&  NoRefineNearBoundary )
                              Aux_Error( ERROR_INFO, "SibPID == %d when NoRefineNearBoundary is on !!\n", SibPID );
#                          endif

                           if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
                        }
                     }
                  } // if ( Flag )
               } // if ( OPT__FLAG_NPAR_PATCH != 0 )
#              endif // #ifdef PARTICLE

            } // if ( ProperNesting )
         } // for (int LocalID=0; LocalID<8; LocalID++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


      if ( Pres     != NULL )    delete [] Pres;
      if ( ParCount != NULL )    delete [] ParCount;

      if ( Lohner_NVar > 0 )    
      {
         delete [] Lohner_Var;
         delete [] Lohner_Ave;
         delete [] Lohner_Slope;
      }

   } // OpenMP parallel region


// free variables of descendant particles
#  ifdef PARTICLE
   if ( OPT__FLAG_NPAR_CELL )    Par_CollectParticleFromDescendant_FreeMemory( lv );
#  endif


// apply the proper-nesting constraint again (should also apply to the buffer patches)
#  pragma omp parallel for
   for (int PID=0; PID<amr->num[lv]; PID++)
   {
      for (int sib=0; sib<26; sib++)
      {
         if (  (  NoRefineNearBoundary && amr->patch[0][lv][PID]->sibling[sib] < 0   )  ||
               ( !NoRefineNearBoundary && amr->patch[0][lv][PID]->sibling[sib] == -1 )     )
         {
            amr->patch[0][lv][PID]->flag = false;
            break;
         }
      }

//    check further --> necessary for OPT__UM_START_DOWNGRADE to avoid refinement near the boundary
      if ( amr->patch[0][lv][PID]->flag  &&  NoRefineNearBoundary )
      {
        for (int d=0; d<3; d++)
        {
           int CornerL = amr->patch[0][lv][PID]->corner[d];
           int CornerR = CornerL + Mis_Cell2Scale( PS1, lv );

           if ( CornerL <= 0                + NoRefineBoundaryRegion  ||
                CornerR >= amr->BoxScale[d] - NoRefineBoundaryRegion     )
            {
               amr->patch[0][lv][PID]->flag = false;
               break;
            }
         }
      }
   } // for (int PID=0; PID<amr->num[lv]; PID++)


// invoke the load-balance functions 
#  ifdef LOAD_BALANCE
   if ( UseLBFunc == USELB_YES )
   {
//    grandson check
      if ( lv < NLEVEL-2 )    LB_GrandsonCheck( lv );

//    exchange the flagged buffer patches 
      LB_ExchangeFlaggedBuffer( lv );    
    
      return;
   }
#  endif


   int SonPID;

// grandson check
   if ( lv < NLEVEL-2 )
   {
//###ISSUE: use atomic ??      
//#     pragma omp parallel for private( SonPID )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         SonPID = amr->patch[0][lv][PID]->son;

         if ( SonPID != -1 )
         {
            for (int LocalID=0; LocalID<8; LocalID++)
            {
               if ( amr->patch[0][lv+1][SonPID+LocalID]->son != -1 )    // if grandson exists
               {
//                flag the corresponding siblings of patch PID
                  Flag_Grandson( lv, PID, LocalID );
   
//                flag the patch PID
                  amr->patch[0][lv][PID]->flag = true;
               }
            }
         }
      }
   }  // if ( lv < NLEVEL-2 )


// set up the BounFlag_NList and BounFlag_PosList 
   Buf_RecordBoundaryFlag( lv );

} // FUNCTION : Flag_Real



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Grandson
// Description :  Properly flag the siblings of patch "PID" at level "lv" to ensure that the proper-nesting
//                condition is satisfied at level "lv+2"
//
// Note        :  Properly take care of the non-periodic BC where there are no buffer patches outside the simulation box
//                and hence the proper-nesting constraint cannot be applied to the real patches adajacent to the simulation 
//                boundary (sibling index <= SIB_OFFSET_NONPERIODIC)
//
// Parameter   :  lv      : Targeted level to be flagged
//                PID     : Targeted patch ID at level "lv"
//                LocalID : Index of son (0~7) which has its own son (grandson of patch "PID" at level "lv")
//-------------------------------------------------------------------------------------------------------
void Flag_Grandson( const int lv, const int PID, const int LocalID )
{

   int SibPID;

   switch ( LocalID )
   {
      case 0:
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[18];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[18];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;
      

      case 1:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[19];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[19];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[10];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;
      

      case 2:
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[20];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[20];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[14];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;

      
      case 3:
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[22];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[22];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 6];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;

      
      case 4:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[21];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[21];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[16];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[11];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 4];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;
      

      case 5:
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[24];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 0];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[24];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 8];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[15];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;
      

      case 6:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[23];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[23];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 7];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 2];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[12];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;
       

      case 7:
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[25];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID >= 0 )   amr->patch[0][lv][SibPID]->flag = true;

#        ifdef GAMER_DEBUG
         SibPID = amr->patch[0][lv][PID]->sibling[ 1];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );  
         SibPID = amr->patch[0][lv][PID]->sibling[ 3];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 5];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[25];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[ 9];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[17];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
         SibPID = amr->patch[0][lv][PID]->sibling[13];   if ( SibPID == -1 )  Aux_Error( ERROR_INFO, "SibPID == -1 !!\n" );
#        endif
         break;

      
      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LocalID", LocalID );

   } // switch ( LocalID )

} // FUNCTION : Flag_Grandson
