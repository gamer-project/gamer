#include "GAMER.h"

#if ( MODEL == ELBDM  &&  defined GAMER_DEBUG )
void ELBDM_GetPhase_DebugOnly( real *CData, const int CSize );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Refine
// Description :  Construct patches at level "lv+1" according to the flagging result at level "lv"
//
// Note        :  1. This function will also construct buffer patches at level "lv+1" by calling the function
//                   "Refine_Buffer"
//                2. Data of all sibling-buffer patches must be prepared in advance for creating new
//                   fine-grid patches by spatial interpolation
//                3. If LOAD_BALANCE is turned on, this function will invoke "LB_Refine" and then return
//
// Parameter   :  lv          : Target refinement level to be refined
//                UseLBFunc   : Use the load-balance alternative functions for the grid refinement
//                              --> USELB_YES : use the load-balance alternative functions
//                                  USELB_NO  : do not use the load-balance alternative functions
//-------------------------------------------------------------------------------------------------------
void Refine( const int lv, const UseLBFunc_t UseLBFunc )
{

// invoke the load-balance refine function
#  ifdef LOAD_BALANCE
   if ( UseLBFunc == USELB_YES )
   {
      LB_Refine( lv );
      return;
   }
#  endif


// check
   if ( lv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : function <%s> should NOT be applied to the finest level !!\n",
                   __FUNCTION__ );
      return;
   }

// check the synchronization
   if ( NPatchTotal[lv+1] != 0 )    Mis_CompareRealValue( Time[lv], Time[lv+1], __FUNCTION__, true );


   const int  Width       = PATCH_SIZE * amr->scale[lv+1];  // scale of a single patch      at level "lv+1"
   const int CFluSg       = amr->FluSg[lv  ];               // sandglass of fluid variables at level "lv"
   const int FFluSg       = amr->FluSg[lv+1];               // sandglass of fluid variables at level "lv+1"
#  ifdef GRAVITY
   const int  CPotSg      = amr->PotSg[lv  ];               // sandglass of potential       at level "lv"
   const int  FPotSg      = amr->PotSg[lv+1];               // sandglass of potential       at level "lv+1"
   const bool SelfGravity = ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH );
#  endif
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   const real  Gamma_m1   = GAMMA - (real)1.0;
   const real _Gamma_m1   = (real)1.0 / Gamma_m1;
#  endif

   int *Cr            = NULL;    // corner coordinates
   int *BufGrandTable = NULL;    // table recording the patch IDs of grandson buffer patches
   int *BufSonTable   = NULL;    // table recording the linking index of each buffer father patch to BufGrandTable
   patch_t *Pedigree  = NULL;    // pointer storing the relation of the target patch at level "lv"


// parameters for spatial interpolation
   const int CRange[3]     = { PATCH_SIZE, PATCH_SIZE, PATCH_SIZE };
   const int FSize         = 2*PATCH_SIZE;
   const int FStart[3]     = { 0, 0, 0 };
   const int FSize3[3]     = { FSize, FSize, FSize };

   int NSide_Flu, CGhost_Flu;
   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, CGhost_Flu );

   const int CSize_Flu     = PATCH_SIZE + 2*CGhost_Flu;
   const int CStart_Flu[3] = { CGhost_Flu, CGhost_Flu, CGhost_Flu };
   const int CSize_Flu3[3] = { CSize_Flu, CSize_Flu, CSize_Flu };

   real Flu_CData[NCOMP_TOTAL][CSize_Flu][CSize_Flu][CSize_Flu];  // coarse-grid fluid array for interpolation
   real Flu_FData[NCOMP_TOTAL][FSize][FSize][FSize];              // fine-grid fluid array storing the interpolation result

#  ifdef GRAVITY
   int NSide_Pot, CGhost_Pot;
   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, CGhost_Pot );

   const int CSize_Pot     = PATCH_SIZE + 2*CGhost_Pot;
   const int CStart_Pot[3] = { CGhost_Pot, CGhost_Pot, CGhost_Pot };

   real Pot_CData[CSize_Pot][CSize_Pot][CSize_Pot];   // coarse-grid potential array for interpolation
   real Pot_FData[FSize][FSize][FSize];               // fine-grid potential array storing the interpolation result
#  endif


// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   const int  NDer       = 0;
   const int *DerVarList = NULL;
   int BC_Face[26], BC_Face_tmp[3], FluVarIdxList[NCOMP_TOTAL];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1 )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1 )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1 )   BC_Face[s] = BC_Face_tmp[0];
   }

   for (int v=0; v<NCOMP_TOTAL; v++)   FluVarIdxList[v] = v;



// a. record the tables "BufGrandTable" and "BufFathTable"
// *** the grandson patch indices must be stored in advance before we remove all buffer patches at level "lv+1"
// ------------------------------------------------------------------------------------------------
   if ( lv < NLEVEL-2 )
   {
      const int NBufSon = amr->NPatchComma[lv+1][27] - amr->NPatchComma[lv+1][1];
      const int NBufFa  = amr->NPatchComma[lv  ][27] - amr->NPatchComma[lv  ][1];

      BufGrandTable = new int [NBufSon];
      BufSonTable   = new int [NBufFa ];

//    initialize the table BufSonTable as -1
//#     pragma omp parallel for
      for (int t=0; t<NBufFa; t++)  BufSonTable[t] = -1;

//#     pragma omp parallel for
      for (int m=0; m<NBufSon; m+=8)
      {
//       record the grandson patch ID
         for (int n=m; n<m+8; n++)
         {
            const int BufSonPID = amr->NPatchComma[lv+1][1] + n;

            BufGrandTable[n] = amr->patch[0][lv+1][BufSonPID]->son;
         }

//       record the index of BufGrandTable array for father patches with son
         const int BufFaID = amr->patch[0][lv+1][ amr->NPatchComma[lv+1][1] + m ]->father
                             - amr->NPatchComma[lv][1];

         BufSonTable[BufFaID] = m;
      }
   } // if ( lv < NLEVEL-2 )



// b. deallocate all buffer patches at level "lv+1"
// ------------------------------------------------------------------------------------------------
//#  pragma omp parallel for
   for (int PID=amr->NPatchComma[lv+1][1]; PID<amr->NPatchComma[lv+1][27]; PID++)
   {
      amr->patch[0][lv+1][PID]->son = -1;
      amr->pdelete( lv+1, PID, OPT__REUSE_MEMORY );
   }

//#  pragma omp parallel for
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][27]; PID++)
      amr->patch[0][lv][PID]->son = -1;



// c. check the refinement flags for all real patches at level "lv"
// ------------------------------------------------------------------------------------------------
//#  pragma omp parallel for private( ??? )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      Pedigree = amr->patch[0][lv][PID];

//    (c1) construct new child patches if they are newly born (one patch group is allocated at a time)
// ================================================================================================
      if ( Pedigree->flag  &&  Pedigree->son == -1 )
      {

//       (c1.1) construct relation : father -> child
         Pedigree->son = amr->num[lv+1];


//       (c1.2) allocate child patches and construct relation : child -> father
         Cr = Pedigree->corner;

         amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2],       PID, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2],       PID, true, true );
         amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2],       PID, true, true );
         amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2]+Width, PID, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2],       PID, true, true );
         amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2]+Width, PID, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2]+Width, PID, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, PID, true, true );


//       (c1.3) assign data to child patches by spatial interpolation
//       (c1.3.1) fill up the central region of CData
         int I, J, K;

//       fluid data
         for (int v=0; v<NCOMP_TOTAL; v++)   {
         for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost_Flu;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost_Flu;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost_Flu;

            Flu_CData[v][K][J][I] = amr->patch[CFluSg][lv][PID]->fluid[v][k][j][i];

         }}}}

//       potential data
#        ifdef GRAVITY
         if ( SelfGravity )
         for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost_Pot;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost_Pot;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost_Pot;

            Pot_CData[K][J][I] = amr->patch[CPotSg][lv][PID]->pot[k][j][i];

         }}}
#        endif


//       (c1.3.2) fill up the ghost zone of CData (no interpolation is required)
         int    SibPID, Loop[3], Disp1[3], Disp2[3], I2, J2, K2, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];
         double xyz[3];

//       calculate the corner coordinates of the coarse-grid data for the user-specified B.C.
         for (int d=0; d<3; d++)    xyz[d] = Pedigree->EdgeL[d] + (0.5-CGhost_Flu)*amr->dh[lv];

//       (c1.3.2.1) prepare the fluid data
         for (int sib=0; sib<NSide_Flu; sib++)
         {
            SibPID   = Pedigree->sibling[sib];

            Loop [0] = TABLE_01( sib, 'x', CGhost_Flu, PATCH_SIZE, CGhost_Flu );
            Loop [1] = TABLE_01( sib, 'y', CGhost_Flu, PATCH_SIZE, CGhost_Flu );
            Loop [2] = TABLE_01( sib, 'z', CGhost_Flu, PATCH_SIZE, CGhost_Flu );
            Disp1[0] = TABLE_01( sib, 'x', 0, CGhost_Flu, CGhost_Flu+PATCH_SIZE );
            Disp1[1] = TABLE_01( sib, 'y', 0, CGhost_Flu, CGhost_Flu+PATCH_SIZE );
            Disp1[2] = TABLE_01( sib, 'z', 0, CGhost_Flu, CGhost_Flu+PATCH_SIZE );

//          (c1.3.2.1-1) if the target sibling patch exists --> just copy data from the nearby patches at the same level
            if ( SibPID >= 0 )
            {
               Disp2[0] = TABLE_01( sib, 'x', PATCH_SIZE-CGhost_Flu, 0, 0 );
               Disp2[1] = TABLE_01( sib, 'y', PATCH_SIZE-CGhost_Flu, 0, 0 );
               Disp2[2] = TABLE_01( sib, 'z', PATCH_SIZE-CGhost_Flu, 0, 0 );

               for (int v=0; v<NCOMP_TOTAL; v++) {
               for (int k=0; k<Loop[2]; k++) {  K = k + Disp1[2];    K2 = k + Disp2[2];
               for (int j=0; j<Loop[1]; j++) {  J = j + Disp1[1];    J2 = j + Disp2[1];
               for (int i=0; i<Loop[0]; i++) {  I = i + Disp1[0];    I2 = i + Disp2[0];

                  Flu_CData[v][K][J][I] = amr->patch[CFluSg][lv][SibPID]->fluid[v][K2][J2][I2];

               }}}}
            } // if ( SibPID >= 0 )


//          (c1.3.2.1-2) if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
            {
               for (int d=0; d<3; d++)
               {
                  BC_Idx_Start[d] = Disp1[d];
                  BC_Idx_End  [d] = Loop[d] + BC_Idx_Start[d] - 1;
               }

               BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

               switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
               {
#                 if ( MODEL == HYDRO  ||  MODEL == MHD )
                  case BC_FLU_OUTFLOW:
                     Hydro_BoundaryCondition_Outflow   ( Flu_CData[0][0][0], BC_Face[BC_Sibling], NCOMP_TOTAL, CGhost_Flu,
                                                         CSize_Flu, CSize_Flu, CSize_Flu, BC_Idx_Start, BC_Idx_End );
                  break;

                  case BC_FLU_REFLECTING:
                     Hydro_BoundaryCondition_Reflecting( Flu_CData[0][0][0], BC_Face[BC_Sibling], NCOMP_TOTAL, CGhost_Flu,
                                                         CSize_Flu, CSize_Flu, CSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                         FluVarIdxList, NDer, DerVarList );
                  break;
#                 if ( MODEL == MHD )
#                 warning : WAIT MHD !!!
#                 endif
#                 endif

                  case BC_FLU_USER:
                     Flu_BoundaryCondition_User        ( Flu_CData[0][0][0],                      NCOMP_TOTAL,
                                                         CSize_Flu, CSize_Flu, CSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                         FluVarIdxList, Time[lv], amr->dh[lv], xyz, _TOTAL, lv );
                  break;

                  default:
                     Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", OPT__BC_FLU[ BC_Face[BC_Sibling] ] );

               } // switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
            } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )


//          (c1.3.2.1-3) it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
            else if ( SibPID == -1 )
               Aux_Error( ERROR_INFO, "no sibling patch is found for FaLv %d, FaPID %d, sib %d !!\n", lv, PID, sib );

            else
               Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, Sib %d) !!\n", SibPID, PID, sib );

         } // for (int sib=0; sib<NSide_Flu; sib++)


//       (c1.3.2.2) prepare the potential data
#        ifdef GRAVITY
         if ( SelfGravity )
         for (int sib=0; sib<NSide_Pot; sib++)
         {
            SibPID = Pedigree->sibling[sib];

//          it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
#           ifdef GAMER_DEBUG
            if ( SibPID < 0 )
               Aux_Error( ERROR_INFO, "no sibling patch is found for FaLv %d, FaPID %d, sib %d !!\n", lv, PID, sib );
#           endif

            Loop [0] = TABLE_01( sib, 'x', CGhost_Pot, PATCH_SIZE, CGhost_Pot );
            Loop [1] = TABLE_01( sib, 'y', CGhost_Pot, PATCH_SIZE, CGhost_Pot );
            Loop [2] = TABLE_01( sib, 'z', CGhost_Pot, PATCH_SIZE, CGhost_Pot );
            Disp1[0] = TABLE_01( sib, 'x', 0, CGhost_Pot, CGhost_Pot+PATCH_SIZE );
            Disp1[1] = TABLE_01( sib, 'y', 0, CGhost_Pot, CGhost_Pot+PATCH_SIZE );
            Disp1[2] = TABLE_01( sib, 'z', 0, CGhost_Pot, CGhost_Pot+PATCH_SIZE );
            Disp2[0] = TABLE_01( sib, 'x', PATCH_SIZE-CGhost_Pot, 0, 0 );
            Disp2[1] = TABLE_01( sib, 'y', PATCH_SIZE-CGhost_Pot, 0, 0 );
            Disp2[2] = TABLE_01( sib, 'z', PATCH_SIZE-CGhost_Pot, 0, 0 );

            for (int k=0; k<Loop[2]; k++) {  K = k + Disp1[2];    K2 = k + Disp2[2];
            for (int j=0; j<Loop[1]; j++) {  J = j + Disp1[1];    J2 = j + Disp2[1];
            for (int i=0; i<Loop[0]; i++) {  I = i + Disp1[0];    I2 = i + Disp2[0];

               Pot_CData[K][J][I] = amr->patch[CPotSg][lv][SibPID]->pot[K2][J2][I2];

            }}}

         } // for (int sib=0; sib<NSide_Pot; sib++)
#        endif // #ifdef GRAVITY


//       (c1.3.3) perform spatial interpolation
         const bool PhaseUnwrapping_Yes    = true;
         const bool PhaseUnwrapping_No     = false;
         const bool EnsureMonotonicity_Yes = true;
         const bool EnsureMonotonicity_No  = false;

//       (c1.3.3.1) determine which variables require **monotonic** interpolation
         bool Monotonicity[NCOMP_TOTAL];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
#           if ( MODEL == HYDRO )
//          we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
            /*
            if ( v == DENS  ||  v == ENGY  ||  v >= NCOMP_FLUID )
                                             Monotonicity[v] = EnsureMonotonicity_Yes;
            else                             Monotonicity[v] = EnsureMonotonicity_No;
            */
                                             Monotonicity[v] = EnsureMonotonicity_Yes;

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            if ( v != REAL  &&  v != IMAG )  Monotonicity[v] = EnsureMonotonicity_Yes;
            else                             Monotonicity[v] = EnsureMonotonicity_No;

#           else
#           error : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION IN THIS NEW MODEL ??
#           endif // MODEL
         }

//       (c1.3.3.2) interpolation
#        if ( MODEL == ELBDM )
         if ( OPT__INT_PHASE )
         {
//          get the wrapped phase (store in the REAL component)
#           ifdef GAMER_DEBUG
            ELBDM_GetPhase_DebugOnly( &Flu_CData[0][0][0][0], CSize_Flu );
#           else
            for (int k=0; k<CSize_Flu; k++)
            for (int j=0; j<CSize_Flu; j++)
            for (int i=0; i<CSize_Flu; i++)
               Flu_CData[REAL][k][j][i] = ATAN2( Flu_CData[IMAG][k][j][i], Flu_CData[REAL][k][j][i] );
#           endif

//          interpolate density
            Interpolate( &Flu_CData[DENS][0][0][0], CSize_Flu3, CStart_Flu, CRange, &Flu_FData[DENS][0][0][0],
                         FSize3, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No,
                         &EnsureMonotonicity_Yes );

//          interpolate phase
            Interpolate( &Flu_CData[REAL][0][0][0], CSize_Flu3, CStart_Flu, CRange, &Flu_FData[REAL][0][0][0],
                         FSize3, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_Yes,
                         &EnsureMonotonicity_No );
         }

         else // if ( OPT__INT_PHASE )
         {
            for (int v=0; v<NCOMP_TOTAL; v++)
            Interpolate( &Flu_CData[v][0][0][0], CSize_Flu3, CStart_Flu, CRange, &Flu_FData[v][0][0][0],
                         FSize3, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No,
                         Monotonicity );
         }

         if ( OPT__INT_PHASE )
         {
//          retrieve real and imaginary parts
            real Amp, Phase, Rho;

            for (int k=0; k<FSize; k++)
            for (int j=0; j<FSize; j++)
            for (int i=0; i<FSize; i++)
            {
               Phase = Flu_FData[REAL][k][j][i];
               Rho   = Flu_FData[DENS][k][j][i];

//             be careful about the negative density introduced from the round-off errors
               if ( Rho < (real)0.0 )
               {
                  Flu_FData[DENS][k][j][i] = (real)0.0;
                  Rho                      = (real)0.0;
               }

               Amp                      = SQRT( Rho );
               Flu_FData[REAL][k][j][i] = Amp*COS( Phase );
               Flu_FData[IMAG][k][j][i] = Amp*SIN( Phase );
            }
         }

#        else // #if ( MODEL == ELBDM )

         for (int v=0; v<NCOMP_TOTAL; v++)
         Interpolate( &Flu_CData[v][0][0][0], CSize_Flu3, CStart_Flu, CRange, &Flu_FData[v][0][0][0],
                      FSize3, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No,
                      Monotonicity );

#        endif // #if ( MODEL == ELBDM ) ... else


#        ifdef GRAVITY
         const int CSize_Pot_Temp[3] = { CSize_Pot, CSize_Pot, CSize_Pot };

         if ( SelfGravity )
         Interpolate( &Pot_CData[0][0][0], CSize_Pot_Temp, CStart_Pot, CRange, &Pot_FData[0][0][0],
                      FSize3, FStart, 1, OPT__REF_POT_INT_SCHEME, PhaseUnwrapping_No,
                      &EnsureMonotonicity_No );
#        endif

//       (c1.3.3.3) check minimum density and pressure
//       --> note that it's unnecessary to check negative passive scalars thanks to the monotonic interpolation
//       --> but we do renormalize passive scalars here
#        if ( MODEL == HYDRO  ||  MODEL == MHD  ||  MODEL == ELBDM  ||  (defined DENS && NCOMP_PASSIVE>0) )
         for (int k=0; k<FSize; k++)
         for (int j=0; j<FSize; j++)
         for (int i=0; i<FSize; i++)
         {
//          check minimum density
            const real DensOld = Flu_FData[DENS][k][j][i];

            if ( DensOld < MIN_DENS )
            {
//             rescale wave function (unnecessary if OPT__INT_PHASE if off, in which case we will rescale all wave functions later)
#              if ( MODEL == ELBDM )
               if ( OPT__INT_PHASE )
               {
                  const real Rescale = SQRT( (real)MIN_DENS / DensOld );

                  Flu_FData[REAL][k][j][i] *= Rescale;
                  Flu_FData[IMAG][k][j][i] *= Rescale;
               }
#              endif

//             apply minimum density
               Flu_FData[DENS][k][j][i] = MIN_DENS;
            }

#           if ( MODEL == HYDRO  ||  MODEL == MHD )
#           ifdef DUAL_ENERGY
//          ensure consistency between pressure, total energy density, and the dual-energy variable
//          --> here we ALWAYS use the dual-energy variable to correct the total energy density
//          --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//              the runtime parameter DUAL_ENERGY_SWITCH here
            const bool CheckMinPres_Yes = true;
            const real UseEnpy2FixEngy  = HUGE_NUMBER;
            char dummy;    // we do not record the dual-energy status here

            CPU_DualEnergyFix( Flu_FData[DENS][k][j][i], Flu_FData[MOMX][k][j][i], Flu_FData[MOMY][k][j][i],
                               Flu_FData[MOMZ][k][j][i], Flu_FData[ENGY][k][j][i], Flu_FData[ENPY][k][j][i],
                               dummy, Gamma_m1, _Gamma_m1, CheckMinPres_Yes, MIN_PRES, UseEnpy2FixEngy );

#           else
//          check minimum pressure
            Flu_FData[ENGY][k][j][i]
               = CPU_CheckMinPresInEngy( Flu_FData[DENS][k][j][i], Flu_FData[MOMX][k][j][i], Flu_FData[MOMY][k][j][i],
                                         Flu_FData[MOMZ][k][j][i], Flu_FData[ENGY][k][j][i],
                                         Gamma_m1, _Gamma_m1, MIN_PRES );
#           endif // #ifdef DUAL_ENERGY ... else ...
#           endif // #if ( MODEL == HYDRO  ||  MODEL == MHD )

//          normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            if ( OPT__NORMALIZE_PASSIVE )
            {
               real Passive[NCOMP_PASSIVE];

               for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Flu_FData[ NCOMP_FLUID + v ][k][j][i];

               CPU_NormalizePassive( Flu_FData[DENS][k][j][i], Passive, PassiveNorm_NVar, PassiveNorm_VarIdx );

               for (int v=0; v<NCOMP_PASSIVE; v++)    Flu_FData[ NCOMP_FLUID + v ][k][j][i] = Passive[v];
            }
#           endif

         } // i,j,k
#        endif // #if ( MODEL == HYDRO  ||  MODEL == MHD  ||  MODEL == ELBDM )


//       (c1.3.4) copy data from IntData to patch pointers
         int SonPID;

         for (int LocalID=0; LocalID<8; LocalID++)
         {
            SonPID = amr->num[lv+1] - 8 + LocalID;
            Disp1[0] = TABLE_02( LocalID, 'x', 0, PATCH_SIZE );
            Disp1[1] = TABLE_02( LocalID, 'y', 0, PATCH_SIZE );
            Disp1[2] = TABLE_02( LocalID, 'z', 0, PATCH_SIZE );

//          fluid data
            for (int v=0; v<NCOMP_TOTAL; v++)   {
            for (int k=0; k<PATCH_SIZE; k++)    {  K = k + Disp1[2];
            for (int j=0; j<PATCH_SIZE; j++)    {  J = j + Disp1[1];
            for (int i=0; i<PATCH_SIZE; i++)    {  I = i + Disp1[0];

               amr->patch[FFluSg][lv+1][SonPID]->fluid[v][k][j][i] = Flu_FData[v][K][J][I];

            }}}}

//          potential data
#           ifdef GRAVITY
            if ( SelfGravity )
            for (int k=0; k<PATCH_SIZE; k++)    {  K = k + Disp1[2];
            for (int j=0; j<PATCH_SIZE; j++)    {  J = j + Disp1[1];
            for (int i=0; i<PATCH_SIZE; i++)    {  I = i + Disp1[0];

               amr->patch[FPotSg][lv+1][SonPID]->pot[k][j][i] = Pot_FData[K][J][I];

            }}}
#           endif

//          rescale real and imaginary parts to get the correct density in ELBDM if OPT__INT_PHASE is off
#           if ( MODEL == ELBDM )
            real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

            if ( !OPT__INT_PHASE )
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
               Real      = amr->patch[FFluSg][lv+1][SonPID]->fluid[REAL][k][j][i];
               Imag      = amr->patch[FFluSg][lv+1][SonPID]->fluid[IMAG][k][j][i];
               Rho_Wrong = Real*Real + Imag*Imag;
               Rho_Corr  = amr->patch[FFluSg][lv+1][SonPID]->fluid[DENS][k][j][i];

//             be careful about the negative density introduced from the round-off errors
               if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
               {
                  amr->patch[FFluSg][lv+1][SonPID]->fluid[DENS][k][j][i] = (real)0.0;
                  Rescale = (real)0.0;
               }
               else
                  Rescale = SQRT( Rho_Corr/Rho_Wrong );

               amr->patch[FFluSg][lv+1][SonPID]->fluid[REAL][k][j][i] *= Rescale;
               amr->patch[FFluSg][lv+1][SonPID]->fluid[IMAG][k][j][i] *= Rescale;
            }
#           endif
         } // for (int LocalID=0; LocalID<8; LocalID++)


//       (c1.4) pass particles from father to son
#        ifdef PARTICLE
         Par_PassParticle2Son( lv, PID );
#        endif

      } // if ( Pedigree->flag  &&  Pedigree->son == -1 )


//    (c2) remove unflagged child patches if they originally existed (one patch group is removed at a time)
// ================================================================================================
      else if ( !Pedigree->flag  &&  Pedigree->son != -1 )
      {

//       (c2.0) pass particles from sons to father
#        ifdef PARTICLE
         Par_PassParticle2Father( lv, PID );
#        endif


//       (c2.1) deallocate the unflagged child patches
         const int SonPID0 = Pedigree->son;
         const int NewPID0 = SonPID0;
         const int OldPID0 = amr->num[lv+1] - 8;

         for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)     amr->pdelete( lv+1, SonPID, OPT__REUSE_MEMORY );


//       (c2.2) construct relation : father -> son
         Pedigree->son = -1;


//       (c2.3) relink the child patch pointers so that no patch indices are skipped
         if ( NewPID0 != OldPID0 )
         {
            int NewPID, OldPID, GrandPID0, FaPID;

            for (int t=0; t<8; t++)
            {
               NewPID = NewPID0 + t;
               OldPID = OldPID0 + t;

//             swap pointers between the old and new PID
//             --> works no matter OPT__REUSE_MEMORY is on or off
//             --> note that when OPT__REUSE_MEMORY is on, we don't want to set pointers of OldPID as NULL
               for (int Sg=0; Sg<2; Sg++)
                  Aux_SwapPointer( (void**)&amr->patch[Sg][lv+1][OldPID], (void**)&amr->patch[Sg][lv+1][NewPID] );

//             re-construct relation : grandson -> son
               GrandPID0 = amr->patch[0][lv+1][NewPID]->son;
               if ( GrandPID0 != -1 )
               {
                  for (int GrandPID=GrandPID0; GrandPID<GrandPID0+8; GrandPID++)
                     amr->patch[0][lv+2][GrandPID]->father = NewPID;
               }
            }

//          re-construct relation : father -> son
            FaPID = amr->patch[0][lv+1][NewPID0]->father;
            amr->patch[0][lv][FaPID]->son = NewPID0;

         } // if ( NewPID0 != OldPID0 )

      } // else if ( !Pedigree->flag  &&  Pedigree->son != -1 )
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// initialize the amr->NPatchComma list for the buffer patches
   for (int m=1; m<28; m++)   amr->NPatchComma[lv+1][m] = amr->num[lv+1];


// d. refine buffer patches
// ------------------------------------------------------------------------------------------------
   Refine_Buffer( lv, BufSonTable, BufGrandTable );


// deallocate tables
   if ( lv < NLEVEL-2 )
   {
      delete [] BufGrandTable;
      delete [] BufSonTable;
   }



// e. re-construct tables and sibling relations
// ------------------------------------------------------------------------------------------------

// set up the BounP_IDMap for the level just created
   Buf_RecordBoundaryPatch( lv+1 );


// construct relation : siblings
   SiblingSearch( lv+1 );


// allocate flux arrays for level "lv"
   if ( amr->WithFlux )
   Flu_AllocateFluxArray( lv );


// allocate flux arrays for level "lv+1"
   if ( lv < NLEVEL-2  &&  amr->WithFlux )
      Flu_AllocateFluxArray( lv+1 );


// get the IDs of patches for sending and receiving data between neighbor ranks
   Buf_RecordExchangeDataPatchID( lv+1 );


// get the total number of patches at lv+1
   Mis_GetTotalPatchNumber( lv+1 );

} // FUNCTION : Refine



#if ( MODEL == ELBDM  &&  defined GAMER_DEBUG )
//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetPhase_DebugOnly
// Description :  Alternative function to calculate phase in the debug mode so that the functions "Refine" and
//                "LB_Refine_AllocateNewPatch" will give EXACTLY THE SAME RESULTS (even the round-off errors
//                are the same)
//
// Note        :  It is found that it is necessary to use this alternative function to calculate phase
//                in order to make runs with "SERIAL", "NO LOAD_BALANCE", and "LOAD_BALANCE" have exactly
//                the same results
//
// Parameter   :  CData : Coarse-grid array
//                CSize : Size of CData in each direction
//-------------------------------------------------------------------------------------------------------
void ELBDM_GetPhase_DebugOnly( real *CData, const int CSize )
{

   const int CSize_1v = CSize*CSize*CSize;

   real *const CData_Dens = CData + DENS*CSize_1v;
   real *const CData_Real = CData + REAL*CSize_1v;
   real *const CData_Imag = CData + IMAG*CSize_1v;

   for (int t=0; t<CSize_1v; t++)   CData_Real[t] = ATAN2( CData_Imag[t], CData_Real[t] );

} // FUNCTION :
#endif // #if ( MODEL == ELBDM  &&  defined GAMER_DEBUG )
