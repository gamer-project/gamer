#include "GAMER.h"
#include "CUFLU.h"


// status of the fluid solver used by AUTO_REDUCE_DT (declared in Flu_AdvanceDt.cpp)
extern int FluStatus_ThisRank;


static void StoreFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                       const int NPG, const int *PID0_List, const real dt );
static void CorrectFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                         const int NPG, const int *PID0_List, const real dt );
#if ( MODEL == HYDRO  ||  MODEL == MHD )
static bool Unphysical( const real Fluid[], const real Gamma_m1, const int CheckMinEngyOrPres );
static void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                               const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                               real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                               char h_DE_Array_F_Out[][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                               real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                               const real dt );
#endif
static int  Table_01( const int lv, const int PID, const int SibID );

extern void CPU_RiemannSolver_Roe ( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );
extern void CPU_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );
extern void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Close
// Description :  1. Save the fluxes across the coarse-fine boundaries at level "lv"
//                2. Correct the fluxes across the coarse-fine boundaries at level "lv-1"
//                3. Copy the data from the "h_Flu_Array_F_Out" and "h_DE_Array_F_Out" arrays to the "amr->patch" pointers
//                4. Get the minimum time-step information of the fluid solver
//
// Parameter   :  lv                : Target refinement level
//                SaveSg            : Sandglass to store the updated data
//                h_Flux_Array      : Host array storing the updated flux data
//                h_Flu_Array_F_Out : Host array storing the updated fluid data
//                h_DE_Array_F_Out  : Host array storing the dual-energy status
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//                dt                : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void Flu_Close( const int lv, const int SaveSg, real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                char h_DE_Array_F_Out[][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                const int NPG, const int *PID0_List, const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                const double dt )
{

// try to correct the unphysical results in h_Flu_Array_F_Out (e.g., negative density)
// --> must be done BEFORE invoking both StoreFlux() and CorrectFlux() since CorrectUnphysical() might modify the flux array
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   CorrectUnphysical( lv, NPG, PID0_List, h_Flu_Array_F_In, h_Flu_Array_F_Out, h_DE_Array_F_Out, h_Flux_Array, dt );
#  endif


// save the flux in the coarse-fine boundary at level "lv"
   if ( OPT__FIXUP_FLUX  &&  lv != TOP_LEVEL )  StoreFlux( lv, h_Flux_Array, NPG, PID0_List, dt );


// correct the flux in the coarse-fine boundary at level "lv-1"
   if ( OPT__FIXUP_FLUX  &&  lv != 0 )    CorrectFlux( lv, h_Flux_Array, NPG, PID0_List, dt );


// copy the updated data from the arrays "h_Flu_Array_F_Out" and "h_DE_Array_F_Out" to each patch pointer
#  if ( FLU_NOUT != NCOMP_TOTAL )
#     error : ERROR : FLU_NOUT != NCOMP_TOTAL (one must specify how to copy data from h_Flu_Array_F_Out to fluid) !!
#  endif

   int I, J, K, KJI, PID0;

#  pragma omp parallel for private( I, J, K, KJI, PID0 ) schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID     = PID0 + LocalID;
         const int Table_x = TABLE_02( LocalID, 'x', 0, PATCH_SIZE );
         const int Table_y = TABLE_02( LocalID, 'y', 0, PATCH_SIZE );
         const int Table_z = TABLE_02( LocalID, 'z', 0, PATCH_SIZE );

         for (int v=0; v<FLU_NOUT; v++)      {
         for (int k=0; k<PATCH_SIZE; k++)    {  K = Table_z + k;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = Table_y + j;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = Table_x + i;

            KJI = K*4*PATCH_SIZE*PATCH_SIZE + J*2*PATCH_SIZE + I;

            amr->patch[SaveSg][lv][PID]->fluid[v][k][j][i] = h_Flu_Array_F_Out[TID][v][KJI];

         }}}}

#        ifdef DUAL_ENERGY
         for (int k=0; k<PATCH_SIZE; k++)    {  K = Table_z + k;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = Table_y + j;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = Table_x + i;

            KJI = K*4*PATCH_SIZE*PATCH_SIZE + J*2*PATCH_SIZE + I;

//          de_status is always stored in Sg=0
            amr->patch[0][lv][PID]->de_status[k][j][i] = h_DE_Array_F_Out[TID][KJI];

         }}}
#        endif
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Flu_Close



//-------------------------------------------------------------------------------------------------------
// Function    :  StoreFlux
// Description :  Save the fluxes across the coarse-fine boundaries for patches at level "lv"
//                (to be fixed later by the function "CorrectFlux")
//
// Parameter   :  lv           : Target refinement level
//                h_Flux_Array : Host array storing the updated flux data
//                NPG          : Number of patch groups to be evaluated
//                PID0_List    : List recording the patch indicies with LocalID==0 to be udpated
//                dt           : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void StoreFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithFlux )
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithFlux is off ??\n", __FUNCTION__ );


   int Table1, Table2, Mapping, mm, nn, PID0;
   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;

#  pragma omp parallel for private( Table1, Table2, Mapping, mm, nn, PID0, FluxPtr ) schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int PID=PID0; PID<PID0+8; PID++)
      for (int s=0; s<6; s++)
      {
//       for bitwise reproducibility, store the fluxes to be corrected in the "flux_bitrep" array
#        ifdef BITWISE_REPRODUCIBILITY
         FluxPtr = amr->patch[0][lv][PID]->flux_bitrep[s];
#        else
         FluxPtr = amr->patch[0][lv][PID]->flux[s];
#        endif

         if ( FluxPtr != NULL )
         {
            switch ( s )
            {
               case 0:  case 1:
                  Mapping = TABLE_02( PID%8, 'x', 0, 1 ) + s;
                  Table1  = TABLE_02( PID%8, 'z', 0, PATCH_SIZE );
                  Table2  = TABLE_02( PID%8, 'y', 0, PATCH_SIZE );
               break;

               case 2:  case 3:
                  Mapping = TABLE_02( PID%8, 'y', 1, 2 ) + s;
                  Table1  = TABLE_02( PID%8, 'z', 0, PATCH_SIZE );
                  Table2  = TABLE_02( PID%8, 'x', 0, PATCH_SIZE );
               break;

               case 4:  case 5:
                  Mapping = TABLE_02( PID%8, 'z', 2, 3 ) + s;
                  Table1  = TABLE_02( PID%8, 'y', 0, PATCH_SIZE );
                  Table2  = TABLE_02( PID%8, 'x', 0, PATCH_SIZE );
               break;

               default:
                  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "s", s );
            }


            for (int v=0; v<NFLUX_TOTAL; v++)   {
            for (int m=0; m<PATCH_SIZE; m++)    {  mm = m + Table1;
            for (int n=0; n<PATCH_SIZE; n++)    {  nn = n + Table2;

               FluxPtr[v][m][n] = dt*h_Flux_Array[TID][Mapping][v][ mm*2*PATCH_SIZE+nn ];

            }}}

         } // if ( FluxPtr != NULL )
      } // for (int PID=PID0; PID<PID0+8; PID++) for (int s=0; s<6; s++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : StoreFlux



//-------------------------------------------------------------------------------------------------------
// Function    :  CorrectFlux
// Description :  Use the fluxes across the coarse-fine boundaries at level "lv" to correct the fluxes
//                at level "lv-1"
//
// Parameter   :  lv           : Target refinement level
//                h_Flux_Array : Array storing the updated flux data
//                NPG          : Number of patch groups to be evaluated
//                PID0_List    : List recording the patch indicies with LocalID==0 to be udpated
//                dt           : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void CorrectFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                  const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithFlux )
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithFlux is off ??\n", __FUNCTION__ );


#  pragma omp parallel
   {
      const int  Mapping  [6] = { 0, 2, 3, 5, 6, 8 };
      const int  MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };
      const real dt_4         = (real)0.25*dt;

      int ID, FaPID, FaSibPID, PID0, mm, nn;
      real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;


#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

         for (int s=0; s<6; s++)
         {
            if ( Table_01( lv, PID0, s ) == -1 )
            {
               FaPID    = amr->patch[0][lv  ][ PID0]->father;
               FaSibPID = amr->patch[0][lv-1][FaPID]->sibling[s];

//             for AUTO_REDUCE_DT, store the updated fluxes in the temporary array "flux_tmp" since
//             we may need to abandon these updated results if the fluid solver fails
               FluxPtr = ( AUTO_REDUCE_DT ) ? amr->patch[0][lv-1][FaSibPID]->flux_tmp[ MirrorSib[s] ] :
                                              amr->patch[0][lv-1][FaSibPID]->flux    [ MirrorSib[s] ];

#              ifdef GAMER_DEBUG
               if ( FluxPtr == NULL )  Aux_Error( ERROR_INFO, "FluxPtr == NULL (PID0 %d, s %d) !!\n", PID0, s );
#              endif

               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<PS2; m++)  { mm = m/2;
               for (int n=0; n<PS2; n++)  { nn = n/2;

                  ID = m*PS2 + n;
                  FluxPtr[v][mm][nn] -= dt_4*h_Flux_Array[TID][ Mapping[s] ][v][ID];

               }}
            } // if ( Table_01( lv, n, s ) == -1 )
         } // for (int s=0; s<6; s++)
      } // for (int TID=0; TID<NPG; TID++)

   } // OpenMP parallel region

} // FUNCTION : CorrectFlux



#if ( MODEL == HYDRO  ||  MODEL == MHD )
//-------------------------------------------------------------------------------------------------------
// Function    :  Unphysical
// Description :  Check whether the input variables are unphysical
//
// Note        :  1. One can put arbitrary criteria here. Cells violating the conditions will be recalculated
//                   in "CorrectUnphysical"
//                2. Currently it is used for MODEL==HYDRO/MHD to check whether the input density and pressure
//                   (or energy density) is smaller than the given thresholds
//                   --> It also checks if any variable is -inf, +inf, or nan
//                   --> It does NOT check if passive scalars are negative
//                       --> We already apply a floor value in CPU_Shared_FullStepUpdate()
//
// Parameter   :  Fluid              : Input fluid variable array with size FLU_NOUT
//                Gamma_m1           : Gamma - 1
//                CheckMinEngyOrPres : (0/1) ==> check (energy/pressure)
//
// Return      :  true/false <==> input Fluid array is unphysical/physical
//-------------------------------------------------------------------------------------------------------
bool Unphysical( const real Fluid[], const real Gamma_m1, const int CheckMinEngyOrPres )
{

   const int  CheckMinEngy = 0;
   const int  CheckMinPres = 1;
   const bool CorrPres_No  = false;


// if any checks below fail, return true
// =================================================
// note that since MIN_DENS and MIN_PRES are declared as double, they must be converted to **real** before the comparison
// --> otherwise LHS in the comparison will be converted from real to double, which is inconsistent with the assignment
//     (e.g., "Update[DENS] = FMAX( Update[DENS], (real)MIN_DENS" )
   if (  !Aux_IsFinite(Fluid[DENS])  ||  !Aux_IsFinite(Fluid[MOMX])  ||  !Aux_IsFinite(Fluid[MOMY])  ||
         !Aux_IsFinite(Fluid[MOMZ])  ||  !Aux_IsFinite(Fluid[ENGY])  ||
         Fluid[DENS] < (real)MIN_DENS  )
      return true;

   if ( CheckMinEngyOrPres == CheckMinEngy  &&  Fluid[ENGY] < (real)MIN_PRES )
      return true;

   if ( CheckMinEngyOrPres == CheckMinPres  &&
        (
#          if ( DUAL_ENERGY == DE_ENPY )
//         when the dual-energy formalism is adopted, do NOT calculate pressure from "Etot-Ekin" since it would suffer
//         from large round-off errors
//         currently we use TINY_NUMBER as the floor value of entropy and hence here we use 2.0*TINY_NUMBER to validate entropy
//         --> in general, for MIN_PRES > 0.0, we expect that unphysical entropy would lead to unphysical pressure
//         --> however, the check "Fluid[ENPY] < (real)2.0*TINY_NUMBER" is necessary when MIN_PRES == 0.0
           CPU_DensEntropy2Pres( Fluid[DENS], Fluid[ENPY], Gamma_m1, CorrPres_No, NULL_REAL ) < (real)MIN_PRES  ||
           Fluid[ENPY] < (real)2.0*TINY_NUMBER

#          else
           CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                            Gamma_m1, CorrPres_No, NULL_REAL ) < (real)MIN_PRES
#          endif
        )
      )
      return true;


// if all checks above pass, return false
// =================================================
   return false;

} // FUNCTION : Unphysical



//-------------------------------------------------------------------------------------------------------
// Function    :  CorrectUnphysical
// Description :  Check if any cell in the output array of Fluid solver contains unphysical retult
//                --> For example, negative density
//
// Note        :  1. Define unphysical values in the function "Unphysical"
//                2. Currently it is only used for MODEL==HYDRO/MHD to check if density or pressure is smaller than
//                   the minimum allowed values (i.e., MIN_DENS and MIN_PRES)
//                   --> It also checks if any variable is -inf, +inf, and nan
//                   --> But one can define arbitrary criteria in "Unphysical" to trigger the correction
//                3. Procedure:
//                   if ( found_unphysical )
//                   {
//                      if ( OPT__1ST_FLUX_CORR ) try to correct unphysical variables using 1st-order fluxes
//                      else                      do nothing
//
//                      apply minimum density and pressure
//
//                      if ( still_found_unphysical )    print debug messages and abort
//                      else                             store corrected results to the output array
//                   }
//
// Parameter   :  lv                : Target refinement level
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated (for debug only)
//                h_Flu_Array_F_In  : Input fluid array
//                h_Flu_Array_F_Out : Output fluid array
//                h_DE_Array_F_Out  : Output dual-energy status array
//                h_Flux_Array      : Output array storing the updated flux data
//                dt                : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                        const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                        real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                        char h_DE_Array_F_Out[][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                        real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                        const real dt )
{

   const real dh           = (real)amr->dh[lv];
   const real dt_dh        = dt/dh;
   const int  didx[3]      = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
   const real  Gamma_m1    = GAMMA - (real)1.0;
   const real _Gamma_m1    = (real)1.0 / Gamma_m1;
   const int  CheckMinEngy = 0;
   const int  CheckMinPres = 1;

   const int LocalID[2][2][2] = { 0, 1, 2, 4, 3, 6, 5, 7 };

// variables shared by all OpenMP threads
   long NCorrThisTime = 0;
   bool CorrectUnphy  = GAMER_SUCCESS;


// OpenMP parallel region
#  pragma omp parallel
   {

// variables private to each OpenMP thread
   real VarL[3][NCOMP_TOTAL], VarC[NCOMP_TOTAL], VarR[3][NCOMP_TOTAL], FluxL[3][NCOMP_TOTAL], FluxR[3][NCOMP_TOTAL];
   real dF[3][NCOMP_TOTAL], Out[NCOMP_TOTAL], Update[NCOMP_TOTAL];
   real FluxL_1D[NCOMP_TOTAL], FluxR_1D[NCOMP_TOTAL];
   int  ijk_out[3];

// variables for OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D
   const int  Corr1D_NBuf          = 1;
   const int  Corr1D_NCell         = 2*Corr1D_NBuf + 1;
// [3][3] = [dimension][i/j/k indices]
// we apply the directionally splitting correction in the order x->y->z
   const int  Corr1D_idx_min[3][3] = { { Corr1D_NBuf,              0,              0 },
                                       { Corr1D_NBuf,    Corr1D_NBuf,              0 },
                                       { Corr1D_NBuf,    Corr1D_NBuf,    Corr1D_NBuf } };
   const int  Corr1D_idx_max[3][3] = { { Corr1D_NBuf, Corr1D_NCell-1, Corr1D_NCell-1 },
                                       { Corr1D_NBuf,    Corr1D_NBuf, Corr1D_NCell-1 },
                                       { Corr1D_NBuf,    Corr1D_NBuf,    Corr1D_NBuf } };
   const int  Corr1D_didx1[3]      = { NCOMP_TOTAL, Corr1D_NCell*NCOMP_TOTAL, SQR(Corr1D_NCell)*NCOMP_TOTAL };
#  if ( DUAL_ENERGY == DE_ENPY )
   const bool CorrPres_Yes         = true;
   const bool CorrPres_No          = false;
#  endif

   int   Corr1D_didx2[3];
   real *Corr1D_InOut_PtrL=NULL, *Corr1D_InOut_PtrC=NULL, *Corr1D_InOut_PtrR=NULL;
   real (*Corr1D_InOut)[Corr1D_NCell][Corr1D_NCell][NCOMP_TOTAL]
         = ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D ) ? new real [Corr1D_NCell][Corr1D_NCell][Corr1D_NCell][NCOMP_TOTAL]
                                                          : NULL;

#  pragma omp for reduction( +:NCorrThisTime ) schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      for (ijk_out[2]=0; ijk_out[2]<PS2; ijk_out[2]++)
      for (ijk_out[1]=0; ijk_out[1]<PS2; ijk_out[1]++)
      for (ijk_out[0]=0; ijk_out[0]<PS2; ijk_out[0]++)
      {
         const int idx_out = (ijk_out[2]*PS2 + ijk_out[1])*PS2 + ijk_out[0];

//       check if the updated values are unphysical
         for (int v=0; v<NCOMP_TOTAL; v++)   Out[v] = h_Flu_Array_F_Out[TID][v][idx_out];

         if ( Unphysical(Out, Gamma_m1, CheckMinPres) )
         {
            const int idx_in = ( (ijk_out[2]+FLU_GHOST_SIZE)*FLU_NXT + (ijk_out[1]+FLU_GHOST_SIZE) )*FLU_NXT + (ijk_out[0]+FLU_GHOST_SIZE);

//          try to correct unphysical results using 1st-order fluxes (which should be more diffusive)
//          ========================================================================================
//          (A) first apply the directionally **unsplitting** correction
//          ========================================================================================
            if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )
            {
//             collect nearby input coserved variables
               for (int v=0; v<NCOMP_TOTAL; v++)
               {
//                here we have assumed that the fluid solver does NOT modify the input fluid array
//                --> not applicable to RTVD and WAF
                  VarC[v] = h_Flu_Array_F_In[TID][v][idx_in];

                  for (int d=0; d<3; d++)
                  {
                     VarL[d][v] = h_Flu_Array_F_In[TID][v][ idx_in - didx[d] ];
                     VarR[d][v] = h_Flu_Array_F_In[TID][v][ idx_in + didx[d] ];
                  }
               }

//             invoke Riemann solver to calculate the fluxes
//             (note that the recalculated flux does NOT include gravity even for UNSPLIT_GRAVITY --> reduce to 1st-order accuracy)
               switch ( OPT__1ST_FLUX_CORR_SCHEME )
               {
                  case RSOLVER_1ST_ROE:
                     for (int d=0; d<3; d++)
                     {
                        CPU_RiemannSolver_Roe( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        CPU_RiemannSolver_Roe( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

                  case RSOLVER_1ST_HLLC:
                     for (int d=0; d<3; d++)
                     {
                        CPU_RiemannSolver_HLLC( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        CPU_RiemannSolver_HLLC( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

                  case RSOLVER_1ST_HLLE:
                     for (int d=0; d<3; d++)
                     {
                        CPU_RiemannSolver_HLLE( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        CPU_RiemannSolver_HLLE( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

                  default:
                     Aux_Error( ERROR_INFO, "unnsupported Riemann solver (%d) !!\n", OPT__1ST_FLUX_CORR_SCHEME );
               }

//             recalculate the first-order solution for a full time-step
               for (int d=0; d<3; d++)
               for (int v=0; v<NCOMP_TOTAL; v++)   dF[d][v] = FluxR[d][v] - FluxL[d][v];

               for (int v=0; v<NCOMP_TOTAL; v++)
                  Update[v] = h_Flu_Array_F_In[TID][v][idx_in] - dt_dh*( dF[0][v] + dF[1][v] + dF[2][v] );

            } // if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )


//          ========================================================================================
//          (B) apply the directionally **splitting** correction if the unsplitting correction failed
//          ========================================================================================
            if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D )
            {
//             apply the dual-energy formalism to correct the internal energy
//             --> gravity solver may update the internal energy and dual-energy variable again when
//                 UNSPLIT_GRAVITY is adopted
//             --> if the corrected internal energy (pressure) passes the test, we don't need to apply the
//                 directionally **splitting** correction
//             --> we do NOT apply the minimum pressure check in CPU_DualEnergyFix() here
//                 --> otherwise the floor value of pressure might disable the 1st-order-flux correction
#              ifdef DUAL_ENERGY
               CPU_DualEnergyFix( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY], Update[ENPY],
                                  h_DE_Array_F_Out[TID][idx_out], Gamma_m1, _Gamma_m1, CorrPres_No, NULL_REAL, DUAL_ENERGY_SWITCH );
#              endif

               if ( Unphysical(Update, Gamma_m1, CheckMinPres) )
               {
//                collect nearby input coserved variables
                  for (int k=0; k<Corr1D_NCell; k++)  { Corr1D_didx2[2] = (k-Corr1D_NBuf)*didx[2];
                  for (int j=0; j<Corr1D_NCell; j++)  { Corr1D_didx2[1] = (j-Corr1D_NBuf)*didx[1];
                  for (int i=0; i<Corr1D_NCell; i++)  { Corr1D_didx2[0] = (i-Corr1D_NBuf)*didx[0];

//                   here we have assumed that the fluid solver does NOT modify the input fluid array
//                   --> not applicable to RTVD and WAF
                     for (int v=0; v<NCOMP_TOTAL; v++)
                        Corr1D_InOut[k][j][i][v] = h_Flu_Array_F_In[TID][v][ idx_in + Corr1D_didx2[0] + Corr1D_didx2[1] + Corr1D_didx2[2] ];
                  }}}

//                apply the 1st-order-flux correction along different spatial directions **successively**
                  for (int d=0; d<3; d++)
                  {
                     for (int k=Corr1D_idx_min[d][2]; k<=Corr1D_idx_max[d][2]; k++)
                     for (int j=Corr1D_idx_min[d][1]; j<=Corr1D_idx_max[d][1]; j++)
                     for (int i=Corr1D_idx_min[d][0]; i<=Corr1D_idx_max[d][0]; i++)
                     {
//                      get the pointers to the left and right states for the Riemann solvers
                        Corr1D_InOut_PtrC = Corr1D_InOut[k][j][i];
                        Corr1D_InOut_PtrL = Corr1D_InOut_PtrC - Corr1D_didx1[d];
                        Corr1D_InOut_PtrR = Corr1D_InOut_PtrC + Corr1D_didx1[d];

//                      invoke Riemann solver to calculate the fluxes
//                      (note that the recalculated flux does NOT include gravity even for UNSPLIT_GRAVITY --> reduce to 1st-order accuracy)
                        switch ( OPT__1ST_FLUX_CORR_SCHEME )
                        {
                           case RSOLVER_1ST_ROE:
                              CPU_RiemannSolver_Roe ( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              CPU_RiemannSolver_Roe ( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                           break;

                           case RSOLVER_1ST_HLLC:
                              CPU_RiemannSolver_HLLC( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              CPU_RiemannSolver_HLLC( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                           break;

                           case RSOLVER_1ST_HLLE:
                              CPU_RiemannSolver_HLLE( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              CPU_RiemannSolver_HLLE( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                           break;

                           default:
                              Aux_Error( ERROR_INFO, "unnsupported Riemann solver (%d) !!\n", OPT__1ST_FLUX_CORR_SCHEME );
                        }

//                      recalculate the first-order solution for a full time-step
                        for (int v=0; v<NCOMP_TOTAL; v++)
                           Corr1D_InOut[k][j][i][v] -= dt_dh*( FluxR_1D[v] - FluxL_1D[v] );

//                      store the 1st-order fluxes used for updating the central cell
                        if ( i == Corr1D_NBuf  &&  j == Corr1D_NBuf  &&  k == Corr1D_NBuf )
                        {
                           for (int v=0; v<NCOMP_TOTAL; v++)
                           {
                              FluxL[d][v] = FluxL_1D[v];
                              FluxR[d][v] = FluxR_1D[v];
                           }
                        }
                     } // i,j,k
                  } // for (int d=0; d<3; d++)

//                store the corrected results in Update[]
                  for (int v=0; v<NCOMP_TOTAL; v++)
                     Update[v] = Corr1D_InOut[Corr1D_NBuf][Corr1D_NBuf][Corr1D_NBuf][v];

               } // if ( Unphysical(Update, Gamma_m1, CheckMinPres) )
            } // if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D )


//          ========================================================================================
//          (C) if OPT__1ST_FLUX_CORR is off, just copy data to Update[] for checking negative density and pressure in the next step
//          ========================================================================================
            if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_NONE )
            {
               for (int v=0; v<NCOMP_TOTAL; v++)   Update[v] = Out[v];
            } // if ( OPT__1ST_FLUX_CORR ) ... else ...


//          ensure positive density
//          --> apply it only when AUTO_REDUCE_DT is disabled
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this density floor
//          --> note that MIN_DENS is declared as double and must be converted to **real** before the comparison
//              --> to be consistent with the check in Unphysical()
//          --> do NOT check the minimum pressure here since we want to apply the dual-energy correction first
            if ( !AUTO_REDUCE_DT )  Update[DENS] = FMAX( Update[DENS], (real)MIN_DENS );


//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Update[v] = FMAX( Update[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               CPU_NormalizePassive( Update[DENS], Update+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif


//          apply the dual-energy formalism to correct the internal energy again
//          --> gravity solver may update the internal energy and dual-energy variable again when
//              UNSPLIT_GRAVITY is adopted
//          --> this might be redundant when OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D
//              --> but it ensures the consistency between all fluid variables since we apply the floor value
//                  of density AFTER the 1st-order-flux correction
//          --> we apply the minimum pressure check in CPU_DualEnergyFix() here only when AUTO_REDUCE_DT is disabled
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this pressure floor
#           ifdef DUAL_ENERGY
            CPU_DualEnergyFix( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY], Update[ENPY],
                               h_DE_Array_F_Out[TID][idx_out], Gamma_m1, _Gamma_m1,
                               (AUTO_REDUCE_DT)?CorrPres_No:CorrPres_Yes, MIN_PRES, DUAL_ENERGY_SWITCH );

//          ensure positive pressure if dual-energy formalism is not adopted
//          --> apply it only when AUTO_REDUCE_DT is disabled
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this pressure floor
#           else
            if ( !AUTO_REDUCE_DT )
            Update[ENGY] = CPU_CheckMinPresInEngy( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                                   Gamma_m1, _Gamma_m1, MIN_PRES );
#           endif


//          check if the newly updated values are still unphysical
//          --> note that, when AUTO_REDUCE_DT is disabled, we check **energy** instead of pressure since even after
//              calling CPU_CheckMinPresInEngy() we may still have pressure < MIN_PRES due to round-off errors
//              (especially when pressure << kinematic energy)
//              --> it will not crash the code since we always apply MIN_PRES when calculating pressure
//          --> when AUTO_REDUCE_DT is enabled, we still check **pressure** instead of energy
            if ( Unphysical(Update, Gamma_m1, (AUTO_REDUCE_DT)?CheckMinPres:CheckMinEngy) )
            {
//             set CorrectUnphy = GAMER_FAILED if any cells fail
//             --> use critical directive to avoid thread racing (may not be necessary here?)
#              pragma omp critical
               CorrectUnphy = GAMER_FAILED;


//             output the debug information (only if AUTO_REDUCE_DT is disabled)
               if ( !AUTO_REDUCE_DT )
               {
                  const int  PID_Failed      = PID0_List[TID] + LocalID[ijk_out[2]/PS1][ijk_out[1]/PS1][ijk_out[0]/PS1];
                  const bool CheckMinPres_No = false;
                  real In[NCOMP_TOTAL], tmp[NCOMP_TOTAL];

                  char FileName[100];
                  sprintf( FileName, "FailedPatchGroup_r%03d_lv%02d_PID0-%05d", MPI_Rank, lv, PID0_List[TID] );

//                use "a" instead of "w" since there may be more than one failed cell in a given patch group
                  FILE *File = fopen( FileName, "a" );

                  for (int v=0; v<NCOMP_TOTAL; v++)   In[v] = h_Flu_Array_F_In[TID][v][idx_in];

//                output information about the failed cell
                  fprintf( File, "PID                              = %5d\n", PID_Failed );
                  fprintf( File, "(i,j,k) in the patch             = (%2d,%2d,%2d)\n", ijk_out[0]%PS1, ijk_out[1]%PS1, ijk_out[2]%PS1 );
                  fprintf( File, "(i,j,k) in the input fluid array = (%2d,%2d,%2d)\n", ijk_out[0], ijk_out[1], ijk_out[2] );
#                 ifdef DUAL_ENERGY
                  fprintf( File, "total energy density update      = %s\n",
                           ( h_DE_Array_F_Out[TID][idx_out] == DE_UPDATED_BY_ETOT     ) ? "Etot" :
                           ( h_DE_Array_F_Out[TID][idx_out] == DE_UPDATED_BY_DUAL     ) ? "Dual" :
                           ( h_DE_Array_F_Out[TID][idx_out] == DE_UPDATED_BY_MIN_PRES ) ? "MinPres" :
                                                                                          "Unknown" );
#                 endif
                  fprintf( File, "\n" );
                  fprintf( File, "input        = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                           CPU_GetPressure(In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                                           Gamma_m1, CheckMinPres_No, NULL_REAL) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", In[ENPY] );
#                 endif
                  fprintf( File, ")\n" );
                  fprintf( File, "ouptut (old) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                           CPU_GetPressure(Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                                           Gamma_m1, CheckMinPres_No, NULL_REAL) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", Out[ENPY] );
#                 endif
                  fprintf( File, ")\n" );
                  fprintf( File, "output (new) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                           CPU_GetPressure(Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                           Gamma_m1, CheckMinPres_No, NULL_REAL) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", Update[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

//                output all data in the input fluid array (including ghost zones)
                  fprintf( File, "\n===============================================================================================\n" );
                  fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s%14s%14s%14s",
                           "i", "j", "k", "Density", "Px", "Py", "Pz", "Energy", "Pressure" );
#                 if ( NCOMP_PASSIVE > 0 )
                  for (int v=0; v<NCOMP_PASSIVE; v++)    fprintf( File, "%14s", PassiveFieldName_Grid[v] );
#                 endif
                  fprintf( File, "\n" );

                  for (int k=0; k<FLU_NXT; k++)
                  for (int j=0; j<FLU_NXT; j++)
                  for (int i=0; i<FLU_NXT; i++)
                  {
                     fprintf( File, "(%2d,%2d,%2d)", i-FLU_GHOST_SIZE, j-FLU_GHOST_SIZE, k-FLU_GHOST_SIZE );

                     for (int v=0; v<NCOMP_TOTAL; v++)
                     {
                        tmp[v] = h_Flu_Array_F_In[TID][v][ ((k*FLU_NXT)+j)*FLU_NXT+i ];
                        fprintf( File, " %13.6e", tmp[v] );
                     }

                     fprintf( File, " %13.6e\n", CPU_GetPressure(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                                                 Gamma_m1, CheckMinPres_No, NULL_REAL) );
                  }

                  fclose( File );

//                output the failed patch (mainly for recording the sibling information)
#                 ifdef GRAVITY
                  Output_Patch( lv, PID_Failed, amr->FluSg[lv], amr->PotSg[lv], "Unphy" );
#                 else
                  Output_Patch( lv, PID_Failed, amr->FluSg[lv], NULL_INT,       "Unphy" );
#                 endif
               } // if ( !AUTO_REDUCE_DT )
            } // if ( Unphysical(Update) )

            else
            {
//             store the corrected solution
               for (int v=0; v<NCOMP_TOTAL; v++)   h_Flu_Array_F_Out[TID][v][idx_out] = Update[v];


//             replace the original coarse-fine boundary fluxes with the 1st-order fluxes for the flux fix-up
               if ( OPT__FIXUP_FLUX  &&  OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )
               {
                  const int dim_map[3][2] = { {1, 2}, {0, 2}, {0, 1} };

                  int   idx_m, idx_n, idx_store, FaceIdx;
                  bool  Store1stFlux;
                  real *FluxPtr=NULL;

                  for (int d=0; d<3; d++)
                  {
                     Store1stFlux = false;

                     if ( ijk_out[d] == 0  ||  ijk_out[d] == PS1 )
                     {
                        FluxPtr      = FluxL[d];
                        FaceIdx      = 3*d + ijk_out[d]/PS1;
                        Store1stFlux = true;
                     }

                     else if ( ijk_out[d] == PS1-1  ||  ijk_out[d] == PS2-1 )
                     {
                        FluxPtr      = FluxR[d];
                        FaceIdx      = 3*d + (ijk_out[d]+1)/PS1;
                        Store1stFlux = true;
                     }

                     if ( Store1stFlux )
                     {
                        idx_m     = ijk_out[ dim_map[d][1] ];
                        idx_n     = ijk_out[ dim_map[d][0] ];
                        idx_store = idx_m*PS2 + idx_n;

                        for (int v=0; v<NCOMP_TOTAL; v++)   h_Flux_Array[TID][FaceIdx][v][idx_store] = FluxPtr[v];
                     }
                  } // for (int d=0; d<3; d++)
               } // if ( OPT__FIXUP_FLUX  &&  OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )

//             record the number of corrected cells
               NCorrThisTime ++;

            } // if ( Unphysical(Update) ) ... else ...
         } // if need correction
      } // i,j,k
   } // for (int TID=0; TID<NPG; TID++)

// "delete" applies to NULL as well
   delete [] Corr1D_InOut;

   } // end of OpenMP parallel region


// operations when CorrectUnphysical() fails
   if ( CorrectUnphy == GAMER_FAILED )
   {
//    if AUTO_REDUCE_DT is enabled, set FluStatus_ThisRank as GAMER_FAILED to rerun Flu_AdvancedDt() with a smaller dt
      if ( AUTO_REDUCE_DT )
         FluStatus_ThisRank = GAMER_FAILED;

//    otherwise, terminate the program
      else
         Aux_Error( ERROR_INFO, "first-order correction failed at Rank %d, lv %d, Time %20.14e, Step %ld, Counter %ld ...\n",
                    MPI_Rank, lv, Time[lv], Step, AdvanceCounter[lv] );
   }

// accumulate the total number of corrected cells in one global time-step if CorrectUnphysical() works
   else
   {
      NCorrUnphy[lv] += NCorrThisTime;
   }

} // FUNCTION : CorrectUnphysical
#endif // #if ( MODEL == HYDRO  ||  MODEL == MHD )



// ============
// |  Tables  |
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  Return the sibling patch ID for the function "CorrectFlux"
//
// Parameter   :  lv    : Target refinement level
//                PID   : Target patch index
//                SibID : Target sibling index (0~5)
//-------------------------------------------------------------------------------------------------------
int Table_01( const int lv, const int PID, const int SibID )
{

   switch ( SibID )
   {
      case 0: case 2: case 4:    return amr->patch[0][lv][PID  ]->sibling[SibID];
      case 1: case 3: case 5:    return amr->patch[0][lv][PID+7]->sibling[SibID];
      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );
         return EXIT_FAILURE;
   }

} // FUNCTION : Table_01
