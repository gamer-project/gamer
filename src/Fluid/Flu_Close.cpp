#include "GAMER.h"
#include "CUFLU.h"


// status of the fluid solver used by AUTO_REDUCE_DT (declared in Flu_AdvanceDt.cpp)
extern int FluStatus_ThisRank;


static void StoreFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                       const int NPG, const int *PID0_List, const real dt );
static void CorrectFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                         const int NPG, const int *PID0_List, const real dt );
#if ( MODEL == HYDRO )
static bool Unphysical( const real Fluid[], const real Gamma_m1, const int CheckMinEngyOrPres );
static void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                               const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                               real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                               char h_DE_Array_F_Out[][ CUBE(PS2) ],
                               real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                               const real dt );
#ifdef MHD
void StoreElectric( const int lv, const real h_Ele_Array[][9][NCOMP_ELE][ PS2_P1*PS2 ],
                    const int NPG, const int *PID0_List, const real dt );
#endif
#endif // #if ( MODEL == HYDRO )
extern void Hydro_RiemannSolver_Roe ( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real Gamma, const real MinPres );
extern void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real Gamma, const real MinPres );
extern void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real Gamma, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Close
// Description :  1. Save the fluxes across the coarse-fine boundaries at level "lv"
//                2. Correct the fluxes across the coarse-fine boundaries at level "lv-1"
//                3. Copy the data from the "h_Flu_Array_F_Out" and "h_DE_Array_F_Out" arrays to the "amr->patch" pointers
//                4. Get the minimum time-step information of the fluid solver
//
// Parameter   :  lv                : Target refinement level
//                SaveSg_Flu        : Sandglass to store the updated fluid data
//                SaveSg_Mag        : Sandglass to store the updated B field (for MHD only)
//                h_Flux_Array      : Host array storing the updated flux data
//                h_Ele_Array       : Host array storing the updated electric field (for MHD only)
//                h_Flu_Array_F_Out : Host array storing the updated fluid data
//                h_Mag_Array_F_Out : Host array storing the updated B field (for MHD only)
//                h_DE_Array_F_Out  : Host array storing the dual-energy status (for DUAL_ENERGY only)
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//                dt                : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void Flu_Close( const int lv, const int SaveSg_Flu, const int SaveSg_Mag,
                real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                real h_Ele_Array[][9][NCOMP_ELE][ PS2_P1*PS2 ],
                real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                real h_Mag_Array_F_Out[][NCOMP_MAG][ PS2_P1*SQR(PS2) ],
                char h_DE_Array_F_Out[][ CUBE(PS2) ],
                const int NPG, const int *PID0_List,
                const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                const double dt )
{

// try to correct the unphysical results in h_Flu_Array_F_Out (e.g., negative density)
// --> must be done BEFORE invoking both StoreFlux() and CorrectFlux() since CorrectUnphysical() might modify the flux array
#  if ( MODEL == HYDRO )
   CorrectUnphysical( lv, NPG, PID0_List, h_Flu_Array_F_In, h_Flu_Array_F_Out, h_DE_Array_F_Out, h_Flux_Array, dt );
#  endif


// operations related to flux fix-up
   if ( OPT__FIXUP_FLUX )
   {
//    save the fluxes on the coarse-fine boundaries at level "lv"
      StoreFlux( lv, h_Flux_Array, NPG, PID0_List, dt );

//    correct the fluxes on the coarse-fine boundaries at level "lv-1"
      CorrectFlux( lv, h_Flux_Array, NPG, PID0_List, dt );
   }


// operations related to electric field fix-up
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC )
   {
//    save the electric field on the coarse-fine boundaries at level "lv"
      if ( lv != TOP_LEVEL )  StoreElectric( lv, h_Ele_Array, NPG, PID0_List, dt );
   }
#  endif


// copy the updated data from output arrays to the corresponding patch pointers
#  if ( FLU_NOUT != NCOMP_TOTAL )
#     error : ERROR : FLU_NOUT != NCOMP_TOTAL (one must specify how to copy data from h_Flu_Array_F_Out to fluid) !!
#  endif


#  pragma omp parallel for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID     = PID0 + LocalID;
         const int Table_x = TABLE_02( LocalID, 'x', 0, PATCH_SIZE );
         const int Table_y = TABLE_02( LocalID, 'y', 0, PATCH_SIZE );
         const int Table_z = TABLE_02( LocalID, 'z', 0, PATCH_SIZE );

         int I, J, K, KJI;

//       fluid variables
         for (int v=0; v<FLU_NOUT; v++)      {
         for (int k=0; k<PATCH_SIZE; k++)    {  K = Table_z + k;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = Table_y + j;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = Table_x + i;

            KJI = IDX321( I, J, K, PS2, PS2 );

            amr->patch[SaveSg_Flu][lv][PID]->fluid[v][k][j][i] = h_Flu_Array_F_Out[TID][v][KJI];

         }}}}

//       dual-energy status
#        ifdef DUAL_ENERGY
         for (int k=0; k<PATCH_SIZE; k++)    {  K = Table_z + k;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = Table_y + j;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = Table_x + i;

            KJI = IDX321( I, J, K, PS2, PS2 );

//          de_status is always stored in Sg=0
            amr->patch[0][lv][PID]->de_status[k][j][i] = h_DE_Array_F_Out[TID][KJI];

         }}}
#        endif

//       magnetic field
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG; v++)
         {
            int ijk_end[3], idx=0;

            for (int d=0; d<3; d++)    ijk_end[d] = ( d == v ) ? PS1+1 : PS1;

            for (int k=0; k<ijk_end[2]; k++)    {  K = Table_z + k;
            for (int j=0; j<ijk_end[1]; j++)    {  J = Table_y + j;
            for (int i=0; i<ijk_end[0]; i++)    {  I = Table_x + i;

               KJI = IDX321( I, J, K, ijk_end[0]+PS1, ijk_end[1]+PS1 );

               amr->patch[SaveSg_Mag][lv][PID]->magnetic[v][ idx ++ ] = h_Mag_Array_F_Out[TID][v][KJI];

            }}}
         } // for (int v=0; v<NCOMP_MAG; v++)
#        endif // #ifdef MHD

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Flu_Close



//-------------------------------------------------------------------------------------------------------
// Function    :  StoreFlux
// Description :  Save the coarse-grid fluxes across the coarse-fine boundaries for patches at level "lv"
//                --> to be fixed later by CorrectFlux()
//
// Parameter   :  lv           : Target refinement level
//                h_Flux_Array : Host array storing the updated flux data
//                NPG          : Number of patch groups to be evaluated
//                PID0_List    : List recording the patch indicies with LocalID==0 to be udpated
//                dt           : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void StoreFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithFlux )   Aux_Error( ERROR_INFO, "amr->WithFlux is off !!\n" );


// nothing to do on the highest level
   if ( lv == TOP_LEVEL )  return;


#  pragma omp parallel for schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;

         for (int s=0; s<6; s++)
         {
//          for bitwise reproducibility, store the fluxes to be corrected in flux_bitrep[]
#           ifdef BITWISE_REPRODUCIBILITY
            real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][PID]->flux_bitrep[s];
#           else
            real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][PID]->flux[s];
#           endif

            if ( FluxPtr != NULL )
            {
               int face_idx, disp_m, disp_n;

               switch ( s )
               {
                  case 0:  case 1:
                     face_idx = TABLE_02( LocalID, 'x', 0, 1 ) + s;
                     disp_m   = TABLE_02( LocalID, 'z', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'y', 0, PS1 );
                     break;

                  case 2:  case 3:
                     face_idx = TABLE_02( LocalID, 'y', 1, 2 ) + s;
                     disp_m   = TABLE_02( LocalID, 'z', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'x', 0, PS1 );
                     break;

                  case 4:  case 5:
                     face_idx = TABLE_02( LocalID, 'z', 2, 3 ) + s;
                     disp_m   = TABLE_02( LocalID, 'y', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'x', 0, PS1 );
                     break;

                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "s", s );
               }


               for (int v=0; v<NFLUX_TOTAL; v++)   {
               for (int m=0; m<PS1; m++)           {  const int mm = m + disp_m;
               for (int n=0; n<PS1; n++)           {  const int nn = n + disp_n;

                  FluxPtr[v][m][n] = dt*h_Flux_Array[TID][face_idx][v][ mm*PS2 + nn ];

               }}}

            } // if ( FluxPtr != NULL )
         } // for (int s=0; s<6; s++)
      } // for (int LocalID=0; LocalID<8; LocalID++)
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
void CorrectFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                  const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithFlux )   Aux_Error( ERROR_INFO, "amr->WithFlux is off !!\n" );


// nothing to do on the root level
   if ( lv == 0 )    return;


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
            FaPID = amr->patch[0][lv][PID0]->father;

#           ifdef GAMER_DEBUG
            if ( FaPID < 0 )
               Aux_Error( ERROR_INFO, "FaPID = %d < 0 (lv %d, PID0 %d) !!\n", FaPID, lv, PID0 );
#           endif

            FaSibPID = amr->patch[0][lv-1][FaPID]->sibling[s];

#           ifdef GAMER_DEBUG
            if ( FaSibPID == -1 )
               Aux_Error( ERROR_INFO, "FaSibPID == -1 (lv %d, FaPID %d, s %d, PID0 %d) !!\n", lv-1, FaPID, s, PID0 );
#           endif

//          skip patches adjacent to non-periodic boundaries
            if ( FaSibPID < -1 )    continue;

//          for AUTO_REDUCE_DT, store the updated fluxes in the temporary array flux_tmp[] since
//          we may need to abandon them if the fluid solver fails
            FluxPtr = ( AUTO_REDUCE_DT ) ? amr->patch[0][lv-1][FaSibPID]->flux_tmp[ MirrorSib[s] ] :
                                           amr->patch[0][lv-1][FaSibPID]->flux    [ MirrorSib[s] ];

//          skip patches not adjacent to coarse-fine boundaries
            if ( FluxPtr == NULL )  continue;

//          store fluxes
            for (int v=0; v<NFLUX_TOTAL; v++)
            for (int m=0; m<PS2; m++)  { mm = m/2;
            for (int n=0; n<PS2; n++)  { nn = n/2;

               ID = m*PS2 + n;
               FluxPtr[v][mm][nn] -= dt_4*h_Flux_Array[TID][ Mapping[s] ][v][ID];

            }}
         } // for (int s=0; s<6; s++)
      } // for (int TID=0; TID<NPG; TID++)

   } // OpenMP parallel region

} // FUNCTION : CorrectFlux



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  Unphysical
// Description :  Check whether the input variables are unphysical
//
// Note        :  1. One can put arbitrary criteria here. Cells violating the conditions will be recalculated
//                   in CorrectUnphysical()
//                2. Currently it is used for MODEL==HYDRO to check whether the input density and pressure
//                   (or energy density) is smaller than the given thresholds
//                   --> It also checks if any variable is -inf, +inf, or nan
//                   --> It does NOT check if passive scalars are negative
//                       --> We already apply a floor value in Hydro_Shared_FullStepUpdate()
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

#  ifndef DUAL_ENERGY
#  ifdef MHD
#  warning : WAIT MHD !!!
   const real EngyB = NULL_REAL;
#  else
   const real EngyB = NULL_REAL;
#  endif
#  endif // #indef DUAL_ENERGY

   if ( CheckMinEngyOrPres == CheckMinPres  &&
         (
//          when the dual-energy formalism is adopted, do NOT calculate pressure from "Etot-Ekin" since it would suffer
//          from large round-off errors
//          currently we use TINY_NUMBER as the floor value of entropy and hence here we use 2.0*TINY_NUMBER to validate entropy
//          --> in general, for MIN_PRES > 0.0, we expect that unphysical entropy would lead to unphysical pressure
//          --> however, the check "Fluid[ENPY] < (real)2.0*TINY_NUMBER" is necessary when MIN_PRES == 0.0
#           if   ( DUAL_ENERGY == DE_ENPY )
            Hydro_DensEntropy2Pres( Fluid[DENS], Fluid[ENPY], Gamma_m1, CorrPres_No, NULL_REAL ) < (real)MIN_PRES  ||
            Fluid[ENPY] < (real)2.0*TINY_NUMBER

#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!

#           else // without DUAL_ENERGY
            Hydro_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                               Gamma_m1, CorrPres_No, NULL_REAL, EngyB ) < (real)MIN_PRES

#           endif // DUAL_ENERGY
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
//                2. Currently it is only used for MODEL==HYDRO to check if density or pressure is smaller than
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
                        const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                        real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                        char h_DE_Array_F_Out[][ CUBE(PS2) ],
                        real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
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
//                --> not applicable to RTVD
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
                        Hydro_RiemannSolver_Roe( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        Hydro_RiemannSolver_Roe( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

#                 ifndef MHD
                  case RSOLVER_1ST_HLLC:
                     for (int d=0; d<3; d++)
                     {
                        Hydro_RiemannSolver_HLLC( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        Hydro_RiemannSolver_HLLC( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;
#                 endif

                  case RSOLVER_1ST_HLLE:
                     for (int d=0; d<3; d++)
                     {
                        Hydro_RiemannSolver_HLLE( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        Hydro_RiemannSolver_HLLE( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

#                 ifdef MHD
#                 warning : WAIT MHD !!!
                  case RSOLVER_1ST_HLLD:
                     /*
                     for (int d=0; d<3; d++)
                     {
                        CPU_RiemannSolver_HLLD( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        CPU_RiemannSolver_HLLD( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     */
                     break;
#                 endif

                  default:
                     Aux_Error( ERROR_INFO, "unnsupported Riemann solver (%d) !!\n", OPT__1ST_FLUX_CORR_SCHEME );
               } // switch ( OPT__1ST_FLUX_CORR_SCHEME )

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
//             --> we do NOT apply the minimum pressure check in Hydro_DualEnergyFix() here
//                 --> otherwise the floor value of pressure might disable the 1st-order-flux correction
#              ifdef DUAL_ENERGY
               Hydro_DualEnergyFix( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY], Update[ENPY],
                                    h_DE_Array_F_Out[TID][idx_out], Gamma_m1, _Gamma_m1, CorrPres_No, NULL_REAL, DUAL_ENERGY_SWITCH );
#              endif

               if ( Unphysical(Update, Gamma_m1, CheckMinPres) )
               {
//                collect nearby input coserved variables
                  for (int k=0; k<Corr1D_NCell; k++)  { Corr1D_didx2[2] = (k-Corr1D_NBuf)*didx[2];
                  for (int j=0; j<Corr1D_NCell; j++)  { Corr1D_didx2[1] = (j-Corr1D_NBuf)*didx[1];
                  for (int i=0; i<Corr1D_NCell; i++)  { Corr1D_didx2[0] = (i-Corr1D_NBuf)*didx[0];

//                   here we have assumed that the fluid solver does NOT modify the input fluid array
//                   --> not applicable to RTVD
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
                              Hydro_RiemannSolver_Roe ( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              Hydro_RiemannSolver_Roe ( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                           break;

#                          ifndef MHD
                           case RSOLVER_1ST_HLLC:
                              Hydro_RiemannSolver_HLLC( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              Hydro_RiemannSolver_HLLC( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                           break;
#                          endif

                           case RSOLVER_1ST_HLLE:
                              Hydro_RiemannSolver_HLLE( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              Hydro_RiemannSolver_HLLE( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                           break;

#                          ifdef MHD
#                          warning : WAIT MHD !!!
                           case RSOLVER_1ST_HLLD:
                              /*
                              CPU_RiemannSolver_HLLD( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, GAMMA, MIN_PRES );
                              CPU_RiemannSolver_HLLD( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, GAMMA, MIN_PRES );
                              */
                           break;
#                          endif

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
               Hydro_NormalizePassive( Update[DENS], Update+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif


//          apply the dual-energy formalism to correct the internal energy again
//          --> gravity solver may update the internal energy and dual-energy variable again when
//              UNSPLIT_GRAVITY is adopted
//          --> this might be redundant when OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D
//              --> but it ensures the consistency between all fluid variables since we apply the floor value
//                  of density AFTER the 1st-order-flux correction
//          --> we apply the minimum pressure check in Hydro_DualEnergyFix() here only when AUTO_REDUCE_DT is disabled
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this pressure floor
#           ifdef DUAL_ENERGY
            Hydro_DualEnergyFix( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY], Update[ENPY],
                                 h_DE_Array_F_Out[TID][idx_out], Gamma_m1, _Gamma_m1,
                                 (AUTO_REDUCE_DT)?CorrPres_No:CorrPres_Yes, MIN_PRES, DUAL_ENERGY_SWITCH );

//          ensure positive pressure if dual-energy formalism is not adopted
//          --> apply it only when AUTO_REDUCE_DT is disabled
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this pressure floor
#           else
            if ( !AUTO_REDUCE_DT )
            {
#              ifdef MHD
#              warning : WAIT MHD !!!
               const real EngyB = NULL_REAL;
#              else
               const real EngyB = NULL_REAL;
#              endif
               Update[ENGY] = Hydro_CheckMinPresInEngy( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                                        Gamma_m1, _Gamma_m1, MIN_PRES, EngyB );
            }
#           endif


//          check if the newly updated values are still unphysical
//          --> note that, when AUTO_REDUCE_DT is disabled, we check **energy** instead of pressure since even after
//              calling Hydro_CheckMinPresInEngy() we may still have pressure < MIN_PRES due to round-off errors
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

#                 ifdef MHD
#                 warning : WAIT MHD !!!
                  const real EngyB_In     = NULL_REAL;
                  const real EngyB_Out    = NULL_REAL;
                  const real EngyB_Update = NULL_REAL;
#                 else
                  const real EngyB_In     = NULL_REAL;
                  const real EngyB_Out    = NULL_REAL;
                  const real EngyB_Update = NULL_REAL;
#                 endif

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

                  fprintf( File, "               (%14s, %14s, %14s, %14s, %14s, %14s",
                           FieldLabel[DENS], FieldLabel[MOMX], FieldLabel[MOMY], FieldLabel[MOMZ], FieldLabel[ENGY], "Pressure" );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14s", FieldLabel[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

                  fprintf( File, "input        = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                           Hydro_GetPressure(In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                                             Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB_In) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", In[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

                  fprintf( File, "ouptut (old) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                           Hydro_GetPressure(Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                                             Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB_Out) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", Out[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

                  fprintf( File, "output (new) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                           Hydro_GetPressure(Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                             Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB_Update) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", Update[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

//                output all data in the input fluid array (including ghost zones)
                  fprintf( File, "\nFull input array including ghost zones\n" );
                  fprintf( File, "===============================================================================================\n" );
                  fprintf( File, "(%2s,%2s,%2s)", "i", "j", "k" );
                  for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, "%14s", FieldLabel[v] );
                  fprintf( File, "%14s", "Pressure" );
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

#                    ifdef MHD
#                    warning : WAIT MHD !!!
                     const real EngyB_tmp = NULL_REAL;
#                    else
                     const real EngyB_tmp = NULL_REAL;
#                    endif

                     fprintf( File, " %13.6e\n", Hydro_GetPressure(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                                                   Gamma_m1, CheckMinPres_No, NULL_REAL, EngyB_tmp) );
                  }

                  fclose( File );

//                output the failed patch (mainly for recording the sibling information)
#                 ifdef MHD
                  const int MagSg = amr->MagSg[lv];
#                 else
                  const int MagSg = NULL_INT;
#                 endif
#                 ifdef GRAVITY
                  const int PotSg = amr->PotSg[lv];
#                 else
                  const int PotSg = NULL_INT;
#                 endif
                  Output_Patch( lv, PID_Failed, amr->FluSg[lv], MagSg, PotSg, "Unphy" );
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



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  StoreElectric
// Description :  Save the coarse-grid electric field on the coarse-fine boundaries for patches at level "lv"
//                --> to be fixed later by CorrectElectric()
//
// Parameter   :  lv          : Target refinement level
//                h_Ele_Array : Host array storing the updated electric field
//                NPG         : Number of patch groups to be evaluated
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//                dt          : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void StoreElectric( const int lv, const real h_Ele_Array[][9][NCOMP_ELE][ PS2_P1*PS2 ],
                    const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithElectric )  Aux_Error( ERROR_INFO, "amr->WithElectric is off !!\n" );


#  pragma omp parallel for schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;

//       (1) store the electric field on faces 0 ~ 5
         for (int s=0; s<6; s++)
         {
//          for bitwise reproducibility, store the electric field to be corrected in electric_bitrep[]
            const int EleSize = PS1_M1*PS1;
#           ifdef BITWISE_REPRODUCIBILITY
            real (*ElePtr)[EleSize] = ( real (*)[EleSize] )amr->patch[0][lv][PID]->electric_bitrep[s];
#           else
            real (*ElePtr)[EleSize] = ( real (*)[EleSize] )amr->patch[0][lv][PID]->electric[s];
#           endif

            if ( ElePtr != NULL )
            {
               int face_idx, disp_m, disp_n;

               switch ( s )
               {
                  case 0:  case 1:
                     face_idx = TABLE_02( LocalID, 'x', 0, 1 ) + s;
                     disp_m   = TABLE_02( LocalID, 'z', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'y', 0, PS1 );
                     break;

                  case 2:  case 3:
                     face_idx = TABLE_02( LocalID, 'y', 1, 2 ) + s;
                     disp_m   = TABLE_02( LocalID, 'x', 0, PS1 ); // (x,z) instead of (z,x)
                     disp_n   = TABLE_02( LocalID, 'z', 0, PS1 );
                     break;

                  case 4:  case 5:
                     face_idx = TABLE_02( LocalID, 'z', 2, 3 ) + s;
                     disp_m   = TABLE_02( LocalID, 'y', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'x', 0, PS1 );
                     break;

                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "s", s );
               }

               const int disp_m_p1 = disp_m + 1;
               const int disp_n_p1 = disp_n + 1;

               for (int m=0; m<PS1_M1; m++)  {  const int mm = m + disp_m_p1;
               for (int n=0; n<PS1;    n++)  {  const int nn = n + disp_n;

                  ElePtr[0][ m*PS1    + n ] = dt*h_Ele_Array[TID][face_idx][0][ mm*PS2    + nn ];
               }}

               for (int m=0; m<PS1;    m++)  {  const int mm = m + disp_m;
               for (int n=0; n<PS1_M1; n++)  {  const int nn = n + disp_n_p1;

                  ElePtr[1][ m*PS1_M1 + n ] = dt*h_Ele_Array[TID][face_idx][1][ mm*PS2_P1 + nn ];
               }}
            } // if ( ElePtr != NULL )
         } // for (int s=0; s<6; s++)


//       (2) store the electric field on edges 6 ~ 17
         real *ElePtr[12];

         for (int s=6; s<18; s++)
         {
//          for bitwise reproducibility, store the electric field to be corrected in electric_bitrep[]
#           ifdef BITWISE_REPRODUCIBILITY
            ElePtr[s-6] = amr->patch[0][lv][PID]->electric_bitrep[s];
#           else
            ElePtr[s-6] = amr->patch[0][lv][PID]->electric[s];
#           endif
         }

         if ( ElePtr[0] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x',   0, PS1 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 3, 4 );
            for (int n=0; n<PS1; n++)  ElePtr[ 0][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[1] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 3, 4 );
            for (int n=0; n<PS1; n++)  ElePtr[ 1][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[2] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x',   0, PS1 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 4, 5 );
            for (int n=0; n<PS1; n++)  ElePtr[ 2][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[3] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 4, 5 );
            for (int n=0; n<PS1; n++)  ElePtr[ 3][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[4] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y',   0, PS1 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 6, 7 );
            for (int n=0; n<PS1; n++)  ElePtr[ 4][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[5] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 6, 7 );
            for (int n=0; n<PS1; n++)  ElePtr[ 5][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[6] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y',   0, PS1 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 7, 8 );
            for (int n=0; n<PS1; n++)  ElePtr[ 6][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[7] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 7, 8 );
            for (int n=0; n<PS1; n++)  ElePtr[ 7][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[8] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z',   0, PS1 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 0, 1 );
            for (int n=0; n<PS1; n++)  ElePtr[ 8][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[9] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 0, 1 );
            for (int n=0; n<PS1; n++)  ElePtr[ 9][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[10] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z',   0, PS1 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 1, 2 );
            for (int n=0; n<PS1; n++)  ElePtr[10][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[11] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 1, 2 );
            for (int n=0; n<PS1; n++)  ElePtr[11][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : StoreElectric
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )
