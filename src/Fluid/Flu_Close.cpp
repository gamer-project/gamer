#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

static void StoreFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                       const int NPG, const int *PID0_List );
static void CorrectFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                         const int NPG, const int *PID0_List );
#if ( MODEL == HYDRO  ||  MODEL == MHD )
static bool Unphysical( const real Fluid[], const real Gamma_m1, const int CheckMinEngyOrPres );
static void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                               const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                               real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE], const real dt );
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
//                3. Copy the data from the "h_Flu_Array_F_Out" array to the "amr->patch" pointers
//                4. Get the minimum time-step information when the option "OPT__ADAPTIVE_DT" is turned on
//
// Parameter   :  lv                : Targeted refinement level
//                SaveSg            : Sandglass to store the updated data
//                h_Flux_Array      : Host array storing the updated flux data
//                h_MinDtInfo_Array : Host array storing the minimum time-step information in each patch group
//                                    --> useful only if "GetMinDtInfo == true"
//                                    --> NOT supported yet
//                h_Flu_Array_F_Out : Host array storing the updated fluid data
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//                GetMinDtInfo      : true --> Gather the minimum time-step information (the CFL condition in
//                                             HYDRO) among all input patch group
//                                         --> NOT supported yet
//-------------------------------------------------------------------------------------------------------
void Flu_Close( const int lv, const int SaveSg, const real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE],
                const real h_MinDtInfo_Array[], const int NPG, const int *PID0_List, const bool GetMinDtInfo,
                const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT], const double dt )
{

// save the flux in the coarse-fine boundary at level "lv"
   if ( OPT__FIXUP_FLUX  &&  lv != TOP_LEVEL )  StoreFlux( lv, h_Flux_Array, NPG, PID0_List );


// correct the flux in the coarse-fine boundary at level "lv-1"
   if ( OPT__FIXUP_FLUX  &&  lv != 0 )    CorrectFlux( lv, h_Flux_Array, NPG, PID0_List );


// try to correct the unphysical results in h_Flu_Array_F_Out (e.g., negative density)
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   CorrectUnphysical( lv, NPG, PID0_List, h_Flu_Array_F_In, h_Flu_Array_F_Out, dt );
#  endif


// copy the updated data from the array "h_Flu_Array_F_Out" to each patch pointer
   int I, J, K, KJI, PID0;

#  pragma omp parallel for private( I, J, K, KJI, PID0 ) schedule( runtime )
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
      }
   } // for (int TID=0; TID<NPG; TID++)


// store the minimum time-step estimated by the hydrodynamical CFL condition
   /*
   if ( GetMinDtInfo )
   {
      if ( PID_Start == 0 )   MinDtInfo_Fluid[lv] = 0.0;

      for (int t=0; t<NPG; t++)
         MinDtInfo_Fluid[lv] = ( h_MinDtInfo_Array[t] > MinDtInfo_Fluid[lv] ) ? h_MinDtInfo_Fluid_Array[t] :
                                                                                MinDtInfo_Fluid[lv];
   }
   */

} // FUNCTION : Flu_Close



//-------------------------------------------------------------------------------------------------------
// Function    :  StoreFlux
// Description :  Save the fluxes across the coarse-fine boundaries for patches at level "lv"
//                (to be fixed later by the function "CorrectFlux")
//
// Parameter   :  lv             : Targeted refinement level
//                h_Flux_Array   : Host array storing the updated flux data
//                NPG            : Number of patch groups to be evaluated
//                PID0_List      : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void StoreFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                const int NPG, const int *PID0_List )
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
//       for debug mode, store the fluxes to be corrected in the "flux_debug" array
#        ifdef GAMER_DEBUG
         FluxPtr = amr->patch[0][lv][PID]->flux_debug[s];
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

               FluxPtr[v][m][n] = h_Flux_Array[TID][Mapping][v][ mm*2*PATCH_SIZE+nn ];

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
// Parameter   :  lv             : Targeted refinement level
//                h_Flux_Array   : Array storing the updated flux data
//                NPG            : Number of patch groups to be evaluated
//                PID0_List      : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void CorrectFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE],
                  const int NPG, const int *PID0_List )
{

// check
   if ( !amr->WithFlux )
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithFlux is off ??\n", __FUNCTION__ );


#  pragma omp parallel
   {
      const int Mapping  [6] = { 0, 2, 3, 5, 6, 8 };
      const int MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };

      int ID, FaPID, FaSibPID, PID0;
      real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE]   = NULL;
      real (*Flux_Temp)[PATCH_SIZE][PATCH_SIZE] = new real [NFLUX_TOTAL][PATCH_SIZE][PATCH_SIZE];


#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

         for (int s=0; s<6; s++)
         {
            if ( Table_01( lv, PID0, s ) == -1 )
            {
               FaPID    = amr->patch[0][lv  ][    PID0]->father;
               FaSibPID = amr->patch[0][lv-1][   FaPID]->sibling[s];
               FluxPtr  = amr->patch[0][lv-1][FaSibPID]->flux[ MirrorSib[s] ];

#              ifdef GAMER_DEBUG
               if ( FluxPtr == NULL )  Aux_Error( ERROR_INFO, "FluxPtr == NULL (PID0 %d, s %d) !!\n", PID0, s );
#              endif

               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<PATCH_SIZE; m++)
               for (int n=0; n<PATCH_SIZE; n++)    Flux_Temp[v][m][n] = (real)0.0;

               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<2*PATCH_SIZE; m++)
               for (int n=0; n<2*PATCH_SIZE; n++)
               {
                  ID = m*2*PATCH_SIZE + n;

                  Flux_Temp[v][m/2][n/2] -= h_Flux_Array[TID][ Mapping[s] ][v][ID];
               }

               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<PATCH_SIZE; m++)
               for (int n=0; n<PATCH_SIZE; n++)
               {
#                 ifdef INDIVIDUAL_TIMESTEP
                  FluxPtr[v][m][n] += (real)0.125*Flux_Temp[v][m][n];
#                 else
                  FluxPtr[v][m][n] += (real)0.250*Flux_Temp[v][m][n];
#                 endif
               }

            } // if ( Table_01( lv, n, s ) == -1 )
         } // for (int s=0; s<6; s++)
      } // for (int TID=0; TID<NPG; TID++)

      delete [] Flux_Temp;

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
//
// Parameter   :  Fluid              : Input fluid variable array with size FLU_NOUT
//                Gamma_m1           : Gamma - 1
//                CheckMinEngyOrPres : (0/1) ==> check (energy/pressure)
//
// Return      :  true/false <==> input Fluid array is unphysical/physical
//-------------------------------------------------------------------------------------------------------
bool Unphysical( const real Fluid[], const real Gamma_m1, const int CheckMinEngyOrPres )
{

   const bool CorrPres_No  = false;
   const int  CheckMinEngy = 0;
   const int  CheckMinPres = 1;

// note that MIN_DENS and MIN_PRES are declared as double
// --> must convert them to **real** before the comparison
// --> otherwise LHS in the comparison will be converted from real to double, which is inconsistent with the assignment
//     (e.g., "Update[DENS] = FMAX( Update[DENS], (real)MIN_DENS" )
   if (  !isfinite(Fluid[DENS])  ||  !isfinite(Fluid[MOMX])  ||  !isfinite(Fluid[MOMY])  ||
         !isfinite(Fluid[MOMZ])  ||  !isfinite(Fluid[ENGY])  ||
         Fluid[DENS] < (real)MIN_DENS  ||
         ( CheckMinEngyOrPres == CheckMinEngy && Fluid[ENGY] < (real)MIN_PRES )  ||
         ( CheckMinEngyOrPres == CheckMinPres && CPU_GetPressure(Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                                                 Gamma_m1, CorrPres_No, NULL_REAL) < (real)MIN_PRES )  )
      return true;
   else
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
// Parameter   :  lv                : Targeted refinement level
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated (for debug only)
//                h_Flu_Array_F_In  : Input array
//                h_Flu_Array_F_Out : Output array
//                dt                : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                        const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                        real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE], const real dt )
{

   const real dh           = (real)amr->dh[lv];
   const real dt_dh        = dt/dh;
   const int  didx[3]      = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
   const real Gamma_m1     = GAMMA - (real)1.0;
   const real _Gamma_m1    = (real)1.0 / Gamma_m1;
   const int  CheckMinEngy = 0;
   const int  CheckMinPres = 1;

   const int LocalID[2][2][2] = { 0, 1, 2, 4, 3, 6, 5, 7 };

   real VarL[3][NCOMP_TOTAL], VarC[NCOMP_TOTAL], VarR[3][NCOMP_TOTAL], FluxL[3][NCOMP_TOTAL], FluxR[3][NCOMP_TOTAL];
   real dF[3][NCOMP_TOTAL], Out[NCOMP_TOTAL], Update[NCOMP_TOTAL];
   int  idx_out, idx_in;
   long NCorrThisTime=0;
   bool CorrectUnphy = GAMER_SUCCESS;  // this variable is shared by all OpenMP threads


#  pragma omp parallel for private( VarL, VarC, VarR, FluxL, FluxR, dF, Out, Update, idx_out, idx_in ) schedule( runtime ) \
   reduction( +:NCorrThisTime )
   for (int TID=0; TID<NPG; TID++)
   {
      for (int k_out=0; k_out<PS2; k_out++)
      for (int j_out=0; j_out<PS2; j_out++)
      for (int i_out=0; i_out<PS2; i_out++)
      {
         idx_out = (k_out*PS2 + j_out)*PS2 + i_out;

//       check if the updated values are unphysical
         for (int v=0; v<NCOMP_TOTAL; v++)   Out[v] = h_Flu_Array_F_Out[TID][v][idx_out];

         if ( Unphysical(Out, Gamma_m1, CheckMinPres) )
         {
//          try to correct unphysical results using 1st-order fluxes (which should be more diffusive)
            if ( OPT__1ST_FLUX_CORR )
            {
               idx_in = ( (k_out+FLU_GHOST_SIZE)*FLU_NXT + (j_out+FLU_GHOST_SIZE) )*FLU_NXT + (i_out+FLU_GHOST_SIZE);

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
                  case RSOLVER_ROE:
                     for (int d=0; d<3; d++)
                     {
                        CPU_RiemannSolver_Roe( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        CPU_RiemannSolver_Roe( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

                  case RSOLVER_HLLC:
                     for (int d=0; d<3; d++)
                     {
                        CPU_RiemannSolver_HLLC( d, FluxL[d], VarL[d], VarC,    GAMMA, MIN_PRES );
                        CPU_RiemannSolver_HLLC( d, FluxR[d], VarC,    VarR[d], GAMMA, MIN_PRES );
                     }
                     break;

                  case RSOLVER_HLLE:
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
            } // if ( OPT__1ST_FLUX_CORR )

            else
            {
//             if OPT__1ST_FLUX_CORR is off, just copy data to Update for checking negative density and pressure in the next step
               for (int v=0; v<NCOMP_TOTAL; v++)   Update[v] = Out[v];
            } // if ( OPT__1ST_FLUX_CORR ) ... else ...


//          ensure positive density and pressure
//          --> note that MIN_DENS is declared as double and must be converted to **real** before the comparison
//              --> to be consistent with the check in Unphysical()
            Update[DENS] = FMAX( Update[DENS], (real)MIN_DENS );
            Update[ENGY] = CPU_CheckMinPresInEngy( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                                   Gamma_m1, _Gamma_m1, MIN_PRES );


//          check if the newly updated values are still unphysical
//          --> note that here we check **energy** instead of pressure since even after calling CPU_CheckMinPresInEngy()
//              we can still have pressure < MIN_PRES due to round-off errors (especially when pressure << kinematic energy)
//              --> it will not crash the code since we always apply MIN_PRES when calculating pressure
            if ( Unphysical(Update, Gamma_m1, CheckMinEngy) )
            {
//             output debug information
               const int  PID_Failed      = PID0_List[TID] + LocalID[k_out/PS1][j_out/PS1][i_out/PS1];
               const bool CheckMinPres_No = false;
               real In[NCOMP_TOTAL], tmp[NCOMP_TOTAL];

               char FileName[100];
               sprintf( FileName, "FailedPatchGroup_r%03d_lv%02d_PID0-%05d", MPI_Rank, lv, PID0_List[TID] );

//             use "a" instead of "w" since there may be more than one failed cell in a given patch group
               FILE *File = fopen( FileName, "a" );

               for (int v=0; v<NCOMP_TOTAL; v++)   In[v] = h_Flu_Array_F_In[TID][v][idx_in];

//             output information about the failed cell
               fprintf( File, "PID                              = %5d\n", PID_Failed );
               fprintf( File, "(i,j,k) in the patch             = (%2d,%2d,%2d)\n", i_out%PS1, j_out%PS1, k_out%PS1 );
               fprintf( File, "(i,j,k) in the input fluid array = (%2d,%2d,%2d)\n", i_out, j_out, k_out );
               fprintf( File, "\n" );
               fprintf( File, "input        = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e)\n",
                        In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                        CPU_GetPressure(In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                                        Gamma_m1, CheckMinPres_No, NULL_REAL) );
               fprintf( File, "ouptut (old) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e)\n",
                        Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                        CPU_GetPressure(Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                                        Gamma_m1, CheckMinPres_No, NULL_REAL) );
               fprintf( File, "output (new) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e)\n",
                        Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                        CPU_GetPressure(Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                        Gamma_m1, CheckMinPres_No, NULL_REAL) );

//             output all data in the input fluid array (including ghost zones)
               fprintf( File, "\n===============================================================================================\n" );
               fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s%14s%14s%14s\n",
                        "i", "j", "k", "Density", "Px", "Py", "Pz", "Energy", "Pressure" );

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

//             output the failed patch (mainly for recording the sibling information)
#              ifdef GRAVITY
               Output_Patch( lv, PID_Failed, amr->FluSg[lv], amr->PotSg[lv], "Unphy" );
#              else
               Output_Patch( lv, PID_Failed, amr->FluSg[lv], NULL_INT,       "Unphy" );
#              endif

//             use critical directive to avoid thread racing (may not be necessary here?)
#              pragma omp critical
               CorrectUnphy = GAMER_FAILED;
            } // if ( Unphysical(Update) )

            else
            {
//             store the corrected solution
               for (int v=0; v<NCOMP_TOTAL; v++)   h_Flu_Array_F_Out[TID][v][idx_out] = Update[v];

//             record the number of corrected cells
               NCorrThisTime ++;
            } // if ( Unphysical(Update) ) ... else ...
         } // if need correction
      } // i,j,k
   } // for (int TID=0; TID<NPG; TID++)

// terminate the program if CorrectUnphysical failed
   if ( CorrectUnphy == GAMER_FAILED )
      Aux_Error( ERROR_INFO, "first-order correction failed at Rank %d, lv %d, Time %20.14e, Step %ld, Counter %ld ...\n",
                 MPI_Rank, lv, Time[lv], Step, AdvanceCounter[lv] );

   NCorrUnphy[lv] += NCorrThisTime;

} // FUNCTION : CorrectUnphysical
#endif // #if ( MODEL == HYDRO  ||  MODEL == MHD )



// ============
// |  Tables  |
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  Return the sibling patch ID for the function "CorrectFlux"
//
// Parameter   :  lv    : Targeted refinement level
//                PID   : Targeted patch index
//                SibID : Targeted sibling index (0~5)
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
