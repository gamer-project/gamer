#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

static void StoreFlux( const int lv, const real Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
                       const int NPG, const int *PID0_List );
static void CorrectFlux( const int lv, const real Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
                         const int NPG, const int *PID0_List );
static bool Unphysical( const real Fluid[] );
static void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                               const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT],
                               real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE], const real dt );
static int  Table_01( const int lv, const int PID, const int SibID );

extern void CPU_RiemannSolver_Roe ( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                    const real Gamma );
extern void CPU_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                    const real Gamma );
extern void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                    const real Gamma );
#if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
extern real CPU_PositivePres_In_Engy( const real ConVar[], const real Gamma_m1, const real _Gamma_m1 );
#endif




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
void Flu_Close( const int lv, const int SaveSg, const real h_Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
                real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE], 
                const real h_MinDtInfo_Array[], const int NPG, const int *PID0_List, const bool GetMinDtInfo,
                const real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT], const double dt )
{

// save the flux in the coarse-fine boundary at level "lv"
   if ( OPT__FIXUP_FLUX  &&  lv != TOP_LEVEL )  StoreFlux( lv, h_Flux_Array, NPG, PID0_List );


// correct the flux in the coarse-fine boundary at level "lv-1"
   if ( OPT__FIXUP_FLUX  &&  lv != 0 )    CorrectFlux( lv, h_Flux_Array, NPG, PID0_List );


// correct the unphysical results in h_Flu_Array_F_Out (e.g., negative density in HYDRO/MHD)
   if ( OPT__CORR_UNPHY )     CorrectUnphysical( lv, NPG, PID0_List, h_Flu_Array_F_In, h_Flu_Array_F_Out, dt );


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
void StoreFlux( const int lv, const real h_Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
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


            for (int v=0; v<NFLUX; v++)         {
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
void CorrectFlux( const int lv, const real h_Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
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
      real (*Flux_Temp)[PATCH_SIZE][PATCH_SIZE] = new real [NFLUX][PATCH_SIZE][PATCH_SIZE];


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

               for (int v=0; v<NFLUX; v++)
               for (int m=0; m<PATCH_SIZE; m++)
               for (int n=0; n<PATCH_SIZE; n++)    Flux_Temp[v][m][n] = (real)0.0;

               for (int v=0; v<NFLUX; v++)
               for (int m=0; m<2*PATCH_SIZE; m++)
               for (int n=0; n<2*PATCH_SIZE; n++)
               {  
                  ID = m*2*PATCH_SIZE + n;

                  Flux_Temp[v][m/2][n/2] -= h_Flux_Array[TID][ Mapping[s] ][v][ID];
               }

               for (int v=0; v<NFLUX; v++)
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



//-------------------------------------------------------------------------------------------------------
// Function    :  Unphysical
// Description :  Check whether the input variables are unphysical 
//
// Note        :  1. One can put arbitriry criterion here. Cell violating the conditions will be recalculated
//                   by "CorrectUnphysical"
//                2. Currently it is used for MODEL==HYDRO/MHD to check whether the input density is negative
//                   (and if any variable is -inf, +inf, and nan)
//
// Parameter   :  Input :  Input fluid variable array with size FLU_NOUT 
//
// Return      :  true/false <==> input Fluid array is unphysical/physical
//-------------------------------------------------------------------------------------------------------
bool Unphysical( const real Fluid[] )
{

#  if ( MODEL == HYDRO || MODEL == MHD )
   if ( !isfinite(Fluid[DENS]) || !isfinite(Fluid[MOMX]) || !isfinite(Fluid[MOMY]) ||
        !isfinite(Fluid[MOMZ]) || !isfinite(Fluid[ENGY]) ||
        Fluid[DENS]<=(real)0.0 || Fluid[ENGY]<=(real)0.0 )
      return true;
   else
      return false;

#  else

   return false;
#  endif

} // FUNCTION : Unphysical



//-------------------------------------------------------------------------------------------------------
// Function    :  CorrectUnphysical 
// Description :  Check if any cell in the output array of Fluid solver contains unphysical retult
//                --> For example, negative density
//
// Note        :  1. Define unphysical values in the  function "Unphysical"
//                2. Currently it is used for MODEL==HYDRO/MHD to check whether the input density is negative
//                   (and if any variable is -inf, +inf, and nan)
//                3. But in principle one can define arbitrary criterion to trigger the correction
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

#  if ( MODEL == HYDRO || MODEL == MHD )
   const real dh        = (real)amr->dh[lv];
   const real dt_dh     = dt/dh;
   const int  didx[3]   = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   const real  Gamma_m1 = GAMMA - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif

   real VarL[3][NCOMP], VarC[NCOMP], VarR[3][NCOMP], FluxL[3][NCOMP], FluxR[3][NCOMP], dF[3][NCOMP], Out[NCOMP], Update[NCOMP];
   int  idx_out, idx_in;
   long NCorrThisTime=0;


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
         for (int v=0; v<NCOMP; v++)   Out[v] = h_Flu_Array_F_Out[TID][v][idx_out];

         if ( Unphysical(Out) )
         {
            idx_in = ( (k_out+FLU_GHOST_SIZE)*FLU_NXT + (j_out+FLU_GHOST_SIZE) )*FLU_NXT + (i_out+FLU_GHOST_SIZE);

//          collect nearby input coserved variables
            for (int v=0; v<NCOMP; v++)
            {
               VarC[v] = h_Flu_Array_F_In[TID][v][idx_in];

               for (int d=0; d<3; d++)
               {
                  VarL[d][v] = h_Flu_Array_F_In[TID][v][ idx_in - didx[d] ];
                  VarR[d][v] = h_Flu_Array_F_In[TID][v][ idx_in + didx[d] ];
               }
            }

//          invoke Riemann solver to calculate the fluxes
//          (note that the recalculated flux does NOT include gravity even for UNSPLIT_GRAVITY --> reduce to 1st-order accuracy)
            switch ( OPT__CORR_UNPHY_SCHEME )
            {
               case RSOLVER_ROE:
                  for (int d=0; d<3; d++)
                  {
                     CPU_RiemannSolver_Roe( d, FluxL[d], VarL[d], VarC,    GAMMA );
                     CPU_RiemannSolver_Roe( d, FluxR[d], VarC,    VarR[d], GAMMA );
                  }
                  break;

               case RSOLVER_HLLC:
                  for (int d=0; d<3; d++)
                  {
                     CPU_RiemannSolver_HLLC( d, FluxL[d], VarL[d], VarC,    GAMMA );
                     CPU_RiemannSolver_HLLC( d, FluxR[d], VarC,    VarR[d], GAMMA );
                  }
                  break;

               case RSOLVER_HLLE:
                  for (int d=0; d<3; d++)
                  {
                     CPU_RiemannSolver_HLLE( d, FluxL[d], VarL[d], VarC,    GAMMA );
                     CPU_RiemannSolver_HLLE( d, FluxR[d], VarC,    VarR[d], GAMMA );
                  }
                  break;

               default:
                  Aux_Error( ERROR_INFO, "unnsupported Riemann solver (%d) !!\n", OPT__CORR_UNPHY_SCHEME );
            }

//          recalculate the first-order solution for a full time-step
            for (int d=0; d<3; d++)
            for (int v=0; v<NCOMP; v++)   dF[d][v] = FluxR[d][v] - FluxL[d][v];

            for (int v=0; v<NCOMP; v++)
               Update[v] = h_Flu_Array_F_In[TID][v][idx_in] - dt_dh*( dF[0][v] + dF[1][v] + dF[2][v] );

//          ensure the positive pressure
#           if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
            Update[ENGY] = CPU_PositivePres_In_Engy( Update, Gamma_m1, _Gamma_m1 );
#           endif

//          check if the newly updated values are still unphysical
            if ( Unphysical(Update) )
            {
               Aux_Message( stderr, "ERROR : first-order correction failed at Rank %d, lv %d, Time %20.14e, Step %ld, Counter %ld ...\n",
                            MPI_Rank, lv, Time[lv], Step, AdvanceCounter[lv] );
               Aux_Message( stderr, "   PID0 = %5d, (i,j,k) = (%2d,%2d,%2d)\n", PID0_List[TID], i_out, j_out, k_out );
               Aux_Message( stderr, "   input        = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e) !!\n",
                            h_Flu_Array_F_In [TID][DENS][idx_in],
                            h_Flu_Array_F_In [TID][MOMX][idx_in],
                            h_Flu_Array_F_In [TID][MOMY][idx_in],
                            h_Flu_Array_F_In [TID][MOMZ][idx_in],
                            h_Flu_Array_F_In [TID][ENGY][idx_in] );
               Aux_Message( stderr, "   ouptut (old) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e) !!\n",
                            h_Flu_Array_F_Out[TID][DENS][idx_out],
                            h_Flu_Array_F_Out[TID][MOMX][idx_out],
                            h_Flu_Array_F_Out[TID][MOMY][idx_out],
                            h_Flu_Array_F_Out[TID][MOMZ][idx_out],
                            h_Flu_Array_F_Out[TID][ENGY][idx_out] );
               Aux_Message( stderr, "   output (new) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e) !!\n",
                           Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY] );

//             output all patches in this patch group
#              ifdef GRAVITY
               for (int PID=PID0_List[TID]; PID<PID0_List[TID]+8; PID++)   Output_Patch( lv, PID, amr->FluSg[lv], amr->PotSg[lv], "Unphy" );
#              else
               for (int PID=PID0_List[TID]; PID<PID0_List[TID]+8; PID++)   Output_Patch( lv, PID, amr->FluSg[lv],       NULL_INT, "Unphy" );
#              endif

               MPI_Exit();
            } // if ( Unphysical(Update) )

//          store the update solution
            for (int v=0; v<NCOMP; v++)   h_Flu_Array_F_Out[TID][v][idx_out] = Update[v];

//          record the number of corrected cells
            NCorrThisTime ++;
         } // if need correction
      } // i,j,k
   } // for (int TID=0; TID<NPG; TID++)

   NCorrUnphy[lv] += NCorrThisTime;
#  endif // HYDRO/MHD

} // FUNCTION : CorrectUnphysical



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
