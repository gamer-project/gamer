#include "GAMER.h"

static void StoreFlux( const int lv, const real Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
                       const int NPG, const int *PID0_List );
static void CorrectFlux( const int lv, const real Flux_Array[][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE],
                         const int NPG, const int *PID0_List );
static int  Table_01( const int lv, const int PID, const int SibID );




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
                const real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE], 
                const real h_MinDtInfo_Array[], const int NPG, const int *PID0_List, const bool GetMinDtInfo )
{

// save the flux in the coarse-fine boundary at level "lv"
   if ( OPT__FIXUP_FLUX  &&  lv != TOP_LEVEL )  StoreFlux( lv, h_Flux_Array, NPG, PID0_List );


// correct the flux in the coarse-fine boundary at level "lv-1"
   if ( OPT__FIXUP_FLUX  &&  lv != 0 )          CorrectFlux( lv, h_Flux_Array, NPG, PID0_List );


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

//          record the unwrapped phase update (= unwrapped phase - initial phase) 
            double Phase_Old, Phase_New;

            Phase_Old = ATAN2( amr->patch[1-SaveSg][lv][PID]->fluid[IMAG][k][j][i], 
                               amr->patch[1-SaveSg][lv][PID]->fluid[REAL][k][j][i] );
            Phase_New = ATAN2( amr->patch[  SaveSg][lv][PID]->fluid[IMAG][k][j][i], 
                               amr->patch[  SaveSg][lv][PID]->fluid[REAL][k][j][i] );

            if ( Phase_New > Phase_Old )  Phase_New -= 2.0*M_PI;

            amr->patch[SaveSg][lv][PID]->fluid[3][k][j][i] = amr->patch[1-SaveSg][lv][PID]->fluid[3][k][j][i] 
                                                             + Phase_New - Phase_Old;
         }}}}
      }
   } // for (int TID=0; TID<NPG; TID++)


// store the minimum time-step estimated by the hydrodynamical CFL condition
   /*
   if ( GetMinDtInfo )
   {
#     ifdef OOC
      if ( ooc.Rank == 0 )
#     endif
      if ( PID_Start == 0 )   MinDtInfo_Fluid[lv] = (real)0.0;

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
