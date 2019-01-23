#include "GAMER.h"

#ifdef LOAD_BALANCE


const real BufSizeFactor = 1.2;  // Send/RecvBufSize = int(NSend/NRecv*BufSizeFactor) --> must be >= 1.0

// MPI buffers are shared by some particle routines
static real *MPI_SendBuf_Shared = NULL;
static real *MPI_RecvBuf_Shared = NULL;
static int   SendBufSize        = -1;
static int   RecvBufSize        = -1;

#ifdef TIMING
extern Timer_t *Timer_MPI[3];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetBufferData
// Description :  Fill up the data(flux) of buffer(real) patches
//
// Note        :  1. This function is the alternative to "Buf_GetBufferData" in LOAD_BALANCE, and is
//                   invoked by "Buf_GetBufferData" and "Flu_FixUp"
//                2. The modes "POT_FOR_POISSON" and "POT_AFTER_REFINE" can be applied to the potential
//                   data only. For others modes, the number of variables to be exchanged depends on
//                   the input parameter "TVar".
//
// Parameter   :  lv          : Target refinement level to exchage data
//                FluSg       : Sandglass of the requested fluid data (useless in POT_FOR_POISSON,
//                              POT_AFTER_REFINEPOT, COARSE_FINE_FLUX )
//                PotSg       : Sandglass of the requested potential data (useless in COARSE_FINE_FLUX)
//                GetBufMode  : Target mode. Each mode has its own MPI lists, by which the amount of data
//                              to be transferred can be minimized.
//                              --> DATA_GENERAL      : data for general-purpose (sibling and coarse-grid data)
//                                  DATA_AFTER_REFINE : subset of DATA_GENERAL after refine
//                                  DATA_AFTER_FIXUP  : subset of DATA_GENERAL after fix-up
//                                  DATA_RESTRICT     : restricted data of the father patches with sons not home
//                                  POT_FOR_POISSON   : potential for the Poisson solver
//                                  POT_AFTER_REFINE  : potential after refine for the Poisson solver
//                                  COARSE_FINE_FLUX  : fluxes across the coarse-fine boundaries (HYDRO ONLY)
//                TVar        : Target variables to exchange
//                              --> Supported variables in different models:
//                                  HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY,[, _POTE]
//                                  MHD   :
//                                  ELBDM : _DENS, _REAL, _IMAG, [, _POTE]
//                              --> _FLUID, _PASSIVE, and _TOTAL apply to all models
//                              --> In addition, the flux variables (e.g., _FLUX_DENS) are also supported
//                              Restrictions :
//                              --> a. DATA_XXX works with all components in (_TOTAL | _POTE)
//                                  b. COARSE_FINE_FLUX works with all components in (_FLUX_TOTAL)
//                                  c. _POTE has no effect on the flux fix-up in DATA_AFTER_FIXUP
//                                  d. POT_FOR_POISSON and POT_AFTER_REFINE only work with _POTE
//                ParaBuf     : Number of ghost zones to exchange (useless in DATA_RESTRICT and COARSE_FINE_FLUX )
//-------------------------------------------------------------------------------------------------------
void LB_GetBufferData( const int lv, const int FluSg, const int PotSg, const GetBufMode_t GetBufMode,
                       const int TVar, const int ParaBuf )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if (  GetBufMode != COARSE_FINE_FLUX  &&  ( TVar & _TOTAL )  &&  ( FluSg != 0 && FluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FluSg", FluSg );

#  ifdef GRAVITY
   if (  GetBufMode != COARSE_FINE_FLUX  &&  ( TVar & _POTE )  &&  ( PotSg != 0 && PotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PotSg", PotSg );

   if (  ( GetBufMode == DATA_GENERAL || GetBufMode == DATA_AFTER_FIXUP || GetBufMode == DATA_AFTER_REFINE ||
           GetBufMode == DATA_RESTRICT )  &&  !( TVar & (_TOTAL|_POTE) )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing (_TOTAL|_POTE) !!\n" );

   if (  ( GetBufMode == POT_FOR_POISSON || GetBufMode == POT_AFTER_REFINE )  &&  !( TVar & _POTE )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing _POTE !!\n" );

   if (  ( GetBufMode == POT_FOR_POISSON || GetBufMode == POT_AFTER_REFINE )  &&  ( TVar & ~_POTE )  )
      Aux_Error( ERROR_INFO, "modes \"%s\" only accept \"%s\" as the target variable !!\n",
                 "POT_FOR_POISSON and POT_AFTER_REFINE", "_POTE" );

   if (  ( GetBufMode == DATA_GENERAL || GetBufMode == DATA_AFTER_FIXUP || GetBufMode == DATA_AFTER_REFINE ||
           GetBufMode == POT_FOR_POISSON || GetBufMode == POT_AFTER_REFINE )  &&
         ( ParaBuf < 0 || ParaBuf > PATCH_SIZE )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d --> accepted range = [0 ... PATCH_SIZE] !!\n",
                 "ParaBuf", ParaBuf );

#  else // #ifdef GRAVITY ... else ...
   if (  ( GetBufMode == DATA_GENERAL || GetBufMode == DATA_AFTER_FIXUP || GetBufMode == DATA_AFTER_REFINE ||
           GetBufMode == DATA_RESTRICT )  &&  !( TVar & _TOTAL )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing _TOTAL !!\n" );

   if (  ( GetBufMode == DATA_GENERAL || GetBufMode == DATA_AFTER_FIXUP || GetBufMode == DATA_AFTER_REFINE )  &&
         ( ParaBuf < 0 || ParaBuf > PATCH_SIZE )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d --> range=[0 ... PATCH_SIZE] !!\n",
                 "ParaBuf", ParaBuf );
#  endif // #ifdef GRAVITY ... else ...

   if (  GetBufMode == COARSE_FINE_FLUX  &&  !( TVar & _FLUX_TOTAL )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing _FLUX_TOTAL !!\n" );

   if ( GetBufMode == COARSE_FINE_FLUX  &&  !amr->WithFlux )
   {
      Aux_Message( stderr, "WARNING : mode COARSE_FINE_FLUX is useless since no flux is required !!\n" );
      return;
   }


// determine the components to be prepared (TFluVarIdx : target fluid variable indices ( = [0 ... NCOMP_TOTAL-1/NFLUX_TOTAL-1] )
   bool ExchangeFlu = ( GetBufMode == COARSE_FINE_FLUX ) ?
                      TVar & _FLUX_TOTAL : TVar & _TOTAL;                        // whether or not to exchage the fluid data
#  ifdef GRAVITY
   bool ExchangePot = (  GetBufMode != COARSE_FINE_FLUX  &&  (TVar & _POTE)  );  // whether or not to exchange the potential data
#  endif

   const int NFluid_Max = ( GetBufMode == COARSE_FINE_FLUX ) ? NFLUX_TOTAL : NCOMP_TOTAL;
   int NVar_Flu, NVar_Tot, TFluVarIdx, *TFluVarIdxList;
   TFluVarIdxList = new int [NFluid_Max];
   NVar_Flu = 0;

   for (int v=0; v<NFluid_Max; v++)
      if ( TVar & (1<<v) )    TFluVarIdxList[ NVar_Flu++ ] = v;

   NVar_Tot = NVar_Flu;
#  ifdef GRAVITY
   if ( ExchangePot )   NVar_Tot ++;

   if ( GetBufMode == POT_FOR_POISSON  ||  GetBufMode == POT_AFTER_REFINE )
   {
      ExchangeFlu = false;
      ExchangePot = true;
      NVar_Flu    = 0;
      NVar_Tot    = 1;
   }
#  endif

// check again
   if ( NVar_Tot == 0  ||  ( GetBufMode == COARSE_FINE_FLUX && NVar_Flu == 0 )  )
   {
      Aux_Message( stderr, "WARNING : no target variable is found !!\n" );
      return;
   }


// _X : for exchanging data after the FLUX fix-up, which has ParaBuf=1
   const int DataUnit_Flux = PS1*PS1*NVar_Flu;

   int   NSend_Total, NRecv_Total, SPID, RPID, SSib, RSib, Counter;
   int   DataUnit_Buf[27], LoopStart[27][3], LoopEnd[27][3];
   int   LoopStart_X[6][3], LoopEnd_X[6][3];
   real *SendPtr=NULL, *RecvPtr=NULL;
   int  *Send_NList=NULL, *Recv_NList=NULL, *Send_NResList=NULL, *Recv_NResList=NULL;
   int **Send_IDList=NULL, **Recv_IDList=NULL, **Send_IDList_IdxTable=NULL, **Recv_IDList_IdxTable=NULL;
   int **Send_SibList=NULL, **Recv_SibList=NULL;
   real (*FluxPtr)[PS1][PS1]=NULL;

   int *Send_NCount = new int [MPI_NRank];
   int *Recv_NCount = new int [MPI_NRank];
   int *Send_NDisp  = new int [MPI_NRank];
   int *Recv_NDisp  = new int [MPI_NRank];


// 1. set up the number of elements to be sent and received in each cell and the send/recv lists
// ============================================================================================================
   switch ( GetBufMode )
   {
      case DATA_GENERAL :        // general-purpose
         Send_NList           = amr->LB->SendH_NList          [lv];
         Send_IDList          = amr->LB->SendH_IDList         [lv];
         Send_SibList         = amr->LB->SendH_SibList        [lv];
         Recv_NList           = amr->LB->RecvH_NList          [lv];
         Recv_IDList          = amr->LB->RecvH_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvH_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvH_SibList        [lv];
         break;

      case DATA_AFTER_REFINE :   // data after refine
         Send_NList           = amr->LB->SendH_NList          [lv];
         Send_IDList          = amr->LB->SendH_IDList         [lv];
         Send_SibList         = amr->LB->SendH_SibDiffList    [lv];
         Recv_NList           = amr->LB->RecvH_NList          [lv];
         Recv_IDList          = amr->LB->RecvH_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvH_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvH_SibDiffList    [lv];
         break;

      case DATA_AFTER_FIXUP :    // data after fix-up
         Send_NList           = amr->LB->SendX_NList          [lv];
         Send_NResList        = amr->LB->SendX_NResList       [lv];
         Send_IDList          = amr->LB->SendX_IDList         [lv];
         Send_SibList         = amr->LB->SendX_SibList        [lv];
         Recv_NList           = amr->LB->RecvX_NList          [lv];
         Recv_NResList        = amr->LB->RecvX_NResList       [lv];
         Recv_IDList          = amr->LB->RecvX_IDList         [lv];
         Recv_SibList         = amr->LB->RecvX_SibList        [lv];
         break;

      case DATA_RESTRICT :       // restricted data
         Send_NList           = amr->LB->SendR_NList          [lv];
         Send_IDList          = amr->LB->SendR_IDList         [lv];
         Send_IDList_IdxTable = amr->LB->SendR_IDList_IdxTable[lv];
         Recv_NList           = amr->LB->RecvR_NList          [lv];
         Recv_IDList          = amr->LB->RecvR_IDList         [lv];
         break;

#     ifdef GRAVITY
      case POT_FOR_POISSON :     // potential for the Poisson solver
         Send_NList           = amr->LB->SendG_NList          [lv];
         Send_IDList          = amr->LB->SendG_IDList         [lv];
         Send_SibList         = amr->LB->SendG_SibList        [lv];
         Recv_NList           = amr->LB->RecvG_NList          [lv];
         Recv_IDList          = amr->LB->RecvG_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvG_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvG_SibList        [lv];
         break;

      case POT_AFTER_REFINE :    // potential after refine for the Poisson solver
         Send_NList           = amr->LB->SendG_NList          [lv];
         Send_IDList          = amr->LB->SendG_IDList         [lv];
         Send_SibList         = amr->LB->SendG_SibDiffList    [lv];
         Recv_NList           = amr->LB->RecvG_NList          [lv];
         Recv_IDList          = amr->LB->RecvG_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvG_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvG_SibDiffList    [lv];
         break;
#     endif // #ifdef GRAVITY

      case COARSE_FINE_FLUX :    // hydro fluxes
         Send_NList           = amr->LB->SendF_NList          [lv];
         Send_IDList          = amr->LB->SendF_IDList         [lv];
         Send_SibList         = amr->LB->SendF_SibList        [lv];
         Recv_NList           = amr->LB->RecvF_NList          [lv];
         Recv_IDList          = amr->LB->RecvF_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvF_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvF_SibList        [lv];
         break;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GetBufMode", GetBufMode );
   } // switch ( GetBufMode )



// 2. set up the loop range, MPI count and displacement arrays, and allocate data for send and recv buffers
// ============================================================================================================
   switch ( GetBufMode )
   {
      case DATA_GENERAL: case DATA_AFTER_REFINE:
#     ifdef GRAVITY
      case POT_FOR_POISSON : case POT_AFTER_REFINE:
#     endif
//    ----------------------------------------------
//       loop range
         for (int s=0; s<26; s++)
         for (int d=0; d<3; d++)
         {
            LoopStart[s][d] = TABLE_01( s, 'x'+d, 0,       0,   PS1-ParaBuf );
            LoopEnd  [s][d] = TABLE_01( s, 'x'+d, ParaBuf, PS1, PS1         );
         }

         LoopStart[26][0] = LoopStart[26][1] = LoopStart[26][2] = 0;
         LoopEnd  [26][0] = LoopEnd  [26][1] = LoopEnd  [26][2] = PS1;

//       MPI count array
         for (int s= 0; s< 6; s++)  DataUnit_Buf[ s] = NVar_Tot*PS1    *PS1    *ParaBuf;  // sibling 0 ~ 5
         for (int s= 6; s<18; s++)  DataUnit_Buf[ s] = NVar_Tot*PS1    *ParaBuf*ParaBuf;  // sibling 6 ~ 17
         for (int s=18; s<26; s++)  DataUnit_Buf[ s] = NVar_Tot*ParaBuf*ParaBuf*ParaBuf;  // sibling 18 ~ 25
                                    DataUnit_Buf[26] = NVar_Tot*PS1    *PS1    *PS1;      // entire patch

         for (int r=0; r<MPI_NRank; r++)
         {
            Send_NCount[r] = 0;
            Recv_NCount[r] = 0;

            for (int t=0; t<Send_NList[r]; t++)
            {
               if ( Send_SibList[r][t] != 0 )
               for (int s=0; s<27; s++)
                  if ( Send_SibList[r][t] & (1<<s) )  Send_NCount[r] += DataUnit_Buf[s];
            }

            for (int t=0; t<Recv_NList[r]; t++)
            {
               if ( Recv_SibList[r][t] != 0 )
               for (int s=0; s<27; s++)
                  if ( Recv_SibList[r][t] & (1<<s) )  Recv_NCount[r] += DataUnit_Buf[s];
            }
         }
         break; // case DATA_GENERAL: case DATA_AFTER_REFINE: case POT_FOR_POISSON : case POT_AFTER_REFINE:


      case DATA_AFTER_FIXUP :
//    ----------------------------------------------
//       loop range for restriction fix-up
         for (int s=0; s<26; s++)
         for (int d=0; d<3; d++)
         {
            LoopStart[s][d] = TABLE_01( s, 'x'+d, 0,       0,   PS1-ParaBuf );
            LoopEnd  [s][d] = TABLE_01( s, 'x'+d, ParaBuf, PS1, PS1         );
         }

         LoopStart[26][0] = LoopStart[26][1] = LoopStart[26][2] = 0;
         LoopEnd  [26][0] = LoopEnd  [26][1] = LoopEnd  [26][2] = PS1;

//       loop range for flux fix-up (with ParaBuf always == 1)
         for (int s=0; s<6; s++)
         for (int d=0; d<3; d++)
         {
            LoopStart_X[s][d] = TABLE_01( s, 'x'+d, 0, 0,   PS1-1 );
            LoopEnd_X  [s][d] = TABLE_01( s, 'x'+d, 1, PS1, PS1   );
         }

//       MPI count array
         for (int s= 0; s< 6; s++)  DataUnit_Buf[ s] = NVar_Tot*PS1    *PS1    *ParaBuf;  // sibling 0 ~ 5
         for (int s= 6; s<18; s++)  DataUnit_Buf[ s] = NVar_Tot*PS1    *ParaBuf*ParaBuf;  // sibling 6 ~ 17
         for (int s=18; s<26; s++)  DataUnit_Buf[ s] = NVar_Tot*ParaBuf*ParaBuf*ParaBuf;  // sibling 18 ~ 25
                                    DataUnit_Buf[26] = NVar_Tot*PS1    *PS1    *PS1;      // entire patch

         for (int r=0; r<MPI_NRank; r++)
         {
            Send_NCount[r] = 0;
            Recv_NCount[r] = 0;

//          for restriction fix-up
            for (int t=0; t<Send_NResList[r]; t++)
            for (int s=0; s<27; s++)
               if ( Send_SibList[r][t] & (1<<s) )  Send_NCount[r] += DataUnit_Buf[s];

            for (int t=0; t<Recv_NResList[r]; t++)
            for (int s=0; s<27; s++)
               if ( Recv_SibList[r][t] & (1<<s) )  Recv_NCount[r] += DataUnit_Buf[s];

//          for flux fix-up
            for (int t=Send_NResList[r]; t<Send_NList[r]; t++)
            for (int s=0; s<6; s++)
               if ( Send_SibList[r][t] & (1<<s) )  Send_NCount[r] += DataUnit_Flux;

            for (int t=Recv_NResList[r]; t<Recv_NList[r]; t++)
            for (int s=0; s<6; s++)
               if ( Recv_SibList[r][t] & (1<<s) )  Recv_NCount[r] += DataUnit_Flux;
         }
         break; // case DATA_AFTER_FIXUP :


      case DATA_RESTRICT :
//    ----------------------------------------------
         for (int r=0; r<MPI_NRank; r++)
         {
            Send_NCount[r] = Send_NList[r]*PS1*PS1*PS1*NVar_Tot;
            Recv_NCount[r] = Recv_NList[r]*PS1*PS1*PS1*NVar_Tot;
         }
         break; // case DATA_RESTRICT :


      case COARSE_FINE_FLUX :
//    ----------------------------------------------
         for (int r=0; r<MPI_NRank; r++)
         {
            Send_NCount[r] = Send_NList[r]*DataUnit_Flux;
            Recv_NCount[r] = Recv_NList[r]*DataUnit_Flux;
         }
         break; // COARSE_FINE_FLUX
   } // switch ( GetBufMode )


// MPI displacement array
   Send_NDisp[0] = 0;
   Recv_NDisp[0] = 0;

   for (int r=1; r<MPI_NRank; r++)
   {
      Send_NDisp[r] = Send_NDisp[r-1] + Send_NCount[r-1];
      Recv_NDisp[r] = Recv_NDisp[r-1] + Recv_NCount[r-1];
   }

   NSend_Total = Send_NDisp[ MPI_NRank-1 ] + Send_NCount[ MPI_NRank-1 ];
   NRecv_Total = Recv_NDisp[ MPI_NRank-1 ] + Recv_NCount[ MPI_NRank-1 ];


// allocate send/recv buffers (only when the current buffer size is not large enough --> improve performance)
   real *SendBuf = LB_GetBufferData_MemAllocate_Send( NSend_Total );
   real *RecvBuf = LB_GetBufferData_MemAllocate_Recv( NRecv_Total );



// 3. prepare the send array
// ============================================================================================================
#  ifdef TIMING
   if ( OPT__TIMING_MPI )  Timer_MPI[0]->Start();
#  endif

   switch ( GetBufMode )
   {
      case DATA_GENERAL: case DATA_AFTER_REFINE:
#     ifdef GRAVITY
      case POT_FOR_POISSON : case POT_AFTER_REFINE:
#     endif
//    ----------------------------------------------
#        pragma omp parallel for private( SendPtr, Counter, SPID, SSib, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            SendPtr = SendBuf + Send_NDisp[r];
            Counter = 0;

            for (int t=0; t<Send_NList[r]; t++)
            {
               SPID = Send_IDList [r][t];    // both SPID and SSib are sorted
               SSib = Send_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu )
               {
                  if ( amr->patch[FluSg][lv][SPID]->fluid == NULL )
                     Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                                GetBufMode, FluSg, lv, SPID );

                  if ( SSib == 0  &&  GetBufMode != DATA_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                                GetBufMode, t, r, SPID );
               }

#              ifdef GRAVITY
               if ( ExchangePot )
               {
                  if ( amr->patch[PotSg][lv][SPID]->pot == NULL )
                     Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->pot has not been allocated !!\n",
                                GetBufMode, PotSg, lv, SPID );

                  if ( SSib == 0  &&  GetBufMode != POT_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                                GetBufMode, t, r, SPID );
               }
#              endif // #ifdef GRAVITY
#              endif // #ifdef GAMER_DEBUG

               if ( SSib )
               for (int s=0; s<27; s++)
               {
                  if ( SSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVar_Flu; v++)
                     {
                        TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][k][j][i];
                     }

#                    ifdef GRAVITY
//                   potential data
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[PotSg][lv][SPID]->pot[k][j][i];
                     }
#                    endif
                  } // if ( SSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_GENERAL: case DATA_AFTER_REFINE: case POT_FOR_POISSON : case POT_AFTER_REFINE:


      case DATA_AFTER_FIXUP :
//    ----------------------------------------------
#        pragma omp parallel for private( SendPtr, Counter, SPID, SSib, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            SendPtr = SendBuf + Send_NDisp[r];
            Counter = 0;

//          for restriction fix-up
            for (int t=0; t<Send_NResList[r]; t++)
            {
               SPID = Send_IDList [r][t];
               SSib = Send_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][SPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, SPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][SPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, SPID );
#              endif

               if ( SSib == 0 )
                  Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                             GetBufMode, t, r, SPID );
#              endif // #ifdef GAMER_DEBUG

               for (int s=0; s<27; s++)
               {
                  if ( SSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVar_Flu; v++)
                     {
                        TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][k][j][i];
                     }

#                    ifdef GRAVITY
//                   potential data
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[PotSg][lv][SPID]->pot[k][j][i];
                     }
#                    endif
                  } // if ( SSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Send_NResList[r]; t++)


//          for flux fix-up
            for (int t=Send_NResList[r]; t<Send_NList[r]; t++)
            {
               SPID = Send_IDList [r][t];    // both SPID and SSib are sorted
               SSib = Send_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][SPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, SPID );

               if ( SSib == 0 )
                  Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                             GetBufMode, t, r, SPID );
#              endif // #ifdef GAMER_DEBUG

               for (int s=0; s<6; s++)
               {
                  if ( SSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVar_Flu; v++)
                     {
                        TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart_X[s][2]; k<LoopEnd_X[s][2]; k++)
                        for (int j=LoopStart_X[s][1]; j<LoopEnd_X[s][1]; j++)
                        for (int i=LoopStart_X[s][0]; i<LoopEnd_X[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][k][j][i];
                     }
                  } // if ( SSib & (1<<s) )
               } // for (int s=0; s<6; s++)
            } //for (int t=Send_NResList[r]; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_AFTER_FIXUP :


      case DATA_RESTRICT :
//    ----------------------------------------------
#        pragma omp parallel for private( SendPtr, SPID, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            SendPtr = SendBuf + Send_NDisp[r];

            for (int t=0; t<Send_NList[r]; t++)
            {
               SPID = Send_IDList[r][ Send_IDList_IdxTable[r][t] ];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][SPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, SPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][SPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, ptr[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, SPID );
#              endif
#              endif // #ifdef GAMER_DEBUG

//             fluid data
               if ( ExchangeFlu )
               for (int v=0; v<NVar_Flu; v++)
               {
                  TFluVarIdx = TFluVarIdxList[v];

                  memcpy( SendPtr, &amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][0][0][0],
                          PS1*PS1*PS1*sizeof(real) );

                  SendPtr += PS1*PS1*PS1;
               }

#              ifdef GRAVITY
//             potential data
               if ( ExchangePot )
               {
                  memcpy( SendPtr, &amr->patch[PotSg][lv][SPID]->pot[0][0][0],
                          PS1*PS1*PS1*sizeof(real) );

                  SendPtr += PS1*PS1*PS1;
               }
#              endif
            } // for (int t=0; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_RESTRICT :


      case COARSE_FINE_FLUX :
//    ----------------------------------------------
#        pragma omp parallel for private( SendPtr, Counter, SPID, SSib, FluxPtr, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            SendPtr = SendBuf + Send_NDisp[r];
            Counter = 0;

            for (int t=0; t<Send_NList[r]; t++)
            {
               SPID    = Send_IDList [r][t];
               SSib    = Send_SibList[r][t];
               FluxPtr = amr->patch[0][lv][SPID]->flux[SSib];

#              ifdef GAMER_DEBUG
               if ( FluxPtr == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, ptr[0][%d][%d]->flux[%d] has not been allocated !!\n",
                             GetBufMode, lv, SPID, SSib );
#              endif

               for (int v=0; v<NVar_Flu; v++)
               {
                  TFluVarIdx = TFluVarIdxList[v];

                  memcpy( SendPtr, FluxPtr[TFluVarIdx], PS1*PS1*sizeof(real) );

                  SendPtr += PS1*PS1;
               }
            } // for (int t=0; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case COARSE_FINE_FLUX :


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GetBufMode", GetBufMode );
   } // switch ( GetBufMode )

#  ifdef TIMING
   if ( OPT__TIMING_MPI )  Timer_MPI[0]->Stop();
#  endif



// 4. transfer data by MPI_Alltoallv
// ============================================================================================================
#  ifdef TIMING
// it's better to add barrier before timing transferring data through MPI
// --> so that the timing results (i.e., the MPI bandwidth reported by OPT__TIMING_MPI ) does NOT include
//     the time waiting for other ranks to reach here
// --> make the MPI bandwidth measured here more accurate
   if ( OPT__TIMING_BARRIER )    MPI_Barrier( MPI_COMM_WORLD );

   if ( OPT__TIMING_MPI )  Timer_MPI[1]->Start();
#  endif

#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf, Send_NCount, Send_NDisp, MPI_DOUBLE,
                  RecvBuf, Recv_NCount, Recv_NDisp, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf, Send_NCount, Send_NDisp, MPI_FLOAT,
                  RecvBuf, Recv_NCount, Recv_NDisp, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif

#  ifdef TIMING
   if ( OPT__TIMING_MPI )  Timer_MPI[1]->Stop();
#  endif



// 5. store the received data to their corresponding patches
// ============================================================================================================
#  ifdef TIMING
   if ( OPT__TIMING_MPI )  Timer_MPI[2]->Start();
#  endif

   switch ( GetBufMode )
   {
      case DATA_GENERAL: case DATA_AFTER_REFINE:
#     ifdef GRAVITY
      case POT_FOR_POISSON : case POT_AFTER_REFINE:
#     endif
//    ----------------------------------------------
#        pragma omp parallel for private( RecvPtr, Counter, RPID, RSib, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            RecvPtr = RecvBuf + Recv_NDisp[r];
            Counter = 0;

            for (int t=0; t<Recv_NList[r]; t++)
            {
               RPID = Recv_IDList [r][ Recv_IDList_IdxTable[r][t] ]; // Recv_IDList is unsorted
               RSib = Recv_SibList[r][t];                            // Recv_SibList is sorted

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu )
               {
                  if ( amr->patch[FluSg][lv][RPID]->fluid == NULL )
                     Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                                GetBufMode, FluSg, lv, RPID );

                  if ( RSib == 0  &&  GetBufMode != DATA_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                                GetBufMode, t, r, RPID );
               }

#              ifdef GRAVITY
               if ( ExchangePot )
               {
                  if ( amr->patch[PotSg][lv][RPID]->pot == NULL )
                     Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->pot has not been allocated !!\n",
                                GetBufMode, PotSg, lv, RPID );

                  if ( RSib == 0  &&  GetBufMode != POT_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                                GetBufMode, t, r, RPID );
               }
#              endif // #ifdef GRAVITY
#              endif // #ifdef GAMER_DEBUG

               if ( RSib )
               for (int s=0; s<27; s++)
               {
                  if ( RSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVar_Flu; v++)
                     {
                        TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][k][j][i] = RecvPtr[ Counter ++ ];
                     }

#                    ifdef GRAVITY
//                   potential data
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[PotSg][lv][RPID]->pot[k][j][i] = RecvPtr[ Counter ++ ];
                     }
#                    endif
                  } // if ( RSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_GENERAL: case DATA_AFTER_REFINE: case POT_FOR_POISSON : case POT_AFTER_REFINE:


      case DATA_AFTER_FIXUP :
//    ----------------------------------------------
#        pragma omp parallel for private( RecvPtr, Counter, RPID, RSib, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            RecvPtr = RecvBuf + Recv_NDisp[r];
            Counter = 0;

//          for restriction fix-up
            for (int t=0; t<Recv_NResList[r]; t++)
            {
               RPID = Recv_IDList [r][t];
               RSib = Recv_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][RPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, RPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][RPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, RPID );
#              endif

               if ( RSib == 0 )
                  Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                             GetBufMode, t, r, RPID );
#              endif // #ifdef GAMER_DEBUG

               for (int s=0; s<27; s++)
               {
                  if ( RSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVar_Flu; v++)
                     {
                        TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][k][j][i] = RecvPtr[ Counter ++ ];
                     }

#                    ifdef GRAVITY
//                   potential data
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[PotSg][lv][RPID]->pot[k][j][i] = RecvPtr[ Counter ++ ];
                     }
#                    endif
                  } // if ( RSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Recv_NResList[r]; t++)


//          for flux fix-up
            for (int t=Recv_NResList[r]; t<Recv_NList[r]; t++)
            {
               RPID = Recv_IDList [r][t];
               RSib = Recv_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][RPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, RPID );

               if ( RSib == 0 )
                  Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                             GetBufMode, t, r, RPID );
#              endif // #ifdef GAMER_DEBUG

               for (int s=0; s<6; s++)
               {
                  if ( RSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVar_Flu; v++)
                     {
                        TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart_X[s][2]; k<LoopEnd_X[s][2]; k++)
                        for (int j=LoopStart_X[s][1]; j<LoopEnd_X[s][1]; j++)
                        for (int i=LoopStart_X[s][0]; i<LoopEnd_X[s][0]; i++)
                           amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][k][j][i] = RecvPtr[ Counter ++ ];
                     }
                  } // if ( RSib & (1<<s) )
               } // for (int s=0; s<6; s++)
            } // for (int t=Recv_NResList[r]; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_AFTER_FIXUP :


      case DATA_RESTRICT :
//    ----------------------------------------------
#        pragma omp parallel for private( RecvPtr, RPID, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            RecvPtr = RecvBuf + Recv_NDisp[r];

            for (int t=0; t<Recv_NList[r]; t++)
            {
               RPID = Recv_IDList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][RPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, RPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][RPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, ptr[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, RPID );
#              endif
#              endif // #ifdef GAMER_DEBUG

//             fluid data
               if ( ExchangeFlu )
               for (int v=0; v<NVar_Flu; v++)
               {
                  TFluVarIdx = TFluVarIdxList[v];

                  memcpy( &amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][0][0][0], RecvPtr,
                          PS1*PS1*PS1*sizeof(real) );

                  RecvPtr += PS1*PS1*PS1;
               }

#              ifdef GRAVITY
//             potential data
               if ( ExchangePot )
               {
                  memcpy( &amr->patch[PotSg][lv][RPID]->pot[0][0][0], RecvPtr,
                          PS1*PS1*PS1*sizeof(real) );

                  RecvPtr += PS1*PS1*PS1;
               }
#              endif
            } // for (int t=0; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_RESTRICT :


      case COARSE_FINE_FLUX :
//    ----------------------------------------------
#        pragma omp parallel for private( RecvPtr, Counter, RPID, RSib, FluxPtr, TFluVarIdx ) schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            RecvPtr = RecvBuf + Recv_NDisp[r];
            Counter = 0;

            for (int t=0; t<Recv_NList[r]; t++)
            {
               RPID    = Recv_IDList [r][ Recv_IDList_IdxTable[r][t] ];
               RSib    = Recv_SibList[r][t];
               FluxPtr = amr->patch[0][lv][RPID]->flux[RSib];

#              ifdef GAMER_DEBUG
               if ( FluxPtr == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, ptr[0][%d][%d]->flux[%d] has not been allocated !!\n",
                             GetBufMode, lv, RPID, RSib );
#              endif

//             add (not replace) flux array with the received flux
               for (int v=0; v<NVar_Flu; v++)
               {
                  TFluVarIdx = TFluVarIdxList[v];

                  for (int m=0; m<PS1; m++)
                  for (int n=0; n<PS1; n++)
                     FluxPtr[TFluVarIdx][m][n] += RecvPtr[ Counter ++ ];
               }
            } // for (int t=0; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case COARSE_FINE_FLUX :


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GetBufMode", GetBufMode );
   } // switch ( GetBufMode )

#  ifdef TIMING
   if ( OPT__TIMING_MPI )  Timer_MPI[2]->Stop();
#  endif



// 6. record the achieved MPI bandwidth
// ============================================================================================================
#  ifdef TIMING
   if ( OPT__TIMING_MPI )
   {
      char FileName[100];
      sprintf( FileName, "Record__TimingMPI_Rank%05d", MPI_Rank );

      char ModeName[100];
      switch ( GetBufMode )
      {
         case DATA_GENERAL :        sprintf( ModeName, "%s", (NVar_Tot==1)?"Dens":"Flu" );   break;
         case DATA_AFTER_REFINE :   sprintf( ModeName, "%s", "FluAfRef"                 );   break;
         case DATA_AFTER_FIXUP :    sprintf( ModeName, "%s", "FluAfFix"                 );   break;
         case DATA_RESTRICT :       sprintf( ModeName, "%s", "FluRes"                   );   break;
#        ifdef GRAVITY
         case POT_FOR_POISSON :     sprintf( ModeName, "%s", "Pot4Poi"                  );   break;
         case POT_AFTER_REFINE :    sprintf( ModeName, "%s", "PotAfRef"                 );   break;
#        endif
         case COARSE_FINE_FLUX :    sprintf( ModeName, "%s", "Flux"                     );   break;
         default:                   Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GetBufMode", GetBufMode );
      }

      FILE *File = fopen( FileName, "a" );

      static bool FirstTime = true;
      if ( FirstTime )  fprintf( File, "%3s %15s %4s %4s %10s %10s %10s %8s %8s %10s %10s\n",
                                 "Lv", "Mode", "NVar", "NBuf", "Prep(s)", "Close(s)", "MPI(s)", "Send(MB)", "Recv(MB)",
                                 "Send(MB/s)", "Recv(MB/s)" );
      FirstTime = false;

      const double SendMB = NSend_Total*sizeof(real)*1.0e-6;
      const double RecvMB = NRecv_Total*sizeof(real)*1.0e-6;

      fprintf( File, "%3d %15s %4d %4d %10.5f %10.5f %10.5f %8.3f %8.3f %10.3f %10.3f\n",
               lv, ModeName, NVar_Tot, (GetBufMode==DATA_RESTRICT || GetBufMode==COARSE_FINE_FLUX)?-1:ParaBuf,
               Timer_MPI[0]->GetValue(), Timer_MPI[2]->GetValue(), Timer_MPI[1]->GetValue(),
               SendMB, RecvMB, SendMB/Timer_MPI[1]->GetValue(), RecvMB/Timer_MPI[1]->GetValue() );

      fclose( File );

      for (int t=0; t<3; t++)    Timer_MPI[t]->Reset();

   } // if ( OPT__TIMING_MPI )
#  endif // #ifdef TIMING


// free memory
   delete [] Send_NCount;
   delete [] Recv_NCount;
   delete [] Send_NDisp;
   delete [] Recv_NDisp;
   delete [] TFluVarIdxList;

} // FUNCTION : LB_GetBufferData



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetBufferData_MemAllocate_Send
// Description :  Allocate the MPI send buffer used by LG_GetBufferData (and Par_LB_SendParticleData)
//
// Note        :  1. Call LB_GetBufferData_MemFree to free memory
//                2. This function is used by some particle routines as well
//                3. We reallocate send/recv buffers only when the current buffer size is not large enough
//                   --> It greatly improves MPI performance
//
// Parameter   :  NSend : Number of elements (with the type "real) to be sent
//
// Return      :  Pointer to the MPI send buffer
//-------------------------------------------------------------------------------------------------------
real *LB_GetBufferData_MemAllocate_Send( const int NSend )
{

   if ( NSend > SendBufSize )
   {
      if ( MPI_SendBuf_Shared != NULL )   delete [] MPI_SendBuf_Shared;

//    allocate BufSizeFactor more memory to sustain longer
      SendBufSize        = int(NSend*BufSizeFactor);
      MPI_SendBuf_Shared = new real [SendBufSize];
   }

   return MPI_SendBuf_Shared;

} // FUNCTION : LB_GetBufferData_MemAllocate_Send



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetBufferData_MemAllocate_Recv
// Description :  Allocate the MPI recv buffer used by LG_GetBufferData (and Par_LB_SendParticleData)
//
// Note        :  1. Call LB_GetBufferData_MemFree to free memory
//                2. This function is used by some particle routines as well
//                3. We reallocate send/recv buffers only when the current buffer size is not large enough
//                   --> It greatly improves MPI performance
//
// Parameter   :  NRecv : Number of elements (with the type "real) to be sent
//
// Return      :  Pointer to the MPI recv buffer
//-------------------------------------------------------------------------------------------------------
real *LB_GetBufferData_MemAllocate_Recv( const int NRecv )
{

   if ( NRecv > RecvBufSize )
   {
      if ( MPI_RecvBuf_Shared != NULL )   delete [] MPI_RecvBuf_Shared;

//    allocate BufSizeFactor more memory to sustain longer
      RecvBufSize        = int(NRecv*BufSizeFactor);
      MPI_RecvBuf_Shared = new real [RecvBufSize];
   }

   return MPI_RecvBuf_Shared;

} // FUNCTION : LB_GetBufferData_MemAllocate_Recv



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GetBufferData_MemFree
// Description :  Free the MPI send and recv buffers
//
// Note        :  This function is invoked by "End_MemFree"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void LB_GetBufferData_MemFree()
{

   if ( MPI_SendBuf_Shared != NULL )
   {
      delete [] MPI_SendBuf_Shared;
      MPI_SendBuf_Shared = NULL;
   }

   if ( MPI_RecvBuf_Shared != NULL )
   {
      delete [] MPI_RecvBuf_Shared;
      MPI_RecvBuf_Shared = NULL;
   }

} // FUNCTION : LB_GetBufferData_MemFree



#endif // #ifdef LOAD_BALANCE
