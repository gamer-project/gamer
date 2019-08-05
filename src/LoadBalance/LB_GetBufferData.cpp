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
// Description :  Fill up the data of buffer patches (and fluxes and electric field of real patches as well)
//
// Note        :  1. Alternative function to Buf_GetBufferData() for LOAD_BALANCE
//                2. Invoked by Buf_GetBufferData(), EvolveLevel(), Flu_CorrAfterAllSyn(), LB_Init_ByFunction()
//                   Init_ByFile(), Init_ByRestart_v1/2(), and End_MemFree()
//                3. The modes "POT_FOR_POISSON" and "POT_AFTER_REFINE" will exchange the potential data only.
//                   The mode "COARSE_FINE_ELECTRIC" will exchange all electric field components.
//                   For others modes, the variables to be exchanged depend on the input parameters "TVarCC" and "TVarFC".
//
// Parameter   :  lv         : Target refinement level to exchage data
//                FluSg      : Sandglass of the requested fluid data
//                             --> Useless in POT_FOR_POISSON, POT_AFTER_REFINEPOT, COARSE_FINE_FLUX, COARSE_FINE_ELECTRIC
//                MagSg      : Sandglass of the requested B field data
//                             --> Useless in POT_FOR_POISSON, POT_AFTER_REFINEPOT, COARSE_FINE_FLUX, COARSE_FINE_ELECTRIC
//                PotSg      : Sandglass of the requested potential data
//                             --> Useless in COARSE_FINE_FLUX and COARSE_FINE_ELECTRIC
//                GetBufMode : Target mode. Each mode has its own MPI lists, by which the amount of data
//                             to be transferred can be minimized.
//                             --> DATA_GENERAL         : data for general-purpose (sibling and coarse-grid data)
//                                 DATA_AFTER_REFINE    : subset of DATA_GENERAL after refine
//                                 DATA_AFTER_FIXUP     : subset of DATA_GENERAL after fix-up
//                                 DATA_RESTRICT        : restricted data of the father patches with sons not home
//                                 POT_FOR_POISSON      : potential for the Poisson solver
//                                 POT_AFTER_REFINE     : potential after refine for the Poisson solver
//                                 COARSE_FINE_FLUX     : fluxes across the coarse-fine boundaries (HYDRO ONLY)
//                                 COARSE_FINE_ELECTRIC : electric field across the coarse-fine boundaries (MHD ONLY)
//                TVarCC     : Target cell-centered variables to exchange
//                             --> Supported variables in different models:
//                                 HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY [, _POTE]
//                                 ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                             --> _FLUID, _PASSIVE, and _TOTAL apply to all models
//                             --> In addition, the flux variables (e.g., _FLUX_DENS) are also supported
//                             Restrictions :
//                             --> a. DATA_XXX works with all components in (_TOTAL | _POTE)
//                                 b. COARSE_FINE_FLUX works with all components in (_FLUX_TOTAL)
//                                 c. _POTE has no effect on the flux fix-up in DATA_AFTER_FIXUP
//                                 d. POT_FOR_POISSON and POT_AFTER_REFINE only work with _POTE
//                TVarFC     : Target face-centered variables to exchange
//                             --> Supported variables in different models:
//                                 HYDRO+MHD : _MAGX, _MAGY, _MAGZ, _MAG
//                                 ELBDM     : none
//                ParaBuf    : Number of ghost zones to exchange
//                             --> Useless in DATA_RESTRICT, COARSE_FINE_FLUX, and COARSE_FINE_ELECTRIC
//-------------------------------------------------------------------------------------------------------
void LB_GetBufferData( const int lv, const int FluSg, const int MagSg, const int PotSg, const GetBufMode_t GetBufMode,
                       const int TVarCC, const int TVarFC, const int ParaBuf )
{

   bool ExchangeFlu = ( GetBufMode == COARSE_FINE_FLUX ) ?
                      TVarCC & _FLUX_TOTAL : TVarCC & _TOTAL;  // whether or not to exchage the fluid data
#  ifdef GRAVITY
   bool ExchangePot = TVarCC & _POTE;                          // whether or not to exchange the potential data
#  endif
#  ifdef MHD
   bool ExchangeMag = TVarFC & _MAG;                           // whether or not to exchange the B field data
#  endif

// TFluVarIdxList: target fluid variable indices ( = [0 ... NCOMP_TOTAL-1/NFLUX_TOTAL-1] )
   const int NFluid_Max = ( GetBufMode == COARSE_FINE_FLUX ) ? NFLUX_TOTAL : NCOMP_TOTAL;
   int NVarCC_Flu, NVarCC_Tot, *TFluVarIdxList;
   TFluVarIdxList = new int [NFluid_Max];
   NVarCC_Flu = 0;

   for (int v=0; v<NFluid_Max; v++)
      if ( TVarCC & (1<<v) )  TFluVarIdxList[ NVarCC_Flu ++ ] = v;

   NVarCC_Tot = NVarCC_Flu;
#  ifdef GRAVITY
   if ( ExchangePot )   NVarCC_Tot ++;
#  endif

#  ifdef MHD
// TMagVarIdxList: target B field indices ( = [0 ... NCOMP_MAG-1] )
   int NVarFC_Mag, *TMagVarIdxList;
   TMagVarIdxList = new int [NCOMP_MAG];
   NVarFC_Mag = 0;

   for (int v=0; v<NCOMP_MAG; v++)
      if ( TVarFC & (1<<v) )  TMagVarIdxList[ NVarFC_Mag ++ ] = v;
#  else
   const int NVarFC_Mag = 0;
#  endif

// overwrite parameters in the special modes POT_FOR_POISSON and POT_AFTER_REFINE
#  ifdef GRAVITY
   if ( GetBufMode == POT_FOR_POISSON  ||  GetBufMode == POT_AFTER_REFINE )
   {
      ExchangeFlu = false;
      ExchangePot = true;
#     ifdef MHD
      ExchangeMag = false;
#     endif

      NVarCC_Flu  = 0;
      NVarCC_Tot  = 1;
#     ifdef MHD
      NVarFC_Mag  = 0;
#     endif
   }
#  endif


// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if (  ( GetBufMode == DATA_GENERAL || GetBufMode == DATA_AFTER_FIXUP || GetBufMode == DATA_AFTER_REFINE ||
           GetBufMode == DATA_RESTRICT )  &&  NVarCC_Tot + NVarFC_Mag == 0  )
      Aux_Error( ERROR_INFO, "no target variable is found !!\n" );

   if ( GetBufMode == COARSE_FINE_FLUX  &&  NVarCC_Flu == 0 )
      Aux_Error( ERROR_INFO, "no target flux is found !!\n" );

   if ( ExchangeFlu  &&  ( FluSg != 0 && FluSg != 1 )  &&  GetBufMode != COARSE_FINE_FLUX )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FluSg", FluSg );

#  ifdef MHD
   if ( ExchangeMag  &&  ( MagSg != 0 && MagSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "MagSg", MagSg );
#  endif

#  ifdef GRAVITY
   if ( ExchangePot  &&  ( PotSg != 0 && PotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PotSg", PotSg );

   if (  ( GetBufMode == POT_FOR_POISSON || GetBufMode == POT_AFTER_REFINE )  &&  !( TVarCC & _POTE )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing _POTE !!\n" );

   if (  ( GetBufMode == POT_FOR_POISSON || GetBufMode == POT_AFTER_REFINE )  &&  ( TVarCC & ~_POTE )  )
      Aux_Error( ERROR_INFO, "modes \"%s\" only accept \"%s\" as the target variable !!\n",
                 "POT_FOR_POISSON and POT_AFTER_REFINE", "_POTE" );
#  endif

   if (  ( GetBufMode == DATA_GENERAL || GetBufMode == DATA_AFTER_FIXUP || GetBufMode == DATA_AFTER_REFINE
#          ifdef GRAVITY
           || GetBufMode == POT_FOR_POISSON || GetBufMode == POT_AFTER_REFINE
#          endif
          )  &&  ( ParaBuf < 0 || ParaBuf > PATCH_SIZE )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d --> accepted range = [0 ... PATCH_SIZE] !!\n",
                 "ParaBuf", ParaBuf );

   if ( GetBufMode == COARSE_FINE_FLUX  &&  !amr->WithFlux )
      Aux_Error( ERROR_INFO, "COARSE_FINE_FLUX failed since flux arrays are not allocated !!\n" );

#  ifdef MHD
   if ( GetBufMode == COARSE_FINE_ELECTRIC  &&  !amr->WithElectric )
      Aux_Error( ERROR_INFO, "WARNING : COARSE_FINE_ELECTRIC failed since electric field arrays are not allocated !!\n" );
#  endif


// _X : for exchanging data after the flux fix-up, which has ParaBuf=1
   const int DataUnit_Flux = SQR( PS1 )*NVarCC_Flu;
#  ifdef MHD
   const int ParaBufP1     = ParaBuf + 1;
#  endif

   int   NSend_Total, NRecv_Total;
   int   DataUnit_Buf[27], LoopStart[27][3], LoopEnd[27][3];
   int   LoopStart_X[6][3], LoopEnd_X[6][3];
   int  *Send_NList=NULL, *Recv_NList=NULL, *Send_NResList=NULL, *Recv_NResList=NULL;
   int **Send_IDList=NULL, **Recv_IDList=NULL, **Send_IDList_IdxTable=NULL, **Recv_IDList_IdxTable=NULL;
   int **Send_SibList=NULL, **Recv_SibList=NULL;

   int *Send_NCount = new int [MPI_NRank];
   int *Recv_NCount = new int [MPI_NRank];
   int *Send_NDisp  = new int [MPI_NRank];
   int *Recv_NDisp  = new int [MPI_NRank];


// 1. set up the number of elements to be sent and received in each cell and the send/recv lists
// ============================================================================================================
   switch ( GetBufMode )
   {
      case DATA_GENERAL :           // general-purpose
         Send_NList           = amr->LB->SendH_NList          [lv];
         Send_IDList          = amr->LB->SendH_IDList         [lv];
         Send_SibList         = amr->LB->SendH_SibList        [lv];
         Recv_NList           = amr->LB->RecvH_NList          [lv];
         Recv_IDList          = amr->LB->RecvH_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvH_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvH_SibList        [lv];
         break;

      case DATA_AFTER_REFINE :      // data after refine
         Send_NList           = amr->LB->SendH_NList          [lv];
         Send_IDList          = amr->LB->SendH_IDList         [lv];
         Send_SibList         = amr->LB->SendH_SibDiffList    [lv];
         Recv_NList           = amr->LB->RecvH_NList          [lv];
         Recv_IDList          = amr->LB->RecvH_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvH_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvH_SibDiffList    [lv];
         break;

      case DATA_AFTER_FIXUP :       // data after fix-up
         Send_NList           = amr->LB->SendX_NList          [lv];
         Send_NResList        = amr->LB->SendX_NResList       [lv];
         Send_IDList          = amr->LB->SendX_IDList         [lv];
         Send_SibList         = amr->LB->SendX_SibList        [lv];
         Recv_NList           = amr->LB->RecvX_NList          [lv];
         Recv_NResList        = amr->LB->RecvX_NResList       [lv];
         Recv_IDList          = amr->LB->RecvX_IDList         [lv];
         Recv_SibList         = amr->LB->RecvX_SibList        [lv];
         break;

      case DATA_RESTRICT :          // restricted data
         Send_NList           = amr->LB->SendR_NList          [lv];
         Send_IDList          = amr->LB->SendR_IDList         [lv];
         Send_IDList_IdxTable = amr->LB->SendR_IDList_IdxTable[lv];
         Recv_NList           = amr->LB->RecvR_NList          [lv];
         Recv_IDList          = amr->LB->RecvR_IDList         [lv];
         break;

#     ifdef GRAVITY
      case POT_FOR_POISSON :        // potential for the Poisson solver
         Send_NList           = amr->LB->SendG_NList          [lv];
         Send_IDList          = amr->LB->SendG_IDList         [lv];
         Send_SibList         = amr->LB->SendG_SibList        [lv];
         Recv_NList           = amr->LB->RecvG_NList          [lv];
         Recv_IDList          = amr->LB->RecvG_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvG_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvG_SibList        [lv];
         break;

      case POT_AFTER_REFINE :       // potential after refine for the Poisson solver
         Send_NList           = amr->LB->SendG_NList          [lv];
         Send_IDList          = amr->LB->SendG_IDList         [lv];
         Send_SibList         = amr->LB->SendG_SibDiffList    [lv];
         Recv_NList           = amr->LB->RecvG_NList          [lv];
         Recv_IDList          = amr->LB->RecvG_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvG_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvG_SibDiffList    [lv];
         break;
#     endif // #ifdef GRAVITY

      case COARSE_FINE_FLUX :       // hydro fluxes
         Send_NList           = amr->LB->SendF_NList          [lv];
         Send_IDList          = amr->LB->SendF_IDList         [lv];
         Send_SibList         = amr->LB->SendF_SibList        [lv];
         Recv_NList           = amr->LB->RecvF_NList          [lv];
         Recv_IDList          = amr->LB->RecvF_IDList         [lv];
         Recv_IDList_IdxTable = amr->LB->RecvF_IDList_IdxTable[lv];
         Recv_SibList         = amr->LB->RecvF_SibList        [lv];
         break;

#     ifdef MHD
      case COARSE_FINE_ELECTRIC :   // MHD electric field
         Aux_Error( ERROR_INFO, "NOT SUPPORTED YET !!\n" );
         break;
#     endif

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GetBufMode", GetBufMode );
   } // switch ( GetBufMode )



// 2. set up the loop range, MPI count and displacement arrays, and allocate data for send and recv buffers
// ============================================================================================================
   switch ( GetBufMode )
   {
      case DATA_GENERAL : case DATA_AFTER_REFINE :
#     ifdef GRAVITY
      case POT_FOR_POISSON : case POT_AFTER_REFINE :
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
         for (int s= 0; s< 6; s++)  DataUnit_Buf[ s] = NVarCC_Tot * PS1     * PS1     * ParaBuf;   // sibling 0 ~ 5
         for (int s= 6; s<18; s++)  DataUnit_Buf[ s] = NVarCC_Tot * PS1     * ParaBuf * ParaBuf;   // sibling 6 ~ 17
         for (int s=18; s<26; s++)  DataUnit_Buf[ s] = NVarCC_Tot * ParaBuf * ParaBuf * ParaBuf;   // sibling 18 ~ 25
                                    DataUnit_Buf[26] = NVarCC_Tot * PS1     * PS1     * PS1;       // entire patch

#        ifdef MHD
         for (int v=0; v<NVarFC_Mag; v++)
         {
            const int TMagVarIdx = TMagVarIdxList[v];

            switch ( TMagVarIdx )
            {
               case MAGX :
                  for (int s= 0; s< 2; s++)  DataUnit_Buf[s] += PS1   * PS1     * ParaBufP1;
                  for (int s= 2; s< 6; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 6; s<10; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  for (int s=10; s<14; s++)  DataUnit_Buf[s] += PS1P1 * ParaBuf * ParaBuf;
                  for (int s=14; s<18; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  break;

               case MAGY :
                  for (int s= 0; s< 2; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 2; s< 4; s++)  DataUnit_Buf[s] += PS1   * PS1     * ParaBufP1;
                  for (int s= 4; s< 6; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 6; s<14; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  for (int s=14; s<18; s++)  DataUnit_Buf[s] += PS1P1 * ParaBuf * ParaBuf;
                  break;

               case MAGZ :
                  for (int s= 0; s< 4; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 4; s< 6; s++)  DataUnit_Buf[s] += PS1   * PS1     * ParaBufP1;
                  for (int s= 6; s<10; s++)  DataUnit_Buf[s] += PS1P1 * ParaBuf * ParaBuf;
                  for (int s=10; s<18; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  break;

               default:
                  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TMagVarIdx", TMagVarIdx );
                  break;
            } // switch ( TMagVarIdx )
         } // for (int v=0; v<NVarFC_Mag; v++)

         for (int s=18; s<26; s++)  DataUnit_Buf[ s] += NVarFC_Mag*ParaBufP1*SQR( ParaBuf );
                                    DataUnit_Buf[26] += NVarFC_Mag*PS1P1*SQR( PS1 );
#        endif // #ifdef MHD

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
         break; // cases DATA_GENERAL, DATA_AFTER_REFINE, POT_FOR_POISSON, POT_AFTER_REFINE


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
         for (int s= 0; s< 6; s++)  DataUnit_Buf[ s] = NVarCC_Tot * PS1     * PS1     * ParaBuf;   // sibling 0 ~ 5
         for (int s= 6; s<18; s++)  DataUnit_Buf[ s] = NVarCC_Tot * PS1     * ParaBuf * ParaBuf;   // sibling 6 ~ 17
         for (int s=18; s<26; s++)  DataUnit_Buf[ s] = NVarCC_Tot * ParaBuf * ParaBuf * ParaBuf;   // sibling 18 ~ 25
                                    DataUnit_Buf[26] = NVarCC_Tot * PS1     * PS1     * PS1;       // entire patch

#        ifdef MHD
         for (int v=0; v<NVarFC_Mag; v++)
         {
            const int TMagVarIdx = TMagVarIdxList[v];

            switch ( TMagVarIdx )
            {
               case MAGX :
                  for (int s= 0; s< 2; s++)  DataUnit_Buf[s] += PS1   * PS1     * ParaBufP1;
                  for (int s= 2; s< 6; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 6; s<10; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  for (int s=10; s<14; s++)  DataUnit_Buf[s] += PS1P1 * ParaBuf * ParaBuf;
                  for (int s=14; s<18; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  break;

               case MAGY :
                  for (int s= 0; s< 2; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 2; s< 4; s++)  DataUnit_Buf[s] += PS1   * PS1     * ParaBufP1;
                  for (int s= 4; s< 6; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 6; s<14; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  for (int s=14; s<18; s++)  DataUnit_Buf[s] += PS1P1 * ParaBuf * ParaBuf;
                  break;

               case MAGZ :
                  for (int s= 0; s< 4; s++)  DataUnit_Buf[s] += PS1P1 * PS1     * ParaBuf;
                  for (int s= 4; s< 6; s++)  DataUnit_Buf[s] += PS1   * PS1     * ParaBufP1;
                  for (int s= 6; s<10; s++)  DataUnit_Buf[s] += PS1P1 * ParaBuf * ParaBuf;
                  for (int s=10; s<18; s++)  DataUnit_Buf[s] += PS1   * ParaBuf * ParaBufP1;
                  break;

               default:
                  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TMagVarIdx", TMagVarIdx );
                  break;
            } // switch ( TMagVarIdx )
         } // for (int v=0; v<NVarFC_Mag; v++)

         for (int s=18; s<26; s++)  DataUnit_Buf[ s] += NVarFC_Mag*ParaBufP1*SQR( ParaBuf );
                                    DataUnit_Buf[26] += NVarFC_Mag*PS1P1*SQR( PS1 );
#        endif // #ifdef MHD

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
         break; // case DATA_AFTER_FIXUP


      case DATA_RESTRICT :
//    ----------------------------------------------
         for (int r=0; r<MPI_NRank; r++)
         {
            Send_NCount[r]  = Send_NList[r]*CUBE( PS1 )*NVarCC_Tot;
            Recv_NCount[r]  = Recv_NList[r]*CUBE( PS1 )*NVarCC_Tot;
#           ifdef MHD
            Send_NCount[r] += Send_NList[r]*SQR( PS1 )*PS1P1*NVarFC_Mag;
            Recv_NCount[r] += Recv_NList[r]*SQR( PS1 )*PS1P1*NVarFC_Mag;
#           endif
         }
         break; // case DATA_RESTRICT


      case COARSE_FINE_FLUX :
//    ----------------------------------------------
         for (int r=0; r<MPI_NRank; r++)
         {
            Send_NCount[r] = Send_NList[r]*DataUnit_Flux;
            Recv_NCount[r] = Recv_NList[r]*DataUnit_Flux;
         }
         break; // case COARSE_FINE_FLUX
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
      case DATA_GENERAL : case DATA_AFTER_REFINE :
#     ifdef GRAVITY
      case POT_FOR_POISSON : case POT_AFTER_REFINE :
#     endif
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *SendPtr = SendBuf + Send_NDisp[r];
            int   Counter = 0;

            for (int t=0; t<Send_NList[r]; t++)
            {
               const int SPID = Send_IDList [r][t];   // both SPID and SSib are sorted
               const int SSib = Send_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu )
               {
                  if ( amr->patch[FluSg][lv][SPID]->fluid == NULL )
                     Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
                                GetBufMode, FluSg, lv, SPID );

                  if ( SSib == 0  &&  GetBufMode != DATA_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                                GetBufMode, t, r, SPID );
               }

#              ifdef GRAVITY
               if ( ExchangePot )
               {
                  if ( amr->patch[PotSg][lv][SPID]->pot == NULL )
                     Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->pot has not been allocated !!\n",
                                GetBufMode, PotSg, lv, SPID );

                  if ( SSib == 0  &&  GetBufMode != POT_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                                GetBufMode, t, r, SPID );
               }
#              endif // #ifdef GRAVITY

#              ifdef MHD
               if ( ExchangeMag )
               {
                  if ( amr->patch[MagSg][lv][SPID]->magnetic == NULL )
                     Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->magnetic has not been allocated !!\n",
                                GetBufMode, MagSg, lv, SPID );

                  if ( SSib == 0  &&  GetBufMode != DATA_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Send mode %d, t %d, TRank %d, SPID %d, SSib == 0 !!\n",
                                GetBufMode, t, r, SPID );
               }
#              endif // #ifdef MHD
#              endif // #ifdef GAMER_DEBUG

               if ( SSib )
               for (int s=0; s<27; s++)
               {
                  if ( SSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVarCC_Flu; v++)
                     {
                        const int TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][k][j][i];
                     }

//                   potential data
#                    ifdef GRAVITY
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[PotSg][lv][SPID]->pot[k][j][i];
                     }
#                    endif

//                   magnetic field data
#                    ifdef MHD
                     if ( ExchangeMag )
                     for (int v=0; v<NVarFC_Mag; v++)
                     {
                        const int TMagVarIdx = TMagVarIdxList[v];

                        switch ( TMagVarIdx )
                        {
                           case MAGX :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i<=LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BX( i, j, k, PS1, PS1 );
                                 SendPtr[ Counter ++ ] = amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][idxB];
                              }
                              break;

                           case MAGY :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j<=LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BY( i, j, k, PS1, PS1 );
                                 SendPtr[ Counter ++ ] = amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][idxB];
                              }
                              break;

                            case MAGZ :
                              for (int k=LoopStart[s][2]; k<=LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BZ( i, j, k, PS1, PS1 );
                                 SendPtr[ Counter ++ ] = amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][idxB];
                              }
                              break;

                            default:
                              Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TMagVarIdx", TMagVarIdx );
                              break;
                        } // switch ( TMagVarIdx )
                     } // for (int v=0; v<NVarFC_Mag; v++)
#                    endif // #ifdef MHD

                  } // if ( SSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // cases DATA_GENERAL, DATA_AFTER_REFINE, POT_FOR_POISSON, POT_AFTER_REFINE


      case DATA_AFTER_FIXUP :
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *SendPtr = SendBuf + Send_NDisp[r];
            int   Counter = 0;

//          for restriction fix-up
            for (int t=0; t<Send_NResList[r]; t++)
            {
               const int SPID = Send_IDList [r][t];
               const int SSib = Send_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][SPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, SPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][SPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, SPID );
#              endif

#              ifdef MHD
               if ( ExchangeMag  &&  amr->patch[MagSg][lv][SPID]->magnetic == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->magnetic has not been allocated !!\n",
                             GetBufMode, MagSg, lv, SPID );
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
                     for (int v=0; v<NVarCC_Flu; v++)
                     {
                        const int TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][k][j][i];
                     }

//                   potential data
#                    ifdef GRAVITY
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[PotSg][lv][SPID]->pot[k][j][i];
                     }
#                    endif

//                   magnetic field data
#                    ifdef MHD
                     if ( ExchangeMag )
                     for (int v=0; v<NVarFC_Mag; v++)
                     {
                        const int TMagVarIdx = TMagVarIdxList[v];

                        switch ( TMagVarIdx )
                        {
                           case MAGX :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i<=LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BX( i, j, k, PS1, PS1 );
                                 SendPtr[ Counter ++ ] = amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][idxB];
                              }
                              break;

                           case MAGY :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j<=LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BY( i, j, k, PS1, PS1 );
                                 SendPtr[ Counter ++ ] = amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][idxB];
                              }
                              break;

                            case MAGZ :
                              for (int k=LoopStart[s][2]; k<=LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BZ( i, j, k, PS1, PS1 );
                                 SendPtr[ Counter ++ ] = amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][idxB];
                              }
                              break;

                            default:
                              Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TMagVarIdx", TMagVarIdx );
                              break;
                        } // switch ( TMagVarIdx )
                     } // for (int v=0; v<NVarFC_Mag; v++)
#                    endif // #ifdef MHD

                  } // if ( SSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Send_NResList[r]; t++)


//          for flux fix-up
            for (int t=Send_NResList[r]; t<Send_NList[r]; t++)
            {
               const int SPID = Send_IDList [r][t];   // both SPID and SSib are sorted
               const int SSib = Send_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][SPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
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
                     for (int v=0; v<NVarCC_Flu; v++)
                     {
                        const int TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart_X[s][2]; k<LoopEnd_X[s][2]; k++)
                        for (int j=LoopStart_X[s][1]; j<LoopEnd_X[s][1]; j++)
                        for (int i=LoopStart_X[s][0]; i<LoopEnd_X[s][0]; i++)
                           SendPtr[ Counter ++ ] = amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][k][j][i];
                     }
                  } // if ( SSib & (1<<s) )
               } // for (int s=0; s<6; s++)
            } //for (int t=Send_NResList[r]; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_AFTER_FIXUP


      case DATA_RESTRICT :
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *SendPtr = SendBuf + Send_NDisp[r];

            for (int t=0; t<Send_NList[r]; t++)
            {
               const int SPID = Send_IDList[r][ Send_IDList_IdxTable[r][t] ];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][SPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, SPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][SPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, SPID );
#              endif
#              ifdef MHD
               if ( ExchangeMag  &&  amr->patch[MagSg][lv][SPID]->magnetic == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[%d][%d][%d]->magnetic has not been allocated !!\n",
                             GetBufMode, MagSg, lv, SPID );
#              endif
#              endif // #ifdef GAMER_DEBUG

//             fluid data
               if ( ExchangeFlu )
               for (int v=0; v<NVarCC_Flu; v++)
               {
                  const int TFluVarIdx = TFluVarIdxList[v];

                  memcpy( SendPtr, &amr->patch[FluSg][lv][SPID]->fluid[TFluVarIdx][0][0][0],
                          PS1*PS1*PS1*sizeof(real) );

                  SendPtr += CUBE( PS1 );
               }

//             potential data
#              ifdef GRAVITY
               if ( ExchangePot )
               {
                  memcpy( SendPtr, &amr->patch[PotSg][lv][SPID]->pot[0][0][0],
                          PS1*PS1*PS1*sizeof(real) );

                  SendPtr += CUBE( PS1 );
               }
#              endif

//             magnetic field data
#              ifdef MHD
               if ( ExchangeMag )
               for (int v=0; v<NVarFC_Mag; v++)
               {
                  const int TMagVarIdx = TMagVarIdxList[v];

                  memcpy( SendPtr, &amr->patch[MagSg][lv][SPID]->magnetic[TMagVarIdx][0],
                          SQR(PS1)*PS1P1*sizeof(real) );

                  SendPtr += SQR( PS1 )*PS1P1;
               }
#              endif
            } // for (int t=0; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_RESTRICT


      case COARSE_FINE_FLUX :
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *SendPtr = SendBuf + Send_NDisp[r];
            int   Counter = 0;

            for (int t=0; t<Send_NList[r]; t++)
            {
               const int SPID = Send_IDList [r][t];
               const int SSib = Send_SibList[r][t];
               const real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][SPID]->flux[SSib];

#              ifdef GAMER_DEBUG
               if ( FluxPtr == NULL )
                  Aux_Error( ERROR_INFO, "Send mode %d, patch[0][%d][%d]->flux[%d] has not been allocated !!\n",
                             GetBufMode, lv, SPID, SSib );
#              endif

               for (int v=0; v<NVarCC_Flu; v++)
               {
                  const int TFluVarIdx = TFluVarIdxList[v];

                  memcpy( SendPtr, FluxPtr[TFluVarIdx], PS1*PS1*sizeof(real) );

                  SendPtr += SQR( PS1 );
               }
            } // for (int t=0; t<Send_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case COARSE_FINE_FLUX


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
      case DATA_GENERAL : case DATA_AFTER_REFINE :
#     ifdef GRAVITY
      case POT_FOR_POISSON : case POT_AFTER_REFINE :
#     endif
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *RecvPtr = RecvBuf + Recv_NDisp[r];
            int   Counter = 0;

            for (int t=0; t<Recv_NList[r]; t++)
            {
               const int RPID = Recv_IDList [r][ Recv_IDList_IdxTable[r][t] ];   // Recv_IDList is unsorted
               const int RSib = Recv_SibList[r][t];                              // Recv_SibList is sorted

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu )
               {
                  if ( amr->patch[FluSg][lv][RPID]->fluid == NULL )
                     Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
                                GetBufMode, FluSg, lv, RPID );

                  if ( RSib == 0  &&  GetBufMode != DATA_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                                GetBufMode, t, r, RPID );
               }

#              ifdef GRAVITY
               if ( ExchangePot )
               {
                  if ( amr->patch[PotSg][lv][RPID]->pot == NULL )
                     Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->pot has not been allocated !!\n",
                                GetBufMode, PotSg, lv, RPID );

                  if ( RSib == 0  &&  GetBufMode != POT_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                                GetBufMode, t, r, RPID );
               }
#              endif // #ifdef GRAVITY

#              ifdef MHD
               if ( ExchangeMag )
               {
                  if ( amr->patch[MagSg][lv][RPID]->magnetic == NULL )
                     Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->magnetic has not been allocated !!\n",
                                GetBufMode, MagSg, lv, RPID );

                  if ( RSib == 0  &&  GetBufMode != DATA_AFTER_REFINE )
                     Aux_Error( ERROR_INFO, "Recv mode %d, t %d, TRank %d, RPID %d, RSib == 0 !!\n",
                                GetBufMode, t, r, RPID );
               }
#              endif // #ifdef MHD
#              endif // #ifdef GAMER_DEBUG

               if ( RSib )
               for (int s=0; s<27; s++)
               {
                  if ( RSib & (1<<s) )
                  {
//                   fluid data
                     if ( ExchangeFlu )
                     for (int v=0; v<NVarCC_Flu; v++)
                     {
                        const int TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][k][j][i] = RecvPtr[ Counter ++ ];
                     }

//                   potential data
#                    ifdef GRAVITY
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[PotSg][lv][RPID]->pot[k][j][i] = RecvPtr[ Counter ++ ];
                     }
#                    endif

//                   magnetic field data
#                    ifdef MHD
                     if ( ExchangeMag )
                     for (int v=0; v<NVarFC_Mag; v++)
                     {
                        const int TMagVarIdx = TMagVarIdxList[v];

                        switch ( TMagVarIdx )
                        {
                           case MAGX :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i<=LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BX( i, j, k, PS1, PS1 );
                                 amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][idxB] = RecvPtr[ Counter ++ ];
                              }
                              break;

                           case MAGY :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j<=LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BY( i, j, k, PS1, PS1 );
                                 amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][idxB] = RecvPtr[ Counter ++ ];
                              }
                              break;

                           case MAGZ :
                              for (int k=LoopStart[s][2]; k<=LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BZ( i, j, k, PS1, PS1 );
                                 amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][idxB] = RecvPtr[ Counter ++ ];
                              }
                              break;

                           default:
                              Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TMagVarIdx", TMagVarIdx );
                              break;
                         } // switch ( TMagVarIdx )
                     } //for (int v=0; v<NVarFC_Mag; v++)
#                    endif // #ifdef MHD

                  } // if ( RSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // cases DATA_GENERAL, DATA_AFTER_REFINE, POT_FOR_POISSON, POT_AFTER_REFINE


      case DATA_AFTER_FIXUP :
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *RecvPtr = RecvBuf + Recv_NDisp[r];
            int   Counter = 0;

//          for restriction fix-up
            for (int t=0; t<Recv_NResList[r]; t++)
            {
               const int RPID = Recv_IDList [r][t];
               const int RSib = Recv_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][RPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, RPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][RPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, RPID );
#              endif

#              ifdef MHD
               if ( ExchangeMag  &&  amr->patch[MagSg][lv][RPID]->magnetic == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->magnetic has not been allocated !!\n",
                             GetBufMode, MagSg, lv, RPID );
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
                     for (int v=0; v<NVarCC_Flu; v++)
                     {
                        const int TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][k][j][i] = RecvPtr[ Counter ++ ];
                     }

//                   potential data
#                    ifdef GRAVITY
                     if ( ExchangePot )
                     {
                        for (int k=LoopStart[s][2]; k<LoopEnd[s][2]; k++)
                        for (int j=LoopStart[s][1]; j<LoopEnd[s][1]; j++)
                        for (int i=LoopStart[s][0]; i<LoopEnd[s][0]; i++)
                           amr->patch[PotSg][lv][RPID]->pot[k][j][i] = RecvPtr[ Counter ++ ];
                     }
#                    endif

//                   magnetic field data
#                    ifdef MHD
                     if ( ExchangeMag )
                     for (int v=0; v<NVarFC_Mag; v++)
                     {
                        const int TMagVarIdx = TMagVarIdxList[v];

                        switch ( TMagVarIdx )
                        {
                           case MAGX :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i<=LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BX( i, j, k, PS1, PS1 );
                                 amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][idxB] = RecvPtr[ Counter ++ ];
                              }
                              break;

                           case MAGY :
                              for (int k=LoopStart[s][2]; k< LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j<=LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BY( i, j, k, PS1, PS1 );
                                 amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][idxB] = RecvPtr[ Counter ++ ];
                              }
                              break;

                           case MAGZ :
                              for (int k=LoopStart[s][2]; k<=LoopEnd[s][2]; k++)
                              for (int j=LoopStart[s][1]; j< LoopEnd[s][1]; j++)
                              for (int i=LoopStart[s][0]; i< LoopEnd[s][0]; i++)
                              {
                                 const int idxB = IDX321_BZ( i, j, k, PS1, PS1 );
                                 amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][idxB] = RecvPtr[ Counter ++ ];
                              }
                              break;

                           default:
                              Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TMagVarIdx", TMagVarIdx );
                              break;
                         } // switch ( TMagVarIdx )
                     } //for (int v=0; v<NVarFC_Mag; v++)
#                    endif // #ifdef MHD

                  } // if ( RSib & (1<<s) )
               } // for (int s=0; s<27; s++)
            } // for (int t=0; t<Recv_NResList[r]; t++)


//          for flux fix-up
            for (int t=Recv_NResList[r]; t<Recv_NList[r]; t++)
            {
               const int RPID = Recv_IDList [r][t];
               const int RSib = Recv_SibList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][RPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
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
                     for (int v=0; v<NVarCC_Flu; v++)
                     {
                        const int TFluVarIdx = TFluVarIdxList[v];

                        for (int k=LoopStart_X[s][2]; k<LoopEnd_X[s][2]; k++)
                        for (int j=LoopStart_X[s][1]; j<LoopEnd_X[s][1]; j++)
                        for (int i=LoopStart_X[s][0]; i<LoopEnd_X[s][0]; i++)
                           amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][k][j][i] = RecvPtr[ Counter ++ ];
                     }
                  } // if ( RSib & (1<<s) )
               } // for (int s=0; s<6; s++)
            } // for (int t=Recv_NResList[r]; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_AFTER_FIXUP


      case DATA_RESTRICT :
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *RecvPtr = RecvBuf + Recv_NDisp[r];

            for (int t=0; t<Recv_NList[r]; t++)
            {
               const int RPID = Recv_IDList[r][t];

#              ifdef GAMER_DEBUG
               if ( ExchangeFlu  &&  amr->patch[FluSg][lv][RPID]->fluid == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->fluid has not been allocated !!\n",
                             GetBufMode, FluSg, lv, RPID );

#              ifdef GRAVITY
               if ( ExchangePot  &&  amr->patch[PotSg][lv][RPID]->pot == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->pot has not been allocated !!\n",
                             GetBufMode, PotSg, lv, RPID );
#              endif

#              ifdef MHD
               if ( ExchangeMag  &&  amr->patch[MagSg][lv][RPID]->magnetic == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[%d][%d][%d]->magnetic has not been allocated !!\n",
                             GetBufMode, MagSg, lv, RPID );
#              endif
#              endif // #ifdef GAMER_DEBUG

//             fluid data
               if ( ExchangeFlu )
               for (int v=0; v<NVarCC_Flu; v++)
               {
                  const int TFluVarIdx = TFluVarIdxList[v];
                  memcpy( &amr->patch[FluSg][lv][RPID]->fluid[TFluVarIdx][0][0][0], RecvPtr, CUBE(PS1)*sizeof(real) );
                  RecvPtr += CUBE( PS1 );
               }

//             potential data
#              ifdef GRAVITY
               if ( ExchangePot )
               {
                  memcpy( &amr->patch[PotSg][lv][RPID]->pot[0][0][0], RecvPtr, CUBE(PS1)*sizeof(real) );
                  RecvPtr += CUBE( PS1 );
               }
#              endif

//             magnetic field data
#              ifdef MHD
               if ( ExchangeMag )
               for (int v=0; v<NVarFC_Mag; v++)
               {
                  const int TMagVarIdx = TMagVarIdxList[v];
                  memcpy( &amr->patch[MagSg][lv][RPID]->magnetic[TMagVarIdx][0], RecvPtr, SQR(PS1)*PS1P1*sizeof(real) );
                  RecvPtr += SQR( PS1 )*PS1P1;
               }
#              endif
            } // for (int t=0; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case DATA_RESTRICT


      case COARSE_FINE_FLUX :
//    ----------------------------------------------
#        pragma omp parallel for schedule( runtime )
         for (int r=0; r<MPI_NRank; r++)
         {
            real *RecvPtr = RecvBuf + Recv_NDisp[r];
            int   Counter = 0;

            for (int t=0; t<Recv_NList[r]; t++)
            {
               const int RPID = Recv_IDList [r][ Recv_IDList_IdxTable[r][t] ];
               const int RSib = Recv_SibList[r][t];
               real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][RPID]->flux[RSib];

#              ifdef GAMER_DEBUG
               if ( FluxPtr == NULL )
                  Aux_Error( ERROR_INFO, "Recv mode %d, patch[0][%d][%d]->flux[%d] has not been allocated !!\n",
                             GetBufMode, lv, RPID, RSib );
#              endif

//             add (not replace) flux array with the received flux
               for (int v=0; v<NVarCC_Flu; v++)
               {
                  const int TFluVarIdx = TFluVarIdxList[v];

                  for (int m=0; m<PS1; m++)
                  for (int n=0; n<PS1; n++)
                     FluxPtr[TFluVarIdx][m][n] += RecvPtr[ Counter ++ ];
               }
            } // for (int t=0; t<Recv_NList[r]; t++)
         } // for (int r=0; r<MPI_NRank; r++)
         break; // case COARSE_FINE_FLUX


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
         case DATA_GENERAL :        sprintf( ModeName, "%s", (NVarCC_Tot==1)?"Dens":"Flu" );   break;
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
               lv, ModeName, NVarCC_Tot, (GetBufMode==DATA_RESTRICT || GetBufMode==COARSE_FINE_FLUX)?-1:ParaBuf,
               Timer_MPI[0]->GetValue(), Timer_MPI[2]->GetValue(), Timer_MPI[1]->GetValue(),
               SendMB, RecvMB, SendMB/Timer_MPI[1]->GetValue(), RecvMB/Timer_MPI[1]->GetValue() );

      fclose( File );

      for (int t=0; t<3; t++)    Timer_MPI[t]->Reset();

   } // if ( OPT__TIMING_MPI )
#  endif // #ifdef TIMING



// 7. free memory
// ============================================================================================================
   delete [] Send_NCount;
   delete [] Recv_NCount;
   delete [] Send_NDisp;
   delete [] Recv_NDisp;
   delete [] TFluVarIdxList;
#  ifdef MHD
   delete [] TMagVarIdxList;
#  endif



#  ifdef MHD
// 8. ensure consistency of B field on the common interfaces between nearby **real** patches
//    --> do this after recording MPI bandwidth to avoid overwriting Timer_MPI[]
//    --> this operation is necessary because nearby real patches may reside on different ranks
//        --> for nearby real patches in the same rank, it is done directly in Flu_FixUp_Restrcit()
//        --> only necessary on coarse-fine interfaces
//    --> procedure: real son -> buffer father in the same rank -> real father in another rank
//                   -> buffer father in all ranks -> real sibling father
//    --> we will further invoke MHD_LB_EnsureBFieldConsistencyAfterRestrict() in EvolveLevel()
//        to ensure consistency of B field on the common interfaces between nearby **buffer** patches
// ============================================================================================================
   if ( GetBufMode == DATA_RESTRICT )
   {
      const int MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };

//    trnasfer the restricted B field from real to buffer patches
//    --> set the number of ghost zones to zero to only transfer B field on the patch boundaries
//###OPTIMIZATION: (a) only need to transfer data in between **non-leaf** patches
//                 (b) only buffer patches with **leaf** sibling real patches need to receive data
      const int ParaBufZero = 0;
      LB_GetBufferData( lv, NULL_INT, MagSg, NULL_INT, DATA_GENERAL, _NONE, TVarFC, ParaBufZero );

//    copy data from non-leaf buffer patches to leaf real patches
//    --> for simplicity, it always works on all three components regardless of TVarFC
      for (int RealPID=0; RealPID<amr->NPatchComma[lv][1]; RealPID++)
      {
         if ( amr->patch[0][lv][RealPID]->son == -1 )
         {
            for (int s=0; s<6; s++)
            {
               const int SibBufPID = amr->patch[0][lv][RealPID]->sibling[s];

               if ( SibBufPID >= amr->NPatchComma[lv][1]  &&  amr->patch[0][lv][SibBufPID]->son != -1 )
                  MHD_CopyPatchInterfaceBField( lv, SibBufPID, MirrorSib[s], MagSg );
            }
         }
      }
   } // if ( GetBufMode == DATA_RESTRICT )



// 9. update the B field on the coarse-fine interfaces of **leaf buffer** patches after applying the restrict operation
//    --> to ensure consistency of B field between leaf buffer patches and their sibling real/buffer non-leaf patches
//    --> necessary because the mode DATA_AFTER_FIXUP will only re-transfer data for **non-leaf** buffer patches
//    --> note that even when OPT__FIXUP_RESTRICT is off we still need to do data restriction in several places
//        (e.g., restart and OPT__CORR_AFTER_ALL_SYNC)
//    --> for simplicity and sustainability, we always perform the following operation even when OPT__FIXUP_RESTRICT is off)
//    --> for simplicity, it always works on all three components regardless of TVarFC
// ============================================================================================================
   if ( GetBufMode == DATA_AFTER_FIXUP )
      MHD_LB_EnsureBFieldConsistencyAfterRestrict( lv );
#  endif // #ifdef MHD

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
// Parameter   :  NSend : Number of elements (with the type "real") to be sent
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
// Parameter   :  NRecv : Number of elements (with the type "real") to be sent
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
