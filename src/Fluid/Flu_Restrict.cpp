#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Restrict
// Description :  Replace the data at level "FaLv" by the average data at level "FaLv+1"
//
// Note        :  Use the input parameter "TVar" to determine the target variables, which can be any
//                subset of (_FLUID | _POTE | _PASSIVE)
//
// Parameter   :  FaLv     : Target refinement level at which the data are going to be replaced
//                SonFluSg : Fluid sandglass at level "FaLv+1"
//                FaFluSg  : Fluid sandglass at level "FaLv"
//                SonPotSg : Potential sandglass at level "FaLv+1"
//                FaPotSg  : Potential sandglass at level "FaLv"
//                TVar     : Target variables
//                           --> Supported variables in different models:
//                               HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY,[, _POTE]
//                               MHD   :
//                               ELBDM : _DENS, _REAL, _IMAG, [, _POTE]
//                           --> _FLUID, _PASSIVE, and _TOTAL apply to all models
//-------------------------------------------------------------------------------------------------------
void Flu_Restrict( const int FaLv, const int SonFluSg, const int FaFluSg, const int SonPotSg, const int FaPotSg,
                   const int TVar )
{

   const int SonLv = FaLv + 1;

// check
   if ( FaLv < 0  ||  FaLv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );

   if ( FaLv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : applying \"%s\" to the maximum level is meaningless !! \n", __FUNCTION__ );
      return;
   }

   if (  ( TVar & (_TOTAL) )  &&  ( SonFluSg != 0 && SonFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonFluSg", SonFluSg );

   if (  ( TVar & (_TOTAL) )  &&  ( FaFluSg != 0 && FaFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaFluSg", FaFluSg );

#  ifdef GRAVITY
   if (  ( TVar & _POTE )  &&  ( SonPotSg != 0 && SonPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonPotSg", SonPotSg );

   if (  ( TVar & _POTE )  &&  ( FaPotSg != 0 && FaPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPotSg", FaPotSg );

   if (  !( TVar & (_TOTAL|_POTE) )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing (_TOTAL|_POTE) !!\n" );
#  else
   if (  !( TVar & (_TOTAL) )  )
      Aux_Error( ERROR_INFO, "no suitable target variable is found --> missing (_TOTAL) !!\n" );
#  endif

// nothing to do if there are no real patches at lv+1
   if ( amr->NPatchComma[SonLv][1] == 0 )    return;

// check the synchronization
   Mis_CompareRealValue( Time[FaLv], Time[SonLv], __FUNCTION__, true );


   const bool ResFlu    = TVar & _TOTAL;
#  ifdef GRAVITY
   const bool ResPot    = TVar & _POTE;
#  endif
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   const real  Gamma_m1 = GAMMA - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif

   int SonPID, FaPID, Disp_i, Disp_j, Disp_k, ii, jj, kk, I, J, K, Ip, Jp, Kp;
   int NVar_Flu, NVar_Tot, TFluVarIdx, TFluVarIdxList[NCOMP_TOTAL];


// determine the components to be restricted (TFluVarIdx : target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
   NVar_Flu = 0;

   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVar & (1<<v) )    TFluVarIdxList[ NVar_Flu++ ] = v;

// check again
   NVar_Tot = NVar_Flu;
#  ifdef GRAVITY
   if ( ResPot )  NVar_Tot ++;
#  endif
   if ( NVar_Tot == 0 )
   {
      Aux_Message( stderr, "WARNING : no target variable is found !!\n" );
      return;
   }


// restrict
#  pragma omp parallel for private( SonPID, FaPID, Disp_i, Disp_j, Disp_k, ii, jj, kk, I, J, K, Ip, Jp, Kp, \
                                    TFluVarIdx ) schedule( static )
   for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)
   {
      FaPID = amr->patch[0][SonLv][SonPID0]->father;

//    check
#     ifdef GAMER_DEBUG
      if ( FaPID < 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch (FaPID = %d) !!\n",
                    SonLv, SonPID0, FaPID );

      if ( ResFlu  &&  amr->patch[FaFluSg][FaLv][FaPID]->fluid == NULL )
         Aux_Error( ERROR_INFO, "FaFluSg %d, FaLv %d, FaPID %d has no fluid array allocated !!\n",
                    FaFluSg, FaLv, FaPID );

#     ifdef GRAVITY
      if ( ResPot  &&  amr->patch[FaPotSg][FaLv][FaPID]->pot == NULL )
         Aux_Error( ERROR_INFO, "FaPotSg %d, FaLv %d, FaPID %d has no potential array allocated !!\n",
                    FaPotSg, FaLv, FaPID );
#     endif
#     endif // #ifdef GAMER_DEBUG


//    loop over eight sons
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         SonPID = SonPID0 + LocalID;
         Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE/2 );
         Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE/2 );
         Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE/2 );

//       check
#        ifdef GAMER_DEBUG
         if ( ResFlu  &&  amr->patch[SonFluSg][SonLv][SonPID]->fluid == NULL )
            Aux_Error( ERROR_INFO, "SonFluSg %d, SonLv %d, SonPID %d has no fluid array allocated !!\n",
                       SonFluSg, SonLv, SonPID );

#        ifdef GRAVITY
         if ( ResPot  &&  amr->patch[SonPotSg][SonLv][SonPID]->pot == NULL )
            Aux_Error( ERROR_INFO, "SonPotSg %d, SonLv %d, SonPID %d has no potential array allocated !!\n",
                       SonPotSg, SonLv, SonPID );
#        endif
#        endif // #ifdef GAMER_DEBUG


//       restrict the fluid data
         if ( ResFlu )
         for (int v=0; v<NVar_Flu; v++)
         {
            TFluVarIdx = TFluVarIdxList[v];

            for (int k=0; k<PATCH_SIZE/2; k++)  {  K = k*2;    Kp = K+1;   kk = k + Disp_k;
            for (int j=0; j<PATCH_SIZE/2; j++)  {  J = j*2;    Jp = J+1;   jj = j + Disp_j;
            for (int i=0; i<PATCH_SIZE/2; i++)  {  I = i*2;    Ip = I+1;   ii = i + Disp_i;

               amr->patch[FaFluSg][FaLv][FaPID]->fluid[TFluVarIdx][kk][jj][ii]
                  = 0.125 * ( amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][J ][I ] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][J ][Ip] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][Jp][I ] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][J ][I ] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][Jp][Ip] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][Jp][I ] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][J ][Ip] +
                              amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][Jp][Ip]   );
            }}}
         }


#        ifdef GRAVITY
//       restrict the potential data
         if ( ResPot )
         {
            for (int k=0; k<PATCH_SIZE/2; k++)  {  K = k*2;    Kp = K+1;   kk = k + Disp_k;
            for (int j=0; j<PATCH_SIZE/2; j++)  {  J = j*2;    Jp = J+1;   jj = j + Disp_j;
            for (int i=0; i<PATCH_SIZE/2; i++)  {  I = i*2;    Ip = I+1;   ii = i + Disp_i;

               amr->patch[FaPotSg][FaLv][FaPID]->pot[kk][jj][ii]
                  = 0.125 * ( amr->patch[SonPotSg][SonLv][SonPID]->pot[K ][J ][I ] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[K ][J ][Ip] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[K ][Jp][I ] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[Kp][J ][I ] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[K ][Jp][Ip] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[Kp][Jp][I ] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[Kp][J ][Ip] +
                              amr->patch[SonPotSg][SonLv][SonPID]->pot[Kp][Jp][Ip]   );
            }}}
         }
#        endif // #ifdef GRAVITY
      } // for (int LocalID=0; LocalID<8; LocalID++)


//    check the minimum pressure and, when the dual-energy formalism is adopted, ensure the consistency between
//    pressure, total energy density, and the dual-energy variable
#     if ( MODEL == HYDRO  ||  MODEL == MHD )
//    apply this correction only when preparing all fluid variables
      if (  ( TVar & _TOTAL ) == _TOTAL  )
      for (int k=0; k<PATCH_SIZE; k++)
      for (int j=0; j<PATCH_SIZE; j++)
      for (int i=0; i<PATCH_SIZE; i++)
      {
#        ifdef DUAL_ENERGY
//       here we ALWAYS use the dual-energy variable to correct the total energy density
//       --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//           the runtime parameter DUAL_ENERGY_SWITCH here
         const bool CheckMinPres_Yes = true;
         const real UseEnpy2FixEngy  = HUGE_NUMBER;
         char dummy;    // we do not record the dual-energy status here

         CPU_DualEnergyFix( amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i],
                            amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMX][k][j][i],
                            amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMY][k][j][i],
                            amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMZ][k][j][i],
                            amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i],
                            amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENPY][k][j][i],
                            dummy, Gamma_m1, _Gamma_m1, CheckMinPres_Yes, MIN_PRES, UseEnpy2FixEngy );

#        else
//       actually it might not be necessary to check the minimum pressure here
         amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i]
            = CPU_CheckMinPresInEngy( amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i],
                                      amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMX][k][j][i],
                                      amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMY][k][j][i],
                                      amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMZ][k][j][i],
                                      amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i],
                                      Gamma_m1, _Gamma_m1, MIN_PRES );
#        endif // #ifdef DUAL_ENERGY ... else ...
      } // i,j,k
#     endif // #if ( MODEL == HYDRO  ||  MODEL == MHD )


//    rescale real and imaginary parts to get the correct density in ELBDM
#     if ( MODEL == ELBDM )
      real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

      if (  ( TVar & _DENS )  &&  ( TVar & _REAL )  &&  (TVar & _IMAG )  )
      for (int k=0; k<PATCH_SIZE; k++)
      for (int j=0; j<PATCH_SIZE; j++)
      for (int i=0; i<PATCH_SIZE; i++)
      {
         Real      = amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i];
         Imag      = amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i];
         Rho_Wrong = Real*Real + Imag*Imag;
         Rho_Corr  = amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
         {
            amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i] = (real)0.0;
            Rescale = (real)0.0;
         }
         else
            Rescale = SQRT( Rho_Corr/Rho_Wrong );

         amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i] *= Rescale;
         amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i] *= Rescale;
      }
#     endif

   } // for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)

} // FUNCTION : Flu_Restrict
