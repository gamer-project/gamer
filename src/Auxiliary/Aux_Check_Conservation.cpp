#include "Copyright.h"
#include "GAMER.h"

#if ( MODEL == MHD )
#warning : WAIT MHD !!!
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Conservation
// Description :  Verify the conservation laws
//                --> HYDRO : check mass, momenum, and energy
//                    MHD   : check mass, momenum, and energy
//                    ELBDM : check mass (and energy in some cases)
//
// Note        :  1. This check only works with the models HYDRO, MHD, ELBDM, and PAR_ONLY
//                2. The values measured during the first time this function is invoked will be taken as the
//                   reference values to estimate errors
//
// Parameter   :  comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Conservation( const char *comment )
{

   static bool FirstTime = true;
   const char *FileName  = "Record__Conservation";


// check
#  if ( MODEL != HYDRO  &&  MODEL != MHD  &&  MODEL != ELBDM  &&  MODEL != PAR_ONLY )
   Aux_Message( stderr, "Warning : function \"%s\" is supported only in the models HYDRO, MHD, ELBDM, and PAR_ONLY !!\n", 
                __FUNCTION__ );
   OPT__CK_CONSERVATION = false;
   return;
#  endif

   if ( FirstTime  &&  MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
   }


#  if   ( MODEL == HYDRO )
#  ifdef GRAVITY
   const int NVar = NCOMP+3;  // 8: mass, momentum (x/y/z), total/kinematic/thermal/potential energies
#  else
   const int NVar = NCOMP+2;
#  endif

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
#  if ( defined GRAVITY  &&  !defined COMOVING )
   const bool ELBDM_GetEngy    = ( OPT__BC_POT == BC_POT_ISOLATED );
#  else
   const bool ELBDM_GetEngy    = false;
#  endif
   const int  NVar             = ( ELBDM_GetEngy ) ? 5 : 1; // 5: mass, kinematic, gravitational, self-interaction, 
                                                            //    total energy

// useful only if ELBDM_GetEngy is true
   const bool IntPhase_No      = false;
   const bool GetTotDens_No    = false;
   const IntScheme_t IntScheme = INT_CQUAR;
   const int  NGhost           = 1;             // number of ghost zones for calculating the gradient of the wave function
   const int  Size_Flu         = PS1 + 2*NGhost;
   const int  NPG              = 1;
   const double _2Eta2         = 0.5/SQR(ELBDM_ETA);

   real (*Flu_Array)[2][Size_Flu][Size_Flu][Size_Flu] = NULL;
   real   GradR[3], GradI[3], _dh2;
   int    ip, im, jp, jm, kp, km;

   if ( ELBDM_GetEngy ) Flu_Array = new real [NPG*8][2][Size_Flu][Size_Flu][Size_Flu];

#  else
#  error : ERROR : unsupported MODEL !!
#  endif


   double dV, Total_local[NVar], Total_sum[NVar], Total_lv[NVar]; // dV : cell volume at each level
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   double EK, ET;
#  ifdef GRAVITY
   double EP;
#  endif
#  endif
   int    FluSg;
#  ifdef GRAVITY
   int    PotSg;
#  endif
   FILE  *File = NULL;


// initialize accumulative variables as zero
   for (int v=0; v<NVar; v++)    Total_local[v] = 0.0;


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {  
      for (int v=0; v<NVar; v++)    Total_lv[v] = 0.0;

      dV    = CUBE( amr->dh[lv] );
      FluSg = amr->FluSg[lv];
#     ifdef GRAVITY
      PotSg = amr->PotSg[lv];
#     endif
#     if ( MODEL == ELBDM )
      if ( ELBDM_GetEngy )
      _dh2  = 0.5/amr->dh[lv];
#     endif


      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {  
#        if ( MODEL == ELBDM)
         if ( ELBDM_GetEngy )
         Prepare_PatchData( lv, Time[lv], Flu_Array[0][0][0][0], NGhost, NPG, &PID0, _REAL|_IMAG,
                            IntScheme, UNIT_PATCH, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE, GetTotDens_No );
#        endif

         for (int PID=PID0; PID<PID0+8; PID++)
         {
//          only check the leaf patches
            if ( amr->patch[0][lv][PID]->son == -1 )
            {
#              if   ( MODEL == HYDRO )
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
                  for (int v=0; v<NCOMP; v++)
                  Total_lv[v]    += (double)amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

                  EK              = 0.5*( SQR(amr->patch[FluSg][lv][PID]->fluid[MOMX][k][j][i]) 
                                         +SQR(amr->patch[FluSg][lv][PID]->fluid[MOMY][k][j][i]) 
                                         +SQR(amr->patch[FluSg][lv][PID]->fluid[MOMZ][k][j][i]) )
                                       / amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];
                  ET              = amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i] - EK;
                  Total_lv[5]    += EK;
                  Total_lv[6]    += ET;

#                 ifdef GRAVITY   
                  EP              = 0.5*amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i]
                                       *amr->patch[PotSg][lv][PID]->pot        [k][j][i];
                  Total_lv[7]    += EP;
                  Total_lv[ENGY] += EP;
#                 endif
               } // i,j,k

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
//                [0] mass
                  Total_lv[0] += (double)amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];

#                 ifdef GRAVITY
                  if ( ELBDM_GetEngy )
                  {
//                   [2] potential energy in ELBDM
                     Total_lv[2] += 0.5*amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i]
                                       *amr->patch[PotSg][lv][PID]->pot        [k][j][i];

//                   [3] quartic self-interaction potential
#                    ifdef QUARTIC_SELF_INTERACTION
                     Total_lv[3] += 0.5*ELBDM_LAMBDA*SQR( amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i] );
#                    endif
                  }
#                 endif
               }

//             [1] kinematic energy in ELBDM
               if ( ELBDM_GetEngy )
               {
                  const int t = PID - PID0;

                  for (int k=NGhost; k<Size_Flu-NGhost; k++)   { kp = k+1; km = k-1;
                  for (int j=NGhost; j<Size_Flu-NGhost; j++)   { jp = j+1; jm = j-1;
                  for (int i=NGhost; i<Size_Flu-NGhost; i++)   { ip = i+1; im = i-1;

                     GradR[0] = _dh2*( Flu_Array[t][0][k ][j ][ip] - Flu_Array[t][0][k ][j ][im] );
                     GradR[1] = _dh2*( Flu_Array[t][0][k ][jp][i ] - Flu_Array[t][0][k ][jm][i ] );
                     GradR[2] = _dh2*( Flu_Array[t][0][kp][j ][i ] - Flu_Array[t][0][km][j ][i ] );

                     GradI[0] = _dh2*( Flu_Array[t][1][k ][j ][ip] - Flu_Array[t][1][k ][j ][im] );
                     GradI[1] = _dh2*( Flu_Array[t][1][k ][jp][i ] - Flu_Array[t][1][k ][jm][i ] );
                     GradI[2] = _dh2*( Flu_Array[t][1][kp][j ][i ] - Flu_Array[t][1][km][j ][i ] );

                     Total_lv[1] += _2Eta2*( SQR(GradR[0]) + SQR(GradR[1]) + SQR(GradR[2]) + 
                                             SQR(GradI[0]) + SQR(GradI[1]) + SQR(GradI[2])   );
                  }}}
               }

#              else
#              error : ERROR : unsupported MODEL !!

#              endif // MODEL

            } // if ( amr->patch[0][lv][PID]->son == -1 )
         } // for (int PID=PID0; PID<PID0+8; PID++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

#     if ( MODEL == ELBDM )
      if ( ELBDM_GetEngy )    Total_lv[4] = Total_lv[1] + Total_lv[2] + Total_lv[3];
#     endif


//    sum over NLEVEL levels
      for (int v=0; v<NVar; v++)    Total_local[v] += Total_lv[v]*dV;
   } // for (int lv=0; lv<NLEVEL; lv++)


// sum over all ranks
   MPI_Reduce( Total_local, Total_sum, NVar, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// calculate particles kinematic and potential energy
#  ifdef PARTICLE
   double Ek_Par, Ep_Par;

   Par_GetEnergy( Ek_Par, Ep_Par );
#  endif


// output
   if ( MPI_Rank == 0 )
   { 
//    a variable length array cannot have static storage duration
#     if ( MODEL == HYDRO  ||  MODEL == MHD )
      static double Total_ref[5]; // maximum size of "NVar" == 5 in ELBDM
#     elif ( MODEL == ELBDM )
      static double Total_ref[NCOMP+3];
#     endif
#     ifdef PARTICLE
      static double Ek_Par_ref, Ep_Par_ref;
#     endif
      double AbsErr[NVar], RelErr[NVar];
   
      if ( FirstTime )
      {  
//       record the reference values
         for (int v=0; v<NVar; v++)    Total_ref[v] = Total_sum[v];

#        ifdef PARTICLE
         Ek_Par_ref = Ek_Par;
         Ep_Par_ref = Ep_Par;
#        endif


//       output header
         FILE *File = fopen( FileName, "a" );

         Aux_Message( File, "# Ref time : %13.7e\n", Time[0] );
         Aux_Message( File, "# Ref step : %ld\n",    Step    );

#        if   ( MODEL == HYDRO  ||  MODEL == MHD )
         Aux_Message( File, "# Mass     : total HYDRO/MHD mass\n" );
         Aux_Message( File, "# MomX/Y/Z : total HYDRO/MHD momentum\n" );
         Aux_Message( File, "# E_Tot    : total HYDRO/MHD energy\n" );
         Aux_Message( File, "# E_K      : total HYDRO/MHD kinematic energy\n" );
         Aux_Message( File, "# E_T      : total HYDRO/MHD thermal energy\n" );
         Aux_Message( File, "# E_P      : total HYDRO/MHD potential energy\n" );

#        elif ( MODEL == ELBDM )
         Aux_Message( File, "# Mass     : total ELBDM mass\n" );
         Aux_Message( File, "# E_K      : total ELBDM kinematic energy\n" );
         Aux_Message( File, "# E_P      : total ELBDM potential energy\n" );
         Aux_Message( File, "# E_S      : total ELBDM self-interaction energy\n" );
         Aux_Message( File, "# E_Tot    : total ELBDM energy\n" );
#        endif

#        ifdef PARTICLE
         Aux_Message( File, "# E_K_Par  : total PARTICLE kinematic energy\n" );
         Aux_Message( File, "# E_P_Par  : total PARTICLE potential energy\n" );
         Aux_Message( File, "# E_All    : total energy in all models (HYDRO/MHD/ELBDM + PARTICLES)\n" );
#        endif

         Aux_Message( File, "\n" );
         Aux_Message( File, "# AbsErr   : absolute error --> (now - ref)\n" );
         Aux_Message( File, "# RelErr   : relative error --> (now - ref)/abs(ref)\n" );

         Aux_Message( File, "#-----------------------------------------------------------------------------" );
         Aux_Message( File, "------------------------------------------------------------------------------\n\n" );

         Aux_Message( File, "#%12s  %10s", "Time", "Step" );

#        if   ( MODEL == HYDRO  ||  MODEL == MHD )
         Aux_Message( File, "  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s",
                      "Mass", "Mass_AbsErr", "Mass_RelErr", "MomX", "MomX_AbsErr", "MomX_RelErr",
                      "MomY", "MomY_AbsErr", "MomY_RelErr", "MomZ", "MomZ_AbsErr", "MomZ_RelErr" );
         Aux_Message( File, "  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s",
                      "E_Tot", "E_Tot_AbsErr", "E_Tot_RelErr", "E_K", "E_K_AbsErr", "E_K_RelErr",
                      "E_T", "E_T_AbsErr", "E_T_RelErr" );
#        ifdef GRAVITY
         Aux_Message( File, "  %14s  %14s  %14s", "E_P", "E_P_AbsErr", "E_P_RelErr" );
#        endif

#        elif ( MODEL == ELBDM )
         Aux_Message( File, "  %14s  %14s  %14s", "Mass", "Mass_AbsErr", "Mass_RelErr" );
         if ( ELBDM_GetEngy ) {
         Aux_Message( File, "  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s",
                      "E_K", "E_K_AbsErr", "E_K_RelErr", "E_P", "E_P_AbsErr", "E_P_RelErr", "E_S", "E_S_AbsErr", "E_S_RelErr",
                      "E_Tot", "E_Tot_AbsErr", "E_Tot_RelErr" ); }
#        endif

#        ifdef PARTICLE
         Aux_Message( File, "  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s",
                      "E_K_Par", "E_K_Par_AbsErr", "E_K_Par_RelErr", "E_P_Par", "E_P_Par_AbsErr", "E_P_Par_RelErr",
                      "E_All", "E_All_AbsErr", "E_All_RelErr" );
#        endif

         Aux_Message( File, "\n" );

         fclose( File );
      } // if ( FirstTime )


//    calculate errors
      for (int v=0; v<NVar; v++)
      {
         AbsErr[v] = Total_sum[v] - Total_ref[v];
         RelErr[v] = AbsErr[v] / fabs(Total_ref[v]);
      }

#     ifdef PARTICLE
      double E_Total     = Ek_Par     + Ep_Par;
      double E_Total_ref = Ek_Par_ref + Ep_Par_ref;

#     if ( MODEL == ELBDM )
      if ( ELBDM_GetEngy )
#     endif
      {
         E_Total     += Total_sum[4];
         E_Total_ref += Total_ref[4];
      }
#     endif // #ifdef PARTICLE

   
//    output
      File = fopen( FileName, "a" );

      Aux_Message( File, "%13.7e  %10ld", Time[0], Step );

      for (int v=0; v<NVar; v++)
         Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Total_sum[v], AbsErr[v], RelErr[v] );

#     ifdef PARTICLE
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
                   Ek_Par, Ek_Par-Ek_Par_ref, (Ek_Par-Ek_Par_ref)/fabs(Ek_Par_ref),
                   Ep_Par, Ep_Par-Ep_Par_ref, (Ep_Par-Ep_Par_ref)/fabs(Ep_Par_ref) );

#     if ( MODEL == ELBDM )
      if ( ELBDM_GetEngy )
#     endif
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e",
                   E_Total, E_Total-E_Total_ref, (E_Total-E_Total_ref)/fabs(E_Total_ref) );
#     endif // #ifdef PARTICLE

      Aux_Message( File, "\n" );

      fclose( File ); 
   } // if ( MPI_Rank == 0 )


   if ( FirstTime )  FirstTime = false;


#  if ( MODEL == ELBDM )
   if ( Flu_Array != NULL )   delete [] Flu_Array;
#  endif

} // FUNCTION : Aux_Check_Conservation
