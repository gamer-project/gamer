#include "GAMER.h"


// global variable to store the ELBDM center-of-mass velocity
// --> declared in "Model_ELBDM/ELBDM_RemoveMotionCM.cpp"
#if ( MODEL == ELBDM )
extern double ELBDM_Vcm[3];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Conservation
// Description :  Verify the conservation laws
//                --> Mass, center of mass, momentum, angular momentum, energy, passive scalars, ...
//
// Note        :  1. This check only works with the models HYDRO, ELBDM, and PAR_ONLY
//                2. The values measured during the first function call will be taken as the reference values
//                   to estimate errors
//                3. For simulations with particles (i.e., when PARTICLE is on), the total conserved variables
//                   (e.g., total energy of gas and particles) will also be recorded
//
// Parameter   :  comment : You can put the location where this function is invoked in this string
//                          (not used currently)
//
// Return      :  Log file "Record__Conservation"
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Conservation( const char *comment )
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Conservation", OUTPUT_DIR );


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM  &&  MODEL != PAR_ONLY )
   Aux_Message( stderr, "WARNING : function \"%s\" is supported only in the models HYDRO, ELBDM, and PAR_ONLY !!\n",
                __FUNCTION__ );
   OPT__CK_CONSERVATION = false;
   return;
#  endif


   if ( FirstTime  &&  MPI_Rank == 0 )
   {
//    check
#     ifdef COMOVING
      Aux_Message( stderr, "WARNING : function \"%s\" is NOT fully supported in COMOVING !! Only mass conservation check works !!\n", __FUNCTION__ );
#     endif

      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
   }


#  if   ( MODEL == HYDRO )
#  ifdef MHD
   const int    NVar_NoPassive    = 12;   // 12: mass, momentum (x/y/z), angular momentum (x/y/z), kinetic/internal/potential/magnetic/total energies
                                          // --> note that **total energy** is put in the last element
   const char   FluLabel[NVar_NoPassive][MAX_STRING] = { "Mass_Gas", "MomX_Gas", "MomY_Gas", "MomZ_Gas",
                                                         "AngMomX_Gas", "AngMomY_Gas", "AngMomZ_Gas",
                                                         "Ekin_Gas", "Eint_Gas", "Epot_Gas", "Emag_Gas",
                                                         "Etot_Gas"
                                                       };
#  else
   const int    NVar_NoPassive    = 11;   // 11: mass, momentum (x/y/z), angular momentum (x/y/z), kinetic/internal/potential/total energies
                                          // --> note that **total energy** is put in the last element
   const char   FluLabel[NVar_NoPassive][MAX_STRING] = { "Mass_Gas", "MomX_Gas", "MomY_Gas", "MomZ_Gas",
                                                         "AngMomX_Gas", "AngMomY_Gas", "AngMomZ_Gas",
                                                         "Ekin_Gas", "Eint_Gas", "Epot_Gas", "Etot_Gas"
                                                       };
#  endif
   const char   FluCoMLabel[3][MAX_STRING] = { "CoMX_Gas", "CoMY_Gas", "CoMZ_Gas" };
   const int    idx_etot_flu      = NVar_NoPassive - 1;
   const bool   CheckMinEint_No   = false;

#  elif ( MODEL == ELBDM )
   const int    NVar_NoPassive    = 11;   // 11: mass, momentum (x/y/z), angular momentum (x/y/z), kinetic/gravitational/self-interaction/total energies
   const int    idx_etot_flu      = NVar_NoPassive - 1;
   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const int    NGhost            = 1;    // number of ghost zones for calculating the gradient of wave function
   const int    Size_Flu          = PS1 + 2*NGhost;
   const int    NPG               = 1;
   const double _Eta              = 1.0/ELBDM_ETA;
   const double _2Eta2            = 0.5/SQR(ELBDM_ETA);
   const IntScheme_t IntScheme    = INT_CQUAR;
   const char   FluLabel[NVar_NoPassive][MAX_STRING] = { "Mass_Psi", "MomX_Psi", "MomY_Psi", "MomZ_Psi",
                                                         "AngMomX_Psi", "AngMomY_Psi", "AngMomZ_Psi",
                                                         "Ekin_Psi", "Epot_Psi", "Esel_Psi", "Etot_Psi"
                                                       };
   const char   FluCoMLabel[3][MAX_STRING] = { "CoMX_Psi", "CoMY_Psi", "CoMZ_Psi" };

   real (*Flu_ELBDM)[2][Size_Flu][Size_Flu][Size_Flu] = new real [NPG*8][2][Size_Flu][Size_Flu][Size_Flu];

#  else
#  error : ERROR : unsupported MODEL !!
#  endif


// get the sum of passive scalars to be normalized
   const bool GetPassiveSum      = ( PassiveNorm_NVar > 0 );
   const int  NVar_Flu           = NVar_NoPassive + NCOMP_PASSIVE + ( (GetPassiveSum)?1:0 );
   const int  idx_offset_flu     = 1;
   const int  idx_offset_flu_com = idx_offset_flu + NVar_Flu;

   int NStoredConRef_noTime = 0;
   NStoredConRef_noTime += NVar_Flu + 3; // +3: center-of-mass position

   double dh, dv, Fluid_ThisRank[NVar_Flu], Fluid_AllRank[NVar_Flu], Fluid_lv[NVar_Flu];   // dv : cell volume at each level
   int    FluSg;
#  ifdef GRAVITY
   int    PotSg;
#  endif
#  ifdef MHD
   int    MagSg;
#  endif
   FILE  *File = NULL;


// initialize accumulative variables as zero
   for (int v=0; v<NVar_Flu; v++)    Fluid_ThisRank[v] = 0.0;


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int v=0; v<NVar_Flu; v++)    Fluid_lv[v] = 0.0;

      dh    = amr->dh[lv];
      dv    = CUBE( amr->dh[lv] );
      FluSg = amr->FluSg[lv];
#     ifdef GRAVITY
      PotSg = amr->PotSg[lv];
#     endif
#     ifdef MHD
      MagSg = amr->MagSg[lv];
#     endif
#     if ( MODEL == ELBDM )
      const real _dh2  = 0.5/amr->dh[lv];
#     endif


      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
#        if ( MODEL == ELBDM )
         const real MinDens_No = -1.0;
         const real MinPres_No = -1.0;
         const real MinTemp_No = -1.0;
         const real MinEntr_No = -1.0;
         long  TVar = -1.0;

#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         if ( amr->use_wave_flag[lv] ) {
#        endif
         TVar = _REAL | _IMAG;
#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         } else {
         TVar = _DENS | _PHAS;
         }
#        endif

         Prepare_PatchData( lv, Time[lv], Flu_ELBDM[0][0][0][0], NULL, NGhost, NPG, &PID0, TVar, _NONE,
                            IntScheme, INT_NONE, UNIT_PATCH, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
#        endif // #if ( MODEL == ELBDM )

         for (int PID=PID0; PID<PID0+8; PID++)
         {
            if ( amr->patch[0][lv][PID]->son != -1 )   continue; // only check the leaf patches

            const double x0  = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
            const double y0  = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
            const double z0  = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

#           if   ( MODEL == HYDRO )
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
               double Dens, MomX, MomY, MomZ, AngMomX, AngMomY, AngMomZ, Etot, Ekin, Eint;
#              ifdef GRAVITY
               double Epot;
#              endif
               double Emag = NULL_REAL;

               Dens = amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];
               MomX = amr->patch[FluSg][lv][PID]->fluid[MOMX][k][j][i];
               MomY = amr->patch[FluSg][lv][PID]->fluid[MOMY][k][j][i];
               MomZ = amr->patch[FluSg][lv][PID]->fluid[MOMZ][k][j][i];
               Etot = amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i];

//             calculate the angular momentum
               const double x  = x0 + i*dh;
               const double y  = y0 + j*dh;
               const double z  = z0 + k*dh;

               const double dX = x - ANGMOM_ORIGIN_X;
               const double dY = y - ANGMOM_ORIGIN_Y;
               const double dZ = z - ANGMOM_ORIGIN_Z;

               AngMomX = dY*MomZ - dZ*MomY;
               AngMomY = dZ*MomX - dX*MomZ;
               AngMomZ = dX*MomY - dY*MomX;

#              ifdef SRHD
//             total energy density also includes rest mass energy density in relativistic hydro
               Etot += Dens;
#              endif

               Fluid_lv[0] += Dens;
               Fluid_lv[1] += MomX;
               Fluid_lv[2] += MomY;
               Fluid_lv[3] += MomZ;

               Fluid_lv[4] += AngMomX;
               Fluid_lv[5] += AngMomY;
               Fluid_lv[6] += AngMomZ;

#              ifdef MHD
               Emag          = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
               Fluid_lv[10] += Emag;
#              endif

#              ifdef GRAVITY
//             set potential energy to zero when enabling both OPT__SELF_GRAVITY and OPT__EXT_POT
//             since the potential energy obtained here would be wrong anyway
//             --> to avoid possible misinterpretation
               if      (  OPT__SELF_GRAVITY  &&  !OPT__EXT_POT )  Epot = 0.5*Dens*amr->patch[PotSg][lv][PID]->pot[k][j][i];
               else if ( !OPT__SELF_GRAVITY  &&   OPT__EXT_POT )  Epot =     Dens*amr->patch[PotSg][lv][PID]->pot[k][j][i];
               else                                               Epot = 0.0;
               Fluid_lv[9] += Epot;
#              endif
#              ifndef SRHD
//             Hydro_Con2Eint() calculates Eint for both HD and SRHD but we disable SRHD for now
               Eint         = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Etot, CheckMinEint_No, NULL_REAL, Emag,
                                              EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                              EoS_AuxArray_Int, h_EoS_Table );
#              else
               Eint = 0.0;
#              endif
               Fluid_lv[8] += Eint;

#              ifdef SRHD
//             For now we disable the calculation of Ekin for SRHD
//             Also, note that the following is equivalent to "Etot - Dens - Lrtz*Eint"
               /*
               real HTilde, Prim[NCOMP_TOTAL], Cons[NCOMP_TOTAL], Lrtz, Lrtz_m1;
               Cons[0]      = Dens;
               Cons[1]      = MomX;
               Cons[2]      = MomY;
               Cons[3]      = MomZ;
               Cons[4]      = Etot;
               for ( int v = NCOMP_FLUID; v < NCOMP_TOTAL; v++ ) Cons[v] = 0.0;
               Hydro_Con2Pri( Cons, Prim, (real)-HUGE_NUMBER, NULL_BOOL, NULL_INT, NULL,
                              NULL_BOOL, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                              EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL, &Lrtz );
               HTilde       = Hydro_Con2HTilde( Cons, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

//             Compute gamma - 1 this way to avoid catastrophic cancellation
               Lrtz_m1      = ( SQR(Prim[1]) + SQR(Prim[2]) + SQR(Prim[3]) ) / ( Lrtz + 1.0 );
               Ekin         = Lrtz_m1*( Dens*(HTilde+1.0) + Prim[4] );
               */
               Ekin = 0.0;
#              else
//###NOTE: assuming Etot = Eint + Ekin + Emag
               Ekin         = Etot - Eint;
#              ifdef MHD
               Ekin        -= Emag;
#              endif
#              endif
               Fluid_lv[7] += Ekin;
            } // i,j,k


#           elif ( MODEL == ELBDM )
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
               double Dens, Esel;
#              ifdef GRAVITY
               double Epot;
#              endif

//             [0] mass
               Dens         = amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];
               Fluid_lv[0] += Dens;

//             [8] potential energy in ELBDM
#              ifdef GRAVITY
//             set potential energy to zero when enabling both OPT__SELF_GRAVITY and OPT__EXT_POT
//             since the potential energy obtained here would be wrong anyway
//             --> to avoid possible misinterpretation
               if      (  OPT__SELF_GRAVITY  &&  !OPT__EXT_POT )  Epot = 0.5*Dens*amr->patch[PotSg][lv][PID]->pot[k][j][i];
               else if ( !OPT__SELF_GRAVITY  &&   OPT__EXT_POT )  Epot =     Dens*amr->patch[PotSg][lv][PID]->pot[k][j][i];
               else                                               Epot = 0.0;
               Fluid_lv[8] += Epot;
#              endif

//             [9] quartic self-interaction potential in ELBDM
#              ifdef QUARTIC_SELF_INTERACTION
               Esel         = 0.5*ELBDM_LAMBDA*SQR( amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i] );
               Fluid_lv[9] += Esel;
#              endif
            }

            const int t = PID - PID0;

            for (int k=NGhost; k<Size_Flu-NGhost; k++)   { const int kp = k+1; const int km = k-1;
            for (int j=NGhost; j<Size_Flu-NGhost; j++)   { const int jp = j+1; const int jm = j-1;
            for (int i=NGhost; i<Size_Flu-NGhost; i++)   { const int ip = i+1; const int im = i-1;

               double MomX, MomY, MomZ, AngMomX, AngMomY, AngMomZ, Ekin;

#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               if ( amr->use_wave_flag[lv] ) {
#              endif
               const real R       = Flu_ELBDM[t][0][k][j][i];
               const real I       = Flu_ELBDM[t][1][k][j][i];

//             compute gradient of real part dR/dx
               const real GradR_X = _dh2*( Flu_ELBDM[t][0][k ][j ][ip] - Flu_ELBDM[t][0][k ][j ][im] );
               const real GradR_Y = _dh2*( Flu_ELBDM[t][0][k ][jp][i ] - Flu_ELBDM[t][0][k ][jm][i ] );
               const real GradR_Z = _dh2*( Flu_ELBDM[t][0][kp][j ][i ] - Flu_ELBDM[t][0][km][j ][i ] );

//             compute gradient of imaginary part dI/dx
               const real GradI_X = _dh2*( Flu_ELBDM[t][1][k ][j ][ip] - Flu_ELBDM[t][1][k ][j ][im] );
               const real GradI_Y = _dh2*( Flu_ELBDM[t][1][k ][jp][i ] - Flu_ELBDM[t][1][k ][jm][i ] );
               const real GradI_Z = _dh2*( Flu_ELBDM[t][1][kp][j ][i ] - Flu_ELBDM[t][1][km][j ][i ] );

//             compute momentum in ELBDM wave scheme
               MomX = _Eta*( R*GradI_X - I*GradR_X );
               MomY = _Eta*( R*GradI_Y - I*GradR_Y );
               MomZ = _Eta*( R*GradI_Z - I*GradR_Z );

//             compute kinetic energy in ELBDM wave scheme
               Ekin = _2Eta2*( SQR(GradR_X) + SQR(GradR_Y) + SQR(GradR_Z) +
                               SQR(GradI_X) + SQR(GradI_Y) + SQR(GradI_Z)   );

#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               } else {
               const double Dens    = Flu_ELBDM[t][DENS][k][j][i];

//             compute bulk velocities v_i = (1/Eta)*dS/dx
               const double Vbulk_X =  _Eta * ( _dh2*( Flu_ELBDM[t][PHAS][k ][j ][ip] - Flu_ELBDM[t][PHAS][k ][j ][im] ) );
               const double Vbulk_Y =  _Eta * ( _dh2*( Flu_ELBDM[t][PHAS][k ][jp][i ] - Flu_ELBDM[t][PHAS][k ][jm][i ] ) );
               const double Vbulk_Z =  _Eta * ( _dh2*( Flu_ELBDM[t][PHAS][kp][j ][i ] - Flu_ELBDM[t][PHAS][km][j ][i ] ) );

//             compute thermal velocities v_r = (1/Eta)*dln(sqrt(rho))/dx = (1/Eta)*0.5*dln(rho)/dx
               const double Vther_X =  _Eta * ( (real)0.5*_dh2*( LOG(Flu_ELBDM[t][DENS][k ][j ][ip]) - LOG(Flu_ELBDM[t][DENS][k ][j ][im]) ) );
               const double Vther_Y =  _Eta * ( (real)0.5*_dh2*( LOG(Flu_ELBDM[t][DENS][k ][jp][i ]) - LOG(Flu_ELBDM[t][DENS][k ][jm][i ]) ) );
               const double Vther_Z =  _Eta * ( (real)0.5*_dh2*( LOG(Flu_ELBDM[t][DENS][kp][j ][i ]) - LOG(Flu_ELBDM[t][DENS][km][j ][i ]) ) );

//             compute momentum in ELBDM fluid scheme
               MomX = Dens * Vbulk_X;
               MomY = Dens * Vbulk_Y;
               MomZ = Dens * Vbulk_Z;

//             compute kinetic energy in ELBDM fluid scheme
               Ekin = 0.5 * Dens * ( SQR(Vbulk_X) + SQR(Vbulk_Y) + SQR(Vbulk_Z) +
                                     SQR(Vther_X) + SQR(Vther_Y) + SQR(Vther_Z)   );

               } // if ( amr->use_wave_flag[lv] ) ... else ...
#              endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )

//             compute angular momentum in ELBDM
               const double x  = x0 + (i-NGhost)*dh;
               const double y  = y0 + (j-NGhost)*dh;
               const double z  = z0 + (k-NGhost)*dh;

               const double dX = x - ANGMOM_ORIGIN_X;
               const double dY = y - ANGMOM_ORIGIN_Y;
               const double dZ = z - ANGMOM_ORIGIN_Z;

               AngMomX = dY*MomZ - dZ*MomY;
               AngMomY = dZ*MomX - dX*MomZ;
               AngMomZ = dX*MomY - dY*MomX;

//             [1-3] momentum in ELBDM
               Fluid_lv[1] += MomX;
               Fluid_lv[2] += MomY;
               Fluid_lv[3] += MomZ;

//             [4-6] angular momentum in ELBDM
               Fluid_lv[4] += AngMomX;
               Fluid_lv[5] += AngMomY;
               Fluid_lv[6] += AngMomZ;

//             [7] kinetic energy in ELBDM
               Fluid_lv[7] += Ekin;

            }}} // i,j,k

#           else
#           error : ERROR : unsupported MODEL !!

#           endif // MODEL


//          individual passive scalars
            for (int v=0; v<NCOMP_PASSIVE; v++)
            {
               const int v1 = NVar_NoPassive + v;
               const int v2 = NCOMP_FLUID    + v;

               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
                  Fluid_lv[v1] += amr->patch[FluSg][lv][PID]->fluid[v2][k][j][i];
            }
         } // for (int PID=PID0; PID<PID0+8; PID++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

//    get the total energy
#     if   ( MODEL == HYDRO )
#     ifdef MHD
      Fluid_lv[idx_etot_flu] = Fluid_lv[7] + Fluid_lv[8] + Fluid_lv[9] + Fluid_lv[10];
#     else
      Fluid_lv[idx_etot_flu] = Fluid_lv[7] + Fluid_lv[8] + Fluid_lv[9];
#     endif
#     elif ( MODEL == ELBDM )
      Fluid_lv[idx_etot_flu] = Fluid_lv[7] + Fluid_lv[8] + Fluid_lv[9];
#     else
#     error : ERROR : unsupported MODEL !!
#     endif

//    sum of passive scalars to be normalized
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<PassiveNorm_NVar; v++)
         Fluid_lv[ NVar_Flu - 1 ] += Fluid_lv[ NVar_NoPassive + PassiveNorm_VarIdx[v] ];
#     endif

//    multiply by the cell volume and sum over all levels
      for (int v=0; v<NVar_Flu; v++)    Fluid_ThisRank[v] += Fluid_lv[v]*dv;
   } // for (int lv=0; lv<NLEVEL; lv++)


// sum over all ranks
   MPI_Reduce( Fluid_ThisRank, Fluid_AllRank, NVar_Flu, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// compute the center of mass
   double CoM_Flu[3];
   double FinaldR_Flu;
   int    FinalNIter_Flu;
   Aux_FindWeightedAverageCenter( CoM_Flu, amr->BoxCenter, __FLT_MAX__, 0.0, _DENS, __FLT_MAX__, 1, &FinaldR_Flu, &FinalNIter_Flu );


// calculate conserved quantities for particles
#  ifdef MASSIVE_PARTICLES
   const int  NVar_Par           = 10; // 10: mass, momentum (x/y/z), angular momentum (x/y/z), kinetic/potential/total energies
   const int  idx_etot_par       = NVar_Par - 1;
   const int  idx_offset_par     = idx_offset_flu_com + 3;
   const int  idx_offset_par_com = idx_offset_par + NVar_Par;
   const char ParLabel[NVar_Par][MAX_STRING] = { "Mass_Par", "MomX_Par", "MomY_Par", "MomZ_Par", "AngMomX_Par",
                                                 "AngMomY_Par", "AngMomZ_Par", "Ekin_Par", "Epot_Par",
                                                 "Etot_Par"
                                               };
   const char ParCoMLabel[3][MAX_STRING] = { "CoMX_Par", "CoMY_Par", "CoMZ_Par" };
   double Par_AllRank[NVar_Par];
   double CoM_Par[3];

   NStoredConRef_noTime += NVar_Par + 3; // +3: center-of-mass position

   Par_Aux_GetConservedQuantity( Par_AllRank[0], CoM_Par[0], CoM_Par[1], CoM_Par[2],
                                 Par_AllRank[1], Par_AllRank[2], Par_AllRank[3],
                                 Par_AllRank[4], Par_AllRank[5], Par_AllRank[6], Par_AllRank[7], Par_AllRank[8] );

   Par_AllRank[idx_etot_par] = Par_AllRank[7] + Par_AllRank[8];
#  endif

// All = fluid + particles
#  if ( defined MASSIVE_PARTICLES  &&  MODEL != PAR_ONLY )
   const int  NVar_All           = 8; // 8: mass, momentum (x/y/z), angular momentum (x/y/z), total energy
   const int  idx_etot_all       = NVar_All - 1;
   const int  idx_offset_all     = idx_offset_par_com + 3;
   const int  idx_offset_all_com = idx_offset_all + NVar_All;
   const char AllLabel[NVar_All][MAX_STRING] = { "Mass_All", "MomX_All", "MomY_All", "MomZ_All", "AngMomX_All",
                                                 "AngMomY_All", "AngMomZ_All", "Etot_All"
                                               };
   const char AllCoMLabel[3][MAX_STRING] = { "CoMX_All", "CoMY_All", "CoMZ_All" };
   double All_AllRank[NVar_Par];
   double CoM_All[3];

   NStoredConRef_noTime += NVar_All + 3; // +3: center-of-mass position
#  endif // if ( defined MASSIVE_PARTICLES  &&  MODEL != PAR_ONLY )


// record the reference values of conserved variables
   if ( MPI_Rank == 0 )
   {
//    calculate the sum of conserved quantities in different models
#     if ( defined MASSIVE_PARTICLES  &&  MODEL != PAR_ONLY )
      for (int v=0; v<7; v++)   All_AllRank[v] = Fluid_AllRank[v] + Par_AllRank[v]; // 0-6: mass, momentum x/y/z, angular momentum x/y/z
      All_AllRank[idx_etot_all] = Fluid_AllRank[idx_etot_flu] + Par_AllRank[idx_etot_par]; // for HYDRO/ELBDM, total energy is stored in the last element

      for (int d=0; d<3; d++)
         CoM_All[d] = ( Fluid_AllRank[0]*CoM_Flu[d] + Par_AllRank[0]*CoM_Par[d] )/All_AllRank[0];
#     endif // if ( defined MASSIVE_PARTICLES  &&  MODEL != PAR_ONLY )

//    record the reference values if not initialized, e.g., first time or restart from an HDF5 snapshot with version < 2502
      if ( ! ConRefInitialized )
      {
         if ( NStoredConRef_noTime > NCONREF_MAX )
            Aux_Error( ERROR_INFO, "exceed NCOMREF_MAX (%d) !!\n", NCONREF_MAX );

         for (int v=0; v<1+NCONREF_MAX; v++)   ConRef[v] = NULL_REAL;

         ConRef[0] = Time[0];
         for (int v=0; v<NVar_Flu; v++)   ConRef[idx_offset_flu    +v] = Fluid_AllRank[v];
         for (int d=0; d<3; d++)          ConRef[idx_offset_flu_com+d] = CoM_Flu[d];

#        ifdef MASSIVE_PARTICLES
         for (int v=0; v<NVar_Par; v++)   ConRef[idx_offset_par    +v] = Par_AllRank[v];
         for (int d=0; d<3; d++)          ConRef[idx_offset_par_com+d] = CoM_Par[d];
#        if ( MODEL != PAR_ONLY )
         for (int v=0; v<NVar_All; v++)   ConRef[idx_offset_all    +v] = All_AllRank[v];
         for (int d=0; d<3; d++)          ConRef[idx_offset_all_com+d] = CoM_All[d];
#        endif // #if ( MODEL != PAR_ONLY )
#        endif // #ifdef MASSIVE_PARTICLES
      } // if ( ! ConRefInitialized )
   } // if ( MPI_Rank == 0 )

   ConRefInitialized = true;


// only record the reference values when conservation check is disabled
   if ( ! OPT__CK_CONSERVATION )
   {
#     if ( MODEL == ELBDM )
      delete [] Flu_ELBDM;
#     endif
      return;
   }


// output
   if ( MPI_Rank == 0 )
   {
      const int index_before_column_CoM = 0;

      double AbsErr_Flu[NVar_Flu], RelErr_Flu[NVar_Flu], AbsErr_CoM_Flu[3], AveVel_CoM_Flu[3];
#     ifdef MASSIVE_PARTICLES
      double AbsErr_Par[NVar_Par], RelErr_Par[NVar_Par], AbsErr_CoM_Par[3], AveVel_CoM_Par[3];
#     if ( MODEL != PAR_ONLY )
      double AbsErr_All[NVar_All], RelErr_All[NVar_All], AbsErr_CoM_All[3], AveVel_CoM_All[3];
#     endif // #if ( MODEL != PAR_ONLY )
#     endif // #ifdef MASSIVE_PARTICLES

      if ( FirstTime )
      {
//       output header
         FILE *File = fopen( FileName, "a" );

         Aux_Message( File, "# Ref time        : %13.7e\n", ConRef[0] );
         Aux_Message( File, "\n" );

#        if   ( MODEL == HYDRO )
         Aux_Message( File, "# Mass_Gas        : total HYDRO mass\n" );
         Aux_Message( File, "# CoMX/Y/Z_Gas    : total HYDRO center of mass\n" );
         Aux_Message( File, "# MomX/Y/Z_Gas    : total HYDRO momentum\n" );
         Aux_Message( File, "# AngMomX/Y/Z_Gas : total HYDRO angular momentum\n" );
         Aux_Message( File, "# Ekin_Gas        : total HYDRO kinetic energy\n" );
         Aux_Message( File, "# Eint_Gas        : total HYDRO internal energy\n" );
         Aux_Message( File, "# Epot_Gas        : total HYDRO potential energy\n" );
#        ifdef MHD
         Aux_Message( File, "# Emag_Gas        : total HYDRO magnetic energy\n" );
#        endif
         Aux_Message( File, "# Etot_Gas        : total HYDRO energy\n" );

#        elif ( MODEL == ELBDM )
         Aux_Message( File, "# Mass_Psi        : total ELBDM mass\n" );
         Aux_Message( File, "# CoMX/Y/Z_Psi    : total ELBDM center of mass\n" );
         Aux_Message( File, "# MomX/Y/Z_Psi    : total ELBDM momentum\n" );
         Aux_Message( File, "# AngMomX/Y/Z_Psi : total ELBDM angular momentum\n" );
         Aux_Message( File, "# Ekin_Psi        : total ELBDM kinetic energy\n" );
         Aux_Message( File, "# Epot_Psi        : total ELBDM potential energy\n" );
         Aux_Message( File, "# Esel_Psi        : total ELBDM self-interaction energy\n" );
         Aux_Message( File, "# Etot_Psi        : total ELBDM energy\n" );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif

         if ( GetPassiveSum )
         Aux_Message( File, "# PassNorm        : sum of all target passive scalars to be normalized\n" );

#        ifdef MASSIVE_PARTICLES
         Aux_Message( File, "# Mass_Par        : total PARTICLE mass\n" );
         Aux_Message( File, "# CoMX/Y/Z_Par    : total PARTICLE center of mass\n" );
         Aux_Message( File, "# MomX/Y/Z_Par    : total PARTICLE momentum\n" );
         Aux_Message( File, "# AngMomX/Y/Z_Par : total PARTICLE angular momentum\n" );
         Aux_Message( File, "# Ekin_Par        : total PARTICLE kinetic energy\n" );
         Aux_Message( File, "# Epot_Par        : total PARTICLE potential energy\n" );
         Aux_Message( File, "# Etot_Par        : total PARTICLE energy\n" );

#        if ( MODEL != PAR_ONLY )
         Aux_Message( File, "# Mass_All        : sum of the total HYDRO/ELBDM + PARTICLE mass\n" );
         Aux_Message( File, "# CoMX/Y/Z_All    : total HYDRO/ELBDM + PARTICLE center of mass\n" );
         Aux_Message( File, "# MomX_All        : sum of the total HYDRO/ELBDM + PARTICLE momentum x\n" );
         Aux_Message( File, "# MomY_All        : sum of the total HYDRO/ELBDM + PARTICLE momentum y\n" );
         Aux_Message( File, "# MomZ_All        : sum of the total HYDRO/ELBDM + PARTICLE momentum z\n" );
         Aux_Message( File, "# AngMomX_All     : sum of the total HYDRO/ELBDM + PARTICLE angular momentum x\n" );
         Aux_Message( File, "# AngMomY_All     : sum of the total HYDRO/ELBDM + PARTICLE angular momentum y\n" );
         Aux_Message( File, "# AngMomZ_All     : sum of the total HYDRO/ELBDM + PARTICLE angular momentum z\n" );
         Aux_Message( File, "# Etot_All        : sum of the total HYDRO/ELBDM + PARTICLE energy\n" );
#        endif // if ( MODEL != PAR_ONLY )
#        endif // #ifdef MASSIVE_PARTICLES

         Aux_Message( File, "\n" );
         Aux_Message( File, "# AErr            : absolute error --> (now - ref)\n" );
         Aux_Message( File, "# RErr            : relative error --> (now - ref) / abs(ref)\n" );
         Aux_Message( File, "# AveV            : average velocity of CoM --> (now - ref) / (time - time_ref)\n" );

         Aux_Message( File, "#-------------------------------------------------------------------------------------------" );
         Aux_Message( File, "--------------------------------------------------------------------------------------------\n\n" );

         const int NColumnMax = 140;
         Aux_Message( File, "#%12s  %10s", "[  1]", "[  2]" );
         for (int c=2; c<NColumnMax; c++)    Aux_Message( File, "  %12s[%3d]", "", c+1 );
         Aux_Message( File, "\n" );

         Aux_Message( File, "#%12s  %10s", "Time", "Step" );

         for (int v=0; v<NVar_NoPassive; v++)
         {
         Aux_Message( File, "  %17s  %12s_AErr  %12s_RErr", FluLabel[v], FluLabel[v], FluLabel[v] );

         if ( v == index_before_column_CoM )
         {
         for (int d=0; d<3; d++)
         Aux_Message( File, "  %17s  %12s_AErr  %12s_AveV", FluCoMLabel[d], FluCoMLabel[d], FluCoMLabel[d] );
         }
         } // for (int v=0; v<NVar_NoPassive; v++)

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         Aux_Message( File, "  %17s  %12s_AErr  %12s_RErr", FieldLabel[v], FieldLabel[v], FieldLabel[v] );

         if ( GetPassiveSum )
         Aux_Message( File, "  %17s  %17s  %17s",    "PassNorm",    "PassNorm_AErr",    "PassNorm_RErr" );

#        ifdef MASSIVE_PARTICLES
         for (int v=0; v<NVar_Par; v++)
         {
         Aux_Message( File, "  %17s  %12s_AErr  %12s_RErr", ParLabel[v], ParLabel[v], ParLabel[v] );

         if ( v == index_before_column_CoM )
         {
         for (int d=0; d<3; d++)
         Aux_Message( File, "  %17s  %12s_AErr  %12s_AveV", ParCoMLabel[d], ParCoMLabel[d], ParCoMLabel[d] );
         }
         } // for (int v=0; v<NVar_Par; v++)

#        if ( MODEL != PAR_ONLY )
         for (int v=0; v<NVar_All; v++)
         {
         Aux_Message( File, "  %17s  %12s_AErr  %12s_RErr", AllLabel[v], AllLabel[v], AllLabel[v] );

         if ( v == index_before_column_CoM )
         {
         for (int d=0; d<3; d++)
         Aux_Message( File, "  %17s  %12s_AErr  %12s_AveV", AllCoMLabel[d], AllCoMLabel[d], AllCoMLabel[d] );
         }
         } // for (int v=0; v<NVar_All; v++)
#        endif // if ( MODEL != PAR_ONLY )
#        endif // #ifdef PARTICLE

         Aux_Message( File, "\n" );

         fclose( File );
      } // if ( FirstTime )


//    calculate errors
      for (int v=0; v<NVar_Flu; v++)
      {
         AbsErr_Flu[v] = Fluid_AllRank[v] - ConRef[idx_offset_flu+v];
         RelErr_Flu[v] = AbsErr_Flu[v] / fabs(ConRef[idx_offset_flu+v]);
      }
      for (int d=0; d<3; d++)
      {
         AbsErr_CoM_Flu[d] = CoM_Flu[d] - ConRef[idx_offset_flu_com+d];
         AveVel_CoM_Flu[d] = AbsErr_CoM_Flu[d] / (Time[0]-ConRef[0]);
      }

#     ifdef MASSIVE_PARTICLES
      for (int v=0; v<NVar_Par; v++)
      {
         AbsErr_Par[v] = Par_AllRank[v] - ConRef[idx_offset_par+v];
         RelErr_Par[v] = AbsErr_Par[v] / fabs(ConRef[idx_offset_par+v]);
      }
      for (int d=0; d<3; d++)
      {
         AbsErr_CoM_Par[d] = CoM_Par[d] - ConRef[idx_offset_par_com+d];
         AveVel_CoM_Par[d] = AbsErr_CoM_Par[d] / (Time[0]-ConRef[0]);
      }

#     if ( MODEL != PAR_ONLY )
      for (int v=0; v<NVar_All; v++)
      {
         AbsErr_All[v] = All_AllRank[v] - ConRef[idx_offset_all+v];
         RelErr_All[v] = AbsErr_All[v] / fabs(ConRef[idx_offset_all+v]);
      }
      for (int d=0; d<3; d++)
      {
         AbsErr_CoM_All[d] = CoM_All[d] - ConRef[idx_offset_all_com+d];
         AveVel_CoM_All[d] = AbsErr_CoM_All[d] / (Time[0]-ConRef[0]);
      }
#     endif // #if ( MODEL != PAR_ONLY )
#     endif // #ifdef MASSIVE_PARTICLES


//    output
      File = fopen( FileName, "a" );

      Aux_Message( File, "%13.7e  %10ld", Time[0], Step );

      for (int v=0; v<NVar_Flu; v++)
      {

         Aux_Message( File, "  %17.7e  %17.7e  %17.7e", Fluid_AllRank[v], AbsErr_Flu[v], RelErr_Flu[v] );

         if ( v == index_before_column_CoM ) {
         for (int d=0; d<3; d++)
         Aux_Message( File, "  %17.7e  %17.7e  %17.7e", CoM_Flu[d], AbsErr_CoM_Flu[d], AveVel_CoM_Flu[d] );
         }
      } // for (int v=0; v<NVar_Flu; v++)

#     ifdef MASSIVE_PARTICLES
      for (int v=0; v<NVar_Par; v++)
      {
         Aux_Message( File, "  %17.7e  %17.7e  %17.7e", Par_AllRank[v], AbsErr_Par[v], RelErr_Par[v] );

         if ( v == index_before_column_CoM ) {
         for (int d=0; d<3; d++)
         Aux_Message( File, "  %17.7e  %17.7e  %17.7e", CoM_Par[d], AbsErr_CoM_Par[d], AveVel_CoM_Par[d] );
         }
      } // for (int v=0; v<NVar_Par; v++)

#     if ( MODEL != PAR_ONLY )
      for (int v=0; v<NVar_All; v++)
      {

         Aux_Message( File, "  %17.7e  %17.7e  %17.7e", All_AllRank[v], AbsErr_All[v], RelErr_All[v] );

         if ( v == index_before_column_CoM ) {
         for (int d=0; d<3; d++)
         Aux_Message( File, "  %17.7e  %17.7e  %17.7e", CoM_All[d], AbsErr_CoM_All[d], AveVel_CoM_All[d] );
         }
      } // for (int v=0; v<NVar_All; v++)
#     endif // if ( MODEL != PAR_ONLY )
#     endif // #ifdef MASSIVE_PARTICLES

      Aux_Message( File, "\n" );

      fclose( File );
   } // if ( MPI_Rank == 0 )


   if ( FirstTime )  FirstTime = false;


#  if ( MODEL == ELBDM )
   delete [] Flu_ELBDM;
#  endif


// calculate the ELBDM center-of-mass velocity for ELBDM_RemoveMotionCM()
#  if ( MODEL == ELBDM )
   if ( ELBDM_REMOVE_MOTION_CM != ELBDM_REMOVE_MOTION_CM_NONE )
   {
//    momentum --> velocity
      if ( MPI_Rank == 0 )
      {
         for (int d=0; d<3; d++)
            ELBDM_Vcm[d] = Fluid_AllRank[d+1]/Fluid_AllRank[0];
      }

//    broadcast
      MPI_Bcast( ELBDM_Vcm, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }
#  endif

} // FUNCTION : Aux_Check_Conservation
