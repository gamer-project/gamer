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
//                   --> Note that during RESTART the reference values will be recalculated since they are NOT
//                       recorded in the output files currently
//                       --> Error estimation will be incorrect ...
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
   const char *FileName  = "Record__Conservation";


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
#  else
   const int    NVar_NoPassive    = 11;   // 11: mass, momentum (x/y/z), angular momentum (x/y/z), kinetic/internal/potential/total energies
                                          // --> note that **total energy** is put in the last element
#  endif
   const int    idx_etot          = NVar_NoPassive - 1;
   const bool   CheckMinEint_No   = false;

#  elif ( MODEL == ELBDM )
   const int    NVar_NoPassive    = 11;   // 11: mass, momentum (x/y/z), angular momentum (x/y/z), kinetic/gravitational/self-interaction/total energies
   const int    idx_etot          = NVar_NoPassive - 1;
   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const int    NGhost            = 1;    // number of ghost zones for calculating the gradient of wave function
   const int    Size_Flu          = PS1 + 2*NGhost;
   const int    NPG               = 1;
   const double _Eta              = 1.0/ELBDM_ETA;
   const double _2Eta2            = 0.5/SQR(ELBDM_ETA);
   const IntScheme_t IntScheme    = INT_CQUAR;

   real (*Flu_ELBDM)[2][Size_Flu][Size_Flu][Size_Flu] = new real [NPG*8][2][Size_Flu][Size_Flu][Size_Flu];

#  else
#  error : ERROR : unsupported MODEL !!
#  endif


// get the sum of passive scalars to be normalized
   const bool GetPassiveSum = ( PassiveNorm_NVar > 0 );
   const int  NVar_Max      = NVar_NoPassive + NCOMP_PASSIVE + 1; // for declaring the static variable Fluid_Ref
   const int  NVar          = NVar_NoPassive + NCOMP_PASSIVE + ( (GetPassiveSum)?1:0 );

   double dh, dv, Fluid_ThisRank[NVar], Fluid_AllRank[NVar], Fluid_lv[NVar];   // dv : cell volume at each level
   int    FluSg;
#  ifdef GRAVITY
   int    PotSg;
#  endif
#  ifdef MHD
   int    MagSg;
#  endif
   FILE  *File = NULL;


// initialize accumulative variables as zero
   for (int v=0; v<NVar; v++)    Fluid_ThisRank[v] = 0.0;


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int v=0; v<NVar; v++)    Fluid_lv[v] = 0.0;

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
      Fluid_lv[idx_etot] = Fluid_lv[7] + Fluid_lv[8] + Fluid_lv[9] + Fluid_lv[10];
#     else
      Fluid_lv[idx_etot] = Fluid_lv[7] + Fluid_lv[8] + Fluid_lv[9];
#     endif
#     elif ( MODEL == ELBDM )
      Fluid_lv[idx_etot] = Fluid_lv[7] + Fluid_lv[8] + Fluid_lv[9];
#     else
#     error : ERROR : unsupported MODEL !!
#     endif

//    sum of passive scalars to be normalized
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<PassiveNorm_NVar; v++)
         Fluid_lv[ NVar - 1 ] += Fluid_lv[ NVar_NoPassive + PassiveNorm_VarIdx[v] ];
#     endif

//    multiply by the cell volume and sum over all levels
      for (int v=0; v<NVar; v++)    Fluid_ThisRank[v] += Fluid_lv[v]*dv;
   } // for (int lv=0; lv<NLEVEL; lv++)


// sum over all ranks
   MPI_Reduce( Fluid_ThisRank, Fluid_AllRank, NVar, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// compute the center of mass
   double CoM_Gas[3];
   double FinaldR_Gas;
   int    FinalNIter_Gas;
   Aux_FindWeightedAverageCenter( CoM_Gas, amr->BoxCenter, __FLT_MAX__, 0.0, _DENS, __FLT_MAX__, 1, &FinaldR_Gas, &FinalNIter_Gas );


// calculate conserved quantities for particles
#  ifdef MASSIVE_PARTICLES
   double Mass_Par, CoMX_Par, CoMY_Par, CoMZ_Par;
   double MomX_Par, MomY_Par, MomZ_Par;
   double AngMomX_Par, AngMomY_Par, AngMomZ_Par, Ekin_Par, Epot_Par, Etot_Par;

   Par_Aux_GetConservedQuantity( Mass_Par, CoMX_Par, CoMY_Par, CoMZ_Par,
                                 MomX_Par, MomY_Par, MomZ_Par,
                                 AngMomX_Par, AngMomY_Par, AngMomZ_Par, Ekin_Par, Epot_Par );

   Etot_Par = Ekin_Par + Epot_Par;
#  endif


// output
   if ( MPI_Rank == 0 )
   {
//    calculate the sum of conserved quantities in different models
#     if ( defined MASSIVE_PARTICLES  &&  MODEL != PAR_ONLY )
      const double Mass_All    = Fluid_AllRank[       0] + Mass_Par;
      const double MomX_All    = Fluid_AllRank[       1] + MomX_Par;
      const double MomY_All    = Fluid_AllRank[       2] + MomY_Par;
      const double MomZ_All    = Fluid_AllRank[       3] + MomZ_Par;
      const double AngMomX_All = Fluid_AllRank[       4] + AngMomX_Par;
      const double AngMomY_All = Fluid_AllRank[       5] + AngMomY_Par;
      const double AngMomZ_All = Fluid_AllRank[       6] + AngMomZ_Par;
      const double Etot_All    = Fluid_AllRank[idx_etot] + Etot_Par;    // for HYDRO/ELBDM, total energy is stored in the last element

      const double CoMX_All    = ( Fluid_AllRank[0]*CoM_Gas[0] + Mass_Par*CoMX_Par )/Mass_All;
      const double CoMY_All    = ( Fluid_AllRank[0]*CoM_Gas[1] + Mass_Par*CoMY_Par )/Mass_All;
      const double CoMZ_All    = ( Fluid_AllRank[0]*CoM_Gas[2] + Mass_Par*CoMZ_Par )/Mass_All;
#     endif // if ( defined MASSIVE_PARTICLES  &&  MODEL != PAR_ONLY )

//    note that a variable length array cannot have static storage duration
      static double Time_Ref;
      static double Fluid_Ref[NVar_Max];
      static double CoM_Gas_Ref[3];
#     ifdef MASSIVE_PARTICLES
      static double Mass_Par_Ref, CoMX_Par_Ref, CoMY_Par_Ref, CoMZ_Par_Ref;
      static double MomX_Par_Ref, MomY_Par_Ref, MomZ_Par_Ref;
      static double AngMomX_Par_Ref, AngMomY_Par_Ref, AngMomZ_Par_Ref, Ekin_Par_Ref, Epot_Par_Ref, Etot_Par_Ref;
#     if ( MODEL != PAR_ONLY )
      static double Mass_All_Ref, CoMX_All_Ref, CoMY_All_Ref, CoMZ_All_Ref;
      static double MomX_All_Ref, MomY_All_Ref, MomZ_All_Ref;
      static double AngMomX_All_Ref, AngMomY_All_Ref, AngMomZ_All_Ref, Etot_All_Ref;
#     endif
#     endif // #ifdef PARTICLE
      double AbsErr[NVar], RelErr[NVar];

      if ( FirstTime )
      {
//       record the reference values
         Time_Ref = Time[0];
         for (int v=0; v<NVar; v++)    Fluid_Ref[v]   = Fluid_AllRank[v];
         for (int d=0; d<3; d++)       CoM_Gas_Ref[d] = CoM_Gas[d];

#        ifdef MASSIVE_PARTICLES
         Mass_Par_Ref    =    Mass_Par;
         CoMX_Par_Ref    =    CoMX_Par;
         CoMY_Par_Ref    =    CoMY_Par;
         CoMZ_Par_Ref    =    CoMZ_Par;
         MomX_Par_Ref    =    MomX_Par;
         MomY_Par_Ref    =    MomY_Par;
         MomZ_Par_Ref    =    MomZ_Par;
         AngMomX_Par_Ref = AngMomX_Par;
         AngMomY_Par_Ref = AngMomY_Par;
         AngMomZ_Par_Ref = AngMomZ_Par;
         Ekin_Par_Ref    =    Ekin_Par;
         Epot_Par_Ref    =    Epot_Par;
         Etot_Par_Ref    =    Etot_Par;

#        if ( MODEL != PAR_ONLY )
         Mass_All_Ref    =    Mass_All;
         CoMX_All_Ref    =    CoMX_All;
         CoMY_All_Ref    =    CoMY_All;
         CoMZ_All_Ref    =    CoMZ_All;
         MomX_All_Ref    =    MomX_All;
         MomY_All_Ref    =    MomY_All;
         MomZ_All_Ref    =    MomZ_All;
         AngMomX_All_Ref = AngMomX_All;
         AngMomY_All_Ref = AngMomY_All;
         AngMomZ_All_Ref = AngMomZ_All;
         Etot_All_Ref    =    Etot_All;
#        endif // #if ( MODEL != PAR_ONLY )

#        endif // #ifdef PARTICLE


//       output header
         FILE *File = fopen( FileName, "a" );

         Aux_Message( File, "# Ref time        : %13.7e\n", Time_Ref );
         Aux_Message( File, "# Ref step        : %ld\n",    Step     );
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

#        if   ( MODEL == HYDRO )
         Aux_Message( File, "  %17s  %17s  %17s",    "Mass_Gas",    "Mass_Gas_AErr",    "Mass_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMX_Gas",    "CoMX_Gas_AErr",    "CoMX_Gas_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMY_Gas",    "CoMY_Gas_AErr",    "CoMY_Gas_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMZ_Gas",    "CoMZ_Gas_AErr",    "CoMZ_Gas_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomX_Gas",    "MomX_Gas_AErr",    "MomX_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomY_Gas",    "MomY_Gas_AErr",    "MomY_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomZ_Gas",    "MomZ_Gas_AErr",    "MomZ_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomX_Gas", "AngMomX_Gas_AErr", "AngMomX_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomY_Gas", "AngMomY_Gas_AErr", "AngMomY_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomZ_Gas", "AngMomZ_Gas_AErr", "AngMomZ_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Ekin_Gas",    "Ekin_Gas_AErr",    "Ekin_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Eint_Gas",    "Eint_Gas_AErr",    "Eint_Gas_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Epot_Gas",    "Epot_Gas_AErr",    "Epot_Gas_RErr" );
#        ifdef MHD
         Aux_Message( File, "  %17s  %17s  %17s",    "Emag_Gas",    "Emag_Gas_AErr",    "Emag_Gas_RErr" );
#        endif
         Aux_Message( File, "  %17s  %17s  %17s",    "Etot_Gas",    "Etot_Gas_AErr",    "Etot_Gas_RErr" );

#        elif ( MODEL == ELBDM )
         Aux_Message( File, "  %17s  %17s  %17s",    "Mass_Psi",    "Mass_Psi_AErr",    "Mass_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMX_Psi",    "CoMX_Psi_AErr",    "CoMX_Psi_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMY_Psi",    "CoMY_Psi_AErr",    "CoMY_Psi_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMZ_Psi",    "CoMZ_Psi_AErr",    "CoMZ_Psi_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomX_Psi",    "MomX_Psi_AErr",    "MomX_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomY_Psi",    "MomY_Psi_AErr",    "MomY_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomZ_Psi",    "MomZ_Psi_AErr",    "MomZ_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomX_Psi", "AngMomX_Psi_AErr", "AngMomX_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomY_Psi", "AngMomY_Psi_AErr", "AngMomY_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomZ_Psi", "AngMomZ_Psi_AErr", "AngMomZ_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Ekin_Psi",    "Ekin_Psi_AErr",    "Ekin_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Epot_Psi",    "Epot_Psi_AErr",    "Epot_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Esel_Psi",    "Esel_Psi_AErr",    "Esel_Psi_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Etot_Psi",    "Etot_Psi_AErr",    "Etot_Psi_RErr" );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         Aux_Message( File, "  %17s  %12s_AErr  %12s_RErr", FieldLabel[v], FieldLabel[v], FieldLabel[v] );

         if ( GetPassiveSum )
         Aux_Message( File, "  %17s  %17s  %17s",    "PassNorm",    "PassNorm_AErr",    "PassNorm_RErr" );

#        ifdef MASSIVE_PARTICLES
         Aux_Message( File, "  %17s  %17s  %17s",    "Mass_Par",    "Mass_Par_AErr",    "Mass_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMX_Par",    "CoMX_Par_AErr",    "CoMX_Par_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMY_Par",    "CoMY_Par_AErr",    "CoMY_Par_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMZ_Par",    "CoMZ_Par_AErr",    "CoMZ_Par_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomX_Par",    "MomX_Par_AErr",    "MomX_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomY_Par",    "MomY_Par_AErr",    "MomY_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomZ_Par",    "MomZ_Par_AErr",    "MomZ_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomX_Par", "AngMomX_Par_AErr", "AngMomX_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomY_Par", "AngMomY_Par_AErr", "AngMomY_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomZ_Par", "AngMomZ_Par_AErr", "AngMomZ_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Ekin_Par",    "Ekin_Par_AErr",    "Ekin_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Epot_Par",    "Epot_Par_AErr",    "Epot_Par_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Etot_Par",    "Etot_Par_AErr",    "Etot_Par_RErr" );

#        if ( MODEL != PAR_ONLY )
         Aux_Message( File, "  %17s  %17s  %17s",    "Mass_All",    "Mass_All_AErr",    "Mass_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMX_All",    "CoMX_All_AErr",    "CoMX_All_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMY_All",    "CoMY_All_AErr",    "CoMY_All_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "CoMZ_All",    "CoMZ_All_AErr",    "CoMZ_All_AveV" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomX_All",    "MomX_All_AErr",    "MomX_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomY_All",    "MomY_All_AErr",    "MomY_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "MomZ_All",    "MomZ_All_AErr",    "MomZ_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomX_All", "AngMomX_All_AErr", "AngMomX_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomY_All", "AngMomY_All_AErr", "AngMomY_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s", "AngMomZ_All", "AngMomZ_All_AErr", "AngMomZ_All_RErr" );
         Aux_Message( File, "  %17s  %17s  %17s",    "Etot_All",    "Etot_All_AErr",    "Etot_All_RErr" );
#        endif // if ( MODEL != PAR_ONLY )
#        endif // #ifdef PARTICLE

         Aux_Message( File, "\n" );

         fclose( File );
      } // if ( FirstTime )


//    calculate errors
      for (int v=0; v<NVar; v++)
      {
         AbsErr[v] = Fluid_AllRank[v] - Fluid_Ref[v];
         RelErr[v] = AbsErr[v] / fabs(Fluid_Ref[v]);
      }


//    output
      File = fopen( FileName, "a" );

      Aux_Message( File, "%13.7e  %10ld", Time[0], Step );

      const int index_before_column_CoM = 0;

      for (int v=0; v<NVar; v++)
      {

      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", Fluid_AllRank[v], AbsErr[v], RelErr[v] );

      if ( v == index_before_column_CoM )
      {
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",  CoM_Gas[0],   CoM_Gas[0]-CoM_Gas_Ref[0],      (CoM_Gas[0]-CoM_Gas_Ref[0])/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",  CoM_Gas[1],   CoM_Gas[1]-CoM_Gas_Ref[1],      (CoM_Gas[1]-CoM_Gas_Ref[1])/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",  CoM_Gas[2],   CoM_Gas[2]-CoM_Gas_Ref[2],      (CoM_Gas[2]-CoM_Gas_Ref[2])/(Time[0]-Time_Ref) );
      }

      }

#     ifdef MASSIVE_PARTICLES
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    Mass_Par,       Mass_Par-Mass_Par_Ref,          (Mass_Par-Mass_Par_Ref)/fabs(Mass_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    CoMX_Par,       CoMX_Par-CoMX_Par_Ref,          (CoMX_Par-CoMX_Par_Ref)/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    CoMY_Par,       CoMY_Par-CoMY_Par_Ref,          (CoMY_Par-CoMY_Par_Ref)/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    CoMZ_Par,       CoMZ_Par-CoMZ_Par_Ref,          (CoMZ_Par-CoMZ_Par_Ref)/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    MomX_Par,       MomX_Par-MomX_Par_Ref,          (MomX_Par-MomX_Par_Ref)/fabs(MomX_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    MomY_Par,       MomY_Par-MomY_Par_Ref,          (MomY_Par-MomY_Par_Ref)/fabs(MomY_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    MomZ_Par,       MomZ_Par-MomZ_Par_Ref,          (MomZ_Par-MomZ_Par_Ref)/fabs(MomZ_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", AngMomX_Par, AngMomX_Par-AngMomX_Par_Ref, (AngMomX_Par-AngMomX_Par_Ref)/fabs(AngMomX_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", AngMomY_Par, AngMomY_Par-AngMomY_Par_Ref, (AngMomY_Par-AngMomY_Par_Ref)/fabs(AngMomY_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", AngMomZ_Par, AngMomZ_Par-AngMomZ_Par_Ref, (AngMomZ_Par-AngMomZ_Par_Ref)/fabs(AngMomZ_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    Ekin_Par,       Ekin_Par-Ekin_Par_Ref,          (Ekin_Par-Ekin_Par_Ref)/fabs(Ekin_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    Epot_Par,       Epot_Par-Epot_Par_Ref,          (Epot_Par-Epot_Par_Ref)/fabs(Epot_Par_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    Etot_Par,       Etot_Par-Etot_Par_Ref,          (Etot_Par-Etot_Par_Ref)/fabs(Etot_Par_Ref) );

#     if ( MODEL != PAR_ONLY )
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    Mass_All,       Mass_All-Mass_All_Ref,          (Mass_All-Mass_All_Ref)/fabs(Mass_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    CoMX_All,       CoMX_All-CoMX_All_Ref,          (CoMX_All-CoMX_All_Ref)/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    CoMY_All,       CoMY_All-CoMY_All_Ref,          (CoMY_All-CoMY_All_Ref)/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    CoMZ_All,       CoMZ_All-CoMZ_All_Ref,          (CoMZ_All-CoMZ_All_Ref)/(Time[0]-Time_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    MomX_All,       MomX_All-MomX_All_Ref,          (MomX_All-MomX_All_Ref)/fabs(MomX_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    MomY_All,       MomY_All-MomY_All_Ref,          (MomY_All-MomY_All_Ref)/fabs(MomY_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    MomZ_All,       MomZ_All-MomZ_All_Ref,          (MomZ_All-MomZ_All_Ref)/fabs(MomZ_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", AngMomX_All, AngMomX_All-AngMomX_All_Ref, (AngMomX_All-AngMomX_All_Ref)/fabs(AngMomX_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", AngMomY_All, AngMomY_All-AngMomY_All_Ref, (AngMomY_All-AngMomY_All_Ref)/fabs(AngMomY_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e", AngMomZ_All, AngMomZ_All-AngMomZ_All_Ref, (AngMomZ_All-AngMomZ_All_Ref)/fabs(AngMomZ_All_Ref) );
      Aux_Message( File, "  %17.7e  %17.7e  %17.7e",    Etot_All,       Etot_All-Etot_All_Ref,          (Etot_All-Etot_All_Ref)/fabs(Etot_All_Ref) );
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
