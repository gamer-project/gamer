#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_GeneralGalaxy
// Description :  Create new star particles in the general methods where ther star formation is independent
//
// Note        :  1. One must turn on STORE_POT_GHOST when adopting STORE_PAR_ACC
//                   --> It is because, currently, this function always uses the pot_ext[] array of each patch
//                       to calculate the gravitationally acceleration of the new star particles
//                2. One must invoke Buf_GetBufferData( ..., _TOTAL, ... ) after calling this function
//
// Parameter   :  lv             : Target refinement level
//                TimeNew        : Current physical time (after advancing solution by dt)
//                dt             : Time interval to advance solution
//                                 --> Currently this function does not distinguish dt and the physical time interval (dTime)
//                                 --> Does NOT support COMOVING yet
//                RNG            : Random number generator
//                UseMetal       : Store the metal mass fraction in star particles
//
// Return      :  1. Particle repository will be updated
//                2. fluid[] array of gas will be updated
//-------------------------------------------------------------------------------------------------------
void SF_CreateStar_GeneralGalaxy( const int lv, const real TimeNew, const real dt, RandomNumber_t *RNG, const bool UseMetal )
{

// check
#  if ( defined STORE_PAR_ACC  &&  !defined STORE_POT_GHOST )
#     error : STAR_FORMATION + STORE_PAR_ACC must work with STORE_POT_GHOST !!
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "Must turn on GRAVITY for %s() !!\n", __FUNCTION__ );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "%s() does not support COMOVING yet !!\n", __FUNCTION__ );
#  endif

#  ifdef GAMER_DEBUG
   if ( UseMetal  &&  Idx_ParMetalFrac == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParMetalFrac is undefined for \"UseMetal\" !!\n" );

   if ( UseMetal  &&  Idx_Metal == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_Metal is undefined for \"UseMetal\" !!\n" );

   if ( Idx_ParCreTime == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParCreTime is undefined !!\n" );
#  endif // #ifdef GAMER_DEBUG


// constant parameters
   const double dh             = amr->dh[lv];
   const real   dv             = CUBE( dh );
   const int    FluSg          = amr->FluSg[lv];
   const int    PotSg          = amr->PotSg[lv];
// const real   GraConst       = ( OPT__GRA_P5_GRADIENT ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh);
   const real   GraConst       = ( false                ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // P5 is NOT supported yet


// start of OpenMP parallel region
#  pragma omp parallel
   {

// thread-private variables
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif

   double x0, y0, z0, x, y, z;
   real   GasDens, _GasDens, GasMass, StarMFrac, StarMass, GasMFracLeft;
   real   (*fluid)[PS1][PS1][PS1]      = NULL;
   real   (*MagCC)[PS1][PS1][PS1]      = NULL;
   real   (*Pres)[PS1][PS1]            = NULL;
   real   (*Cs2)[PS1][PS1]             = NULL;
#  ifdef STORE_POT_GHOST
   real   (*pot_ext)[GRA_NXT][GRA_NXT] = NULL;
#  endif

   bool NeedPres = false;
   bool NeedCs2  = false;

   if ( SF_CREATE_STAR_SCHEME == SF_CREATE_STAR_SCHEME_DWARFGALAXY )   NeedCs2 = true;
   if ( NeedCs2 )   NeedPres = true;

#  ifdef MHD
   if ( NeedPres )   MagCC = new real [3][PS1][PS1][PS1];
#  endif
   if ( NeedPres )   Pres  = new real    [PS1][PS1][PS1];
   if ( NeedCs2  )   Cs2   = new real    [PS1][PS1][PS1];

   const int    MaxNewParPerPatch = CUBE(PS1);
   real_par   (*NewParAttFlt)[PAR_NATT_FLT_TOTAL] = new real_par [MaxNewParPerPatch][PAR_NATT_FLT_TOTAL];
   long_par   (*NewParAttInt)[PAR_NATT_INT_TOTAL] = new long_par [MaxNewParPerPatch][PAR_NATT_INT_TOTAL];
   long        *NewParID                          = new long     [MaxNewParPerPatch];

   int NNewPar;


// loop over all real patches
// use static schedule to ensure bitwise reproducibility when running with the same numbers of OpenMP threads and MPI ranks
// --> bitwise reproducibility will still break when running with different numbers of OpenMP threads and/or MPI ranks
//     unless both BITWISE_REPRODUCIBILITY and SF_CREATE_STAR_DET_RANDOM are enabled
#  pragma omp for schedule( static )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;


//    to get deterministic and different random numbers for all patches, reset the random seed of each patch according to
//    its location and time
//    --> patches at different time and/or AMR levels may still have the same random seeds...
      if ( SF_CREATE_STAR_DET_RANDOM )
      {
//       the factor "1.0e6" in the end is just to make random seeds at different times more different, especially for
//       extremely small time-step
         const long RSeed = SF_CREATE_STAR_RSEED + amr->patch[0][lv][PID]->LB_Idx + long(TimeNew*UNIT_T/Const_yr*1.0e6);
         RNG->SetSeed( TID, RSeed );
      }


      fluid   = amr->patch[FluSg][lv][PID]->fluid;
#     ifdef STORE_POT_GHOST
      pot_ext = amr->patch[PotSg][lv][PID]->pot_ext;
#     endif

#     ifdef MHD
//    evaluate cell-centered B field
      if ( NeedPres )
      {
         real MagCC_1Cell[NCOMP_MAG];

         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            MHD_GetCellCenteredBFieldInPatch( MagCC_1Cell, lv, PID, i, j, k, amr->MagSg[lv] );

            for (int v=0; v<NCOMP_MAG; v++)  MagCC[v][k][j][i] = MagCC_1Cell[v];
         }
      } // if ( NeedPres )
#     endif // #ifdef MHD

//    evaluate pressure
      if ( NeedPres )
      {
         const bool CheckMinPres_Yes = true;

         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
//          if applicable, compute pressure from the dual-energy variable to reduce the round-off errors
#           ifdef DUAL_ENERGY

#           if   ( DUAL_ENERGY == DE_ENPY )
            Pres[k][j][i] = Hydro_DensDual2Pres( fluid[DENS][k][j][i], fluid[DUAL][k][j][i],
                                                 EoS_AuxArray_Flt[1], CheckMinPres_Yes, MIN_PRES );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

#           else // #ifdef DUAL_ENERGY

#           ifdef MHD
            const real Emag = (real)0.5*(  SQR( MagCC[MAGX][k][j][i] )
                                         + SQR( MagCC[MAGY][k][j][i] )
                                         + SQR( MagCC[MAGZ][k][j][i] )  );
#           else
            const real Emag = NULL_REAL;
#           endif

#           if ( EOS != EOS_GAMMA  &&  EOS != EOS_ISOTHERMAL  &&  NCOMP_PASSIVE > 0 )
            real Passive[NCOMP_PASSIVE];
            for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = fluid[ NCOMP_FLUID + v ][k][j][i];
#           else
            const real *Passive = NULL;
#           endif

            Pres[k][j][i] = Hydro_Con2Pres( fluid[DENS][k][j][i], fluid[MOMX][k][j][i], fluid[MOMY][k][j][i],
                                            fluid[MOMZ][k][j][i], fluid[ENGY][k][j][i], Passive,
                                            CheckMinPres_Yes, MIN_PRES, PassiveFloorMask, Emag,
                                            EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                            NULL );
#           endif // #ifdef DUAL_ENERGY ... else ...
         } // k,j,i
      } // if ( NeedPres )

//    evaluate sound speed squared
      if ( NeedCs2 )
      {
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
#           if ( EOS != EOS_GAMMA  &&  EOS != EOS_ISOTHERMAL  &&  NCOMP_PASSIVE > 0 )
            real Passive[NCOMP_PASSIVE];
            for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = fluid[ NCOMP_FLUID + v ][k][j][i];
#           else
            const real *Passive = NULL;
#           endif

            Cs2[k][j][i] = EoS_DensPres2CSqr_CPUPtr( fluid[DENS][k][j][i], Pres[k][j][i], Passive,
                                                     EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         } // k,j,i
      } // if ( NeedCs2 )

      x0      = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0      = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0      = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;
      NNewPar = 0;

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {

         GasDens = fluid[DENS][k][j][i];
         GasMass = GasDens*dv;

//       1. star formation criteria and star mass based on various star formation prescriptions
//       ===========================================================================================================

//       1-1. check star formation criteria
         if ( !SF_CreateStar_Check( lv, PID, i, j, k, dh, fluid, Pres, Cs2 ) )   continue;

//       1-2. get the star mass
         StarMass = SF_CreateStar_GetStarMass( GasDens, dv, dt, RNG, TID );
         if ( StarMass <= 0.0 )   continue;

//       check the maximum gas mass fraction allowed to convert to stars
         StarMFrac = MIN( StarMass/GasMass, SF_CREATE_STAR_MAX_STAR_MFRAC );
         StarMass  = GasMass*StarMFrac;



//       2. store the information of new star particles
//       --> we will not create these new particles until looping over all cells in a patch in order to reduce
//           the OpenMP synchronization overhead
//       ===========================================================================================================
//       check
#        ifdef GAMER_DEBUG
         if ( NNewPar >= MaxNewParPerPatch )
            Aux_Error( ERROR_INFO, "NNewPar (%d) >= MaxNewParPerPatch (%d) !!\n", NNewPar, MaxNewParPerPatch );
#        endif

//       2-1. intrinsic attributes
         _GasDens = (real)1.0 / GasDens;
         x        = x0 + i*dh;
         y        = y0 + j*dh;
         z        = z0 + k*dh;

         NewParAttFlt[NNewPar][PAR_MASS] = StarMass;
         NewParAttFlt[NNewPar][PAR_POSX] = x;
         NewParAttFlt[NNewPar][PAR_POSY] = y;
         NewParAttFlt[NNewPar][PAR_POSZ] = z;
         NewParAttFlt[NNewPar][PAR_VELX] = fluid[MOMX][k][j][i]*_GasDens;
         NewParAttFlt[NNewPar][PAR_VELY] = fluid[MOMY][k][j][i]*_GasDens;
         NewParAttFlt[NNewPar][PAR_VELZ] = fluid[MOMZ][k][j][i]*_GasDens;
         NewParAttFlt[NNewPar][PAR_TIME] = TimeNew;
         NewParAttInt[NNewPar][PAR_TYPE] = PTYPE_STAR;

//       particle acceleration
#        ifdef STORE_PAR_ACC
         real GasAcc[3] = { (real)0.0, (real)0.0, (real)0.0 };

//       external acceleration
         if ( OPT__EXT_ACC )  CPUExtAcc_Ptr( GasAcc, x, y, z, TimeNew, ExtAcc_AuxArray );

//       self-gravity and external potential
         if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
         {
            const int ii = i + GRA_GHOST_SIZE;
            const int jj = j + GRA_GHOST_SIZE;
            const int kk = k + GRA_GHOST_SIZE;

#           ifdef STORE_POT_GHOST
            const real pot_xm = pot_ext[kk  ][jj  ][ii-1];
            const real pot_xp = pot_ext[kk  ][jj  ][ii+1];
            const real pot_ym = pot_ext[kk  ][jj-1][ii  ];
            const real pot_yp = pot_ext[kk  ][jj+1][ii  ];
            const real pot_zm = pot_ext[kk-1][jj  ][ii  ];
            const real pot_zp = pot_ext[kk+1][jj  ][ii  ];
#           endif

            GasAcc[0] += GraConst*( pot_xp - pot_xm );
            GasAcc[1] += GraConst*( pot_yp - pot_ym );
            GasAcc[2] += GraConst*( pot_zp - pot_zm );
         }

         NewParAttFlt[NNewPar][PAR_ACCX] = GasAcc[0];
         NewParAttFlt[NNewPar][PAR_ACCY] = GasAcc[1];
         NewParAttFlt[NNewPar][PAR_ACCZ] = GasAcc[2];
#        endif // ifdef STORE_PAR_ACC


//       2-2. extrinsic attributes
//       note that we store the metal mass **fraction** instead of density in particles
         if ( UseMetal )
         NewParAttFlt[NNewPar][Idx_ParMetalFrac] = fluid[Idx_Metal][k][j][i] * _GasDens;

         NewParAttFlt[NNewPar][Idx_ParCreTime  ] = TimeNew;

         NNewPar ++;



//       3. remove the gas that has been converted to stars
//       ===========================================================================================================
         GasMFracLeft = (real)1.0 - StarMFrac;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v][k][j][i] *= GasMFracLeft;
      } // i,j,k



//    4. create new star particles
//    ===========================================================================================================
//    use OpenMP critical construct since both amr->Par->AddOneParticle() and amr->patch[0][lv][PID]->AddParticle()
//    will modify some global variables
//    --> note that the order of which thread calls amr->Par->AddOneParticle() is non-deterministic and may change from run to run
//        --> order of particles stored in the particle repository (i.e., their particle ID) may change from run to run
//        --> particle text file may change from run to run since it's dumped according to the order of particle ID
//    --> but it's not an issue since the actual data of each particle will not be affected
#     pragma omp critical
      {
//       4-1. add particles to the particle repository
         for (int p=0; p<NNewPar; p++)
            NewParID[p] = amr->Par->AddOneParticle( NewParAttFlt[p], NewParAttInt[p] );


//       4-2. add particles to the patch
         const long_par *PType = amr->Par->Type;
#        ifdef DEBUG_PARTICLE
//       do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//       may change after calling amr->Par->AddOneParticle()
         const real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         char Comment[100];
         sprintf( Comment, "%s", __FUNCTION__ );

         amr->patch[0][lv][PID]->AddParticle( NNewPar, NewParID, &amr->Par->NPar_Lv[lv],
                                              PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#        else
         amr->patch[0][lv][PID]->AddParticle( NNewPar, NewParID, &amr->Par->NPar_Lv[lv], PType );
#        endif
      } // pragma omp critical

   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// free memory
#  ifdef MHD
   if ( NeedPres )   delete [] MagCC;
#  endif
   if ( NeedPres )   delete [] Pres;
   if ( NeedCs2  )   delete [] Cs2;
   delete [] NewParAttFlt;
   delete [] NewParAttInt;
   delete [] NewParID;

   } // end of OpenMP parallel region


// get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : SF_CreateStar_GeneralGalaxy



#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )
