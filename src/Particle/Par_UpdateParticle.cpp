#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE

#ifndef GRAVITY
#   error : ERROR : GRAVITY is not defined !!
#endif

#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_UpdateParticle
// Description :  Update particle position and velocity at the target level
//
// Note        :  1. Does not work with "P5_Gradient"
//                2. Does not take into account the "periodic B.C." when updating particles
//                   --> After update, the particles' position may lie outside the simulation box
//                   --> It will be corrected in the function "Par_PassParticle2Sibling"
//                   (however, periodic B.C. is taken into account when calculating the gravitational force)
//                3. For the K-D-K scheme, this function performs either prediction (K-D) or correction (last K) operation
//                   --> Use the input parameter "UpdateStep" to control
//                4. For the Euler scheme, this function completes the full update
//                   --> UpdateStep==PAR_UPSTEP_CORR is meaningless
//                5. For KDK, particle time is set to -0.5*dt after the K-D operation to indicate that they require
//                   velocity correction (the last K operation) later. Otherwise particles just cross from fine to coarse
//                   grids cannot be distinguished from those already in the coarse grid, while we only want to apply
//                   velocity correction to the former. After the velocity correction, particle time is set to TimeNew
//                   For Euler, particle time is set to TimeNew in the PAR_UPSTEP_PRED step
//                6. Particle time may not be synchronized (so different particles may have different time).
//                   --> For example, particles just cross from coarse (lv) to fine (lv+1) grids may have time greater than
//                       other particles at lv+1. Also, particles just cross from fine (lv) to coarse (lv-1) grids may have
//                       time less than other particles at lv-1.
//                   --> Only update particles with time < TimeNew
//
// Parameter   :  lv          : Target refinement level
//                TimeNew     : Targeted physical time to reach
//                TimeOld     : Physical time before update
//                UpdateStep  : PAR_UPSTEP_PRED/PAR_UPSTEP_CORR (CORR is for PAR_INTEG_KDK only)
//-------------------------------------------------------------------------------------------------------
void Par_UpdateParticle( const int lv, const double TimeNew, const double TimeOld, const ParUpStep_t UpdateStep )
{

   const ParInterp_t IntScheme  = amr->Par->Interp;
   const OptFluBC_t *FluBC_None = NULL;
   const bool IntPhase_No       = false;
   const bool GetTotDens_No     = false;
   const double dh              = amr->dh[lv];
   const double _dh             = 1.0/dh;
   const double PrepPotTime     = ( UpdateStep == PAR_UPSTEP_PRED ) ? TimeOld : TimeNew;

   const int  GraGhost_Par      = 1;                     // always set to 1 for particles (P5 gradient is not supported yet)
   const int  ParGhost          = amr->Par->GhostSize;
   const int  PotGhost          = GraGhost_Par + ParGhost;
   const int  PotSize           = PS1 + 2*PotGhost;
   const int  AccSize           = PS1 + 2*ParGhost;
   const real Const_8           = (real)8.0;

   real *Pos[3]  = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3]  = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   real *ParTime = amr->Par->Time;

   real *Pot = new real [ 8*CUBE(PotSize) ];             // 8: number of patches per patch group
   real *Acc = new real [ 3*CUBE(AccSize) ];             // 3: three dimension

   real (*Pot3D)[PotSize][PotSize][PotSize] = ( real (*)[PotSize][PotSize][PotSize] )Pot;
   real (*Acc3D)[AccSize][AccSize][AccSize] = ( real (*)[AccSize][AccSize][AccSize] )Acc;

   bool   GotYou;
   long   ParID;
   real   Acc_Temp[3];
   double PhyCorner_ExtAcc[3], PhyCorner_ExtPot[3], x, y, z, dt, dt_half;


// determine the coefficient for calculating acceleration
   real GraConst = ( OPT__GRA_P5_GRADIENT ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // but P5 is NOT supported yet


// determine PotSg for STORE_POT_GHOST
   int PotSg;

#  ifdef STORE_POT_GHOST
   if      (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][   amr->PotSg[lv] ], NULL, false )  )
      PotSg =   amr->PotSg[lv];

   else if (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][ 1-amr->PotSg[lv] ], NULL, false )  )
      PotSg = 1-amr->PotSg[lv];

   else
      Aux_Error( ERROR_INFO, "Cannot determine PotSg (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                 lv, PrepPotTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );
#  endif


// check
#  ifdef COMOVING
#  error : ERROR : does not support COMOVING because time-step has not been converted to comoving !!
#  endif

#  ifdef DEBUG_PARTICLE
// Par->ImproveAcc only works particle interpolation schemes with ParGhost == 1 (CIC & TSC)
   if ( amr->Par->ImproveAcc  &&  PotGhost != GRA_GHOST_SIZE )
      Aux_Error( ERROR_INFO, "PotGhost (%d) != GRA_GHOST_SIZE (%d) for amr->Par->ImproveAcc !!\n",
                 PotGhost, GRA_GHOST_SIZE );

// Par->ImproveAcc must work with STORE_POT_GHOST
#  ifndef STORE_POT_GHOST
   if ( amr->Par->ImproveAcc )
      Aux_Error( ERROR_INFO, "amr->Par->ImproveAcc must work with STORE_POT_GHOST !!\n" );
#  endif

   if ( amr->Par->Integ == PAR_INTEG_EULER  &&  UpdateStep == PAR_UPSTEP_CORR )
      Aux_Error( ERROR_INFO, "UpdateStep == PAR_UPSTEP_CORR is meaningless for PAR_INTEG_EULER !!\n" );

   if ( UpdateStep != PAR_UPSTEP_PRED  &&  UpdateStep != PAR_UPSTEP_CORR )
      Aux_Error( ERROR_INFO, "unsupported UpdateStep (%d) !!\n", UpdateStep );
#  endif // #ifdef DEBUG_PARTICLE


   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    1. find the patch groups with particles 
//    --> use patch group as the calculation unit since "Prepare_PatchData" only work with patch group
//    --> some patches may not have particles ...
      GotYou = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         if ( amr->patch[0][lv][PID]->NPar > 0 )
         {
            GotYou = true;
            break;
         }
      }


//    nothing to do if there are no particles in the target patch group
      if ( !GotYou )    continue;


//    2. prepare the potential data for the patch group with particles (need NSIDE_26 for ParGhost>0 )
      if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
      {
#        ifdef STORE_POT_GHOST
         if ( amr->Par->ImproveAcc  &&  lv > 0 )
         {
            for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
            {
               if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;   // skip patches with no particles

               memcpy( Pot3D[P], amr->patch[PotSg][lv][PID]->pot_ext, CUBE(PotSize)*sizeof(real) );
            }
         }

         else
#        endif
            Prepare_PatchData( lv, PrepPotTime, Pot, PotGhost, 1, &PID0, _POTE,
                               OPT__GRA_INT_SCHEME, UNIT_PATCH, NSIDE_26, IntPhase_No, FluBC_None, OPT__BC_POT, GetTotDens_No );
      }

      for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      {
         if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;   // skip patches with no particles

         for (int d=0; d<3; d++)    
         {
            PhyCorner_ExtAcc[d] = amr->patch[0][lv][PID]->EdgeL[d] + (0.5-ParGhost)*dh;
            PhyCorner_ExtPot[d] = amr->patch[0][lv][PID]->EdgeL[d] + (0.5-PotGhost)*dh;
         }

//       3. external potential (currently useful only for ELBDM)
         if ( OPT__EXTERNAL_POT )
         {
            for (int k=0; k<PotSize; k++)    {  z = PhyCorner_ExtPot[2] + (double)k*dh;
            for (int j=0; j<PotSize; j++)    {  y = PhyCorner_ExtPot[1] + (double)j*dh;
            for (int i=0; i<PotSize; i++)    {  x = PhyCorner_ExtPot[0] + (double)i*dh;

               Pot3D[P][k][j][i] += CPU_ExternalPot( x, y, z, PrepPotTime, ExtPot_AuxArray );

            }}}
         }

         for (int k=GraGhost_Par, kk=0; k<PotSize-GraGhost_Par; k++, kk++)  {  z = PhyCorner_ExtAcc[2] + (double)kk*dh;
         for (int j=GraGhost_Par, jj=0; j<PotSize-GraGhost_Par; j++, jj++)  {  y = PhyCorner_ExtAcc[1] + (double)jj*dh;
         for (int i=GraGhost_Par, ii=0; i<PotSize-GraGhost_Par; i++, ii++)  {  x = PhyCorner_ExtAcc[0] + (double)ii*dh;

            Acc_Temp[0] = (real)0.0;
            Acc_Temp[1] = (real)0.0;
            Acc_Temp[2] = (real)0.0;

//          4. external gravity (currently useful only for HYDRO)
            if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
               CPU_ExternalAcc( Acc_Temp, x, y, z, PrepPotTime, ExtAcc_AuxArray );

//          5. self-gravity
            if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
            {
//             OPT__GRA_P5_GRADIENT is not supported yet
//             if ( OPT__GRA_P5_GRADIENT )
               if ( false )
               {
                  Acc_Temp[0] += GraConst * ( -         Pot3D[P][k  ][j  ][i+2] +         Pot3D[P][k  ][j  ][i-2]
                                              + Const_8*Pot3D[P][k  ][j  ][i+1] - Const_8*Pot3D[P][k  ][j  ][i-1] );
                  Acc_Temp[1] += GraConst * ( -         Pot3D[P][k  ][j+2][i  ] +         Pot3D[P][k  ][j-2][i  ]
                                              + Const_8*Pot3D[P][k  ][j+1][i  ] - Const_8*Pot3D[P][k  ][j-1][i  ] );
                  Acc_Temp[2] += GraConst * ( -         Pot3D[P][k+2][j  ][i  ] +         Pot3D[P][k-2][j  ][i  ]
                                              + Const_8*Pot3D[P][k+1][j  ][i  ] - Const_8*Pot3D[P][k-1][j  ][i  ] );
               }

               else
               {
                  Acc_Temp[0] += GraConst * ( Pot3D[P][k  ][j  ][i+1] - Pot3D[P][k  ][j  ][i-1] );
                  Acc_Temp[1] += GraConst * ( Pot3D[P][k  ][j+1][i  ] - Pot3D[P][k  ][j-1][i  ] );
                  Acc_Temp[2] += GraConst * ( Pot3D[P][k+1][j  ][i  ] - Pot3D[P][k-1][j  ][i  ] );
               }
            }

            for (int d=0; d<3; d++)    Acc3D[d][kk][jj][ii] = Acc_Temp[d];

         }}} // i,j,k


//       6. update particles
         switch ( IntScheme )
         {
//          6.1 NGP
            case ( PAR_INTERP_NGP ):
            {
               int idx[3];

               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
//                6.1.1 calculate the nearest grid index
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  for (int d=0; d<3; d++)    
                  {
                     idx[d] = int( ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh );

//                   prevent from round-off errors (especially for NGP and TSC)
                     if ( idx[d] < 0 )
                     {
#                       ifdef DEBUG_PARTICLE
                        if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                                   d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idx[d] );
#                       endif

                        idx[d] = 0;
                     }

                     else if ( idx[d] >= AccSize )
                     {
#                       ifdef DEBUG_PARTICLE
                        if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                                      d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idx[d] );
#                       endif

                        idx[d] = AccSize - 1;
                     }
                  } // for (int d=0; d<3; d++)

//                6.1.2 update particles position and velocity
//                6.1.2-1: Euler method
                  if ( amr->Par->Integ == PAR_INTEG_EULER )
                  {
                     dt = TimeNew - ParTime[ParID];

                     for (int d=0; d<3; d++)
                     {
                        Pos[d][ParID] += Vel  [d][ParID]*dt;   // update position first for a full dt
                        Vel[d][ParID] += Acc3D[d][ idx[2] ][ idx[1] ][ idx[0] ]*dt;
                     }

                     ParTime[ParID] = TimeNew;
                  }

//                6.1.2-2: KDK scheme
                  else if ( amr->Par->Integ == PAR_INTEG_KDK )
                  {
//                   KDK prediction
                     if ( UpdateStep == PAR_UPSTEP_PRED )
                     {
                        dt      = TimeNew - ParTime[ParID];
                        dt_half = 0.5*dt;

                        for (int d=0; d<3; d++)
                        {
                           Vel[d][ParID] += Acc3D[d][ idx[2] ][ idx[1] ][ idx[0] ]*dt_half; // predict velocity for 0.5*dt
                           Pos[d][ParID] += Vel  [d][ParID]*dt;   // update position by the half-step velocity for a full dt
                        }

                        ParTime[ParID] = -dt_half;   // negative --> indicating that it requires velocity correction
                     }

//                   KDK correction for velocity
                     else // UpdateStep == PAR_UPSTEP_CORR
                     {
                        dt_half = -ParTime[ParID];

                        for (int d=0; d<3; d++)
                           Vel[d][ParID] += Acc3D[d][ idx[2] ][ idx[1] ][ idx[0] ]*dt_half; // correct velocity for 0.5*dt

                        ParTime[ParID] = TimeNew;
                     }
                  } // amr->Par->Integ
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // PAR_INTERP_NGP
            break;


//          6.2 CIC
            case ( PAR_INTERP_CIC ):
            {
               int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
               double dr      [3];     // distance to the center of the left cell
               double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  for (int d=0; d<3; d++)
                  {
//                   6.2.1 calculate the array index of the left and right cells
                     dr      [d] = ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost - 0.5;
                     idxLR[0][d] = int( dr[d] );
                     idxLR[1][d] = idxLR[0][d] + 1;

//                   prevent from round-off errors
//                   (CIC should be clear off this issue unless round-off erros are comparable to dh)
                     if ( idxLR[0][d] < 0 )  
                     {
#                       ifdef DEBUG_PARTICLE
                        if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                   d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLR[0][d], idxLR[1][d] );
#                       endif

                        idxLR[0][d] = 0;
                        idxLR[1][d] = 1;
                     }

                     else if ( idxLR[1][d] >= AccSize )
                     {
#                       ifdef DEBUG_PARTICLE
                        if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                   d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLR[0][d], idxLR[1][d] );
#                       endif

                        idxLR[0][d] = AccSize - 2;
                        idxLR[1][d] = AccSize - 1;
                     }

//                   6.2.2 get the weighting of the nearby 8 cells
                     dr     [d] -= (double)idxLR[0][d];
                     Frac[0][d]  = 1.0 - dr[d];
                     Frac[1][d]  =       dr[d];
                  } // for (int d=0; d<3; d++)

//                6.2.3 update particles position and velocity
//                6.2.3-1: Euler method
                  if ( amr->Par->Integ == PAR_INTEG_EULER )
                  {
                     dt = TimeNew - ParTime[ParID];

                     for (int d=0; d<3; d++)
                     {
                        Pos[d][ParID] += Vel[d][ParID]*dt;  // update position first for a full dt

                        for (int k=0; k<2; k++)
                        for (int j=0; j<2; j++)
                        for (int i=0; i<2; i++)
                        Vel[d][ParID] += Acc3D[d][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                        *Frac[i][0]*Frac[j][1]*Frac[k][2]*dt;
                     }

                     ParTime[ParID] = TimeNew;
                  }

//                6.1.2-2: KDK scheme
                  else if ( amr->Par->Integ == PAR_INTEG_KDK )
                  {
//                   KDK prediction
                     if ( UpdateStep == PAR_UPSTEP_PRED )
                     {
                        dt      = TimeNew - ParTime[ParID];
                        dt_half = 0.5*dt;

                        for (int d=0; d<3; d++)
                        {
                           for (int k=0; k<2; k++)
                           for (int j=0; j<2; j++)
                           for (int i=0; i<2; i++)
                           Vel[d][ParID] += Acc3D[d][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                           *Frac[i][0]*Frac[j][1]*Frac[k][2]*dt_half;  // predict velocity for 0.5*dt

                           Pos[d][ParID] += Vel[d][ParID]*dt;  // update position by the half-step velocity for a full dt
                        }

                        ParTime[ParID] = -dt_half;   // negative --> indicating that it requires velocity correction
                     }

//                   KDK correction for velocity
                     else // UpdateStep == PAR_UPSTEP_CORR
                     {
                        dt_half = -ParTime[ParID];

                        for (int d=0; d<3; d++)
                        {
                           for (int k=0; k<2; k++)
                           for (int j=0; j<2; j++)
                           for (int i=0; i<2; i++)
                           Vel[d][ParID] += Acc3D[d][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                           *Frac[i][0]*Frac[j][1]*Frac[k][2]*dt_half; // correct velocity for 0.5*dt
                        }

                        ParTime[ParID] = TimeNew;
                     }
                  } // amr->Par->Integ
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // PAR_INTERP_CIC
            break;


//          6.3 TSC
            case ( PAR_INTERP_TSC ):
            {
               int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
               double dr       [3];    // distance to the left edge of the central cell
               double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  for (int d=0; d<3; d++)
                  {
//                   6.3.1 calculate the array index of the left, central, and right cells
                     dr       [d] = ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost;
                     idxLCR[1][d] = int( dr[d] );
                     idxLCR[0][d] = idxLCR[1][d] - 1;
                     idxLCR[2][d] = idxLCR[1][d] + 1;

//                   prevent from round-off errors (especially for NGP and TSC)
                     if ( idxLCR[0][d] < 0 )
                     {
#                       ifdef DEBUG_PARTICLE
                        if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                   d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#                       endif

                        idxLCR[0][d] = 0;
                        idxLCR[1][d] = 1;
                        idxLCR[2][d] = 2;
                     }

                     else if ( idxLCR[2][d] >= AccSize )
                     {
#                       ifdef DEBUG_PARTICLE
                        if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                   d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#                       endif

                        idxLCR[0][d] = AccSize - 3;
                        idxLCR[1][d] = AccSize - 2;
                        idxLCR[2][d] = AccSize - 1;
                     }

//                   6.3.2 get the weighting of the nearby 27 cells
                     dr     [d] -= (double)idxLCR[1][d];
                     Frac[0][d]  = 0.5*SQR( 1.0 - dr[d] );
                     Frac[1][d]  = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
                     Frac[2][d]  = 0.5*SQR( dr[d] );
                  } // for (int d=0; d<3; d++)

//                6.3.3 update particles position and velocity
//                6.3.3-1: Euler method
                  if ( amr->Par->Integ == PAR_INTEG_EULER )
                  {
                     dt = TimeNew - ParTime[ParID];

                     for (int d=0; d<3; d++)
                     {
                        Pos[d][ParID] += Vel[d][ParID]*dt;  // update position first for a full dt 

                        for (int k=0; k<3; k++)
                        for (int j=0; j<3; j++)
                        for (int i=0; i<3; i++)
                        Vel[d][ParID] += Acc3D[d][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                        *Frac[i][0]*Frac[j][1]*Frac[k][2]*dt;
                     }

                     ParTime[ParID] = TimeNew;
                  }

//                6.3.3-2: KDK scheme
                  else if ( amr->Par->Integ == PAR_INTEG_KDK )
                  {
//                   KDK prediction
                     if ( UpdateStep == PAR_UPSTEP_PRED )
                     {
                        dt      = TimeNew - ParTime[ParID];
                        dt_half = 0.5*dt;

                        for (int d=0; d<3; d++)
                        {
                           for (int k=0; k<3; k++)
                           for (int j=0; j<3; j++)
                           for (int i=0; i<3; i++)
                           Vel[d][ParID] += Acc3D[d][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                           *Frac[i][0]*Frac[j][1]*Frac[k][2]*dt_half;  // predict velocity for 0.5*dt

                           Pos[d][ParID] += Vel[d][ParID]*dt;  // update position by the half-step velocity for a full dt
                        }

                        ParTime[ParID] = -dt_half;   // negative --> indicating that it requires velocity correction
                     }

//                   KDK correction for velocity
                     else // UpdateStep == PAR_UPSTEP_CORR
                     {
                        dt_half = -ParTime[ParID];

                        for (int d=0; d<3; d++)
                        {
                           for (int k=0; k<3; k++)
                           for (int j=0; j<3; j++)
                           for (int i=0; i<3; i++)
                           Vel[d][ParID] += Acc3D[d][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                           *Frac[i][0]*Frac[j][1]*Frac[k][2]*dt_half; // correct velocity for 0.5*dt
                        }

                        ParTime[ParID] = TimeNew;
                     }
                  } // amr->Par->Integ
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // PAR_INTERP_TSC
            break;

            default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
         } // switch ( IntScheme )
      } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


// 7. free memory
   delete [] Pot;
   delete [] Acc;

} // FUNCTION : Par_UpdateParticle



#endif // #ifdef PARTICLE
