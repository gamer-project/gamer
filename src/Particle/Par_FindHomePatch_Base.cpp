#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_FindHomePatch_Base
// Description :  Find the particle's home patch at the base level
//
// Note        :  1. This function uses BaseP[] to find the home patch
//                   --> It should NOT be used for LOAD_BALANCE except during the initialization
//                2. Invoked by Init_ByFunction()
//                3. This function assumes that all particles stored in "amr->Par->XXX" are located
//                   in the sub-domain calculated by this rank
//
// Parameter   :  BaseP : Array storing indices of all base-level patches
//-------------------------------------------------------------------------------------------------------
void Par_FindHomePatch_Base( const int *BaseP )
{

// check
   if ( BaseP == NULL )    Aux_Error( ERROR_INFO, "BaseP == NULL !!\n" );

#  ifndef SERIAL
   if ( amr->ParaVar == NULL )   Aux_Error( ERROR_INFO, "amr->ParaVar == NULL !!\n" );
#  endif


   const double PS0  = PS1*amr->dh[0];                               // base-level patch size
   const int  NP0[3] = { NX0[0]/PS1+4, NX0[1]/PS1+4, NX0[2]/PS1+4 }; // number of base-level patches including
                                                                     // 2 buffer patches on each sides

// 1D -> 3D array
   const int (*BaseP3D)[ NP0[1] ][ NP0[0] ] = ( const int (*)[ NP0[1] ][ NP0[0] ] )BaseP;

   real   *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   double dr[3], *EdgeL, *EdgeR;
   int    ijk[3], TPID, SibPID;
   bool   Pass;


// initialize the array for counting particles
   int *NParList = new int [ amr->NPatchComma[0][1] ];

   for (int PID=0; PID<amr->NPatchComma[0][1]; PID++)    NParList[PID] = 0;


// Stage 1: determine the number of particles at each patch and allocate the particle list
// Stage 2: store the particle indices into the particle list of each patch
   for (int Stage=0; Stage<2; Stage++)
   {
      for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
      {
//       0. skip inactive particles
         if ( amr->Par->Mass[ParID] < 0.0 )  continue;


//       1. evaluate the relative particle position and get the corresponding array indices in BaseP3D
         for (int d=0; d<3; d++)
         {
#           ifdef SERIAL
            dr [d] = Pos[d][ParID];
#           else
            dr [d] = Pos[d][ParID] - amr->ParaVar->SubDomain_EdgeL[d];
#           endif
            ijk[d] = 2 + int( dr[d]/PS0 );

#           ifdef DEBUG_PARTICLE
//          check: whether the array index is within the correct range
            if ( ijk[d] < 2  ||  ijk[d] > NP0[d]-3 )
               Aux_Error( ERROR_INFO, "particle (%ld) lies outside the target sub-domain (ijk[%d] = %d) !!\n",
                          ParID, d, ijk[d] );
#           endif // #ifdef DEBUG_PARTICLE

         } // for (int d=0; d<3; d++)


//       2. get the home patch ID
         TPID = BaseP3D[ ijk[2] ][ ijk[1] ][ ijk[0] ];

//       verify that the home patch exists
#        ifdef DEBUG_PARTICLE
         if ( TPID <= -1 )    Aux_Error( ERROR_INFO, "particle %ld, TPID (%d) <= -1 !!\n", ParID, TPID );
#        endif


//       3. check the home patch carefully to prevent from any issue resulting from round-off errors
         Pass = true;

         for (int d=0; d<3; d++)
         {
            if ( Pos[d][ParID] < amr->patch[0][0][TPID]->EdgeL[d]  ||  Pos[d][ParID] >= amr->patch[0][0][TPID]->EdgeR[d] )
            {
               Aux_Message( stderr, "WARNING : wrong home patch (PID %d, particle %ld)\n", TPID, ParID );
               for (int dd=0; dd<3; dd++)
               Aux_Message( stderr, "          L/R edge[%d] %13.7e/%13.7e, particle pos[%d] %13.7e\n",
                            dd, amr->patch[0][0][TPID]->EdgeL[dd], amr->patch[0][0][TPID]->EdgeR[dd], dd, Pos[dd][ParID] );
               Aux_Message( stderr, "          --> trying to find the correct home patch ...\n" );

               Pass     = false;
               break;
            }
         }

//       loop over sibling patches to try to find the correct home patch
         if ( !Pass )
         {
            for (int s=0; s<26; s++)
            {
               SibPID = amr->patch[0][0][TPID]->sibling[s];
               EdgeL  = amr->patch[0][0][SibPID]->EdgeL;
               EdgeR  = amr->patch[0][0][SibPID]->EdgeR;

               if ( SibPID < amr->NPatchComma[0][1]  &&  SibPID >= 0 )  // make sure that the sibling patch is a real patch
               if ( Pos[0][ParID] >= EdgeL[0]  &&  Pos[0][ParID] < EdgeR[0]  &&
                    Pos[1][ParID] >= EdgeL[1]  &&  Pos[1][ParID] < EdgeR[1]  &&
                    Pos[2][ParID] >= EdgeL[2]  &&  Pos[2][ParID] < EdgeR[2]    )
               {
                  Aux_Message( stderr, "WARNING : correct home patch has been found (PID %d)\n", SibPID );
                  for (int dd=0; dd<3; dd++)
                  Aux_Message( stderr, "          L/R edge[%d] %13.7e/%13.7e, particle pos[%d] %13.7e\n",
                               dd, amr->patch[0][0][SibPID]->EdgeL[dd], amr->patch[0][0][SibPID]->EdgeR[dd],
                               dd, Pos[dd][ParID] );

                  TPID = SibPID;
                  Pass = true;
                  break;
               }
            }
         } // if ( !Pass );

         if ( !Pass )   Aux_Error( ERROR_INFO, "home patch for particle %ld is lost !!\n", ParID );


#        ifdef DEBUG_PARTICLE
//       check: whether the particle indeed lies in the target home patch
         for (int d=0; d<3; d++)
         {
            if ( Pos[d][ParID] <  amr->patch[0][0][TPID]->EdgeL[d]  ||
                 Pos[d][ParID] >= amr->patch[0][0][TPID]->EdgeR[d]     )
               Aux_Error( ERROR_INFO, "wrong home patch (PID %d, L/R edge %13.7e/%13.7e, particle pos[%d][%d] %13.7e) !!\n",
                          TPID, amr->patch[0][0][TPID]->EdgeL[d], amr->patch[0][0][TPID]->EdgeR[d], d, ParID, Pos[d][ParID] );
         }

//       check: whether the home patch is a real patch
         if ( TPID >= amr->NPatchComma[0][1]  ||  TPID <= -1 )
            Aux_Error( ERROR_INFO, "particle %ld's home patch (PID %d) is a buffer patch !!\n", ParID, TPID );
#        endif // #ifdef DEBUG_PARTICLE


//       4-Stage0. count the accumulated number of particles
         if ( Stage == 0 )    NParList[TPID] ++;

//       4-Stage1. store the particle index into the particle list
         else
         {
//###NOTE : No OpenMP since AddParticle will modify amr->Par->NPar_Lv[]
#           ifdef DEBUG_PARTICLE
            const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
            amr->patch[0][0][TPID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[0], ParPos, amr->Par->NPar_AcPlusInac, __FUNCTION__ );
#           else
            amr->patch[0][0][TPID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[0] );
#           endif
         }
      } // for (long p=0; p<amr->Par->NPar; p++)


//    5-Stage0. allocate particle list for each patch (we set ParListSize==NPar during the initialization)
      if ( Stage == 0 )
      {
         for (int PID=0; PID<amr->NPatchComma[0][1]; PID++)
         {
            if ( NParList[PID] > 0 )
            {
               amr->patch[0][0][PID]->ParListSize = NParList[PID];
               amr->patch[0][0][PID]->ParList     = (long*)malloc( NParList[PID]*sizeof(long) );
            }
         }
      }
   } // for (int Stage=0; Stage<2; Stage++)


// free memory
   delete [] NParList;

} // FUNCTION : Par_FindHomePatch_Base



#endif // #ifdef PARTICLE
