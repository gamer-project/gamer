#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Aux_Check_Particle
// Description :  Verify the following.
//                1.   Particles reside in their home patches
//                2.   Particles always reside in the leaf patches
//                3.   There are no missing or redundant particles
//                4.   No active particles have mass=-1.0
//                5/6. Each particle has one and only one home patch
//                7.   NPar_AcPlusInac == NPar_Active + NPar_Inactive
//                8.   NPar_Active_AllRank = sum(NPar_Active, All ranks)
//                9.   NPar_Active = sum(NPar_Lv, all levels)
//                10.  All patches have NPar_Copy == -1
//                11.  Particle types are recognizable
//
// Note        :  None
//
// Parameter   :  comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Par_Aux_Check_Particle( const char *comment )
{

   const int       NCheck    = 11;
   const real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };

   int     PassAll    = true;
   long    NParInLeaf = 0;
   bool   *ParHome    = new bool [amr->Par->NPar_AcPlusInac];  // true/false --> particle has home/is homeless
   int     NParThisPatch;
   long    ParID, NPar_Active_AllRank_Expect;
   double *EdgeL, *EdgeR;
   int     PassCheck[NCheck];


// initialize the check list
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   ParHome[p] = false;

   for (int t=0; t<NCheck; t++)  PassCheck[t] = true;


// get the total number of active particles
   MPI_Allreduce( &amr->Par->NPar_Active, &NPar_Active_AllRank_Expect, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );


// start checking
   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
//       loop over all "real" patches
         for (int lv=0; lv<NLEVEL; lv++)
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            NParThisPatch = amr->patch[0][lv][PID]->NPar;

            if ( amr->patch[0][lv][PID]->son == -1 )
            {
//             count the number of particles in the leaf patches
               NParInLeaf += NParThisPatch;


               EdgeL = amr->patch[0][lv][PID]->EdgeL;
               EdgeR = amr->patch[0][lv][PID]->EdgeR;

               for (int p=0; p<NParThisPatch; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

//                Check 5: one particle should has only one home patch
                  if ( ParHome[ParID] == true )
                  {
                     if ( PassAll )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, Time[lv], Step );
                        PassAll = false;
                     }

                     if ( PassCheck[4] )
                     {
                        Aux_Message( stderr, "Check 5: %4s  %2s  %7s  %10s\n", "Rank", "Lv", "PID", "ParID" );
                        PassCheck[4] = false;
                     }

                     Aux_Message( stderr, "Check 5: %4d  %2d  %7d  %10ld\n", MPI_Rank, lv, PID, ParID );
                  }

                  else
                     ParHome[ParID] = true;


//                Check 1: check whether or not all particles find their home patches
                  for (int d=0; d<3; d++)
                  {
                     if ( ParPos[d][ParID] < EdgeL[d]  ||  ParPos[d][ParID] >= EdgeR[d] )
                     {
                        if ( PassAll )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, Time[lv], Step );
                           PassAll = false;
                        }

                        if ( PassCheck[0] )
                        {
                           Aux_Message( stderr, "Check 1: %4s  %2s  %7s  %10s  %3s  %20s  %20s  %20s\n",
                                        "Rank", "Lv", "PID", "ParID", "Dim", "EdgeL", "EdgeR", "ParPos"  );
                           PassCheck[0] = false;
                        }

                        Aux_Message( stderr, "Check 1: %4d  %2d  %7d  %10ld  %3d  %20.13e  %20.13e  %20.13e\n",
                                     MPI_Rank, lv, PID, ParID, d, EdgeL[d], EdgeR[d], ParPos[d][ParID] );
                     }
                  } // for (int d=0; d<3; d++)


//                Check 4: no active particles should have mass < 0.0
                  if ( amr->Par->Mass[ParID] < 0.0 )
                  {
                     if ( PassAll )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, Time[lv], Step );
                        PassAll = false;
                     }

                     if ( PassCheck[3] )
                     {
                        Aux_Message( stderr, "Check 4: %4s  %2s  %7s  %10s  %20s  %20s  %20s  %20s\n",
                                     "Rank", "Lv", "PID", "ParID", "PosX", "PosY", "PosZ", "Mass"  );
                        PassCheck[3] = false;
                     }

                     Aux_Message( stderr, "Check 4: %4d  %2d  %7d  %10ld  %20.13e  %20.13e  %20.13e  %20.13e\n",
                                  MPI_Rank, lv, PID, ParID, ParPos[0][ParID], ParPos[1][ParID], ParPos[2][ParID],
                                  amr->Par->Mass[ParID] );
                  }


//                Check 11: particle types
                  bool CheckTypePass = true;

//                particle types must be recognizable
                  if ( amr->Par->Type[ParID] < (long_par)0  ||  amr->Par->Type[ParID] >= (long_par)PAR_NTYPE )
                     CheckTypePass = false;

//                only support tracer particles when disabling GRAVITY
#                 ifndef GRAVITY
                  if ( amr->Par->Type[ParID] != PTYPE_TRACER )
                     CheckTypePass = false;
#                 endif

//                must enable TRACER for tracer particles
#                 ifndef TRACER
                  if ( amr->Par->Type[ParID] == PTYPE_TRACER )
                     CheckTypePass = false;
#                 endif

//                tracer particles must be massless
#                 ifdef TRACER
                  if ( amr->Par->Type[ParID] == PTYPE_TRACER  &&  amr->Par->Mass[ParID] != (real_par)0.0 )
                     CheckTypePass = false;
#                 endif

                  if ( !CheckTypePass )
                  {
                     if ( PassAll )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, Time[lv], Step );
                        PassAll = false;
                     }

                     if ( PassCheck[10] )
                     {
                        Aux_Message( stderr, "Check 11: %4s  %2s  %7s  %10s  %10s  %20s\n",
                                     "Rank", "Lv", "PID", "ParID", "Type", "Mass" );
                        PassCheck[10] = false;
                     }

                     Aux_Message( stderr, "Check 11: %4d  %2d  %7d  %10ld  %10d  %20.13e\n",
                                  MPI_Rank, lv, PID, ParID, (int)amr->Par->Type[ParID], amr->Par->Mass[ParID] );
                  } // if ( !CheckTypePass )

               } // for (int p=0; p<NParThisPatch; p++)
            } // if ( amr->patch[0][lv][PID]->son == -1 )


//          Check 2: particles should only reside in the leaf patches
            else if ( NParThisPatch != 0 )
            {
               if ( PassAll )
               {
                  Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                               comment, __FUNCTION__, Time[lv], Step );
                  PassAll = false;
               }

               if ( PassCheck[1] )
               {
                  Aux_Message( stderr, "Check 2: %4s  %2s  %7s  %7s  %7s\n", "Rank", "Lv", "PID", "SonPID", "NPar" );
                  PassCheck[1] = false;
               }

               Aux_Message( stderr, "Check 2: %4d  %2d  %7d  %7d  %7d\n", MPI_Rank, lv, PID, amr->patch[0][lv][PID]->son, NParThisPatch );
            } // if ( amr->patch[0][lv][PID]->son == -1 ) ... else ...
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++); for (int lv=0; lv<NLEVEL; lv++)


//       Check 10: all patches (both real and buffer) should have NPar_Copy == -1
         for (int lv=0; lv<NLEVEL; lv++)
         for (int PID=0; PID<amr->num[lv]; PID++)
         {
            if ( amr->patch[0][lv][PID]->NPar_Copy != -1 )
            {
               if ( PassAll )
               {
                  Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                               comment, __FUNCTION__, Time[0], Step, MPI_Rank );
                  PassAll = false;
               }

               if ( PassCheck[9] )
               {
                  Aux_Message( stderr, "Check 9: lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                               lv, PID, amr->patch[0][lv][PID]->NPar_Copy );
                  PassCheck[9] = false;
               }
            }
         } // for lv, PID


//       Check 3: check the total number of active particles
         if ( NParInLeaf != amr->Par->NPar_Active )
         {
            if ( PassAll )
            {
               Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                            comment, __FUNCTION__, Time[0], Step, MPI_Rank );
               PassAll = false;
            }

            if ( PassCheck[2] )
            {
               Aux_Message( stderr, "Check 3: total number of active particles found in the leaf patches (%ld) != expect (%ld) !!\n",
                            NParInLeaf, amr->Par->NPar_Active );
               Aux_Message( stderr, "         (inactive + active particles = %ld)\n", amr->Par->NPar_AcPlusInac );
               PassCheck[2] = false;
            }
         }


//       Check 6: check if any active particle is homeless
         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
            if ( amr->Par->Mass[p] >= 0.0  &&  ParHome[p] == false )
            {
               if ( PassAll )
               {
                  Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                               comment, __FUNCTION__, Time[0], Step );
                  PassAll = false;
               }

               if ( PassCheck[5] )
               {
                  Aux_Message( stderr, "Check 6: %4s  %10s\n", "Rank", "ParID" );
                  PassCheck[5] = false;
               }

               Aux_Message( stderr, "Check 6: %4d  %10ld\n", MPI_Rank, p );
            }
         }


//       Check 7/8/9: consistency of the number of particles
         if ( amr->Par->NPar_Active + amr->Par->NPar_Inactive != amr->Par->NPar_AcPlusInac )
         {
            if ( PassAll )
            {
               Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                            comment, __FUNCTION__, Time[0], Step, MPI_Rank );
               PassAll = false;
            }

            if ( PassCheck[6] )
            {
               Aux_Message( stderr, "Check 7: NPar_Active (%ld) + NPar_Inactive (%ld) = %ld != NPar_AcPlusInac (%ld) !!\n",
                            amr->Par->NPar_Active, amr->Par->NPar_Inactive,
                            amr->Par->NPar_Active+amr->Par->NPar_Inactive, amr->Par->NPar_AcPlusInac );
               PassCheck[6] = false;
            }
         }


         if ( amr->Par->NPar_Active_AllRank != NPar_Active_AllRank_Expect )
         {
            if ( PassAll )
            {
               Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                            comment, __FUNCTION__, Time[0], Step, MPI_Rank );
               PassAll = false;
            }

            if ( PassCheck[7] )
            {
               Aux_Message( stderr, "Check 8: NPar_Active_AllRank (%ld) != expect (%ld) !!\n",
                            amr->Par->NPar_Active_AllRank, NPar_Active_AllRank_Expect );
               PassCheck[7] = false;
            }
         }


         long NPar_Lv_Sum=0;
         for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv_Sum += amr->Par->NPar_Lv[lv];

         if ( NPar_Lv_Sum != amr->Par->NPar_Active )
         {
            if ( PassAll )
            {
               Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                            comment, __FUNCTION__, Time[0], Step, MPI_Rank );
               PassAll = false;
            }

            if ( PassCheck[8] )
            {
               Aux_Message( stderr, "Check 9: NPar_Lv_Sum (%ld) != expect (%ld) !!\n",
                            NPar_Lv_Sum, amr->Par->NPar_Active );
               PassCheck[8] = false;
            }
         }
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &PassAll,        1, MPI_INT, TargetRank, MPI_COMM_WORLD );
      MPI_Bcast(  PassCheck, NCheck, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)


   if ( PassAll )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at Time = %13.7e, Step = %ld\n", comment, __FUNCTION__, Time[0], Step );
   }

   else
      MPI_Exit();


   delete [] ParHome;

} // FUNCTION : Par_Aux_Check_Particle



#endif // #ifdef PARTICLE
