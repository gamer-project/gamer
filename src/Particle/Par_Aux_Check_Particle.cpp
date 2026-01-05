#include "GAMER.h"

#ifdef PARTICLE


static void Check_FindHomePatch      ( int &PassAll, int &PassOne, const char *comment, const int lv,
                                       const int PID, const int ParID, double *EdgeL, double *EdgeR,
                                       const real_par *ParPos[] );
static void Check_InLeafPatch        ( int &PassAll, int &PassOne, const char *comment, const int lv,
                                       const int PID, const int NParThisPatch );
static void Check_NActive            ( int &PassAll, int &PassOne, const char *comment, const long NParInLeaf );
static void Check_InactiveMass       ( int &PassAll, int &PassOne, const char *comment, const int lv,
                                       const int PID, const int ParID, const real_par *ParPos[] );
static void Check_OneHomePatch       ( int &PassAll, int &PassOne, const char *comment, const int lv,
                                       const int PID, const int ParID, bool *ParHome );
static void Check_ActiveHome         ( int &PassAll, int &PassOne, const char *comment, const bool *ParHome );
static void Check_NPar_AcPlusInac    ( int &PassAll, int &PassOne, const char *comment );
static void Check_NPar_Active_AllRank( int &PassAll, int &PassOne, const char *comment, const long NPar_Active_AllRank_Expect );
static void Check_NPar_Lv_Sum        ( int &PassAll, int &PassOne, const char *comment );
static void Check_NPar_Copy          ( int &PassAll, int &PassOne, const char *comment );
static void Check_ValidType          ( int &PassAll, int &PassOne, const char *comment, const int lv,
                                       const int PID, const int ParID );
static void Check_ValidUID           ( int &PassAll, int &PassOne, const char *comment );
static void Check_UniqueUID          ( int &PassAll, int &PassOne, const char *comment, int *Par_CountUID );



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
//                12.  Particle UID is valid
//                13.  Particle UID is unique
//
// Note        :  None
//
// Parameter   :  comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Par_Aux_Check_Particle( const char *comment )
{

   const int       NCheck    = 13;
   const real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };

   int     PassAll    = true;
   long    NParInLeaf = 0;
   bool   *ParHome    = new bool [amr->Par->NPar_AcPlusInac];  // true/false --> particle has home/is homeless
   int     NParThisPatch;
   long    ParID, NPar_Active_AllRank_Expect;
   double *EdgeL, *EdgeR;
   int     PassCheck[NCheck];
   int    *Par_CountUID = new int [amr->Par->NextUID];


// initialize the check list
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   ParHome[p] = false;

   for (int t=0; t<NCheck; t++)  PassCheck[t] = true;

   for (long uid=0; uid<amr->Par->NextUID; uid++)   Par_CountUID[uid] = 0;


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

                  Check_OneHomePatch( PassAll, PassCheck[4], comment, lv, PID, ParID, ParHome );

                  Check_FindHomePatch( PassAll, PassCheck[0], comment, lv, PID, ParID, EdgeL, EdgeR, ParPos );

                  Check_InactiveMass( PassAll, PassCheck[3], comment, lv, PID, ParID, ParPos );

                  Check_ValidType( PassAll, PassCheck[10], comment, lv, PID, ParID );

               } // for (int p=0; p<NParThisPatch; p++)
            } // if ( amr->patch[0][lv][PID]->son == -1 )
            else
            {
               Check_InLeafPatch( PassAll, PassCheck[1], comment, lv, PID, NParThisPatch );
            } // if ( amr->patch[0][lv][PID]->son == -1 ) ... else ...
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++); for (int lv=0; lv<NLEVEL; lv++)


         Check_NPar_Copy( PassAll, PassCheck[9], comment );

         Check_NActive( PassAll, PassCheck[2], comment, NParInLeaf );

         Check_ActiveHome( PassAll, PassCheck[5], comment, ParHome );

         Check_NPar_AcPlusInac( PassAll, PassCheck[6], comment );

         Check_NPar_Active_AllRank( PassAll, PassCheck[7], comment, NPar_Active_AllRank_Expect );

         Check_NPar_Lv_Sum( PassAll, PassCheck[8], comment );

         Check_ValidUID( PassAll, PassCheck[11], comment );

         Check_UniqueUID( PassAll, PassCheck[12], comment, Par_CountUID );
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast(  Par_CountUID, amr->Par->NextUID, MPI_INT, TargetRank, MPI_COMM_WORLD );
      MPI_Bcast( &PassAll,      1,                 MPI_INT, TargetRank, MPI_COMM_WORLD );
      MPI_Bcast(  PassCheck,    NCheck,            MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)


   if ( ! PassAll )   MPI_Exit();

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "\"%s\" : <%s> PASSED at Time = %13.7e, Step = %ld\n", comment, __FUNCTION__, Time[0], Step );


   delete [] ParHome;
   delete [] Par_CountUID;

} // FUNCTION : Par_Aux_Check_Particle



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_FindHomePatch
// Description :  Check if the particle reside in its home patch
//
// Note        :  None
//
// Parameter   :  PassAll : Pass all the checkes or not
//                PassOne : Pass this check or not
//                comment : You can put the location where this function is invoked in this string
//                lv      : Target refinement level
//                PID     : Target patch ID
//                ParID   : Target particle ID
//                EdgeL   : Target patch left edges
//                EdgeR   : Target patch right edges
//                ParPos  : Particle postion
//-------------------------------------------------------------------------------------------------------
void Check_FindHomePatch( int &PassAll, int &PassOne, const char *comment, const int lv, const int PID,
                          const int ParID, double *EdgeL, double *EdgeR, const real_par *ParPos[] )
{

   for (int d=0; d<3; d++)
   {
      if ( ParPos[d][ParID] >= EdgeL[d]  &&  ParPos[d][ParID] < EdgeR[d] )   continue;

      if ( PassAll )
         Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                      comment, __FUNCTION__, Time[lv], Step );

      if ( PassOne )
         Aux_Message( stderr, "Check 1: %4s  %2s  %7s  %10s  %3s  %20s  %20s  %20s\n",
                      "Rank", "Lv", "PID", "ParID", "Dim", "EdgeL", "EdgeR", "ParPos"  );

      Aux_Message( stderr, "Check 1: %4d  %2d  %7d  %10ld  %3d  %20.13e  %20.13e  %20.13e\n",
                   MPI_Rank, lv, PID, ParID, d, EdgeL[d], EdgeR[d], ParPos[d][ParID] );

      PassAll = false;
      PassOne = false;
   } // for (int d=0; d<3; d++)

} // FUNCTION : Check_FindHomePatch



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_InLeafPatch
// Description :  Check if the particle reside in the leaf patches
//
// Note        :  None
//
// Parameter   :  PassAll       : Pass all the checkes or not
//                PassOne       : Pass this check or not
//                comment       : You can put the location where this function is invoked in this string
//                lv            : Target refinement level
//                PID           : Target patch ID
//                NParThisPatch : Number of particle in target patch
//-------------------------------------------------------------------------------------------------------
void Check_InLeafPatch( int &PassAll, int &PassOne, const char *comment, const int lv, const int PID,
                        const int NParThisPatch )
{

   if ( NParThisPatch == 0 )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                   comment, __FUNCTION__, Time[lv], Step );

   if ( PassOne )
      Aux_Message( stderr, "Check 2: %4s  %2s  %7s  %7s  %7s\n", "Rank", "Lv", "PID", "SonPID", "NPar" );

   Aux_Message( stderr, "Check 2: %4d  %2d  %7d  %7d  %7d\n", MPI_Rank, lv, PID, amr->patch[0][lv][PID]->son, NParThisPatch );

   PassAll = false;
   PassOne = false;

} // FUNCTION : CheckInLeafPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_NActive
// Description :  Check if there are missing or redundant particles
//
// Note        :  None
//
// Parameter   :  PassAll    : Pass all the checkes or not
//                PassOne    : Pass this check or not
//                comment    : You can put the location where this function is invoked in this string
//                NParInLeaf : Number of particles in leaf patches
//-------------------------------------------------------------------------------------------------------
void Check_NActive( int &PassAll, int &PassOne, const char *comment, const long NParInLeaf )
{

   if ( NParInLeaf == amr->Par->NPar_Active )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                   comment, __FUNCTION__, Time[0], Step, MPI_Rank );

   if ( PassOne )
   {
      Aux_Message( stderr, "Check 3: total number of active particles found in the leaf patches (%ld) != expect (%ld) !!\n",
                   NParInLeaf, amr->Par->NPar_Active );
      Aux_Message( stderr, "         (inactive + active particles = %ld)\n", amr->Par->NPar_AcPlusInac );
   }

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_NActive



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_InactiveMass
// Description :  Check if the inactive particle has negative mass
//
// Note        :  None
//
// Parameter   :  PassAll : Pass all the checkes or not
//                PassOne : Pass this check or not
//                comment : You can put the location where this function is invoked in this string
//                lv      : Target refinement level
//                PID     : Target patch ID
//                ParID   : Target particle ID
//                ParPos  : Particle postion
//-------------------------------------------------------------------------------------------------------
void Check_InactiveMass( int &PassAll, int &PassOne, const char *comment, const int lv, const int PID,
                         const int ParID, const real_par *ParPos[] )
{

   if ( amr->Par->Mass[ParID] >= 0.0 )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                   comment, __FUNCTION__, Time[lv], Step );

   if ( PassOne )
      Aux_Message( stderr, "Check 4: %4s  %2s  %7s  %10s  %20s  %20s  %20s  %20s\n",
                   "Rank", "Lv", "PID", "ParID", "PosX", "PosY", "PosZ", "Mass"  );

   Aux_Message( stderr, "Check 4: %4d  %2d  %7d  %10ld  %20.13e  %20.13e  %20.13e  %20.13e\n",
                MPI_Rank, lv, PID, ParID, ParPos[0][ParID], ParPos[1][ParID], ParPos[2][ParID],
                amr->Par->Mass[ParID] );

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_InactiveMass



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_OneHomePatch
// Description :  Check if particle has only one home patch
//
// Note        :  None
//
// Parameter   :  PassAll : Pass all the checkes or not
//                PassOne : Pass this check or not
//                comment : You can put the location where this function is invoked in this string
//                lv      : Target refinement level
//                PID     : Target patch ID
//                ParID   : Target particle ID
//                ParHome : Array of particle has home or not
//-------------------------------------------------------------------------------------------------------
void Check_OneHomePatch( int &PassAll, int &PassOne, const char *comment, const int lv, const int PID,
                         const int ParID, bool *ParHome )
{

   if ( ParHome[ParID] == false ) { ParHome[ParID] = true; return; }

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                   comment, __FUNCTION__, Time[lv], Step );

   if ( PassOne )
      Aux_Message( stderr, "Check 5: %4s  %2s  %7s  %10s\n", "Rank", "Lv", "PID", "ParID" );

   Aux_Message( stderr, "Check 5: %4d  %2d  %7d  %10ld\n", MPI_Rank, lv, PID, ParID );

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_OneHomePatch



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_ActiveHome
// Description :  Check if the active particle has only one home patch
//
// Note        :  1. Must be executed after Check_OneHomePatch
//
// Parameter   :  PassAll : Pass all the checkes or not
//                PassOne : Pass this check or not
//                comment : You can put the location where this function is invoked in this string
//                ParHome : Array of particle has home or not
//-------------------------------------------------------------------------------------------------------
void Check_ActiveHome( int &PassAll, int &PassOne, const char *comment, const bool *ParHome )
{

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
      if ( amr->Par->Mass[p] < 0.0  ||  ParHome[p] == true )   continue;

      if ( PassAll )
         Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                      comment, __FUNCTION__, Time[0], Step );

      if ( PassOne )
         Aux_Message( stderr, "Check 6: %4s  %10s\n", "Rank", "ParID" );

      Aux_Message( stderr, "Check 6: %4d  %10ld\n", MPI_Rank, p );

      PassAll = false;
      PassOne = false;
   } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)

} // FUNCTION : Check_ActiveHome



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_AcPlusInac
// Description :  Check if NPar_AcPlusInac == NPar_Active + NPar_Inactive
//
// Note        :  None
//
// Parameter   :  PassAll : Pass all the checkes or not
//                PassOne : Pass this check or not
//                comment : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Check_NPar_AcPlusInac( int &PassAll, int &PassOne, const char *comment )
{

   if ( amr->Par->NPar_Active + amr->Par->NPar_Inactive == amr->Par->NPar_AcPlusInac )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                   comment, __FUNCTION__, Time[0], Step, MPI_Rank );

   if ( PassOne )
      Aux_Message( stderr, "Check 7: NPar_Active (%ld) + NPar_Inactive (%ld) = %ld != NPar_AcPlusInac (%ld) !!\n",
                   amr->Par->NPar_Active, amr->Par->NPar_Inactive,
                   amr->Par->NPar_Active+amr->Par->NPar_Inactive, amr->Par->NPar_AcPlusInac );

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_NPar_AcPlusInac



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_NPar_Active_AllRank
// Description :  Check if NPar_Active_AllRank = sum(NPar_Active, All ranks)
//
// Note        :  None
//
// Parameter   :  PassAll                    : Pass all the checkes or not
//                PassOne                    : Pass this check or not
//                comment                    : You can put the location where this function is invoked in this string
//                NPar_Active_AllRank_Expect : Expected number of active particles
//-------------------------------------------------------------------------------------------------------
void Check_NPar_Active_AllRank( int &PassAll, int &PassOne, const char *comment, const long NPar_Active_AllRank_Expect )
{

   if ( amr->Par->NPar_Active_AllRank == NPar_Active_AllRank_Expect )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                   comment, __FUNCTION__, Time[0], Step, MPI_Rank );

   if ( PassOne )
      Aux_Message( stderr, "Check 8: NPar_Active_AllRank (%ld) != expect (%ld) !!\n",
                   amr->Par->NPar_Active_AllRank, NPar_Active_AllRank_Expect );

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_NPar_Active_AllRank



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_NPar_Lv_Sum
// Description :  Check if NPar_Active = sum(NPar_Lv, all levels)
//
// Note        :  None
//
// Parameter   :  PassAll  : Pass all the checkes or not
//                PassOne  : Pass this check or not
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Check_NPar_Lv_Sum( int &PassAll, int &PassOne, const char *comment )
{

   long NPar_Lv_Sum=0;
   for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv_Sum += amr->Par->NPar_Lv[lv];

   if ( NPar_Lv_Sum == amr->Par->NPar_Active )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                   comment, __FUNCTION__, Time[0], Step, MPI_Rank );

   if ( PassOne )
      Aux_Message( stderr, "Check 9: NPar_Lv_Sum (%ld) != expect (%ld) !!\n",
                   NPar_Lv_Sum, amr->Par->NPar_Active );

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_NPar_Lv_Sum



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_NPar_Copy
// Description :  Check if all patches have NPar_Copy == -1
//
// Note        :  None
//
// Parameter   :  PassAll  : Pass all the checkes or not
//                PassOne  : Pass this check or not
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Check_NPar_Copy( int &PassAll, int &PassOne, const char *comment )
{

   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->num[lv]; PID++)
   {
      if ( amr->patch[0][lv][PID]->NPar_Copy == -1 )   continue;

      if ( PassAll )
         Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld, Rank = %d !!\n",
                      comment, __FUNCTION__, Time[0], Step, MPI_Rank );

      if ( PassOne )
         Aux_Message( stderr, "Check 10: lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                      lv, PID, amr->patch[0][lv][PID]->NPar_Copy );

      PassAll = false;
      PassOne = false;
   } // for lv, PID

} // FUNCTION : Check_NPar_Copy



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_ValidType
// Description :  Check if particle type is valid
//
// Note        :  None
//
// Parameter   :  PassAll : Pass all the checkes or not
//                PassOne : Pass this check or not
//                comment : You can put the location where this function is invoked in this string
//                lv      : Target refinement level
//                PID     : Target patch ID
//                ParID   : Target particle ID
//-------------------------------------------------------------------------------------------------------
void Check_ValidType( int &PassAll, int &PassOne, const char *comment, const int lv, const int PID, const int ParID )
{

   bool CheckTypePass = true;

// particle types must be recognizable
   if ( amr->Par->Type[ParID] < (long_par)0  ||  amr->Par->Type[ParID] >= (long_par)PAR_NTYPE )
      CheckTypePass = false;

// only support tracer particles when disabling GRAVITY
#  ifndef GRAVITY
   if ( amr->Par->Type[ParID] != PTYPE_TRACER )
      CheckTypePass = false;
#  endif

// must enable TRACER for tracer particles
#  ifndef TRACER
   if ( amr->Par->Type[ParID] == PTYPE_TRACER )
      CheckTypePass = false;
#  endif

// tracer particles must be massless
#  ifdef TRACER
   if ( amr->Par->Type[ParID] == PTYPE_TRACER  &&  amr->Par->Mass[ParID] != (real_par)0.0 )
      CheckTypePass = false;
#  endif

   if ( CheckTypePass )   return;

   if ( PassAll )
      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                   comment, __FUNCTION__, Time[lv], Step );

   if ( PassOne )
      Aux_Message( stderr, "Check 11: %4s  %2s  %7s  %10s  %10s  %20s\n",
                   "Rank", "Lv", "PID", "ParID", "Type", "Mass" );

   Aux_Message( stderr, "Check 11: %4d  %2d  %7d  %10ld  %10d  %20.13e\n",
                MPI_Rank, lv, PID, ParID, (int)amr->Par->Type[ParID], amr->Par->Mass[ParID] );

   PassAll = false;
   PassOne = false;

} // FUNCTION : Check_ValidType



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_ValidUID
// Description :  Check if particle has valid UID
//
// Note        :  None
//
// Parameter   :  PassAll  : Pass all the checkes or not
//                PassOne  : Pass this check or not
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Check_ValidUID( int &PassAll, int &PassOne, const char *comment )
{

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
      if ( amr->Par->PUid[p] > (long_par)0  &&  amr->Par->PUid[p] < amr->Par->NextUID )  continue;

//    exclude inactive particles, which could have been skipped during UID assignment
      if ( amr->Par->Mass[p] < (real_par)0.0 )   continue;

      if ( PassAll )
         Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                      comment, __FUNCTION__, Time[0], Step );

      if ( PassOne )
         Aux_Message( stderr, "Check 12: %4s  %10s  %10s  %10s\n", "Rank", "ParID", "PUid", "NextUID" );

      Aux_Message( stderr, "Check 12: %4d  %10ld  %10ld  %10ld\n", MPI_Rank, p, (long)amr->Par->PUid[p], amr->Par->NextUID );

      PassAll = false;
      PassOne = false;
   } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)

} // FUNCTION : Check_ValidUID



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_UniqueUID
// Description :  Check if particle has unique UID
//
// Note        :  None
//
// Parameter   :  PassAll      : Pass all the checkes or not
//                PassOne      : Pass this check or not
//                comment      : You can put the location where this function is invoked in this string
//                Par_CountUID : Array of the counted number of each particle UID
//-------------------------------------------------------------------------------------------------------
void Check_UniqueUID( int &PassAll, int &PassOne, const char *comment, int *Par_CountUID )
{

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
//    exclude inactive particles, which could have duplicated UIDs
      if ( amr->Par->Mass[p] < (real_par)0.0 )   continue;

      const long uid = amr->Par->PUid[p];
      Par_CountUID[uid]++;

      if ( Par_CountUID[uid] <= 1 )   continue;

      if ( PassAll )
         Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n",
                      comment, __FUNCTION__, Time[0], Step );

      if ( PassOne )
         Aux_Message( stderr, "Check 13: %4s  %10s  %10s  %10s\n", "Rank", "ParUID", "Number", "NextUID" );

      Aux_Message( stderr, "Check 13: %4d  %10ld  %10d  %10ld\n", MPI_Rank, uid, Par_CountUID[uid], amr->Par->NextUID );

      PassAll = false;
      PassOne = false;
   } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)

} // FUNCTION : Check_UniqueUID



#endif // #ifdef PARTICLE
