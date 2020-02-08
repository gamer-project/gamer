#include "GAMER.h"

#ifdef MHD


// tolerated errors
#ifdef FLOAT8
#  define ERROR_TOLERANCE  1.0e-11
#else
#  define ERROR_TOLERANCE  1.0e-5f
#endif


static void CheckError( int &Pass, const real B, const real BSib, const int lv, const int PID,
                        const int Sib, const int idx1, const int idx2, const char *comment );




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Aux_Check_InterfaceB
// Description :  Verify that adjacent patches have the same B field on their shared interfaces
//
// Note        :  1. Tolerated error is set by the symbolic constant ERROR_TOLERANCE
//                   --> Should be close to the round-off errors
//                2. After data restriction, Flu_FixUp_Restrict() will also copy B field between nearby coarse patches
//                   to ensure this B field consistency
//                3. Invoked by Aux_Check()
//
// Parameter   :  comment : You can put the location where this function is invoked in this string
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void MHD_Aux_Check_InterfaceB( const char *comment )
{

   int Pass = true;
   real B, BSib;

// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      const int MagSg = amr->MagSg[lv];

      for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
      {
         if ( MPI_Rank == TargetRank )
         {
            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            for (int s=0; s<6; s++)
            {
               const int SibPID = amr->patch[0][lv][PID]->sibling[s];

               if ( SibPID >= 0 )
               {
                  const int Offset    = ( s   %2)*PS1;
                  const int OffsetSib = ((s+1)%2)*PS1;

                  switch ( s )
                  {
                     case 0: case 1:
                        for (int k=0; k<PS1; k++)
                        for (int j=0; j<PS1; j++)
                        {
                           B     = amr->patch[MagSg][lv][   PID]->magnetic[MAGX][ IDX321_BX(Offset,    j, k, PS1, PS1) ];
                           BSib  = amr->patch[MagSg][lv][SibPID]->magnetic[MAGX][ IDX321_BX(OffsetSib, j, k, PS1, PS1) ];

                           CheckError( Pass, B, BSib, lv, PID, s, j, k, comment );
                        }
                        break;

                     case 2: case 3:
                        for (int k=0; k<PS1; k++)
                        for (int i=0; i<PS1; i++)
                        {
                           B     = amr->patch[MagSg][lv][   PID]->magnetic[MAGY][ IDX321_BY(i, Offset,    k, PS1, PS1) ];
                           BSib  = amr->patch[MagSg][lv][SibPID]->magnetic[MAGY][ IDX321_BY(i, OffsetSib, k, PS1, PS1) ];

                           CheckError( Pass, B, BSib, lv, PID, s, i, k, comment );
                        }
                        break;

                     case 4: case 5:
                        for (int j=0; j<PS1; j++)
                        for (int i=0; i<PS1; i++)
                        {
                           B     = amr->patch[MagSg][lv][   PID]->magnetic[MAGZ][ IDX321_BZ(i, j, Offset,    PS1, PS1) ];
                           BSib  = amr->patch[MagSg][lv][SibPID]->magnetic[MAGZ][ IDX321_BZ(i, j, OffsetSib, PS1, PS1) ];

                           CheckError( Pass, B, BSib, lv, PID, s, i, j, comment );
                        }
                        break;

                     default:
                        Aux_Error( ERROR_INFO, "s = %d !!\n" );
                        break;
                  } // switch( s )
               } // if ( SibPID >= 0 )
            } // for (int s=0; s<6; s++) for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // if ( MPI_Rank == TargetRank )

         MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );
         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( Pass  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "\"%s\" : <%s> PASSED, Time = %13.7e, Step = %ld \n",
                   comment, __FUNCTION__, Time[0], Step );

} // FUNCTION : MHD_Aux_Check_InterfaceB



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckError
// Description :  Check if error exceeds a given tolerance
//
// Note        :  1. Invoked by MHD_Aux_Check_InterfaceB()
//                2. Tolerated error is set by the symbolic constant ERROR_TOLERANCE
//
// Parameter   :  Pass    : Whether this test has failed or not (call-by-reference)
//                B       : B field on the target patch
//                BSib    : B field on the sibling patch
//                lv      : AMR level
//                PID     : Patch ID
//                Sib     : Sibling direction
//                idx1/2  : Cell index on the target face
//                comment : Output error message
//
// Return      :  Pass
//-------------------------------------------------------------------------------------------------------
void CheckError( int &Pass, const real B, const real BSib, const int lv, const int PID,
                 const int Sib, const int idx1, const int idx2, const char *comment )
{

   const real Error = (real)2.0*FABS(B - BSib)/( FABS(B) + FABS(BSib) );

// be careful about the case where both B and BSib are zero
   if (  Error > ERROR_TOLERANCE  &&  !( B == 0.0 && BSib == 0.0 )  )
   {
      if ( Pass )
      {
         Aux_Message( stderr, "\"%s\" : <%s> FAILED, Time = %13.7e, Step = %ld !!\n",
                      comment, "MHD_Aux_Check_InterfaceB", Time[0], Step );
         Aux_Message( stderr, "%4s   %2s   %7s   %3s   %4s   %4s   %21s   %21s   %21s\n",
                      "Rank", "Lv", "PID", "Sib", "Idx1", "Idx2", "B", "BSib", "Error" );
         Pass = false;
      }

      Aux_Message( stderr, "%4d   %2d   %7d   %3d   %4d   %4d   %21.14e   %21.14e   %21.14e\n",
                   MPI_Rank, lv, PID, Sib, idx1, idx2, B, BSib, Error );
   }

} // FUNCTION : CheckError



#endif // #ifdef MHD
