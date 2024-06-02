#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Restrict
// Description :  Verify the restriction condition at the input level
//                --> Output grids at level lv that violate the restriction condition
//
// Note        :  1. Restriction condition --> coarse-grid data = average of fine-grid data
//                2. Not supported in the LOAD_BALANCE mode
//
// Parameter   :  lv       : Target refinement level
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Restrict( const int lv, const char *comment )
{

// check
#  ifdef LOAD_BALANCE
   if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : check \"%s\" is not supported in LOAD_BALANCE and has been disabled !!\n",
                   __FUNCTION__ );

   OPT__CK_RESTRICT = false;
   return;
#  endif

   if ( lv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : function \"%s\" should NOT be applied to the finest level !!\n",
                   __FUNCTION__ );
      return;
   }


   const int FSg = amr->FluSg[lv+1];
   const int CSg = amr->FluSg[lv  ];

#  ifdef FLOAT8
   const double TolErr = 1.0e-13;
#  else
   const double TolErr = 1.0e-5;
#  endif
   int Pass = true;

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
// convert between phase/dens and re/im
   const bool ConvertWaveToFluid = ( !amr->use_wave_flag[lv] && amr->use_wave_flag[lv+1] );
#  endif

   int    SonPID0, SonPID, ii0, jj0, kk0, ii, jj, kk;
   double ResData[NCOMP_TOTAL][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            SonPID0 = amr->patch[0][lv][PID]->son;

            if ( SonPID0 != -1 )
            {
//             initialize the restrict-data array as zero
               for (int v=0; v<NCOMP_TOTAL; v++)
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
                  ResData[v][k][j][i] = 0.0;


//             fill up the restrict-data array
               for (int LocalID=0; LocalID<8; LocalID++)
               {
                  SonPID   = SonPID0 + LocalID;
                  ii0      = TABLE_02( LocalID, 'x', 0, PATCH_SIZE/2 );
                  jj0      = TABLE_02( LocalID, 'y', 0, PATCH_SIZE/2 );
                  kk0      = TABLE_02( LocalID, 'z', 0, PATCH_SIZE/2 );

                  for (int v=0; v<NCOMP_TOTAL; v++)   {
                  for (int k=0; k<PATCH_SIZE; k++)    {  kk = kk0 + k/2;
                  for (int j=0; j<PATCH_SIZE; j++)    {  jj = jj0 + j/2;
                  for (int i=0; i<PATCH_SIZE; i++)    {  ii = ii0 + i/2;

                     ResData[v][kk][jj][ii] += 0.125*amr->patch[FSg][lv+1][SonPID]->fluid[v][k][j][i];

                  }}}}
               }


//             compare the data of the restrict-data array and the data stored in the patch pointers
               double Err, u;

               for (int v=0; v<NCOMP_TOTAL; v++)
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
                  u = amr->patch[CSg][lv][PID]->fluid[v][k][j][i];

#                 if ( ELBDM_SCHEME == ELBDM_HYBRID )
//                to convert from wave to fluid, store the phase in the REAL component and ignore the imaginary part
                  if ( ConvertWaveToFluid  &&  v == REAL  &&  v == PHAS ) {
                     ResData[v][k][j][i] = ELBDM_UnwrapPhase( u, SATAN2(ResData[IMAG][k][j][i], ResData[REAL][k][j][i]) );
                  }
#                 endif

                  Err = fabs(  ( u - ResData[v][k][j][i] ) / ResData[v][k][j][i]  );

#                 if ( ELBDM_SCHEME == ELBDM_HYBRID )
//                skip stub component
                  if ( ConvertWaveToFluid  &&  v == IMAG  &&  v == STUB )
                  {
                     Err = 0.0;
                  }
#                 endif

                  if ( Err > TolErr )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Aux_Message( stderr, "%4s\t%7s\t%19s\t%10s\t%7s\t%14s\t%14s\t%14s\n",
                                     "Rank", "PID", "Patch Corner", "Grid ID", "Var", "Patch Data",
                                     "Restrict Data", "Error" );

                        Pass = false;
                     }

                     Aux_Message( stderr,"%4d\t%7d\t(%10d,%10d,%10d)\t(%2d,%2d,%2d)\t%7d\t%14.7e\t%14.7e\t%14.7e\n",
                                  MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                                 amr->patch[0][lv][PID]->corner[1],
                                                 amr->patch[0][lv][PID]->corner[2],
                                  i, j, k, v, u, ResData[v][k][j][i], Err );

                  } // if ( Err > TolErr )
               } // i,j,k,v

            } // if ( SonPID0 != -1 )
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)


   if ( Pass )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld\n",
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

} // Aux_Check_Restrict
