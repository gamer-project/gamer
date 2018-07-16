#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_NormalizePassive
// Description :  Verify that the passive scalars are properly normalized
//                --> sum(passive scalar mass density) = gas_mass_density
//
// Note        :  1. Only check variables stored in PassiveNorm_VarIdx[]
//                2. This check will likely fail if OPT__NORMALIZE_PASSIVE is disabled
//                3. Target passive scalars to be normalized are determined when invoking AddField()
//
// Parameter   :  lv       : Target refinement level
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_NormalizePassive( const int lv, const char *comment )
{


// check
#  ifndef DENS
   OPT__CK_NORMALIZE_PASSIVE = false;

   if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : \"OPT__CK_NORMALIZE_PASSIVE\" is supported only when DENS is defined !!\n" );

   return;
#  endif

   if ( PassiveNorm_NVar <= 0 )
   {
      OPT__CK_NORMALIZE_PASSIVE = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : \"OPT__CK_NORMALIZE_PASSIVE\" is disabled since there are no passive scalars to be normalized !!\n" );

      return;
   }


#  ifdef FLOAT8
   const double TolErr = 1.0e-13;
#  else
   const double TolErr = 1.0e-5;
#  endif
   const int    FluSg  = amr->FluSg[lv];

   double GasDens, Sum, RelErr;
   int    Pass = true;


   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
//             avoid compilation error for models without DENS
#              ifdef DENS
               GasDens = amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];
#              endif
               Sum     = 0;

               for (int v=0; v<PassiveNorm_NVar; v++)
               Sum += amr->patch[FluSg][lv][PID]->fluid[ NCOMP_FLUID + PassiveNorm_VarIdx[v] ][k][j][i];

               RelErr = fabs( (Sum-GasDens)/GasDens );

               if ( RelErr > TolErr )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s  %7s  %7s", "Rank", "PID", "(i,j,k)" );
                     for (int v=0; v<PassiveNorm_NVar; v++)
                     Aux_Message( stderr, "  %13s", FieldLabel[ NCOMP_FLUID + PassiveNorm_VarIdx[v] ] );
                     Aux_Message( stderr, "  %13s  %13s  %13s", "Sum", "GasDens", "RelErr" );
                     Aux_Message( stderr, "\n" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d  %7d  (%d,%d,%d)", MPI_Rank, PID, i, j, k );
                  for (int v=0; v<PassiveNorm_NVar; v++)
                  Aux_Message( stderr, "  %13.7e", amr->patch[FluSg][lv][PID]->fluid[ NCOMP_FLUID + PassiveNorm_VarIdx[v] ][k][j][i] );
                  Aux_Message( stderr, "  %13.7e  %13.7e  %13.7e", Sum, GasDens, RelErr );
                  Aux_Message( stderr, "\n" );
               } // if ( RelErr > TolErr )
            } // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)


   if ( Pass )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld \n",
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

} // FUNCTION : Aux_Check_NormalizePassive
