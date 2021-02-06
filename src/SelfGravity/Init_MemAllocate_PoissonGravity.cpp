#include "GAMER.h"

#if ( !defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_PoissonGravity
// Description :  Allocate memory for the Poisson and Gravity solvers
//
// Note        :  Only work when using CPUs only
//
// Parameter   :  Pot_NPG  : Number of patch groups to be evaluated at a time
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_PoissonGravity( const int Pot_NPG )
{

   const int Pot_NP = 8*Pot_NPG;

   for (int t=0; t<2; t++)
   {
      h_Rho_Array_P     [t] = new real   [Pot_NP][RHO_NXT][RHO_NXT][RHO_NXT];
      h_Pot_Array_P_In  [t] = new real   [Pot_NP][POT_NXT][POT_NXT][POT_NXT];
      h_Pot_Array_P_Out [t] = new real   [Pot_NP][GRA_NXT][GRA_NXT][GRA_NXT];
#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_G [t] = new real   [Pot_NP][USG_NXT_G][USG_NXT_G][USG_NXT_G];
      h_Flu_Array_USG_G [t] = new real   [Pot_NP][GRA_NIN-1][PS1][PS1][PS1];
#     endif
      h_Flu_Array_G     [t] = new real   [Pot_NP][GRA_NIN  ][PS1][PS1][PS1];

      if ( OPT__EXT_ACC  ||  OPT__EXT_POT )
      h_Corner_Array_PGT[t] = new double [Pot_NP][3];

#     ifdef DUAL_ENERGY
      h_DE_Array_G      [t] = new char   [Pot_NP][PS1][PS1][PS1];
#     endif

#     ifdef MHD
      h_Emag_Array_G    [t] = new real   [Pot_NP][PS1][PS1][PS1];
#     endif

      h_Pot_Array_T     [t] = new real   [Pot_NP][ CUBE(GRA_NXT) ];
   }


// external potential table
   if ( OPT__EXT_POT == EXT_POT_TABLE ) {

      const long TableSize = (long)EXT_POT_TABLE_NPOINT[0]*EXT_POT_TABLE_NPOINT[1]*EXT_POT_TABLE_NPOINT[2];

      if ( TableSize <= 0L )
         Aux_Error( ERROR_INFO, "External potential table size = %ld <= 0 (EXT_POT_TABLE_NPOINT = [%d, %d, %d]) !!\n",
                    TableSize, EXT_POT_TABLE_NPOINT[0], EXT_POT_TABLE_NPOINT[1], EXT_POT_TABLE_NPOINT[2] );

      h_ExtPotTable = new real [TableSize];
   }

} // FUNCTION : Init_MemAllocate_PoissonGravity



#endif // #if ( !defined GPU  &&  defined GRAVITY )
