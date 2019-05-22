#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_FixUp_Electric
// Description :  Use the fine-grid electric field on the edges of coarse-fine boundaries to correct the
//                coarse-grid magnetic field
//
// Note        :  1. Invoked by Flu_FixUp()
//                2. Only correct B field on the coarse-coarse interfaces
//                   --> Does NOT correct B field on the coarse-fine interfaces since that has been done
//                       by Flu_FixUp_Restrict()
//
// Parameter   :  lv : Targeted coarse level
//-------------------------------------------------------------------------------------------------------
void MHD_FixUp_Electric( const int lv )
{

   const real _dh             = (real)1.0 / amr->dh[lv];
   const int  MagSg           = amr->MagSg[lv];

   const int  PS1M1_PS1       = PS1M1*PS1;
   const int  PS1P1_PS1       = PS1P1*PS1;
   const int  PS1_PS1         = SQR( PS1 );
   const int  PS1_PS1_PS1     = CUBE( PS1 );
   const int  PS1M1_PS1P1     = PS1M1*PS1P1;
   const int  PS1_PS1M1_PS1P1 = PS1*PS1M1_PS1P1;


// check
#  ifdef GAMER_DEBUG
   if ( !amr->WithElectric )
      Aux_Error( ERROR_INFO, "amr->WithElectric is off -> no electric field array is allocated for OPT__FIXUP_ELECTRIC !!\n" );
#  endif


#  pragma omp parallel for
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int s=6; s<18; s++)
      {
         const real *EPtr = amr->patch[0][lv][PID]->electric[s];

//       skip the edges not adjacent to the coarse-fine boundaries
         if ( EPtr == NULL )  continue;

//       set up array indices
         const int e     = s - 6;         // 0 ~ 11 (edge index)
         const int xyz   = (e/4+2)%3;     // (2,2,2,2,0,0,0,0,1,1,1,1): edge direction
         const int LR1   = e%2;           // (0,1,0,1,0,1,0,1,0,1,0,1): left(0)/right(1) edge along the 1st transverse direction
         const int LR2   = e%4/2;         // (0,0,1,1,0,0,1,1,0,0,1,1): left(0)/right(1) egee along the 2nd transverse direction
         const int TDir1 = ( xyz + 1 )%3; // 1st transverse direction
         const int TDir2 = ( xyz + 2 )%3; // 2nd transverse direction

  //     array offsets and strides along the 1st/2nd transverse directions
         int offset1=-1, offset2=-1, stride1=-1, stride2=-1;

         switch ( xyz )
         {
//          Ex --> transverse 1/2 = By/Bz
            case 0 :
               offset1 = LR1*PS1_PS1   + LR2*PS1_PS1M1_PS1P1;
               offset2 = LR1*PS1M1_PS1 + LR2*PS1_PS1_PS1;
               stride1 = 1;
               stride2 = 1;
               break;

//          Ey --> transverse 1/2 = Bz/Bx
            case 1 :
               offset1 = LR1*PS1_PS1_PS1     + LR2*PS1M1;
               offset2 = LR1*PS1_PS1M1_PS1P1 + LR2*PS1;
               stride1 = PS1;
               stride2 = PS1P1;
               break;

//          Ez --> transverse 1/2 = Bx/By
            case 2 :
               offset1 = LR1*PS1   + LR2*PS1M1_PS1P1;
               offset2 = LR1*PS1M1 + LR2*PS1_PS1;
               stride1 = PS1P1_PS1;
               stride2 = PS1P1_PS1;
               break;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "xyz", xyz );
               break;
         } // switch ( xyz )


//       correct B field
         for (int t=0; t<PS1; t++)
         {
            const real CorrB   = EPtr[t]*_dh;
            const real Sign1   = (real)+2.0*(real)LR2 - (real)1.0;
            const real Sign2   = (real)-2.0*(real)LR1 + (real)1.0;
            const int  idx1    = offset1 + t*stride1;
            const int  idx2    = offset2 + t*stride2;
            const int  SibPID1 = amr->patch[0][lv][PID]->sibling[ 2*TDir1 + LR1 ];
            const int  SibPID2 = amr->patch[0][lv][PID]->sibling[ 2*TDir2 + LR2 ];

//          "SibPID == -1" will violate the proper-nesting constraint
//          --> because "EPtr != NULL" only happens around coarse-fine boundaries
#           ifdef GAMER_DEBUG
            if ( SibPID1 == -1 )
               Aux_Error( ERROR_INFO, "SibPID1 == -1 (lv %d, PID %d, sib %d) !!\n", lv, PID, 2*TDir1+LR1 );

            if ( SibPID2 == -1 )
               Aux_Error( ERROR_INFO, "SibPID2 == -1 (lv %d, PID %d, sib %d) !!\n", lv, PID, 2*TDir2+LR2 );
#           endif

//          only need to correct B field on the **coarse-coarse** interfaces
//          --> still need to check if SibPID>=0 for the non-periodic BC
            if ( SibPID1 >= 0  &&  amr->patch[0][lv][SibPID1]->son == -1 )
               amr->patch[MagSg][lv][PID]->magnetic[TDir1][idx1] += Sign1*CorrB;

            if ( SibPID2 >= 0  &&  amr->patch[0][lv][SibPID2]->son == -1 )
               amr->patch[MagSg][lv][PID]->magnetic[TDir2][idx2] += Sign2*CorrB;
         } // for (int t=0; t<PS1; t++)
      } // for (int s=6; s<18; s++)
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : MHD_FixUp_Electric



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
