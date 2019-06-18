#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_FixUp_Electric
// Description :  Use the fine-grid electric field on the coarse-fine boundaries to correct the
//                coarse-grid magnetic field
//
// Note        :  1. Invoked by Flu_FixUp()
//                2  Because the longitudinal (normal) B component on the coarse-fine interfaces has been
//                   corrected by Flu_FixUp_Restrict(), this function only corrects the two **transverse**
//                   B components on such interfaces (i.e., the two B components perpendicular to the normal
//                   vectors of C-F interfaces)
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
//    1. E field on the patch 6 faces
      for (int s=0; s<6; s++)
      {
         const real (*EPtr)[PS1M1_PS1] = ( real (*)[PS1M1_PS1] )amr->patch[0][lv][PID]->electric[s];

//       skip the faces not adjacent to the coarse-fine boundaries
         if ( EPtr == NULL )  continue;


//       1-1. set array indices
         const int xyz = s / 2;           // (0,0,1,1,2,2): face direction
         const int LR  = s % 2;           // (0,1,0,1,0,1): left/right face along xyz
         const int B1  = ( xyz + 2 )%3;   // B component to be fixed by E_((xyz+1)%3)
         const int B2  = ( xyz + 1 )%3;   // B component to be fixed by E_((xyz+2)%3)


//       1-2. set array offsets and strides for B1/B2
         int offset1=-1, offset2=-1, stride1n=-1, stride1m=-1, stride2n=-1, stride2m=-1;

         switch ( xyz )
         {
//          B1: Ey->Bz, (m,n)=(z,y)
//          B2: Ez->By, (m,n)=(z,y)
            case 0 :
               offset1  = LR*PS1M1;
               offset2  = offset1;
               stride1n = PS1;
               stride1m = PS1_PS1;
               stride2n = PS1;
               stride2m = PS1P1_PS1;
               break;

//          B1: Ez->Bx, (m,n)=(x,z)
//          B2: Ex->Bz, (m,n)=(x,z)
            case 1 :
               offset1  = LR*PS1M1_PS1P1;
               offset2  = LR*PS1M1_PS1;
               stride1n = PS1P1_PS1;
               stride1m = 1;
               stride2n = PS1_PS1;
               stride2m = 1;
               break;

//          B1: Ex->By, (m,n)=(y,x)
//          B2: Ey->Bx, (m,n)=(y,x)
            case 2 :
               offset1  = LR*PS1_PS1M1_PS1P1;
               offset2  = offset1;
               stride1n = 1;
               stride1m = PS1;
               stride2n = 1;
               stride2m = PS1P1;
               break;

            default :
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "xyz", xyz );
               break;
         } // switch ( xyz )


//       1-3. correct B field
         const real Coeff = _dh*( (real)2.0*(real)LR - (real)1.0 );  // correction coefficient
         int idx_E, idx_B;

//       B1
         idx_E = 0;
         for (int m=1; m<PS1; m++)  {  idx_B = m*stride1m + 0*stride1n + offset1;
         for (int n=0; n<PS1; n++)  {

            amr->patch[MagSg][lv][PID]->magnetic[B1][idx_B] += Coeff*EPtr[0][idx_E];

            idx_B += stride1n;
            idx_E ++;
         }}

//       B2
         idx_E = 0;
         for (int m=0; m<PS1; m++)  {  idx_B = m*stride2m + 1*stride2n + offset2;
         for (int n=1; n<PS1; n++)  {

            amr->patch[MagSg][lv][PID]->magnetic[B2][idx_B] -= Coeff*EPtr[1][idx_E];

            idx_B += stride2n;
            idx_E ++;
         }}
      } // for (int s=0; s<6; s++)


//    2. E field on the patch 12 edges
      for (int s=6; s<18; s++)
      {
         const real *EPtr = amr->patch[0][lv][PID]->electric[s];

//       skip the edges not adjacent to the coarse-fine boundaries
         if ( EPtr == NULL )  continue;


//       2-1. set array indices
         const int e   = s - 6;           // 0 ~ 11 (edge index)
         const int xyz = (e/4+2)%3;       // (2,2,2,2,0,0,0,0,1,1,1,1): edge direction
         const int B1  = ( xyz + 1 )%3;   // 1st B component to be fixed
         const int B2  = ( xyz + 2 )%3;   // 2nd B component to be fixed
         const int LR1 = e%2;             // (0,1,0,1,0,1,0,1,0,1,0,1): left(0)/right(1) edge along B1
         const int LR2 = e%4/2;           // (0,0,1,1,0,0,1,1,0,0,1,1): left(0)/right(1) egee along B2


//       2-2. set array offsets and strides for B1/B2
         int offset1=-1, offset2=-1, stride1=-1, stride2=-1;

         switch ( xyz )
         {
//          Ex -> B1/B2=By/Bz
            case 0 :
               offset1 = LR1*PS1_PS1   + LR2*PS1_PS1M1_PS1P1;
               offset2 = LR1*PS1M1_PS1 + LR2*PS1_PS1_PS1;
               stride1 = 1;
               stride2 = 1;
               break;

//          Ey -> B1/B2=Bz/Bx
            case 1 :
               offset1 = LR1*PS1_PS1_PS1     + LR2*PS1M1;
               offset2 = LR1*PS1_PS1M1_PS1P1 + LR2*PS1;
               stride1 = PS1;
               stride2 = PS1P1;
               break;

//          Ez -> B1/B2=Bx/By
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


//       2-3. correct B field along B1/2
//       --> only need to correct B field on the (i) coarse-coarse interfaces and (ii) simulation boundaries
//       --> skip **coarse-fine** interfaces since B field on that has been corrected by Flu_Restrict()
         const int  SibPID1 = amr->patch[0][lv][PID]->sibling[ 2*B1 + LR1 ];  // sibling direction along B1/B2
         const int  SibPID2 = amr->patch[0][lv][PID]->sibling[ 2*B2 + LR2 ];
         const real Coeff1  = _dh*( (real)+2.0*(real)LR2 - (real)1.0 );       // correction coefficient for B1/B2
         const real Coeff2  = _dh*( (real)-2.0*(real)LR1 + (real)1.0 );

//       "SibPID == -1" will violate the proper-nesting constraint
//       --> because "EPtr != NULL" only happens around coarse-fine boundaries
#        ifdef GAMER_DEBUG
         if ( SibPID1 == -1 )
            Aux_Error( ERROR_INFO, "SibPID1 == -1 (lv %d, PID %d, sib %d) !!\n", lv, PID, 2*B1+LR1 );

         if ( SibPID2 == -1 )
            Aux_Error( ERROR_INFO, "SibPID2 == -1 (lv %d, PID %d, sib %d) !!\n", lv, PID, 2*B2+LR2 );
#        endif

//       B1
         if (  ( SibPID1 >= 0 && amr->patch[0][lv][SibPID1]->son == -1 )  ||  SibPID1 <= SIB_OFFSET_NONPERIODIC  )
            for (int t=0; t<PS1; t++)
               amr->patch[MagSg][lv][PID]->magnetic[B1][ offset1 + t*stride1 ] += Coeff1*EPtr[t];

//       B2
         if (  ( SibPID2 >= 0 && amr->patch[0][lv][SibPID2]->son == -1 )  ||  SibPID2 <= SIB_OFFSET_NONPERIODIC  )
            for (int t=0; t<PS1; t++)
               amr->patch[MagSg][lv][PID]->magnetic[B2][ offset2 + t*stride2 ] += Coeff2*EPtr[t];
      } // for (int s=6; s<18; s++)
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : MHD_FixUp_Electric



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
