#include "ExtractProfile.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  FindFather
// Description :  Construct the relation between levels "lv" and "lv-1"
//
// Note        :  a. This function assumes that the relations between levels "0, 1, 2, ..., lv-1" have already
//                   been constructed correctly
//                b. Currently this function only works for the function "Init_Reload"
//
// Parameter   :  The targeted refinement level
//-------------------------------------------------------------------------------------------------------
void FindFather( const int lv )
{

   if ( lv <= 0  ||  lv > NLEVEL-1 )
   {
      fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "lv", lv );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit( -1 );
   }


   const int scale0        = amr.scale[0];
   const int PS0           = PATCH_SIZE*scale0;
   const int BaseP_Disp    = 2;
   const int NBaseP[3]     = {  NX0_TOT[0]/PATCH_SIZE + 4,
                                NX0_TOT[1]/PATCH_SIZE + 4,
                                NX0_TOT[2]/PATCH_SIZE + 4  };

   int FaPID, BaseP_ID, BaseP_xyz[3], GrandPaLv, GrandPaPID, FaPS, LocalPos[3], LocalID_1D, LocalID=0;
   int *Corner, *GrandPaCorner;


   for (int PID0=0; PID0<amr.num[lv]; PID0+=8)
   {
      Corner = amr.patch[lv][PID0]->corner;

//    a. find the ancestor patch at the base level
      for (int d=0; d<3; d++)    BaseP_xyz[d] = Corner[d] / PS0  +  BaseP_Disp;

      BaseP_ID = BaseP_xyz[2]*NBaseP[1]*NBaseP[0] + BaseP_xyz[1]*NBaseP[0] + BaseP_xyz[0];
      FaPID    = BaseP[ BaseP_ID ];


//    b. find the father patch
      for (int FaLv=1; FaLv<lv; FaLv++)
      {
         GrandPaLv      = FaLv - 1;
         GrandPaPID     = FaPID;
         GrandPaCorner  = amr.patch[GrandPaLv][GrandPaPID]->corner;
         FaPS           = PATCH_SIZE*amr.scale[FaLv];

         for (int d=0; d<3; d++)
            LocalPos[d] = ( Corner[d] - GrandPaCorner[d] ) / FaPS;

         LocalID_1D     = LocalPos[2]*4 + LocalPos[1]*2 + LocalPos[0];

         switch ( LocalID_1D )
         {
            case 0:  LocalID = 0;   break;
            case 1:  LocalID = 1;   break;
            case 2:  LocalID = 2;   break;
            case 3:  LocalID = 4;   break;
            case 4:  LocalID = 3;   break;
            case 5:  LocalID = 6;   break;
            case 6:  LocalID = 5;   break;
            case 7:  LocalID = 7;   break;
            default:
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "LocalID_1D", LocalID_1D );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n",
                        __FILE__, __LINE__,  __FUNCTION__  );
               exit( -1 );
         }

         FaPID = amr.patch[GrandPaLv][GrandPaPID]->son + LocalID;

      } // for (int FaLv=1; FaLv<lv; FaLv++)


//    c. store the relation
      amr.patch[lv-1][FaPID]->son = PID0;

      for (int PID=PID0; PID<PID0+8; PID++)        amr.patch[lv][PID]->father = FaPID;

   } // for (int PID0=0; PID0<amr.num[lv]; PID0+=8)

} // FUNCTION : FindFather



