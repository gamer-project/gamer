#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Refine
// Description :  Perform the refine operation during initialization (including the buffer patches)
//
// Note        :  Allocate patches at level "lv+1"
//
// Parameter   :  lv : Target level to be refined (lv --> lv+1)
//-------------------------------------------------------------------------------------------------------
void Init_Refine( const int lv )
{

   if ( lv == NLEVEL-1 )   Aux_Error( ERROR_INFO, "refine the maximum level !!\n" );


   const int Width = PATCH_SIZE*amr->scale[lv+1];
   bool AllocData[8];         // allocate data or not
   int *Cr;

   for (int m=0; m<27; m++)
   {
//    all real patches must store physical data
      if ( m == 0 )
         for (int LocalID=0; LocalID<8; LocalID++)    AllocData[LocalID] = true;

//    only the outer buffer patches do NOT need to store physical data
      else
      {
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            if (  TABLE_01( m-1, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
                  TABLE_01( m-1, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
                  TABLE_01( m-1, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )
               AllocData[LocalID] = false;

            else
               AllocData[LocalID] = true;
         }
      }


      for (int PID=amr->NPatchComma[lv][m]; PID<amr->NPatchComma[lv][m+1]; PID++)
      {
         if ( amr->patch[0][lv][PID]->flag )
         {
//          construct relation : father -> child
            amr->patch[0][lv][PID]->son = amr->num[lv+1];


//          allocate child patches and construct relation : child -> father
            Cr = amr->patch[0][lv][PID]->corner;

            amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2],       PID, AllocData[0], AllocData[0], AllocData[0] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2],       PID, AllocData[1], AllocData[1], AllocData[1] );
            amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2],       PID, AllocData[2], AllocData[2], AllocData[2] );
            amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2]+Width, PID, AllocData[3], AllocData[3], AllocData[3] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2],       PID, AllocData[4], AllocData[4], AllocData[4] );
            amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2]+Width, PID, AllocData[5], AllocData[5], AllocData[5] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2]+Width, PID, AllocData[6], AllocData[6], AllocData[6] );
            amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, PID, AllocData[7], AllocData[7], AllocData[7] );


//          record the number of real/buffer patches along each sibling direction
            amr->NPatchComma[lv+1][m+1] += 8;


//          pass particles from father to son
#           ifdef PARTICLE
            Par_PassParticle2Son_SinglePatch( lv, PID );
#           endif
         } // if ( amr->patch[0][lv][PID]->flag )
      } // for (int PID=amr->NPatchComma[lv][s+1]; PID<amr->NPatchComma[lv][s+2]; PID++)

      for (int n=m+2; n<28; n++)    amr->NPatchComma[lv+1][n] = amr->num[lv+1];

   } // for (int m=0; m<27; m++)


// set up the BounP_IDMap for the level just created
   Buf_RecordBoundaryPatch( lv+1 );

// construct the sibling relation for the level just created (including the buffer patches)
   SiblingSearch( lv+1 );

// get the patch IDs for sending and receiving data between neighboring ranks
   Buf_RecordExchangeDataPatchID( lv+1 );

// allocate flux arrays on level "lv"
   if ( amr->WithFlux )
   Flu_AllocateFluxArray( lv );

// allocate electric arrays on level "lv"
#  ifdef MHD
   if ( amr->WithElectric )
   MHD_AllocateElectricArray( lv );
#  endif

} // FUNCTION : Init_Refine
