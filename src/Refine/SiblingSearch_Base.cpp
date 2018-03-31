#include "GAMER.h"

static void FindBaseSibling( const int PID );




//-------------------------------------------------------------------------------------------------------
// Function    :  SiblingSearch_Base
// Description :  Construct the sibling relation for the base level
//
// Note        :  For the base-level sibling relation, it uses the BaseP[] list to locate the sibling patches
//-------------------------------------------------------------------------------------------------------
void SiblingSearch_Base()
{

// initialize all siblings as -1
   for (int PID=0; PID<amr->num[0]; PID++)
   for (int s=0; s<26; s++)
      amr->patch[0][0][PID]->sibling[s] = -1;

// construct the sibling relation for the real patches
   for (int PID=0; PID<amr->NPatchComma[0][1]; PID++)
      FindBaseSibling( PID );

// construct the sibling relation for the buffer patches
#  ifndef SERIAL
   for (int s=0; s<26; s++)
   for (int PID0=amr->NPatchComma[0][s+1]; PID0<amr->NPatchComma[0][s+2]; PID0+=8)
   for (int t=0; t<TABLE_04(s); t++)
      FindBaseSibling( PID0+TABLE_03(s,t) );
#  endif

} // FUNCTION : SiblingSearch_Base



//-------------------------------------------------------------------------------------------------------
// Function    :  FindBaseSibling
// Description :  Construct the sibling relation for the input base-level patch
//
// Note        :  For the base-level sibling relation, it uses the BaseP[] list to locate sibling patches
//
// Parameter   :  PID : Index of the target patch
//-------------------------------------------------------------------------------------------------------
void FindBaseSibling( const int PID )
{

   const int NPatch1D[3] = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };
   const int scale0      = amr->scale[0];

   int order[3], Sib[26];     // order[3] : (i,j,k)th patch in (x,y,z) direction

   for (int d=0; d<3; d++)
      order[d] = (amr->patch[0][0][PID]->corner[d] - MPI_Rank_X[d]*NX0[d]*scale0) / (PATCH_SIZE*scale0) + 2;

   Sib[ 0] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]-1);
   Sib[ 1] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]+1);
   Sib[ 2] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]  );
   Sib[ 3] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]  );
   Sib[ 4] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]  );
   Sib[ 5] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]  );
   Sib[ 6] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]-1);
   Sib[ 7] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]+1);
   Sib[ 8] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]-1);
   Sib[ 9] = (order[2]  )*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);
   Sib[10] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]  );
   Sib[11] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]  );
   Sib[12] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]  );
   Sib[13] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]  );
   Sib[14] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]-1);
   Sib[15] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]-1);
   Sib[16] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]+1);
   Sib[17] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]  )*NPatch1D[0] + (order[0]+1);
   Sib[18] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]-1);
   Sib[19] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]+1);
   Sib[20] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]-1);
   Sib[21] = (order[2]-1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);
   Sib[22] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]-1);
   Sib[23] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]-1)*NPatch1D[0] + (order[0]+1);
   Sib[24] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]-1);
   Sib[25] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);

// record the sibling relation
   for (int s=0; s<26; s++)      amr->patch[0][0][PID]->sibling[s] = BaseP[ Sib[s] ];

} // FUNCTION : FindBaseSibling
