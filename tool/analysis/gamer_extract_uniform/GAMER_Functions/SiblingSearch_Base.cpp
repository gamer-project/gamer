#include "ExtractUniform.h"

static void FindBaseSibling( const int PID );




//-------------------------------------------------------------------------------------------------------
// Function    :  SiblingSearch_Base
// Description :  Construct the sibling relation for the base level
//
// Note        :  For the base-level sibling relation, it uses the BaseP list to locate the sibling patches
//-------------------------------------------------------------------------------------------------------
void SiblingSearch_Base()
{

// initialize all siblings as -1
   for (int PID=0; PID<amr.num[0]; PID++)
   for (int s=0; s<26; s++)
      amr.patch[0][PID]->sibling[s] = -1;

// construct the sibling relation for the real patches
   for (int PID=0; PID<NPatchComma[0][1]; PID++)
      FindBaseSibling( PID );

// construct the sibling relation for the buffer patches
   for (int s=0; s<26; s++)
   for (int PID0=NPatchComma[0][s+1]; PID0<NPatchComma[0][s+2]; PID0+=8)
   for (int t=0; t<TABLE_04(s); t++)
      FindBaseSibling( PID0+TABLE_03(s,t) );

}



//-------------------------------------------------------------------------------------------------------
// Function    :  FindBaseSibling
// Description :  Construct the sibling relation for the input base-level patch
//
// Note        :  For the base-level sibling relation, it uses the BaseP list to locate sibling patches
//
// Parameter   :  PID : The index of the targeted patch
//-------------------------------------------------------------------------------------------------------
void FindBaseSibling( const int PID )
{

   const int NPatch1D[3]   = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };
   const int scale         = amr.scale[0];

   int ip, jp, kp, Sib[26];     // (i,j,k)p : (i,j,k)th patch in (x,y,z) direction

   ip = (amr.patch[0][PID]->corner[0] - MyRank_X[0]*NX0[0]*scale) / (PATCH_SIZE*scale) + 2;
   jp = (amr.patch[0][PID]->corner[1] - MyRank_X[1]*NX0[1]*scale) / (PATCH_SIZE*scale) + 2;
   kp = (amr.patch[0][PID]->corner[2] - MyRank_X[2]*NX0[2]*scale) / (PATCH_SIZE*scale) + 2;

   Sib[ 0]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip-1);
   Sib[ 1]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip+1);
   Sib[ 2]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip  );
   Sib[ 3]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip  );
   Sib[ 4]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip  );
   Sib[ 5]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip  );
   Sib[ 6]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip-1);
   Sib[ 7]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip+1);
   Sib[ 8]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip-1);
   Sib[ 9]  =  (kp  ) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip+1);
   Sib[10]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip  );
   Sib[11]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip  );
   Sib[12]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip  );
   Sib[13]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip  );
   Sib[14]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip-1);
   Sib[15]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip-1);
   Sib[16]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip+1);
   Sib[17]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp  ) * NPatch1D[0]  +  (ip+1);
   Sib[18]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip-1);
   Sib[19]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip+1);
   Sib[20]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip-1);
   Sib[21]  =  (kp-1) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip+1);
   Sib[22]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip-1);
   Sib[23]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp-1) * NPatch1D[0]  +  (ip+1);
   Sib[24]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip-1);
   Sib[25]  =  (kp+1) * NPatch1D[1] * NPatch1D[0]  +  (jp+1) * NPatch1D[0]  +  (ip+1);

// record the sibling relation
   for (int s=0; s<26; s++)      amr.patch[0][PID]->sibling[s] = BaseP[ Sib[s] ];

}



