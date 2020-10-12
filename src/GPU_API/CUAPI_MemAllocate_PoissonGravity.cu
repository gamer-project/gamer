#include "CUAPI.h"

#if ( defined GPU  &&  defined GRAVITY )



extern real (*d_Rho_Array_P    )[ CUBE(RHO_NXT) ];
extern real (*d_Pot_Array_P_In )[ CUBE(POT_NXT) ];
extern real (*d_Pot_Array_P_Out)[ CUBE(GRA_NXT) ];
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_G)[ CUBE(USG_NXT_G) ];
extern real (*d_Flu_Array_USG_G)[GRA_NIN-1][ CUBE(PS1) ];
#endif
extern real (*d_Flu_Array_G    )[GRA_NIN  ][ CUBE(PS1) ];
extern double (*d_Corner_Array_PGT)[3];
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_G     )[ CUBE(PS1) ];
#endif
#ifdef MHD
extern real (*d_Emag_Array_G   )[ CUBE(PS1) ];
#endif
extern real (*d_Pot_Array_T)    [ CUBE(GRA_NXT) ];
extern real  *d_ExtPotTable;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_PoissonGravity
// Description :  Allocate device and host memory for the Poisson and Gravity solvers
//
// Parameter   :  Pot_NPG  : Number of patch groups evaluated simultaneously by GPU
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate_PoissonGravity( const int Pot_NPG )
{

   const long Pot_NP            = 8*Pot_NPG;
   const long Rho_MemSize_P     = sizeof(real  )*Pot_NP*CUBE(RHO_NXT);
   const long Pot_MemSize_P_In  = sizeof(real  )*Pot_NP*CUBE(POT_NXT);
   const long Pot_MemSize_P_Out = sizeof(real  )*Pot_NP*CUBE(GRA_NXT);
#  ifdef UNSPLIT_GRAVITY
   const long Pot_MemSize_USG_G = sizeof(real  )*Pot_NP*CUBE(USG_NXT_G);
   const long Flu_MemSize_USG_G = sizeof(real  )*Pot_NP*CUBE(PS1)*(GRA_NIN-1);
#  endif
   const long Flu_MemSize_G     = sizeof(real  )*Pot_NP*CUBE(PS1)*(GRA_NIN  );
   const long Corner_MemSize    = sizeof(double)*Pot_NP*3;
#  ifdef DUAL_ENERGY
   const long DE_MemSize_G      = sizeof(char  )*Pot_NP*CUBE(PS1);
#  endif
#  ifdef MHD
   const long Emag_MemSize_G    = sizeof(real  )*Pot_NP*CUBE(PS1);
#  endif
   const long Pot_MemSize_T     = sizeof(real  )*Pot_NP*CUBE(GRA_NXT);
   const long ExtPot_MemSize    = (long)sizeof(real)*EXT_POT_TABLE_NPOINT[0]*EXT_POT_TABLE_NPOINT[1]*EXT_POT_TABLE_NPOINT[2];


// output the total memory requirement
   long TotalSize = Rho_MemSize_P + Pot_MemSize_P_In + Pot_MemSize_P_Out + Flu_MemSize_G + Pot_MemSize_T;
#  ifdef UNSPLIT_GRAVITY
   TotalSize += Pot_MemSize_USG_G + Flu_MemSize_USG_G;
#  endif
#  ifdef DUAL_ENERGY
   TotalSize += DE_MemSize_G;
#  endif
#  ifdef MHD
   TotalSize += Emag_MemSize_G;
#  endif
   if ( OPT__EXT_POT == EXT_POT_TABLE )
   TotalSize += ExtPot_MemSize;

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : total memory requirement in GPU Poisson and gravity solver = %ld MB\n",
                   TotalSize/(1<<20) );


// allocate the device memory
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Rho_Array_P,      Rho_MemSize_P     )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_P_In,   Pot_MemSize_P_In  )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_P_Out,  Pot_MemSize_P_Out )  );
#  ifdef UNSPLIT_GRAVITY
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_USG_G,  Pot_MemSize_USG_G )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_USG_G,  Flu_MemSize_USG_G )  );
#  endif
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_G,      Flu_MemSize_G     )  );

   if ( OPT__EXT_ACC  ||  OPT__EXT_POT )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Corner_Array_PGT, Corner_MemSize    )  );

#  ifdef DUAL_ENERGY
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_DE_Array_G,       DE_MemSize_G      )  );
#  endif

#  ifdef MHD
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Emag_Array_G,     Emag_MemSize_G   )  );
#  endif

   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_T,      Pot_MemSize_T     )  );

   if ( OPT__EXT_POT == EXT_POT_TABLE )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_ExtPotTable,      ExtPot_MemSize    )  );


// allocate the host memory by CUDA
   for (int t=0; t<2; t++)
   {
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Rho_Array_P     [t], Rho_MemSize_P     )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_P_In  [t], Pot_MemSize_P_In  )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_P_Out [t], Pot_MemSize_P_Out )  );
#     ifdef UNSPLIT_GRAVITY
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_USG_G [t], Pot_MemSize_USG_G )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_USG_G [t], Flu_MemSize_USG_G )  );
#     endif
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_G     [t], Flu_MemSize_G     )  );

      if ( OPT__EXT_ACC  ||  OPT__EXT_POT )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Corner_Array_PGT[t], Corner_MemSize    )  );

#     ifdef DUAL_ENERGY
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_DE_Array_G      [t], DE_MemSize_G      )  );
#     endif

#     ifdef MHD
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Emag_Array_G    [t], Emag_MemSize_G    )  );
#     endif

      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_T     [t], Pot_MemSize_T     )  );
   } // for (int t=0; t<2; t++)

   if ( OPT__EXT_POT == EXT_POT_TABLE )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_ExtPotTable,         ExtPot_MemSize    )  );

} // FUNCTION : CUAPI_MemAllocate_PoissonGravity



#endif // #if ( defined GPU  &&  defined GRAVITY )
