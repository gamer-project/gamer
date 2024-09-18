#include "CUFLU.h"



// external functions and GPU-related set-up
#ifdef __CUDACC__

#if ( MODEL == HYDRO )
#include "CUFLU_Shared_FluUtility.cu"
#endif
#include "CUDA_ConstMemory.h"

#endif // #ifdef __CUDACC__




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/GPU_SrcSolver_IterateAllCells
// Description :  Iterate over all cells to add each source term
//
// Note        :  1. Invoked by CPU_SrcSolver() and CUAPI_Asyn_SrcSolver()
//                2. No ghost zones
//                   --> Should support ghost zones in the future
//
// Parameter   :  g_Flu_Array_In    : Array storing the input fluid variables
//                g_Flu_Array_Out   : Array to store the output fluid variables
//                g_Mag_Array_In    : Array storing the input B field (for MHD only)
//                g_Corner_Array    : Array storing the physical corner coordinates of each patch
//                SrcTerms          : Structure storing all source-term variables
//                NPatchGroup       : Number of patch groups to be evaluated
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                EoS               : EoS object
//
// Return      : fluid[] in all patches
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUSRC_SrcSolver_IterateAllCells(
   const real g_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
         real g_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
   const real g_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
   const double g_Corner_Array[][3],
   const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
   const double TimeNew, const double TimeOld,
   const real MinDens, const real MinPres, const real MinEint, const EoS_t EoS )
#else
void CPU_SrcSolver_IterateAllCells(
   const real g_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
         real g_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
   const real g_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
   const double g_Corner_Array[][3],
   const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
   const double TimeNew, const double TimeOld,
   const real MinDens, const real MinPres, const real MinEint, const EoS_t EoS )
#endif
{

// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//     to work on different patches
#  ifdef __CUDACC__
   const int p = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int p=0; p<8*NPatchGroup; p++)
#  endif
   {
      const double x0 = g_Corner_Array[p][0] - SRC_GHOST_SIZE*dh;
      const double y0 = g_Corner_Array[p][1] - SRC_GHOST_SIZE*dh;
      const double z0 = g_Corner_Array[p][2] - SRC_GHOST_SIZE*dh;

//###REVISE: support ghost zones
      CGPU_LOOP( idx_out, CUBE(PS1) )
      {
//       compute the cell-centered coordinates
         double x, y, z;
         int    i_in, j_in, k_in, idx_in;

         i_in   = SRC_GHOST_SIZE + idx_out % PS1;
         j_in   = SRC_GHOST_SIZE + idx_out % SQR(PS1) / PS1;
         k_in   = SRC_GHOST_SIZE + idx_out / SQR(PS1);
         idx_in = IDX321( i_in, j_in, k_in, SRC_NXT, SRC_NXT );

         x      = x0 + double(i_in*dh);
         y      = y0 + double(j_in*dh);
         z      = z0 + double(k_in*dh);


//       get the input arrays
         real fluid[FLU_NIN_S];

         for (int v=0; v<FLU_NIN_S; v++)  fluid[v] = g_Flu_Array_In[p][v][idx_in];

#        ifdef MHD
         real B[NCOMP_MAG];
         MHD_GetCellCenteredBField( B, g_Mag_Array_In[p][MAGX], g_Mag_Array_In[p][MAGY], g_Mag_Array_In[p][MAGZ],
                                    SRC_NXT, SRC_NXT, SRC_NXT, i_in, j_in, k_in );
#        else
         real *B = NULL;
#        endif


//       add all source terms one by one
//       (1) deleptonization
#        if ( MODEL == HYDRO )
         if ( SrcTerms.Deleptonization )
            SrcTerms.Dlep_FuncPtr( fluid, B, &SrcTerms, dt, dh, x, y, z, TimeNew, TimeOld, MinDens, MinPres, MinEint, &EoS,
                                   SrcTerms.Dlep_AuxArrayDevPtr_Flt, SrcTerms.Dlep_AuxArrayDevPtr_Int );
//       (2) exact cooling
         if ( SrcTerms.ExactCooling )
            SrcTerms.EC_FuncPtr( fluid, B, &SrcTerms, dt, dh, x, y, z, TimeNew, TimeOld, MinDens, MinPres, MinEint, &EoS,
                                 SrcTerms.EC_AuxArrayDevPtr_Flt, SrcTerms.EC_AuxArrayDevPtr_Int );
#        endif

//       (3) user-defined
         if ( SrcTerms.User )
            SrcTerms.User_FuncPtr( fluid, B, &SrcTerms, dt, dh, x, y, z, TimeNew, TimeOld, MinDens, MinPres, MinEint, &EoS,
                                   SrcTerms.User_AuxArrayDevPtr_Flt, SrcTerms.User_AuxArrayDevPtr_Int );

//       store the updated results
         for (int v=0; v<FLU_NOUT_S; v++)   g_Flu_Array_Out[p][v][idx_out] = fluid[v];

      } // CGPU_LOOP( idx_out, CUBE(PS1) )
   } // for (int p=0; p<8*NPatchGroup; p++)

} // FUNCTION : CPU/GPU_SrcSolver_IterateAllCells
