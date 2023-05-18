#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Init_BField_ByVecPot_Function
// Description :  Use the function pointer "Init_BField_ByVecPot_User_Ptr"
//                to assign a magnetic vector potential to all real patches on level
//                "B_lv" and take the curl of this potential to compute the magnetic field
//
// Note        :
//
// Parameter   :  B_lv : Target AMR level
//
// Return      :  amr->patch->magnetic
//-------------------------------------------------------------------------------------------------------
void MHD_Init_BField_ByVecPot_Function( const int B_lv )
{

// check
#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

   if ( Init_BField_ByVecPot_User_Ptr == NULL  &&  OPT__INIT_BFIELD_BYVECPOT == 2 )
      Aux_Error( ERROR_INFO, "Init_BField_ByVecPot_User_Ptr == NULL !!\n" );


   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Constructing the magnetic field from vector potential ...\n" );


// set the number of OpenMP threads
#  ifdef OPENMP
   const int OMP_NThread = ( OPT__INIT_GRID_WITH_OMP ) ? OMP_NTHREAD : 1;
#  else
   const int OMP_NThread = 1;
#  endif


// allocate memory for per-thread arrays
   double **Ax = NULL, **Ay = NULL, **Az = NULL;

   Aux_AllocateArray2D( Ax, OMP_NThread, CUBE(PS1+1) );
   Aux_AllocateArray2D( Ay, OMP_NThread, CUBE(PS1+1) );
   Aux_AllocateArray2D( Az, OMP_NThread, CUBE(PS1+1) );


// loop over the patches on the target AMR level
   const double dh = amr->dh[B_lv];

#  pragma omp parallel for schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)
   {

#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

      const double EdgeL[3] = { amr->patch[0][B_lv][PID]->EdgeL[0],
                                amr->patch[0][B_lv][PID]->EdgeL[1],
                                amr->patch[0][B_lv][PID]->EdgeL[2] };


//    set the magnetic vector potential from function
      for (int k=0; k<PS1+1; k++) {  const double z0 = EdgeL[2] + k*dh;  const double zc = z0 + 0.5*dh;
      for (int j=0; j<PS1+1; j++) {  const double y0 = EdgeL[1] + j*dh;  const double yc = y0 + 0.5*dh;
      for (int i=0; i<PS1+1; i++) {  const double x0 = EdgeL[0] + i*dh;  const double xc = x0 + 0.5*dh;

         int idx = IDX321( i, j, k, PS1+1, PS1+1 );

         Ax[TID][idx] = 0.0;
         Ay[TID][idx] = 0.0;
         Az[TID][idx] = 0.0;

         if ( i != PS1 )   Init_BField_ByVecPot_User_Ptr( &Ax[TID][idx], xc, y0, z0, Time[B_lv], B_lv, 'x', NULL );
         if ( j != PS1 )   Init_BField_ByVecPot_User_Ptr( &Ay[TID][idx], x0, yc, z0, Time[B_lv], B_lv, 'y', NULL );
         if ( k != PS1 )   Init_BField_ByVecPot_User_Ptr( &Az[TID][idx], x0, y0, zc, Time[B_lv], B_lv, 'z', NULL );

      }}}


//    calculate Bx from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1+1; i++) {
         int idx  = IDX321   ( i, j,   k,   PS1+1, PS1+1 );
         int idxj = IDX321   ( i, j+1, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i, j,   k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BX( i, j,   k,   PS1,   PS1   );
         real Bx = ( Az[TID][idxj] - Az[TID][idx] - Ay[TID][idxk] + Ay[TID][idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[0][idxB] = Bx;
      }}}

//    calculate By from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1+1; j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j, k,   PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i,   j, k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BY( i,   j, k,   PS1,   PS1   );
         real By = ( Ax[TID][idxk] - Ax[TID][idx] - Az[TID][idxi] + Az[TID][idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[1][idxB] = By;
      }}}

//    calculate Bz from vector potential
      for (int k=0; k<PS1+1; k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j,   k, PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j,   k, PS1+1, PS1+1 );
         int idxj = IDX321   ( i,   j+1, k, PS1+1, PS1+1 );
         int idxB = IDX321_BZ( i,   j,   k, PS1,   PS1   );
         real Bz = ( Ay[TID][idxi] - Ay[TID][idx] - Ax[TID][idxj] + Ax[TID][idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[2][idxB] = Bz;
      }}}

   } // for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)


// free memory
   Aux_DeallocateArray2D( Ax );
   Aux_DeallocateArray2D( Ay );
   Aux_DeallocateArray2D( Az );


   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Constructing the magnetic field from vector potential ... done\n" );

} // FUNCTION : MHD_Init_BField_ByVecPot_Function
