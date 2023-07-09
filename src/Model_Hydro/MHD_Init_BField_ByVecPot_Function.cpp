#include "GAMER.h"

#ifdef MHD


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static double Init_BField_ByVecPot_User_Template( const double x, const double y, const double z, const double Time,
                                                  const int lv, const char Component, double AuxArray[] );

// this function pointer must be set by a test problem initializer
double (*Init_BField_ByVecPot_User_Ptr)( const double x, const double y, const double z, const double Time,
                                         const int lv, const char Component, double AuxArray[] ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_BField_ByVecPot_User_Template
// Description :  Function template to initialize the magnetic vector potential
//
// Note        :  1. Invoked by MHD_Init_BField_ByVecPot_Function() using the function pointer
//                   "Init_BField_ByVecPot_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target refinement level
//                Component : Component of the output magnetic vector potential
//                            --> Supported components: 'x', 'y', 'z'
//                AuxArray  : Auxiliary array
//                            --> Useless since it is currently fixed to NULL
//
// Return      :  "Component"-component of the magnetic vector potential at (x, y, z, Time)
//-------------------------------------------------------------------------------------------------------
double Init_BField_ByVecPot_User_Template( const double x, const double y, const double z, const double Time,
                                           const int lv, const char Component, double AuxArray[] )
{

   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };

   const double x0    = x - BoxCenter[0];
   const double y0    = y - BoxCenter[1];
   const double varpi = 1.0e-4*sqrt( SQR(x0) + SQR(y0) + SQR(amr->dh[MAX_LEVEL]) );

   double mag_vecpot;


   switch ( Component )
   {
      case 'x' :   mag_vecpot = 0.0;           break;
      case 'y' :   mag_vecpot = 0.0;           break;
      case 'z' :   mag_vecpot = 1.0 / varpi;   break;
      default  :   Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
   }


   return mag_vecpot;

} // FUNCTION : Init_BField_ByVecPot_User_Template



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

   if ( OPT__INIT_BFIELD_BYVECPOT != INIT_MAG_BYVECPOT_FUNC )
      Aux_Error( ERROR_INFO, "OPT__INIT_BFIELD_BYVECPOT != %d !!\n", INIT_MAG_BYVECPOT_FUNC );

   if ( Init_BField_ByVecPot_User_Ptr == NULL )
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
   const int    sample_res  = 1 << (MAX_LEVEL - B_lv);
   const double sample_fact = 1.0 / sample_res;
   const double dh          = amr->dh[B_lv];
   const double dh_sample   = dh * sample_fact;

#  pragma omp parallel for schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif


//    set the magnetic vector potential from function
      for (int k=0; k<PS1+1; k++) {  const double z0 = amr->patch[0][B_lv][PID]->EdgeL[2] + k*dh;
      for (int j=0; j<PS1+1; j++) {  const double y0 = amr->patch[0][B_lv][PID]->EdgeL[1] + j*dh;
      for (int i=0; i<PS1+1; i++) {  const double x0 = amr->patch[0][B_lv][PID]->EdgeL[0] + i*dh;

         const int idx = IDX321( i, j, k, PS1+1, PS1+1 );

         Ax[TID][idx] = 0.0;
         Ay[TID][idx] = 0.0;
         Az[TID][idx] = 0.0;


         if ( i != PS1 ) {
            for (int ii=0; ii<sample_res; ii++) {
               const double x = x0 + (ii+0.5)*dh_sample;
               Ax[TID][idx] += Init_BField_ByVecPot_User_Ptr( x,  y0, z0, Time[B_lv], B_lv, 'x', NULL );
            }
         }

         if ( j != PS1 ) {
            for (int jj=0; jj<sample_res; jj++) {
               const double y = y0 + (jj+0.5)*dh_sample;
               Ay[TID][idx] += Init_BField_ByVecPot_User_Ptr( x0, y,  z0, Time[B_lv], B_lv, 'y', NULL );
            }
         }

         if ( k != PS1 ) {
            for (int kk=0; kk<sample_res; kk++) {
               const double z = z0 + (kk+0.5)*dh_sample;
               Az[TID][idx] += Init_BField_ByVecPot_User_Ptr( x0, y0, z,  Time[B_lv], B_lv, 'z', NULL );
            }
         }

         Ax[TID][idx] *= sample_fact;
         Ay[TID][idx] *= sample_fact;
         Az[TID][idx] *= sample_fact;
      }}}


//    calculate Bx from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1+1; i++) {
         const int  idx  = IDX321   ( i, j,   k,   PS1+1, PS1+1 );
         const int  idxj = IDX321   ( i, j+1, k,   PS1+1, PS1+1 );
         const int  idxk = IDX321   ( i, j,   k+1, PS1+1, PS1+1 );
         const int  idxB = IDX321_BX( i, j,   k,   PS1,   PS1   );
         const real Bx   = ( Az[TID][idxj] - Az[TID][idx] - Ay[TID][idxk] + Ay[TID][idx] ) / dh;

         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[MAGX][idxB] = Bx;
      }}}

//    calculate By from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1+1; j++) {
      for (int i=0; i<PS1;   i++) {
         const int  idx  = IDX321   ( i,   j, k,   PS1+1, PS1+1 );
         const int  idxi = IDX321   ( i+1, j, k,   PS1+1, PS1+1 );
         const int  idxk = IDX321   ( i,   j, k+1, PS1+1, PS1+1 );
         const int  idxB = IDX321_BY( i,   j, k,   PS1,   PS1   );
         const real By   = ( Ax[TID][idxk] - Ax[TID][idx] - Az[TID][idxi] + Az[TID][idx] ) / dh;

         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[MAGY][idxB] = By;
      }}}

//    calculate Bz from vector potential
      for (int k=0; k<PS1+1; k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1;   i++) {
         const int  idx  = IDX321   ( i,   j,   k, PS1+1, PS1+1 );
         const int  idxi = IDX321   ( i+1, j,   k, PS1+1, PS1+1 );
         const int  idxj = IDX321   ( i,   j+1, k, PS1+1, PS1+1 );
         const int  idxB = IDX321_BZ( i,   j,   k, PS1,   PS1   );
         const real Bz   = ( Ay[TID][idxi] - Ay[TID][idx] - Ax[TID][idxj] + Ax[TID][idx] ) / dh;

         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[MAGZ][idxB] = Bz;
      }}}

   } // for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)


// free memory
   Aux_DeallocateArray2D( Ax );
   Aux_DeallocateArray2D( Ay );
   Aux_DeallocateArray2D( Az );


   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Constructing the magnetic field from vector potential ... done\n" );

} // FUNCTION : MHD_Init_BField_ByVecPot_Function



#endif // #ifdef MHD
