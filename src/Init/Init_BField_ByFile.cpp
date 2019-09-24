#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static int nAx, nAy, nAz;
static double Axmin, Aymin, Azmin;
static double Adx, Ady, Adz;
static double *Axcoord, *Aycoord, *Azcoord;


double TSC_Weight( const double x );
void VecPot_ReadField( const int ibegin, const int jbegin, const int kbegin,
                       const int iend, const int jend, const int kend, 
                       double Ax[], double Ay[], double Az[] );

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_VecPot_BField
// Description :  Use the input uniform-mesh array stored in the file "UM_Filename" to assign data to all
//                real patches on level "UM_lv"
//
// Note        :  1. The function pointer Init_ByFile_User_Ptr() points to Init_ByFile_Default() by default
//                   but may be overwritten by various test problem initializers
//
// Parameter   :  UM_Filename  : Target file name
//                UM_lv        : Target AMR level
//                UM_NVar      : Number of variables
//                UM_LoadNRank : Number of parallel I/O
//                UM_Format    : Data format of the file "UM_Filename"
//                               --> UM_IC_FORMAT_VZYX: [field][z][y][x] in a row-major order
//                                   UM_IC_FORMAT_ZYXV: [z][y][x][field] in a row-major order
//
// Return      :  amr->patch->fluid
//-------------------------------------------------------------------------------------------------------
void Init_BField_ByFile( const char B_Filename[], const int B_lv )
{

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  if SUPPORT_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading the magnetic field from the input file ...\n" );

   const double dh       = amr->dh[B_lv];

   double *Axf, *Ayf, *Azf;

   herr_t status;
   hid_t dataset, dataspace;
   hsize_t dims[3], maxdims[3];

   hid_t mag_file_id = H5Fopen(B_Filename, H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset = H5Dopen(file_id, "magnetic_vector_potential_z");

   dataspace = H5Dget_space(dataset);

   H5Sget_simple_extent_dims(dataspace, dims, maxdims);

   H5Sclose(dataspace);
   H5Dclose(dataset);

   nAx = dims[0];
   nAy = dims[1];
   nAz = dims[2];

   Axcoord = new double [ nAx ];
   Aycoord = new double [ nAy ];
   Azcoord = new double [ nAz ];

   dataset = H5Dopen(mag_file_id, "/x");
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
            H5S_ALL, H5P_DEFAULT, Axcoord);
   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "/y");
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, ycoord);
   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "/z");
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, zcoord);
   H5Dclose(dataset);
   
   double *Ax = new double [ CUBE(PS1+1) ];
   double *Ay = new double [ CUBE(PS1+1) ];
   double *Az = new double [ CUBE(PS1+1) ];

   double sample_res = 2**(sim_lrefineMax-B_lv);
   double sample_fact = 1.0/((double)sample_res);

   Adx = Axcoord[1]-Axcoord[0];
   Ady = Aycoord[1]-Aycoord[0];
   Adz = Azcoord[1]-Azcoord[0];

   Axmin = Axcoord[0]-0.5*Adx;
   Aymin = Aycoord[0]-0.5*Ady;
   Azmin = Azcoord[0]-0.5*Adz;

#  pragma omp parallel for schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[VP_lv][1]; PID++)
   {

      int ibegin = max(int((amr->patch[0][lv][PID]->EdgeL[0]-Axmin)/Adx), 0)
      int jbegin = max(int((amr->patch[0][lv][PID]->EdgeL[1]-Aymin)/Ady), 0)
      int kbegin = max(int((amr->patch[0][lv][PID]->EdgeL[2]-Azmin)/Adz), 0)
      
      int iend = min((int)((amr->patch[0][lv][PID]->EdgeR[0]-Axmin)/Adx)+2, nAx-1);
      int jend = min((int)((amr->patch[0][lv][PID]->EdgeR[1]-Aymin)/Ady)+2, nAy-1);
      int kend = min((int)((amr->patch[0][lv][PID]->EdgeR[2]-Azmin)/Adz)+2, nAz-1);
      
      int nlocx = iend-ibegin+1;
      int nlocy = jend-jbegin+1;
      int nlocz = kend-kbegin+1;

      int fdims[3] = { nlocx, nlocy, nlocz };
      int fbegin[3] = { ibegin, jbegin, kbegin };
      int nloc = nlocx*nlocy*nlocz;
      
      Axf = new double [nloc];
      Ayf = new double [nloc];
      Azf = new double [nloc];

      VecPot_ReadField( ibegin, jbegin, kbegin, iend, jend, kend, 
                        Axf, Ayf, Azf );

      for (int k=0; k<PS1+1; k++)    {  const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh;
      for (int j=0; j<PS1+1; j++)    {  const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh;
      for (int i=0; i<PS1+1; i++)    {  const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh;

         int idx = IDX321( i, j, k, PS1+1, PS1+1 );

         for (int ii=0; ii<sample_res; ii++) {  
            const double x = x0 + (ii+0.5)*dh*sample_fact;  
            Ax[idx] += VecPot_Interp( Axf, x, y0, z0, fdims, fbegin );
         }

         for (int jj=0; jj<sample_res; jj++) {  
            const double y = y0 + (jj+0.5)*dh*sample_fact;  
            Ay[idx] += VecPot_Interp( Ayf, x0, y, z0, fdims, fbegin );
         }

         for (int kk=0; kk<sample_res; kk++) {
            const double z = z0 + (kk+0.5)*dh*sample_fact;  
            Az[idx] += VecPot_Interp( Azf, x0, y0, z, fdims, fbegin );
         }

         Ax[idx] *= sample_fact;
         Ay[idx] *= sample_fact;
         Az[idx] *= sample_fact;

      }}}

//    Calculate Bx from vector potential
      for (int k=0; k<PS1; k++)    { 
      for (int j=0; j<PS1; j++)    {  
      for (int i=0; i<PS1+1; i++)  { 
         int idx = IDX321( i, j, k, PS1+1, PS1+1 );
         int idxj = IDX321( i, j+1, k, PS1+1, PS1+1 );
         int idxk = IDX321( i, j, k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BX( i, j, k, PS1, PS1 );
         double Bx = ( Az[idxj] - Az[idx] - Ay[idxk] + Ay[idx] ) / dh;
         amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[0][idxB] = Bx;
      }}}

//    Calculate By from vector potential
      for (int k=0; k<PS1; k++)    { 
      for (int j=0; j<PS1+1; j++)  {  
      for (int i=0; i<PS1; i++)    {
         int idx = IDX321( i, j, k, PS1+1, PS1+1 );
         int idxi = IDX321( i+1, j, k, PS1+1, PS1+1 );
         int idxk = IDX321( i, j, k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BY( i, j, k, PS1, PS1 );
         double By = ( Ax[idxk] - Ax[idx] - Az[idxi] + Az[idx] ) / dh;
         amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[1][idxB] = By;
      }}}

//    Calculate Bz from vector potential
      for (int k=0; k<PS1+1; k++)    { 
      for (int j=0; j<PS1; j++)    {  
      for (int i=0; i<PS1; i++)    { 
         int idx = IDX321( i, j, k, PS1+1, PS1+1 );
         int idxi = IDX321( i+1, j, k, PS1+1, PS1+1 );
         int idxj = IDX321( i, j+1, k, PS1+1, PS1+1 );
         int idxB = IDX321_BY( i, j, k, PS1, PS1 );
         double Bz = ( Ay[idxi] - Ay[idx] - Ax[idxj] + Ax[idx] ) / dh;
         amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[2][idxB] = Bz;
      }}}
   
      delete [] Axf;
      delete [] Ayf;
      delete [] Azf;

   } // for (int PID=0; PID<amr->NPatchComma[VP_lv][1]; PID++)

// Close the magnetic field file
   H5Fclose(mag_file_id);

   delete [] Ax;
   delete [] Ay;
   delete [] Az;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading the magnetic field from the input file ... done\n" );

#  else

   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be set !!\n" );
   
#  endif

} // FUNCTION : Init_BField_ByFile

double VecPot_Interp( const double field[], const double xx, const double yy, const double zz 
                      const int fdims[], const int fbegin[] )
{

   const int ii = (int)((xx-Axmin)/Adx);
   const int jj = (int)((yy-Aymin)/Ady);
   const int kk = (int)((zz-Azmin)/Adz);

   const int ib = ii - fbegin[0];
   const int jb = jj - fbegin[1];
   const int kb = kk - fbegin[2];

   double pot = 0.0;

   if ( ib > 0 && ib < fdims[0]-1 && 
        jb > 0 && jb < fdims[1]-1 && 
        kb > 0 && kb < fdims[2]-1 ) {

      for (int i = -1; i <= 1; i++) { double dx = (xx-Axcoord[ii+i])/Adx;
      for (int j = -1; j <= 1; j++) { double dy = (yy-Aycoord[jj+j])/Ady;
      for (int k = -1; k <= 1; k++) { double dz = (zz-Azcoord[kk+k])/Adz;
         int idx = IDX321( ib+i, jb+j, kb+k, fdims[0], fdims[1] ); 
         pot += field[idx]*TSC_Weight(dx)*TSC_Weight(dy)*TSC_Weight(dz);
      }}}

   }

   return pot;

} // FUNCTION : VecPot_Interp

//-------------------------------------------------------------------------------------------------------
// Function    :  TSC_Weight
// Description :  Function to compute second-order interpolation on a uniform
//                grid to arbitrary points
//
// Parameter   :  x : Scaled coordinate for the TSC kernel
//
// Return      :  weight
//-------------------------------------------------------------------------------------------------------
double TSC_Weight( const double x )
{

   double weight;
   const double xx = fabs(x);

   if ( xx <= 0.5 ) {
      weight = 0.75 - SQR( xx );
   } else if ( xx >= 0.5 && xx <= 1.5 ) {
      weight = 0.5*SQR( 1.5-xx );
   } else {
      weight = 0.0;
   }

   return weight;

} // FUNCTION : TSC_Weight

void VecPot_ReadField( const int ibegin, const int jbegin, const int kbegin,
                       const int iend, const int jend, const int kend, 
                       double Ax[], double Ay[], double Az[] )
{
   hid_t dataset, dataspace, memspace, dxfer_template;

   herr_t status;

   hsize_t start[3], stride[3], count[3], dims[3];

   int rank, ierr;
   
   rank = 3;

   start[0] = ibegin;
   start[1] = jbegin;
   start[2] = kbegin;

   stride[0] = 1;
   stride[1] = 1;
   stride[2] = 1;

   count[0] = iend-ibegin+1;
   count[1] = jend-jbegin+1;
   count[2] = kend-kbegin+1;

   dims[0] = count[0];
   dims[1] = count[1];
   dims[2] = count[2];

// Read Ax
   dataset = H5Dopen(file_id, "magnetic_vector_potential_x");
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Ax);
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);
  
// Read Ay
   dataset = H5Dopen(file_id, "magnetic_vector_potential_y");
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Ay);
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Az
   dataset = H5Dopen(file_id, "magnetic_vector_potential_z");
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Az);
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   return;

} // FUNCTION : VecPot_ReadField
