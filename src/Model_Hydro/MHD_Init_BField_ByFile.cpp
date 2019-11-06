#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static int nAx, nAy, nAz;
static double Axmin, Aymin, Azmin;
static double Adx, Ady, Adz;
static double *Axcoord, *Aycoord, *Azcoord;

double TSC_Weight( const double x );
double VecPot_Interp( const double field[], const double xx, const double yy,
		      const double zz, const int fdims[], const int fbegin[] );
#ifdef SUPPORT_HDF5
void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
                       const int kbegin, const int iend, const int jend,
                       const int kend, double Ax[], double Ay[], double Az[] );
#endif

//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Init_BField_ByFile
// Description :  Use the input uniform-mesh array stored in the file "B_Filename"
//                to assign a magnetic vector potential to all real patches on level
//                "B_lv" and take the curl of this potential to compute the magnetic field
//
// Note        :
//
// Parameter   :  B_lv         : Target AMR level
//
// Return      :  amr->patch->magnetic
//-------------------------------------------------------------------------------------------------------
void MHD_Init_BField_ByFile( const int B_lv )
{

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_HDF5
   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be set to load a vector potential from a file !!\n" );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading the magnetic field from the input file ...\n" );

   const char B_Filename[] = "B_IC";

   const double dh       = amr->dh[B_lv];

   double *Axf, *Ayf, *Azf;

   if ( !Aux_CheckFileExist(B_Filename) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", B_Filename );

// Open the magnetic field file and determine the dimensionality of the vector
// potential grid

#  ifdef SUPPORT_HDF5

   herr_t status;
   hid_t dataset, dataspace;
   hsize_t dims[3], maxdims[3];
   int ndim;

   hid_t mag_file_id = H5Fopen(B_Filename, H5F_ACC_RDONLY, H5P_DEFAULT);

   if ( B_lv == 0 ) {

     dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_z", H5P_DEFAULT);

     dataspace = H5Dget_space(dataset);

     ndim = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

     if ( ndim != 3 ) Aux_Error( ERROR_INFO, "Incorrect dimensionality of vector potential ndim=%d !!\n", ndim );

     H5Sclose(dataspace);
     H5Dclose(dataset);

//   NOTE: Magnetic vector potential arrays are stored in column-major order,
//   i.e., Ax[nAx][nAy][nAz]

     nAx = dims[0];
     nAy = dims[1];
     nAz = dims[2];

   }

#  endif

// Read the coordinate information from the vector potential grid
   Axcoord = new double [ nAx ];
   Aycoord = new double [ nAy ];
   Azcoord = new double [ nAz ];

#  ifdef SUPPORT_HDF5

   dataset = H5Dopen(mag_file_id, "x", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Axcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load x-coordinate !!\n" );

   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "y", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Aycoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load y-coordinate !!\n" );
   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "z", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Azcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load z-coordinate !!\n" );
   H5Dclose(dataset);

#  endif

// Cell spacing and left edge of vector potential grid

   Adx = Axcoord[1]-Axcoord[0];
   Ady = Aycoord[1]-Aycoord[0];
   Adz = Azcoord[1]-Azcoord[0];

   Axmin = Axcoord[0]-0.5*Adx;
   Aymin = Aycoord[0]-0.5*Ady;
   Azmin = Azcoord[0]-0.5*Adz;

   double Axmax = Axcoord[nAx-1]+0.5*Adx;
   double Aymax = Aycoord[nAy-1]+0.5*Ady;
   double Azmax = Azcoord[nAz-1]+0.5*Adz;

   if ( amr->BoxEdgeL[0] < Axmin+2*Adx || amr->BoxEdgeR[0] >= Axmax-2*Adx ||
        amr->BoxEdgeL[1] < Aymin+2*Ady || amr->BoxEdgeR[1] >= Aymax-2*Ady ||
        amr->BoxEdgeL[2] < Azmin+2*Adz || amr->BoxEdgeR[2] >= Azmax-2*Adz )
      Aux_Error( ERROR_INFO, "Input grid is smaller than the simulation domain !!" );

   double *Ax = new double [ CUBE(PS1+1) ];
   double *Ay = new double [ CUBE(PS1+1) ];
   double *Az = new double [ CUBE(PS1+1) ];

   double sample_res = pow(2, MAX_LEVEL-B_lv);
   double sample_fact = 1.0/((double)sample_res);

   for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++) {

      double EdgeL[3];
      double EdgeR[3];

      for (int i=0; i<3; i++) {
         EdgeL[i] = amr->patch[0][B_lv][PID]->EdgeL[i];
         EdgeR[i] = amr->patch[0][B_lv][PID]->EdgeR[i];
      }

//    Compute the beginning and ending indices on the vector potential grid
//    +/- 1 are necessary because we will be computing derivatives

      int ibegin = (int)((EdgeL[0]-Axmin)/Adx)-2;
      int jbegin = (int)((EdgeL[1]-Aymin)/Ady)-2;
      int kbegin = (int)((EdgeL[2]-Azmin)/Adz)-2;

      int iend   = (int)((EdgeR[0]-Axmin)/Adx)+2;
      int jend   = (int)((EdgeR[1]-Aymin)/Ady)+2;
      int kend   = (int)((EdgeR[2]-Azmin)/Adz)+2;

      int nlocx  = iend-ibegin+1;
      int nlocy  = jend-jbegin+1;
      int nlocz  = kend-kbegin+1;

      int fdims[3] = { nlocx, nlocy, nlocz };
      int fbegin[3] = { ibegin, jbegin, kbegin };
      int nloc = nlocx*nlocy*nlocz;

//    Allocate for the data on the vector potential grid local to this patch and
//    read it from the file

      Axf = new double [nloc];
      Ayf = new double [nloc];
      Azf = new double [nloc];

#     ifdef SUPPORT_HDF5
      VecPot_ReadField( mag_file_id, ibegin, jbegin, kbegin,
			               iend, jend, kend, Axf, Ayf, Azf );
#     endif

//    Loop over the indices in this patch and interpolate the vector potential
//    to the current refinement level's resolution

      for (int k=0; k<PS1+1; k++) {  const double z0 = EdgeL[2] + k*dh;
      for (int j=0; j<PS1+1; j++) {  const double y0 = EdgeL[1] + j*dh;
      for (int i=0; i<PS1+1; i++) {  const double x0 = EdgeL[0] + i*dh;

         int idx = IDX321( i, j, k, PS1+1, PS1+1 );

         Ax[idx] = 0.0;
         Ay[idx] = 0.0;
         Az[idx] = 0.0;

         if ( i != PS1 ) {

            for ( int ii=0; ii<sample_res; ii++ ) {
               const double x = x0 + (ii+0.5)*dh*sample_fact;
               Ax[idx] += VecPot_Interp( Axf, x, y0, z0, fdims, fbegin );
            }

         }

         if ( j != PS1 ) {

            for ( int jj=0; jj<sample_res; jj++ ) {
               const double y = y0 + (jj+0.5)*dh*sample_fact;
               Ay[idx] += VecPot_Interp( Ayf, x0, y, z0, fdims, fbegin );
            }

         }

         if ( k != PS1 ) {

            for ( int kk=0; kk<sample_res; kk++ ) {
               const double z = z0 + (kk+0.5)*dh*sample_fact;
               Az[idx] += VecPot_Interp( Azf, x0, y0, z, fdims, fbegin );
            }

         }

         Ax[idx] *= sample_fact;
         Ay[idx] *= sample_fact;
         Az[idx] *= sample_fact;

      }}}

//    Calculate Bx from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1+1; i++) {
         int idx  = IDX321   ( i, j,   k,   PS1+1, PS1+1 );
         int idxj = IDX321   ( i, j+1, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i, j,   k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BX( i, j,   k,   PS1,   PS1   );
         real Bx = ( Az[idxj] - Az[idx] - Ay[idxk] + Ay[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[0][idxB] = Bx;
      }}}

//    Calculate By from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1+1; j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j, k,   PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i,   j, k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BY( i,   j, k,   PS1,   PS1   );
         real By = ( Ax[idxk] - Ax[idx] - Az[idxi] + Az[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[1][idxB] = By;
      }}}

//    Calculate Bz from vector potential
      for (int k=0; k<PS1+1; k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j,   k, PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j,   k, PS1+1, PS1+1 );
         int idxj = IDX321   ( i,   j+1, k, PS1+1, PS1+1 );
         int idxB = IDX321_BZ( i,   j,   k, PS1,   PS1   );
         real Bz = ( Ay[idxi] - Ay[idx] - Ax[idxj] + Ax[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[2][idxB] = Bz;
      }}}

      delete [] Axf;
      delete [] Ayf;
      delete [] Azf;

   } // for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)

// Close the magnetic field file

#  ifdef SUPPORT_HDF5
   H5Fclose(mag_file_id);
#  endif

   delete [] Ax;
   delete [] Ay;
   delete [] Az;

   delete [] Axcoord;
   delete [] Aycoord;
   delete [] Azcoord;

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Loading the magnetic field from the input file ... done\n" );

} // FUNCTION : Hydro_Init_BField_ByFile

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_Interp
// Description :  Use triangle-shaped cloud interpolation to interpolate
//                a vector potential from the input grid to the AMR grid
//
// Parameter   :  field  : Local patch of input vector potential, one component
//                xx     : coordinate along the x-axis
//                yy     : coordinate along the y-axis
//                zz     : coordinate along the z-axis
//                fdims  : size of the input vector potential patch
//                fbegin : index location of the input vector potential patch
//
// Return      :  vector potential component on AMR grid at (xx, yy, zz)
//-------------------------------------------------------------------------------------------------------
double VecPot_Interp( const double field[], const double xx, const double yy,
		                const double zz, const int fdims[], const int fbegin[] )
{

   // Indices into the coordinate vectors
   const int ii = (int)((xx-Axmin)/Adx);
   const int jj = (int)((yy-Aymin)/Ady);
   const int kk = (int)((zz-Azmin)/Adz);

   // Indices into the local vector potential patch
   const int ib = ii - fbegin[0];
   const int jb = jj - fbegin[1];
   const int kb = kk - fbegin[2];

   double pot = 0.0;

   if ( ib == 0 || ib == fdims[0]-1 ||
        jb == 0 || jb == fdims[1]-1 ||
        kb == 0 || kb == fdims[2]-1 ) {
      Aux_Error( ERROR_INFO, "VecPot_Interp: An invalid index was entered!!\n" );
   }

   for (int i = -1; i <= 1; i++) { double dx = (xx-Axcoord[ii+i])/Adx;
   for (int j = -1; j <= 1; j++) { double dy = (yy-Aycoord[jj+j])/Ady;
   for (int k = -1; k <= 1; k++) { double dz = (zz-Azcoord[kk+k])/Adz;
      int idx = (ib+i)*fdims[2]*fdims[1] + (jb+j)*fdims[2] + (kb+k);
      pot += field[idx]*TSC_Weight(dx)*TSC_Weight(dy)*TSC_Weight(dz);
   }}}

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

#ifdef SUPPORT_HDF5

void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
		       const int kbegin, const int iend, const int jend,
		       const int kend, double Ax[], double Ay[], double Az[] )
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
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_x", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Ax);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_x !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Ay
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_y", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                    H5P_DEFAULT, Ay);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_y !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Az
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_z", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Az);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_z !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   return;

} // FUNCTION : VecPot_ReadField

#endif // #ifdef SUPPORT_HDF5
