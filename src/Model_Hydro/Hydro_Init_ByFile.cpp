#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static int Nx, Ny, Nz;
static double fluid_xmin, fluid_ymin, fluid_zmin;
static double fluid_dx, fluid_dy, fluid_dz;
static double *fluid_xcoord, *fluid_ycoord, *fluid_zcoord;

#ifdef SUPPORT_HDF5
void ReadField( hid_t hydro_file_id, const int ibegin, const int jbegin,
                const int kbegin, const int iend, const int jend,
                const int kend, double Dens[], double MomX[], double MomY[], double MomZ[], double Engy[] );
#endif

//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Init_BField_ByFile
// Description :  Use the input uniform-mesh array stored in the file "Hydro_Filename"
//                to assign a magnetic vector potential to all real patches on level
//                "lv" and take the curl of this potential to compute the magnetic field
//
// Note        :
//
// Parameter   :  lv         : Target AMR level
//
// Return      :  amr->patch->magnetic
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_ByFile( const int lv )
{

#  ifndef SUPPORT_HDF5
   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be set to load a vector potential from a file !!\n" );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading the hydrodynamics variables from the input file ...\n" );

   const char Hydro_Filename[] = "Hydro_IC";

   const double dh             = amr->dh[lv];

   double *Densf, *MomXf, *MomYf, *MomZf, *Engyf;

   if ( !Aux_CheckFileExist(Hydro_Filename) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Hydro_Filename );


#  ifdef SUPPORT_HDF5

   herr_t status;
   hid_t dataset, dataspace;
   hsize_t dims[3], maxdims[3];
   int ndim;

   hid_t hydro_file_id = H5Fopen(Hydro_Filename, H5F_ACC_RDONLY, H5P_DEFAULT);

   if ( lv == 0 ) {

     dataset = H5Dopen(hydro_file_id, "Dens", H5P_DEFAULT);

     dataspace = H5Dget_space(dataset);

     ndim = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

     if ( ndim != 3 ) Aux_Error( ERROR_INFO, "Incorrect dimensionality of initial consition ndim=%d !!\n", ndim );

     H5Sclose(dataspace);
     H5Dclose(dataset);

//   NOTE: Initial condition arrays are stored in column-major order,
//   i.e., Dens[Nx][Ny][Nz]

     Nx = dims[0];
     Ny = dims[1];
     Nz = dims[2];
   }

   
   //hydro_file_id = H5Fopen(Hydro_Filename, H5F_ACC_RDONLY, H5P_DEFAULT);
   //dataset = H5Dopen(hydro_file_id, "Dens", H5P_DEFAULT);
   //dataspace = H5Dget_space(dataset);
   //double *MomZZ = new double [ Nx*Ny*Nz ];
   //H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MomZZ);

   //for (int i=0;i<Nx*Ny*Nz;i++)
   //{
   //    if (MomZZ[i]!=0.0)
   //    {
   //     printf("MomZZ[i]=%f\n", MomZZ[i]);
   //     exit(0);
   //    }
   //}


#  endif

// Read the coordinate information from the vector potential grid
   fluid_xcoord = new double [ Nx ];
   fluid_ycoord = new double [ Ny ];
   fluid_zcoord = new double [ Nz ];

#  ifdef SUPPORT_HDF5

   dataset = H5Dopen(hydro_file_id, "x", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, fluid_xcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load x-coordinate !!\n" );

   H5Dclose(dataset);

   dataset = H5Dopen(hydro_file_id, "y", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, fluid_ycoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load y-coordinate !!\n" );
   H5Dclose(dataset);

   dataset = H5Dopen(hydro_file_id, "z", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, fluid_zcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load z-coordinate !!\n" );
   H5Dclose(dataset);

#  endif

// Cell spacing
   fluid_dx = fluid_xcoord[1]-fluid_xcoord[0];
   fluid_dy = fluid_ycoord[1]-fluid_ycoord[0];
   fluid_dz = fluid_zcoord[1]-fluid_zcoord[0];

// Left edge of grid
   fluid_xmin = fluid_xcoord[0]-0.5*fluid_dx;
   fluid_ymin = fluid_ycoord[0]-0.5*fluid_dy;
   fluid_zmin = fluid_zcoord[0]-0.5*fluid_dz;

// Right edge of grid
   double fluid_xmax = fluid_xcoord[Nx-1]+0.5*fluid_dx;
   double fluid_ymax = fluid_ycoord[Ny-1]+0.5*fluid_dy;
   double fluid_zmax = fluid_zcoord[Nz-1]+0.5*fluid_dz;

   //printf("fluid_xcoord[0]=%f, fluid_xcoord[1]=%f, fluid_xcoord[127]=%f\n", fluid_xcoord[0], fluid_xcoord[1], fluid_xcoord[127]);
   //printf("fluid_ycoord[0]=%f, fluid_ycoord[1]=%f, fluid_ycoord[127]=%f\n", fluid_ycoord[0], fluid_ycoord[1], fluid_ycoord[127]);
   //printf("fluid_zcoord[0]=%f, fluid_zcoord[1]=%f, fluid_zcoord[127]=%f\n", fluid_zcoord[0], fluid_zcoord[1], fluid_zcoord[127]);
   //printf("fluid_dx=%f, fluid_dy=%f, fluid_dz=%f\n", fluid_dx, fluid_dy, fluid_dz);
   //printf("fluid_xmin=%f, fluid_ymin=%f, fluid_zmin=%f\n", fluid_xmin, fluid_ymin, fluid_zmin);
   //printf("fluid_xmax=%f, fluid_ymax=%f, fluid_zmax=%f\n", fluid_xmax, fluid_ymax, fluid_zmax);


   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++) {

      double EdgeL[3];
      double EdgeR[3];

      for (int i=0; i<3; i++) {
         EdgeL[i] = amr->patch[0][lv][PID]->EdgeL[i];
         EdgeR[i] = amr->patch[0][lv][PID]->EdgeR[i];
      }


      int ibegin = (int)((EdgeL[0]-fluid_xmin)/fluid_dx);
      int jbegin = (int)((EdgeL[1]-fluid_ymin)/fluid_dy);
      int kbegin = (int)((EdgeL[2]-fluid_zmin)/fluid_dz);

      int iend   = (int)((EdgeR[0]-fluid_xmin)/fluid_dx)-1;
      int jend   = (int)((EdgeR[1]-fluid_ymin)/fluid_dy)-1;
      int kend   = (int)((EdgeR[2]-fluid_zmin)/fluid_dz)-1;

      int nlocx  = iend-ibegin+1;
      int nlocy  = jend-jbegin+1;
      int nlocz  = kend-kbegin+1;

      //printf("%03d %03d %03d\n", ibegin, jbegin, kbegin);
      //printf("%03d %03d %03d\n\n", iend, jend, kend);
      int nloc = nlocx*nlocy*nlocz;
      //printf("ibegin=%d, jbegin=%d, kbegin=%d\n", ibegin, jbegin, kbegin);
      //printf("iend=%d,     jend=%d,   kend=%d\n",   iend,   jend,   kend);
      //printf("nlocx=%d, nlocy=%d, nlocz=%d, PID=%d\n", nlocx, nlocy, nlocz, PID);

//    Allocate for the data on the vector potential grid local to this patch and
//    read it from the file

      Densf = new double [nloc];
      MomXf = new double [nloc];
      MomYf = new double [nloc];
      MomZf = new double [nloc];
      Engyf = new double [nloc];

#     ifdef SUPPORT_HDF5
      ReadField( hydro_file_id, ibegin, jbegin, kbegin,
	             iend, jend, kend, Densf, MomXf, MomYf, MomZf, Engyf );
#     endif

      for (int k=0; k<PS1; k++) {
      for (int j=0; j<PS1; j++) {
      for (int i=0; i<PS1; i++) {

         int idx = IDX321( i, j, k, PS1, PS1 ); 

         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] = Densf[idx];
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i] = MomXf[idx];
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i] = MomYf[idx];
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i] = MomZf[idx];
         amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i] = Engyf[idx];

      }}}

      delete [] Densf;
      delete [] MomXf;
      delete [] MomYf;
      delete [] MomZf;
      delete [] Engyf;

   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// Close the magnetic field file

#  ifdef SUPPORT_HDF5
   H5Fclose(hydro_file_id);
#  endif


   delete [] fluid_xcoord;
   delete [] fluid_ycoord;
   delete [] fluid_zcoord;

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Loading the magnetic field from the input file ... done\n" );

} // FUNCTION : Hydro_Init_BField_ByFile


#ifdef SUPPORT_HDF5

void ReadField( hid_t hydro_file_id, const int ibegin, const int jbegin,
		       const int kbegin, const int iend, const int jend,
		       const int kend, double Dens[], double MomX[], double MomY[], double MomZ[], double Engy[]  )
{
   hid_t dataset, dataspace, memspace, dxfer_template;

   herr_t status;

   hsize_t start[3], stride[3], count[3], dims[3], block[3];

   int rank, ierr;

   rank = 3;

   start[0] = kbegin;
   start[1] = jbegin;
   start[2] = ibegin;

   stride[0] = 1;
   stride[1] = 1;
   stride[2] = 1;

   block[0] = 1;
   block[1] = 1;
   block[2] = 1;

   count[0] = kend-kbegin+1;
   count[1] = jend-jbegin+1;
   count[2] = iend-ibegin+1;

   dims[0] = count[0];
   dims[1] = count[1];
   dims[2] = count[2];

// Read Dens
// Opens an existing dataset.
// Returns a dataset identifier if successful; otherwise returns a negative value.
   dataset = H5Dopen(hydro_file_id, "Dens", H5P_DEFAULT);

// Returns an identifier for a copy of the dataspace for a dataset.
// Returns a dataspace identifier if successful; otherwise returns a negative value.
   dataspace = H5Dget_space(dataset);


// Selects a hyperslab region to add to the current selected region. (in-place to replace)
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);

// Creates a new simple dataspace and opens it for access
// Returns a dataspace identifier if successful; otherwise returns a negative value.
   memspace = H5Screate_simple(rank, dims, NULL);

// Reads raw data from a dataset into a buffer.
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, Dens);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load Dens !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read MomX
   dataset = H5Dopen(hydro_file_id, "MomX", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                    H5P_DEFAULT, MomX);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load MomX !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read MomY
   dataset = H5Dopen(hydro_file_id, "MomY", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, MomY);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load MomY !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read MomZ
   dataset = H5Dopen(hydro_file_id, "MomZ", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                    H5P_DEFAULT, MomZ);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load MomZ !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Engy
   dataset = H5Dopen(hydro_file_id, "Engy", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                    H5P_DEFAULT, Engy);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load Engy !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   return;

} // FUNCTION : ReadField

#endif // #ifdef SUPPORT_HDF5
