#include "GAMER.h"
#include "hdf5.h"

typedef double real_par_in;

int Read_Num_Points(string filename)
{

  hid_t   file_id, dataset, dataspace;
  herr_t  status;
  hsize_t dims[1], maxdims[1];

  int rank;

  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset   = H5Dopen(file_id, "/fields/radius", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  H5Fclose(file_id);

  return (int)dims[0];

}

void Read_Profile(string filename, string fieldname, double field[])
{

  hid_t   file_id, dataset;
  herr_t  status;

  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = H5Dopen(file_id, fieldname.c_str(), H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                   H5S_ALL, H5P_DEFAULT, field);
  H5Dclose(dataset);

  H5Fclose(file_id);

  return;

}

long Read_Particle_Number(string filename)
{

  hid_t   file_id, dataset, dataspace;
  herr_t  status;
  hsize_t dims[1], maxdims[1];

  int rank;

  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset   = H5Dopen(file_id, "/dm/particle_mass", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  H5Fclose(file_id);

  return (long)dims[0];

}

void Read_Particles(string filename, long offset, long num, 
                    real_par_in xpos[], real_par_in ypos[], real_par_in zpos[], 
                    real_par_in xvel[], real_par_in yvel[], real_par_in zvel[], 
                    real_par_in mass[])
{

  hid_t   file_id, dataset, dataspace, memspace;
  herr_t  status;

  hsize_t start[2], stride[2], count[2], dims[2], maxdims[2];
  hsize_t start1d[1], stride1d[1], count1d[1], dims1d[1], maxdims1d[1];

  int rank;

  stride[0] = 1;
  stride[1] = 1;
  start[0] = (hsize_t)offset;
  
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset   = H5Dopen(file_id, "/dm/particle_position", H5P_DEFAULT);

  dataspace = H5Dget_space(dataset);

  rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  count[0] = (hsize_t)num;
  count[1] = 1;

  dims[1] = 1;
 
  start[1] = 0;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, xpos);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 1;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, ypos);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 2;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, zpos);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);

  H5Dclose(dataset);

  dataset   = H5Dopen(file_id, "/dm/particle_velocity", H5P_DEFAULT);

  dataspace = H5Dget_space(dataset);

  start[1] = 0;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, xvel);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 1;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, yvel);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 2;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, zvel);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  H5Dclose(dataset);

  dataset   = H5Dopen(file_id, "/dm/particle_mass", H5P_DEFAULT);

  dataspace = H5Dget_space(dataset);

  rank      = H5Sget_simple_extent_dims(dataspace, dims1d, maxdims1d);

  count1d[0] = dims1d[0];
  stride1d[0] = 1;
  start1d[0] = 0;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start1d,
                               stride1d, count1d, NULL);
  memspace = H5Screate_simple(rank, dims1d, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, mass);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  H5Dclose(dataset);

  H5Fclose(file_id);

  return;

}
