#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static int nAx, nAy, nAz;
static double Axmin, Aymin, Azmin;
static double Adx, Ady, Adz;
static double *Axcoord, *Aycoord, *Azcoord;

// Grid metadata (coordinate arrays, spacing, lower edge) needed by VecPot_Interp.
// A GridMeta can describe either the root-level coarse grid (Axcoord/Adx/Axmin, etc.)
// or a single refinement patch (BFieldPatch, below) -- VecPot_Interp itself doesn't
// care which.
struct GridMeta
{
   const double *xcoord, *ycoord, *zcoord;
   double        dx, dy, dz;
   double        xmin, ymin, zmin;
};

// A single locally-refined B-field vector-potential patch, as written by
// cluster_generator's RandomClusterField.write_file when refinement_regions is used
// (HDF5 group "patch_00", "patch_01", ...). Optional: files with no such groups
// (old files, or ones written without refinement_regions) behave exactly as before --
// see NumBPatches below.
//
// "children" holds any further-refined patches nested inside this one (HDF5
// groups "child_00", "child_01", ... -- cluster_generator's refinement_regions
// "num_levels" > 1), read/applied recursively by VecPot_ReadPatchGroup/
// VecPot_ApplyPatch; it is NULL (num_children == 0) for a patch with no nested
// children, which is also how old (single-level) patch files behave.
struct BFieldPatch
{
   int         nx, ny, nz;
   double      dx, dy, dz;
   double      xmin, ymin, zmin;
   double      frac_low;
   double     *xcoord, *ycoord, *zcoord;
   double     *Ax, *Ay, *Az;
   double     *window;
   int         num_children;
   BFieldPatch *children;
};

static int         NumBPatches = 0;
static BFieldPatch *BPatches   = NULL;

double TSC_Weight( const double x );
double VecPot_Interp( const double field[], const double xx, const double yy,
                      const double zz, const int fdims[], const int fbegin[],
                      const GridMeta &grid );
bool   VecPot_PointInPatch( const BFieldPatch &patch, const double xx, const double yy, const double zz );
void   VecPot_ApplyPatch( double &pot, const BFieldPatch &patch, const double xx, const double yy,
                          const double zz, const int comp );
double VecPot_EvalComponent( const double coarse_field[], const int fdims[], const int fbegin[],
                             const double xx, const double yy, const double zz, const int comp );
void   VecPot_FreePatch( BFieldPatch &patch );
#ifdef SUPPORT_HDF5
void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
                       const int kbegin, const int iend, const int jend,
                       const int kend, double Ax[], double Ay[], double Az[] );
void VecPot_ReadPatches( hid_t mag_file_id );
void VecPot_ReadPatchGroup( hid_t group_id, BFieldPatch &patch );
#endif

//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Init_BField_ByVecPot_File
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
void MHD_Init_BField_ByVecPot_File( const int B_lv )
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

//   Optionally read any locally-refined vector-potential patches
//   (cluster_generator's RandomClusterField refinement_regions). Old files, or
//   files written without refinement_regions, simply have no "num_patches"
//   attribute -- NumBPatches stays 0 and everything below behaves exactly as
//   it did before this feature existed.
     VecPot_ReadPatches( mag_file_id );

     if ( MPI_Rank == 0 && NumBPatches > 0 )
        Aux_Message( stdout, "   Loading %d B-field refinement patch(es) from \"%s\" ...\n", NumBPatches, B_Filename );

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
               Ax[idx] += VecPot_EvalComponent( Axf, fdims, fbegin, x, y0, z0, 0 );
            }

         }

         if ( j != PS1 ) {

            for ( int jj=0; jj<sample_res; jj++ ) {
               const double y = y0 + (jj+0.5)*dh*sample_fact;
               Ay[idx] += VecPot_EvalComponent( Ayf, fdims, fbegin, x0, y, z0, 1 );
            }

         }

         if ( k != PS1 ) {

            for ( int kk=0; kk<sample_res; kk++ ) {
               const double z = z0 + (kk+0.5)*dh*sample_fact;
               Az[idx] += VecPot_EvalComponent( Azf, fdims, fbegin, x0, y0, z, 2 );
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

// The patch data read once in VecPot_ReadPatches() (guarded by B_lv==0) is reused on
// every level, so it's only freed after the last level has been initialized.
   if ( B_lv == TOP_LEVEL ) {

      for (int p=0; p<NumBPatches; p++)  VecPot_FreePatch( BPatches[p] );

      delete [] BPatches;
      BPatches    = NULL;
      NumBPatches = 0;

   }

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Loading the magnetic field from the input file ... done\n" );

} // FUNCTION : MHD_Init_BField_ByVecPot_File

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
//                grid   : coordinate arrays/spacing/lower-edge of the grid "field"
//                         lives on -- either the root coarse grid or a single
//                         refinement patch (see GridMeta)
//
// Return      :  vector potential component on AMR grid at (xx, yy, zz)
//-------------------------------------------------------------------------------------------------------
double VecPot_Interp( const double field[], const double xx, const double yy,
                                const double zz, const int fdims[], const int fbegin[],
                                const GridMeta &grid )
{

   // Indices into the coordinate vectors
   const int ii = (int)((xx-grid.xmin)/grid.dx);
   const int jj = (int)((yy-grid.ymin)/grid.dy);
   const int kk = (int)((zz-grid.zmin)/grid.dz);

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

   for (int i = -1; i <= 1; i++) { double dx = (xx-grid.xcoord[ii+i])/grid.dx;
   for (int j = -1; j <= 1; j++) { double dy = (yy-grid.ycoord[jj+j])/grid.dy;
   for (int k = -1; k <= 1; k++) { double dz = (zz-grid.zcoord[kk+k])/grid.dz;
      int idx = (ib+i)*fdims[2]*fdims[1] + (jb+j)*fdims[2] + (kb+k);
      pot += field[idx]*TSC_Weight(dx)*TSC_Weight(dy)*TSC_Weight(dz);
   }}}

   return pot;

} // FUNCTION : VecPot_Interp

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_PointInPatch
// Description :  Determine whether (xx,yy,zz) is far enough inside the given refinement
//                patch's array to be safely interpolated by VecPot_Interp (which needs
//                a +/-1 cell margin around the nearest grid index). Unlike the coarse
//                grid (guaranteed by construction to cover the whole simulation domain
//                with margin, so an out-of-range index there is a fatal input error),
//                being outside this margin is the normal case for most query points --
//                refinement patches are small, local sub-boxes -- so this fails
//                gracefully (the caller falls back to the coarse-only value) instead of
//                raising an error.
//
// Parameter   :  patch : The refinement patch to test against
//                xx    : coordinate along the x-axis
//                yy    : coordinate along the y-axis
//                zz    : coordinate along the z-axis
//
// Return      :  true if (xx,yy,zz) can be safely interpolated on "patch"
//-------------------------------------------------------------------------------------------------------
bool VecPot_PointInPatch( const BFieldPatch &patch, const double xx, const double yy, const double zz )
{

   const int ii = (int)((xx-patch.xmin)/patch.dx);
   const int jj = (int)((yy-patch.ymin)/patch.dy);
   const int kk = (int)((zz-patch.zmin)/patch.dz);

   return ( ii >= 1 && ii <= patch.nx-2 &&
            jj >= 1 && jj <= patch.ny-2 &&
            kk >= 1 && kk <= patch.nz-2 );

} // FUNCTION : VecPot_PointInPatch

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_EvalComponent
// Description :  Evaluate one component of the vector potential at (xx,yy,zz),
//                combining the coarse (root-level) grid with any overlapping
//                refinement patch. Outside every patch -- or when there are no patches
//                at all (NumBPatches == 0: old files, or files written without
//                refinement_regions) -- this returns exactly what VecPot_Interp on the
//                coarse grid alone would, so the uniform-grid-only case is unaffected.
//                Inside a patch, the coarse and patch contributions are blended with
//                the same taper-based correction used on the cluster_generator side
//                (RandomClusterField._interpolate_at_points); patches are assumed
//                non-overlapping.
//
// Parameter   :  coarse_field : Local buffer of the coarse grid's vector potential,
//                                one component (as read by VecPot_ReadField)
//                fdims        : size of the coarse_field buffer
//                fbegin       : index location of the coarse_field buffer
//                xx           : coordinate along the x-axis
//                yy           : coordinate along the y-axis
//                zz           : coordinate along the z-axis
//                comp         : which component: 0=x, 1=y, 2=z
//
// Return      :  combined vector potential component at (xx, yy, zz)
//-------------------------------------------------------------------------------------------------------
double VecPot_EvalComponent( const double coarse_field[], const int fdims[], const int fbegin[],
                             const double xx, const double yy, const double zz, const int comp )
{

   const GridMeta coarse_grid = { Axcoord, Aycoord, Azcoord, Adx, Ady, Adz, Axmin, Aymin, Azmin };

   double pot = VecPot_Interp( coarse_field, xx, yy, zz, fdims, fbegin, coarse_grid );

   for (int p=0; p<NumBPatches; p++) {

      if ( VecPot_PointInPatch( BPatches[p], xx, yy, zz ) ) {

         VecPot_ApplyPatch( pot, BPatches[p], xx, yy, zz, comp );
         break;   // top-level patches are assumed non-overlapping

      }

   }

   return pot;

} // FUNCTION : VecPot_EvalComponent

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_ApplyPatch
// Description :  Blend a single refinement patch into "pot" -- the vector potential
//                component already evaluated at (xx,yy,zz) from "patch"'s parent (the
//                coarse grid, or an enclosing, coarser patch) -- then recurse into
//                "patch.children" to blend in any further refinement nested inside it
//                (cluster_generator's refinement_regions "num_levels" > 1). Mirrors
//                cluster_generator's ClusterField._apply_patch exactly, one telescoping
//                chain of patches per cluster (see RandomClusterField.refinement_regions).
//
// Parameter   :  pot   : Vector potential component evaluated at (xx,yy,zz) so far;
//                        updated in place
//                patch : The refinement patch to blend in -- must already be known to
//                        contain (xx,yy,zz) (see VecPot_PointInPatch); not re-checked
//                        here
//                xx    : coordinate along the x-axis
//                yy    : coordinate along the y-axis
//                zz    : coordinate along the z-axis
//                comp  : which component: 0=x, 1=y, 2=z
//
// Return      :  "pot", updated in place
//-------------------------------------------------------------------------------------------------------
void VecPot_ApplyPatch( double &pot, const BFieldPatch &patch, const double xx, const double yy,
                        const double zz, const int comp )
{

   const int      pfdims[3]  = { patch.nx, patch.ny, patch.nz };
   const int      pfbegin[3] = { 0, 0, 0 };
   const GridMeta patch_grid = { patch.xcoord, patch.ycoord, patch.zcoord,
                                  patch.dx, patch.dy, patch.dz,
                                  patch.xmin, patch.ymin, patch.zmin };

   const double  w          = VecPot_Interp( patch.window, xx, yy, zz, pfdims, pfbegin, patch_grid );
   const double  correction = 1.0 - ( 1.0 - sqrt(patch.frac_low) )*w;
   const double *comp_field = ( comp == 0 ) ? patch.Ax : ( comp == 1 ) ? patch.Ay : patch.Az;
   const double  pot_patch  = VecPot_Interp( comp_field, xx, yy, zz, pfdims, pfbegin, patch_grid );

   pot = pot*correction + pot_patch;

   for (int c=0; c<patch.num_children; c++) {

      if ( VecPot_PointInPatch( patch.children[c], xx, yy, zz ) ) {

         VecPot_ApplyPatch( pot, patch.children[c], xx, yy, zz, comp );
         break;   // siblings (children of the same parent) are assumed non-overlapping

      }

   }

} // FUNCTION : VecPot_ApplyPatch

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_FreePatch
// Description :  Recursively free a BFieldPatch's own arrays and every patch nested
//                inside it (patch.children), as populated by VecPot_ReadPatchGroup.
//
// Parameter   :  patch : The patch (and its whole nested chain) to free
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void VecPot_FreePatch( BFieldPatch &patch )
{

   delete [] patch.xcoord;
   delete [] patch.ycoord;
   delete [] patch.zcoord;
   delete [] patch.Ax;
   delete [] patch.Ay;
   delete [] patch.Az;
   delete [] patch.window;

   for (int c=0; c<patch.num_children; c++)  VecPot_FreePatch( patch.children[c] );

   delete [] patch.children;

} // FUNCTION : VecPot_FreePatch

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

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_ReadPatches
// Description :  Detect and read any locally-refined vector-potential patches
//                (HDF5 groups "patch_00", "patch_01", ... as written by
//                cluster_generator's RandomClusterField.write_file when
//                refinement_regions is used), populating the static BPatches array.
//
// Note        :  1. If the root-level "num_patches" attribute is absent -- an old file,
//                   or one written without refinement_regions -- NumBPatches is left at
//                   0 and BPatches at NULL, so VecPot_EvalComponent falls back to
//                   exactly the coarse-grid-only behavior.
//                2. Each patch is small enough (by construction, a local sub-box near a
//                   single cluster) that its full arrays -- not just a local hyperslab,
//                   unlike the coarse grid -- are read whole into every MPI rank.
//                3. Called once (guarded by B_lv==0 at the call site); the data is
//                   reused across all levels and freed only when B_lv==TOP_LEVEL.
//                4. Each top-level patch's own "child_00", "child_01", ... subgroups
//                   (cluster_generator's refinement_regions "num_levels" > 1) are read
//                   recursively by VecPot_ReadPatchGroup.
//
// Parameter   :  mag_file_id : Open HDF5 file identifier for the B-field input file
//
// Return      :  NumBPatches, BPatches (static)
//-------------------------------------------------------------------------------------------------------
void VecPot_ReadPatches( hid_t mag_file_id )
{

   const htri_t has_patches = H5Aexists( mag_file_id, "num_patches" );

   if ( has_patches <= 0 ) {
      NumBPatches = 0;
      BPatches    = NULL;
      return;
   }

   hid_t  attr_id = H5Aopen( mag_file_id, "num_patches", H5P_DEFAULT );
   herr_t status  = H5Aread( attr_id, H5T_NATIVE_INT, &NumBPatches );
   H5Aclose( attr_id );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load the \"num_patches\" attribute !!\n" );

   if ( NumBPatches <= 0 ) {
      NumBPatches = 0;
      BPatches    = NULL;
      return;
   }

   BPatches = new BFieldPatch [ NumBPatches ];

   char group_name[32];

   for (int p=0; p<NumBPatches; p++) {

      snprintf( group_name, sizeof(group_name), "patch_%02d", p );

      hid_t group_id = H5Gopen( mag_file_id, group_name, H5P_DEFAULT );
      if ( group_id < 0 ) Aux_Error( ERROR_INFO, "Failed to open group \"%s\" !!\n", group_name );

      VecPot_ReadPatchGroup( group_id, BPatches[p] );

      H5Gclose( group_id );

   } // for (int p=0; p<NumBPatches; p++)

   return;

} // FUNCTION : VecPot_ReadPatches

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_ReadPatchGroup
// Description :  Read a single vector-potential patch group's own data (coordinates,
//                Ax/Ay/Az, taper window, frac_low) into "patch", then recurse into any
//                "child_00", "child_01", ... subgroups nested inside it
//                (cluster_generator's refinement_regions "num_levels" > 1), populating
//                "patch.children"/"patch.num_children". Used both for each top-level
//                "patch_NN" group (from VecPot_ReadPatches) and, recursively, for every
//                "child_NN" group nested inside one.
//
// Note        :  An absent (or zero) "num_children" attribute -- an old (single-level)
//                patch file, or the finest level of this region's chain -- leaves this
//                patch with no children, exactly as before this feature existed.
//
// Parameter   :  group_id : Open HDF5 group identifier for this patch (or child)
//                patch    : Filled in place
//
// Return      :  "patch", filled in place
//-------------------------------------------------------------------------------------------------------
void VecPot_ReadPatchGroup( hid_t group_id, BFieldPatch &patch )
{

   hid_t   dataset, dataspace;
   hsize_t dims[3], maxdims[3];
   herr_t  status;

   char group_name[256];
   H5Iget_name( group_id, group_name, sizeof(group_name) );

// Dimensions, from the x/y/z coordinate datasets
   dataset   = H5Dopen( group_id, "x", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   H5Sget_simple_extent_dims( dataspace, dims, maxdims );
   patch.nx  = (int)dims[0];
   H5Sclose( dataspace );
   H5Dclose( dataset );

   dataset   = H5Dopen( group_id, "y", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   H5Sget_simple_extent_dims( dataspace, dims, maxdims );
   patch.ny  = (int)dims[0];
   H5Sclose( dataspace );
   H5Dclose( dataset );

   dataset   = H5Dopen( group_id, "z", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   H5Sget_simple_extent_dims( dataspace, dims, maxdims );
   patch.nz  = (int)dims[0];
   H5Sclose( dataspace );
   H5Dclose( dataset );

// Coordinate arrays
   patch.xcoord = new double [ patch.nx ];
   patch.ycoord = new double [ patch.ny ];
   patch.zcoord = new double [ patch.nz ];

   dataset = H5Dopen( group_id, "x", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.xcoord );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/x\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "y", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.ycoord );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/y\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "z", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.zcoord );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/z\" !!\n", group_name );
   H5Dclose( dataset );

   patch.dx = patch.xcoord[1] - patch.xcoord[0];
   patch.dy = patch.ycoord[1] - patch.ycoord[0];
   patch.dz = patch.zcoord[1] - patch.zcoord[0];

   patch.xmin = patch.xcoord[0] - 0.5*patch.dx;
   patch.ymin = patch.ycoord[0] - 0.5*patch.dy;
   patch.zmin = patch.zcoord[0] - 0.5*patch.dz;

// Vector potential components and taper window, read whole (see Note 2 on VecPot_ReadPatches)
   const int npts = patch.nx*patch.ny*patch.nz;

   patch.Ax     = new double [ npts ];
   patch.Ay     = new double [ npts ];
   patch.Az     = new double [ npts ];
   patch.window = new double [ npts ];

   dataset = H5Dopen( group_id, "magnetic_vector_potential_x", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.Ax );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/magnetic_vector_potential_x\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "magnetic_vector_potential_y", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.Ay );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/magnetic_vector_potential_y\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "magnetic_vector_potential_z", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.Az );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/magnetic_vector_potential_z\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "window", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.window );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/window\" !!\n", group_name );
   H5Dclose( dataset );

// frac_low attribute
   hid_t attr = H5Aopen( group_id, "frac_low", H5P_DEFAULT );
   status = H5Aread( attr, H5T_NATIVE_DOUBLE, &patch.frac_low );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/frac_low\" !!\n", group_name );
   H5Aclose( attr );

// Recurse into any further-refined children (see Note above)
   patch.num_children = 0;
   patch.children     = NULL;

   if ( H5Aexists( group_id, "num_children" ) > 0 ) {

      hid_t  nc_attr = H5Aopen( group_id, "num_children", H5P_DEFAULT );
      status = H5Aread( nc_attr, H5T_NATIVE_INT, &patch.num_children );
      H5Aclose( nc_attr );
      if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/num_children\" !!\n", group_name );

   }

   if ( patch.num_children > 0 ) {

      patch.children = new BFieldPatch [ patch.num_children ];

      char child_name[32];

      for (int c=0; c<patch.num_children; c++) {

         snprintf( child_name, sizeof(child_name), "child_%02d", c );

         hid_t child_group = H5Gopen( group_id, child_name, H5P_DEFAULT );
         if ( child_group < 0 ) Aux_Error( ERROR_INFO, "Failed to open group \"%s/%s\" !!\n", group_name, child_name );

         VecPot_ReadPatchGroup( child_group, patch.children[c] );

         H5Gclose( child_group );

      }

   }

} // FUNCTION : VecPot_ReadPatchGroup

#endif // #ifdef SUPPORT_HDF5
