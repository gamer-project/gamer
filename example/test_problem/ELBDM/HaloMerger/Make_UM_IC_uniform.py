import argparse
import numpy as np
import scipy.interpolate

# =============================================================================================================
# Table of Contents
# =============================================================================================================
# Step 01. Load the command-line parameters
# Step 02. Set input parameters
# Step 03. Load the input UM_IC from file
# Step 04. Load the Input__UM_IC_RefineRegion
# Step 05. Set the input UM_IC structure information
# Step 06. Print the input UM_IC information
# Step 07. Set the information of the output single-level UM_IC
# Step 08. Print the output UM_IC information
# Step 09. Define the functions for data construction
# Step 10. Construct the output UM_IC
# Step 11. Write the output UM_IC to the file
# =============================================================================================================

# =============================================================================================================
# Step 01. Load the command-line parameters
parser = argparse.ArgumentParser( description='Make multi-AMR-level UM_IC uniform \n\
                                  Example usage: python3 Make_UM_IC_uniform.py -input UM_IC_input -output UM_IC_output -n0_in 256 -l0_x_in 1.0 -level_out 0')

parser.add_argument( '-input', action='store', required=False, type=str, dest='input_file',
                     help='input file [%(default)s]', default='./UM_IC_input')
parser.add_argument( '-output', action='store', required=False, type=str, dest='output_file',
                     help='output file [%(default)s]', default='./UM_IC_output')
parser.add_argument( '-n0_in', action='store', required=True, type=int, dest='n0_in',
                     help='number of input UM_IC base level grids', nargs='+' )
parser.add_argument( '-l0_x_in', action='store', required=False, type=float, dest='l0_x_in',
                     help='box size of input UM_IC', default=1.0 )
parser.add_argument( '-level_out', action='store', required=False, type=int, dest='level_out',
                     help='target level to output [%(default)d]', default=0 )

args=parser.parse_args()

# =============================================================================================================

# =============================================================================================================
# Step 02. Set input parameters
Target_lv                           = args.level_out              # target level in the input UM_IC to be made uniform to
Input_filename                      = args.input_file             # filename of the input UM_IC to be load
Output_filename                     = args.output_file            # filename of the output UM_IC to be saved
Input__UM_IC_RefineRegion_filename  = 'Input__UM_IC_RefineRegion' # the file to specify the refinement region of the input multi-level UM_IC

UM_IC_Input_BoxSize_x               = args.l0_x_in                # box size of the input UM_IC in the x direction

if len(args.n0_in) == 1:                                          # box is cubic
    UM_IC_Input_N_x_base            = args.n0_in[0]               # number of base-level cells in the input UM_IC in the x direction
    UM_IC_Input_N_y_base            = args.n0_in[0]               # ...                                                  y direction
    UM_IC_Input_N_z_base            = args.n0_in[0]               # ...                                                  z direction
elif len(args.n0_in) == 3:                                        # box is rectangular
    UM_IC_Input_N_x_base            = args.n0_in[0]               # number of base-level cells in the input UM_IC in the x direction
    UM_IC_Input_N_y_base            = args.n0_in[1]               # ...                                                  y direction
    UM_IC_Input_N_z_base            = args.n0_in[2]               # ...                                                  z direction
else:
    raise RuntimeError('n0_in only accepts 1 or 3 paramters !!')

PatchSize                           = 8                           # PATCH_SIZE in GAMER
Float8                              = False                       # whether the UM_IC is in double precision

# Data construction method
Method_Lv_LtoH                      = 2                           # from the lower levels to the target level
                                                                  # (1=repeat, 2=interpolate)
Method_Lv_Same                      = 1                           # from the same level as the target level
                                                                  # (1=paste,  2=interpolate)
Method_Lv_HtoL                      = 1                           # from the higher levels to the target level
                                                                  # (1=pass,   2=interpolate, 3=average)
# Set the data type for UM_IC
if Float8:
    dtype_UM_IC = np.double
else:
    dtype_UM_IC = np.single

# =============================================================================================================

# =============================================================================================================
# Step 03. Load the input UM_IC from file
print( '' )
print( 'Loading input data %s ... '%Input_filename )

UM_IC_Input = np.fromfile( Input_filename, dtype=dtype_UM_IC )

print( 'done!' )

# =============================================================================================================

# =============================================================================================================
# Step 04. Load the Input__UM_IC_RefineRegion
print( '' )
print( 'Loading %s ... '%Input__UM_IC_RefineRegion_filename )

with open( Input__UM_IC_RefineRegion_filename, 'r' ) as f:
    for line in f:
        if line.startswith( '#' ):
            header = line
        elif line.startswith( '\n' ):
            continue
        else:
            break #stop when there are no more #

    f.close()

Input__UM_IC_RefineRegion_header = header[1:].strip().split()
Input__UM_IC_RefineRegion_table  = np.genfromtxt( Input__UM_IC_RefineRegion_filename, delimiter=None, comments='#',
                                                  names=Input__UM_IC_RefineRegion_header, dtype=None, encoding=None )

print( 'done!' )

# =============================================================================================================

# =============================================================================================================
# Step 05. Set the input UM_IC structure information
UM_IC_Input_dh_base       = UM_IC_Input_BoxSize_x/UM_IC_Input_N_x_base       # base-level cell size of the input UM_IC

UM_IC_Input_BoxSize_y     = UM_IC_Input_dh_base*UM_IC_Input_N_y_base         # box size of the input UM_IC in the y direction
UM_IC_Input_BoxSize_z     = UM_IC_Input_dh_base*UM_IC_Input_N_z_base         # ...                                z direction

UM_IC_Input_NLEVEL        = 1+np.max(Input__UM_IC_RefineRegion_table['dLv']) # number of AMR levels of the input UM_IC

# Initialize the arrays
UM_IC_Input_NP_Skip_xL    = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # number of patches on the parent level to be skipped in the x direction from the left edge of the parent refinement region
UM_IC_Input_NP_Skip_xR    = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                                        x ...                right ...
UM_IC_Input_NP_Skip_yL    = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                                        y ...                left  ...
UM_IC_Input_NP_Skip_yR    = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                                        y ...                right ...
UM_IC_Input_NP_Skip_zL    = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                                        z ...                left ...
UM_IC_Input_NP_Skip_zR    = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                                        z ...                right ...
UM_IC_Input_NP_x          = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # number of patches on each level in the x direction
UM_IC_Input_NP_y          = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                    y direction
UM_IC_Input_NP_z          = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                    z direction
UM_IC_Input_N_x           = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # number of cells on each level in the x direction
UM_IC_Input_N_y           = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                  y direction
UM_IC_Input_N_z           = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # ...                                  z direction
UM_IC_Input_index0        = np.zeros( UM_IC_Input_NLEVEL, dtype=np.uint32 )  # starting index of data in the input UM_IC for each level
UM_IC_Input_x0            = np.zeros( UM_IC_Input_NLEVEL )                   # left edge of the refinement region for each level in the x direction
UM_IC_Input_y0            = np.zeros( UM_IC_Input_NLEVEL )                   # left ...                                                 y direction
UM_IC_Input_z0            = np.zeros( UM_IC_Input_NLEVEL )                   # left ...                                                 z direction
UM_IC_Input_x1            = np.zeros( UM_IC_Input_NLEVEL )                   # right ...                                                x direction
UM_IC_Input_y1            = np.zeros( UM_IC_Input_NLEVEL )                   # right ...                                                y direction
UM_IC_Input_z1            = np.zeros( UM_IC_Input_NLEVEL )                   # right ...                                                z direction
UM_IC_Input_dh            = np.zeros( UM_IC_Input_NLEVEL )                   # cell size of the input UM_IC for each level

# Loop for each level to set the input UM_IC structure information
for lv in range( 0, UM_IC_Input_NLEVEL, 1 ):

    UM_IC_Input_NP_Skip_xL[lv] = 0                              if lv == 0 else Input__UM_IC_RefineRegion_table['NP_Skip_xL'][(Input__UM_IC_RefineRegion_table['dLv'] == lv)]
    UM_IC_Input_NP_Skip_xR[lv] = 0                              if lv == 0 else Input__UM_IC_RefineRegion_table['NP_Skip_xR'][(Input__UM_IC_RefineRegion_table['dLv'] == lv)]
    UM_IC_Input_NP_Skip_yL[lv] = 0                              if lv == 0 else Input__UM_IC_RefineRegion_table['NP_Skip_yL'][(Input__UM_IC_RefineRegion_table['dLv'] == lv)]
    UM_IC_Input_NP_Skip_yR[lv] = 0                              if lv == 0 else Input__UM_IC_RefineRegion_table['NP_Skip_yR'][(Input__UM_IC_RefineRegion_table['dLv'] == lv)]
    UM_IC_Input_NP_Skip_zL[lv] = 0                              if lv == 0 else Input__UM_IC_RefineRegion_table['NP_Skip_zL'][(Input__UM_IC_RefineRegion_table['dLv'] == lv)]
    UM_IC_Input_NP_Skip_zR[lv] = 0                              if lv == 0 else Input__UM_IC_RefineRegion_table['NP_Skip_zR'][(Input__UM_IC_RefineRegion_table['dLv'] == lv)]
    UM_IC_Input_NP_x      [lv] = UM_IC_Input_N_x_base/PatchSize if lv == 0 else 2*(UM_IC_Input_NP_x[lv-1]-UM_IC_Input_NP_Skip_xL[lv]-UM_IC_Input_NP_Skip_xR[lv])
    UM_IC_Input_NP_y      [lv] = UM_IC_Input_N_y_base/PatchSize if lv == 0 else 2*(UM_IC_Input_NP_y[lv-1]-UM_IC_Input_NP_Skip_yL[lv]-UM_IC_Input_NP_Skip_yR[lv])
    UM_IC_Input_NP_z      [lv] = UM_IC_Input_N_z_base/PatchSize if lv == 0 else 2*(UM_IC_Input_NP_z[lv-1]-UM_IC_Input_NP_Skip_zL[lv]-UM_IC_Input_NP_Skip_zR[lv])
    UM_IC_Input_N_x       [lv] = UM_IC_Input_N_x_base           if lv == 0 else UM_IC_Input_NP_x[lv]*PatchSize
    UM_IC_Input_N_y       [lv] = UM_IC_Input_N_y_base           if lv == 0 else UM_IC_Input_NP_y[lv]*PatchSize
    UM_IC_Input_N_z       [lv] = UM_IC_Input_N_z_base           if lv == 0 else UM_IC_Input_NP_z[lv]*PatchSize
    UM_IC_Input_index0    [lv] = 0                              if lv == 0 else UM_IC_Input_index0[lv-1]+2*UM_IC_Input_N_z[lv-1]*UM_IC_Input_N_y[lv-1]*UM_IC_Input_N_x[lv-1]
    UM_IC_Input_x0        [lv] = 0.0                            if lv == 0 else UM_IC_Input_x0[lv-1]+UM_IC_Input_NP_Skip_xL[lv]*PatchSize*UM_IC_Input_dh[lv-1]
    UM_IC_Input_y0        [lv] = 0.0                            if lv == 0 else UM_IC_Input_y0[lv-1]+UM_IC_Input_NP_Skip_yL[lv]*PatchSize*UM_IC_Input_dh[lv-1]
    UM_IC_Input_z0        [lv] = 0.0                            if lv == 0 else UM_IC_Input_z0[lv-1]+UM_IC_Input_NP_Skip_zL[lv]*PatchSize*UM_IC_Input_dh[lv-1]
    UM_IC_Input_x1        [lv] = UM_IC_Input_BoxSize_x          if lv == 0 else UM_IC_Input_x1[lv-1]-UM_IC_Input_NP_Skip_xR[lv]*PatchSize*UM_IC_Input_dh[lv-1]
    UM_IC_Input_y1        [lv] = UM_IC_Input_BoxSize_y          if lv == 0 else UM_IC_Input_y1[lv-1]-UM_IC_Input_NP_Skip_yR[lv]*PatchSize*UM_IC_Input_dh[lv-1]
    UM_IC_Input_z1        [lv] = UM_IC_Input_BoxSize_z          if lv == 0 else UM_IC_Input_z1[lv-1]-UM_IC_Input_NP_Skip_zR[lv]*PatchSize*UM_IC_Input_dh[lv-1]
    UM_IC_Input_dh        [lv] = UM_IC_Input_dh_base            if lv == 0 else UM_IC_Input_dh_base/(2**lv)

# =============================================================================================================
# Step 06. Print the input UM_IC information
print( '' )
print( '------------------------------------------------------------------------------------------------' )
print( 'Input UM_IC information' )
print( '------------------------------------------------------------------------------------------------' )
print( f'{PatchSize               = }' )
print( f'{Float8                  = }' )
print( '' )
print( f'{UM_IC_Input_BoxSize_x   = }' )
print( f'{UM_IC_Input_BoxSize_y   = }' )
print( f'{UM_IC_Input_BoxSize_z   = }' )
print( '' )
print( f'{UM_IC_Input_N_x_base    = }' )
print( f'{UM_IC_Input_N_y_base    = }' )
print( f'{UM_IC_Input_N_z_base    = }' )
print( '' )
print( f'{UM_IC_Input_dh_base     = }' )
print( '' )
print( f'{UM_IC_Input_NLEVEL      = }' )
print( '' )
print( f'{UM_IC_Input_NP_Skip_xL  = }' )
print( f'{UM_IC_Input_NP_Skip_xR  = }' )
print( f'{UM_IC_Input_x0          = }' )
print( f'{UM_IC_Input_x1          = }' )
print( '' )
print( f'{UM_IC_Input_NP_Skip_yL  = }' )
print( f'{UM_IC_Input_NP_Skip_yR  = }' )
print( f'{UM_IC_Input_y0          = }' )
print( f'{UM_IC_Input_y1          = }' )
print( '' )
print( f'{UM_IC_Input_NP_Skip_zL  = }' )
print( f'{UM_IC_Input_NP_Skip_zR  = }' )
print( f'{UM_IC_Input_z0          = }' )
print( f'{UM_IC_Input_z1          = }' )
print( '' )
print( f'{UM_IC_Input_NP_x        = }' )
print( f'{UM_IC_Input_NP_y        = }' )
print( f'{UM_IC_Input_NP_z        = }' )
print( '' )
print( f'{UM_IC_Input_N_x         = }' )
print( f'{UM_IC_Input_N_y         = }' )
print( f'{UM_IC_Input_N_z         = }' )
print( '' )
print( f'{UM_IC_Input_dh          = }' )
print( '' )
print( f'{UM_IC_Input_index0      = }' )
print( '------------------------------------------------------------------------------------------------' )
print( '' )

# =============================================================================================================

# =============================================================================================================
# Step 07. Set the information of the output single-level UM_IC
UM_IC_Output_BoxSize_x = UM_IC_Input_BoxSize_x                    # box size of the output UM_IC in the x direction
UM_IC_Output_BoxSize_y = UM_IC_Input_BoxSize_y                    # ...                                 y direction
UM_IC_Output_BoxSize_z = UM_IC_Input_BoxSize_z                    # ...                                 z direction
UM_IC_Output_N_x       = UM_IC_Input_N_x_base*(2**Target_lv)      # number of cells of the output UM_IC in the x direction
UM_IC_Output_N_y       = UM_IC_Input_N_y_base*(2**Target_lv)      # ...                                        y direction
UM_IC_Output_N_z       = UM_IC_Input_N_z_base*(2**Target_lv)      # ...                                        z direction
UM_IC_Output_dh        = UM_IC_Output_BoxSize_x/UM_IC_Output_N_x  # cell size of the output UM_IC

# =============================================================================================================

# =============================================================================================================
# Step 08. Print the output UM_IC information
print( '' )
print( '------------------------------------------------------------------------------------------------' )
print( 'Output UM_IC information' )
print( '------------------------------------------------------------------------------------------------' )
print( f'{Target_lv               = }' )
print( '' )
print( f'{UM_IC_Output_N_x        = }' )
print( f'{UM_IC_Output_N_y        = }' )
print( f'{UM_IC_Output_N_z        = }' )
print( '' )
print( f'{UM_IC_Output_BoxSize_x  = }' )
print( f'{UM_IC_Output_BoxSize_y  = }' )
print( f'{UM_IC_Output_BoxSize_z  = }' )
print( '' )
print( f'{UM_IC_Output_dh         = }' )
print( '' )
print( f'{Method_Lv_LtoH          = }'+' (1=repeat, 2=interpolate           )' )
print( f'{Method_Lv_Same          = }'+' (1=paste,  2=interpolate           )' )
print( f'{Method_Lv_HtoL          = }'+' (1=pass,   2=interpolate, 3=average)' )
print( '------------------------------------------------------------------------------------------------' )
print( '' )

# =============================================================================================================

# =============================================================================================================
# Step 09. Define the functions for data construction

# Use the single level input data as the interpolation table and perform interpolation for the grid points on the output data
def Interpolated_Data( UM_IC_Input_SingleLevel, lv ):

    # difference of levels
    delta_lv                 = lv - Target_lv

    if delta_lv > np.log2( 2*PatchSize ):
        raise RuntimeError( 'Interpolated_Data() does not work when (lv - Target_lv) > log_2( 2*PatchSize ) !!' )

    # Interpolation Table
    InterpolationTable_Real  = UM_IC_Input_SingleLevel[ 0, :, :, : ]
    InterpolationTable_Imag  = UM_IC_Input_SingleLevel[ 1, :, :, : ]

    # Interpolation coordinates
    InterpolationTable_Z     = np.array( [ UM_IC_Input_z0[lv] + (k+0.5)*UM_IC_Input_dh[lv] for k in range( 0, UM_IC_Input_N_z[lv], 1 ) ] )
    InterpolationTable_Y     = np.array( [ UM_IC_Input_y0[lv] + (j+0.5)*UM_IC_Input_dh[lv] for j in range( 0, UM_IC_Input_N_y[lv], 1 ) ] )
    InterpolationTable_X     = np.array( [ UM_IC_Input_x0[lv] + (i+0.5)*UM_IC_Input_dh[lv] for i in range( 0, UM_IC_Input_N_x[lv], 1 ) ] )

    # Interpolator
    Interpolator_Real        = scipy.interpolate.RegularGridInterpolator( (InterpolationTable_Z, InterpolationTable_Y, InterpolationTable_X),
                                                                           InterpolationTable_Real, bounds_error=False, fill_value=None )
    Interpolator_Imag        = scipy.interpolate.RegularGridInterpolator( (InterpolationTable_Z, InterpolationTable_Y, InterpolationTable_X),
                                                                           InterpolationTable_Imag, bounds_error=False, fill_value=None )

    # Interpolated data size
    Interpolated_Data_N_z    = np.around( UM_IC_Input_N_z[lv]/(2**delta_lv) ).astype(int)
    Interpolated_Data_N_y    = np.around( UM_IC_Input_N_y[lv]/(2**delta_lv) ).astype(int)
    Interpolated_Data_N_x    = np.around( UM_IC_Input_N_x[lv]/(2**delta_lv) ).astype(int)

    # Interpolation points
    Interpolation_Pts_Z      = np.array( [ UM_IC_Input_z0[lv] + (k+0.5)*UM_IC_Output_dh for k in range( 0, Interpolated_Data_N_z, 1 ) ], dtype=np.single )
    Interpolation_Pts_Y      = np.array( [ UM_IC_Input_y0[lv] + (j+0.5)*UM_IC_Output_dh for j in range( 0, Interpolated_Data_N_y, 1 ) ], dtype=np.single )
    Interpolation_Pts_X      = np.array( [ UM_IC_Input_x0[lv] + (i+0.5)*UM_IC_Output_dh for i in range( 0, Interpolated_Data_N_x, 1 ) ], dtype=np.single )
    Interpolation_Grid_Z, Interpolation_Grid_Y, Interpolation_Grid_X = np.meshgrid( Interpolation_Pts_Z, Interpolation_Pts_Y, Interpolation_Pts_X, indexing='ij' )
    Interpolation_Points     = np.array( [Interpolation_Grid_Z.ravel(), Interpolation_Grid_Y.ravel(), Interpolation_Grid_X.ravel()], dtype=np.single ).T

    # Interpolation
    Interpolated_Real        = Interpolator_Real( Interpolation_Points ).reshape( (Interpolated_Data_N_z, Interpolated_Data_N_y, Interpolated_Data_N_x) )
    Interpolated_Imag        = Interpolator_Imag( Interpolation_Points ).reshape( (Interpolated_Data_N_z, Interpolated_Data_N_y, Interpolated_Data_N_x) )
    Interpolated_Data        = np.array( [Interpolated_Real, Interpolated_Imag], dtype=dtype_UM_IC )

    return Interpolated_Data

# Repeat the lower level input data to fill the grids of the higher level output data (e.g. 2x2x2 cube -> 4x4x4 cube)
def Repeated_Data( UM_IC_Input_LowerLevel, lv ):

    # difference of levels
    delta_lv                 = lv - Target_lv

    if delta_lv > 0:
        raise RuntimeError( 'lv must <= Target_lv for repeating !!' )

    # Repeating
    Repeated_Data            = np.repeat( np.repeat( np.repeat( UM_IC_Input_LowerLevel, 2**(-delta_lv), axis=3 ), 2**(-delta_lv), axis=2 ), 2**(-delta_lv), axis=1 )

    return Repeated_Data

# Average the higher level input data to fill the grids of the lower level output data (e.g. 2x2x2 cube -> 1x1x1 cube)
def Averaged_Data( UM_IC_Input_HigherLevel, lv ):

    # difference of levels
    delta_lv                 = lv - Target_lv

    if delta_lv < 0:
        raise RuntimeError( 'lv must >= Target_lv for averaging !!' )

    if delta_lv > np.log2( 2*PatchSize ):
        raise RuntimeError( 'Averaged_Data() does not work when (lv - Target_lv) > log_2( 2*PatchSize ) !!' )

    # Averaged data size
    Averaged_Data_N_z        = np.around( UM_IC_Input_N_z[lv]/(2**delta_lv) ).astype(int)
    Averaged_Data_N_y        = np.around( UM_IC_Input_N_y[lv]/(2**delta_lv) ).astype(int)
    Averaged_Data_N_x        = np.around( UM_IC_Input_N_x[lv]/(2**delta_lv) ).astype(int)

    # Averaging
    Averaged_Data            = ( 1.0/(2**delta_lv)**3 )*UM_IC_Input_HigherLevel.reshape( 2, Averaged_Data_N_z, 2**delta_lv,
                                                                                            Averaged_Data_N_y, 2**delta_lv,
                                                                                            Averaged_Data_N_x, 2**delta_lv ).sum(axis=2).sum(axis=3).sum(axis=4)

    return Averaged_Data

# =============================================================================================================

# =============================================================================================================
# Step 10. Construct the output UM_IC
print( '' )
print( 'Constructing the output UM_IC ...' )

UM_IC_Output = np.zeros( ( 2, UM_IC_Output_N_z, UM_IC_Output_N_y, UM_IC_Output_N_x ), dtype=dtype_UM_IC )

for lv in range( 0, UM_IC_Input_NLEVEL, 1 ):

    print( '    lv %d ...'%lv )

    # information for this level
    index0 = UM_IC_Input_index0[lv]  # starting index of data in the UM_IC_Input for this level
    N_z    = UM_IC_Input_N_z[lv]     # number of cells in the z direction for this level
    N_y    = UM_IC_Input_N_y[lv]     # number of cells in the y direction for this level
    N_x    = UM_IC_Input_N_x[lv]     # number of cells in the x direction for this level
    index1 = index0 + 2*N_z*N_y*N_x  # ending index of data in the UM_IC_Input for this level

    # Input UM_IC data for this level
    UM_IC_Input_thislevel = UM_IC_Input[index0:index1].reshape( ( 2, N_z, N_y, N_x ) )

    # Construct the data according to the level and the methods
    if lv < Target_lv:    # Lower levels

        if Method_Lv_LtoH == 1:    # Repeat
            OutputRegion_Data = Repeated_Data( UM_IC_Input_thislevel, lv )

        elif Method_Lv_LtoH == 2:  # Interpolate
            OutputRegion_Data = Interpolated_Data( UM_IC_Input_thislevel, lv )

        else:
            raise RuntimeError( 'Unsported Method_Lv_LtoH !!' )

    elif lv == Target_lv: # Target level

        if Method_Lv_Same == 1:    # Paste
            OutputRegion_Data = UM_IC_Input_thislevel

        elif Method_Lv_Same == 2:  # Interpolate
            OutputRegion_Data = Interpolated_Data( UM_IC_Input_thislevel, lv )

        else:
            raise RuntimeError( 'Unsported Method_Lv_Same !!' )

    elif lv > Target_lv and lv <= Target_lv + np.log2( 2*PatchSize ):  # Higher levels

        if Method_Lv_HtoL == 1:    # Pass
            continue

        elif Method_Lv_HtoL == 2:  # Interpolate
            OutputRegion_Data = Interpolated_Data( UM_IC_Input_thislevel, lv )

        elif Method_Lv_HtoL == 3:  # Average
            OutputRegion_Data = Averaged_Data( UM_IC_Input_thislevel, lv )

        else:
            raise RuntimeError( 'Unsported Method_Lv_HtoL !!' )

    elif lv > Target_lv + np.log2( 2*PatchSize ):  # level is too high that no method can be applied yet
        print( 'levels higher than Target_lv+log_2( 2*PATCH_SIZE ) cannot be handled and will be ignored' )
        break

    else:
        raise RuntimeError('Unknown lv !!')

    # Find the indices of where to put the constructed data
    OutputRegion_index_k0 = np.around( UM_IC_Input_z0[lv]/UM_IC_Output_dh ).astype(int)
    OutputRegion_index_j0 = np.around( UM_IC_Input_y0[lv]/UM_IC_Output_dh ).astype(int)
    OutputRegion_index_i0 = np.around( UM_IC_Input_x0[lv]/UM_IC_Output_dh ).astype(int)

    OutputRegion_index_k1 = OutputRegion_index_k0+OutputRegion_Data.shape[1]
    OutputRegion_index_j1 = OutputRegion_index_j0+OutputRegion_Data.shape[2]
    OutputRegion_index_i1 = OutputRegion_index_i0+OutputRegion_Data.shape[3]

    # Put the constructed data into the output data
    UM_IC_Output[ 0, OutputRegion_index_k0:OutputRegion_index_k1, OutputRegion_index_j0:OutputRegion_index_j1, OutputRegion_index_i0:OutputRegion_index_i1 ] = OutputRegion_Data[ 0, :, :, : ]
    UM_IC_Output[ 1, OutputRegion_index_k0:OutputRegion_index_k1, OutputRegion_index_j0:OutputRegion_index_j1, OutputRegion_index_i0:OutputRegion_index_i1 ] = OutputRegion_Data[ 1, :, :, : ]

print( 'done!' )

# =============================================================================================================

# =============================================================================================================
# Step 11. Write the output UM_IC to the file
print( '' )
print( 'Writing output file %s ...'%Output_filename )
with open( Output_filename, 'wb' ) as f:
    UM_IC_Output.tofile( f )
    f.close()

print( 'done!' )

# =============================================================================================================
