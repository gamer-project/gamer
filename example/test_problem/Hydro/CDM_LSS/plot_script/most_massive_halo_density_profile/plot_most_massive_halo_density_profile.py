import sys
import argparse
import yt
import matplotlib.pyplot as plt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot profile and out put Halo_parameter' )
parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
args=parser.parse_args()
idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------' )

# load input data frames
yt.enable_parallelism()
ts      = yt.DatasetSeries( [ prefix+'../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
ts0     = ts[0]
storage = {}    # storage dictionary

# "GAMER" physical constants (/include/PhysicalConstant.h) and input parameters
BoxSize_3D               = ts0.domain_width                    # [BoxSize, BoxSize, BoxSize]
BoxSize                  = BoxSize_3D[1]                       # [Mpc/h] = [code_length]
base_level_cell_num_3D   = ts0.domain_dimensions               # [base_level_cell_num, base_level_cell_num, base_level_cell_num]
base_level_cell_num      = base_level_cell_num_3D[1]           # total number of base-level cells
max_AMR_level            = ts0.parameters["MaxLevel"]
input_resolution_limit   = ts0.domain_width.in_units("kpc")[0].d/2**max_AMR_level/ts0.domain_dimensions[0] # [Mpc/h] = [code_length]
length_per_grid          = BoxSize/base_level_cell_num
halo_sphere_radius       = 2.0                                 # [Mpc/h]; expected maximum halo length scale
nbin                     = 128                                 # number of spherical shell bins in the output density profile
print("Simulation BoxSize [Mpc/h]:             "+str(BoxSize_3D))
print("Simulation Base-level Cell Number:      "+str(base_level_cell_num_3D))
print("Smoothed Grid Resolution Limit [Mpc/h]: "+str(length_per_grid))


for sto, ds in ts.piter(storage=storage):
    # cosmological parameters in the current data frame
    Data_ID_Now    = ds.parameters["DumpID"]
    density_filter = ds.all_data()
    current_time_a = ds.current_time                          # scale factor
    current_time_z = (1.0/float(current_time_a)) - 1.0        # redshift

    # [Part I] Identify the Soliton-halo CM Coordinates
    # ------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------

    # CM coordinates for the most massive halo in the input data snapshot
    center_find         = density_filter.quantities.max_location(('gamer','ParDens'))
    # append the correct length unit to the coordinate entries
    if str(center_find[1].units) == "cm":
        CM_actual_coord = ds.arr([center_find[1].d, center_find[2].d, center_find[3].d], 'cm')
    elif str(center_find[1].units) == "code_length":
        CM_actual_coord = ds.arr([center_find[1].d, center_find[2].d, center_find[3].d], 'code_length')
    else:
        print("ERROR CONVERSION ERROR: Check Point #2")
        sys.exit()


    # [Part II] Compute Soliton-halo Density Profile
    # ------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------
    # compute the shell-averaged density and enclosed mass profiles
    sp                   = ds.sphere(CM_actual_coord, (halo_sphere_radius, 'code_length'))
    prof_dens            = yt.create_profile(sp, 'radius', fields= ('gamer','ParDens'), weight_field = 'cell_volume', n_bins = nbin ,
                           units={'radius': 'kpc',('gamer','ParDens'): 'Msun/kpc**3'}, extrema = {'radius': (input_resolution_limit,1e3)})
    radius_raw           = prof_dens.x.value                        # radius grid of the binned spherical shells
    density_raw          = prof_dens[('gamer','ParDens')].value   # local density on the shell at radius

    # remove zero entries in "prof_dens" and "prof_mass_accumulate"
    radius_processed          = []
    density_processed         = []
    for i in range(len(radius_raw)):
        if(density_raw[i]!=0):
            radius_processed.append(radius_raw[i])
            density_processed.append(density_raw[i])
    radius_processed = np.array(radius_processed)
    density_processed = np.array(density_processed)

    # output the halo density profile
    print("Output Density Profile.")
    with open('Data_%06d_Most_Massive_Halo_Density_Profile'%Data_ID_Now,'w') as writer:
        writer.write('Data_%06d, Redshift: %.6f, '%(Data_ID_Now,current_time_z))
        writer.write("CM Coordinates [Mpc/h]: (%.6f, %.6f, %.6f) \n"%(float(CM_actual_coord[0].in_units("code_length")), float(CM_actual_coord[1].in_units("code_length")), float(CM_actual_coord[2].in_units("code_length"))))
        writer.write("---------------------------------------------------------------------------------------------------\n")
        writer.write('Radius [kpc/a]           Density [Msun/kpc^3]     \n')
        for radius_density_pair in zip(radius_processed,density_processed):
            [writer.write(format("%e"%value,"<25")) for value in radius_density_pair]
            writer.write("\n")


    # [Part III] Plot Soliton-halo Density Profile
    # ------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------
    # import the density profiles
    input_halo_density_profile_filename = 'Data_%06d_Most_Massive_Halo_Density_Profile'%Data_ID_Now

    with open("%s"%input_halo_density_profile_filename, "r") as writer:
        header_line    = writer.readline().split(",")
        current_time_z = float(header_line[1].split(" ")[-1])
        coord_1        = ''.join(str(x) for x in header_line[2:]).split(" ")[4:7]
        x_coord_1      = float(coord_1[0][1:])
        y_coord_1      = float(coord_1[1])
        z_coord_1      = float(coord_1[2][:-1])

    file_input_1 = np.genfromtxt("%s"%input_halo_density_profile_filename, skip_header=3)
    radius_1     = file_input_1[:,0]
    density_1    = file_input_1[:,1]

    # generate the density profile plot
    plt.title("Most Massive Halo Density Profile at z = %.3g"%current_time_z)
    plt.loglog(radius_1, density_1, '-o', markersize=2.7, label="Data_%06d[CDM:128$^3$,NX0_TOT:%i$^3$,CM:(%.6g,%.6g,%.6g))]"%(Data_ID_Now,base_level_cell_num,x_coord_1,y_coord_1,z_coord_1))
    plt.xlabel("r [kpc/a]")
    plt.ylabel(r"$\rho_h$ "+"[M$_\odot$/kpc$^3$]")
    plt.legend(prop={'size': 6})
    plt.savefig("Data_%06d_Most_Massive_Halo_Density_Profile.png"%Data_ID_Now, dpi=300)
    plt.clf()