import argparse
import sys
import yt
import math
import matplotlib.image
import numpy as np
import matplotlib as mpl
import os
from matplotlib import pyplot as plt
from numpy import loadtxt
from scipy.optimize import fsolve
mpl.rcParams['agg.path.chunksize'] = 10000

# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters (consistent with "Input__TestProb")
Gas_Par_Setup    = 1             # Gas-only [1] or Particle-only [2] setup
n_Pert_Wave_Len  = 1             # Pert_Wave_Len = BOXSIZE_x/n_Pert_Wave_Len
zc_Collapse      = 1
Temp_gas_init    = 100           # [Kelvin]

# -------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Slice of mass density' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]) )
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

# -------------------------------------------------------------------------------------------------------------------------
# set up basic parameters
colormap_dens    = 'algae'
center_mode      = 'c'
dpi              = 150
if (Gas_Par_Setup == 1): # gas-only setup
    field_dens       = ("gas", "density")
    field_xvel       = ('gas', 'velocity_x')
    field_temp       = ('gas', 'temperature')
    field_pressure   = ('gas', 'pressure')
    field_sound_spd  = ('gas', 'sound_speed')
    Gas_Par_Setup_T  = "Gas-only"
elif (Gas_Par_Setup == 2): # particle-only setup
    field_dens       = ("all", "particle_mass")
    field_xvel       = ('all', 'particle_velocity_x')
    Gas_Par_Setup_T  = "Particle-only"
    Analysis_Par_Bin = 20
else:
    print("Error: Gas_Par_Setup != (1) or (2)")
    sys.exit()

# -------------------------------------------------------------------------------------------------------------------------
# load simulation data with yt
yt.enable_parallelism()
ts                     = yt.DatasetSeries( [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
ts0                    = ts[0]
base_level_cell_num    = int(ts0.domain_dimensions[1])
BoxSize_x              = float(ts0.domain_width[0])        # [Mpc/h] = [code_length]
A_Init                 = ts0.parameters["A_Init"]
Z_Init                 = (1.0/float(A_Init)) - 1.0
Unit_E                 = ts0.parameters["Unit_E"]
Unit_M                 = ts0.parameters["Unit_M"]
Pert_Wave_Len          = BoxSize_x/n_Pert_Wave_Len
if (Gas_Par_Setup == 1): # gas-only setup
    Molecular_Weight       = ts0.parameters["MolecularWeight"]
    base_level_cell_num_3D = ts0.domain_dimensions             # [base_level_cell_num, base_level_cell_num, base_level_cell_num]
    base_level_mid_index   = base_level_cell_num_3D/2          # mid-point (x,y,z) grid indices
    Const_amu               = 1.660539040e-24                  # in cgs
    Const_kB                = 1.38064852e-16                   # in cgs
    # import data columns from GAMER "Output_L1Error()" output files
    def Output_L1Error_Read_In(input_filename):
        data                 = loadtxt("%s"%input_filename, skiprows=1, dtype=float)
        Input_Coordinates    = data[:,0].astype(float)
        Input_Numerical      = data[:,1].astype(float)
        Input_Analytical     = data[:,2].astype(float)
        Input_Absolute_Error = data[:,3].astype(float)
        Input_Relative_Error = Input_Absolute_Error/abs(Input_Analytical)
        return [Input_Coordinates,Input_Numerical,Input_Analytical,Input_Relative_Error]
    # combine two subplots vertically
    def Subplots_Vertical_Combine(input_filename_1,input_filename_2,output_filename):
        subplot_merge_1 = matplotlib.image.imread(input_filename_1)
        subplot_merge_2 = matplotlib.image.imread(input_filename_2)
        merged_image    = np.concatenate((subplot_merge_1, subplot_merge_2))
        matplotlib.image.imsave(output_filename, merged_image)
        os.remove("./"+input_filename_1)
        os.remove("./"+input_filename_2)
elif (Gas_Par_Setup == 2): # particle-only setup
    BoxSize_y              = float(ts0.domain_width[1])        # [Mpc/h] = [code_length]
    BoxSize_z              = float(ts0.domain_width[2])        # [Mpc/h] = [code_length]
    x_Len_per_Bin          = BoxSize_x/Analysis_Par_Bin
    AveDens_Init           = ts0.parameters["AveDens_Init"]
    BoxScale               = ts0.parameters["BoxScale"]
    Unit_T                 = ts0.parameters["Unit_T"]
    Hubble0                = ts0.parameters["Hubble0"]
    Hubble0_Units          = Hubble0*(100*3.2407792894458e-20)*Unit_T

# -------------------------------------------------------------------------------------------------------------------------
# data analysis with yt
for ds in ts.piter():
    # load data frame cosmological parameters
    DumpID_now       = ds.parameters["DumpID"]
    current_time_a   = ds.current_time                          # scale factor
    current_time_z   = (1.0/float(current_time_a)) - 1.0        # redshift

    # -------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------
    if (Gas_Par_Setup == 1): # gas-only setup
        if current_time_z > zc_Collapse:
            # load numerical errors in perturbed density/x-momentum/pressure profiles
            [x_line_Coordinates,Dens_Numerical,Dens_Analytical,Dens_Relative_Error] = Output_L1Error_Read_In("../Zeldovich_Dens_%06d"%DumpID_now)
            [x_line_Coordinates,MomX_Numerical,MomX_Analytical,MomX_Relative_Error] = Output_L1Error_Read_In("../Zeldovich_MomX_%06d"%DumpID_now)
            [x_line_Coordinates,Pres_Numerical,Pres_Analytical,Pres_Relative_Error] = Output_L1Error_Read_In("../Zeldovich_Pres_%06d"%DumpID_now)
            [x_line_Coordinates,Temp_Numerical,Temp_Analytical,Temp_Relative_Error] = Output_L1Error_Read_In("../Zeldovich_Temp_%06d"%DumpID_now)
            Temp_Numerical = Temp_Numerical*(Molecular_Weight*Const_amu/Const_kB*(Unit_E/Unit_M)/(current_time_a**2))            # [Kelvin]
            # Zeldovich pancake collapse temperature profile analytical solution at current_time_z
            def Zeldovich_Temp_Analytical(density_profile):
                return Temp_gas_init*np.abs(((1+current_time_z)/(1+Z_Init))**3*density_profile)**(2.0/3.0)

        # extract on-grid physical quantities
        smoothed_grid           = ds.smoothed_covering_grid(0,[0.0, 0.0, 0.0], dims=base_level_cell_num_3D)
        smoothed_grid_den       = smoothed_grid[field_dens].in_units("code_mass/code_length**3")
        smoothed_grid_temp      = smoothed_grid[field_temp]/(current_time_a**2)              # [Kelvin]
        smoothed_grid_sspd      = smoothed_grid[field_sound_spd].in_units("km/s")
        smoothed_grid_x         = smoothed_grid["gamer","x"].in_units("code_length")
        smoothed_grid_x_mid     = smoothed_grid_x[-1,base_level_mid_index[1],base_level_mid_index[2]]/2
        smoothed_grid_den_line  = np.array(smoothed_grid_den[:,base_level_mid_index[1],base_level_mid_index[2]])
        smoothed_grid_den_line  = list(smoothed_grid_den_line/np.mean(smoothed_grid_den_line))
        smoothed_grid_temp_line = list(np.array(smoothed_grid_temp[:,base_level_mid_index[1],base_level_mid_index[2]]))
        smoothed_grid_sspd_line = list(np.array(smoothed_grid_sspd[:,base_level_mid_index[1],base_level_mid_index[2]]))
        smoothed_grid_x_line    = list(np.array(smoothed_grid_x[:,base_level_mid_index[1],base_level_mid_index[2]] - smoothed_grid_x_mid))
        mean_cov_den            = np.mean(smoothed_grid[field_dens].in_units("Msun/kpc**3")[:,base_level_mid_index[1],base_level_mid_index[2]])
        if current_time_z <= zc_Collapse:
            smoothed_grid_vel       = smoothed_grid[field_xvel].in_units("km/s")/current_time_a
            smoothed_grid_vel_line  = list(np.array(smoothed_grid_vel[:,base_level_mid_index[1],base_level_mid_index[2]]))
            smoothed_grid_pres      = smoothed_grid[field_pressure].in_units("code_pressure")
            smoothed_grid_pres_line = list(np.array(smoothed_grid_pres[:,base_level_mid_index[1],base_level_mid_index[2]]))

        #-----------------------------------------------------------------------------------------------------------------
        # [yt] xy-plane slice plot
        px_dens_sliceOG = yt.SlicePlot(ds, 'z', field_dens, center = center_mode)
        px_dens_sliceOG.set_unit(field_dens, 'Msun/kpc**3')
        px_dens_sliceOG.set_axes_unit("Mpc/h")
        px_dens_sliceOG.set_cmap(field_dens, colormap_dens)
        px_dens_sliceOG.annotate_grids()
        px_dens_sliceOG.annotate_title(r"Slice Comoving Density: $z_\mathrm{c}$ = %.2f, $\rho_\mathrm{mean,comoving}$ = %.2e [Msun/kpc$^3$]"%(zc_Collapse, mean_cov_den))
        px_dens_sliceOG.annotate_timestamp(corner='upper_right', redshift=True, time=False, text_args={'color':'k'})
        px_dens_sliceOG.save("Data_%06d_Density_Slice.png"%DumpID_now, mpl_kwargs={"dpi":dpi})

        # [yt] xy-plane slice plot
        px_vel_sliceOG = yt.SlicePlot(ds, 'z', field_xvel, center = center_mode)
        px_vel_sliceOG.set_unit(field_xvel, 'km/s')
        px_vel_sliceOG.set_axes_unit("Mpc/h")
        px_vel_sliceOG.set_zlim(field_xvel, -2.3e3, 2.3e3)
        px_vel_sliceOG.set_cmap(field_xvel, colormap_dens)
        px_vel_sliceOG.annotate_grids()
        px_vel_sliceOG.annotate_title("Velocity: %s setup, $\lambda_\mathrm{pert} =$ %.1f [Mpc/h], $z_\mathrm{c}$ = %.2f"%(Gas_Par_Setup_T,Pert_Wave_Len,zc_Collapse))
        px_vel_sliceOG.annotate_timestamp(corner='upper_right', redshift=True, time=False, text_args={'color':'k'})
        px_vel_sliceOG.save("Data_%06d_Velocity_Slice.png"%DumpID_now, mpl_kwargs={"dpi":dpi})

        # merge subplots
        Subplots_Vertical_Combine("Data_%06d_Density_Slice.png"%DumpID_now,"Data_%06d_Velocity_Slice.png"%DumpID_now
                                ,"Data_%06d_Density_Velocity_Slice.png"%DumpID_now)

        #-----------------------------------------------------------------------------------------------------------------
        # x-line density plot
        plt.axhline(y=1, color='0.9', linestyle='--')
        if current_time_z > zc_Collapse:
            plt.plot(x_line_Coordinates - BoxSize_x/2, Dens_Numerical, label = "GAMER")
            plt.plot(x_line_Coordinates - BoxSize_x/2, Dens_Analytical, color="r", linestyle='-.', label = "Analytical")
        else:
            px_dens_line = plt.plot(smoothed_grid_x_line,smoothed_grid_den_line, label = "GAMER")
        if max(smoothed_grid_den_line)-min(smoothed_grid_den_line) < 0.18:
            plt.ylim(0.9,1.1)
        elif max(smoothed_grid_den_line)-min(smoothed_grid_den_line) < 0.5:
            plt.ylim(0.7,1.3)
        else:
            plt.yscale('log')
            plt.ylim(10**math.floor(math.log10(min(smoothed_grid_den_line))), 1.5*10**math.ceil(math.log10(max(smoothed_grid_den_line))))
        plt.xlabel('x [Mpc/h]')
        plt.ylabel(r'$\rho_\mathrm{gas}/\rho_\mathrm{ave}$')
        plt.title("%s density profile along x-axis at z = %.2f"%(Gas_Par_Setup_T,current_time_z))
        plt.tick_params(top=True, right=True)
        plt.legend(loc="upper right")
        plt.savefig("Data_%06d_x_Density_Line.png"%DumpID_now, dpi = dpi)
        plt.clf()

        # x-component velocity plot
        plt.axhline(y=0, color='0.9', linestyle='--')
        if current_time_z > zc_Collapse:
            plt.plot(x_line_Coordinates - BoxSize_x/2, (100/current_time_a)*MomX_Numerical/Dens_Numerical, label = "GAMER")
            plt.plot(x_line_Coordinates - BoxSize_x/2, (100/current_time_a)*MomX_Analytical/Dens_Analytical, color="r", linestyle='-.', label = "Analytical")
            print("Maximum x-velocity [km/sec]: %f"%np.amax((100/current_time_a)*MomX_Numerical/Dens_Numerical))
        else:
            px_dens_line = plt.plot(smoothed_grid_x_line,smoothed_grid_vel_line, label = "GAMER")
            print("Maximum x-velocity [km/sec]: %f"%np.amax(smoothed_grid_vel_line))
        plt.ylim(-2500,2500)
        plt.xlabel('x [Mpc/h]')
        plt.ylabel('$v_x$ [km/s]')
        plt.title("%s x-component velocity at z = %.2f"%(Gas_Par_Setup_T,current_time_z))
        plt.tick_params(top=True, right=True)
        plt.legend(loc="upper right")
        plt.savefig("Data_%06d_x_Velocity_Line.png"%DumpID_now, dpi = dpi)
        plt.clf()

        # x-line temperature plot
        if current_time_z > zc_Collapse:
            Temp_Analytical = np.array(list(map(Zeldovich_Temp_Analytical, smoothed_grid_den_line)))
            plt.plot(x_line_Coordinates - BoxSize_x/2, Temp_Numerical, label = "GAMER-Output_L1Error()")
            plt.plot(smoothed_grid_x_line,smoothed_grid_temp_line, label = "GAMER-yt")
            plt.plot(smoothed_grid_x_line, Temp_Analytical, color="r", linestyle='-.', label = "Analytical")
        else:
            plt.plot(smoothed_grid_x_line,smoothed_grid_temp_line, label = "GAMER")
        plt.yscale('log')
        plt.ylim(0.9*min(smoothed_grid_temp_line), 1.1*max(smoothed_grid_temp_line))
        plt.xlabel('x [Mpc/h]')
        plt.ylabel('$T$ [K]')
        plt.title("%s Temperature at z = %.2f"%(Gas_Par_Setup_T,current_time_z))
        plt.tick_params(top=True, right=True)
        plt.legend(loc="upper right")
        plt.savefig("Data_%06d_Temperature.png"%DumpID_now, dpi = dpi)
        plt.clf()
        print("Mean Gas Temperature [K]: %f"%np.mean(smoothed_grid_temp_line))

        # merge subplots
        Subplots_Vertical_Combine("Data_%06d_x_Density_Line.png"%DumpID_now,"Data_%06d_x_Velocity_Line.png"%DumpID_now
                                ,"Data_%06d_x_Density_Velocity_Line.png"%DumpID_now)
        Subplots_Vertical_Combine("Data_%06d_x_Density_Velocity_Line.png"%DumpID_now,"Data_%06d_Temperature.png"%DumpID_now
                                ,"Data_%06d_x_Density_Velocity_Temperature_Line.png"%DumpID_now)

        # x-line pressure plot
        if current_time_z > zc_Collapse:
            plt.plot(x_line_Coordinates - BoxSize_x/2, Pres_Numerical, label = "GAMER")
            plt.plot(x_line_Coordinates - BoxSize_x/2, Pres_Analytical, color="r", linestyle='-.', label = "Analytical")
        else:
            px_dens_line = plt.plot(smoothed_grid_x_line,smoothed_grid_pres_line, label = "GAMER")
        plt.xlabel('x [Mpc/h]')
        plt.ylabel('Pressure [UNIT_P]')
        plt.title("%s Gas Pressure at z = %.2f"%(Gas_Par_Setup_T,current_time_z))
        plt.tick_params(top=True, right=True)
        plt.legend(loc="upper right")
        plt.savefig("Data_%06d_Pressure.png"%DumpID_now, dpi = dpi)
        plt.clf()

        # x-line sound speed
        plt.plot(smoothed_grid_x_line,smoothed_grid_sspd_line, label = "GAMER")
        plt.xlabel('x [Mpc/h]')
        plt.ylabel('Sound Speed [km/sec]')
        plt.title("%s Sound Speed [km/sec] at z = %.2f"%(Gas_Par_Setup_T,current_time_z))
        plt.tick_params(top=True, right=True)
        plt.legend(loc="upper right")
        plt.savefig("Data_%06d_Sound_Speed.png"%DumpID_now, dpi = dpi)
        plt.clf()
        print("Mean Sound Speed [km/sec]: %f"%np.mean(smoothed_grid_sspd_line))

        if current_time_z > zc_Collapse:
            # compilation of relative errors
            Temp_Relative_Error = abs(smoothed_grid_temp_line-Temp_Analytical)/Temp_Analytical
            plt.plot(x_line_Coordinates - BoxSize_x/2, Dens_Relative_Error, label = "Density")
            plt.plot(x_line_Coordinates - BoxSize_x/2, MomX_Relative_Error, label = "x-Momentum")
            plt.plot(smoothed_grid_x_line, Temp_Relative_Error, "-.", label = "Temperature")
            plt.plot(x_line_Coordinates - BoxSize_x/2, Pres_Relative_Error, ":", label = "Pressure")
            plt.xlabel('x [Mpc/h]')
            plt.ylabel('Relative Error: |Numerical-Aanalytical|/Aanalytical')
            plt.yscale('log')
            plt.title("%s Relative Errors at z = %.2f"%(Gas_Par_Setup_T,current_time_z))
            plt.tick_params(top=True, right=True)
            plt.legend(loc="upper right")
            plt.savefig("Data_%06d_Relative_Errors.png"%DumpID_now, dpi = dpi)
            plt.clf()

            # merge subplots
            Subplots_Vertical_Combine("Data_%06d_Pressure.png"%DumpID_now,"Data_%06d_Sound_Speed.png"%DumpID_now
                                    ,"Data_%06d_Pressure_Sound_Speed_Line.png"%DumpID_now)
            Subplots_Vertical_Combine("Data_%06d_Pressure_Sound_Speed_Line.png"%DumpID_now,"Data_%06d_Relative_Errors.png"%DumpID_now
                                    ,"Data_%06d_x_Pressure_Sound_Speed_Relative_Errors.png"%DumpID_now)

            #-----------------------------------------------------------------------------------------------------------------
            # output the merged master plot
            pz_dens_1 = matplotlib.image.imread("Data_%06d_Density_Velocity_Slice.png"%DumpID_now)
            pz_dens_2 = matplotlib.image.imread("Data_%06d_x_Density_Velocity_Temperature_Line.png"%DumpID_now)
            pz_dens_3 = matplotlib.image.imread("Data_%06d_x_Pressure_Sound_Speed_Relative_Errors.png"%DumpID_now)

            xdim_1,ydim_1,zdim_1 = pz_dens_2.shape
            xdim_2,ydim_2,zdim_2 = pz_dens_3.shape
            pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_1_zero[0:xdim_1, 0:ydim_1, 0:zdim_1] = pz_dens_2
            pz_dens_2_zero[0:xdim_2, 0:ydim_2, 0:zdim_2] = pz_dens_3
            merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)

            xdim_1,ydim_1,zdim_1 = pz_dens_1.shape
            xdim_2,ydim_2,zdim_2 = merged_image.shape
            pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_1_zero[0:xdim_1, 0:ydim_1, 0:zdim_1] = pz_dens_1
            pz_dens_2_zero[0:xdim_2, 0:ydim_2, 0:zdim_2] = merged_image
            merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)
            matplotlib.image.imsave("Combine_Data_%06d.png"%DumpID_now, merged_image)

            os.remove("./Data_%06d_Density_Velocity_Slice.png"%DumpID_now)
            os.remove("./Data_%06d_x_Pressure_Sound_Speed_Relative_Errors.png"%DumpID_now)
            os.remove("./Data_%06d_x_Density_Velocity_Temperature_Line.png"%DumpID_now)
        else:
            # merge subplots
            Subplots_Vertical_Combine("Data_%06d_Pressure.png"%DumpID_now,"Data_%06d_Sound_Speed.png"%DumpID_now
                                    ,"Data_%06d_Pressure_Sound_Speed_Line.png"%DumpID_now)

            #-----------------------------------------------------------------------------------------------------------------
            # output the merged master plot
            pz_dens_1 = matplotlib.image.imread("Data_%06d_Density_Velocity_Slice.png"%DumpID_now)
            pz_dens_2 = matplotlib.image.imread("Data_%06d_x_Density_Velocity_Temperature_Line.png"%DumpID_now)
            pz_dens_3 = matplotlib.image.imread("Data_%06d_Pressure_Sound_Speed_Line.png"%DumpID_now)

            xdim_1,ydim_1,zdim_1 = pz_dens_2.shape
            xdim_2,ydim_2,zdim_2 = pz_dens_3.shape
            pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_1_zero[0:xdim_1, 0:ydim_1, 0:zdim_1] = pz_dens_2
            pz_dens_2_zero[0:xdim_2, 0:ydim_2, 0:zdim_2] = pz_dens_3
            merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)

            xdim_1,ydim_1,zdim_1 = pz_dens_1.shape
            xdim_2,ydim_2,zdim_2 = merged_image.shape
            pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
            pz_dens_1_zero[0:xdim_1, 0:ydim_1, 0:zdim_1] = pz_dens_1
            pz_dens_2_zero[0:xdim_2, 0:ydim_2, 0:zdim_2] = merged_image
            merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)
            matplotlib.image.imsave("Combine_Data_%06d.png"%DumpID_now, merged_image)

            os.remove("./Data_%06d_Density_Velocity_Slice.png"%DumpID_now)
            os.remove("./Data_%06d_Pressure_Sound_Speed_Line.png"%DumpID_now)
            os.remove("./Data_%06d_x_Density_Velocity_Temperature_Line.png"%DumpID_now)
    # -------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------
    elif (Gas_Par_Setup == 2): # particle-only setup
        # import the particle attribute table
        data                 = ds.all_data()
        ParMass              = data[('all', 'ParMass')].in_units("code_mass").value
        ParPosX              = data[('all', 'ParPosX')].in_units("code_length").value                          # [Mpc/h]
        ParVelX              = data[('all', 'ParVelX')].in_units("code_velocity").value*(100/current_time_a)   # [km/s]
        ParticleNumber       = len(ParMass)
        ParDenBin            = [[] for _ in xrange(Analysis_Par_Bin)]
        ParDenListx          = []
        ParDenListValue      = []
        for i_now in range(ParticleNumber):
            BinIndexNow = int(math.floor(ParPosX[i_now]/x_Len_per_Bin))
            ParDenBin[BinIndexNow].append(i_now)
        for i_now, list_now in enumerate(ParDenBin):
            ParDenListx.append((i_now+0.5)*x_Len_per_Bin)
            ParDenListValue.append(ParMass[0]*len(list_now)/(x_Len_per_Bin*BoxSize_y*BoxSize_z))

        # Zeldovich pancake collapse analytical solutions at current_time_z
        def Zeldovich_coord_transformation_solver(x_code_length):
            x_Lagrangian = []
            for x_code_now in x_code_length:
                f = lambda x_Lagrangian_Trial: (x_Lagrangian_Trial - ((1+zc_Collapse)/(1+current_time_z))*math.sin(2*math.pi*(x_Lagrangian_Trial/Pert_Wave_Len))/(2*math.pi/Pert_Wave_Len)) - x_code_now
                x_Lagrangian.append(fsolve(f, x_code_now))
            return x_Lagrangian
        def Zeldovich_Dens_Analytical(x_Lagrangian):
            return 1/(1-((1+zc_Collapse)/(1+current_time_z))*math.cos(2*math.pi*(x_Lagrangian/Pert_Wave_Len)))
        def Zeldovich_Vel_Analytical(x_Lagrangian):
            return -Hubble0_Units*((1+zc_Collapse)/(1+current_time_z)**0.5)*math.sin(2*math.pi*(x_Lagrangian/Pert_Wave_Len))/(2*math.pi/Pert_Wave_Len)*100

        # -------------------------------------------------------------------------------------------------------------------------
        # [yt] particle projection plot
        projplot = yt.ParticleProjectionPlot(ds, "z", ("all", "particle_mass"))#, center="c", width=(width_value, "Mpc/h"))
        projplot.annotate_title("%s Setup: Data_%06d, Proj_Axis = %s"%(Gas_Par_Setup_T, ds.parameters["DumpID"],"z"))
        projplot.annotate_timestamp(corner='upper_right', redshift=True, time=False, text_args={'color':'k'})
        projplot.set_colorbar_label(("all", "particle_mass"),"Projected Mass [M$_\odot$]")
        projplot.set_axes_unit ("code_length")
        projplot.set_unit(("all", "particle_mass"), "Msun")
        projplot.save()

        # output interpolated halo density and dispersion profiles
        uniform_grid_x_line = np.linspace(0,BoxSize_x,1000)
        radius_sample = np.logspace(1.3,5,100)
        fig, axs      = plt.subplots(2, figsize=(10.5,12.5))
        # x-axis density profile
        axs[0].axhline(y=1, color='0.9', linestyle='--')
        if current_time_z > zc_Collapse:
            axs[0].plot(uniform_grid_x_line, np.array(list(map(Zeldovich_Dens_Analytical, Zeldovich_coord_transformation_solver(uniform_grid_x_line-BoxSize_x/2)))), color="r", linestyle='-.', label = "Analytical")
        axs[0].plot(ParDenListx, ParDenListValue, alpha=0.175)
        axs[0].scatter(ParDenListx, ParDenListValue, s=12, facecolors='none', edgecolors='C0', label = "GAMER")
        axs[0].set_title("NPar_AllRank: %i, Scaled Density at z = %.2f"%(ParticleNumber,current_time_z), fontsize=18)
        axs[0].set_xlabel('x [Mpc/h]', fontsize=14)
        axs[0].set_ylabel(r'$\rho(x)$/$\rho_\mathrm{crit}$', fontsize=14)
        axs[0].tick_params(top=True, right=True, labelsize=14)
        axs[0].legend(loc="upper right", fontsize=14)
        # x-component velocity plot
        axs[1].axhline(y=0, color='0.9', linestyle='--')
        if current_time_z > zc_Collapse:
            axs[1].plot(uniform_grid_x_line, np.array(list(map(Zeldovich_Vel_Analytical, Zeldovich_coord_transformation_solver(uniform_grid_x_line-BoxSize_x/2)))), color="r", linestyle='-.', label = "Analytical")
        axs[1].scatter(ParPosX, ParVelX, s=12, facecolors='none', edgecolors='C0', label = "GAMER")
        axs[1].set_ylim(-2500,2500)
        axs[1].set_title("NPar_AllRank: %i, $v_x$ at z = %.2f"%(ParticleNumber,current_time_z), fontsize=18)
        axs[1].set_xlabel('x [Mpc/h]', fontsize=14)
        axs[1].set_ylabel('$v_x$ [km/s]', fontsize=14)
        axs[1].tick_params(top=True, right=True, labelsize=14)
        axs[1].legend(loc="upper right", fontsize=14)
        fig.tight_layout()
        plt.savefig("./Data_%06d_x_Density_Velocity_Line.png"%DumpID_now)
        plt.close()
        print("Maximum x-velocity [km/sec]: %f"%np.amax(ParVelX))

        # -------------------------------------------------------------------------------------------------------------------------
        # merge subplots
        pz_dens_1 = matplotlib.image.imread("Data_%06d_Particle_z_particle_mass.png"%DumpID_now)
        pz_dens_2 = matplotlib.image.imread("Data_%06d_x_Density_Velocity_Line.png"%DumpID_now)
        xdim_1,ydim_1,zdim_1 = pz_dens_1.shape
        xdim_2,ydim_2,zdim_2 = pz_dens_2.shape
        pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
        pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
        pz_dens_1_zero[0:xdim_1, 0:ydim_1, 0:zdim_1] = pz_dens_1
        pz_dens_2_zero[0:xdim_2, 0:ydim_2, 0:zdim_2] = pz_dens_2
        merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)
        matplotlib.image.imsave("Combine_Data_%06d.png"%DumpID_now, merged_image)
        os.remove("./Data_%06d_Particle_z_particle_mass.png"%DumpID_now)
        os.remove("./Data_%06d_x_Density_Velocity_Line.png"%DumpID_now)
