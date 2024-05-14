#!/bin/python3
import yt
import sys
import numpy as np
import matplotlib.pylab as plt


sys.dont_write_bytecode = True
import ana_solution as ana



#====================================================================================================
# Test probelm parameters
#====================================================================================================
# 1. gaussian ball
K_PARA   = 0.05                         # parallel diffusion coefficient
K_PERP   = 0.5                          # perpendicular diffusion coefficient
MAG      = np.array([5., 5., 5.])       # magnetic field
CR_BG    = 1.e-1                        # cosmic ray energy density background value
CR_E0    = 1.e+2                        # cosmic ray energy density amplitude
CR_R2    = 40                           # inverse of variance of Gaussian distribution
ANA_FUNC = ana.ana_ball                 # analytical solution function
all_args = {"K_PARA":K_PARA, "K_PERP":K_PERP, "MAG":MAG, "CR_BG":CR_BG, "CR_E0":CR_E0, "CR_R2":CR_R2}

# 2. step function ring
# K_PARA   = 0.05                         # parallel diffusion coefficient
# CR_BG    = 1.e-1                        # cosmic ray energy density background value
# CR_E0    = 1.e-1                        # cosmic ray energy density amplitude
# R_IN     = 0.5                          # inner radius of the ring
# R_OUT    = 0.7                          # outer radius of the ring
# PLANE    = "xy"                         # the plane where the ring exist [xy/xz/yz]
# ANA_FUNC = ana.ana_step_ring            # analytical solution function
# all_args = {"K_PARA":K_PARA, "CR_BG":CR_BG, "CR_E0":CR_E0, "R_IN":R_IN, "R_OUT":R_OUT, "PLANE":PLANE}

# 3. gaussian distribution ring
# K_PARA   = 0.05                         # parallel diffusion coefficient
# CR_BG    = 0.1                          # cosmic ray energy density background value
# CR_E0    = 0.1                          # cosmic ray energy density amplitude
# R        = 0.6                          # center radius of the ring
# DEL_R    = 0.05                         # the standard deviation on r direction
# DEL_PHI  = 0.5                          # the standard deviation on phi direction
# PLANE    = "xy"                         # the plane where the ring exist [xy/xz/yz]
# ANA_FUNC = ana.ana_gaussian_ring        # analytical solution function
# all_args = {"K_PARA":K_PARA, "CR_BG":CR_BG, "CR_E0":CR_E0, "R":R, "DEL_R":DEL_R, "DEL_PHI":DEL_PHI, "PLANE":PLANE}

# 4. gaussian distribution plane
# K_PARA   = 0.05                         # parallel diffusion coefficient
# CR_BG    = 0.1                          # cosmic ray energy density background value
# CR_E0    = 0.1                          # cosmic ray energy density amplitude
# CR_R2    = 40                           # inverse of variance of Gaussian distribution
# DIR      = "x"                          # the direction of gaussian and magnetic field. [x/y/z/xy/xz/yz/xyz]
# ANA_FUNC = ana.ana_plane                # analytical solution function
# all_args = {"K_PARA":K_PARA, "CR_BG":CR_BG, "CR_E0":CR_E0, "CR_R2":CR_R2, "DIR":DIR}



#====================================================================================================
# Main
#====================================================================================================
N_res      = 4          # number of resolutions
base_res   = 6          # base resolution N=2**base_res
center     = 0.5        # center coordinate of the box
data_index = 10         # the data dump index to be analyzed

all_res    = np.array([ int(2**i) for i in range(base_res, base_res+N_res) ], dtype=np.int32)
all_err    = np.ones(N_res)

for i in range(N_res):
    ds = yt.load("./res_%04d/Data_%06d"%(all_res[i], data_index))
    sl = ds.r[:, :, :]

    cray   = np.array(sl["CRay"])
    plot_x = np.array(sl["x"]) - center
    plot_y = np.array(sl["y"]) - center
    plot_z = np.array(sl["z"]) - center
    time   = ds.current_time / ds.time_unit

    ana    = ANA_FUNC( plot_x, plot_y, plot_z, time, **all_args )

    all_err[i] = np.mean( np.abs(1. - cray / ana) )

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
coeff = np.polyfit( np.log(all_res), np.log(all_err), 1 )
y_N1 = np.exp(coeff[1]) * all_res**(-1.)
y_N2 = np.exp(coeff[1]) * all_res**(-2.)
y_N1 = y_N1 / y_N1[-1] * all_err[-1]
y_N2 = y_N2 / y_N2[-1] * all_err[-1]

ax.plot(all_res, ( np.exp(coeff[1]) * all_res**coeff[0] ), label = "$N^{%.2f}$"%coeff[0])
ax.scatter(all_res, all_err, label = "L1error")
#ax.plot(all_res, y_N2, '-', label = "$N^{%.2f}$"%(-2.))
ax.set(xlabel = "N", ylabel = "L1Err")
ax.set(xscale = "log", yscale = "log")
ax.legend(loc = 1)

#plt.suptitle("anisotropic diffusion, t_GAMER = %.2f"%time)
save_name = "L1error.png"
if save_name != '':
    plt.savefig(save_name, dpi = 200)
plt.close()
#plt.show()
