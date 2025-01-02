#!/usr/bin/env python3
#
# Check the mesh attributes stored in the HDF5 snapshots
# at the locations of the tracer particles
#

import matplotlib
matplotlib.use("Agg")


import re
from glob import glob

import h5py
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 16


### helper function
def comp_dens(ParX, ParY, Center, Dens_Bg, BoxSize):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the density at these particles' locations
    Radius = np.hypot(ParX - Center[0], ParY - Center[1])

    return Dens_Bg * (1.0 + 5.0 * Radius / BoxSize)


def comp_pres(ParX, ParY, Center, Pres_Bg, BoxSize):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the pressure at these particles' locations
    Radius = np.hypot(ParX - Center[0], ParY - Center[1])

    return Pres_Bg * (1.0 + 5.0 * Radius / BoxSize)


def comp_velx(ParX, ParY, Center, Ang_Freq):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the velocity in the x direction at these particles' locations
    Radius    = np.hypot(ParX - Center[0], ParY - Center[1])
    Sin_Theta = (ParY - Center[1]) / Radius
    Velocity  = Ang_Freq * Radius

    return -1.0 * Velocity * Sin_Theta


### retrieve runtime parameters in Input__Parameter and Input__TestProb
regex_num = r"\s*([-+]?\d+\.?\d*[eE]?[-+]?\d*)"

param = dict()

# Input__Parameter
with open("Input__Parameter", "r") as f:
    param_in = f.read()

    key_list = "BOX_SIZE",
    for key in key_list:
        value = re.findall(key + regex_num, param_in)

        # assume the value is a float
        param[key] = float(value[0])

# Input__TestProb
with open("Input__TestProb", "r") as f:
    param_in = f.read()

    key_list = "ParTest_Dens_Bg", "ParTest_Pres_Bg", "ParTest_Ang_Freq"
    for key in key_list:
        value = re.findall(key + regex_num, param_in)

        # assume the value is a float
        param[key] = float(value[0])


### load data
fn_list = glob("Data_" + "[0-9]" * 6)
fn_in = sorted(fn_list)[-1]  # use the final snapshot
print("Loading data in {}".format(fn_in))

hdf_data = h5py.File(fn_in, "r")

field_par  = ["ParPosX", "ParPosY", "ParType"]
field_mesh = ["Mesh{}".format(label)
              for label in np.genfromtxt("Input__Par_Mesh", dtype = "|U")]

for key in field_par + field_mesh:
    cmd = "{} = hdf_data['Particle']['{}'][:]".format(key, key)
    exec(cmd)

time_phys = hdf_data["Info"]["KeyInfo"]["Time"][0]

hdf_data.close()


### post-process
cond_tracer = (ParType == 0)

print("Number of tracer particles: {}/{}".format(cond_tracer.sum(), cond_tracer.size))

## check non-tracer particles
for field in field_mesh:
    cmd = "set_nontracer = set({}[~cond_tracer])".format(field)
    exec(cmd)

    if len(set_nontracer) > 1:
        print("mesh quantities on non-tracer particles have more than one value: ",
              " ".join(str(i) for i in set_nontracer))
    else:
        # the value for non-tracer particles should be __FLT_MAX__
        print("{} on non-tracer particles: {}".format(field, set_nontracer))

## check tracer particles
for key in field_par + field_mesh:
    cmd = "{}_tracer = {}[cond_tracer]".format(key, key)
    exec(cmd)

# compute the linear momentum in the x direction
Center_Bg  = 0.25 * param["BOX_SIZE"], 0.25 * param["BOX_SIZE"]
Center_Mom = 0.50 * param["BOX_SIZE"], 0.50 * param["BOX_SIZE"]

Radius   = np.hypot(ParPosX_tracer - Center_Bg[0], ParPosY_tracer - Center_Bg[1])
Dens_ref = comp_dens(ParPosX_tracer, ParPosY_tracer, Center_Bg,  param["ParTest_Dens_Bg"], param["BOX_SIZE"])
Pres_ref = comp_pres(ParPosX_tracer, ParPosY_tracer, Center_Bg,  param["ParTest_Pres_Bg"], param["BOX_SIZE"])
VelX_ref = comp_velx(ParPosX_tracer, ParPosY_tracer, Center_Mom, param["ParTest_Ang_Freq"])

# create a data set and sort the data based on the x-coordinate position
dataset = zip(Radius, MeshDens_tracer, Dens_ref, MeshPres_tracer, Pres_ref, MeshVelX_tracer, VelX_ref)
dataset = sorted(dataset)
dataset = [np.array(data) for data in zip(*dataset)]

# compute the relative difference
reldiff_list = [np.abs(data_interp / data_ref - 1.0)
                for data_interp, data_ref in zip(dataset[1::2], dataset[2::2])]

# print information
for reldiff, key in zip(reldiff_list, field_mesh):
    obj_tracer = "{}_tracer".format(key)
    cmd = "minval = {}.min(); maxval = {}.max()".format(obj_tracer, obj_tracer)
    exec(cmd)

    print("Min/Max         of {} on tracer particles: {:14.7e} / {:14.7e}".format(key, minval, maxval))
    print("Min/Max RelDiff of {} on tracer particles: {:14.7e} / {:14.7e}".format(key, reldiff.min(), reldiff.max()))

## visualization
fig = plt.figure(figsize = (14, 8))
ax1 = plt.subplot2grid((5, 3), (0, 0), rowspan = 3)
ax2 = plt.subplot2grid((5, 3), (3, 0), rowspan = 2)
ax3 = plt.subplot2grid((5, 3), (0, 1), rowspan = 3)
ax4 = plt.subplot2grid((5, 3), (3, 1), rowspan = 2)
ax5 = plt.subplot2grid((5, 3), (0, 2), rowspan = 3)
ax6 = plt.subplot2grid((5, 3), (3, 2), rowspan = 2)

axes_list_data    = ax1, ax3, ax5
axes_list_reldiff = ax2, ax4, ax6

# field quantity
label_list  = "Interpolation", "Reference"
marker_list = "o", "x"
ylabel_list = "Dens", "Pres", "VelX"

for idx_1, (ax, ylabel) in enumerate(zip(axes_list_data, ylabel_list)):
    for idx_2, (label, marker) in enumerate(zip(label_list, marker_list), 2 * idx_1 + 1):
        ax.scatter(dataset[0], dataset[idx_2], label = label, marker = marker)

    ax.set_ylabel(ylabel)
    ax.legend(framealpha = 0)

# relative difference
for reldiff, ax in zip(reldiff_list, axes_list_reldiff):
    ax.scatter(dataset[0], reldiff)

    ax.set_yscale("log")
    ax.set_xlabel("Radius of Tracer Particles")
    ax.set_ylabel("Relative Difference")

fig.suptitle("Time = {}".format(time_phys))
fig.tight_layout()
plt.savefig("check_mesh2tracer.png", bbox_inches = "tight")
plt.close()
