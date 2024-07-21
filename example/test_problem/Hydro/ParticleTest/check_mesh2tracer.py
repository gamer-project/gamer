#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")


import re

import h5py
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 16


### helper function
def comp_momx(ParX, ParY, Dens_Bg, Ang_Freq, BoxSize):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the linear momentum in the x direction at these particles' locations
    BoxSize_Half = 0.5 * BoxSize
    Radius       = np.hypot(ParX - BoxSize_Half, ParY - BoxSize_Half)
    Sin_Theta    = (ParY - BoxSize_Half) / Radius
    Velocity     = Ang_Freq * Radius

    return -Dens_Bg * Velocity * Sin_Theta


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
fn_in = "Data_000000"
hdf_data = h5py.File(fn_in, "r")

field_par  = ["ParPosX", "ParPosY", "ParType"]
field_mesh = ["Mesh{}".format(label)
              for label in np.genfromtxt("Input__Par_Mesh", dtype = "|U")]

for key in field_par + field_mesh:
    cmd = "{} = hdf_data['Particle']['{}'][:]".format(key, key)
    exec(cmd)

hdf_data.close()


### postprocess
cond_tracer = (ParType == 0.0)

print("Number of tracer particles: {}/{}".format(cond_tracer.sum(), cond_tracer.size))

## check non-tracer particles
for field in field_mesh:
    cmd = "set_nontracer = set({}[~cond_tracer])".format(field)
    exec(cmd)

    if len(set_nontracer) > 1:
        print("mesh quantities on non-tracer particles have more than two values: ",
              " ".join(str(i) for i in set_nontracer))

    # the value for non-tracer particles should be __FLT_MAX__
    print("{} on non-tracer particles: {}".format(field, set_nontracer))

## check tracer particles
for key in field_par + field_mesh:
    cmd = "{}_tracer = {}[cond_tracer]".format(key, key)
    exec(cmd)

for key in field_mesh:
    obj_tracer = "{}_tracer".format(key)
    cmd = "minval = {}.min(); maxval = {}.max()".format(obj_tracer, obj_tracer)
    exec(cmd)

    print("Min/Max of {} on tracer particles: {} / {}".format(key, minval, maxval))

# compute the linear momentum in the x direction
MomX_Par = MeshDens_tracer * MeshVelX_tracer
MomX_IC  = comp_momx(ParPosX_tracer, ParPosY_tracer,
                     param["ParTest_Dens_Bg"], param["ParTest_Ang_Freq"], param["BOX_SIZE"])
Radius   = np.hypot(ParPosX_tracer - 0.5 * param["BOX_SIZE"],
                    ParPosY_tracer - 0.5 * param["BOX_SIZE"])

# create a dataset for sorting the data based on the x-coordinate position
dataset = zip(Radius, MeshDens_tracer, MeshPres_tracer, MomX_Par, MomX_IC)
dataset = sorted(dataset)
dataset = [np.array(data) for data in zip(*dataset)]

reldiff_dens = np.abs(dataset[1] - param  ["ParTest_Dens_Bg"]) / param  ["ParTest_Dens_Bg"]
reldiff_pres = np.abs(dataset[2] - param  ["ParTest_Pres_Bg"]) / param  ["ParTest_Pres_Bg"]
reldiff_momx = np.abs(dataset[3] - dataset[4]                ) / dataset[4]

# visualization
fig, axes = plt.subplots(figsize = (8, 10), nrows = 3, sharex = True)

axes[0].scatter(dataset[0], reldiff_dens)
axes[1].scatter(dataset[0], reldiff_pres)
axes[2].scatter(dataset[0], reldiff_momx)

axes[2].set_xlabel("Radius of Tracer Particles")

axes[0].set_ylabel("Relative Difference\nin Dens")
axes[1].set_ylabel("Relative Difference\nin Pres")
axes[2].set_ylabel("Relative Difference\nin MomX")

fig.tight_layout()
plt.savefig("check_mesh2tracer.png", bbox_inches = "tight")
plt.close()
