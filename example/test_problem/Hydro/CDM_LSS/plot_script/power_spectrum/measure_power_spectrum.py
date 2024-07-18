# measure high resolution power spectrum using nbodykit
# to install nbodykit, follow the instruction from https://nbodykit.readthedocs.io/en/latest/getting-started/install.html
# hint: install cython version 0.29 instead of the latest version to avoid gil error (by pip install Cython==0.29.37)
from nbodykit.lab import RandomCatalog, FFTPower
import numpy as np
import h5py

target = "Data_000072" # snapshot at z=0
Nmesh = 512  # mesh resolution

snapshot = h5py.File("../../"+target, "r")

npar = snapshot['Info']['KeyInfo']['Par_NPar']
boxsize = snapshot['Info']['KeyInfo']['BoxSize']

mass = np.array(snapshot["Particle/ParMass"])
posx = np.array(snapshot["Particle/ParPosX"])
posy = np.array(snapshot["Particle/ParPosY"])
posz = np.array(snapshot["Particle/ParPosZ"])

pos = np.vstack((posx, posy, posz)).T

# create a catalog object
particle_catalog = RandomCatalog(csize=npar) # overwritten by the actual values below

particle_catalog['Mass'] = mass
particle_catalog['Position'] = pos

# create a mesh object
mesh = particle_catalog.to_mesh(Nmesh=Nmesh, BoxSize=boxsize)#, value='Mass')

# compute the power spectrum
boxlen = boxsize[0]
r = FFTPower(mesh, mode='1d', kmin=1/boxlen, kmax=Nmesh/boxlen, dk=0.1)

# the result is stored at "power" attribute
Pk = r.power

# print out the meta-data
for k in r.attrs:
    print("%s = %s" %(k, str(r.attrs[k])))

# save P(k)
out = [Pk['k'], Pk['power'].real]
npout = np.array(out).T
np.savetxt("powerspec_highres_"+target+".dat", npout)


