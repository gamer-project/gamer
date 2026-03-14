# measure high resolution power spectrum using nbodykit
# to install nbodykit, follow the instruction from https://nbodykit.readthedocs.io/en/latest/getting-started/install.html
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
# ref: https://nbodykit.readthedocs.io/en/latest/catalogs/mock-data.html#RandomCatalog
particle_catalog = RandomCatalog(csize=npar) # overwritten by the actual values below

particle_catalog['Mass'] = mass
particle_catalog['Position'] = pos

# create a mesh object 
# ref: https://nbodykit.readthedocs.io/en/latest/api/_autosummary/nbodykit.base.catalog.html#nbodykit.base.catalog.CatalogSource.to_mesh
mesh = particle_catalog.to_mesh(Nmesh=Nmesh, BoxSize=boxsize)

# compute the power spectrum
# ref: https://nbodykit.readthedocs.io/en/latest/api/_autosummary/nbodykit.algorithms.fftpower.html#nbodykit.algorithms.fftpower.FFTPower
#      https://nbodykit.readthedocs.io/en/latest/results/algorithms/fftpower.html#fftpower
boxlen = boxsize[0]
r = FFTPower(mesh, mode='1d', kmin=2*np.pi/boxlen, kmax=np.pi*Nmesh/boxlen, dk=2*np.pi/boxlen)

# the result is stored at "power" attribute
Pk = r.power

# print out the meta-data
for k in r.attrs:
    print("%s = %s" %(k, str(r.attrs[k])))

# save P(k)
out = [Pk['k'], Pk['power'].real]
npout = np.array(out).T
np.savetxt("powerspec_highres_"+target+".dat", npout)


