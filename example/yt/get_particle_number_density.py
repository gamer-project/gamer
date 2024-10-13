import sys, yt
import numpy as np

#### define derived field for particle number density ####
def number_density(field, data):
   number_count = data["deposit","all_count"]
   number_density = number_count / data["index","cell_volume"]
   return number_density
####

#### Input check ####
if len(sys.argv) !=5:
   print("Wrong input number! Usage: ./EXE DATA_PATH START_ID END_ID SEPARATION")
   sys.exit(1)
####

yt.enable_parallelism()
radius_range     = 45.                         # radius range for sphere covering the particles, in unit of kpc
N_BIN            = 16
field            = ("gamer","number_density")
data_saving_path = '.'
r_min, r_max     = 0.4, 40.0                   # minimum and maximum for computing the particle number density profile

PREFIX = sys.argv[1]                           # path for snapshot
START_ID = np.int(sys.argv[2])                 # start ID for snapshot
END_ID = np.int(sys.argv[3])                   # end ID for snapshot
SEPARATION = np.int(sys.argv[4])               # separatioon between two snapshots

#ts = yt.load( [ '%s/Data_%06d'%(PREFIX,idx) for idx in np.arange(START_ID, END_ID+1, SEPARATION) ] )
ts = yt.DatasetSeries( [ '%s/Data_%06d'%(PREFIX,idx) for idx in np.arange(START_ID, END_ID+1, SEPARATION) ] )
for ds in ts.piter():
    print("Working on %s..."%str(ds))
#### add number desnity field ####
    ds.add_field(
    field,
    function=number_density,
    sampling_type="cell",
    units="1/cm**3",)
####

#### compute profile ####
    source = ds.sphere( "m", (radius_range, "kpc"))
    prof = yt.create_profile(source, "radius", field, n_bins=N_BIN, units = {'radius': 'kpc',field: '1/cm**3'}, extrema = {'radius': ((r_min, 'kpc'), (r_max, 'kpc'))}, weight_field='cell_volume')
    r = prof.x
    density = prof["number_density"]
####
#### save data ####
    print("Start to write file %s/profile_%s.txt ..."%(data_saving_path,str(ds)))
    with open("%s/profile_%s.txt"%(data_saving_path,str(ds)),"w") as f:
        np.savetxt(f, np.vstack((r,density)).T, fmt='%.8e', header="Radius(kpc)\tDensity(1/cm^3)")
####
