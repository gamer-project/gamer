import argparse
import sys
import yt
import matplotlib.pyplot as plt
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='./' )
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
#   print str(sys.argv[t]), # for python2.X
   print(str(sys.argv[t]),end=', ') # for python3.X
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start    = args.idx_start
idx_end      = args.idx_end
didx         = args.didx
prefix       = args.prefix

dpi          = 150
field        = ('gas','particle_density_on_grid')
r_sphere     = 1.5

yt.enable_parallelism()

d0           = yt.load(  prefix+'Data_000000' )
center_pos_0 = d0.find_max(field)[1].in_units("code_length")  # return maximum value and position, here use [1] to take peak position only
sp0          = d0.sphere(center_pos_0, (r_sphere, "code_length"))
rp0          = yt.create_profile(sp0,"radius",field,units={("radius"): "code_length"},logs={("radius"): True})

ts = yt.DatasetSeries([ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ])
for ds in ts.piter():
#   print(str(ds)+".png")
   center_pos_1 = ds.find_max(field)[1].in_units("code_length")  # return maximum value and position, here use [1] to take peak position only
   sp1          = ds.sphere(center_pos_1, (r_sphere, "code_length"))
   rp1          = yt.create_profile(sp1,"radius",field,units={("radius"): "code_length"},logs={("radius"): True})

   fig          = plt.figure()
   ax           = fig.add_subplot(111)
   ax.plot(rp0.x.value,rp0[field].in_units("code_mass/code_length**3").value)
   ax.plot(rp1.x.value,rp1[field].in_units("code_mass/code_length**3").value)
   ax.set_xscale('log')
   ax.set_yscale('log')
   ax.legend(["Later", "Original"])
   fig.savefig(str(ds)+"_profile.png")
