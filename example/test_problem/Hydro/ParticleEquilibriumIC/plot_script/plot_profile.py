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
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

center_mode = 'max'
dpi         = 150


yt.enable_parallelism()

d0= yt.load(  prefix+'Data_000000' )
#ts = yt.load( [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
files = [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ]
for f in files:
   print(f+".png")
   ds  = yt.load(f)
   sp1 = ds.sphere("c", (1.5, "code_length")) #ds.sphere( center_mode, 0.5*ds.domain_width.to_value().max() )
   sp0 = d0.sphere("c", (1.5, "code_length"))
   field ='particle_density_on_grid'
   rp0 = yt.create_profile(sp0,"radius",field,units={("radius"): "code_length"},logs={("radius"): True})
   rp1 = yt.create_profile(sp1,"radius",field,units={("radius"): "code_length"},logs={("radius"): True})
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.plot(rp0.x.value,rp0[field].in_units("code_mass/code_length**3").value)
   ax.plot(rp1.x.value,rp1[field].in_units("code_mass/code_length**3").value)
   ax.set_xscale('log')
   ax.set_yscale('log')
   ax.legend(["Later", "Original"])
   fig.savefig(f+"_profile.png")
