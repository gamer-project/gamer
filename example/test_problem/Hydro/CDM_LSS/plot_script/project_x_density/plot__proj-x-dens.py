import argparse
import sys
import numpy as np
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Projection of mass density' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
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
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

field           = 'density'
colormap_dens   = 'algae'
center_mode     = 'c'
dpi             = 150
projection_axis = "x"
width_value     = 30

yt.enable_parallelism()
ts                   = yt.DatasetSeries( [ prefix+'../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
ts0                  = ts[0]
base_level_cell_num  = int(ts0.domain_dimensions[1])
max_AMR_level        = int(ts0.parameters["MaxLevel"])
NPar_base            = int(round(np.power(int(ts0["Par_NPar"]), 0.3333333)))

for ds in ts.piter():
   p = yt.ParticleProjectionPlot(ds, projection_axis, ("all", "particle_mass"), center="c", width=(width_value, "Mpc/h"))
   p.annotate_title("GAMER-%i$^3$: Data_%06d, Proj_Axis = %s"%(NPar_base,ds.parameters["DumpID"],projection_axis))
   p.annotate_timestamp(corner='upper_right', redshift=True, time=False, text_args={'color':'k'})
   p.set_colorbar_label(("all", "particle_mass"),"Projected Particle Mass [M$_\odot$]")
   p.set_zlim(("all", "particle_mass"), 1e8, 6e12)
   p.set_unit(("all", "particle_mass"), "Msun")
   # see https://yt-project.org/doc/visualizing/callbacks.html
   p.annotate_text((0,0.8,0.8),"[NX0_TOT:%i$^3$][AMR_MAX:%i]"%(base_level_cell_num,max_AMR_level),coord_system="data",text_args={"color": "black"},
   inset_box_args={
        "boxstyle": "square,pad=0.2",
        "facecolor": "white",
        "linewidth": 3,
        "edgecolor": "white",
        "alpha": 0.5,})
   p.save()
