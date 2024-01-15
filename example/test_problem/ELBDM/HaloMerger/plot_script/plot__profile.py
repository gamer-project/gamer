import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Density profile' )

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

center_mode = 'max'
r_sphere    = 0.06
dpi         = 150
nbin        = 32

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   for field in ['density']:

      center_pos  = ds.find_max(field)[1].in_units("code_length")  # return maximum value and position, here use [1] to take peak position only
      sp = ds.sphere(center_pos, (r_sphere, "code_length") )

      prof = yt.ProfilePlot( sp, 'radius', field, weight_field='cell_volume', n_bins=nbin )
      prof.set_unit( 'radius', 'kpc' )
      prof.set_xlim( 1.0e-1, 1.0e2 )
      prof.set_ylim( field, 1.0e0, 1.0e6 )
      prof.set_unit(  field, 'code_density' )
      prof.annotate_title("t = %13.7e Gyr"%(ds.current_time.in_units("Gyr")) )

      prof.save( mpl_kwargs={"dpi":dpi} )

      prof_dens = yt.create_profile( sp, 'radius', fields=field,
                                     weight_field='cell_volume', n_bins=nbin ,
                                     units={'radius': 'kpc',field: 'code_density'})

      np.savetxt('%s_%s_profile'%(ds,field), np.column_stack( (prof_dens.x.in_units("kpc").d, prof_dens[field].in_units("code_density").d)              ), fmt='          %9.8e',  header='                r (kpc)  density (code_density)'              )
