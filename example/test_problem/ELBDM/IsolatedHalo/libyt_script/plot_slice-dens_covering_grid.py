#!/usr/bin/env python3.7

#reference 1 uniform grid: https://yt-project.org/doc/examining/generic_array_data.html#Generic-Unigrid-Data
#reference 2 parallel    : https://yt-project.org/doc/analyzing/parallel_computation.html#parallelizing-over-multiple-objects

import argparse
import sys
import yt
import glob
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Slice of mass density' )

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

num_procs     = 4
lv            = 1
colormap_dens = 'jet'
center_mode   = 'c'
dpi           = 512
buff_size     = 1024
h_0           = 0.6955
box_size      = 0.0875 # Mpc/h
proj_axes     = ['x']
suffix        = '_central_halo'
field         = 'Dens'

yt.enable_parallelism()
idx_range = np.arange(idx_start, idx_end+1, didx)
file_list = ( [ '%s/covering-grid_test_Data_%06d_lv=%d.npz'%(prefix,idx,lv) for idx in idx_range ] )
my_storage = {}

for sto, fn in yt.parallel_objects(file_list, num_procs, storage=my_storage):
   dens = np.load(fn)[field]
   idx  = int(fn.split("_")[-2])   # get the data index
   data = dict(Dens = (dens,"code_mass/code_length**3"))  # in unit of background density
   bbox = np.array([ [-box_size/2./h_0*1e3,box_size/2./h_0*1e3],\
                     [-box_size/2./h_0*1e3,box_size/2./h_0*1e3],\
                     [-box_size/2./h_0*1e3,box_size/2./h_0*1e3]])  # Mpc/h -> kpc
   ds   = yt.load_uniform_grid(data, dens.shape, length_unit = "kpc", bbox=bbox, nprocs = num_procs)

   for proj_axis in proj_axes:
       sz_dens = yt.SlicePlot( ds, proj_axis, field, width = (50, 'kpc'), center=center_mode)
       sz_dens.set_buff_size((buff_size,buff_size))  #for resolution increasing
       sz_dens.set_zlim( field, 1.0e-2, 1e6 )
       sz_dens.set_cmap( field, colormap_dens )
       sz_dens.annotate_scale( pos=(0.1,0.9), coeff=5., corner='upper_left', unit = 'kpc', text_args={'size':16.0}, size_bar_args={'color':'black'})
       sz_dens.save( "Fig%06d_Slice_%s_Dens_mode_%s%s.png"%(idx,proj_axis,center_mode,suffix), mpl_kwargs={"dpi":dpi}) #for resolution increasing
