import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Projection of mass density' )

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

colormap_dens = 'algae'
dpi           = 150

yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
for ds in ts.piter():

   fields_list = ['particle_density_on_grid']

   for center_mode in ['c']:
      for direction in ['z']:
         for field in fields_list:

            # decide the center
            if center_mode == 'm':
               center = ds.all_data().quantities.max_location(field)[1:]
            elif center_mode == 'c':
               center = ds.domain_center

            # ProjectionPlot
            p_dens = yt.ProjectionPlot( ds, direction, field, center=center, buff_size=(1024, 1024) )

            # setting for the figure
            p_dens.set_axes_unit( 'code_length' )
            p_dens.set_unit( field, 'code_mass/code_length**2' )
            p_dens.set_zlim( field, 1.0e-5, 1.0e0 )
            p_dens.set_cmap( field, colormap_dens )
            p_dens.set_background_color( field )
            p_dens.annotate_timestamp( time_unit='code_time', corner='upper_right' )

            # zoom in
            if center_mode == 'm':
               p_dens.zoom(4)

            # save the figure
            p_dens.save( 'fig_%s_%s_Projection_%s_%s.png'%(ds, center_mode, direction, field), mpl_kwargs={'dpi':dpi} )

            # annotate the grids and save again
            p_dens.annotate_grids()
            p_dens.save( 'fig_%s_%s_Projection_%s_%s_withgrids.png'%(ds, center_mode, direction, field), mpl_kwargs={'dpi':dpi} )

            # additionally, plot the particles
            if field == 'particle_density_on_grid':
               # decide the direction
               if direction == 'z':
                  p_par = yt.ParticlePlot( ds, 'particle_position_x', 'particle_position_y', 'particle_mass', center=center )
               elif direction == 'y':
                  p_par = yt.ParticlePlot( ds, 'particle_position_z', 'particle_position_x', 'particle_mass', center=center )

               # setting for the figure
               p_par.set_axes_unit( 'code_length' )
               p_par.set_unit( 'particle_mass', 'code_mass'    )
               p_par.set_zlim( 'particle_mass', 1.0e-8, 1.0e-4 )
               p_par.set_cmap( 'particle_mass', colormap_dens  )
               p_par.set_background_color('particle_mass' )
               p_par.annotate_timestamp( time_unit='code_time', corner='upper_right' )

               # zoom in
               if center_mode == 'm':
                  p_par.zoom(4)

               # save the figure
               p_par.save( 'fig_%s_%s_Particle_%s_particle_mass.png'%(ds, center_mode, direction), mpl_kwargs={'dpi':dpi} )
