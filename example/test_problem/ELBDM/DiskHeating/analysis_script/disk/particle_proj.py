import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
from mpl_toolkits.axes_grid1 import AxesGrid
import argparse
import sys


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
FONT_SIZE   = 24.0
LINE_WIDTH  = 0.5
colormap    = 'arbre'


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the face-on and edge-on projection of stellar disk' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


# -------------------------------------------------------------------------------------------------------------------------
# specify script unit and output figure settings
cm       = 1/2.54  # centimeters in inches
field = ('Disk', 'ParMass')
Center   = np.loadtxt('../../Record__Center', skiprows=1, dtype=float)
if Center.ndim == 1:
   Center = Center.reshape(1,len(Center)) # reshape the array if there is only one row


plt.rcParams['font.size']         = FONT_SIZE
plt.rcParams['figure.titlesize']  = 2*FONT_SIZE

plt.rcParams['axes.titlesize']    = 2*FONT_SIZE
plt.rcParams['axes.labelsize']    = FONT_SIZE
plt.rcParams['axes.labelpad']     = 0.05
plt.rcParams['axes.linewidth']    = LINE_WIDTH

plt.rcParams['legend.fontsize']   = 6.0
plt.rcParams['lines.linewidth']   = LINE_WIDTH

plt.rcParams['xtick.major.size']  = 2
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.minor.size']  = 1
plt.rcParams['xtick.minor.width'] = 0.25
plt.rcParams['xtick.labelsize']   = FONT_SIZE

plt.rcParams['ytick.major.size']  = 2
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['ytick.minor.size']  = 1
plt.rcParams['ytick.minor.width'] = 0.25
plt.rcParams['ytick.labelsize']   = FONT_SIZE

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'

# add yt particle filter
def Disk(pfilter, data):
   filter = data['all', 'ParType'] == 2
   return filter
yt.add_particle_filter('Disk', function = Disk, filtered_type = 'all', requires = ['ParType'])

# -------------------------------------------------------------------------------------------------------------------------
# output figures
fig = plt.figure()
fig.dpi = 150

grid = AxesGrid( fig, (0.1, 0.05, 3.2, 2.7), nrows_ncols=(1, 2), axes_pad=(1.2,0.5), label_mode="all", share_all=True, cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="2%")

for idx in range(idx_start, idx_end+1, didx):
   ds             = yt.load( '../../Data_%06d'%idx )
   ds.add_particle_filter('Disk')
   if sys.version_info[0] == 2:
      ds.periodicity=(True,True,True)
   current_step   = ds.parameters["Step"]
   print("Current Simulation Time = %.5e [code units]"%ds.parameters["Time"][0])
   print("Current Simulation Step = %i"%current_step)

   parx = yt.ParticleProjectionPlot( ds, 'x', fields = field, center = Center[current_step,3:6], width =( (60, 'kpc'),(60,'kpc')))
   parx.set_background_color( field )
   parx.set_cmap( field, colormap )
   parx.set_colorbar_label(field, "Projected stellar mass (M$_\odot$)")
   parx.set_font( {'size':FONT_SIZE} )
   parx.set_axes_unit( 'kpc' )
   parx.set_unit( field, 'Msun' )
   parx.set_zlim( field, 5.5e+5, 5.5e+2, dynamic_range=None)
   parx.annotate_text([ 0.05, 0.92], 'edge-on', coord_system="axis",text_args={"size":FONT_SIZE,"color":"black","weight":"normal","bbox":dict(boxstyle="round",ec='white',fc='white',alpha=0.7)})
   plot = parx.plots[field]
   plot.figure = fig
   plot.axes = grid[0].axes
   plot.cax = grid.cbar_axes[0]
   parx._setup_plots()

   parz = yt.ParticleProjectionPlot( ds, 'z', fields = field, center = Center[current_step,3:6], width =( (60, 'kpc'),(60,'kpc')))
   parz.set_background_color( field )
   parz.set_cmap( field, colormap )
   parz.set_colorbar_label(field, "Projected stellar mass (M$_\odot$)")
   parz.set_font( {'size':FONT_SIZE} )
   parz.set_axes_unit( 'kpc' )
   parz.set_unit( field, 'Msun' )
   parz.set_zlim( field, 5.5e+5, 5.5e+2, dynamic_range=None)
   parz.annotate_text([ 0.05, 0.92], 'face-on', coord_system="axis",text_args={"size":FONT_SIZE,"color":"black","weight":"normal","bbox":dict(boxstyle="round",ec='white',fc='white',alpha=0.7)})
   parz.annotate_text([ 0.75, 0.92], r'$t_{\rm rel}$ = %2.1f Gyr'%(140.59*idx/1000.), coord_system="axis", text_args={"size":FONT_SIZE,"color":"white"})
   plot = parz.plots[field]
   plot.figure = fig
   plot.axes = grid[1].axes
   plot.cax = grid.cbar_axes[1]
   parz._setup_plots()

   fig.set_size_inches(18*cm, 8*cm)
   fig.savefig("particle_proj_%06d.png"%idx, bbox_inches='tight',pad_inches=0.02)
   #fig.savefig("particle_proj_%06d.pdf"%idx, bbox_inches='tight',pad_inches=0.02)

   print('\nparticle_proj_%06d.png completed\n'%idx)
