import numpy as np
import h5py
import yt
import argparse
import sys
import gc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Recreate turbulence acceleration field from snapshot' )

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

dpi       = 150
fontsize  = 24
titlepad  = 16


yt.enable_parallelism()

ts = yt.DatasetSeries( [ '../Data_%06d' % idx for idx in range(idx_start, idx_end + 1, didx) ] )

for ds in ts.piter():
   filename = ds.parameter_filename

   f = h5py.File(filename, "r")
   NMode     = f["Turbulence/NMode"    ][()]
   TimeLast  = f["Turbulence/TimeLast" ][()]
   TimeNext  = f["Turbulence/TimeNext" ][()]
   OUArrLast = np.array( f["Turbulence/OUArrLast"], dtype=np.float64 )
   OUArrNext = np.array( f["Turbulence/OUArrNext"], dtype=np.float64 )
   OUArr     = np.stack([ OUArrLast.reshape(NMode, 3, 2), OUArrNext.reshape(NMode, 3, 2) ])
   Time      = np.array([ TimeLast,  TimeNext ])
   f.close()

   L    = ds.parameters["BoxSize"            ][0]
   vel  = ds.parameters["Src_Turb_Vel"       ]
   ampl = ds.parameters["Src_Turb_AmplFactor"]
   kdri = ds.parameters["Src_Turb_Kdriv"     ]
   kmin = ds.parameters["Src_Turb_Kmin"      ]
   kmax = ds.parameters["Src_Turb_Kmax"      ]
   zeta = ds.parameters["Src_Turb_Zeta"      ]
   form = ds.parameters["Src_Turb_SpecForm"  ]
   N    = ds.parameters["Src_Turb_TableSize" ]
   dh   = L/N

   zeta_norm = ( 3.0 / (1.0 - 2.0*zeta + 3.0*zeta**2) )**0.5
   tau = L / kdri / vel
   var = ( ( ampl*0.15*vel )**3 / L / tau )**0.5

   kmode1d = 2*np.pi/L * np.arange( -int(kmax), int(kmax) + 1)

   kmin *= 2*np.pi/L
   kmax *= 2*np.pi/L
   kmid  = 0.5*(kmin + kmax)

   kx = np.broadcast_to( kmode1d[None, None, :], (len(kmode1d), len(kmode1d), len(kmode1d)) )
   ky = np.broadcast_to( kmode1d[None, :, None], (len(kmode1d), len(kmode1d), len(kmode1d)) )
   kz = np.broadcast_to( kmode1d[:, None, None], (len(kmode1d), len(kmode1d), len(kmode1d)) )

   kmag = np.sqrt(kx**2 + ky**2 + kz**2)

   mask = (kmag >= kmin) & (kmag <= kmax)

   kmag  = kmag[mask]
   kmode = np.vstack(( kx[mask], ky[mask], kz[mask] ))

   nmode = len(kmag)

   if nmode != NMode:
      print( "number of modes computed from Input__Parameter ( %d )  does not match stored number ( %d )\n"%(nmode, NMode) )
      sys.exit(1)

   if form == 0:
      amplitude = kmin/kmag
   elif form == 1:
      amplitude = np.sqrt( np.abs( -4 * ((kmag - kmid) / (kmax - kmin))**2 + 1) )*kmid/kmag
   elif form == 2:
      power     = ds.parameters["Src_Turb_Pow"]
      amplitude = np.sqrt( (kmag / kmin)**power )*kmin/kmag
   else:
      print( "unknown SRC_TURB_SPEC_FORM\n" )
      sys.exit(1)

   amplitude *= 2*zeta_norm

   x = (np.arange(N) + 0.5) * dh

#  phase[nmode, k, j, i] = k dot x
   phase = ( kmode[0, :, None, None, None] * x[None, None, None, :]
           + kmode[1, :, None, None, None] * x[None, None, :, None]
           + kmode[2, :, None, None, None] * x[None, :, None, None] )

   real = np.cos(phase)
   imag = np.sin(phase)

   del kx, ky, kz, phase
   gc.collect()

   ds_acc = [None]*2
   bbox = np.array([[0.0, L], [0.0, L], [0.0, L]])

   for i in range(2):
      Acc = np.sum( amplitude[:, None, None, None, None]
                 * ( OUArr[i, :, :, 0, None, None, None] * real[:, None, :, :, :]
                   - OUArr[i, :, :, 1, None, None, None] * imag[:, None, :, :, :] ), axis=0 )

      AccMag = np.sqrt( Acc[0]**2 + Acc[1]**2 + Acc[2]**2 )

      data = { ("gas", "AccMag"): AccMag }

      ds_acc[i] = yt.load_uniform_grid( data, domain_dimensions=(N, N, N), time_unit='s', sim_time=Time[i],
                                        bbox=bbox, axis_order=("z", "y", "x") )
      print( "" )
      print( "   time = %13.7e"%Time[i]        )
      print( "   min  = %13.7e"%np.min(AccMag) )
      print( "   max  = %13.7e"%np.max(AccMag) )
      print( "   std  = %13.7e"%np.std(AccMag) )
      print( "")

   del Acc, AccMag, real, imag
   gc.collect()

#  plot
   fig = plt.figure()
   fig.dpi = dpi
   grid = AxesGrid( fig, (0.1, 0.05, 1.8, 1.7), nrows_ncols=(1, 2), axes_pad=(1.5,0.5), label_mode="all",
                          share_all=True, cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="2%" )

   slc = [None]*2
   for i in range(2):
      fieldname = "AccMag"
      slc[i] = yt.SlicePlot( ds_acc[i], 0, fields = fieldname, center = 'c')
      slc[i].set_background_color( fieldname )
      slc[i].set_axes_unit( 'code_length' )
      slc[i].set_zlim( fieldname, 5e-4, 5e-1 )
      slc[i].set_cmap( fieldname, 'viridis' )
      slc[i].set_font( {'size':fontsize} )
      slc[i].annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )

      plot = slc[i].plots[fieldname]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[i]
      slc[i]._setup_plots()
      grid[i].set_title(fieldname, fontsize=fontsize, pad=titlepad)

   fig.savefig("fig_%s_AccTable.png"%(ds), bbox_inches='tight',pad_inches=0.02)


