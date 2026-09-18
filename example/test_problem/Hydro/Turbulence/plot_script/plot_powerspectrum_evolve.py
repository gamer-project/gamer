import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess
from matplotlib.pyplot import cm
from matplotlib.ticker import LogLocator, NullFormatter
from matplotlib.patches import Rectangle


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot power spectrum evolution' )

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

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
dpi       = 150
FONT_SIZE = 24
TICK_SIZE = 20

kmax = subprocess.check_output( ["awk", '$1 == "SRC_TURB_KMAX" {print $2}', "../Record__Note"], text=True ).strip()

# Plot
fig = plt.figure(figsize=(8.5,7.0), dpi = dpi)
axs = fig.add_subplot(111)
label_rows = [ ['time'  , r'$E_{\rm kin}(k)$', r'$E_{\rm mag}(k)$'] ]
extra = Rectangle((0, 0), 1, 1.0, fc="w", fill=False, edgecolor='none', linewidth=0)
legend_handle = [ [extra, extra, extra] ]

# Load data
for idx in range(idx_start, idx_end+1, didx):
   table = np.loadtxt("EnergyPowerSpec_%06d"%idx)
   k    = table[:, 0]
   emag = table[:, 1]
   ekin = table[:, 2]

   color1 = cm.Blues (0.3+0.6*(idx-idx_start)/(idx_end-idx_start)) if idx_start != idx_end else cm.Blues (0.9)
   color2 = cm.Reds  (0.3+0.6*(idx-idx_start)/(idx_end-idx_start)) if idx_start != idx_end else cm.Reds  (0.9)

   l1, = axs.plot(k[1:], ekin[1:], color = color1, label = r'$E_{\rm kin}(k)$')
   l2, = axs.plot(k[1:], emag[1:], color = color2, label = r'$E_{\rm mag}(k)$')
   label_rows.append( ['%3.1f'%(idx*2.5), '', ''])
   legend_handle.append( [extra, l1, l2] )
legend_labels = np.array( label_rows    ).flatten('F')
legend_handle = np.array( legend_handle ).flatten('F')
axs.legend(legend_handle, legend_labels, bbox_to_anchor=(1.35, 0), loc='lower right', borderaxespad=0, fontsize=16, ncol = 3, handletextpad = -2.5, handlelength=2.5, handleheight=1.5, columnspacing=0.6, labelspacing=0.4)
axs.set_xlim(2*np.pi, k[-1])
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel(r"$k$",    fontsize=FONT_SIZE)
axs.set_ylabel(r"$E(k)$", fontsize=FONT_SIZE )
axs.axvline(float(kmax)*2*np.pi, color = '0.8', ls = '--')
plt.xticks(fontsize=TICK_SIZE)
plt.yticks(fontsize=TICK_SIZE)
plt.savefig("fig_PowerSpecEvolve.png", dpi = dpi, bbox_inches="tight", pad_inches=0.05)

plt.close()


