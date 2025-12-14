#**************************************************************************#
# Code:        GetColMap                                                   # 
#                                                                          #
# Author:      Wei-An Chen                                                 #
#                                                                          #
# This scipt is useful for checking the evolution of simulation by         #
# plotting the column density map.                                         #
#                                                                          #
# s_index:   the starting index of the snapshot,                           #
#            e.g. s_index = 10 for Data_000010.                            #
# e_index:   the end index of the snapshot.                                #
# npixel:    pixel number of the map at each dimension.                    #
# ProjPlane: the projection plane [xy/xz].                                 #
#**************************************************************************#

import yt
from astropy.constants import m_p
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import matplotlib as mpl

mpl.rcParams.update({'font.size': 15})
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['patch.linewidth'] = 2
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.width'] = 2

s_index   = 0
e_index   = 300
npixel    = 256
ProjPlane = 'xy'
Output    = "./ColMap"
tff       = 16.560403596485347 # kyr

def ProjCol(ds, c, normal_vector, north_vector, width, npixel):
    field = ('gamer', 'Dens')
    print("Projecting Dens...")
    Col = yt.off_axis_projection(ds, center=c,  normal_vector=normal_vector, north_vector=north_vector, width=width,
                          resolution=npixel, item=field, 
                          no_ghost=False).to('g/cm**2').value
    Col /= 2.3*m_p.cgs.value
    
    return Col.T

fns = ['Data_{}'.format(format(index_id,'06')) for index_id in range(s_index, e_index+1)]
index_list = [format(index_id,'06') for index_id in range(s_index, e_index+1)]

if not os.path.exists(Output):
      os.mkdir(Output)

TimeArr    = []
Par1MArr   = []
Par2MArr   = []
GasMassArr = []

for fn, index_id in zip(fns, index_list):
    try: 
        # Since some snapshots may be deleted, use "try" to skip them
        ds = yt.load(fn)
        u = ds.units
        ad = ds.all_data()
    except:
        continue
        
    print("Current snapshot = %s"%fn)
    print('Simulation time = %.2f kyr'%(ds.current_time.to('kyr')))
    
    # Read the box dimensions
    box_size_lxy = ds.domain_right_edge.value[0] # cm, for x and y
    box_size_lz  = ds.domain_right_edge.value[2] # cm, for z
    
    # Projection box size and center
    width = ([box_size_lxy, box_size_lxy, box_size_lxy]*u.cm)
    c = ds.domain_center.value*u.cm
    
    print("Center = ", c.to("pc"))
    print("Box size for projection = ", width.to("pc"))
    
    # Set normal and north (up) vector
    if ProjPlane == 'xy':
        normal_vector = np.array([0., 0., 1.]) # the direction normal to the plane of sky, which is toward the observer and must be normalized
        north_vector = np.array([0., 1., 0.])

    elif ProjPlane == 'xz':
        normal_vector = np.array([0., -1., 0.])
        north_vector  = np.array([-1., 0., 0.])
        
    img_extend = width.to('pc').value[0]/2 # the boundary for plotting

    # Information for the gas
    Time  = ds.current_time.to('kyr').value/tff
    Gas_M = (ad["gamer", "Dens"]*ad["gamer", "cell_volume"]).to('Msun').value
    
    # Information for sink particle
    Par_M  = ad['all', 'ParMass'].to('Msun').value
    Par_ID = ad['all', 'PAR_ID'].value
    Par_x  = (ad['all', 'ParPosX'].to('pc') - c[0].to('pc')).value
    Par_y  = (ad['all', 'ParPosY'].to('pc') - c[1].to('pc')).value
    Par_z  = (ad['all', 'ParPosZ'].to('pc') - c[2].to('pc')).value
    
    TimeArr.append(Time)
    if Par_M.shape[0] > 0:
        print("Particle number = ", Par_M.shape[0])
        print("Total sink mass = ", np.sum(Par_M), " M_sun")
        print("Most massive sink mass = ", np.max(Par_M), " M_sun")

        M_tot = np.sum(Gas_M) + np.sum(Par_M)
        print("Total mass = %.13e M_sun"%M_tot)

        Par1MArr.append(Par_M[Par_ID==0][0])
        Par2MArr.append(Par_M[Par_ID==1][0])
    
    else:
        print("Total mass = %.13e M_sun"%np.sum(Gas_M))

        Par1MArr.append(0)
        Par2MArr.append(0)

    GasMassArr.append(np.sum(Gas_M))
    
    # Projection
    Col = ProjCol(ds, c, normal_vector, north_vector, width, npixel)
    
    # Plotting
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    im = ax.imshow(Col, origin='lower', cmap='gist_heat', norm=colors.LogNorm(vmin=1e21, vmax=1e27), 
                 extent=[-img_extend, img_extend, -img_extend, img_extend])
    
    img_limt = 0.01
    ax.set_xlim(-img_limt, img_limt)
    ax.set_ylim(-img_limt, img_limt)

    if ProjPlane == 'xy':
        ax.scatter(Par_x, Par_y, color='b', s=1)
    elif ProjPlane == 'xz':
        ax.scatter(Par_z, -Par_x, color='b', s=1)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, aspect=30)
    cbar.set_label(r'Column number density (cm$^{-2}$)', fontsize=13)

    if ProjPlane == 'xy':
        ax.set_xlabel(r'x (pc)')
        ax.set_ylabel(r'y (pc)')
    elif ProjPlane == 'xz':
        ax.set_xlabel(r'z (pc)')
        ax.set_ylabel(r'x (pc)')
        
    ax.set_title(r'%.2f $t_{ff}$'%Time)
    
    plt.savefig('%s/Column_Density_%s_%s.png'%(Output, index_id, ProjPlane), dpi=150,bbox_inches='tight')
    plt.close()

    print("Finished!\n")

# Plot sink particle mass evolution
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

ax.plot(TimeArr, Par1MArr, c="r", label='Sink Particle 1')
ax.plot(TimeArr, Par2MArr, c="g", label='Sink Particle 2')
ax.set_xlabel(r'Time ($t_{ff}$)')
ax.set_ylabel(r'Sink Particle Mass (M$_\odot$)')
ax.legend()

plt.savefig('./SinkMassEvo.png', dpi=150,bbox_inches='tight')
plt.close()

# Plot error in total mass
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

ax.plot(TimeArr, (np.array(GasMassArr) + np.array(Par1MArr) + np.array(Par2MArr) - GasMassArr[0])/GasMassArr[0], c="k")
ax.axhline(y=0, color='k', linestyle='--')

ax.set_xlabel(r'Time ($t_{ff}$)')
ax.set_ylabel(r'$((M_{\rm{gas}} + M_{\rm{par, tot}}) - M_{\rm{gas. ini}})/M_{\rm{gas. ini}}$')

plt.savefig('./MassError.png', dpi=150,bbox_inches='tight')
plt.close()