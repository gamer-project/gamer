#!/bin/python3
import sys
import numpy as np
import matplotlib.pylab as plt
import yt

sys.dont_write_bytecode = True
import shocktube_ana as shock_ana


def get_data(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
        name  = lines[0].split()[1:]
        lines = lines[1:]
        for i in range(len(lines)):
            lines[i] = lines[i].split()
            for j in range(len(lines[i])):
                lines[i][j] = float(lines[i][j])
    data = np.array(lines)
    return data


#===============================
# GAMER data
#===============================
SHOCK_DIR = 0   # shock tube direction [0/1/2]

if SHOCK_DIR == 0:
    gamer_file = './Xline_y0.000_z0.000_000005'
elif SHOCK_DIR == 1:
    gamer_file = './Yline_x0.000_z0.000_000005'
elif SHOCK_DIR == 2:
    gamer_file = './Zline_x0.000_y0.000_000005'
else:
    raise ValueError( "Wrong SHOCK_DIR (%d)"%SHOCK_DIR )

gamer_data = get_data(gamer_file)

#['i', 'j', 'k', 'x', 'y', 'z', 'Dens', 'MomX', 'MomY', 'MomZ', 'Engy', 'CRay', 'MagX', 'MagY', 'MagZ', 'MagEngy']
gamer_x    = gamer_data[:, 3]
gamer_y    = gamer_data[:, 4]
gamer_z    = gamer_data[:, 5]
gamer_rho  = gamer_data[:, 6]
gamer_vx   = gamer_data[:, 7] / gamer_data[:, 6]
gamer_vy   = gamer_data[:, 8] / gamer_data[:, 6]
gamer_vz   = gamer_data[:, 9] / gamer_data[:, 6]
gamer_Ek   = 0.5 * gamer_rho * (gamer_vx**2 + gamer_vy**2 + gamer_vz**2)
gamer_engy = gamer_data[:, 10]
gamer_cray = gamer_data[:, 11]
gamer_pres = (gamer_engy - gamer_Ek - gamer_cray) * (5./3. - 1.)

if SHOCK_DIR == 0:
    gamer_r  = gamer_x
    gamer_vr = gamer_vx
elif SHOCK_DIR == 1:
    gamer_r  = gamer_y
    gamer_vr = gamer_vy
elif SHOCK_DIR == 2:
    gamer_r  = gamer_z
    gamer_vr = gamer_vz


#===============================
# analytical solution
#===============================
ana_rho, ana_pres_cr, ana_pres, ana_velx = shock_ana.shock_sol(1.0, 0.2, 1.3e+5, 2.4e+2, 6.7e+4, 2.4e+2, 0., 0., gamer_r, 0.5, 4.4e-4)
ana_cray = ana_pres_cr / (4./3. - 1.)



#===============================
# plot
#===============================
fig, ax = plt.subplots(2, 2, figsize = (8, 8))
plt.subplots_adjust(hspace=0, wspace=0)

ax[0, 0].set(ylabel = r'$\rho$')
ax[0, 0].plot(gamer_r, ana_rho, label = 'Analytical')
ax[0, 0].plot(gamer_r, gamer_rho, 'o', ms = 3, label = 'GAMER')
ax[0, 0].legend(loc = 0)

ax[0, 1].set(ylabel = r'$v$')
ax[0, 1].plot(gamer_r, ana_velx, label = 'Analytical')
ax[0, 1].plot(gamer_r, gamer_vr, 'o', ms = 3, label = 'GAMER')
ax[0, 1].yaxis.set_label_position("right")
ax[0, 1].yaxis.tick_right()
ax[0, 1].legend(loc = 0)

ax[1, 0].set(ylabel = r'$P/(\gamma-1)$')
ax[1, 0].plot(gamer_r, ana_pres/(5./3. - 1.), label = 'Analytical')
ax[1, 0].plot(gamer_r, gamer_pres/(5./3. - 1.), 'o', ms = 3, label = 'GAMER')
ax[1, 0].ticklabel_format(axis = 'y', style = 'sci', scilimits = (4, 5))
ax[1, 0].legend(loc = 0)

ax[1, 1].set(ylabel = r'$e_{cr}$')
ax[1, 1].plot(gamer_r, ana_cray, label = 'Analytical')
ax[1, 1].plot(gamer_r, gamer_cray, 'o', ms = 3, label = 'GAMER')
ax[1, 1].ticklabel_format(axis = 'y', style = 'sci', scilimits = (4, 5))
ax[1, 1].yaxis.set_label_position("right")
ax[1, 1].yaxis.tick_right()
ax[1, 1].legend(loc = 0)

plt.savefig('Shocktube.png')
plt.show()
