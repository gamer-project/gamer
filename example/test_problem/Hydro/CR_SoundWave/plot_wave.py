#!/bin/python3
import numpy as np
import matplotlib.pylab as plt

time = np.loadtxt("Record__Dump", skiprows=1)

N_dump = len(time)

print("WARNING : Assuming Passive00 to be cosmic ray energy density.")
for i in range(N_dump):
    end_sign = "\r" if i != N_dump-1 else "\n"
    print("Plotting: %06d/%06d"%(i+1, N_dump), end=end_sign)

    dens = np.loadtxt("CosmicRay_Acousticwave_Dens_%06d"%i)
    momx = np.loadtxt("CosmicRay_Acousticwave_MomX_%06d"%i)
    momy = np.loadtxt("CosmicRay_Acousticwave_MomY_%06d"%i)
    momz = np.loadtxt("CosmicRay_Acousticwave_MomZ_%06d"%i)
    pres = np.loadtxt("CosmicRay_Acousticwave_Pres_%06d"%i)
    cray = np.loadtxt("CosmicRay_Acousticwave_Passive00_%06d"%i)

    fig, ax = plt.subplots(3, 2, figsize=(15, 10))
    plt.subplots_adjust(hspace=0, wspace=0)

    ax[0, 0].set(title="Profile")
    ax[0, 0].set(ylabel="Dens")
    ax[1, 0].set(ylabel="MomX")
    ax[2, 0].set(ylabel="MomY")
    ax[0, 1].set(title="Profile")
    ax[0, 1].set(ylabel="MomZ")
    ax[1, 1].set(ylabel="Pres")
    ax[2, 1].set(ylabel="CRay")

    ax[0, 0].plot(dens[:, 0], dens[:, 2], label="Analytical")
    ax[0, 0].scatter(dens[:, 0], dens[:, 1], label="Simulation")
    ax[0, 0].legend()

    ax[1, 0].plot(momx[:, 0], momx[:, 2], label="Analytical")
    ax[1, 0].scatter(momx[:, 0], momx[:, 1], label="Simulation")
    ax[1, 0].legend()

    ax[2, 0].plot(momy[:, 0], momy[:, 2], label="Analytical")
    ax[2, 0].scatter(momy[:, 0], momy[:, 1], label="Simulation")
    ax[2, 0].legend()

    ax[0, 1].plot(momz[:, 0], momz[:, 2], label="Analytical")
    ax[0, 1].scatter(momz[:, 0], momz[:, 1], label="Simulation")
    ax[0, 1].yaxis.set_label_position("right")
    ax[0, 1].yaxis.tick_right()
    ax[0, 1].legend()

    ax[1, 1].plot(pres[:, 0], pres[:, 2], label="Analytical")
    ax[1, 1].scatter(pres[:, 0], pres[:, 1], label="Simulation")
    ax[1, 1].yaxis.set_label_position("right")
    ax[1, 1].yaxis.tick_right()
    ax[1, 1].legend()

    ax[2, 1].plot(cray[:, 0], cray[:, 2], label="Analytical")
    ax[2, 1].scatter(cray[:, 0], cray[:, 1], label="Simulation")
    ax[2, 1].yaxis.set_label_position("right")
    ax[2, 1].yaxis.tick_right()
    ax[2, 1].legend()



    plt.suptitle("t = %.2f"%time[i, 1])
    plt.savefig("CosmicRay_Acousticwave_%06d.png"%i)
    plt.close()
