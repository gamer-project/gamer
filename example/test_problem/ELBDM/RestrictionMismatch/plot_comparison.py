import yt
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable


def make_1d_continuous(f):
    for i in range(len(f) - 1):
        while (f[i] - f[i + 1]) > np.pi:
            f[i + 1 :] += 2 * np.pi
        while (f[i] - f[i + 1]) < -np.pi:
            f[i + 1 :] -= 2 * np.pi
    return f


def make_2d_continuous(f):
    for i in range(f.shape[0]):
        make_1d_continuous(f[i, :])
    for i in range(f.shape[1]):
        make_1d_continuous(f[:, i])
    return f


def getLaplacian(field):
    return (np.abs(np.roll(field,-1, 0) + np.roll(field,+1, 0) +np.roll(field,-1, 1) +np.roll(field,+1, 1) - 4*field))[3:-3, 3:-3]


def getSlice(field):
    return field[:, :, int(field.shape[2]/2)]


ds1  = yt.load("Data_000001_NoRefinement")
ds2  = yt.load("Data_000001_OldRestriction")
ds3  = yt.load("Data_000001_NewRestriction")

level = 0
index = 1
axis  = 2

d1 = ds1.covering_grid(
            level=level, left_edge=[0, 0.0, 0.0], dims=ds1.domain_dimensions * 2 ** level
            )

d2 = ds2.covering_grid(
            level=level, left_edge=[0, 0.0, 0.0], dims=ds2.domain_dimensions * 2 ** level
            )

d3 = ds3.covering_grid(
            level=level, left_edge=[0, 0.0, 0.0], dims=ds3.domain_dimensions * 2 ** level
            )



phase1 = getSlice(np.arctan2(d1["gamer", "Imag"], d1["gamer", "Real"]))
phase2 = getSlice(np.arctan2(d2["gamer", "Imag"], d2["gamer", "Real"]))
phase3 = getSlice(np.arctan2(d3["gamer", "Imag"], d3["gamer", "Real"]))
phase1 = make_2d_continuous(phase1)
phase2 = make_2d_continuous(phase2)
phase3 = make_2d_continuous(phase3)

fig, axes = plt.subplots(2, 3, dpi = 200, figsize=(18, 12))
ax = axes.reshape(6)
fig.suptitle("Comparison of restriction methods (old/new) at timestep %d on level %d" % (index, level))

ax[0].set_title("Phase (old)")
im1 = ax[0].imshow(phase2)

divider = make_axes_locatable(ax[0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')

ax[1].set_title("Laplacian of phase (old)")
im2 = ax[1].imshow(getLaplacian(phase2))

divider = make_axes_locatable(ax[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical')


ax[2].set_title("Difference (With/Without refinement)")
im3 = ax[2].imshow(np.abs(phase1-phase2))

divider = make_axes_locatable(ax[2])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im3, cax=cax, orientation='vertical')

ax[3].set_title("Phase (new)")
im4 = ax[3].imshow(phase3)

divider = make_axes_locatable(ax[3])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im4, cax=cax, orientation='vertical')

ax[4].set_title("Laplacian of phase (new)")
im5 = ax[4].imshow(getLaplacian(phase3))

divider = make_axes_locatable(ax[4])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im5, cax=cax, orientation='vertical')


ax[5].set_title("Difference (With/Without refinement)")
im6 = ax[5].imshow(np.abs(phase1-phase3))

divider = make_axes_locatable(ax[5])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im6, cax=cax, orientation='vertical')


plt.savefig("ComparisonOfRestrictionMethodsBeforeEvolution.png")
plt.close()

ds1  = yt.load("Data_000002_NoRefinement")
ds2  = yt.load("Data_000002_OldRestriction")
ds3  = yt.load("Data_000002_NewRestriction")

level = 0
index = 2
axis  = 2

d1 = ds1.covering_grid(
            level=level, left_edge=[0, 0.0, 0.0], dims=ds1.domain_dimensions * 2 ** level
            )

d2 = ds2.covering_grid(
            level=level, left_edge=[0, 0.0, 0.0], dims=ds2.domain_dimensions * 2 ** level
            )

d3 = ds3.covering_grid(
            level=level, left_edge=[0, 0.0, 0.0], dims=ds3.domain_dimensions * 2 ** level
            )



phase1 = getSlice(np.arctan2(d1["gamer", "Imag"], d1["gamer", "Real"]))
phase2 = getSlice(np.arctan2(d2["gamer", "Imag"], d2["gamer", "Real"]))
phase3 = getSlice(np.arctan2(d3["gamer", "Imag"], d3["gamer", "Real"]))
phase1 = make_2d_continuous(phase1)
phase2 = make_2d_continuous(phase2)
phase3 = make_2d_continuous(phase3)

fig, axes = plt.subplots(2, 3, dpi = 200, figsize=(18, 12))
ax = axes.reshape(6)
fig.suptitle("Comparison of restriction methods (old/new) at timestep %d on level %d" % (index, level))

ax[0].set_title("Phase (old)")
im1 = ax[0].imshow(phase2)

divider = make_axes_locatable(ax[0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')

ax[1].set_title("Laplacian of phase (old)")
im2 = ax[1].imshow(getLaplacian(phase2))

divider = make_axes_locatable(ax[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical')


ax[2].set_title("Difference (With/Without refinement)")
im3 = ax[2].imshow(np.abs(phase1-phase2))

divider = make_axes_locatable(ax[2])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im3, cax=cax, orientation='vertical')

ax[3].set_title("Phase (new)")
im4 = ax[3].imshow(phase3)

divider = make_axes_locatable(ax[3])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im4, cax=cax, orientation='vertical')

ax[4].set_title("Laplacian of phase (new)")
im5 = ax[4].imshow(getLaplacian(phase3))

divider = make_axes_locatable(ax[4])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im5, cax=cax, orientation='vertical')


ax[5].set_title("Difference (With/Without refinement)")
im6 = ax[5].imshow(np.abs(phase1-phase3))

divider = make_axes_locatable(ax[5])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im6, cax=cax, orientation='vertical')


plt.savefig("ComparisonOfRestrictionMethodsAfterEvolution.png")
plt.close()


