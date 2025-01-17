from pylab import *
import sys

'''
Usage: python VG_turb.py [power law index/2] [kmin] [random seed]
'''

nmodes = int(128)

grid = nmodes+1

A0 = 1

Ax = zeros((grid,grid,grid))
Ay = zeros((grid,grid,grid))
Az = zeros((grid,grid,grid))
phasex = zeros((grid,grid,grid))
phasey = zeros((grid,grid,grid))
phasez = zeros((grid,grid,grid))

n = float(sys.argv[1])
kmin = int(sys.argv[2])
sd = int(sys.argv[3])
print(n,kmin,sd)
seed(sd)

half_modes = int(nmodes/2 + 1)

# The Amplitudes of the modes are governed by a power law spectrum of the form A(k) = A0 * k**(-n)

for ii in range(0,half_modes):
	kz = ii
	print(kz)
	if(kz == 0):
		for jj in range(0,half_modes):
			ky = jj
			if(ky==0):
				for ll in range(0,half_modes):
					kx = ll
					k = sqrt(kx**2 + ky**2 + kz**2)
					if(kx==0 or k<kmin):
						Ax[ii,jj,ll] = 0
						Ay[ii,jj,ll] = 0
						Az[ii,jj,ll] = 0
					else:
						Ax[ii,jj,ll] = A0 * normal() * k**n
						Ay[ii,jj,ll] = A0 * normal() * k**n
						Az[ii,jj,ll] = A0 * normal() * k**n

					phasex[ii,jj,ll]  = 2*pi*random()
					phasey[ii,jj,ll]  = 2*pi*random()
					phasez[ii,jj,ll]  = 2*pi*random()
			if(ky>0):
				for ll in range(0,nmodes + 1):
					kx = ll
					if(ll > nmodes/2):
						kx = ll - nmodes - 1
					k = sqrt(kx**2 + ky**2 + kz**2)

					Ax[ii,jj,ll] = A0 * normal() * k**n
					Ay[ii,jj,ll] = A0 * normal() * k**n
					Az[ii,jj,ll] = A0 * normal() * k**n

					if(k<kmin):
						Ax[ii,jj,ll] = 0
						Ay[ii,jj,ll] = 0
						Az[ii,jj,ll] = 0

					phasex[ii,jj,ll]  = 2*pi*random()
					phasey[ii,jj,ll]  = 2*pi*random()
					phasez[ii,jj,ll]  = 2*pi*random()	
	if(kz>0):
		for jj in range(0,nmodes + 1):
			ky= jj
			if(jj > nmodes/2):
				ky = jj - nmodes - 1
			for ll in range(0,nmodes + 1):
				kx = ll
				if(ll > nmodes/2):
					kx = ll - nmodes - 1
				k = sqrt(kx**2 + ky**2 + kz**2)

				Ax[ii,jj,ll] = A0 * normal() * k**n
				Ay[ii,jj,ll] = A0 * normal() * k**n
				Az[ii,jj,ll] = A0 * normal() * k**n

				if(k<kmin):
					Ax[ii,jj,ll] = 0
					Ay[ii,jj,ll] = 0
					Az[ii,jj,ll] = 0

				phasex[ii,jj,ll]  = 2*pi*random()
				phasey[ii,jj,ll]  = 2*pi*random()
				phasez[ii,jj,ll]  = 2*pi*random()



Bx1 = Ax*cos(phasex)
Bx2 = 1j*Ax*sin(phasex)

vx = ifftn(Bx1 + Bx2)
vx = vx.real

By1 = Ay*cos(phasey)
By2 = 1j*Ay*sin(phasey)

vy = ifftn(By1 + By2)
vy = vy.real

Bz1 = Az*cos(phasez)
Bz2 = 1j*Az*sin(phasez)

vz = ifftn(Bz1 + Bz2)
vz = vz.real

x = linspace(-1,1,grid)
y = linspace(-1,1,grid)
z = linspace(-1,1,grid)

f=open("v_field.dat","w")

for ii in range(0,grid):
	for jj in range(0,grid):
		for ll in range(0,grid):
			f.write("{0:.5g} {1:.5g} {2:.5g} {3:.5g} {4:.5g} {5:.5g}\n".format(x[ll],  y[jj],  z[ii],  vx[ii,jj,ll], vy[ii,jj,ll], vz[ii,jj,ll]))
			
f.close()



