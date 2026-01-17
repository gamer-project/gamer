import numpy as np
import argparse

'''
Usage: python generate_TubulentVelocityField.py -n [power law index] -kmin [kmin] -seed [random seed] -nmodes [number of modes]
Author: Hsinhao Huang and S. D. Clarke
Reference: O. Lomax et al., 2015, MNRAS, 449, 662
'''

# load the command-line parameters
parser = argparse.ArgumentParser( description='Generation of turbulent velocity field' )

parser.add_argument( '-n', action='store', required=False, type=float, dest='n',
                     help='power law index [%(default)s]', default=-4 )
parser.add_argument( '-kmin', action='store', required=False, type=int, dest='kmin',
                     help='minimum of k [%(default)d]', default=1 )
parser.add_argument( '-seed', action='store', required=False, type=int, dest='seed',
                     help='random number seed [%(default)d]', default=0 )
parser.add_argument( '-nmodes', action='store', required=False, type=int, dest='nmodes',
                     help='number of modes [%(default)d]', default=128 )

args=parser.parse_args()

# set the input parameters
n          = args.n/2   # this is the power law index of |\hat{v}_k|, the AMPLTUDE of the k-space velocity field
kmin       = args.kmin  # minimum of k
seed       = args.seed  # random seed
nmodes     = args.nmodes

# set the number of modes
grid       = nmodes+1

# compute the k-space coordinates
k1D_index  = np.arange(grid)
k1D        = np.where( k1D_index > nmodes/2, k1D_index-nmodes-1, k1D_index ) # if index > nmodes/2, map to negative frequencies
kz, ky, kx = np.meshgrid( k1D, k1D, k1D, indexing='ij' )
k          = np.sqrt( kx**2 + ky**2 + kz**2 )

# get the random numbers following a specific order
np.random.seed(seed)
sampledRandomNumbers_Normals = np.zeros((3,grid,grid,grid))
sampledRandomNumbers_Randoms = np.zeros((3,grid,grid,grid))
getNormals = ~( (kz < 0) | ((kz == 0) & (ky < 0)) | ((kz == 0) & (ky == 0) & (kx < 0)) | ((kz == 0) & (ky == 0) & (k < kmin)) ).ravel()
getRandoms = ~( (kz < 0) | ((kz == 0) & (ky < 0)) | ((kz == 0) & (ky == 0) & (kx < 0)) ).ravel()
for i in range(0, grid**3):
    sampledRandomNumbers_Normals[0,].ravel()[i] = np.random.normal() if getNormals[i] else 0
    sampledRandomNumbers_Normals[1,].ravel()[i] = np.random.normal() if getNormals[i] else 0
    sampledRandomNumbers_Normals[2,].ravel()[i] = np.random.normal() if getNormals[i] else 0
    sampledRandomNumbers_Randoms[0,].ravel()[i] = np.random.random() if getRandoms[i] else 0
    sampledRandomNumbers_Randoms[1,].ravel()[i] = np.random.random() if getRandoms[i] else 0
    sampledRandomNumbers_Randoms[2,].ravel()[i] = np.random.random() if getRandoms[i] else 0

# compute the amplitudes, a power law spectrum of the form A(k) = A0 * k**n
A0           = 1
A            = np.zeros((3,grid,grid,grid)) # initialize A
highK        = (k >= kmin)                  # high-k filter
A[0,][highK] = A0 * sampledRandomNumbers_Normals[0,][highK] * k[highK]**n
A[1,][highK] = A0 * sampledRandomNumbers_Normals[1,][highK] * k[highK]**n
A[2,][highK] = A0 * sampledRandomNumbers_Normals[2,][highK] * k[highK]**n

# compute the phases
phase        = 2 * np.pi * sampledRandomNumbers_Randoms

# compute the velocities
vx           = np.fft.ifftn( A[0,]*( np.cos(phase[0,]) + 1j*np.sin(phase[0,]) ) ).real
vy           = np.fft.ifftn( A[1,]*( np.cos(phase[1,]) + 1j*np.sin(phase[1,]) ) ).real
vz           = np.fft.ifftn( A[2,]*( np.cos(phase[2,]) + 1j*np.sin(phase[2,]) ) ).real

# compute the x-space coordinates
z, y, x      = np.meshgrid( np.linspace(-1, 1, grid), np.linspace(-1, 1, grid), np.linspace(-1, 1, grid), indexing='ij')

# save to the file
np.savetxt( "Tur_Table.dat",
            np.column_stack( ( x.ravel(), y.ravel(), z.ravel(), vx.ravel(), vy.ravel(), vz.ravel() ) ),
            fmt='%12.4e', delimiter=' ', header='%10s%13s%13s%13s%13s%13s'%("x","y","z","vx","vy","vz") )