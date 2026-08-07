import numpy as np
import argparse

'''
Usage: python generate_TubulentVelocityField.py -n [power law index] -kmin [kmin] -seed [random seed] -nmodes [number of modes]
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
vel_n      = args.n/2   # the power law index of |\hat{v}_k|, the AMPLTUDE of the k-space velocity field
kmin       = args.kmin  # minimum of k
seed       = args.seed  # random seed
nmodes     = args.nmodes

# set the number of modes
grid       = nmodes

# compute the k-space coordinates
k1D_index  = np.fft.fftfreq(grid) * grid
k1D_r      = np.fft.rfftfreq(grid) * grid
kz, ky, kx = np.meshgrid( k1D_index, k1D_index, k1D_r, indexing='ij' )
k          = np.sqrt( kx**2 + ky**2 + kz**2 )

# get the random numbers following a specific order
np.random.seed(seed)
rfft_shape = kx.shape
sampledRandomNumbers_Normals = np.random.standard_normal((3, *rfft_shape))
sampledRandomNumbers_Randoms = np.random.uniform(0, 1, size=(3, *rfft_shape))

# compute the amplitudes, a power law spectrum of the form A(k) = A0 * k**vel_n
A0    = 1.0
A     = np.zeros((3, *rfft_shape)) # initialize A
highK = (k >= kmin)                # high-k filter

A[0][highK] = A0 * sampledRandomNumbers_Normals[0][highK] * k[highK]**vel_n
A[1][highK] = A0 * sampledRandomNumbers_Normals[1][highK] * k[highK]**vel_n
A[2][highK] = A0 * sampledRandomNumbers_Normals[2][highK] * k[highK]**vel_n

# compute the phases
phase = 2 * np.pi * sampledRandomNumbers_Randoms

# compute the velocities using irfftn to enforce Hermitian symmetry
vx = np.fft.irfftn( A[0] * np.exp(1j * phase[0]), s=(grid, grid, grid) )
vy = np.fft.irfftn( A[1] * np.exp(1j * phase[1]), s=(grid, grid, grid) )
vz = np.fft.irfftn( A[2] * np.exp(1j * phase[2]), s=(grid, grid, grid) )

# compute the x-space coordinates
z, y, x      = np.meshgrid( np.linspace(-1, 1, grid, endpoint=False), np.linspace(-1, 1, grid, endpoint=False), np.linspace(-1, 1, grid, endpoint=False), indexing='ij')

# save to the file
np.savetxt( "Tur_Table.dat",
            np.column_stack( ( x.ravel(), y.ravel(), z.ravel(), vx.ravel(), vy.ravel(), vz.ravel() ) ),
            fmt='%12.4e', delimiter=' ', header='%10s%13s%13s%13s%13s%13s'%("x","y","z","vx","vy","vz") )