import numpy as np
import scipy
import sys
import argparse

# load the command-line parameters
parser = argparse.ArgumentParser( description='Upscale periodic GAMER ELBDM initial conditions' )

parser.add_argument( '-float8', action='store', required=False,  type=bool, dest='float8',
                            help='assume double precision for input and output file', default=False )
parser.add_argument( '-n_in', action='store', required=True, type=int, dest='n_in',
                            help='input file resolution')
parser.add_argument( '-n_out', action='store', required=True, type=int, dest='n_out',
                            help='output file resolution')
parser.add_argument( '-n_thread', action='store', required=False, type=int, dest='n_thread',
                            help='number of threads to use', default=-1)
parser.add_argument( '-input', action='store', required=False,  type=str, dest='input',
                      help='input file [%(default)s]', default='./UM_IC_lr' )
parser.add_argument( '-output', action='store', required=False,  type=str, dest='output',
                      help='output file [%(default)s]', default='./UM_IC' )


###
args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
    print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )

# set up interpolation parameters
input_file  = args.input
output_file = args.output
N_in        = args.n_in
N_out       = args.n_out
NThreads    = args.n_thread         # -1 wraps around from os.cpu_count()
FLOAT_8     = False           # enable double precision

if FLOAT_8:
     rprec = np.double
     cprec = np.cdouble
else:
     rprec = np.float32
     cprec = np.complex64


def interp(psi, N_in, N_out, NThreads):
     # compute fft of array
     print("Computing FFT... ", end="")
     psihat = scipy.fft.fftn(psi, workers = NThreads, norm="forward")
     print("done!")

     # free input arrays to save memory
     del psi

     # shift zero frequencies to center of cube
     psihat = np.fft.fftshift(psihat)

     # pad cube from the outside
     N_pad  = int(np.floor(N_out/2 - N_in/2))
     psihat = np.pad(psihat, N_pad, mode="constant")

     # shift zero frequencies back to outside of cube
     print("Computing IFFT... ", end="")
     psihat = np.fft.fftshift(psihat)
     print("done!")

     # inverse FFT
     psihr  = scipy.fft.ifftn(psihat, workers = NThreads, norm="forward")

     # delete large array in frequency domain to save memory
     del psihat

     return psihr


# load binary from file
print("Loading input data... ", end="")
lr   = np.fromfile(input_file, dtype=rprec).reshape((2, N_in, N_in, N_in))
psi  = lr[0, :, :, : ] + 1j * lr[1, :, :, : ]
psi  = psi.astype(cprec)
del lr
print("done!")

# interpolate data
psihr  = interp(psi, N_in, N_out, NThreads).astype(cprec)


# write data to to disk
print("Writing wave function to binary file... ", end="")

# overwrite output file with real part
with open(output_file,"wb") as f:
     np.real(psihr).astype(rprec).tofile(f)
     f.close()

# append imaginary part
with open(output_file,"a") as f:
     np.imag(psihr).astype(rprec).tofile(f)
     f.close()

print("done!")
