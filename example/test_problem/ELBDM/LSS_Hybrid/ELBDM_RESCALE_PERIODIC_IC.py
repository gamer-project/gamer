#!/usr/bin/env python

import numpy as np
import scipy
import sys
import argparse

# load the command-line parameters
parser = argparse.ArgumentParser( description='Upscale periodic GAMER ELBDM initial conditions. Example usage: ./ELBDM_RESCALE_PERIODIC_IC.py -n_in 256 -n_out 64 -input UM_IC_high_resolution -output UM_IC_low_resolution' )

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

#credit to stackoverflow user Bily for this unpadding function
def unpad(x, pad_width):
    slices = []
    for c in pad_width:
        e = None if c[1] == 0 else -c[1]
        slices.append(slice(c[0], e))
    return x[tuple(slices)]

def interp(psi, N_in, N_out, NThreads):
     if N_in == N_out:
         print("N_in == N_out, no rescaling necessary!")
         return psi

     # compute fft of array
     print("Computing FFT... ")
     psihat = scipy.fft.fftn(psi, workers = NThreads, norm="forward")
     print("done!")

     # free input arrays to save memory
     del psi

     # shift zero frequencies to center of cube
     psihat = np.fft.fftshift(psihat)

     # upscale
     if N_out > N_in:
         # pad cube from the outside
         N_pad  = int(np.floor(N_out/2 - N_in/2))
         psihat = np.pad(psihat, ((N_pad, N_pad), (N_pad, N_pad), (N_pad, N_pad)), mode="constant")
     #downscale
     elif N_in > N_out:
         # unpad cube from the outside
         N_pad  = int(np.floor(N_in/2 - N_out/2))
         psihat = unpad(psihat, ((N_pad, N_pad), (N_pad, N_pad), (N_pad, N_pad)))

     # shift zero frequencies back to outside of cube
     print("Computing IFFT... ")
     psihat = np.fft.fftshift(psihat)
     print("done!")

     # inverse FFT
     psihr  = scipy.fft.ifftn(psihat, workers = NThreads, norm="forward")

     # delete large array in frequency domain to save memory
     del psihat

     return psihr


# load binary from file
print("Loading input data... ")
lr   = np.fromfile(input_file, dtype=rprec).reshape((2, N_in, N_in, N_in))
psi  = lr[0, :, :, : ] + 1j * lr[1, :, :, : ]
psi  = psi.astype(cprec)
del lr
print("done!")

# interpolate data
psihr  = interp(psi, N_in, N_out, NThreads).astype(cprec)


# write data to to disk
print("Writing wave function to binary file... ")

# overwrite output file with real part
with open(output_file,"wb") as f:
     np.real(psihr).astype(rprec).tofile(f)
     f.close()

# append imaginary part
with open(output_file,"a") as f:
     np.imag(psihr).astype(rprec).tofile(f)
     f.close()

print("done!")
