from __future__ import print_function, division  # Ensure Python 2/3 compatibility

import numpy as np
import sys
import argparse

# Load the command-line parameters
parser = argparse.ArgumentParser(description='Rescale periodic GAMER ELBDM initial conditions. \n\
                                 Example usage: ./elbdm_rescale_periodic_ic.py -n_in 256 -n_out 64 -input UM_IC_high_resolution -output UM_IC_low_resolution')

# Define command-line arguments
parser.add_argument('-float8', action='store', required=False, type=bool, dest='float8',
                    help='assume double precision for input and output file', default=False)
parser.add_argument('-n_in', action='store', required=True, type=int, dest='n_in',
                    help='input file resolution')
parser.add_argument('-n_out', action='store', required=True, type=int, dest='n_out',
                    help='output file resolution')
parser.add_argument('-input', action='store', required=False, type=str, dest='input',
                    help='input file [%(default)s]', default='./UM_IC_lr')
parser.add_argument('-output', action='store', required=False, type=str, dest='output',
                    help='output file [%(default)s]', default='./UM_IC')

# Parse the command-line arguments
args = parser.parse_args()

# Print the command-line arguments for reference
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
for t in range(len(sys.argv)):
    print(str(sys.argv[t]), end=' ')
print('')
print('-------------------------------------------------------------------\n')

# Set up interpolation parameters
input_file  = args.input
output_file = args.output
n_in        = args.n_in
n_out       = args.n_out
float8      = args.float8  # Enable double precision

if float8:
    rprec = np.double
    cprec = np.cdouble
else:
    rprec = np.float32
    cprec = np.complex64

# Credit to stackoverflow user Bily for this unpadding function
def unpad(x, pad_width):
    slices = []
    for c in pad_width:
        e = None if c[1] == 0 else -c[1]
        slices.append(slice(c[0], e))
    return x[tuple(slices)]

def interp(psi, n_in, n_out):
    if n_in == n_out:
        print("n_in == n_out, no rescaling necessary!")
        return psi

    n_total = np.prod(psi.shape)

    # Compute fft of array
    print("Computing FFT... ", end="")
    # Divide by n_total for forward normalisation
    psihat = np.fft.fftn(psi) / n_total
    print("done!")

    # Free input arrays to save memory
    del psi

    # Shift zero frequencies to center of cube
    psihat = np.fft.fftshift(psihat)

    # Upscale
    if n_out > n_in:
        # Pad cube from the outside
        n_pad = int(np.floor(n_out / 2 - n_in / 2))
        psihat = np.pad(psihat, ((n_pad, n_pad), (n_pad, n_pad), (n_pad, n_pad)), mode="constant")
    # Downscale
    elif n_in > n_out:
        # Unpad cube from the outside
        n_pad = int(np.floor(n_in / 2 - n_out / 2))
        psihat = unpad(psihat, ((n_pad, n_pad), (n_pad, n_pad), (n_pad, n_pad)))

    # Shift zero frequencies back to outside of cube
    psihat = np.fft.fftshift(psihat)

    print("Computing IFFT... ", end="")

    # Inverse FFT
    psihr = np.fft.ifftn(psihat)

    # Multiply by new n_total to undo backward normalisation
    psihr *= np.prod(psihr.shape)

    print("done!")

    # Delete large array in frequency domain to save memory
    del psihat

    return psihr

# Load binary from file
print("Loading input data... ", end="")
lr = np.fromfile(input_file, dtype=rprec).reshape((2, n_in, n_in, n_in))
psi = lr[0, :, :, :] + 1j * lr[1, :, :, :]
psi = psi.astype(cprec)
del lr
print("done!")

# Interpolate data
psihr = interp(psi, n_in, n_out).astype(cprec)

# Write data to disk
print("Writing wave function to binary file... ", end="")

# Overwrite output file with real part
with open(output_file, "wb") as f:
    np.real(psihr).astype(rprec).tofile(f)
    f.close()

# Append imaginary part
with open(output_file, "ab") as f:
    np.imag(psihr).astype(rprec).tofile(f)
    f.close()

print("done!")
