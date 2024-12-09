from mpi4py import MPI

import matplotlib.pyplot as plt
import numpy as np
import scipy
import os

from mpmath import *


# This script computes the tables required for the spectral interpolation in GAMER
#
# No OpenMP parallelisation, launch with as many MPI nodes as desired, e.g. with 3 nodes and
# mpirun -map-by ppr:16:socket:pe=1 python3 compute_interpolation_tables.py
# The script takes a few hours to run and will not output any log messages during that time.
#
# Spectral interpolation for the psidm branch can be enabled
# via the compile time option SUPPORT_SPECTRAL_INT

# The spectral interpolation algorithm requires two kinds of tables
# 1. Interpolation matrices that map an input vector of size N to an interpolated vector at 2 * (N - 2) with high accuracy
#    Given a vector with values f(x_i) with x_0 = 0, x_1 = 1, x_2 = 2 and N = 3
#    Multiplication with the interpolation matrices yields an output vector with the interpolated values
#    f(y_j) with y_0 = 0.75 and y_1 = 1.25
#
#    The respective tables are stored as binary files in double precision in the folder interpolation_tables
#
# 2. Tables that given an input function, compute its periodic Gram-Fourier extension
#    This is because for large N, the interpolation algorithm uses a Gram FE extension at runtime and interpolates the input function using the FFT algorithm
#
#    The respective tables are stored as binary files in double precision in the folder boundary2extension_tables

mp.dps = 256
eps    = 1e-256
comm   = MPI.COMM_WORLD
rank   = comm.Get_rank()
nprocs = comm.Get_size()

base_path = "./"

sub_directories = ["boundary2extension_tables", "interpolation_tables"]


if rank == 0:
    try:
        os.mkdir(base_path)
    except OSError as error:
        print(error)

    for sub_dir in sub_directories:
        try:
            os.mkdir(base_path + "/" + sub_dir)
        except OSError as error:
            print(error)

# Make other ranks wait for rank 0
comm.barrier()

# Implement the Gram Schmidt orthogonalisation algorithm using mpmath
class GramSchmidt:
    def __init__(self, x, m):
        self.x = x
        self.m = m
        self.A = mp.zeros(m, len(x))
        #Linear map for polynomial scalar product
        for i in range(m):
            for j in range(len(x)):
                #Polynomial basis {1, x, x^2, x^3, x^4, ..., x^m}
                self.A[i, j] = x[j]**i

        #Write basis vector as columns of matrix V
        self.V = mp.eye(m)

        self.U = self.modifiedGramSchmidt(self.V)

    def evaluateBasis(self, x, basis_element):
        #Linear map for polynomial scalar product
        A = mp.zeros(self.m, len(x))
        for i in range(self.m):
            for j in range(len(x)):
                #Polynomial basis {1, x, x^2, x^3, x^4, ..., x^m}
                A[i, j] = x[j]**i
        ei = self.U[:, basis_element].T * A

        return ei

    def sp(self, u, v):
        return mp.fsum((u.T * self.A) * (v.T * self.A).T)

    def proj(self, u, v):
        a1 = self.sp(v, u)
        a2 = self.sp(u, u)
        return a1/a2 * u

    def norm(self, u):
        return mp.sqrt(self.sp(u, u))

    def modifiedGramSchmidt(self, V):
        n, k = V.rows, V.cols
        U    = V.copy()
        U[:, 0] = V[:, 0] / self.norm(V[:, 0])

        for i in range(1, k):
            for j in range(i, k):
                U[:, j] = U[:, j] - self.proj(U[:, i - 1], U[:, j])


            U[:, i] = U[:, i] / self.norm(U[:, i])
        return U

    def projectFunction(self, f):
        coeffs = mp.matrix(1, self.m)

        for i in range(self.m):
            basis = (self.U[:, i].T * self.A)
            coeffs[0, i] = mp.fsum(f * basis.T)


        return coeffs

    def reconstructFunction(self, coeffs, x = None):
        if x == None:
            A = self.A
        else:
            A = mp.zeros(self.m, len(x))
            for i in range(self.m):
                for j in range(len(x)):
                    #Polynomial basis {1, x, x^2, x^3, x^4, ..., x^m}
                    self.A[i] = x[j]**i

        frec = mp.matrix(1, A.cols)
        for i in range(self.m):
            frec += coeffs[0, i] * (self.U[:, i].T * A)
        return frec

    def debug(self):
        m = self.m
        u_ij = mp.zeros(m)

        plt.title(f"Unorthogonalised polynomials m = {m}")
        for i in range(m):
            plt.plot(self.x, self.V[:, i].T * self.A, label=f"x^{i}")
        plt.legend()
        plt.show()



        plt.title(f"Orthonormalised polynomials m = {m}")
        for i in range(m):
            plt.plot(self.x, self.U[:, i].T * self.A, label=f"{i}")
        plt.legend()
        plt.show()

        print("The orthonormalised polynomials and their scalar products")
        for i in range(m):
            for j in range(m):
                u_ij[i, j] = self.sp(self.U[:, i], self.U[:, j])
            print(f"i = {i} u_ij = {u_ij[i, :]}")

class SVDFourierExtension:
    # Implement the SVD Fourier continuation algorithm using mpmath
    M_ALL_K  = 0
    M_EVEN_K = 1
    M_ODD_K  = 2

    def __init__(self, m, nDelta, nd, Gamma, g):
        self.m      = m
        self.nDelta = nDelta
        self.nd     = nd
        self.Gamma  = Gamma
        self.g      = g
        # Set up evaluation grid
        self.h      = 1/(nd - 1)
        self.d      = (nd - 1) * self.h
        self.Delta  = (nDelta  - 1) * self.h

        x = mp.linspace(0, 1, nd)

        # Compute left and right Gram Schmidt extensions
        leftBoundary  = x[       :nDelta]
        rightBoundary = x[-nDelta:      ]

        self.lgs = GramSchmidt(leftBoundary, m)
        self.rgs = GramSchmidt(rightBoundary, m)

        dxeval = self.Delta/(Gamma - 1)
        self.xeval  = mp.matrix(1, Gamma)
        for i in range(Gamma):
            self.xeval[0, i] = 1 - self.Delta + i * dxeval

        #Set up extension grid
        self.xext  = mp.linspace(1 - self.Delta, 1 + self.Delta + 2*self.d, 1000)
        mode  = self.M_EVEN_K
        M     = self.getM(g, Gamma, self.Delta, self.d, mode)
        Minv  = self.invertComplexM(M, 0)
        self.evencoeffs = []
        self.evenbasis  = []
        self.evenfrecs  = []
        for i in range(m):
            yeval = self.rgs.evaluateBasis(self.xeval, i)
            a     = self.iterativeRefinement(M, Minv, yeval)
            frec  = self.reconstruct(self.xext, a, g, Gamma, self.Delta, self.d, mode)
            self.evencoeffs.append(a)
            self.evenbasis.append(yeval)
            self.evenfrecs.append(frec)


        mode  = self.M_ODD_K
        M     = self.getM(g, Gamma, self.Delta, self.d, mode)
        Minv  = self.invertComplexM(M, 0)
        self.oddcoeffs = []
        self.oddbasis = []
        self.oddfrecs = []
        for i in range(m):
            yeval = self.rgs.evaluateBasis(self.xeval, i)
            a     = self.iterativeRefinement(M, Minv, yeval)
            frec  = self.reconstruct(self.xext, a, g, Gamma, self.Delta, self.d, mode)
            self.oddcoeffs.append(a)
            self.oddbasis.append(yeval)
            self.oddfrecs.append(frec)

        Next = 2 * nd + 2 * nDelta - 4
        xstore = mp.matrix(1, Next)
        for i in range(Next):
            xstore[i] = 1 - self.Delta + i * self.h

        self.F = mp.matrix(2 * m, Next)

        mode = self.M_EVEN_K

        for i in range(m):
            self.F[i, :] = self.reconstruct(xstore, self.evencoeffs[i], g, Gamma, self.Delta, self.d, mode)

        mode = self.M_ODD_K
        for i in range(m):
            self.F[i+m, :] = self.reconstruct(xstore, self.oddcoeffs[i], g, Gamma, self.Delta, self.d, mode)

        self.Pr = mp.matrix(m, nDelta)
        self.Pl = mp.matrix(m, nDelta)
        for i in range(m):
            self.Pr[i, :] = self.rgs.evaluateBasis(rightBoundary, i)
            self.Pl[i, :] = self.lgs.evaluateBasis(leftBoundary, i)


        #self.numpyF = np.array(self.F.apply(mp.re), dtype=np.float128).reshape(2 * m, Next)
        #self.numpyF.tofile(f"{base_path}/extension_tables/F_nD={nDelta}_nd={nd}_m={m}_g={g}_Gamma={Gamma}.bin")
        #self.numpyPr = np.array(self.Pr, dtype=np.float128).reshape(m, nDelta)
        #self.numpyPr.tofile(f"{base_path}/polynomial_tables/Pright_m={m}_nD={nDelta}.bin")
        #self.numpyPl = np.array(self.Pl, dtype=np.float128).reshape(m, nDelta)
        #self.numpyPl.tofile(f"{base_path}/polynomial_tables/Pleft_m={m}_nD={nDelta}.bin")

    # Pick Fourier modes
    def t(self, g, mode = M_ALL_K):
        if g % 2 == 0:
            k = np.arange(-int(-g/2) + 1, int(g/2) + 1)
        else:
            k = np.arange(-int((g-1)/2), int((g-1)/2) + 1)

        if mode == self.M_EVEN_K:
            k = k[k % 2 == 0]
        elif mode == self.M_ODD_K:
            k = k[k % 2 == 1]

        return k * mp.mpf(1)

    # Return array for evaluation of Gram polynomials
    def getX(self, Delta, Gamma):
        dxeval = Delta/(Gamma - 1)
        xeval  = mp.matrix(1, Gamma)
        for i in range(Gamma):
            xeval[0, i] = 1 - Delta + i * dxeval
        return xeval

    # Return array with values of plane waves at evaluation points
    def getM(self, g, Gamma, Delta, d, mode):
        ks = self.t(g, mode)
        x  = self.getX(Delta, Gamma)
        M  = mp.matrix(Gamma, len(ks))
        for i in range(Gamma):
            for j, k in enumerate(ks):
                M[i, j] = mp.exp(1j * k * np.pi / (d + Delta) * x[0, i])
        return M

    # Invert plane wave array using SVD with truncation of singular values below threshold
    def invertComplexM(self, M, cutoff):
        U, s, Vh = mp.svd(M)
        sinv = mp.diag(s)
        r = M.cols
        if M.rows < M.cols:
            r = M.rows
        for i in range(r):
            if s[i] < cutoff:
                sinv[i, i] = 0
            else:
                sinv[i, i] = 1/s[i]

        Vht = Vh.transpose_conj()
        Ut  = U.transpose_conj()
        f1  = sinv * Ut
        f2  = Vht * f1
        return  f2

    # Evaluate Fourier extension at point x
    def reconstruct(self, x, a, g, Gamma, Delta, d, mode):
        ks = self.t(g, mode)
        rec = mp.matrix(1, len(x))
        for j, coeff in enumerate(a):
            for i in range(len(x)):
                rec[i] += coeff * mp.exp(1j * ks[j] * np.pi / (d + Delta) * x[i])
        return rec

    # Iterative refinement for SVD-based matrix inversion
    def iterativeRefinement(self, M, Minv, f, threshold = 100, maxiter = 3000):
        a       = Minv * f.T
        r       = M * a - f.T
        counter = 0
        while mp.norm(r) > eps * mp.norm(a) and counter < maxiter:
            delta    = Minv * r
            a        = a - delta
            r        = M * a - f.T
            counter += 1
        return a

    def computeExtension(self, x, g, Gamma, Delta, d, mode, f):
        M     = self.getM(g, Gamma, Delta, d, mode)
        Minv  = self.invertComplexM(M, 0)
        a     = self.iterativeRefinement(M, Minv, f)
        frec  = self.reconstruct(x, a, g, Gamma, Delta, d, mode)
        return frec

    def debug(self):
        fig, axs = plt.subplots(self.m, 2, figsize=(3.54*2, 3.54*4), dpi=100)
        fig.suptitle("Reproduce figure 5.2 in Mark Lyon's thesis")

        for i, (ybasis, yrec) in enumerate(zip(self.evenbasis, self.evenfrecs)):
            axs[i,0].plot(self.xeval, ybasis, lw = 5)
            y = np.array([mp.re(yrec[i]) for i in range(len(yrec))])
            axs[i,0].plot(self.xext, y)
        for i, (ybasis, yrec) in enumerate(zip(self.oddbasis, self.oddfrecs)):
            axs[i,1].plot(self.xeval, ybasis, lw = 5)
            y = np.array([mp.re(yrec[i]) for i in range(len(yrec))])
            axs[i,1].plot(self.xext, y)
        plt.show()

        fig, axs = plt.subplots(self.m, 2, figsize=(3.54*2, 3.54*4), dpi=200)
        fig.suptitle("Verify that imaginary parts are zero up to double precision")
        for i, (ybasis, yrec) in enumerate(zip(self.evenbasis, self.evenfrecs)):
            y = [mp.im(yrec[i]) for i in range(len(yrec))]
            axs[i,0].plot(self.xext, y)
        for i, (ybasis, yrec) in enumerate(zip(self.oddbasis, self.oddfrecs)):
            y = [mp.im(yrec[i]) for i in range(len(yrec))]
            axs[i,1].plot(self.xext, y)
        plt.show()

        print("Verify how well the Gram polynomials are approximated in the new basis...\n")

        xtest  = mp.linspace(1 - self.Delta, 1, 900)
        r = self.m
        mode = self.M_EVEN_K
        evenerrors = []
        for i in range(r):
            yeval = self.rgs.evaluateBasis(xtest, i)
            frec  = self.reconstruct(xtest, self.evencoeffs[i], self.g, self.Gamma, self.Delta, self.d, mode)
            y = np.array([mp.re(frec[i]) for i in range(len(frec))])
            evenerrors.append(mp.norm(yeval - frec, p = mp.inf))
            plt.title(f"Even error {np.float(evenerrors[-1]):1.1e}")
            plt.plot(xtest, yeval, label=f"{i} even")
            plt.plot(xtest, y, label=f"{i} even")
            plt.legend()
            plt.show()

        mode = self.M_ODD_K
        odderrors = []
        for i in range(r):
            yeval = self.rgs.evaluateBasis(xtest, i)
            frec  = self.reconstruct(xtest, self.oddcoeffs[i], self.g, self.Gamma, self.Delta, self.d, mode)
            y = np.array([mp.re(frec[i]) for i in range(len(frec))])
            odderrors.append(mp.norm(yeval - frec, p = mp.inf))
            plt.title(f"Odd error {np.float(odderrors[-1]):1.1e}")
            plt.plot(xtest, yeval, label=f"{i} even")
            plt.plot(xtest, y, label=f"{i} even")
            plt.legend()
            plt.show()

        for i in range(r):
            print(f"f{i}_even: {evenerrors[i]} f{i}_odd: {odderrors[i]}")

class GramFEFixedSizeExtension:

    def __init__(self, N, m, nDelta, nd, Gamma, g):
        self.N      = N

        # compute accurate Gram-Fourier extension
        extension = SVDFourierExtension(m, nDelta, nd, Gamma, g)
        self.extension = extension

        # size of extension
        nExt = nd - 2

        # matrix with left and right Gram polynomials
        # Pb  = np.block([[Pl, np.zeros((nDelta, nDelta))], [np.zeros((nDelta, nDelta)), Pr]])
        Pb = mp.zeros(2*m, 2*nDelta)
        for i in range(m):
            for j in range(nDelta):
                Pb[i,     j         ] = extension.Pl[i, j]
                Pb[i + m, j + nDelta] = extension.Pr[i, j]

        # matrix that combines left and right Gram polynomials to get even and odd extensions
        # mix = np.block([[np.identity(nDelta) * 0.5, np.identity(nDelta) * 0.5], [np.identity(nDelta) * (-0.5), np.identity(nDelta) * 0.5]])
        Shuffle = mp.zeros(2*m)
        for i in range(m):
            Shuffle[i    , i    ] = +mp.mpf(1)/mp.mpf(2)
            Shuffle[i    , i + m] = +mp.mpf(1)/mp.mpf(2)
            Shuffle[i + m, i    ] = -mp.mpf(1)/mp.mpf(2)
            Shuffle[i + m, i + m] = +mp.mpf(1)/mp.mpf(2)

        # matrix that combines left and right extensions
        # Fb = np.transpose(np.concatenate([Fe, Fo], axis=0))
        # Fe   = F[:nDelta, nDelta:nDelta + nd - 2]
        # Fo   = F[nDelta:, nDelta:nDelta + nd - 2]

        Fb = mp.matrix(nExt, 2 * m)
        for i in range(nExt):
            for j in range(m):
                Fb[i, j    ] = extension.F[j    , nDelta + i]
                Fb[i, j + m] = extension.F[j + m, nDelta + i]

        # matrix that maps the input function to the extended domain
        #extendWavefunction = np.block([[np.identity(N)],
        #                              [np.identity(nDelta), np.zeros((nDelta, N - nDelta))],
        #                              [np.zeros((nDelta, N - nDelta)), np.identity(nDelta)]]
        #                              )
        extendWavefunction = mp.zeros(N + 2 * nDelta, N)
        for i in range(N):
            extendWavefunction[i,i] = 1

        for i in range(nDelta):
            extendWavefunction[N          + i,              i] = 1
            extendWavefunction[N + nDelta + i, N - nDelta + i] = 1

        # matrix that maps left and right boundary to extension
        self.boundary2Extension      = Fb * Shuffle * Pb
        self.numpyboundary2Extension = np.array(self.boundary2Extension.tolist(), dtype=np.complex256).reshape(nExt, 2*nDelta).astype(np.float128)
        self.numpyboundary2Extension.tofile(f"{base_path}/boundary2extension_tables/nD={nDelta}_nd={nd}_m={m}_g={g}_Gamma={Gamma}.bin")

        self.nExtended = N + nExt

        # matrix that maps the input function to the extended wave function
        #computeExtension = np.block([[np.identity(N), np.zeros((N, 2*nDelta))],
        #                            [np.zeros((nd - 2, N)), transform]])
        computeExtension = mp.zeros( N + nExt, N + 2*nDelta )
        for i in range(N):
            computeExtension[i, i] = 1
        for i in range(nExt):
            for j in range(2 * nDelta):
                computeExtension[N + i, N +j] = self.boundary2Extension[i, j]

        self.computeExtension = computeExtension @ extendWavefunction
        #self.numpyComputeExtension = np.array(self.computeExtension.tolist(), dtype=np.complex256).reshape(self.nExtended, self.N).astype(np.float128)
        #self.numpyComputeExtension.tofile(f"{base_path}/extension_tables/N={self.N}_nD={nDelta}_nd={nd}_m={m}_g={g}_Gamma={Gamma}.bin")

        # matrix that computes FFT
        self.computeFFT  = self.dftmat(self.nExtended)
        #self.numpyFFTMatrix = np.array(self.computeFFT.tolist(), dtype=np.complex256).reshape(self.nExtended, self.nExtended)
        #self.numpyFFTMatrix.tofile(f"{base_path}/fft_tables/N={self.nExtended}.bin")

        # matrix that extends and computes FFT
        self.computeExtendedFFT = self.computeFFT * self.computeExtension
        #self.numpyComputeExtensionFFT = np.array(self.computeExtendedFFT.tolist(), dtype=np.complex256).reshape(self.nExtended, self.N)
        #self.numpyComputeExtensionFFT.tofile(f"{base_path}/extension_tables/FFT_N={self.N}_nD={nDelta}_nd={nd}_m={m}_g={g}_Gamma={Gamma}.bin")

        # matrix that computes inverse FFT
        self.computeIFFT  = self.idftmat(self.nExtended)
        # self.numpyIFFTMatrix = np.array(self.computeIFFT.tolist(), dtype=np.complex256).reshape(self.nExtended, self.nExtended)
        # self.numpyIFFTMatrix.tofile(f"{base_path}/ifft_tables/N={self.nExtended}.bin")

    def dftmat(self, N):
        M = mp.matrix(N, N)
        for i in range(N):
            for j in range(N):
                M[i, j] = mp.exp(-2j * mp.pi * i * j / N)

        return M

    def idftmat(self, N):
        M = mp.matrix(N, N)
        for i in range(N):
            for j in range(N):
                M[i, j] = mp.exp(+2j * mp.pi * i * j / N) / N

        return M

    def debug(self, func = lambda x: np.sin(10 * x)):
        f = func(np.linspace(0, 1, self.N))
        fext = self.numpyComputeExtension @ f
        plt.title("x-space")
        plt.plot(f, label="f")
        plt.plot(fext, label="periodic extension of f")
        plt.legend()
        plt.show()

        fhat1 = scipy.fft.fft(fext)
        fhat2 = self.numpyComputeExtensionFFT @ f

        plt.title("k-space")
        plt.loglog(np.abs(fhat1), label="Scipy FFT")
        plt.loglog(np.abs(fhat2), label="Matrix FFT")
        plt.legend()
        plt.show()

        print(f"Difference between f and inverse f: {np.mean(np.abs((self.numpyIFFTMatrix @ fhat2)[:len(f)] - f))}")



class GramFEInterpolation:

    def __init__(self, N, m, nDelta, nd, Gamma, g):
        self.N      = N
        self.m      = m
        self.nDelta = nDelta
        self.nd     = nd
        self.Gamma  = Gamma
        self.g      = g

        # compute accurate Gram-Fourier extension
        self.extension = GramFEFixedSizeExtension(N, m, nDelta, nd, Gamma, g)

        # matrix that evaluates extended input function in k-space at the interpolation points xx
        xend = 1 * ((N - 1) / (self.extension.nExtended))
        dx   = xend / (N - 1)
        xx   = mp.linspace(mp.mpf(3)/mp.mpf(4) * dx, xend - mp.mpf(3)/mp.mpf(4) * dx , 2 * (N - 2))
        computeInterpolation = self.computeInterpolationMatrix(self.extension.nExtended, xx)


        # matrix that maps the input function of size N to the 2 * (N - 2) interpolated values (requires a ghost boundary of at least one)
        self.interpolationMatrix = computeInterpolation @ self.extension.computeFFT @ self.extension.computeExtension
        self.numpyInterpolationMatrix = np.array(self.interpolationMatrix.tolist(), dtype=complex).reshape(2 * (N - 2), N).astype(float)
        self.numpyInterpolationMatrix.tofile(f"{base_path}/interpolation_tables/N={N}_m={m}_nDelta={nDelta}_nd={nd}_Gamma={Gamma}_g={g}.bin")



    def computeInterpolationMatrix(self, N, xarray):
        N1 = len(xarray)
        M = mp.matrix(N1, N)
        for i in range(N1):
            for j in range(N):
                if j < N/2 + 1:
                    kn = j
                else:
                    kn = j - N
                M[i, j] = mp.exp(2j * mp.pi * xarray[i] * kn) / N
        return M

    def debug(self, func = lambda x: np.sin(10 * x)):

        # Generate the x values at which to sample the function
        x = np.linspace(0, 1, self.N)
        dx = x[1]-x[0]

        # Evaluate the function at the sample points
        y = func(x)

        dx = 1/(self.N - 1)
        xx = np.linspace(0.75 * dx, 1 - 0.75*dx, 2 * (self.N - 2))

        # Compute the FFT of the function
        interp = self.numpyInterpolationMatrix @ y


        # Print the interpolated value
        plt.title("Gram interpolation routine")
        plt.plot(x,func(x),'-', label="Original")
        plt.plot(xx, interp, label="Interpolated")
        plt.legend()
        plt.show()

        plt.title("Interpolation error")
        plt.yscale("log")
        plt.plot(xx, np.abs(interp-func(xx)),'.')
        plt.show()


# Fixed parameters of Gram FE algorithn
Gamma  = 150
g      = 63

# Input sizes for interpolation tables
N_min  = 4
N_max  = 32
interpolation_table_iterations = N_max + 1 - N_min

# Extension sizes for extension tables
nd_min = 24
nd_max = 36
extension_table_iterations = nd_max + 1 - nd_min

total_iterations = interpolation_table_iterations + extension_table_iterations

# Compute the interpolation and extension tables in parallel loop
for i in range(rank, total_iterations + 1, nprocs):
    if i < interpolation_table_iterations:
        N  = i + N_min
        nd = 32
        m  = 8

        print(f"Rank {rank}: Computing interpolation table for N = {N}")
        # Let polynomial order be size of interpolation domain
        if N <= m:
            GramFEInterpolation(N = N, m = N, nDelta = N, nd = nd, Gamma = Gamma, g = g)
        # Keep polynomial order fixed for larger domain sizes for stability
        else:
            GramFEInterpolation(N = N, m = m, nDelta = m, nd = nd, Gamma = Gamma, g = g)
    else:
        # Compute extension tables with different sizes for fast FFTs
        nd = i - interpolation_table_iterations + nd_min
        m  = 8

        print(f"Rank {rank}: Computing fixed extension for nd = {nd}")
        # Keyword N not used for computing extension tables, 16 is placeholder
        GramFEFixedSizeExtension(N = 16, m = m, nDelta = m, nd = nd, Gamma = Gamma, g = g)