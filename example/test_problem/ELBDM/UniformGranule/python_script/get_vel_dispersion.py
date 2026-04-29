import yt
import numpy as np
import argparse
import sys

#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get average velocity dispersion of the entire box' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-m', action='store', required=False, type=str, dest='method',
                     help='gradient method', default='yt' )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
method      = args.method

if not (method == 'yt' or method == 'numpy' or method == 'fft'):
   print('Unknown method!!')
   sys.exit(1)

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )

yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   idx     = ds.parameters["DumpID"]
   time    = ds.parameters["Time"][0]
   N       = ds.parameters["NX0"]
   N_tot   = N[0]*N[1]*N[2]
   BoxSize = ds.parameters["BoxSize"]
   UNIT_L  = ds.parameters["Unit_L"]
   ma      = ds.parameters["ELBDM_Mass"]
   hbar    = ds.parameters["ELBDM_PlanckConst"]
   fac     = hbar/ma
   dh      = BoxSize/N

   if method == 'yt':
      grad_Real = ds.add_gradient_fields(("Real"))
      grad_Imag = ds.add_gradient_fields(("Imag"))

   dd = ds.all_data()

   if ( N_tot != len(dd["Dens"])):
      print('Data_%06d file size not matched!'%idx)
      sys.exit(1)

   if method == 'yt':
      dens = np.array(dd["Dens"])
      real = np.array(dd["Real"])
      imag = np.array(dd["Imag"])
      avedens = np.mean(dens)

   elif method == 'numpy' or method == 'fft':
      # re-arange arrays base on 3d position
      X = np.array( dd['x'] )
      Y = np.array( dd['y'] )
      Z = np.array( dd['z'] )

      if len(np.unique(X)) != N[0] or len(np.unique(Y)) != N[1] or len(np.unique(Z)) != N[2]:
         print('Array size not matched!'%idx)
         sys.exit(1)

      indices = np.lexsort((Z, Y, X))

      # re-arange then reshape to 3d array
      dens = np.reshape( np.array(dd["Dens"][indices]), (N[0], N[1], N[2]) )
      real = np.reshape( np.array(dd["Real"][indices]), (N[0], N[1], N[2]) )
      imag = np.reshape( np.array(dd["Imag"][indices]), (N[0], N[1], N[2]) )

      if method == 'fft':
         kx = 2*np.pi*np.fft.fftfreq(N[0], d=dh[0])
         ky = 2*np.pi*np.fft.fftfreq(N[1], d=dh[1])
         kz = 2*np.pi*np.fft.fftfreq(N[2], d=dh[2])

         KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
         psi_k = np.fft.fftn(real + 1j*imag)

         grad_psi = [ np.fft.ifftn(1j*KX*psi_k),
                      np.fft.ifftn(1j*KY*psi_k),
                      np.fft.ifftn(1j*KZ*psi_k) ]

   sigma_square_bk = [0., 0., 0.]
   sigma_square_qp = [0., 0., 0.]

   for i in range(3):
      if method == 'yt':
         field_r = "Real_gradient_%s"%(chr(ord('x') + i))
         field_i = "Imag_gradient_%s"%(chr(ord('x') + i))
         grad_real = np.array(dd[field_r])*UNIT_L
         grad_imag = np.array(dd[field_i])*UNIT_L

      elif method == 'numpy':
         grad_real = (np.roll(real, -1, axis=i) - np.roll(real, 1, axis=i)) / (2*dh[i])
         grad_imag = (np.roll(imag, -1, axis=i) - np.roll(imag, 1, axis=i)) / (2*dh[i])

      elif method == 'fft':
         grad_real = grad_psi[i].real
         grad_imag = grad_psi[i].imag

      v_bk = fac*(imag*grad_real - real*grad_imag)/dens
      v_qp = fac*(imag*grad_imag + real*grad_real)/dens

      sigma_square_bk[i] = np.average(v_bk**2, weights=dens) - np.average(v_bk, weights=dens)**2
      sigma_square_qp[i] = np.average(v_qp**2, weights=dens) - np.average(v_qp, weights=dens)**2

   sigma_bk = ((sigma_square_bk[0] + sigma_square_bk[1] + sigma_square_bk[2])/3.)**0.5
   sigma_qp = ((sigma_square_qp[0] + sigma_square_qp[1] + sigma_square_qp[2])/3.)**0.5
   sigma_total = (sigma_bk**2+sigma_qp**2)**0.5

   d = 0.35*2*np.pi*hbar/(ma*sigma_total)

   print('\nDumpID = %06d, time = %13.7e\n'%(idx, time) +
         'average density = %13.7e\n'%(np.mean(dens)) +
         'minimum density = %13.7e\n'%(np.min(dens)) +
         'method          = %s\n'%(method) +
         'velocity dispersion (bulk, thermal, total) = (%13.7e, %13.7e, %13.7e)\n'%(sigma_bk, sigma_qp, sigma_total) +
         'estimated granule diameter = %13.7e\n'%d)

