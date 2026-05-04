import numpy as np
import argparse
import sys
import warnings
from scipy.integrate import quad


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
Div_disk = 500        # total number of evenly divided radial bins
r_min    = 5.0e-2
r_max    = 1.2e+2
ratio    = (r_max/r_min)**(1.0/(Div_disk-1.0))
r_log    = np.zeros(Div_disk)

for i in range(Div_disk):
   r_log[i] = r_min * (ratio)**(i)


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Convert halo properties for later heating rate calculation' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


# -------------------------------------------------------------------------------------------------------------------------
# specify script unit setting
Unit_L = 3.08568e+21
Unit_M = 1.989e+43
Unit_V = 100000
G      = 6.6743e-8*Unit_M/Unit_L


#-----------------------------------------------------------------------------------------------------------------
# predefined functions
def interpolation(r, array_radius, array_value):
   count = 0
   log_radius = np.log(array_radius)
   while(r > array_radius[count]):
      count = count+1
   output = array_value[count-1]+(array_value[count]-array_value[count-1])*(np.log(r)-log_radius[count-1])/(log_radius[count]-log_radius[count-1])
   return output

def log_interpolation(r, array_radius, array_value):
   count = 0
   log_value = np.log(array_value)
   log_radius = np.log(array_radius)
   while(r > array_radius[count]):
      count = count+1
   output = log_value[count-1]+(log_value[count]-log_value[count-1])*(np.log(r)-log_radius[count-1])/(log_radius[count]-log_radius[count-1])
   return np.exp(output)
def density(r):
   if (isinstance(r, int) or isinstance(r, float)):
      if r < r_d[0]:
         return dens[0]
      elif r > r_d[-1]:
         return 0
      else:
         return log_interpolation(r, r_d, dens)
   else:
      out = np.zeros(len(r))
      for i in range(len(r)):
         if r[i] < r_d[0]:
            out[i] = dens[0]
         elif r[i] > r_d[-1]:
            out[i] = dens[-1]
         else:
            out[i] = log_interpolation(r[i], r_d, dens)
      return out

def enclosed_mass_int(r):
   warnings.filterwarnings("ignore")
   if (isinstance(r, int) or isinstance(r, float)):
      return quad(lambda x: 4*np.pi*x*x*density(x), 0, r)[0]
   else:
      temp = []
      for i in r:
         temp.append( quad(lambda x: 4*np.pi*x*x*density(x), 0, i)[0] )
      result = np.array(temp)
      return result

def potential(r):
   if (isinstance(r, int) or isinstance(r, float)):
      if r < r_p[0]:
         return pote[0]
      elif r > r_p[-1]:
         return pote[-1]
      else:
         return log_interpolation(r, r_p, pote)
   else:
      out = np.zeros(len(r))
      for i in range(len(r)):
         if r[i] < r_p[0]:
            out[i] = pote[0]
         elif r[i] > r_p[-1]:
            out[i] = pote[-1]
         else:
            out[i] = interpolation(r[i], r_p, pote)
      return out

def jeans_int(r):
   if (isinstance(r, int) or isinstance(r, float)):
      return quad(lambda x: density(x)*jeans_interpolation(x), r, np.inf)[0]
   else:
      temp = []
      for i in r:
         temp.append( quad(lambda x: density(x)*jeans_interpolation(x), i, np.inf)[0] )
      result = np.array(temp)
      return result

def jeans_interpolation(r):
   if (isinstance(r, int) or isinstance(r, float)):
      if r < r_p[0]:
         return DPotDr[0]
      elif r > r_p[-1]:
         return DPotDr[-1]
      else:
         return interpolation(r, r_p, DPotDr)
   else:
      out = np.zeros(len(r))
      for i in range(len(r)):
         if r[i] < r_p[0]:
            out[i] = DPotDr[0]
         elif r[i] > r_p[-1]:
            out[i] = DPotDr[-1]
         else:
            out[i] = interpolation(r[i], r_p, DPotDr)
      return out


#-----------------------------------------------------------------------------------------------------------------
# analyze simulation data and generate processed output files
for idx in range(idx_start, idx_end+1, didx):
   Dens =  np.load('Halo_Dens_Data_%06d.npy'%idx)
   r_d = Dens[0]
   dens_cgs = Dens[1]
   dens = dens_cgs*Unit_L**3/Unit_M
   enmass = enclosed_mass_int(r_log)
   mass_cgs = enmass*Unit_M

   Pote =  np.load('Halo_Pote_Data_%06d.npy'%idx)
   r_p = Pote[0]
   pote = Pote[1]
   pote_cgs = potential(r_log)

   DPotDr = np.gradient(pote, r_p)
   sigma = ((jeans_int(r_p)/density(r_p))*0.5)**0.5

   r_in_kpc = np.zeros(Div_disk)
   rho = np.zeros(Div_disk)
   M_h =  np.zeros(Div_disk)
   disp = np.zeros(Div_disk)
   phi_h = np.zeros(Div_disk)

   r = 0
   num_dens = 0
   num_log = 0
   num_disp = 0

   dr = 15.0/float(Div_disk)
   for j in range(Div_disk):
      r_in_kpc[j] = r + 0.5 * dr
      while r_d[num_dens] < (r + 0.5*dr):
         num_dens += 1
      rho[j] = dens_cgs[num_dens-1] + (dens_cgs[num_dens] - dens_cgs[num_dens-1])*(r_in_kpc[j]-r_d[num_dens-1])/(r_d[num_dens]-r_d[num_dens-1])
      while r_log[num_log] < (r + 0.5*dr):
         num_log += 1
      M_h[j] = mass_cgs[num_log-1] + (mass_cgs[num_log] - mass_cgs[num_log-1])*(r_in_kpc[j]-r_log[num_log-1])/(r_log[num_log]-r_log[num_log-1])
      while r_p[num_disp] < (r + 0.5*dr):
         num_disp += 1
      disp[j] = sigma[num_disp-1] + (sigma[num_disp] - sigma[num_disp-1])*(r_in_kpc[j]-r_p[num_disp-1])/(r_p[num_disp]-r_p[num_disp-1])

      r = r + dr

   Data = np.zeros((6,Div_disk))
   Data[0] = r_in_kpc
   Data[1] = rho       # halo density
   Data[2] = M_h       # enclosed mass
   Data[3] = disp      # Jeans velocity dispersion/sqrt(2)
   Data[4] = r_log     # r in log bins (used for <P2> computation)
   Data[5] = pote_cgs  # Phi(halo potential) in log bins (used for <P2> computation)

   np.save('Data_Halo_%06d'%idx, Data)
   print('Data_Halo_%06d.npy completed'%idx)
