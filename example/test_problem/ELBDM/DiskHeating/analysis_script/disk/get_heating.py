import numpy as np
import argparse
import sys
import scipy.special as sp
from scipy.integrate import quad
import yt


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
Div_disk = 500        # total number of evenly divided radial bins
dr = 15.0/Div_disk

r_min = 5.0e-2
r_max = 1.2e+2
n_bin = 500
ratio = (r_max/r_min)**(1.0/(n_bin-1.0))


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Compute theoretical heating rate of each snapshot' )

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
'''
unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}
'''
G = 6.6743e-8   # gravitational constant (cm**3/(g*s**2))
hb = 1.0546e-27 # reduced planck const cm**2*g/s
kappa1 = 0.526
kappa2 = 1.0
ILambda1 = 4.82
ILambda2 = 3.33


#-----------------------------------------------------------------------------------------------------------------
# predefined functions
def Gx(x):
   return (sp.erf(x)-2*x*np.exp(-x*x)/(np.pi)**0.5)/2/x**2

def interpolation(r, array_radius, array_value):
   log_radius = np.log(array_radius)
   count = int(np.log(r/r_min)/np.log(ratio))
   while(r > array_radius[count]):
      count = count+1
   output = array_value[count-1]+(array_value[count]-array_value[count-1])*(np.log(r)-log_radius[count-1])/(log_radius[count]-log_radius[count-1])
   return output

def phi_bg(R, x, h):
   z = x*h
   radius = (R*R+z*z)**0.5/3.08568e+21
   if radius < r_p[0]:
      out = pote_data[0]
   elif radius > r_p[-2]:
      out = pote_data[-2]
   else:
      out = interpolation(radius, r_p, pote_data)
   return out

def phi_int(r, h, x, x1):
   phi_diff = phi_bg(r, x, h) - phi_bg(r, x1, h)
   if((phi_diff) < 1.0e-50):
      return 0.0
   else:
      return 2.0/(2.0*phi_diff)**0.5

def IP2(r, x, h):
   return quad(lambda x1: phi_int(r, h, x, x1), 0, x)[0]

def IP2_log_IP2_int(r, x, h):
   ip2 = IP2(r, x, h)
   if (ip2 < 1.0e-50):
      return 0.0
   else:
      return ip2*np.log(ip2)*2.0*np.exp(-x*x)/(np.pi**0.5)


#-----------------------------------------------------------------------------------------------------------------
# analyze simulation data and generate processed output files
ds = yt.load( '../../Data_%06d'%idx_start )
ma = ds.parameters["ELBDM_Mass"]*ds.parameters["Unit_M"]      # FDM particle mass [gram]

for idx in range(idx_start, idx_end+1, didx):
   data       = np.load('Data_Disk_%06d.npy'%idx)
   data_h     = np.load('../halo/Data_Halo_%06d.npy'%idx)

   # data in cgs
   r          = data[0]
   r_in_kpc   = r/3.08568e+21
   v_c        = data[1]
   sigma_r    = data[2]
   sigma_z    = data[3]
   Sigma      = data[4]
   v_c2       = data[5]
   enmass     = data[6] + data_h[2]
   dMdR       = np.gradient(enmass, r)
   rho_eff    = dMdR/(4.0*np.pi*r**2)-enmass/(4.0*np.pi*r**3)
   sigma_h    = data_h[3]
   rho_h      = data_h[1]
   rho_bg_eff = rho_h - rho_eff
   rho_bg_eff = rho_bg_eff.clip(min=0)
   nu         = (((v_c2 + sigma_z**2 + sigma_r**2)/3.0)**0.5)/sigma_z
   r_p        = data_h[4]
   pote_data  = data_h[5]
   h          = data[7]
   mean_P1    = sigma_z*4.82/(np.pi*G*Sigma)
   mean_P2    = np.zeros(Div_disk)

   for i in range(Div_disk):
      mean_P2[i] =  2*h[i]*IP2(r[i], 1, h[i])
#   np.save('mean_P2_%06d'%idx, mean_P2)
#   mean_P2 = np.load('mean_P2_%06d.npy'%idx)

   X = nu*sigma_z*3**0.5/2**0.5/sigma_h
   E = (4.0*X)/(3.0*nu**2*np.pi**0.5)*np.exp(-X**2) + (1-nu**(-2))*(sp.erf(X)-Gx(X))
   tau = hb/(2**0.5*ma*sigma_h**2)
   M = (np.pi**(2.5)*G**2*hb**3*rho_h**2)/(6**0.5*ma**3*sigma_h**3*nu*sigma_z)*E

   heat_SGD = kappa1*M*np.log(mean_P1/tau)
   heat_BGD = kappa2*M*np.log(mean_P2/tau)

   heat = np.zeros(Div_disk)
   for i in range(Div_disk):
      k = sigma_z[i]*(8*G*rho_bg_eff[i])**0.5/(np.pi*G*Sigma[i])
      if (k<0.5):
          heat[i] = heat_SGD[i]
      elif (k>1.5):
          heat[i] = heat_BGD[i]
      else:
          heat[i] = (k-0.5)*heat_BGD[i] + (1.5-k)*heat_SGD[i]

   np.savez('Heating_%06d'%idx, a = r, b = heat)
