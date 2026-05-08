import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import yt

# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
sigma  = 6.0                                  # velocity dispersion in km/s
# -------------------------------------------------------------------------------------------------------------------------
ds     = yt.load( '../Data_000000' )
UNIT_L = ds.parameters["Unit_L"]              # code length unit (default is kpc)
UNIT_V = ds.parameters["Unit_V"]              # code velocity unit (default is km/s)
UNIT_T = ds.parameters["Unit_T"]
ma     = ds.parameters["ELBDM_Mass"]          # FDM particle mass in code unit
h_bar  = ds.parameters["ELBDM_PlanckConst"]   # reduced planck constant in code unit
d      = 0.35*2.0*np.pi*h_bar/ma/sigma        # granule diameter in code uniti (kpc)

# C(t) decay time scale, C(t_corr)=(1+2*0.35**2)**(-1.5)*C(t=0) = 0.72*C(t=0)
t_corr = d/(2**0.5*np.pi*sigma)*UNIT_T/(3600.*24.*365.*1e6) # expected correlation time scale in Myr

print("velocity dispersion             = %g km/s"%sigma)
print("estimated granule size          = %g kpc"%d)
print("expected correlation time scale = %g Myr"%t_corr)
print("")

file_dir = '../Record__Correlation/'
files = glob.glob(file_dir + 'correlation_function_t=*.txt')

filename = np.array(files)

r = []
C = []
t = []
for f in filename:
   Corr = np.loadtxt(f, skiprows=1, dtype=float)
   if not r:
      r.append(Corr[:,0])
   else:
      if not np.array_equal(r[0], Corr[:,0]):
         print('radius bin not matched!! filename=\"%s\"'%f)
         sys.exit(1)

   C.append(Corr[:,1])
   match = re.search(r'correlation_function_t=(\d+\.\d+e[+-]\d+)', f)
   if match:
      time = float(match.group(1))
      t.append(time)
   else:
      print('time pattern not matched!! filename=\"%s\"'%f)
      sys.exit(1)

if not r:
   print("no file loaded, check ../Record__Correlation/ !!")

t = np.array(t)*UNIT_T/(3600.*24.*365.*1e6)
C = np.array(C)

color = plt.cm.turbo(np.linspace(0.1, 0.9, len(r[0])))

plt.figure()
for i in range(len(r[0])):
   plt.plot(t, C[:,i], c=color[i], label = r'$r$ = %1.3e kpc'%(r[0][i]))
plt.axvline(x=t_corr, ls='--', lw = '0.9', color = '0.7', label = 'expected correlation time scale')
plt.xlabel('$t$ (Myr)')
plt.ylabel(r'$C(t)$')
plt.xlim(0, t[-1])
plt.legend(bbox_to_anchor=(1.03,0.03), loc='lower left')
plt.savefig('fig_correlation.png', dpi = 150, bbox_inches="tight")
plt.close()


