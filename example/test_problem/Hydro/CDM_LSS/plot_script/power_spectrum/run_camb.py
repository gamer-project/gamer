# %%
import numpy as np
import camb
from camb import model, initialpower

# Set parameters
class param:
  z = 0
  kmax = 1e2
  kmin = 1e-1
  npoints = 1000
  cosmology = "Planck18"
  nonlinear=True

# Set up cosmologal parameters
class Cosmo:
  def __init__(self, cosmology):
    if cosmology == "WMAP7":
      self.H0 = 70.2
      self.ombh2 = 0.02255
      self.omch2 = 0.1126
      self.As = 2.430e-9
      self.ns = 0.968
    elif cosmology == "WMAP9":
      self.H0 = 69.33
      self.ombh2 = 0.02266
      self.omch2 = 0.1157
      self.As = 2.427e-9
      self.ns = 0.971
    elif cosmology == "Planck15":
      self.H0 = 67.74
      self.ombh2 = 0.02230
      self.omch2 = 0.1188
      self.As = 2.142e-9
      self.ns = 0.9667
    elif cosmology == "Planck18":
      self.H0 = 67.66
      self.ombh2 = 0.02242
      self.omch2 = 0.11933  
      self.As = 2.105e-9    # sigma8 = 0.8102
      self.ns = 0.9665
    else:
      pass

cosmo = Cosmo(param.cosmology)

# Set up parameters for CAMB
pars = camb.CAMBparams()
pars.set_cosmology(H0=cosmo.H0, ombh2=cosmo.ombh2, omch2=cosmo.omch2)
pars.InitPower.set_params(As=cosmo.As, ns=cosmo.ns)
pars.set_matter_power(redshifts=[param.z], nonlinear=param.nonlinear, kmax=param.kmax)

# Linear spectra
results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=param.kmin, maxkh=param.kmax, npoints = param.npoints)

# save to file
pkm = pk / (8 * np.pi**3)
powerspec_music = np.vstack((kh, pkm))
header  = "k-Pk table assuming {} cosmology\n".format(param.cosmology)
header += "Power spectrum is computed with CAMB and divided by 8 pi^3 for consistency with MUSIC.\n"
header += "k [h/Mpc]  P(k) [Mpc/h]^3"
np.savetxt("powerspec_z{}_{}.dat".format(param.z, param.cosmology), powerspec_music.T, header=header)

