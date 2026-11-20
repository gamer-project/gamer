import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt
import yt

# expectation
try:
    ds = yt.load('../Data_000000')
    expected_num_SNe_per_par = ds.parameters['FB_ResolvedSNeII_NPerMass']*ds.parameters['SF_CreateStar_MinStarMass']
except:
    SF_CREATE_STAR_MIN_STAR_MASS = 3.2e4
    FB_RESOLVED_SNEII_N_PER_MASS = 1.0912e-2
    expected_num_SNe_per_par = SF_CREATE_STAR_MIN_STAR_MASS * FB_RESOLVED_SNEII_N_PER_MASS
num_SNe_per_par = np.arange(0, 2*int(np.ceil(expected_num_SNe_per_par)))

# simulation results
feedback_record = np.genfromtxt('../Record__FB_Resolved_SNeII',
                                delimiter=None,
                                comments='#',
                                names=['Rank', 'TID', 'lv', 'TimeOld', 'TimeNew', 'SNII_Time_first', 'SNII_Time_last', 'SNII_AccumNum', 'SNII_Energy', 'SNII_Mass', 'SNII_Metal', 'FB_Diameter', 'FB_Flu_Mass', 'Par_PUID', 'Par_NxtE', 'Par_Mass', 'Par_PosX', 'Par_PosY', 'Par_PosZ', 'Par_VelX', 'Par_VelY', 'Par_VelZ', 'Par_MetalMass', 'Flu_Dens', 'Flu_MomX', 'Flu_MomY', 'Flu_MomZ', 'Flu_Engy', 'Flu_Metal'],
                                dtype=None,
                                encoding=None)
num_of_SNe = feedback_record['Par_NxtE'][:].astype('int') // 10000
idx_of_SNe = feedback_record['Par_NxtE'][:].astype('int') %  10000
num_of_completed_SNe = num_of_SNe[idx_of_SNe==(num_of_SNe-1)]

# histogram
hist = np.histogram(num_of_completed_SNe,
                    bins=[i-0.5 for i in range(0,len(num_SNe_per_par)+1)])[0]
prob = hist/np.sum(hist)
plt.plot(num_SNe_per_par, prob,
         'o-', color='C0', label='Simulation')
plt.axvline(np.mean(np.repeat(num_SNe_per_par, hist)),
            linestyle='--', color='C0')

# Poisson distribution
plt.plot(num_SNe_per_par, poisson.pmf(num_SNe_per_par, expected_num_SNe_per_par),
         '-', color='C1', label='Expected Poisson')
plt.axvline(expected_num_SNe_per_par,
            color='C1', linestyle=':')

# formatting
plt.legend()
plt.xlabel('Number of SN Events Each Particle')
plt.ylabel('Probability')
plt.yscale('log')
plt.ylim(bottom=0.5/len(num_of_completed_SNe), top=2.0*np.max(prob))
plt.xlim(np.min(num_of_completed_SNe) +
         0.5*(np.min(num_of_completed_SNe)-np.mean(num_of_completed_SNe)),
         np.max(num_of_completed_SNe) +
         0.5*(np.max(num_of_completed_SNe)-np.mean(num_of_completed_SNe)))
plt.savefig('fig__NumSNeParticleDistribution.png', bbox_inches='tight')
