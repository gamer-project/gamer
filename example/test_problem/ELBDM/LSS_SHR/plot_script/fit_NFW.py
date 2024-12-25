# Fit NFW profile
# calculate potential energy in physical coordinate (without constant G)

from scipy.integrate import quad
from scipy.optimize import curve_fit
import numpy as np

def NFW_dens(r, dens_para):
    rho0 = dens_para[0]            # Msun/kpc^3
    Rs = dens_para[1]              # kpc
    return rho0/(r/Rs*(1+r/Rs)**2) # Msun/kpc^3


def fit_NFW_dens(r, dens, range):
    # fit in log scale
    def NFW_dens_log(log_r, rho0, Rs):
        r = 10**log_r
        return np.log10(NFW_dens(r, (rho0, Rs)))
    
    [rho0,Rs], pcov = curve_fit(NFW_dens_log, np.log10(r[range]),
                                np.log10(dens[range]),bounds=(0, [np.inf,200]))
    return rho0,Rs


def interp_dens(x, dens_para):
    r = dens_para[0]
    dens = dens_para[1]
    return 10**np.interp(np.log10(x), np.log10(r), np.log10(dens))

def shell_mass(r,dens):
    return 4*np.pi*r**2*dens

def potential_NFW(r,dens_para,halo_radius,halo_mass):
    rho0 = dens_para[0]
    Rs = dens_para[1]
    # ignores contributions outside the halo radius
    # set potential = GM/R at r = halo radius
    shift = 4*np.pi*rho0*Rs**3/halo_radius*np.log(1+ halo_radius/Rs) - halo_mass/halo_radius
    return 4*np.pi*rho0*Rs**3/r*np.log(1+ r/Rs) - shift

def potential_r(r,dens_para,halo_radius):
    # ignores contributions outside the halo radius
    def p(s):
        if s<r:
            return shell_mass(s,interp_dens(s,dens_para))/r
        else:
            return shell_mass(s,interp_dens(s,dens_para))/s

    potential_r, abserr = quad(p, 0, halo_radius, epsrel = 0.001)
    return potential_r

def shell_Ep(r, dens_para,halo_radius, halo_mass, typ):
    # Ep  = 1/2 \int_0^{halo_radius} potential(r) dM(r)
    if typ == 'NFW':
        return 0.5*potential_NFW(r, dens_para, halo_radius, halo_mass)*shell_mass(r, NFW_dens(r, dens_para))
    else:
        return 0.5*potential_r(r, dens_para,halo_radius)*shell_mass(r,interp_dens(r,dens_para))

def get_Ep_sph_sym(dens_para, halo_radius, halo_mass, typ='general'):
    # add up Ep of each shell
    ep, abserr = quad(lambda x: shell_Ep(x, dens_para, halo_radius, halo_mass, typ), 0, halo_radius, epsrel = 0.001)
    return ep/halo_mass**2*halo_radius

def NFW_mass(para, halo_radius):
    rho0 = para[0]
    Rs = para[1]
    return 4*np.pi*rho0*Rs**3*(np.log((Rs+halo_radius)/Rs)-halo_radius/(Rs+halo_radius))

def get_Ep_ideal(Mh, Rh, c, type):

    def f_c(c):
        # ignores contributions outside the halo radius
        return (2*c*(1+c)*np.log(1+c)-2*c**2-c**3)/((1+c)*np.log(1+c)-c)**2
    
    if type == 'NFW':
        Ep = -1/2*f_c(c)
    elif type == 'Top Hat':
        Ep = -0.6*Mh**2/Rh
    else:
        # not support yet
        print('Profile type not support yet')

    return Ep

if __name__== '__main__':

    rho0, Rs, halo_radius = 1e5, 10, 100
    radius = np.logspace(-1, 2.4, 128)
    density = NFW_dens(radius, (rho0, Rs))

    mass = NFW_mass((rho0,Rs), halo_radius)
    mass_acc, err = quad(lambda x: shell_mass(x, interp_dens(x,(radius,density))), 0, halo_radius, epsrel = 0.001)
    # check the mass is the same between formula and integral
    print("halo mass (from formula) = ", mass)
    print("halo mass (from integral)= ", mass_acc)

    # get the Ep from profile
    ep1 = get_Ep_sph_sym((radius,density), halo_radius, mass)
    print("Ep (from profile)        = ", ep1)

    # get Ep from formula
    ep2 = get_Ep_sph_sym((rho0, Rs, 10), halo_radius, mass, 'NFW')
    print("Ep (from NFW pontential) = ", ep2)
    ep3 = get_Ep_ideal(mass, halo_radius, halo_radius/Rs, 'NFW')
    print("Ep (from NFW formula)    = ", ep3)


