#######################################

# Computation of Power spectrum within a target halo.
# Run Compute_profile.py to call these functions.

#######################################
import os
import yt
import numpy as np
import csv

rfac  = 10 # the enclosed box has length of 2 * rfac * rc
rmfac = 2 # density within rmfac * rc is zeroed
lv    = 7

def corezero( x, y, z, rho, r_min, center ):

    r = np.sqrt( (x - center[0])**2 +
                 (y - center[1])**2 +
                 (z - center[2])**2)
    mask = np.where( (r < r_min) )[0]
    rho[mask] = 0.0
    return rho

def get_frb_cube( ds, c, radius, lv, dxs ):

    L       = radius*2
    frb_res = np.int64( L/dxs[lv] )
    cedge   = c - radius
    cube    = ds.covering_grid( lv, left_edge=cedge, dims=[frb_res, frb_res, frb_res] )

    return cube

def density_profile( x, y, z, rho, center ):

    r = np.sqrt( (x - center[0])**2 +
                 (y - center[1])**2 +
                 (z - center[2])**2)

    N            = np.int64( 8e1 )
    rbin_edge    = np.logspace( -3, np.log10(r.max()), N )
    rbin         = (rbin_edge[1:] + rbin_edge[:-1])/2
    densprof     = np.zeros( N - 1 )
    counts       = np.zeros( N - 1 )

    inds = np.digitize( r, rbin_edge )
    inds -= 1

    for ind in np.unique( inds ):
        if ( ind > 0 and ind < inds.max() ):
            densprof[ind] += np.sum( (rho)[inds==ind] )
            counts  [ind] += np.sum( inds==ind )

    mask = (counts==0)
    return rbin[~mask], densprof[~mask]/counts[~mask], counts[~mask]

def correlation_function_Pk( rho, rhoave ):

    N     = np.int64( np.ceil( rho.shape[0]**(1./3) ) )
    rho   = rho.reshape( (N, N, N) )
    rho   = (rho - rhoave)/rhoave
    Pk3D  = np.abs( np.fft.fftn( rho ) )**2

    k  = np.fft.fftfreq( N )
    kx,ky,kz = np.meshgrid( k, k, k )
    kx = kx.flatten()
    ky = ky.flatten()
    kz = kz.flatten()
    Pk3D = Pk3D.flatten()
    distance = np.sqrt( kx**2 +  ky**2 + kz**2 )

    N            = np.int64( 8e1 )
    rbin_edge    = np.logspace( -2, np.log10( distance.max() ),N )
    rbin         = (rbin_edge[1:]+rbin_edge[:-1])/2
    corrfunc     = np.zeros( N-1 )
    counts       = np.zeros( N-1 )

    inds  = np.digitize( distance, rbin_edge )
    inds -= 1

    for ind in np.unique( inds ):
        if ( ind > 0 and ind < inds.max() ):
            corrfunc[ind] += np.sum( (Pk3D)[inds==ind] )
            counts  [ind] += np.sum( inds==ind )

    mask = (counts==0)
    return rbin[~mask], corrfunc[~mask]/counts[~mask], counts[~mask]

def normalize_rho( x, y, z, rho, dr, dprof, center ):

    rhonorm = rho.copy()
    r = np.sqrt( (x - center[0])**2 +
                 (y - center[1])**2 +
                 (z - center[2])**2 )

    inds = np.digitize( r, dr )-1

    for i in range( inds.shape[0] ):
        rhonorm[i] /= dprof[inds[i]]

    return rhonorm

def Compute_PowerSpectrum( ds, center, rc, halo_id, path_data ):

    a      = 1/(1+ds.current_redshift)
    h      = ds.hubble_constant
    rc    *= h / 1e3 # from comoving kpc to comoving Mpc /h
    dx_max = (ds.domain_width[0]/ds.domain_dimensions[0] ).d  # in cMpc/h
    dxs    = ds.arr( [dx_max/(2**i) for i in np.arange(ds.max_level+1)], 'code_length' )  # in cMpc/h

    rmax = rfac * rc
    rmin = rmfac * rc
    cube = get_frb_cube( ds, center, rmax, lv, dxs )

    x    = cube["gamer","x"].to( "Mpc/h" ).d.flatten()
    y    = cube["gamer","y"].to( "Mpc/h" ).d.flatten()
    z    = cube["gamer","z"].to( "Mpc/h" ).d.flatten()
    rho  = cube["gas", "density"].to("Msun/kpc**3").d.flatten()

    dr,dprof,_ = density_profile( x, y, z, rho, center )
    rho        = normalize_rho( x, y, z, rho, dr, dprof, center )
    rho        = corezero( x, y, z, rho, rmin, center )
    rhoave     = np.mean( rho )

    k,pk,c  = correlation_function_Pk( rho, rhoave )
    pknorm  = pk * k**3
    pknorm /= pknorm.max()
    k       = k / dxs[-1] / 1e3 / a * h * np.pi # from 1/ pixel to 1/kpc == wave vector (cycle/dx) but not angular wave vector (2pi/dx)

    if not os.path.exists( path_data + '/powerspec_within_halo' ):
        os.makedirs( path_data + '/powerspec_within_halo' )

    with open( '%s/powerspec_within_halo/%s_%d_powerspectrum_within_halo'%(path_data,ds,halo_id) , 'w' ) as file:
        writer = csv.writer( file, delimiter='\t' )
        writer.writerow( [f"{'#wavenumber(2pi/kpc)':<22}", f"{'powerspec':<15}"] )
        for i in range( len(k) ):
            writer.writerow( [f"{k[i].d:<22.8f}", f"{pknorm[i]:<15.8f}"] )
