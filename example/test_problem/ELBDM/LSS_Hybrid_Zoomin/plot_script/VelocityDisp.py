#######################################

# Computation of Velocity Dispersion.
# Run Compute_profile.py to call these functions.

# carefully adjust rmax_for_vd to avoid exceeeding memory limit
#######################################
import os
import yt
import numpy as np
import csv
from   ELBDM_DerivedField import *

rmax_for_vd = 0.2 # in comoving Mpc/h. Maximum enclosing radius for computation of velocity dispersion

def get_frb_cube( ds, c, radius, lv, dxs ):

    L       = radius*2
    frb_res = np.int64( L/dxs[lv] )
    cedge   = c - radius
    cube    = ds.covering_grid( lv, left_edge=cedge, dims=[frb_res, frb_res, frb_res],num_ghost_zones=1 )

    return cube

def density_profile( x, y, z, rho, center ):

    r = np.sqrt( (x - center[0])**2 +
                 (y - center[1])**2 +
                 (z - center[2])**2 )

    N            = np.int64( 128 )
    rbin_edge    = np.logspace( -4, np.log10(r.max()), N )
    rbin         = (rbin_edge[1:]+rbin_edge[:-1])/2
    densprof     = np.zeros( N-1 )
    counts       = np.zeros( N-1 )

    inds  = np.digitize(r, rbin_edge)
    inds -= 1

    for ind in np.unique(inds):
        if (ind > 0 and ind < inds.max()):
            densprof[ind] += np.sum((rho)[inds==ind])
            counts  [ind] += np.sum(inds==ind)

    mask = (counts==0)

    return rbin[~mask], densprof[~mask]/counts[~mask], counts[~mask]

def Hz( z, H0, OmegaL, Omegam ):
    return H0 * np.sqrt( Omegam* (1+z)**3 + OmegaL )

def Compute_VelocityDispersion( ds, center, halo_id, path_data ):

    omegaL = ds.omega_lambda
    omegam = ds.omega_matter
    a      = 1/(1 + ds.current_redshift)
    z      = ds.current_redshift
    h      = ds.hubble_constant
    H      = Hz( z, h*100, omegaL, omegam )
    dx_max = (ds.domain_width[0]/ds.domain_dimensions[0]).d  # in cMpc/h
    dxs    = ds.arr( [dx_max/(2**i) for i in np.arange(ds.max_level+1)], 'code_length' )  # in cMpc/h
    lv     = 7

    Add_derived_fields( ds )

    cube = get_frb_cube( ds, center, rmax_for_vd, lv, dxs )

    x    = cube["gamer","x"].to( "Mpc/h" ).d.flatten()
    y    = cube["gamer","y"].to( "Mpc/h" ).d.flatten()
    z    = cube["gamer","z"].to( "Mpc/h" ).d.flatten()
    rho  = cube["gas", "density"].to( "Msun/kpc**3" ).d.flatten()
    vx   = cube["gamer", "v_x"].to( "km/s" ).d.flatten() / a + H * a * x # physical v = comoving v * a +  Hax_c
    vy   = cube["gamer", "v_y"].to( "km/s" ).d.flatten() / a + H * a * y # remind: v = \grad S = dS/dy = dS/dy_c * a^(-1)
    vz   = cube["gamer", "v_z"].to( "km/s" ).d.flatten() / a + H * a * z

    vx_center = np.sum( rho*vx )/np.sum( rho )
    vy_center = np.sum( rho*vy )/np.sum( rho )
    vz_center = np.sum( rho*vz )/np.sum( rho )

    vx -= vx_center; # rest frame at center of subhalo
    vy -= vy_center;
    vz -= vz_center;

    rhovx2  = rho * vx*vx
    rhovy2  = rho * vy*vy
    rhovz2  = rho * vz*vz
    rhovx   = rho * vx
    rhovy   = rho * vy
    rhovz   = rho * vz

    dr,dprof,_      = density_profile( x, y, z, rho,    center )
    dr,rhovx2prof,_ = density_profile( x, y, z, rhovx2, center )
    dr,rhovy2prof,_ = density_profile( x, y, z, rhovy2, center )
    dr,rhovz2prof,_ = density_profile( x, y, z, rhovz2, center )
    dr,rhovxprof,_  = density_profile( x, y, z, rhovx,  center )
    dr,rhovyprof,_  = density_profile( x, y, z, rhovy,  center )
    dr,rhovzprof,_  = density_profile( x, y, z, rhovz,  center )

    rhovxprof2 = rhovxprof**2
    rhovyprof2 = rhovyprof**2
    rhovzprof2 = rhovzprof**2

    sigmax2 = rhovx2prof/dprof - rhovxprof2/(dprof*dprof)
    sigmay2 = rhovy2prof/dprof - rhovyprof2/(dprof*dprof)
    sigmaz2 = rhovz2prof/dprof - rhovzprof2/(dprof*dprof)

    sigma = np.sqrt( (sigmax2 + sigmay2 + sigmaz2) / 3.0 )

    dr  *= 1e3/h # in ckpc
    if not os.path.exists( path_data + '/prof_veldisp' ):
       os.makedirs( path_data + '/prof_veldisp' )

    with open( '%s/prof_veldisp/%s_%d_veldisp_haloRestFrame'%(path_data,ds,halo_id) , 'w' ) as file:
        writer = csv.writer( file, delimiter='\t' )
        writer.writerow( [f"{'#radius(ckpc)':<15}", f"{'VelDisp(km/s)':<15}"] )
        for i in range( len(dr) ):
            writer.writerow( [f"{dr[i]:<15.8f}", f"{sigma[i]:<15.8f}"] )
