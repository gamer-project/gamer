"""
NOTE: assuming the center is at (0, 0, 0)
"""
import numpy as np


#====================================================================================================
# Functions
#====================================================================================================
def ana_ball( x, y, z, t, **kwargs ):
    K_PARA = kwargs["K_PARA"]   # parallel diffusion coefficient
    K_PERP = kwargs["K_PERP"]   # perpendicular diffusion coefficient
    MAG    = kwargs["MAG"]      # magnetic field
    CR_BG  = kwargs["CR_BG"]    # cosmic ray energy density background value
    CR_E0  = kwargs["CR_E0"]    # cosmic ray energy dendity amplitude
    CR_R2  = kwargs["CR_R2"]    # inverse of variance of Gaussian distribution

    mag     = np.copy(MAG) / np.linalg.norm(MAG) # normalization
    mag_new = np.array([1., 0., 0.])

    # f(x, t) = 1 + 4 * R^2 * K_x * t
    f_x = 1. + (4. * CR_R2 * K_PARA * t)
    f_y = 1. + (4. * CR_R2 * K_PERP * t)
    f_z = 1. + (4. * CR_R2 * K_PERP * t)

    C   = np.cross( mag, mag_new )
    cos = np.dot( mag_new, mag )

    # rotation matrix
    R = np.array([[    1. - (C[1]**2 + C[2]**2) / (1+cos), -C[2] +       (C[0] * C[1]) / (1+cos),  C[1] +       (C[0] * C[2]) / (1+cos) ], \
                  [  C[2] +       (C[0] * C[1]) / (1+cos),    1. - (C[2]**2 + C[0]**2) / (1+cos), -C[0] +       (C[1] * C[2]) / (1+cos) ], \
                  [ -C[1] +       (C[0] * C[2]) / (1+cos),  C[0] +       (C[1] * C[2]) / (1+cos),    1. - (C[0]**2 + C[1]**2) / (1+cos) ]])

    # get the rotated coordinate
    pos     = np.array([x, y, z])
    pos_new = R.dot(pos)
    x_new = pos_new[0, :]
    y_new = pos_new[1, :]
    z_new = pos_new[2, :]

    # e(x, t) = e0 * exp( -R^2 * x^2/ f ) / sqrt(f)
    ex = np.exp( -CR_R2*x_new**2 / f_x ) / np.sqrt(f_x)
    ey = np.exp( -CR_R2*y_new**2 / f_y ) / np.sqrt(f_y)
    ez = np.exp( -CR_R2*z_new**2 / f_z ) / np.sqrt(f_z)

    return CR_BG + CR_E0 * ex * ey * ez


def ana_plane( x, y, z, t, **kwargs ):
    print("NOTE: Assuming uniform magnetic field and the direction is same as the gaussiain direction.")
    K_PARA = kwargs["K_PARA"]   # parallel diffusion coefficient
    CR_BG  = kwargs["CR_BG"]    # cosmic ray energy density background value
    CR_E0  = kwargs["CR_E0"]    # cosmic ray energy dendity amplitude
    CR_R2  = kwargs["CR_R2"]    # inverse of variance of Gaussian distribution
    DIR    = kwargs["DIR"]      # the direction of gaussian and magnetic field.

    if DIR == "x":
        r = x
    elif PLANE == "y":
        r = y
    elif PLANE == "z":
        r = z
    elif PLANE == "xy":
        r = x + y
    elif PLANE == "xz":
        r = x + z
    elif PLANE == "yz":
        r = y + z
    elif PLANE == "xyz":
        r = x + y + z
    else:
        raise ValueError( "DIR can only be one of [x/y/z/xy/xz/yz/xyz]. DIR: %s"%(DIR) )

    first = 1. + ( 4 * CR_R2 * K_PARA * t )

    return CR_BG + CR_E0 * np.exp( -CR_R2*r*r / first ) / np.sqrt(first)


def ana_step_ring( x, y, z, t, **kwargs ):
    print("NOTE: Assuming cicular magnetic field, and only parallel diffusion")
    K_PARA = kwargs["K_PARA"]   # parallel diffusion coefficient
    CR_BG  = kwargs["CR_BG"]    # cosmic ray energy density background value
    CR_E0  = kwargs["CR_E0"]    # cosmic ray energy dendity amplitude
    R_IN   = kwargs["R_IN"]     # inner radius of the ring
    R_OUT  = kwargs["R_OUT"]    # outer radius of the ring
    PLANE  = kwargs["PLANE"]    # the plane where the ring exist

    if PLANE == "xy":
        r   = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        N   = len(x)
    elif PLANE == "xz":
        r   = np.sqrt(z**2 + x**2)
        phi = np.arctan2(x, z)
        N   = len(y)
    elif PLANE == "yz":
        r   = np.sqrt(y**2 + z**2)
        phi = np.arctan2(z, y)
        N   = len(z)
    else:
        raise ValueError( "PLANE can only be one of [xy/xz/yz]. PLANE: %s"%(PLANE) )

    D  = np.sqrt(4. * t * K_PARA)
    t1 = np.empty(N)
    t2 = np.empty(N)
    for i in range(N):
        if r[i] < R_IN or r[i] > R_OUT:
            t1[i] = 0.
            t2[i] = 0.
        else:
            t1[i] = math.erfc( (phi[i] - np.pi/12.) * r[i] / D )
            t2[i] = math.erfc( (phi[i] + np.pi/12.) * r[i] / D )
    return CR_BR + t1 - t2


def ana_gaussian_ring( x, y, t ):
    print("NOTE: Assuming cicular magnetic field, and only parallel diffusion")
    K_PARA  = kwargs["K_PARA"]  # parallel diffusion coefficient
    CR_BG   = kwargs["CR_BG"]   # cosmic ray energy density background value
    CR_E0   = kwargs["CR_E0"]   # cosmic ray energy dendity amplitude
    R       = kwargs["R"]       # center radius of the ring
    DEL_R   = kwargs["DEL_R"]   # the standard deviation on r direction
    DEL_PHI = kwargs["DEL_PHI"] # the standard deviation on phi direction
    PLANE   = kwargs["PLANE"]   # the plane where the ring exist

    if PLANE == "xy":
        r   = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
    elif PLANE == "xz":
        r   = np.sqrt(z**2 + x**2)
        phi = np.arctan2(x, z)
    elif PLANE == "yz":
        r   = np.sqrt(y**2 + z**2)
        phi = np.arctan2(z, y)
    else:
        raise ValueError( "PLANE can only be one of [xy/xz/yz]. PLANE: %s"%(PLANE) )

    del_r2 = DEL_R**2
    del_phi2 = DEL_PHI**2 + 2 * K_PARA * t / r**2
    e0  = CR_E0 * np.sqrt(0.5**2 / del_phi2)

    tr  = np.exp( -0.5 * (r-R)**2 / del_r2)
    tp  = np.exp( -0.5 * phi**2 / del_phi2)
    return CR_BG + e0 * tr * tp
