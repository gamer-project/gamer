import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.special import erf
from tqdm import tqdm
from read_parameters import read_selected_parameters,extract_parameters_from_logfile

print("==========Calculating the dynamical friction using Chandrasekhar's DF formula==========")
version = 'interpolate_table' # 'analytical' or 'interpolate_table' -> 'analytical' means using the analyical result "without" any hard truncation, 'interpolate_table' means using the table generated in GAMER routine

LN_LAMBDA_BMAX_FACTOR = 3  # BMAX fitting factor for ln_Lambda, 3 is the default value meaning b_max = r*3 (r is the distance to the center)
LN_LAMBDA_BMIN_FACTOR = 3  # BMIN fitting factor for ln_Lambda, 3 is the default value meaning b_min = dc*3 (dc is the spatial resolution in simulation)

# Plotting and saving the result
SAVE_PLOT_RESULT = True
SAVE_PLOT_LOG_SCALE = True
WRITE_TO_FILE = True

DATA_SOURE = 'SIMU' # 'USER' or 'SIMU' : ['SIMU']. If you choose 'USER', you need to set the parameters below!!!

if DATA_SOURE == 'SIMU':
    # read the data from simulation : spatial resolution and GC_mass
    params_TB =extract_parameters_from_logfile("../log")
    m = params_TB['GC Mass']     # mass of the massive object in [Msun]
    GC_velocity = params_TB['GC Velocity X'], params_TB['GC Velocity Y'], params_TB['GC Velocity Z'] # [pc/Myr] = [kpc/Gyr]
    GC_position = params_TB['GC Position X']/1000, params_TB['GC Position Y']/1000, params_TB['GC Position Z']/1000 # [kpc]
    Halo_Center = params_TB['central coordinate [0]']/1000, params_TB['central coordinate [1]']/1000, params_TB['central coordinate [2]']/1000 # [kpc]
    params_IP = read_selected_parameters("../Record__Note")
    SPATIAL_RESO = params_IP['BOX_SIZE_X']/params_IP['NX0_TOT[0]'] # [pc]

    DPLD_a, DPLD_b, DPLD_c = params_TB['Halo_Profle_Param_a'], params_TB['Halo_Profle_Param_b'], params_TB['Halo_Profle_Param_c']

    t_end = params_IP['END_T']/1000 # end time in [Gyr]
    Interpolate_table_path = '../Profile_Table'
    halo_rs,rho0_halo, halo_type = params_TB['HALO_Rs']/1000, params_TB['HALO_RHO_0']*1E9,params_TB['Halo Type']# [kpc], [Msun/kpc^3], Type
    stellar_rs, rho0_stellar, stellar_type = 0, 0, "None" # [kpc], [Msun/kpc^3], Type
elif DATA_SOURE == 'USER' :
    m = 0 # mass of the massive object in [Msun]
    SPATIAL_RESO = 0 # [pc]
    t_end = 0 # end time in [Gyr]
    GC_position = 0,0,0 # [kpc]
    GC_velocity = 0,0,0 # [pc/Myr]=[kpc/Gyr]
    Halo_Center = 0,0,0 # [kpc]
    Interpolate_table_path = 'PUT_YOUR_PATH_HERE'
    halo_rs,rho0_halo, halo_type = 0, 0, "0" # [kpc], [Msun/kpc^3], Type
    stellar_rs, rho0_stellar, stellar_type = 0, 0, "0" # [kpc], [Msun/kpc^3], Type
else:
    print('Wrong DATA_SOURE!!! Please set the DATA_SOURE to be "USER" or "SIMU"')




GC_softening = 12.5  # softening length in [pc] -> For softening the gravitational acceleration
dt = 0.01   # time step in [Gyr]
G = 4.492e-6  # gravitational constant in [kpc^3 / (Gyr^2 Msun)]
UPDATE_METHOD = 'RK4'  # ['RK4'] : 'RK4', 'Euler' or 'RK2'


def NFW_dens(x):
    return (1/((1+x)*(1+x)*x))


def Burkert_dens(x):
    return ((1/((1+x)*(1+x*x))))


def Plummer_dens(x):
    return (1+x*x)**(-2.5)


def King_dens(x):
    return ((1+x**2)**(-0.5))

def DoublePowerLaw_dens(x,alpha,beta,gamma):
    return x**(-gamma)*(1+x**alpha)**((gamma-beta)/alpha)


def density(rho_s, r0, r, model_name):
    x = r/r0
    if model_name == "Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name == "NFW":
        return rho_s*NFW_dens(x)
    elif model_name == "Plummer":
        return rho_s*Plummer_dens(x)
    elif model_name == "King":
        return rho_s*King_dens(x)
    elif model_name == "DoublePowerLaw":
        return rho_s*DoublePowerLaw_dens(x,DPLD_a,DPLD_b,DPLD_c)
    else:
        return 0


def clustermass(rho_s, r0, r, model_name):

    x = r/r0

    def massbase(x):
        if model_name == "Burkert":
            return 4*np.pi*x*x*(r0**3) * Burkert_dens(x)
        elif model_name == "NFW":
            return 4*np.pi*x*x*(r0**3) * NFW_dens(x)
        elif model_name == "Plummer":
            return 4*np.pi*x*x*(r0**3) * Plummer_dens(x)
        elif model_name == "King":
            return 4*np.pi*x*x*(r0**3) * King_dens(x)
        elif model_name == "DoublePowerLaw":
            return 4*np.pi*x*x*(r0**3) * DoublePowerLaw_dens(x,DPLD_a,DPLD_b,DPLD_c)
        else:
            return 0
    clustermass_result = integrate.nquad(massbase, [[0, x]])[0] * rho_s
    return clustermass_result


def getrho(mass, r0, rt, model_name):
    if model_name == "None":
        return 0
    rho_0 = 3.0
    mass_init = clustermass(rho_0, r0, rt, model_name)
    rho_0 = rho_0*mass/mass_init
    return rho_0


def accumulated_mass_analytical(r, rs, rho_0):
    return 4*np.pi*rs**3*rho_0*((r+rs)*np.log(r+rs)+(-np.log(rs)-1)*r-rs*np.log(rs))/(r+rs)

def accumulated_mass_from_file(r):
    # find the index of the radius which is closest to the input radius
    index = np.abs(Input_Radius-r).argmin()
    accumulated_mass = Input_Mass_Profile[index]
    return accumulated_mass  # in [Msun]


def density_from_file(r):
    # find the index of the radius which is closest to the input radius
    index = np.abs(Input_Radius-r).argmin()
    density = Input_Density_Scaled[index]
    return density  # in [Msun/kpc^3]

def velocity_dispersion(r):
    if version == 'analytical':
        rho_r = density(rho0_halo, halo_rs, np.linalg.norm(r), halo_type) + density(rho0_stellar, stellar_rs, np.linalg.norm(r), stellar_type)
        integral_part = integrate.quad(lambda r_: (clustermass(rho0_halo, halo_rs, np.linalg.norm(r_), halo_type)+clustermass(rho0_stellar, stellar_rs, np.linalg.norm(r_), stellar_type)) *
                                       (density(rho0_halo, halo_rs, np.linalg.norm(
                                           r_), halo_type)+density(rho0_stellar, stellar_rs, np.linalg.norm(r_), stellar_type))
                                       * G/np.linalg.norm(r_)**2, r, np.infty)
        result = integral_part[0]/rho_r
    elif version == 'interpolate_table':
        rho_r = density_from_file(np.linalg.norm(r))
        integral_part = integrate.quad(lambda r_: accumulated_mass_from_file(np.linalg.norm(
            r_)) * density_from_file(np.linalg.norm(r_)) * G/np.linalg.norm(r_)**2, r, np.infty)
        result = integral_part[0]/rho_r
    else:
        return 0

    return result


def dynamical_friction_Chandrasekhar_erf(r, v, m):
    sigma = np.sqrt(velocity_dispersion(np.linalg.norm(r)))
    x = np.linalg.norm(v) / (np.sqrt(2) * sigma)

    if version == 'analytical':
        rho = density(rho0_halo, halo_rs, np.linalg.norm(r), halo_type) + density(rho0_stellar, stellar_rs, np.linalg.norm(r), stellar_type)
    elif version == 'interpolate_table':
        rho = density_from_file(np.linalg.norm(r))
    else:
        rho = 0
        print('Wrong rho!!!!')

    ln_Lambda = -999
    Lamda = (np.linalg.norm(r)*LN_LAMBDA_BMAX_FACTOR)/(SPATIAL_RESO*0.001*LN_LAMBDA_BMIN_FACTOR)
    ln_Lambda = np.log(Lamda)
    if ln_Lambda < 0:
        ln_Lambda = np.log(1+Lamda**2)
        print('Warning : ln_Lambda < 0!!! Correct to ln_Lambda = ln(1+Lamda^2) : ', ln_Lambda)

    erf_term = erf(x)
    exp_term = 2 * x / np.sqrt(np.pi) * np.exp(-x**2)
    a_df = -4*np.pi*G**2 * m * rho * ln_Lambda / np.linalg.norm(v)**3 * v * (erf_term-exp_term)
    return a_df

r0 = np.array([GC_position[0]-Halo_Center[0],GC_position[1]-Halo_Center[1],GC_position[2]-Halo_Center[2]], dtype=np.float64)  # initial position in [kpc]
v0 = np.array(GC_velocity, dtype=np.float64)  # initial velocity in [kpc/Gyr]

def calculate_trajectory(r0_=r0, v0_=v0, m=m, t_end=t_end, dt_=dt):
    r = r0_
    v = v0_
    dt = dt_
    t = 0
    r_history = [np.linalg.norm(r)]
    t_history = [float(t)]
    # Define the total number of iterations
    total_iterations = int(t_end / dt_)
    with tqdm(total=total_iterations) as pbar:
        while t < t_end and r_history[-1] > 5e-2: # 0.05 is a threshold to stop the calculation
            def a(r, v):
                # Calculate acceleration due to dynamical friction
                df_accel = dynamical_friction_Chandrasekhar_erf(r, v, m)
                if version == 'analytical':
                    # Calculate acceleration due to gravity and consider the softening length
                    g_accel = -G * (clustermass(rho0_halo, halo_rs, np.linalg.norm(r), halo_type) +
                                    clustermass(rho0_stellar, stellar_rs, np.linalg.norm(r), stellar_type))/(np.linalg.norm(r)+GC_softening*0.001)**3 * r
                elif version == 'interpolate_table':
                    g_accel = -G * (accumulated_mass_from_file(np.linalg.norm(r))) / (np.linalg.norm(r)+GC_softening*0.001)**3 * r
                else:
                    print('g_accel wrong!!!')
                    g_accel = 0
                    # return
                return df_accel+g_accel
            if UPDATE_METHOD == 'Euler':
                # Update velocity and position
                v += a(r, v) * dt
                r += v * float(dt)
            # apply RK2 for the position and the velocity
            if UPDATE_METHOD == 'RK2':
                y = [r, v]

                def evaluate(t, y):
                    v = y[1]
                    a_ = a(y[0], y[1])
                    return np.array([v, a_])
                k1 = evaluate(t, y)*dt
                k2 = evaluate(t+dt, y+k1)*dt
                y += (k1+k2)/2
                r = y[0]
                v = y[1]


            if UPDATE_METHOD == 'RK4':
                # using RK4 method to update
                y = [r, v]

                def evaluate(t, y):
                    v = y[1]
                    a_ = a(y[0], y[1])
                    return np.array([v, a_])
                k1 = evaluate(t, y)*dt
                k2 = evaluate(t+0.5*dt, y+0.5*k1)*dt
                k3 = evaluate(t+0.5*dt, y+0.5*k2)*dt
                k4 = evaluate(t+dt, y+k3)*dt
                y += (k1+2*k2+2*k3+k4)/6
                r = y[0]
                v = y[1]

            # Update time and save history
            t += dt

            r_mag = np.linalg.norm(r)

            r_history.append(r_mag)
            t_history.append(t)
            # Update the progress bar
            pbar.update(1)


    return np.array(r_history), np.array(t_history)


if version == 'interpolate_table':
    print(f'You are using the interpolation version, will use the data from {Interpolate_table_path}!')
    data = np.loadtxt("%s" % Interpolate_table_path, skiprows=1, dtype=float)
    Input_Radius = data[:, 0].astype(float)  # in [pc]
    Input_Radius = Input_Radius/1000  # change to [kpc]
    Input_Density_Scaled = data[:, 1].astype(float)  # in Msun/pc^3
    Input_Density_Scaled = Input_Density_Scaled*1e9  # change to [Msun/kpc^3]
    Input_Mass_Profile = data[:, 2].astype(float)  # in Msun
elif version == 'analytical':
    print('You are using the analytical version please make sure you have set the parameters correctly!')
else:
    print('Wrong version!!! Please set the version to be "analytical" or "interpolate_table"')

r_history, t_history = calculate_trajectory(r0_=r0, v0_=v0, m=m, t_end=t_end, dt_=dt)


if WRITE_TO_FILE:
    filename = f'bminFactor={LN_LAMBDA_BMIN_FACTOR}_bmaxFactor={LN_LAMBDA_BMAX_FACTOR}'
    with open(filename, 'w') as f:
        f.write(f'# t_end={t_end}, dt={dt}\n')
        f.write(f'# r0={r0}, v0={v0}\n')
        f.write(f'# t r\n')
        for i in range(len(r_history)):
            f.write(f'{t_history[i]:.10f}\t{r_history[i]:.10f}\n')
    print(f'Write to file Done! : {filename}')

if SAVE_PLOT_RESULT:
    print("5. Plotting and showing the result")
    plt.plot(t_history, r_history, label=rf'Chandrasekhar DF', linewidth=0.8)
    plt.xlim(0, t_end+1)
    plt.gcf().set_size_inches(8, 6)
    plt.legend()
    # open the grid line
    plt.xlabel('Time (Gyr)')
    plt.ylabel('Radius (kpc)')
    plt.grid(True)
    plt.title(f'Chandrasekhar DF, bminFactor={LN_LAMBDA_BMIN_FACTOR}_bmaxFactor={LN_LAMBDA_BMAX_FACTOR} profile')
    if SAVE_PLOT_LOG_SCALE: plt.yscale('log')
    plt.savefig(f'bminFactor={LN_LAMBDA_BMIN_FACTOR}_bmaxFactor={LN_LAMBDA_BMAX_FACTOR}.png', dpi=300)
    # plt.show()
    print('Save the plot result...Done! : '+f'bminFactor={LN_LAMBDA_BMIN_FACTOR}_bmaxFactor={LN_LAMBDA_BMAX_FACTOR}.png')
print('========================================================')
