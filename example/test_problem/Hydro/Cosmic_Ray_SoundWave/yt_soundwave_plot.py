import numpy as np
import matplotlib.pylab as plt
import yt

def plot_wave(num, start, end, waveV, waveLen, r_dir, axis, save_name):
    p_0     = 1.
    rho_0   = 1.
    c_s     = np.sqrt( ( (5./3.)*p_0 + (4./3.)*p_0 ) / rho_0 )
    perturb = 1.e-6
    
    ds = yt.load("./Data_%06d"%num)
    dd = ds.all_data()
    
    ray = ds.r[start:end]
    srt = np.argsort(ray[axis])
    
    dens   = np.array(ray["Dens"][srt])
    cray   = np.array(ray["CRay"][srt])
    cray_p = cray * (4./3. - 1)
    gas_p  = ( np.array(ray["Engy"][srt]) - cray ) * (5./3. - 1.)
    delta_cray_p = cray_p - p_0
    
    cray_sol = (4./3.) * perturb / c_s
    pres_sol = (5./3.) * perturb / c_s
    dens_sol =           perturb / c_s
    
    plot_x = np.array(ray[axis][srt])
    plot_y = np.array(ray[axis][srt])
    plot_z = np.array(ray[axis][srt])
    
    time = ds.current_time / ds.time_unit
    
    plot_r = plot_x + plot_y + plot_z
    axis_num = 3.

    if not r_dir[0]:
        plot_r -= plot_x
        axis_num -= 1.
    
    if not r_dir[1]:
        plot_r -= plot_y
        axis_num -= 1.
    
    if not r_dir[2]:
        plot_r -= plot_z
        axis_num -= 1.

    plot_r = plot_r / np.sqrt(axis_num) - waveV * time
    
    waveK   = 2. * np.pi / waveLen
    waveW   = waveK * c_s
    wave    = np.sin( plot_r*waveK - waveW * time )
    
    fig, ax = plt.subplots(1, 2, figsize = (12, 5))
    
    # Cosmic ray
    ax[0].scatter(plot_r, delta_cray_p / cray_p, s = 5, label = "Sim")
    ax[0].plot(plot_r, cray_sol * wave, '-r', label = "Analytical")
    ax[0].set(xlabel = "r", ylabel = r"$\delta p_{cr}/p_{cr}$", title = "cosmic ray")
    ax[0].set(ylim = [-1.e-6, 1.e-6])
    ax[0].legend(loc = 1)
    ax[0].ticklabel_format(axis = "y", style = "sci", scilimits = (4, 5))
    
    # density
    ax[1].scatter(plot_r, np.array(ray["density"][srt]) - rho_0, label = "Sim")
    ax[1].plot(plot_r, dens_sol * wave, '-r',  label = "Analytical")
    ax[1].set(xlabel = "r", ylabel = r"$\delta \rho /\rho_0$", title = "density")
    ax[1].set(ylim = [-1.e-6, 1.e-6])
    ax[1].legend(loc = 1)
    ax[1].ticklabel_format(axis = "y", style = "sci", scilimits = (4, 5))
    
    plt.suptitle("Sound Wave Test, t = %.3f"%time)
    if save_name != '':
        plt.savefig(save_name, dpi = 200)
    plt.close()
            
    return

#===============================================================================
# plot  wave
#===============================================================================
start, end = [0., 0., 0.], [1., 1., 1.]         # start and end line of plot

r_dir = [True, True, True] # x, y, z
axis_sort = 'x'                                 # the sort direction should be True at "r_dir"

waveV = 0.0                                     # background velocity
waveL = np.sqrt(3.) / 3.        # 0.5           # wavelength

max_num = 2
for i in range(max_num):
    save_name = './CR_soundwave_%04d.png'%i
    plot_wave(i, start, end, waveV, waveL, r_dir, axis_sort, save_name)

