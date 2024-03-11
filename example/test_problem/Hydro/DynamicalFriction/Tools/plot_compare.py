import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from read_parameters import read_selected_parameters

Chandreasekhar_result_path = ''

x_min, x_max = 0, 20
y_min, y_max = 0.1, 13

def plot_rt():
    # 0. Read the data
    # 0.1 set up the desired resolution and particle number
    # 0.2 read the data
    folder = '../'
    df_center = pd.read_csv(folder+'Record__Result_Center',sep='\t', skiprows=1, header=None)
    df_center.columns = ['Time', 'CenterX', 'CenterY', 'CenterZ']
    df_GC = pd.read_csv(folder+'Record__Result_GC_position',sep='\t', skiprows=1, header=None)
    df_GC.columns = ['Time', 'GCX', 'GCY', 'GCZ']

    # 1. calculate the radius
    rx = df_GC['GCX'] - df_center['CenterX']
    ry = df_GC['GCY'] - df_center['CenterY']
    rz = df_GC['GCZ'] - df_center['CenterZ']

    r = np.sqrt(rx**2 + ry**2 + rz**2)
    df_GC['Radius'] = r/1e3 # convert to [kpc]

    # 2. set the coreect unit
    df_GC['Time'] = df_GC['Time']/1e3  # convert from [Myr] to [Gyr]

    params_IP = read_selected_parameters("../Record__Note")
    SPATIAL_RESO = params_IP['BOX_SIZE_X']/params_IP['NX0_TOT[0]'] # [pc]

    return df_GC['Radius'], df_GC['Time'],SPATIAL_RESO


def plot_chandrasekhar(filename):
    df = pd.read_csv(filename, sep='\t', header=None, comment='#')
    df.columns = ['Time', 'Radius']
    # Split the filename at underscores
    parts = filename.split("_")

    # Extract values using string split
    bminFactor = int(parts[0].split("=")[1])
    bmaxFactor = int(parts[1].split("=")[1])

    return df['Radius'], df['Time'],[bminFactor,bmaxFactor]

print("1. Read and plot the simulation result...")
result_simu = plot_rt()
plt.plot(result_simu[1], result_simu[0],label=f'Simulation       : Spatial Resolution={result_simu[2]:.3f} [pc]', linewidth=0.75)
print("1. Read and plot the simulation result...Done")
print("2. Read and plot the Chandrasekhar result...")
result_DF = plot_chandrasekhar(Chandreasekhar_result_path)
plt.plot(result_DF[1], result_DF[0],label=f'Chandrasekhar DF : bmin Factor={result_DF[2][0]}, bmax Factor={result_DF[2][1]}', linewidth=0.75)
print("2. Read and plot the Chandrasekhar result...Done")


plt.xlabel('Time (Gyr)')
plt.ylabel('Radius (kpc)')
plt.ylim(y_min, y_max)
plt.xlim(x_min, x_max)
plt.grid(True, 'both')
plt.yscale('log')
plt.title('Result Comparison between Simulation and Chandrasekhar DF')
plt.legend(loc='best')
plt.savefig('Compare_result.png',dpi=300)
print(f"3. Save the result figure : Compare_result.png")
# plt.show()
