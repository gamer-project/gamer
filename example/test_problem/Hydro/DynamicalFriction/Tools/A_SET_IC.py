import numpy as np
import scipy.integrate as integrate
import pandas as pd
import datetime
from X_some_function import read_parameters
import os

params = read_parameters("../Input__TestProb")

nbin = 5000
# Assign the parameters to variables
Center_Halo = params["Cloud_CenterX"], params["Cloud_CenterY"], params["Cloud_CenterZ"] #Units: pc, pc, pc
Halo = params["HALO_Rs"], params["HALO_RHO_0"], params["HALO_Rt"],params["HALO_TYPE"] #Units: pc, Msun/pc^3, pc, None
Stellar = params["STELLAR_Rs"], params["STELLAR_RHO_0"], params["STELLAR_Rt"],params["STELLAR_TYPE"]
GC_mass = params["GC_MASS"] #Units: Msun
GC_ri = params["GC_initial_R"] #Units: pc
r_cutoff = max(params["HALO_Rt"], params["STELLAR_Rt"]) #Units: pc

PURE_TABLE = params["PURE_TABLE"]
if PURE_TABLE: pure_table_path = params["Chandrasekhar_DF_Table_Name"]


G = 4.302e-3  # gravitational constant in pc^3 / (Myr^2 Msun)
# # print all the assigned parameters to check
# print("Center_Halo:", Center_Halo)
# print("Halo:", Halo)
# print("Stellar:", Stellar)
# print("GC_mass:", GC_mass)
# print("GC_ri:", GC_ri)
# print("r_cutoff:", r_cutoff)


def NFW_dens(x):
    return (1/((1+x)*(1+x)*x))

def Burkert_dens(x):
    return ((1/((1+x)*(1+x*x))))

def Plummer_dens(x):
    return (1+x*x)**(-2.5)

def LC_dens(x):
    gamma_0 = 0.07
    gamma_inf = 4.65
    eta = 3.7
    r0 = 1.4
    return x**(-gamma_0)*(1+x**eta)**((gamma_0-gamma_inf)/eta)

def King_dens(x):
    return ((1+x**2)**(-0.5))

def density(rho_s, r0, r, model_name):
    x = r/r0
    if model_name == "Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name == "NFW":
        return rho_s*NFW_dens(x)
    elif model_name == "Plummer":
        return rho_s*Plummer_dens(x)
    elif model_name == "LC":
        return rho_s*LC_dens(x)
    elif model_name == "King":
        return rho_s*King_dens(x)
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
        elif model_name == "LC":
            return 4*np.pi*x*x*(r0**3) * LC_dens(x)
        elif model_name == "King":
            return 4*np.pi*x*x*(r0**3) * King_dens(x)
        else:
            return 0
    f = integrate.nquad(massbase, [[0, x]])[0] * rho_s
    return f


def getrho(mass, r0, rt, model_name):
    if model_name == "None":
        return 0
    rho_0 = 3.0
    mass_init = clustermass(rho_0, r0, rt, model_name)
    rho_0 = rho_0*mass/mass_init
    return rho_0


def create_table(table_filename, r, dens, mass):
    f = open(table_filename, "w")
    f.write("#r[pc] density[Msun/pc^3] enclosed_mass[Msun]\n")

    for i in range(len(r)):
        if dens[i] != 0.:
            a = str(r[i])+" "+str(dens[i])+" "+str(mass[i]) + "\n"
            f.write(a)

    f.close()

def SetICs(Halo, Stellar, r_orbit, r_cutoff, center_position):
    
    center_x, center_y, center_z = center_position #[pc]
    halo_rs, halo_rho0, halo_rt, halo_type = Halo
    stellar_rs, stellar_rho0, stellar_rt, stellar_type = Stellar
    rho_halo = halo_rho0
    rho_stellar = stellar_rho0


    # create the table
    if PURE_TABLE == False:
        print('A-1-1. Create the table...')
        table_filename = Halo[3]+'_table_'+datetime.datetime.now().strftime("%Y%m%d")
        r = np.logspace(np.log10(halo_rs/100), np.log10(r_cutoff), nbin)
        if Stellar == (0, 0, 0, "None"):
            dens = np.array([density(rho_halo, halo_rs, r[i], halo_type) for i in range(len(r))])
            mass = np.array([clustermass(rho_halo, halo_rs, r[i], halo_type) for i in range(len(r))])
            create_table(table_filename+'.txt', r, dens, mass)
        else:
            dens = np.array([density(rho_halo, halo_rs, r[i], halo_type)+density(
                rho_stellar, stellar_rs, r[i], stellar_type) for i in range(len(r))])
            mass = np.array([clustermass(rho_halo, halo_rs, r[i], halo_type) + clustermass(
                rho_stellar, stellar_rs, r[i], stellar_type) for i in range(len(r))])
            create_table(table_filename+'.txt', r, dens, mass)
        print('A-1-1. Create the table...Done : '+table_filename+' is created')
    elif PURE_TABLE == True:
        print('A-1-2. Pure Table (Still need to specify rho0 and Rs manually in Input__Testprob)...') 
        table_filename, file_extension = os.path.splitext(pure_table_path)
        print('A-1-2. Pure Table...Done : '+table_filename+' is read')

    print('A-2. Extract the table for simulation need...')
    if PURE_TABLE == False:
        # Extract the table for the density profile
        with open('../'+table_filename+'_wo_mass.txt', 'w') as f:
            for i in range(len(r)):
                if i == 0:
                    f.write("#r[pc] density[Msun/pc^3]\n")
                if dens[i] != 0.:
                    f.write(str(r[i])+" "+str(dens[i])+"\n")
    elif PURE_TABLE == True:
        # Extract the table for the density profile
        with open(pure_table_path, 'r') as f:
            lines = f.readlines()

        with open('../'+table_filename+'_wo_mass.txt', 'w') as f:
            for line in lines:
                columns = line.split()
                f.write(f"{columns[0]} {columns[1]}\n")
    print('A-2. Extract the table for simulation need...Done : ../'+table_filename+'_wo_mass.txt is created')


    print("A-3 Calculate the initial conditions from table...")
    # Use the table to caculate the initial conditions
    data = np.loadtxt("%s" % table_filename+'.txt' , comments='#', dtype=float)
    Input_Radius = data[:, 0].astype(float) # in pc
    Input_Density_Scaled = data[:, 1].astype(float) # in Msun/pc^3
    Input_Mass_Profile = data[:, 2].astype(float) # in Msun

    def accumulated_mass_from_file(r):
        # find the index of the radius which is closest to the input radius
        index = np.abs(Input_Radius-r).argmin()
        accumulated_mass = Input_Mass_Profile[index]
        return accumulated_mass

    v_initial = (G*accumulated_mass_from_file(GC_ri)/GC_ri)**0.5
    theta = 0

    # set the GC position
    center_pos = [Center_Halo[0], Center_Halo[1], Center_Halo[2]]
    GC_pos = [center_pos[0]+GC_ri * np.cos(theta), center_pos[1]+GC_ri*np.sin(theta), center_pos[2]]

    # set the GC velocity (counterclockwise)
    halo_cutoff_mass = accumulated_mass_from_file(r_cutoff)
    GC_v = [-v_initial*np.sin(theta), v_initial*np.cos(theta), 0]


    
    if Stellar == (0, 0, 0, "None"):
        print("!!!!!!!!REMINDER : No stellar component!!!!!!!!")

    

    print('--------------------------------')
    print('scaling density:')
    print("DM Halo Density (rho0):{:>15.10f} [Msun/pc^3]".format( rho_halo))
    print("Stellar Density (rho0):{:>15.10f} [Msun/pc^3]".format( rho_stellar))
    print("Total Peak density    :{:>15.10f} [Msun/pc^3]".format( rho_halo + rho_stellar))

    # if Stellar == (0, 0, 0, "None"):
    #     enclosed_mass = clustermass(rho_halo, halo_rs, ri, halo_type)
    #     halo_cutoff_mass = clustermass(rho_halo, halo_rs, r_cutoff, halo_type)
    # else:
    #     enclosed_mass = clustermass(rho_halo, halo_rs, ri, halo_type) +  clustermass(rho_stellar, stellar_rs, ri, stellar_type)
    #     halo_cutoff_mass = clustermass(rho_halo, halo_rs, r_cutoff, halo_type)+clustermass(rho_stellar, stellar_rs, r_cutoff, stellar_type)
        
    
    # v = (G*enclosed_mass/ri)**0.5
    # theta = 0
    # print('--------------------------------')
    # print("Enclosed Mass from ri :{:>15.1f} [Msun]".format(enclosed_mass))
    # print('intial x,y:')
    print("GC mass               :{:>15.1f} [Msun]".format(GC_mass))
    print("GC x                  :{:>15.1f} [pc]".format(GC_pos[0]))
    print("GC y                  :{:>15.1f} [pc]".format(GC_pos[1]))
    print("GC z                  :{:>15.1f} [pc]".format(GC_pos[2]))
    # print('intial velocity (No Consider relative speed between halo and GC):')
    print("GC Vx                 :{:>15.10f} [pc/Myr]".format(GC_v[0]))
    print("GC Vy                 :{:>15.10f} [pc/Myr]".format(GC_v[1]))
    print("Halo Vx               :{:>15.10f} [pc/Myr]".format(0))
    print("Halo Vy               :{:>15.10f} [pc/Myr]".format(0))
    print("Halo,GC Vz            :{:>15.10f} [pc/Myr]".format(0))
    print('--------------------------------')
    print("A-3. Calculate the initial conditions...Done")
    

    
    # update the input file for the result:
    print('A-4. Update Input__TestProb...')
    with open("../Input__TestProb", "r") as f:
        lines = f.readlines()

    with open("../Input__TestProb", "w") as f:
        for line in lines:
            if line.startswith("GC_POSX"):
                f.write("GC_POSX               "+str(GC_pos[0])+"                       #+\n")
            elif line.startswith("GC_POSY"):
                f.write("GC_POSY               "+str(GC_pos[1])+"                       #+\n")
            elif line.startswith("GC_POSZ"):
                f.write("GC_POSZ               "+str(GC_pos[2])+"                       #+\n")
            elif line.startswith("GC_VELX"):
                f.write("GC_VELX               "+str(GC_v[0])+"                          #+\n")
            elif line.startswith("GC_VELY"):
                f.write("GC_VELY               "+str(GC_v[1])+"            #+\n")
            elif line.startswith("GC_VELZ"):
                f.write("GC_VELZ               0.0"+"                           #+\n")
            elif line.startswith("Density_Table_Name"):
                f.write("Density_Table_Name     "+table_filename+'_wo_mass.txt'+" #+ The table only contains radius and density for sumulation need\n")
            elif line.startswith("Chandrasekhar_DF_Table_Name"):
                f.write("Chandrasekhar_DF_Table_Name   "+table_filename+'.txt'+" #+ if PURE_TABLE = 0; *if PURE_TABLE=1 the table must contains radius, density and enclosed mass\n")
            elif line.startswith("Cloud_BulkVelX"):
                f.write("Cloud_BulkVelX        "+str(0)+"                             #+ Bulk Velocity for the whole cloud\n")
            elif line.startswith("Cloud_BulkVelY"):
                f.write("Cloud_BulkVelY        "+str(0)+"                             #+ Bulk Velocity for the whole cloud\n")
            elif line.startswith("Cloud_BulkVelZ"):
                f.write("Cloud_BulkVelZ        0.0"+"                           #+ Bulk Velocity for the whole cloudn\n")
            elif line.startswith("Cloud_Rho0"):
                f.write("Cloud_Rho0            "+str(rho_halo+rho_stellar)+"            #+ Maximum core density for the cloud : set to be HALO_RHO_0+STELLAR_RHO_0 as default\n")
            elif line.startswith("Cloud_MaxR"):
                f.write("Cloud_MaxR            "+str(r_cutoff)+"                        #+Maximum radius for the particle : set to be HALO_Rt as default"+"\n")
            elif line.startswith("Cloud_R0"):
                f.write("Cloud_R0              "+str(Halo[0])+"                         #+ Maximum core radius for the cloud : set to be HALO_Rs as default \n")
            else:
                f.write(line)
    print('A-4. Update Input__TestProb...Done')


SetICs(Halo, Stellar, GC_ri, r_cutoff, Center_Halo)

