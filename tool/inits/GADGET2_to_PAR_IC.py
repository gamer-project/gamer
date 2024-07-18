#!/usr/bin/env python3.9
import h5py
import numpy as np
import math

# user-specified parameters
input_GADGET_file      = "Data_000.hdf5" # GADGET-2 HDF5 file
                                         # --> an example "Data_000.hdf5" (128^3 particles) can be downloaded via "curl -JO http://use.yt/upload/1d18ed50"
output_float_precision = np.float32      # "np.float32" or "np.float64"
Gadget_UnitMass        = 1.989e43        # identical to the values adopted under "% System of units" in GADGET-2's runtime parameter file
Gadget_UnitLength      = 3.085678e24     # ...
Gadget_UnitVelocity    = 100000.0        # ...

# load the input file
with h5py.File(input_GADGET_file, "r") as f:
    Header_list = f["Header"]
    OMEGA_M0    = float(Header_list.attrs["Omega0"])        # omega matter at the present time
    HUBBLE0     = float(Header_list.attrs["HubbleParam"])   # dimensionless Hubble parameter "h"
    print("GADGET input OMEGA_M0: %.16f"%OMEGA_M0)
    print("GADGET input HUBBLE0:  %.16f"%HUBBLE0)
    print("[IMPORTANT] Set OMEGA_M0 and HUBBLE0 in GAMER Input__Parameter accordingly")

# "GAMER" physical constants (/include/PhysicalConstant.h), input parameters (Input__Parameter), and COMOVING units (src/Init/Init_Unit.cpp)
Const_cm      = 1.0
Const_km      = 1.0e5*Const_cm
Const_pc      = 3.08567758149e18*Const_cm                # parsec
Const_Mpc     = 1.0e6*Const_pc
Const_s       = 1.0                                      # second
Const_yr      = 3.15569252e7*Const_s                     # year
Const_g       = 1.0
Const_Msun    = 1.9885e33                                # solar mass
Const_NewtonG = 6.6738e-8                                # gravitational constant in cgs
H0            = 100*HUBBLE0*Const_km/(Const_s*Const_Mpc) # H0 = 100*h*km/(s*Mpc)
# see https://github.com/gamer-project/gamer/wiki/Runtime-Parameters%3A-Units#units-in-cosmological-simulations
GAMER_UNIT_L  = HUBBLE0**(-1.0)*Const_Mpc
GAMER_UNIT_T  = H0**(-1.0)
GAMER_UNIT_D  = 3.0*OMEGA_M0*(H0**2.0)/(8.0*math.pi*Const_NewtonG)

# GADGET-to-GAMER mass, position, and velocity unit conversions
mass_conversion = (Gadget_UnitMass/HUBBLE0)/(GAMER_UNIT_D*GAMER_UNIT_L**3.0)
pos_conversion  = (Gadget_UnitLength/HUBBLE0)/GAMER_UNIT_L
vel_conversion  = Gadget_UnitVelocity/(GAMER_UNIT_L/GAMER_UNIT_T)
print("GADGET-to-GAMER mass     unit conversion: %.16f"%mass_conversion)
print("GADGET-to-GAMER length   unit conversion: %.16f"%pos_conversion)
print("GADGET-to-GAMER velocity unit conversion: %.16f"%vel_conversion)

# load the input file
with h5py.File(input_GADGET_file, "r") as f:
    '''
    # print all root level object names/keys
    print("Keys: %s" % f.keys())
    input_key_list = list(f.keys())
    '''
    # print all attributes associated with "Header"
    Header_list                = f["Header"]
    Header_list_attribute_list = f["Header"].attrs.keys()
    input_a_scale_factor       = Header_list.attrs["Time"]
    print(Header_list_attribute_list)

    # output simulation parameters adopted in the GADGET input file
    with open('GADGET2_Input_Parameters.txt', 'w') as output:

        output.write("GADGET2 Input Parameters\n\n")

        # loop through the complete attribute information
        for Header_key_i in Header_list_attribute_list:

            # output attribute name
            output.write(format("%s:"%Header_key_i,"<25"))

            # output attribute value(s) as string(s)
            if isinstance(Header_list.attrs[Header_key_i],(list,np.ndarray)):
                str_place_holder ='['+",  ".join(map(str, Header_list.attrs[Header_key_i]))+"]"
            else: str_place_holder = str(Header_list.attrs[Header_key_i])

            output.write(str_place_holder)
            output.write("\n")

        output.close()

    print('GADGET2_Input_Parameters.txt complete')

    # start sorting the particle attributes and initial conditions
    # (1) ordinary point masses should have non-zero entry in the "MassTable"
    GADGET_MassTable = list(Header_list.attrs["MassTable"])
    # (2) "NumPart_ThisFile" should be identical to "NumPart_Total" IFF "input_GADGET_file" is NOT split over several smaller files
    GADGET_NumPart_ThisFile = list(Header_list.attrs["NumPart_ThisFile"])
    # (3) record the number of particles for each type
    GADGET_NumPart_Total = list(Header_list.attrs["NumPart_Total"])

    with open('PAR_IC', 'wb') as output:
        mass = np.array([], dtype=output_float_precision)
        posx = np.array([], dtype=output_float_precision)
        posy = np.array([], dtype=output_float_precision)
        posz = np.array([], dtype=output_float_precision)
        velx = np.array([], dtype=output_float_precision)
        vely = np.array([], dtype=output_float_precision)
        velz = np.array([], dtype=output_float_precision)
        partype = np.array([], dtype=output_float_precision)

        for i, NumPart_i in enumerate(GADGET_NumPart_Total):
            # if for "PartType(i)" the total particle number is greater than zero, we import the information to PAR_IC
            if NumPart_i > 0:
                key_NOW = "PartType"+str(i)
                PartType_Coordinates = np.array(f[key_NOW]["Coordinates"], float)
                Velocities_Coordinates = np.array(f[key_NOW]["Velocities"], float)

                if (GADGET_MassTable[i] == 0.0):
                    mass_i = np.array(f[key_NOW]["Masses"], float)*mass_conversion
                else:
                    mass_i = np.full(NumPart_i, GADGET_MassTable[i])*mass_conversion

                # set the origin of Cartesian coordinates to the box left bottom corner
                posx_i = PartType_Coordinates[:,0]*pos_conversion
                posy_i = PartType_Coordinates[:,1]*pos_conversion
                posz_i = PartType_Coordinates[:,2]*pos_conversion

                # For GADGET, v_peculiar = sqrt(a) * v_snapshot 
                # --> See "3. Particle velocities" on p.32 of https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf
                # For GAMER, v_code = a * v_physical
                # --> See Eq. (16) of Schive et al. (2010)
                velx_i = Velocities_Coordinates[:,0]*vel_conversion*(input_a_scale_factor**1.5)
                vely_i = Velocities_Coordinates[:,1]*vel_conversion*(input_a_scale_factor**1.5)
                velz_i = Velocities_Coordinates[:,2]*vel_conversion*(input_a_scale_factor**1.5)

                # For GADGET, see "Table 3" on p.31 of https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf
                # For GAMER, we fix the imported particle type to "1 (PTYPE_GENERIC_MASSIVE)."
                # --> See https://github.com/gamer-project/gamer/wiki/Initial-Conditions#IC-File-Particles
                partype_i = np.full(NumPart_i, 1)

                mass = np.append(mass, mass_i)
                posx = np.append(posx, posx_i)
                posy = np.append(posy, posy_i)
                posz = np.append(posz, posz_i)
                velx = np.append(velx, velx_i)
                vely = np.append(vely, vely_i)
                velz = np.append(velz, velz_i)
                partype = np.append(partype, partype_i)

                print('Number of particles = %i'%NumPart_i)

        output.write(mass.astype(dtype=output_float_precision).tobytes())
        output.write(posx.astype(dtype=output_float_precision).tobytes())
        output.write(posy.astype(dtype=output_float_precision).tobytes())
        output.write(posz.astype(dtype=output_float_precision).tobytes())
        output.write(velx.astype(dtype=output_float_precision).tobytes())
        output.write(vely.astype(dtype=output_float_precision).tobytes())
        output.write(velz.astype(dtype=output_float_precision).tobytes())
        output.write(partype.astype(dtype=output_float_precision).tobytes())

        output.close()

    print('PAR_IC complete')
