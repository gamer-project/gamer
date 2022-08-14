#!/usr/bin/env python3.9
import h5py
import numpy as np
import math

# "GADGET" input hdf5 file
input_GADGET_file          = "Data_000.hdf5" # an example "Data_000.hdf5" (128^3 particles) can be downloaded via "curl -JO http://use.yt/upload/1d18ed50"
conversion_float_percision = np.float32      # "np.float32" or "np.float64"
input_GADGET_Units_Option  = 0               # 0 -> GADGET default; non-zero -> manual input
if (input_GADGET_Units_Option != 0):
    Gadget_UnitMass     = 1
    Gadget_UnitLength   = 1
    Gadget_UnitVelocity = 1

# "GAMER" physical constants (/include/PhysicalConstant.h), input parameters (Input__Parameter), COMOVING units (src/Init/Init_Unit.cpp)
Const_cm            = 1.0
Const_km            = 1.0e5*Const_cm
Const_pc            = 3.08567758149e18*Const_cm                # parsec
Const_Mpc           = 1.0e6*Const_pc
Const_s             = 1.0                                      # second
Const_yr            = 3.15569252e7*Const_s                     # year
Const_g             = 1.0
Const_Msun          = 1.9885e33                                # solar mass
Const_NewtonG       = 6.6738e-8                                # gravitational constant in cgs
GAMER_OMEGA_M0      = 0.3158230904284232                       # omega matter at the present time
GAMER_HUBBLE0       = 0.6732117                                # dimensionless Hubble parameter "h"
GAMER_H0            = 100*GAMER_HUBBLE0*Const_km/(Const_s*Const_Mpc) # H_0 = 100*h*km/(s*Mpc)
# see also https://github.com/gamer-project/gamer/wiki/Runtime-Parameters%3A-Units#units-in-cosmological-simulations
GAMER_UNIT_L        = GAMER_HUBBLE0**(-1)*Const_Mpc
GAMER_UNIT_T        = GAMER_H0**(-1)
GAMER_UNIT_D        = 3*GAMER_OMEGA_M0*(GAMER_H0**2)/(8*math.pi*Const_NewtonG)


# load the input file
with h5py.File(input_GADGET_file, "r") as f:

    # print all root level object names/keys
    print("Keys: %s" % f.keys())
    input_key_list = list(f.keys())

    # print all attributes associated with "Header"
    Header_list = f["Header"]
    Header_list_attribute_list = f["Header"].attrs.keys()
    input_a_scale_factor  = Header_list.attrs["Time"]
    print(Header_list_attribute_list)

    # output simulation parameters adopted in the GADGET input file
    with open('Input_Parameters.txt', 'wb') as output:

        output.write("GADGET Input Parameters\n\n")

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

    print('Input_Parameters.txt complete')

    if (input_GADGET_Units_Option == 0):
        # GADGET default system units (see User guide: https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf)
        Gadget_HUBBLE0      = float(Header_list.attrs["HubbleParam"])  # dimensionless Hubble parameter "h"
        Gadget_UnitMass     = (1e10*Const_Msun/Gadget_HUBBLE0)         # see "UnitMass in g" on p.26 of the User guide
        Gadget_UnitLength   = (Const_Mpc/Gadget_HUBBLE0)               # see "BoxSize" on p.18 of the User guide
        Gadget_UnitVelocity = (Const_km/Const_s)                       # see "Particle velocities" on p.32 of the User guide

    # "GADGET" to "GAMER" position and velocity unit conversions
    mass_conversion     = Gadget_UnitMass/(GAMER_UNIT_D*GAMER_UNIT_L**3)
    pos_conversion      = Gadget_UnitLength/GAMER_UNIT_L
    vel_conversion      = Gadget_UnitVelocity/(GAMER_UNIT_L/GAMER_UNIT_T)
    print("GADGET-to-GAMER Mass Conversion: %.6f"%mass_conversion)
    print("GADGET-to-GAMER Mass Conversion: %.6f"%pos_conversion)
    print("GADGET-to-GAMER Mass Conversion: %.6f"%vel_conversion)

    # start sorting the particle attributes and initial conditions
    # (1) ordinary point masses should have non-zero entry in the "MassTable"
    GADGET_MassTable = list(Header_list.attrs["MassTable"])
    # (2) "NumPart_ThisFile" should be identical to "NumPart_Total" IFF "input_GADGET_file" is NOT split over several smaller files
    GADGET_NumPart_ThisFile = list(Header_list.attrs["NumPart_ThisFile"])
    # (3) record the number of particles for each type
    GADGET_NumPart_Total = list(Header_list.attrs["NumPart_Total"])
    GADGET_NumPart_Total_HighWord = list(Header_list.attrs["NumPart_Total_HighWord"])
    GADGET_BoxSize = Header_list.attrs["BoxSize"]


    with open('PAR_IC', 'wb') as output:

        for i, NumPart_i in enumerate(GADGET_NumPart_Total):
            # if for "PartType(i)" the total particle number is greater than zero, we import the information to PAR_IC
            if NumPart_i > 0:
                key_NOW = "PartType"+str(i)
                PartType_Coordinates = np.array(f[key_NOW]["Coordinates"], float)
                Velocities_Coordinates = np.array(f[key_NOW]["Velocities"], float)

                mass = np.full(NumPart_i, GADGET_MassTable[i])*mass_conversion

                # set the origin of Cartesian coordinates to the box left bottom corner
                posx = PartType_Coordinates[:,0]*pos_conversion
                posy = PartType_Coordinates[:,1]*pos_conversion
                posz = PartType_Coordinates[:,2]*pos_conversion

                velx = Velocities_Coordinates[:,0]*vel_conversion*(input_a_scale_factor**1.5)
                vely = Velocities_Coordinates[:,1]*vel_conversion*(input_a_scale_factor**1.5)
                velz = Velocities_Coordinates[:,2]*vel_conversion*(input_a_scale_factor**1.5)

                # For GADGET, see "Table 3" on p.31 of https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf
                # For GAMER, we fix the imported particle type to "= 1 (PTYPE_GENERIC_MASSIVE)." See https://github.com/gamer-project/gamer/wiki/Initial-Conditions#IC-File-Particles
                partype = np.full(NumPart_i, 1)

                output.write(mass.astype(dtype=conversion_float_percision).tobytes())
                output.write(posx.astype(dtype=conversion_float_percision).tobytes())
                output.write(posy.astype(dtype=conversion_float_percision).tobytes())
                output.write(posz.astype(dtype=conversion_float_percision).tobytes())
                output.write(velx.astype(dtype=conversion_float_percision).tobytes())
                output.write(vely.astype(dtype=conversion_float_percision).tobytes())
                output.write(velz.astype(dtype=conversion_float_percision).tobytes())
                output.write(partype.astype(dtype=conversion_float_percision).tobytes())

                print('number of particles = %i'%NumPart_i)

        output.close()

    print('PAR_IC complete')