import numpy as np
import re



#====================================================================================================
# Function
#====================================================================================================
def LoadValueFromInputFile(filename, key_list):
    regex_num = r"\s*([-+]?\d+\.?\d*[eE]?[-+]?\d*)"
    param = dict()
    with open(filename, "r") as f:
        param_in = f.read()

        for key in key_list:
            value = re.findall(key + regex_num, param_in)

            # assume the value is a float
            param[key] = float(value[0])
    return param


def AnalyticalDens(ParX, ParY, Center, Dens_Bg, BoxSize):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the density at these particles' locations
    Radius = np.hypot(ParX - Center[0], ParY - Center[1])

    return Dens_Bg * (1.0 + 5.0 * Radius / BoxSize)


def AnalyticalPres(ParX, ParY, Center, Pres_Bg, BoxSize):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the pressure at these particles' locations
    Radius = np.hypot(ParX - Center[0], ParY - Center[1])

    return Pres_Bg * (1.0 + 5.0 * Radius / BoxSize)


def AnalyticalVelX(ParX, ParY, Center, Ang_Freq):
    # use the formulae in SetGridIC() from the ParticleTest test problem
    # to compute the velocity in the x direction at these particles' locations
    Radius    = np.hypot(ParX - Center[0], ParY - Center[1])
    Sin_Theta = (ParY - Center[1]) / Radius
    Velocity  = Ang_Freq * Radius

    return -1.0 * Velocity * Sin_Theta
