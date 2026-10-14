#!/usr/bin/env python3

import numpy as np

NEWTON_G  = 1.0000000e+00
Total_M   = 2.3062183e-02
Max_R     = 1.5000000e+00

Ave_Rho   = 3.0*Total_M/(4.0*np.pi*Max_R**3)
t_dyn     = ( NEWTON_G*Ave_Rho )**(-0.5)

print( f'{NEWTON_G =: .8e}' )
print( f'{Total_M  =: .8e}' )
print( f'{Max_R    =: .8e}' )
print( f'{Ave_Rho  =: .8e}' )
print( f'{t_dyn    =: .8e}' )
