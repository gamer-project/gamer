import numpy as np
import random as rd


BOX_SIZE_X = 1.0
BOX_SIZE_Y = 1.0
BOX_SIZE_Z = 1.0


Theta_Min = 0.0
Theta_Max = np.pi

Phi_Min   = 0.0
Phi_Max   = 2.0*np.pi

Temp_Min = 1.0
Temp_Max = 10.0*Temp_Min

Number_Blast_Waves = 64

print("#%19s%20s%20s%20s%20s%20s" %  ("x" ,"y" ,"z" ,"theta" ,"phi" ,"temperature" ))

for i in range(Number_Blast_Waves):
  x     = rd.uniform(0.0, BOX_SIZE_X)
  y     = rd.uniform(0.0, BOX_SIZE_Y)
  z     = rd.uniform(0.0, BOX_SIZE_Z)
  
  theta = rd.uniform(Theta_Min, Theta_Max)
  phi   = rd.uniform(Phi_Min,     Phi_Max)
  
  Temp  = rd.uniform(Temp_Min, Temp_Max)
  
  
  print("%20.12e%20.12e%20.12e%20.12e%20.12e%20.12e" % (x, y, z, theta, phi, Temp))
