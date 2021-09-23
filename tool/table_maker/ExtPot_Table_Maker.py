import numpy as np
def Zero_Pot(x):
    return 0
def Plummer(x):
    rho_0 = 1.0
    r0 = 0.1
    Frac = 1
    G = 1.0
    GM = G*(4.0/3.0)*np.pi*(r0**3)*rho_0*Frac
    print( -GM/(1+x**2)**0.5)
    return -GM/(1+x**2)**0.5

r0   = 0.1
N    = 1000
f = open("external_pot_table.txt", "w")
r = r0*np.logspace(-1,1,N)
Pot = [Plummer(r[i]/r0) for i in range(N)]

f.write("External Potential Table\n")
for i in range(N):
    f.write("{:10.4f} {:10.7f}\n".format(r[i], Pot[i]))
f.close()
