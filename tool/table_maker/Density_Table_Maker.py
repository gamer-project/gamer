import numpy as np
def NFW(x):
    return 1/(x*(1+x)*(1+x))
def Plummer(x):
    return (1+x*x)**(-2.5)

rho0 = 1.0
r0   = 0.1
N    = 1000
f = open("density_table.txt", "w")
r = r0*np.logspace(-1,1,N)
dens = [NFW(r[i]) for i in range(N)]

f.write("Density Talbe\n")
for i in range(N):
    f.write("{:10.4f} {:10.4f}\n".format(r[i], dens[i]))
f.close()
