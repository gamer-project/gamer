import numpy as np
def Zero_Pot(x):
    return 0

r0   = 0.1
N    = 1000
f = open("external_pot_table.txt", "w")
r = r0*np.logspace(-1,1,N)
Pot = [Zero_Pot(r[i]) for i in range(N)]

f.write("External Potential Table\n")
for i in range(N):
    f.write("{:10.4f} {:10.4f}\n".format(r[i], Pot[i]))
f.close()
