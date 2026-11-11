import matplotlib.pyplot as plt
import numpy as np

data       = np.loadtxt('Record__CoolingErr', skiprows=1)

time  = data[:, 0]
temp  = data[:, 2]
tcool = data[:, 3]

time      /= tcool[0]

plt.figure(figsize=(6, 4.5))
plt.scatter(time, temp, s=7, color='black')
plt.yscale('log')
plt.xlabel('time ($t_{cool}$)')
plt.ylabel('temperature (K)')
plt.title('Exact Cooling')
plt.legend()
plt.savefig('evolution_ExactCooling.png', dpi=300)