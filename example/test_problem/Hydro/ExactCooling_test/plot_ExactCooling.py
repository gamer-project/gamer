import matplotlib.pyplot as plt
import numpy as np

data       = np.loadtxt('Record__CoolingErr', skiprows=1)

time       = data[:, 3]
temp_nume  = data[:, 5]
temp_anal  = data[:, 6]
tcool_anal = data[:, 9]

time      /= tcool_anal[0]

plt.figure(figsize=(6, 4.5))
plt.plot(   time, temp_anal, label='Analytical', linestyle='-', color='black')
plt.scatter(time, temp_nume, label='Numerical' , s=7          , color='black')
plt.yscale('log')
plt.xlabel('time ($t_{cool}$)')
plt.ylabel('temperature (K)')
plt.title('Exact Cooling')
plt.legend()
plt.savefig('evolution_ExactCooling.png', dpi=300)