import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('Output__Error', skiprows=1)

time = data[:, 0]
dump_id = data[:, 1]
temp_nume = data[:, 2]
temp_anal = data[:, 3]
err = data[:, 4]
tcool_nume = data[:, 5]
tcool_anal = data[:, 6]

time /= tcool_anal[0]

plt.figure(figsize=(6, 4.5))
plt.plot(time, temp_anal, label='Analytical', linestyle='-', color='black')
plt.scatter(time, temp_nume, label='Numerical', s=7, color='black')
plt.yscale('log')
plt.xlabel(r'time ($t_{cool}$)')
plt.ylabel(r'temperature (K)')
plt.legend()
plt.savefig('evolution_ExactCooling.png', dpi=300)
