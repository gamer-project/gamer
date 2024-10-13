import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pylab as plt


N_res = 5
N_base = 32

cases = []
cases.append(["HighT", "HighT, T = 1.e+10"])
cases.append(["LowT",  "LowT,  T = 1.e-10"])

fig, ax = plt.subplots(1,1, figsize=(10, 10))
order_2 = False
for case in cases:
    all_res = []
    all_err = []
    folder  = case[0]
    label   = case[1]
    for i in range(N_res):
        N = N_base * 2**i
        read_data = False

        try:
            data = np.loadtxt("./%s/res_%03d/Record__L1Err"%(folder, N))
            E = data[-1][2]

            all_res.append(N)
            all_err.append(E)
            read_data = True
        except:
            pass
            #print("No gamer data under %s/res_%04d/"%(folder, N))

        if read_data: continue

        print("Data missing!(%s/res_%03d/)"%(folder, N))

    if len(all_err) < 2: continue

    coeff = np.polyfit( np.log(all_res), np.log(all_err), 1 )
    print(case, coeff[0])

    if "ATHENA" in label:
        ax.scatter(all_res, all_err, s=40, label="%s, $N^{%.2f}$"%(label, coeff[0]))
    else:
        ax.scatter(all_res, all_err, s=20, label="%s, $N^{%.2f}$"%(label, coeff[0]))
    #ax.plot(all_res, (np.exp(coeff[1]) * all_res**coeff[0]), '--', label="$N^{%.2f}$"%coeff[0])

    coeff[0] = -2.
    if not order_2:
        ax.plot(all_res, ( all_err[0] * all_res[0]**(-coeff[0])  * all_res**coeff[0]), 'k--', label="$N^{-2.00}$")
        order_2 = True

ax.set(xlabel="Number of cells (N)")
ax.set(ylabel="L1 error")
ax.set(xscale='log', yscale='log')
ax.set(xlim=[10, 1000])
ax.set(ylim=[1.e-11, 1.e-7])

ax.legend()
plt.tight_layout()
plt.savefig("SRHD_L1error.png")
plt.close()
#plt.show()
