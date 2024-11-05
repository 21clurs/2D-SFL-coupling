# quick and dirty

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

outdir = "./sim/out/oscillation_tests/"
plots = os.listdir(outdir)
start = 0
end = len(plots)

periods = np.empty((len(plots),2))

for i in range(start,end):
    n,_ = plots[i].split(".")
    t = []
    amplitudes = []  # x-positions
    dt = 0
    with open(outdir+plots[i],"r") as curr_plot_file:
        for row in curr_plot_file:
            r = row.split(" ")
            t.append(float(r[0]))
            amplitudes.append(float(r[1]))
    t_subset = t[300:]
    amplitudes_subset = amplitudes[300:]
    periods[i,0] = int(n)
    periods[i,1] = t_subset[amplitudes_subset.index(max(amplitudes[300:]))]
    #periods.append([int(n), max(amplitudes[300:]), t_subset[amplitudes_subset.index(max(amplitudes[300:]))]])

periods = periods[periods[:,0].argsort()]
print(periods)

errs = periods[:,1]-math.pi
log_errs = np.log(errs)
log_n = np.log(periods[:,0])

plt.plot(log_n,log_errs)
plt.scatter(log_n,log_errs)

plt.title(r"log-log plot of $T_{error}$ vs. mesh resolution $N$")
plt.xlabel(r"log $|N|$")
plt.ylabel(r"log $|T-\pi|$")
#plt.grid(True)

plt.gca().spines['top'].set_visible(False) 
plt.gca().spines['right'].set_visible(False) 

#plt.show()
plt.savefig('./plots/loglog_Terr_N.png', format="png")

    
