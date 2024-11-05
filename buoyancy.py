# quick and dirty

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

outdir = "./sim/out/buoyancy_tests/"
test_files = os.listdir(outdir)
start = 0
end = len(test_files)

vels = np.empty((len(test_files),2))

for i in range(start,end):
    n,_ = test_files[i].split(".")
    with open(outdir+test_files[i],"r") as curr_plot_file:
        for row in curr_plot_file:
            r = row.split(" ")
            v = [float(r[0]),float(r[1])]
            vels[i,0] = n
            vels[i,1] = np.linalg.norm(v)


vels = vels[vels[:,0].argsort()]
print(vels)

errs = np.abs(vels[:,1])
log_errs = np.log(errs)
log_n = np.log(vels[:,0])

plt.plot(log_n,log_errs)
plt.scatter(log_n,log_errs)

plt.title(r"log-log plot of $V_{error}$ vs. mesh resolution $N$")
plt.xlabel(r"log $|N|$")
plt.ylabel(r"log $|V|$")

#plt.grid(True)
plt.gca().spines['top'].set_visible(False) 
plt.gca().spines['right'].set_visible(False) 

#plt.show()
plt.savefig('./plots/loglog_Verr_N.png', format="png")

    
