# quick and dirty

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

if len(sys.argv)>1:
    rb_density = float(sys.argv[1])

outdir = f"./sim/out/added_mass_tests/added_mass_{int(rb_density*10)}/"
test_files = os.listdir(outdir)
start = 0
end = len(test_files)

domain_size_list = []
vels_list = []

for i in range(start,end):
    domain_size,_ = test_files[i].split(".")
    if(int(domain_size)>=10):
        with open(outdir+test_files[i],"r") as curr_plot_file:
            for row in curr_plot_file:
                r = row.split(" ")
                v = [float(r[0]),float(r[1])]
                domain_size_list.append(int(domain_size))
                vels_list.append(np.linalg.norm(v))
                #vels[i,0] = domain_size
                #vels[i,1] = np.linalg.norm(v)
vels = np.column_stack((np.array(domain_size_list), np.array(vels_list)))
vels = vels[vels[:,0].argsort()]
print(vels)

domain = vels[:,0]
#domain = np.square(domain)
#domain = 4*domain
#domain = np.pi * np.square(domain/2)

velocity_data =  vels[:,1]

dt = 0.0001
analytic = dt*np.abs((1-rb_density)/(1+rb_density))
plt.axhline(y = analytic, color = 'orange', ls = '--',zorder=0) 

errs = np.abs(vels[:,1] - analytic)

print(domain)
print(errs)

log_errs = np.log(errs)
log_domain = np.log(domain)

#find line of best fit
a, b = np.polyfit(log_domain, log_errs, 1)
x = np.linspace(min(log_domain),max(log_domain))
#plt.plot(x, a*x+b, color="orange", zorder=0)
print(a)

#plot first order convergence
#plt.plot(x, -x+b + ((a*log_domain[0]+b)-(-log_domain[0]+b)), color="orange", ls="--", zorder=0)

#plt.plot(log_n,log_errs)
plt.scatter(domain, velocity_data, zorder=1)
#plt.scatter(log_domain,log_errs, zorder=1)


#plt.title(r"log-log plot of $V_{error}$ vs. domain size $d$, $\rho_{solid}=$"+str(rb_density)+r"$\rho_{liquid}$")
#plt.xlabel(r"log $|d|$")
#plt.ylabel(r"log $|V_{error}|$")
plt.title(r"Plot of $||\mathbf{v}||$ vs. domain size $d$, $\rho_{solid}=$"+str(rb_density)+r"$\rho_{liquid}$")
plt.xlabel(r"$d$")
plt.ylabel(r"$||\mathbf{v}||$")

#plt.grid(True)
plt.gca().spines['top'].set_visible(False) 
plt.gca().spines['right'].set_visible(False) 

#plt.show()
plt.savefig(f'./plots/V_domain_{int(rb_density*10):02d}.png', format="png")