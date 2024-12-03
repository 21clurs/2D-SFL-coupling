# quick and dirty

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

if len(sys.argv)>1:
    target_domain_size = int(sys.argv[1])

rb_density_list = [0.3, 0.5, 0.7, 0.9, 1.1, 1.5]
vels_list = []

for i in range(len(rb_density_list)):
    outdir = f"./sim/out/added_mass_tests_again/added_mass_{int(rb_density_list[i]*10)}/"
    test_files = os.listdir(outdir)
    for j in range(len(test_files)):
        domain_size,_ = test_files[j].split(".")
        if(int(domain_size) == target_domain_size):
            with open(outdir+test_files[j],"r") as curr_plot_file:
                for row in curr_plot_file:
                    r = row.split(" ")
                    v = [float(r[0]),float(r[1])]
                    vels_list.append(np.linalg.norm(v))

vels = np.column_stack((np.array(rb_density_list), np.array(vels_list)))
vels = vels[vels[:,0].argsort()]
print(vels)
plt.scatter(vels[:,0],vels[:,1],zorder=1)

x = np.linspace(min(rb_density_list)-.1, max(rb_density_list)+.1,1000)
dt = 0.0001
y = dt * np.abs((1-x)/(1+x))
plt.plot(x,y,color='orange',ls="--",zorder=0)

plt.title(r"$||\mathbf{v}||$ vs. density relative $\alpha$, $\rho_{solid}=\alpha\rho_{fluid}$")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$||\mathbf{v}||$")

#plt.grid(True)
plt.gca().spines['top'].set_visible(False) 
plt.gca().spines['right'].set_visible(False) 

#plt.show()
plt.savefig(f'./plots/V_density_{target_domain_size:02d}.png', format="png")
plt.savefig(f'./plots/V_density_{target_domain_size:02d}.pdf', format="pdf")



