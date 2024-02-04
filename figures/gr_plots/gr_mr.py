from plot_functions2 import dataReader2,rho_scale,r_scale,pp_mass_conv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mycolorpy import colorlist as mcp

rcParams['font.size'] = 16
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'grey'
rcParams['axes.labelcolor'] = 'black'

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/figures/gr_plots/GR/"
eos1= prefix+ r"eos1/eos_1_beta_0000p000_mass_0000p000"
eos2= prefix+ r"eos2/eos_2_beta_0000p000_mass_0000p000"
eos3= prefix+ r"eos3/eos_3_beta_0000p000_mass_0000p000"
eos4= prefix+ r"eos4/eos_4_beta_0000p000_mass_0000p000"
eos5= prefix+ r"eos5/eos_5_beta_0000p000_mass_0000p000"
plt.figure()
fig = plt.gcf()
fig.set_size_inches(8, 15)
plt.subplot(2, 1, 1)
plt.xlim([7,17])
plt.ylim([0,3])
EOSs = [eos1, eos2, eos3, eos4, eos5]
name = [r'$\lambda=$2H',r'$\lambda=$H',r'$\lambda=$HB', r'$\lambda=$B', r'$\lambda=$2B' ]
colors = ['tab:brown','tab:purple','tab:blue','tab:orange', 'tab:green']
# colors = ['dimgray','teal','forestgreen','tab:blue', 'darkblue']

# colors=mcp.gen_color(cmap="GnBu",n=7)
print(colors)
for i in range(5):
    _, mass_ADMGR, _, radiusGR, _, rho_cGR, _, mass_ADMGR_stable, _, \
    radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
    = dataReader2(EOSs[i])

    plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,color=colors[i],linestyle="--",linewidth=3.5)
    plt.plot (radiusGR_stable*r_scale,mass_ADMGR_stable/pp_mass_conv, color=colors[i], label = name[i],linewidth=3.5)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel (r"R (km)",fontsize=14)
plt.ylabel(r"$M / M_{\odot}$",fontsize=14)
plt.grid(True)


plt.subplot(2, 1, 2)
plt.xlim([0,4e15])

for i in range(5):
    _, mass_ADMGR, _, radiusGR, _, rho_cGR, _, mass_ADMGR_stable, _, \
    radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
    = dataReader2(EOSs[i])

    plt.plot (rho_cGR*rho_scale,mass_ADMGR/pp_mass_conv,color=colors[i],linestyle="--",linewidth=3.5)
    plt.plot (rho_cGR_stable*rho_scale,mass_ADMGR_stable/pp_mass_conv, color=colors[i], label = name[i],linewidth=3.5)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel (r"$\rho_c$ $(g/cm^3)$",fontsize=14)
plt.ylabel(r"$M / M_{\odot}$",fontsize=14)
plt.grid(True)
plt.legend()
plt.show()