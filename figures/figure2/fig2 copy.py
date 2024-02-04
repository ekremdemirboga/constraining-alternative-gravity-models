import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
import numpy as np
from matplotlib import rcParams
rcParams['font.size'] = 16
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'grey'
rcParams['axes.labelcolor'] = 'black'


def readtxt(eostxt):
    with open(eostxt) as f:
        beta = []
        p = []
        for line in f:
            this_line = line.split()
            if 'beta' in this_line:
                continue
            beta.append(float(this_line[0])*2)
            p.append(float(this_line[1]))
    # beta = np.array(beta)
    # p = np.array(p)
    return beta,p

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/figures/figure2/likelihood_quadratic/"

eos3_positive=prefix+ r"likelihood_quadratic_eos3.txt"
eos4_positive=prefix+ r"likelihood_quadratic_eos4.txt"
eos5_positive=prefix+ r"likelihood_quadratic_eos5.txt"
betas = np.arange(0,89*2,2)

beta3_positive,p3_positive = readtxt(eos3_positive)

beta3 = beta3_positive
p3 = p3_positive
p_eos_3 = np.interp(betas,beta3,p3)
# cs3 = CubicSpline(beta3, p3)
norm3 =  simpson(p_eos_3, betas)


beta4_positive,p4_positive = readtxt(eos4_positive)


beta4 = beta4_positive
p4 = p4_positive
p_eos_4 = np.interp(betas,beta4,p4)
# cs4 = CubicSpline(beta4, p4)
norm4 =  simpson(p_eos_4, betas)


beta5_positive,p5_positive = readtxt(eos5_positive)


beta5 = beta5_positive
p5 = p5_positive
p_eos_5 = np.interp(betas,beta5,p5)
# cs5 = CubicSpline(beta5, p5)
norm5 =  simpson(p_eos_5, betas)

plt.figure()
plt.grid()
fig = plt.gcf()
fig.set_size_inches(8, 8)


rcParams['font.size'] = 17
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'black'
rcParams['axes.labelcolor'] = 'black'
# rcParams['ytick.minor.size'] = 0
# rcParams['ytick.minor.width'] = 0

p_eos_3 =np.array(p_eos_3)
p_eos_4 =np.array(p_eos_4)
p_eos_5 =np.array(p_eos_5)
NORM = norm3+norm4+norm5

plt.subplot(1, 2, 1)

plt.step (beta3,p3/NORM,color ='tab:blue', label = r'$\lambda = $HB',linewidth=3.5)
plt.step (beta4,p4/NORM,color ='tab:orange', label = r'$\lambda = $B',linewidth=3.5)
plt.step (beta5,p5/NORM,color = 'tab:green', label =r'$\lambda = $2B',linewidth=3.5)

ax = plt.gca()
# ax.set_yscale('log')
# every_nth =3
# for n, a in enumerate(ax.yaxis.get_ticklabels()):
#     if n % every_nth != 0:
#         a.set_visible(False)

# rcParams['ytick.minor.size'] = 0
# rcParams['ytick.minor.width'] = 0



p_marg = (p_eos_3/NORM+p_eos_4/NORM+p_eos_5/NORM)
norm_marg = simpson(p_marg, betas)
print(norm_marg)
plt.step(betas, p_marg,color ='tab:red', label = r"marginal", linewidth =2.4)
plt.axhline(y = 1/(89*2), color = 'k', linestyle = '--',label = r"prior",linewidth=3.5*0.6)
plt.xticks(fontsize=15)
plt.yticks([1e-3,3e-3,5e-3,7e-3,9e-3],fontsize=15)
plt.ylabel(r"$P(\beta)$",fontsize=22)

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.grid()

plt.xlim([-5,89*2])
plt.ylim([1e-3,1e-2])

plt.xlabel (r"$\beta$",fontsize=22)
# plt.legend(loc='lower right')
plt.grid(True)


plt.subplot(1, 2, 2)


NORM = norm3
plt.step (beta3,p3/NORM,color ='tab:blue', label = r'$\lambda = $HB',linewidth=3.5)
plt.step ([],[],color ='tab:orange', label = r'$\lambda = $B',linewidth=3.5)
# plt.step ([],[],color = 'tab:green', label =r'$\lambda = $2B',linewidth=3.5)

ax = plt.gca()
# ax.set_yscale('log')

p_marg = (p3/NORM)
norm_marg = simpson(p_marg, beta3)
print(norm_marg)
plt.step(beta3, p_marg,color ='tab:red', label = r"marginal", linewidth =2.4)
plt.axhline(y = 1/(89*2), color = 'k', linestyle = '--',label = r"prior",linewidth=3.5*0.6)
plt.xticks(fontsize=15)
plt.yticks([1e-3,3e-3,5e-3,7e-3,9e-3],fontsize=15)
plt.xlabel (r"$\beta$",fontsize=22)
# plt.ylabel(r"$P(\beta)$",fontsize=20)
plt.grid(True)
plt.xlim([-5,89*2])
plt.ylim([1e-3,1e-2])
plt.yticks(visible=False)
frame1 = plt.gca()
frame1.axes.yaxis.set_ticklabels([])

handles, labels = ax.get_legend_handles_labels()
plt.legend( loc = 'lower right', fancybox=True, ncol=1, labelspacing=0.)
# plt.figlegend(handles, labels, loc = 'lower center', fancybox=True, shadow=True, ncol=6, labelspacing=0.)
plt.show()
# plt.savefig(r"./asd.png", bbox_inches ="tight",dpi=300)