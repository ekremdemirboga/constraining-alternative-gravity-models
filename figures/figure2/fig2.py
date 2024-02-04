import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
import numpy as np
from matplotlib import rcParams
rcParams['font.size'] = 14
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'grey'
rcParams['axes.labelcolor'] = 'black'


def readtxt(eostxt):
    print("asd",eostxt)
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

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/figures/figure2/likelihood_linear/"
eos1_negative=prefix+ r"likelihood_negative_eos1.txt"
eos1_positive=prefix+ r"likelihood_positive_eos1.txt"
eos2_negative=prefix+ r"likelihood_negative_eos2.txt"
eos2_positive=prefix+ r"likelihood_positive_eos2.txt"
eos3_negative=prefix+ r"likelihood_negative_eos3.txt"
eos3_positive=prefix+ r"likelihood_positive_eos3.txt"
eos4_negative=prefix+ r"likelihood_negative_eos4.txt"
eos4_positive=prefix+ r"likelihood_positive_eos4.txt"
eos5_negative=prefix+ r"likelihood_negative_eos5.txt"
eos5_positive=prefix+ r"likelihood_positive_eos5.txt"
betas = np.arange(-202,212,1)

beta1_negative,p1_negative = readtxt(eos1_negative)
beta1_positive,p1_positive = readtxt(eos1_positive)

beta1_negative.reverse()
p1_negative.reverse()
beta1 = beta1_negative+beta1_positive
p1 = p1_negative+p1_positive
norm1 =  simpson(p1, beta1)

beta2_negative,p2_negative = readtxt(eos2_negative)
beta2_positive,p2_positive = readtxt(eos2_positive)

beta2_negative.reverse()
p2_negative.reverse()
beta2 = beta2_negative+beta2_positive
p2 = p2_negative+p2_positive
norm2 =  simpson(p2, beta2)


beta3_negative,p3_negative = readtxt(eos3_negative)
beta3_positive,p3_positive = readtxt(eos3_positive)

beta3_negative.reverse()
p3_negative.reverse()
beta3 = beta3_negative+beta3_positive
p3 = p3_negative+p3_positive
# p_eos_3 = np.interp(betas,beta3,p3)
# cs3 = CubicSpline(beta3, p3)
norm3 =  simpson(p3, beta3)


beta4_negative,p4_negative = readtxt(eos4_negative)
beta4_positive,p4_positive = readtxt(eos4_positive)

beta4_negative.reverse()
p4_negative.reverse()
beta4 = beta4_negative+beta4_positive
p4 = p4_negative+p4_positive
# p_eos_4 = np.interp(betas,beta4,p4)
# cs4 = CubicSpline(beta4, p4)
norm4 =  simpson(p4, beta4)


beta5_negative,p5_negative = readtxt(eos5_negative)
beta5_positive,p5_positive = readtxt(eos5_positive)

beta5_negative.reverse()
p5_negative.reverse()
beta5 = beta5_negative+beta5_positive
p5 = p5_negative+p5_positive
# p_eos_5 = np.interp(betas,beta5,p5)
# cs5 = CubicSpline(beta5, p5)
norm5 =  simpson(p5, beta5)

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
p1 =np.array(p1)
p2 =np.array(p2)
p3 =np.array(p3)
p4 =np.array(p4)
p5 =np.array(p5)
NORM = norm1+norm2+norm3+norm4+norm5
print(NORM)

plt.subplot(1, 2, 1)
# plt.semilogy (beta1,p1/NORM,color ='tab:brown', label = r'$\lambda = $2H',linewidth=3.5)
plt.step (beta2,p2/NORM,color ='tab:purple', label = r'$\lambda = $H',linewidth=3.5)
plt.step (beta3,p3/NORM,color ='tab:blue', label = r'$\lambda = $HB',linewidth=3.5)
plt.step (beta4,p4/NORM,color ='tab:orange', label = r'$\lambda = $B',linewidth=3.5)
plt.step (beta5,p5/NORM,color = 'tab:green', label =r'$\lambda = $2B',linewidth=3.5)

ax = plt.gca()
ax.set_yscale('log')
every_nth =2
for n, a in enumerate(ax.yaxis.get_ticklabels()):
    if n % every_nth != 0:
        a.set_visible(False)

# rcParams['ytick.minor.size'] = 0
# rcParams['ytick.minor.width'] = 0

p1 =np.array(p1)
p2 =np.array(p2)
p3 =np.array(p3)
p4 =np.array(p4)
p5 =np.array(p5)
p_marg = (p1/NORM+p2/NORM+p3/NORM+p4/NORM+p5/NORM)
norm_marg = simpson(p_marg, beta3)
print(norm_marg)
plt.semilogy(beta3, p_marg,color ='tab:red', label = r"marginal", linewidth =3.4)
plt.axhline(y = 1/414, color = 'k', linestyle = '--',label = r"prior",linewidth=3.5*0.6)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel(r"$P(\beta)$",fontsize=17)
plt.grid()

plt.xlim([-202,212])
plt.ylim([1e-5,1e-2])

plt.xlabel (r"$\beta$",fontsize=17)
# plt.legend(loc='lower right')
plt.grid(True)


plt.subplot(1, 2, 2)
eos4_positive=prefix+ r"likelihood_positive_eos4_2ndprior.txt"
beta4_positive,p4_positive = readtxt(eos4_positive)
p4 = p4_negative+p4_positive
norm4 =  simpson(p4, beta4)
p4 =np.array(p4)
NORM = norm1+norm2+norm3+norm4
print(NORM)

# plt.step (beta1,p1/NORM,color ='tab:brown', label = r'$\lambda = $2H',linewidth=3.5)
plt.step (beta2,p2/NORM,color ='tab:purple', label = r'$\lambda = $H',linewidth=3.5)
plt.step (beta3,p3/NORM,color ='tab:blue', label = r'$\lambda = $HB',linewidth=3.5)
plt.step (beta4,p4/NORM,color ='tab:orange', label = r'$\lambda = $B',linewidth=3.5)
plt.step ([],[],color = 'tab:green', label =r'$\lambda = $2B',linewidth=3.5)

ax = plt.gca()
ax.set_yscale('log')

p_marg = (p1/NORM+p2/NORM+p3/NORM+p4/NORM)
norm_marg = simpson(p_marg, beta3)
print(norm_marg)
plt.semilogy(beta3, p_marg,color ='tab:red', label = r"marginal", linewidth =3.4)
plt.axhline(y = 1/414, color = 'k', linestyle = '--',label = r"prior",linewidth=3.5*0.6)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel (r"$\beta$",fontsize=17)
# plt.ylabel(r"$P(\beta)$",fontsize=20)
plt.grid(True)
plt.xlim([-202,212])
plt.ylim([1e-5,1e-2])
plt.yticks(visible=False)

handles, labels = ax.get_legend_handles_labels()
plt.legend( loc = 'lower right', fancybox=True, ncol=1, labelspacing=0.,fontsize = 14)

# plt.figlegend(handles, labels, loc = 'lower center', fancybox=True, shadow=True, ncol=6, labelspacing=0.)
plt.show()
# plt.savefig(r"./linearPosterior.png", bbox_inches ="tight",dpi=300)