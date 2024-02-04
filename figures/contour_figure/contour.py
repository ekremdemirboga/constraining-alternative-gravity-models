
from likelihood_funcs import contourplot
from plot_functions2 import dataReader2,dataReader,dataReader3,interpolated,interpolated2,pp_mass_conv,r_scale,dataReaderDEF,rho_scale
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

rcParams['font.size'] = 16
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'grey'
rcParams['axes.labelcolor'] = 'black'

datafile = "MRprob_1724.dat"
prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/figures/"
savedir= r"./"
type=1
beta_st = -8
GR =prefix+ r"gr/eos_3_beta_0000p000_mass_0000p000"
MENDES_negative = prefix+r"Mendes/negative/eos_3_beta_0004p000_mass_0000p000/eos_3_beta_0004p000_mass_0000p000"
MENDES_postive = prefix + r"Mendes/positive/eos_3_beta_0046p000_mass_0000p000/eos_3_beta_0046p000_mass_0000p000"

mass_BaryonGR, mass_ADMGR, _, radiusGR, _, rho_cGR, mass_BaryonGR_stable, mass_ADMGR_stable, _, \
radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
= dataReader2(GR)
_,_, mass_ADM_Mendes_neg, A_r_Mendes_neg, radius_Mendes_neg,_, rho_Mendes_neg, _, _, _, _ = dataReader(MENDES_negative)
fR_stt_MENDES_neg,m_stt_Mendes_neg,fR_stt_Mendes_neg2,m_stt_Mendes_neg2  =interpolated(radius_Mendes_neg,mass_ADM_Mendes_neg,A_r_Mendes_neg,radiusGR_stable,mass_ADMGR_stable)
mass_ADM_Mendes, radius_Mendes, mass_ADM_stable1_Mendes,radius_stable1_Mendes,mass_ADM_stable2_Mendes,radius_stable2_Mendes,A_r,A_r_stable1,A_r_stable2,rho,rho_stable1,rho_stable2 = dataReader3(MENDES_postive)



plt.figure()
fig = plt.gcf()
fig.set_size_inches(10,9)
plt.xlim([6,15])
plt.ylim([0.6,2.4])

Radius,Mass,Prob2D = contourplot(datafile,savedir,type,beta_st)
plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=2.5)
plt.plot (radiusGR_stable[:725]*r_scale,mass_ADMGR_stable[:725]/pp_mass_conv, 'k', label = 'GR',linewidth=2.5)


# plt.plot(fR_stt_MENDES_neg*r_scale,m_stt_Mendes_neg/pp_mass_conv,'red',label=r"Linear ",linewidth=2.5)
# plt.plot(fR_stt_Mendes_neg2*r_scale,m_stt_Mendes_neg2/pp_mass_conv,'red',linewidth=2.5,linestyle='--')


plt.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv, 'red',label = r"Linear DEF",linewidth=2.5)
plt.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'red',linewidth=2.5,linestyle='--')
plt.plot(radius_stable2_Mendes*r_scale*A_r_stable2, mass_ADM_stable2_Mendes/pp_mass_conv, 'red',linewidth=2.5,linestyle='--')




plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel (r"R (km)",fontsize=14)
plt.ylabel(r"$M / M_{\odot}$",fontsize=14)
plt.legend()
ax = plt.gca()

axins = zoomed_inset_axes(ax, 5, loc='lower left')
axins.set_xlim(10.1, 10.7)
axins.set_ylim(2.05, 2.15)

plt.xticks(visible=False)

plt.yticks(visible=False)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
axins.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv,'red',linewidth=2.5)
axins.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'red',linewidth=2.5,linestyle='--')
axins.plot(radiusGR_stable[:728]*r_scale,mass_ADMGR_stable[:728]/pp_mass_conv, 'k',linewidth=2.5)
axins.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=2.5)

axins.contourf(Radius,Mass, Prob2D.T,levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7],cmap ='Greys')

axins.grid(True)




plt.show()