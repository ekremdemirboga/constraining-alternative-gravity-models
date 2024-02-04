import matplotlib.pyplot as plt
from plot_functions2 import dataReader2,dataReader,dataReader3,interpolated,interpolated2,pp_mass_conv,r_scale,dataReaderDEF,rho_scale
from matplotlib import rcParams
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

rcParams['font.size'] = 16
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'black'
rcParams['axes.labelcolor'] = 'black'

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/figures/"
GR = r"gr/eos_3_beta_0000p000_mass_0000p000"
DEF =r"DEF/eos_3_beta_0004p000_mass_0000p000/eos_3_beta_0004p000_mass_0000p000"
DEF_positive = r"DEF/positive/eos_3_beta_0046p000_mass_0000p000/eos_3_beta_0046p000_mass_0000p000"
MENDES_negative = r"Mendes/negative/eos_3_beta_0004p000_mass_0000p000/eos_3_beta_0004p000_mass_0000p000"
MENDES_postive = r"Mendes/positive/eos_3_beta_0046p000_mass_0000p000/eos_3_beta_0046p000_mass_0000p000"

mass_BaryonGR, mass_ADMGR, _, radiusGR, _, rho_cGR, mass_BaryonGR_stable, mass_ADMGR_stable, _, \
radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
= dataReader2(GR)

instability,_, mass_ADM_DEF, A_r_DEF, radius_DEF,_, rho_cDEF, _, _, _, _ = dataReader(DEF)


rhostt1_DEF, m_stt1_DEF, rhostt2_DEF, m_stt2_DEF = interpolated2(rho_cDEF,mass_ADM_DEF)
fR_stt_DEF,m_stt_DEF,fR_stt_DEF2,m_stt_DEF2  =interpolated(radius_DEF,mass_ADM_DEF,A_r_DEF,radiusGR_stable,mass_ADMGR_stable)

mass_ADM_Mendes, radius_Mendes, mass_ADM_stable1_Mendes,radius_stable1_Mendes,mass_ADM_stable2_Mendes,radius_stable2_Mendes,A_r,A_r_stable1,A_r_stable2,rho,rho_stable1,rho_stable2 = dataReader3(MENDES_postive)
mass_ADM_Linear, radius_Linear, mass_ADM_stable1_Linear,radius_stable1_Linear,mass_ADM_stable2_Linear,radius_stable2_Linear,A_r_Linear,A_r_stable1_Linear,A_r_stable2_Linear,rho_Linear,rho_stable1_Linear,rho_stable2_Linear = dataReader3(DEF_positive)



_,_, mass_ADM_Mendes_neg, A_r_Mendes_neg, radius_Mendes_neg,_, rho_Mendes_neg, _, _, _, _ = dataReader(MENDES_negative)
fR_stt_MENDES_neg,m_stt_Mendes_neg,fR_stt_Mendes_neg2,m_stt_Mendes_neg2  =interpolated(radius_Mendes_neg,mass_ADM_Mendes_neg,A_r_Mendes_neg,radiusGR_stable,mass_ADMGR_stable)
rhostt_Mendes, m_stt_Mendes, rhostt_Mendes2, m_stt_Mendes2 = interpolated2(rho_Mendes_neg,mass_ADM_Mendes_neg)
plt.figure()
fig = plt.gcf()
fig.set_size_inches(9,9)

plt.subplot(2, 2, 1)
plt.xlim([7.5,15])
plt.ylim([0,3.5])

plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
plt.plot (radiusGR_stable[:734]*r_scale,mass_ADMGR_stable[:734]/pp_mass_conv, 'k', label = 'GR',linewidth=3)

plt.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv, 'limegreen',label = r"Linear DEF",linewidth=3)
plt.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'limegreen',linewidth=3,linestyle='--')
plt.plot(radius_stable2_Mendes*r_scale*A_r_stable2, mass_ADM_stable2_Mendes/pp_mass_conv, 'limegreen',linewidth=3,linestyle='--')


plt.plot(radius_stable1_Linear*r_scale*A_r_stable1_Linear, mass_ADM_stable1_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
plt.plot(radius_Linear*r_scale*A_r_Linear, mass_ADM_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
plt.plot(radius_stable2_Linear*r_scale*A_r_stable2_Linear, mass_ADM_stable2_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
plt.grid()

plt.xticks(visible=False)
plt.yticks(visible=False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.xlabel (r"R (km)",fontsize=14)
# plt.ylabel(r"$M / M_{\odot}$",fontsize=14)

ax = plt.gca()
plt.text(0.045, 0.95, r"$\beta = 92$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)



axins = zoomed_inset_axes(ax, 5, loc='upper right')
axins.set_xlim(10.1, 10.7)
axins.set_ylim(2.00, 2.17)

plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
axins.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv,'limegreen',linewidth=3)
axins.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'limegreen',linewidth=3,linestyle='--')
axins.plot(radiusGR_stable[:734]*r_scale,mass_ADMGR_stable[:734]/pp_mass_conv, 'k',linewidth=3)
axins.plot(radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
axins.plot(radius_stable1_Linear*r_scale*A_r_stable1_Linear, mass_ADM_stable1_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
axins.grid(True)



plt.subplot(2, 2, 2)
plt.ylim([0,3])
plt.xlim([0,3])
plt.plot (rho_cGR*rho_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
plt.plot (rho_cGR_stable*rho_scale,mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR',linewidth=3)

plt.plot(rho_stable1*rho_scale, mass_ADM_stable1_Mendes/pp_mass_conv, 'limegreen',label = r"Linear DEF",linewidth=3)
plt.plot(rho*rho_scale, mass_ADM_Mendes/pp_mass_conv, 'limegreen',linewidth=3,linestyle='--')
plt.plot(rho_stable2*rho_scale, mass_ADM_stable2_Mendes/pp_mass_conv, 'limegreen',linewidth=3)


plt.plot([],[],'firebrick',linewidth=3,label = r"Quadratic DEF")
plt.plot(rho_stable1_Linear*rho_scale, mass_ADM_stable1_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
plt.plot(rho_Linear*rho_scale, mass_ADM_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
plt.plot(rho_stable2_Linear*rho_scale, mass_ADM_stable2_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')

plt.grid(True)
plt.xticks(visible=False)
plt.yticks(visible=False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.xlabel (r"$\rho_c$ $(g/cm^3)$",fontsize=14)


ax = plt.gca()

plt.text(0.045, 0.95, r"$\beta = 92$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)

axins = zoomed_inset_axes(ax, 3, loc='lower right')
axins.set_xlim(1.4, 1.9)
axins.set_ylim(2.00, 2.17)

plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
axins.plot(rho_stable1*rho_scale, mass_ADM_stable1_Mendes/pp_mass_conv, 'limegreen',label = r"Linear DEF",linewidth=3)
axins.plot(rho*rho_scale, mass_ADM_Mendes/pp_mass_conv, 'limegreen',linewidth=3,linestyle='--')
axins.plot(rho_cGR_stable*rho_scale,mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR',linewidth=3)
axins.plot(rho_stable1_Linear*rho_scale, mass_ADM_stable1_Linear/pp_mass_conv, 'firebrick',linewidth=3,linestyle='--')
axins.grid(True)


plt.subplot(2, 2, 3)
plt.xlim([7.5,15])
plt.ylim([0,3.5])
plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
plt.plot (radiusGR_stable*r_scale,mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR',linewidth=3)

plt.plot(fR_stt_MENDES_neg*r_scale,m_stt_Mendes_neg/pp_mass_conv,'limegreen',label=r"Linear ",linewidth=3)
plt.plot(fR_stt_Mendes_neg2*r_scale,m_stt_Mendes_neg2/pp_mass_conv,'limegreen',linewidth=3,linestyle='--')
plt.plot(fR_stt_Mendes_neg2[71:]*r_scale,m_stt_Mendes_neg2[71:]/pp_mass_conv,'limegreen',linewidth=3,linestyle='--')


plt.plot(fR_stt_DEF*r_scale,m_stt_DEF/pp_mass_conv,'firebrick',linewidth=3,label=r"Quadratic ")
plt.plot(fR_stt_DEF2*r_scale,m_stt_DEF2/pp_mass_conv,'firebrick',linewidth=3,linestyle='--')
ax = plt.gca()
plt.text(0.045, 0.95, r"$\beta = -8$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.xlabel (r"R (km)",fontsize=20)
plt.ylabel(r"$M / M_{\odot}$",fontsize=20)
plt.grid(True)

plt.subplot(2, 2, 4)
plt.xlim([0,3])
plt.ylim([0,3])
plt.plot (rho_cGR*rho_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
plt.plot (rho_cGR_stable*rho_scale,mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR',linewidth=3)

plt.plot(rhostt_Mendes,m_stt_Mendes,'limegreen',label=r"Linear",linewidth=3)
plt.plot(rhostt_Mendes2,m_stt_Mendes2,'limegreen',linewidth=3,linestyle='--')

plt.plot([],[],'firebrick',linewidth=3,label = r"Quadratic DEF")
plt.plot(rhostt1_DEF,m_stt1_DEF,'firebrick',linewidth=3)
plt.plot(rhostt2_DEF,m_stt2_DEF,'firebrick',linewidth=3,linestyle="--")

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.xticks(visible=False)
plt.yticks(visible=False)
plt.xlabel (r"$\rho_c$ $(g/cm^3)$",fontsize=20)
# plt.ylabel(r"$M / M_{\odot}$",fontsize=14)
plt.grid(True)

ax = plt.gca()
plt.text(0.045, 0.95, r"$\beta = -8$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center',ncol=3)
plt.show()






############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

# plt.figure()
# plt.grid()
# fig = plt.gcf()
# fig.set_size_inches(10,7)
# # fig.tight_layout()
# plt.xlim([7.5,18])

# plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--', label = 'unstable',linewidth=3)
# plt.plot (radiusGR_stable*r_scale,mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR',linewidth=3)

# plt.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv, 'b',label = r"Linear $ \beta=92$",linewidth=3)
# plt.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'b--',linewidth=3)
# plt.plot(radius_stable2_Mendes*r_scale*A_r_stable2, mass_ADM_stable2_Mendes/pp_mass_conv, 'b',linewidth=3)

# plt.plot(radius_stable1_Linear*r_scale*A_r_stable1_Linear, mass_ADM_stable1_Linear/pp_mass_conv, 'r',label = r"Quadratic $ \beta=92$",linewidth=3)
# plt.plot(radius_Linear*r_scale*A_r_Linear, mass_ADM_Linear/pp_mass_conv, 'r--',linewidth=3)
# plt.plot(radius_stable2_Linear*r_scale*A_r_stable2_Linear, mass_ADM_stable2_Linear/pp_mass_conv, 'r',linewidth=3)


# plt.legend()

# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.xlabel (r"R (km)",fontsize=14)
# plt.ylabel(r"$M / M_{\odot}$",fontsize=14)


# ax = plt.gca()

# axins = zoomed_inset_axes(ax, 6, loc='right')
# axins.set_xlim(10.1, 10.9)
# axins.set_ylim(2.02, 2.15)
# axins.grid(True)

# plt.xticks(visible=False)

# plt.yticks(visible=False)
# mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
# axins.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv,'b',linewidth=3)
# axins.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'b--',linewidth=3)
# axins.plot(radiusGR_stable*r_scale,mass_ADMGR_stable/pp_mass_conv, 'k',linewidth=3)
# axins.plot(radius_stable1_Linear*r_scale*A_r_stable1_Linear, mass_ADM_stable1_Linear/pp_mass_conv, 'r',linewidth=3)
# plt.show()