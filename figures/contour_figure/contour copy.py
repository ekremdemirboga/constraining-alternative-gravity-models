
from likelihood_funcs import contourplot,arrange_sgb
from plot_functions2 import dataReader2,dataReader,dataReader3,interpolated,interpolated2,pp_mass_conv,r_scale,dataReaderDEF,rho_scale
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
import pandas as pd
import numpy as np

dataframe1 = pd.read_excel('lambda14_beta150.xlsx',index_col=0,converters={'mass':float,'R_surface':float}) 
df1 = pd.DataFrame(dataframe1, columns=['mass'])
df2 = pd.DataFrame(dataframe1, columns=['R_surface'])
Mass_sGB_14 =np.asarray(df1.values.tolist()).T[:][0]
Radius_sGB_14 = np.asarray(df2.values.tolist()).T[:][0]
Mass_sGB_14,Radius_sGB_14,Mass_sGB_14_unstable,Radius_sGB_14_unstable = arrange_sgb(Mass_sGB_14,Radius_sGB_14)


dataframe1 = pd.read_excel('lambda22_beta150.xlsx',index_col=0,converters={'mass':float,'R_surface':float}) 
df3 = pd.DataFrame(dataframe1, columns=['mass'])
df4 = pd.DataFrame(dataframe1, columns=['R_surface'])
Mass_sGB_22 =np.asarray(df3.values.tolist()).T[:][0]
print(len(Mass_sGB_22))
Radius_sGB_22 = np.asarray(df4.values.tolist()).T[:][0]
Mass_sGB_22,Radius_sGB_22,Mass_sGB_22_unstable,Radius_sGB_22_unstable = arrange_sgb(Mass_sGB_22,Radius_sGB_22)

print(len(Mass_sGB_22))
print(len(Mass_sGB_22_unstable))

rcParams['font.size'] = 16
rcParams['font.family'] = 'serif'
rcParams['text.color'] = 'black'
rcParams['axes.labelcolor'] = 'black'

datafile = "MRprob_1745.dat"
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
##################################################################################################
plt.subplot(2, 2, 1)
##################################################################################################

plt.xlim([7.5,12])
plt.ylim([0.0,2.3])
Radius,Mass,Prob2D,cs = contourplot(datafile,savedir,type,beta_st)
plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3.5,dashes=(1, 0.7))
plt.plot (radiusGR_stable[:728]*r_scale,mass_ADMGR_stable[:728]/pp_mass_conv, 'k', label = 'GR',linewidth=3.5)

plt.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv, 'red',label = r"Linear DEF",linewidth=3.5)
plt.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'red',linewidth=3.5,linestyle='--',dashes=(1, 0.7))
plt.plot(radius_stable2_Mendes*r_scale*A_r_stable2, mass_ADM_stable2_Mendes/pp_mass_conv, 'red',linewidth=3.5,linestyle='--',dashes=(1, 0.7))

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(visible=False)
plt.yticks(visible=False)
# plt.xlabel (r"R (km)",fontsize=14)
# plt.ylabel(r"$M / M_{\odot}$",fontsize=14)
plt.grid(True)

ax = plt.gca()
plt.text(0.045, 0.15, r"$\beta = 92$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)

axins = zoomed_inset_axes(ax, 5, loc='center left')
axins.set_xlim(10.1, 10.7)
axins.set_ylim(2.05, 2.15)

plt.xticks(visible=False)

plt.yticks(visible=False)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
axins.plot(radius_stable1_Mendes*r_scale*A_r_stable1, mass_ADM_stable1_Mendes/pp_mass_conv,'red',linewidth=3.5)
axins.plot(radius_Mendes*r_scale*A_r, mass_ADM_Mendes/pp_mass_conv, 'red',linewidth=3.5,linestyle='--',dashes=(1, 0.7))
axins.plot(radiusGR_stable[:728]*r_scale,mass_ADMGR_stable[:728]/pp_mass_conv, 'k',linewidth=3.5)
axins.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3.5,dashes=(1, 0.7))
axins.contourf(Radius,Mass, Prob2D.T,levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7],cmap ='Greys')
axins.grid(True)

##################################################################################################
plt.subplot(2, 2, 2)
##################################################################################################
plt.xlim([7.5,12])
plt.ylim([0.0,2.3])
Radius,Mass,Prob2D,cs = contourplot(datafile,savedir,type,beta_st)
plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3.5,dashes=(1, 0.7))
plt.plot (radiusGR_stable[:231]*r_scale,mass_ADMGR_stable[:231]/pp_mass_conv, 'k', label = 'GR',linewidth=3.5)

plt.plot(fR_stt_MENDES_neg*r_scale,m_stt_Mendes_neg/pp_mass_conv,'red',label=r"Linear ",linewidth=3.5)
plt.plot(fR_stt_Mendes_neg2*r_scale,m_stt_Mendes_neg2/pp_mass_conv,'red',linewidth=3.5,linestyle='--',dashes=(1, 0.7))

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xticks(visible=False)
plt.yticks(visible=False)
ax = plt.gca()
plt.text(0.045, 0.15, r"$\beta = -8$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)
plt.grid(True)



##################################################################################################
plt.subplot(2, 2, 3)
##################################################################################################

plt.xlim([7.5,12])
plt.ylim([0.0,2.3])
Radius,Mass,Prob2D,cs = contourplot(datafile,savedir,type,beta_st)
plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3.5,dashes=(1, 0.7))
plt.plot (radiusGR_stable[:281]*r_scale,mass_ADMGR_stable[:281]/pp_mass_conv, 'k', label = 'GR',linewidth=3.5)
print(radiusGR_stable[281]*r_scale)
print(mass_ADMGR_stable[281])

plt.plot(Radius_sGB_22,Mass_sGB_22/pp_mass_conv,'limegreen',label = 'sGB',linewidth=3.5)
plt.plot(Radius_sGB_22_unstable,Mass_sGB_22_unstable/pp_mass_conv,'limegreen',linewidth=4,linestyle='--',dashes=(1, 0.4))


ax = plt.gca()
plt.text(0.045, 0.15, r"$\beta_{GB} = 484$ $\sigma = 150$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)
plt.xlabel (r"R (km)",fontsize=16)
plt.ylabel(r"$M / M_{\odot}$",fontsize=16)
plt.grid(True)
##################################################################################################
plt.subplot(2, 2, 4)
##################################################################################################
plt.xlim([7.5,12])
plt.ylim([0.0,2.3])
Radius,Mass,Prob2D,cs = contourplot(datafile,savedir,type,beta_st)
plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3.5,dashes=(1, 0.7))
plt.plot (radiusGR_stable[:390]*r_scale,mass_ADMGR_stable[:390]/pp_mass_conv, 'k',linewidth=3.5)

plt.plot(Radius_sGB_14,Mass_sGB_14/pp_mass_conv,'limegreen', label = r'sGB',linewidth=3.5)
plt.plot(Radius_sGB_14_unstable,Mass_sGB_14_unstable/pp_mass_conv,'limegreen',linewidth=3.5,linestyle='--',dashes=(1, 0.7))
plt.plot([],[],'k', label = r'GR',linewidth=3.5)
plt.xticks(visible=False)
plt.yticks(visible=False)
plt.grid(True)
plt.plot([],[],'r',label=r'Linear DEF',linewidth=3.5)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=3)
ax = plt.gca()
plt.text(0.045, 0.15, r"$\beta_{GB} = 196$ $\sigma = 150$", verticalalignment='top', fontsize = 14, 
         bbox = dict(facecolor = 'gray', alpha = 0.2),transform=ax.transAxes)
# plt.colorbar(cs)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.90, 0.15, 0.05, 0.7])
fig.colorbar(cs, cax=cbar_ax)
plt.show()