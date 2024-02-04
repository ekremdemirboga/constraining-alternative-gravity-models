

import numpy as np
import matplotlib.pyplot as plt

#conversion constants
pp_mass_conv =0.01185185
Mass_Sun = 1.989e30
G = 6.674e-11
c = 299792458
Goverc2 = 7.42591549e-28
r_scale = (Goverc2*Mass_Sun/1000)/pp_mass_conv
beta_conv =2
rho0 = 1.66*10**14 #g/cm^3
rho_scale = ((Mass_Sun*1000)/(pp_mass_conv*r_scale**3*10**15))/rho0

def nameArrange(directory):
    ## takes directory as string and gives the name of the datafile. directory = ..../datafile
    datafile = ""
    for i in directory[::-1]:
        if (i == '/'):
            datafile = datafile[::-1]
            break
        datafile += i
    
    return datafile

def Repeated(Mass_array,Mass_value):
    for i in range(len(Mass_array)):
        if (abs(Mass_value - Mass_array[i])<1e-12):
            print("repeat!!!")
            return True
    return False


def dataReader(directory):
    
    #Takes the input as "only" the name of the file and reads the data files 
    #with the names format: eos_x_beta_xxxxpxxx_mass_xxxxpxxx
    
    #Order: isScalarized | massIncreased | isMonotoneInside | isMonotoneOutside | 
    #mass_Baryon | mass_ADM | AST(radius) | radius  |  phi_c | rho_c | ctr |
    
    #Also seperates the stable solutions
    datafile=nameArrange(directory)

    isScalarized = []
    massIncreased = []
    isMonotoneInside = []
    isMonotoneOutside = []
    
    mass_Baryon = []
    mass_Baryon_stable = []
    mass_ADM = []
    mass_ADM_stable = []
    A_r = []
    A_r_stable = []
    radius = []
    radius_stable = []
    phi_c = []
    phi_c_stable =[]
    rho_c = []
    rho_c_stable = []
    ctr = []
    beta_st = float(datafile[11:15]) #+ float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29]) #+ float("0."+ datafile[30:33])
    eos = datafile[4]
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))

        
    
    with open(directory) as f:

        
        for line in f:
            this_line = line.split()

            if 'isScalarized' in this_line:
                continue
            if ((this_line) ==["\n"] or this_line ==[]):
                print("EMPTY")
                return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, \
                A_r_stable, radius_stable,phi_c_stable, rho_c_stable, ctr, beta_st, mass_st, eos

            # if (Repeated(mass_ADM,float(this_line[5]))):
            #     continue

            isScalarized.append(int(this_line[0]))
            massIncreased.append(int(this_line[1]))
            isMonotoneInside.append(int(this_line[2]))
            isMonotoneOutside.append(int(this_line[3]))

            ## categorizing according to instabilities
            if ((int(this_line[1]) == 1)): # if stable 
                mass_Baryon_stable.append(float(this_line[4]))
                mass_ADM_stable.append(float(this_line[5]))
                A_r_stable.append(float(this_line[6]))
                radius_stable.append(float(this_line[7]))
                phi_c_stable.append(abs(float(this_line[8])))
                rho_c_stable.append(float(this_line[9]))
                ctr.append(float(this_line[10]))
            
            mass_Baryon.append(float(this_line[4]))
            mass_ADM.append(float(this_line[5]))
            A_r.append(float(this_line[6]))
            radius.append(float(this_line[7]))
            phi_c.append(abs(float(this_line[8])))
            rho_c.append(float(this_line[9]))
            ctr.append(float(this_line[10]))
            
            

        # turning to arrays
        isScalarized = np.asarray(isScalarized)
        massIncreased = np.asarray(massIncreased)
        isMonotoneInside = np.asarray(isMonotoneInside)
        isMonotoneOutside = np.asarray(isMonotoneOutside)
        
        mass_Baryon = np.asarray(mass_Baryon)
        mass_ADM = np.asarray(mass_ADM)
        A_r = np.asarray(A_r)
        radius = np.asarray(radius)
        phi_c = abs(np.asarray(phi_c))
        rho_c = np.asarray(rho_c)
        mass_Baryon_stable = np.asarray(mass_Baryon_stable)
        mass_ADM_stable = np.asarray(mass_ADM_stable)
        A_r_stable = np.asarray(A_r_stable)
        radius_stable = np.asarray(radius_stable)
        phi_c_stable = abs(np.asarray(phi_c_stable))
        rho_c_stable = np.asarray(rho_c_stable)
        ctr = np.asarray(ctr)
        
        #fix the first star
        # if (mass_ADM[1]>mass_ADM[0]):
        #     massIncreased[0] = 1
        # else:
        #     massIncreased[0] = 0

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, \
            A_r_stable, radius_stable,phi_c_stable, rho_c_stable, ctr, beta_st, mass_st, eos


def dataReader_negative(directory):
    
    #Takes the input as "only" the name of the file and reads the data files 
    #with the names format: eos_x_beta_xxxxpxxx_mass_xxxxpxxx
    
    #Order: isScalarized | massIncreased | isMonotoneInside | isMonotoneOutside | 
    #mass_Baryon | mass_ADM | AST(radius) | radius  |  phi_c | rho_c | ctr |
    
    #Also seperates the stable solutions
    datafile=nameArrange(directory)

    isScalarized = []
    massIncreased = []
    isMonotoneInside = []
    isMonotoneOutside = []
    
    mass_Baryon = []
    mass_Baryon_stable = []
    mass_ADM = []
    mass_ADM_stable = []
    A_r = []
    A_r_stable = []
    radius = []
    radius_stable = []
    phi_c = []
    phi_c_stable =[]
    rho_c = []
    rho_c_stable = []
    ctr = []

    beta_st = float(datafile[11:15]) #+ float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29]) #+ float("0."+ datafile[30:33])
    eos = datafile[4]
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))

        
    
    with open(directory) as f:

        stable =0
        MADM_prev=0
        first_line_flag =0
        for line in f:
            this_line = line.split()

            if 'isScalarized' in this_line:
                continue
            if ((this_line) ==["\n"] or this_line ==[]):
                print("EMPTY")
                return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, \
                A_r_stable, radius_stable,phi_c_stable, rho_c_stable, ctr, beta_st, mass_st, eos

            # if (Repeated(mass_ADM,float(this_line[5]))):
            #     continue

            if (float(this_line[5])>MADM_prev) and first_line_flag:
                stable =1
                first_line_flag = 1
            else:
                stable =0
                first_line_flag =1
            
            isScalarized.append(int(this_line[0]))
            massIncreased.append(int(this_line[1]))
            isMonotoneInside.append(int(this_line[2]))
            isMonotoneOutside.append(int(this_line[3]))

            ## categorizing according to instabilities
            if (stable): # if stable 
                mass_Baryon_stable.append(float(this_line[4]))
                mass_ADM_stable.append(float(this_line[5]))
                A_r_stable.append(float(this_line[6]))
                radius_stable.append(float(this_line[7]))
                phi_c_stable.append(abs(float(this_line[8])))
                rho_c_stable.append(float(this_line[9]))
                ctr.append(float(this_line[10]))
                
            
            mass_Baryon.append(float(this_line[4]))
            mass_ADM.append(float(this_line[5]))
            print(float(this_line[5]))
            A_r.append(float(this_line[6]))
            radius.append(float(this_line[7]))
            phi_c.append(abs(float(this_line[8])))
            rho_c.append(float(this_line[9]))
            ctr.append(float(this_line[10]))
            MADM_prev = float(this_line[5])
            

        # turning to arrays
        isScalarized = np.asarray(isScalarized)
        massIncreased = np.asarray(massIncreased)
        isMonotoneInside = np.asarray(isMonotoneInside)
        isMonotoneOutside = np.asarray(isMonotoneOutside)
        
        mass_Baryon = np.asarray(mass_Baryon)
        mass_ADM = np.asarray(mass_ADM)
        A_r = np.asarray(A_r)
        radius = np.asarray(radius)
        phi_c = abs(np.asarray(phi_c))
        rho_c = np.asarray(rho_c)
        mass_Baryon_stable = np.asarray(mass_Baryon_stable)
        mass_ADM_stable = np.asarray(mass_ADM_stable)
        A_r_stable = np.asarray(A_r_stable)
        radius_stable = np.asarray(radius_stable)
        phi_c_stable = abs(np.asarray(phi_c_stable))
        rho_c_stable = np.asarray(rho_c_stable)
        ctr = np.asarray(ctr)
        
        #fix the first star
        # if (mass_ADM[1]>mass_ADM[0]):
        #     massIncreased[0] = 1
        # else:
        #     massIncreased[0] = 0

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, \
            A_r_stable, radius_stable,phi_c_stable, rho_c_stable, ctr, beta_st, mass_st, eos

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/"
GR_dir_eos3 = prefix + r"GR/eos3/eos_3_beta_0000p000_mass_0000p000"
def plot(directory,savedir, directoryGR = GR_dir_eos3, Compactness=0, GR_on=1, Madm_R = 1, Madm_rhoc=1, Madm_phic=0, phic_rhoc =0, R_rhoc=0, R_phic=0, Binding=0, Mbaryon_rhoc = 0):

    ### plots desired figures.

    print("plotting")
    mass_BaryonGR, mass_ADMGR, _, radiusGR, _, rho_cGR, mass_BaryonGR_stable, mass_ADMGR_stable, _, \
    radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
    = dataReader(directoryGR)


    if (dataReader(directory) != 0):
        mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, A_r_stable, \
        radius_stable,phi_c_stable, rho_c_stable, _, beta_st, mass_st, eos \
        = dataReader(directory)
    else:
        return 0
    

    datafile=nameArrange(directory)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.serif": ["Palatino"],
        "xtick.labelsize": 15,
        "ytick.labelsize": 15
    })
#     plt.xlim([8,radius[-1]*A_r[-1]*r_scale+1])
    
    
    if (Madm_R==1):
        maxR = [9.47651,11.08725,9.691275,9.90604,10.7651,10.65772,11.125,9.125,11.375,9.875,11.875,9.625,10.436,10.832]
        maxMass = [1.592617,1.884564,1.631544,1.67047,1.826174,1.787248,1.1625,0.9375,1.2875,1.0125,0.9875,1.0875,0.5,1.46]
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        plt.xlim([6,55])

        if (GR_on==1):
            plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--', label = 'GR unstable')
            plt.plot (radiusGR_stable*r_scale,mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR stable')
        plt.xlabel ("Radius (km)",fontsize=20)
        plt.ylabel(r"$M_{ADM} / M_{sun}$",fontsize=20)
        title = (r"$M_{ADM}$ vs R \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        plt.plot (radius*A_r*r_scale,mass_ADM/pp_mass_conv, 'g.', label = 'STT unstable')
        plt.plot (radius_stable*A_r_stable*r_scale,mass_ADM_stable/pp_mass_conv, 'r.', label = 'STT stable')

        # plt.scatter(maxR,maxMass)
        plt.legend()
        plt.savefig(savedir+"Madm_R_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()

    if (Compactness==1):
        
        maxR = [9.47651,11.08725,9.691275,9.90604,10.7651,10.65772,11.125,9.125,11.375,9.875,11.875,9.625,10.436,10.832]
        maxMass = [1.592617,1.884564,1.631544,1.67047,1.826174,1.787248,1.1625,0.9375,1.2875,1.0125,0.9875,1.0875,0.5,1.46]
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        # plt.xlim([8,15])

        if (GR_on==1):
            plt.plot ((mass_ADMGR/pp_mass_conv)/(radiusGR*r_scale),(mass_ADMGR/pp_mass_conv),'k--', label = 'GR unstable')
            plt.plot ((mass_ADMGR_stable/pp_mass_conv)/(radiusGR_stable*r_scale),mass_ADMGR_stable/pp_mass_conv, 'k', label = 'GR stable')
        plt.xlabel ("M/R",fontsize=20)
        plt.ylabel(r"$M_{ADM} / M_{sun}$",fontsize=20)
        title = (r"$M/R$ vs M \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        plt.plot ((mass_ADM/pp_mass_conv)/(radius*A_r*r_scale),mass_ADM/pp_mass_conv, 'g.', label = 'STT unstable')
        plt.plot ((mass_ADM_stable/pp_mass_conv)/(radius_stable*A_r_stable*r_scale),mass_ADM_stable/pp_mass_conv, 'ro', label = 'STT stable')

        plt.scatter(np.asarray(maxMass)/np.asarray(maxR),np.asarray(maxMass))
        plt.legend()
        plt.savefig(savedir+"Madm_R_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()

    if (phic_rhoc==1):

        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        # plt.xlim([8,15])

        plt.xlabel(r"$\rho_c / \rho_0$",fontsize=20)
        plt.ylabel (r"$\phi_c$",fontsize=20)
        title = (r"$\rho_c$ vs $\phi_c$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        plt.plot (rho_c*rho_scale,phi_c, 'g.', label = 'STT unstable')
        plt.plot (rho_c_stable*rho_scale, phi_c_stable, 'ro', label = 'STT stable')

        plt.legend()
        plt.savefig(savedir+"phic_rhoc_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()
    if (Madm_rhoc==1):
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        if (GR_on==1):
            plt.plot (rho_cGR*rho_scale,mass_ADMGR/pp_mass_conv,'k--')
            plt.plot (rho_cGR_stable*rho_scale,mass_ADMGR_stable/pp_mass_conv, 'k')
        plt.xlabel (r"$\rho_c / \rho_0$",fontsize=20)
        plt.ylabel(r"$M_{ADM} / M_{sun}$",fontsize=20)
        title = (r"$M_{ADM}$ vs $\rho_c$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        plt.plot (rho_c*rho_scale,mass_ADM/pp_mass_conv, '-', label = 'beta = '+str(beta_st*beta_conv),linewidth = 4)
        plt.plot (rho_c_stable*rho_scale,mass_ADM_stable/pp_mass_conv, 'ro', label = 'STT stable')

        plt.legend()
        plt.savefig(savedir+"Madm_rhoc_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()
    if (Madm_phic==1):
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        plt.xlabel (r"$\phi_c$",fontsize=20)
        plt.ylabel(r"$M_{ADM} / M_{sun}$",fontsize=20)
        title = (r"$M_{ADM}$ vs $\phi_c$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        
        plt.plot (phi_c,mass_ADM/pp_mass_conv, 'g.', label = 'STT unstable')
        plt.plot (phi_c_stable,mass_ADM_stable/pp_mass_conv, 'ro', label = 'STT stable')

        plt.legend()
        plt.savefig(savedir+"Madm_phic_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()
    if (R_rhoc==1):

        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        plt.xlabel(r"$\rho_c / \rho_0$",fontsize=20)
        plt.ylabel ("Radius (km)",fontsize=20)
        if (GR_on==1):
            plt.plot (rho_cGR*rho_scale, radiusGR*r_scale,'k--', label = 'GR unstable')
            plt.plot (rho_cGR_stable*rho_scale, radiusGR_stable*r_scale, 'k', label = 'GR stable')

        title = (r"R vs $\rho_c$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        
        plt.plot (rho_c*rho_scale, radius*A_r*r_scale, 'g.', label = 'STT unstable')
        plt.plot (rho_c_stable*rho_scale, radius_stable*A_r_stable*r_scale, 'ro', label = 'STT stable')

        plt.legend()
        plt.savefig(savedir+"R_rhoc_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()
    if (R_phic==1):

        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        plt.xlabel(r"$\phi_c$",fontsize=20)
        plt.ylabel ("Radius (km)",fontsize=20)

        title = (r"R vs $\phi_c$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        
        plt.plot (phi_c, radius*A_r*r_scale, 'g.', label = 'STT unstable')
        plt.plot (phi_c_stable, radius_stable*A_r_stable*r_scale, 'ro', label = 'STT stable')

        
        plt.legend()
        plt.savefig(savedir+"R_phic_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()

    if (Binding==1):

        bindingEnergy = np.subtract (mass_Baryon,mass_ADM)
        bindingEnergy_stable = np.subtract (mass_Baryon_stable,mass_ADM_stable)
        bindingEnergyGR = np.subtract (mass_BaryonGR, mass_ADMGR)
        bindingEnergyGR_stable = np.subtract (mass_BaryonGR_stable, mass_ADMGR_stable)
        
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        plt.xlabel(r"$M_{baryon}/M_{sun}$",fontsize=20)
        plt.ylabel ("Binding Energy",fontsize=20)
        if (GR_on==1):
            plt.plot (mass_BaryonGR/pp_mass_conv,bindingEnergyGR/pp_mass_conv,'k--', label = 'GR unstable')
            plt.plot (mass_BaryonGR_stable/pp_mass_conv,bindingEnergyGR_stable/pp_mass_conv, 'k', label = 'GR stable')

        title = (r"BindingEnergy($M_{baryon}-M_{ADM}$) vs $M_{baryon}$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        plt.plot (mass_Baryon/pp_mass_conv, bindingEnergy/pp_mass_conv, 'g.', label = 'STT unstable')
        plt.plot (mass_Baryon_stable/pp_mass_conv, bindingEnergy_stable/pp_mass_conv, 'ro', label = 'STT stable')

        plt.legend()
        plt.savefig(savedir+"Binding_Mbaryon_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()
    if (Mbaryon_rhoc==1):
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 9)
        title = (r"$M_{b}$ vs $\rho_c$ \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        if (GR_on==1):
            plt.title(title, fontsize=20)
            plt.plot (rho_cGR*rho_scale,mass_BaryonGR/pp_mass_conv,'k--', label = 'GR unstable')
            plt.plot (rho_cGR_stable*rho_scale,mass_BaryonGR_stable/pp_mass_conv, 'k', label = 'GR stable')
        plt.xlabel (r"$\rho_c / \rho_0$",fontsize=20)
        plt.ylabel(r"$M_{b} / M_{sun}$",fontsize=20)
        
        plt.title(title, fontsize=20)
        plt.plot (rho_c*rho_scale,mass_Baryon/pp_mass_conv, 'g.', label = 'STT unstable')
        plt.plot (rho_c_stable*rho_scale,mass_Baryon_stable/pp_mass_conv, 'r.', label = 'STT stable')
        plt.legend()
        plt.savefig(savedir+"Mbaryon_rhoc_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        plt.close()