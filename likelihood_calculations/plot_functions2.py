

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep,splev,CubicSpline,interp1d
#conversion constants
pp_mass_conv =0.01185185
Mass_Sun = 1.989e30
G = 6.674e-11
c = 299792458
Goverc2 = 7.42591549e-28
r_scale = (Goverc2*Mass_Sun/1000)/pp_mass_conv
beta_conv =-2
rho0 = 1.66*10**14 #g/cm^3
rho_scale = ((Mass_Sun*1000)/(pp_mass_conv*r_scale**3*10**15))

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
            return True
    return False

def dataReader(directory):

    datafile=nameArrange(directory)

    isScalarized = []
    massIncreased = []
    isMonotoneInside = []
    isMonotoneOutside = []
    mass_Baryon = []
    mass_ADM = []
    A_r = []
    radius = []
    phi_c = []
    rho_c = []
    ctr = []
    beta_st = float(datafile[11:15]) + float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29])+ float("0."+ datafile[30:33])
    eos = datafile[4]

    stable_branch = False
    
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))
    with open(directory) as f:
        mas=0.0
        rad=0.0
        for line in f:
            if len(line.strip()) == 0:
                continue
            this_line = line.split()
            if 'isScalarized' in this_line:
                continue

            # if(abs( (float(this_line[5])) -mas)<1e-12 or abs( float(this_line[7])-rad)<1e-15):
            #     continue

            if (Repeated(mass_ADM,float(this_line[5]))):
                continue

            # if (int(this_line[1]) ==1):
            #     stable_branch = True
            # if (int(this_line[1]) ==1 and stable_branch):
            #     isScalarized.append(int(this_line[0]))
            #     massIncreased.append(int(this_line[1]))
            #     isMonotoneInside.append(int(this_line[2]))
            #     isMonotoneOutside.append(int(this_line[3]))
            #     mass_Baryon.append(float(this_line[4]))
            #     mass_ADM.append(float(this_line[5]))
            #     A_r.append(float(this_line[6]))
            #     radius.append(float(this_line[7]))
            #     phi_c.append(float(this_line[8]))
            #     rho_c.append(float(this_line[9]))
            #     ctr.append(float(this_line[10]))
                
            # if(stable_branch and int(this_line[1]) ==0):
            #     if(beta_st>11 and float(this_line[5])<1.80e-02):
            #         massIncreased.append(int(this_line[1]))
            #         mass_Baryon.append(float(this_line[4]))
            #         mass_ADM.append(float(this_line[5]))
            #         A_r.append(float(this_line[6]))
            #         radius.append(float(this_line[7]))
            #         phi_c.append(float(this_line[8]))
            #         rho_c.append(float(this_line[9]))
            #         ctr.append(float(this_line[10]))

            #         # rho_c, mass_ADM, A_r, radius = sort_mass(rho_c,mass_ADM,A_r, radius)
                    
            #         mass_Baryon = np.asarray(mass_Baryon)
            #         mass_ADM = np.asarray(mass_ADM)
            #         A_r = np.asarray(A_r)
            #         radius = np.asarray(radius)
            #         phi_c = np.asarray(phi_c)
            #         rho_c = np.asarray(rho_c)
            #         ctr = np.asarray(ctr)
            #         return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos
            #     else:
            #         continue


            isScalarized.append(int(this_line[0]))
            massIncreased.append(int(this_line[1]))
            isMonotoneInside.append(int(this_line[2]))
            isMonotoneOutside.append(int(this_line[3]))
            mass_Baryon.append(float(this_line[4]))
            mass_ADM.append(float(this_line[5]))
            A_r.append(float(this_line[6]))
            radius.append(float(this_line[7]))
            phi_c.append(float(this_line[8]))
            rho_c.append(float(this_line[9]))
            ctr.append(float(this_line[10]))  
                  
            rad = float(this_line[7])
            mas = float(this_line[5])

        # turning to arrays

        # rho_c, mass_ADM, A_r, radius = sort_mass(rho_c,mass_ADM,A_r, radius)

        
        mass_Baryon = np.asarray(mass_Baryon)
        mass_ADM = np.asarray(mass_ADM)
        A_r = np.asarray(A_r)
        radius = np.asarray(radius)
        phi_c = np.asarray(phi_c)
        rho_c = np.asarray(rho_c)
        ctr = np.asarray(ctr)
        
        

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos

def sort_mass(rho,mass,A_r, radius):
    new_Mass = [x for _,x in sorted(zip(rho,mass),reverse=True)]
    new_radius = [x for _,x in sorted(zip(rho,radius),reverse=True)]
    new_A = [x for _,x in sorted(zip(rho,A_r),reverse=True)]
    rho = sorted(rho,reverse=True)
    return rho,new_Mass, new_A, new_radius
    

def dataReader2(directory):
    
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

    beta_st = float(datafile[11:15]) + float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29])+ float("0."+ datafile[30:33])
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

def split(R_stt,M_stt):
    min_r  =999
    R_stt_1 = []
    R_stt_2 = []
    M_stt_1 = []
    M_stt_2 = []
    for r in range(len(R_stt)):
        if R_stt[r]< min_r:
            min_r = R_stt[r]
            R_stt_1.append(R_stt[r])
            M_stt_1.append(M_stt[r])
        else:
            R_stt_2.append(R_stt[r])
            M_stt_2.append(M_stt[r])
    return R_stt_1, R_stt_2, M_stt_1, M_stt_2

def interpolated(R_stt,M_stt,A_r,R_gr,M_gr):
    #interpolate GR
    fR_gr = CubicSpline(M_gr,R_gr, bc_type='natural')

    fR_stt,m_stt= Interpolated_stt(R_stt*A_r, M_stt)
    # ##interpolate in two parts,
    # R_stt_1, R_stt_2, M_stt_1, M_stt_2 = split(A_r*R_stt,M_stt)
    # ###1###
    # knots = np.linspace(min(M_stt_1),max(M_stt_1), num=200, endpoint = True)
    # tck, u = splprep([M_stt_1, R_stt_1],s=0,t = knots)
    # new_points = splev(u, tck)
    # m_stt_1 = new_points[0]
    # fR_stt_1 = new_points[1]
    # ###2###
    # knots = np.linspace(min(M_stt_2),max(M_stt_2), num=50, endpoint = True)
    # tck, u = splprep([M_stt_2, R_stt_2],s=0, t = knots)
    # new_points = splev(u, tck)
    # m_stt_2 = new_points[0]
    # fR_stt_2 = new_points[1]


    # knots = np.linspace(min(M_stt),max(M_stt), num=500, endpoint = True)
    # tck, u = splprep([M_stt, R_stt*A_r],s=0, t = knots)
    # new_points = splev(u, tck)
    # m_stt = new_points[0]
    # fR_stt = new_points[1]
    # m_stt, fR_stt = arrange(m_stt,fR_stt)
    print("------------------")
    m_gr1 = np.linspace(min(M_gr),min(M_stt), num=200,endpoint = False)
    m_gr2 = np.linspace(max(M_stt),max(M_gr), num=200)

    return fR_stt,m_stt, fR_gr,m_gr1,m_gr2 


def Interpolated_stt(X,Y):
    # Define some points:
    points = np.array([X.tolist(),
                    Y.tolist()]).T  # a (nbre_points x nbre_dim) array
    # Linear length along the line:
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)/distance[-1]
    # Interpolation for different methods:
    alpha = np.linspace(0, 1, 1000)


    interpolator =  interp1d(distance, points, kind='cubic', axis=0)

    # Graph:
    # plt.figure(figsize=(7,7))
    
    # plt.plot(*interpolator(alpha).T, '-', label='cubic')
    X_values = interpolator(alpha)[:,0]
    Y_values = interpolator(alpha)[:,1]
    # plt.plot(*points.T, 'ok', label='original points')
    # plt.axis('equal'); plt.legend(); plt.xlabel('x'); plt.ylabel('y');
    # plt.show()
    return X_values, Y_values

def arrange(m_stt,fR_stt):
    new_m_stt = []
    new_fR_stt = []
    for i in range(len(m_stt)-1):
        if m_stt[i]<1.8e-02:   
            if (m_stt[i+1]<=m_stt[i]):
                new_m_stt.append(m_stt[i])
                new_fR_stt.append(fR_stt[i])
            else:
                new_m_stt.append(m_stt[i])
                new_fR_stt.append(fR_stt[i])
                return np.asarray(new_m_stt), np.asarray(new_fR_stt)
        else:
            new_m_stt.append(m_stt[i])
            new_fR_stt.append(fR_stt[i])
    return np.asarray(new_m_stt), np.asarray(new_fR_stt)


prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/"
GR_dir_eos3 = prefix + r"GR/eos4/eos_4_beta_0000p000_mass_0000p000"
def plot2(directory,savedir, directoryGR, GR_on=1, Madm_R = 1,rhoc_Madm=0,clr='black',):

    ### plots desired figures.
    mass_BaryonGR, mass_ADMGR, _, radiusGR, _, rho_cGR, mass_BaryonGR_stable, mass_ADMGR_stable, _, \
    radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
    = dataReader2(directoryGR)

    
    mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, _, beta_st, mass_st, eos = dataReader(directory)
    
    datafile=nameArrange(directory)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.serif": ["Palatino"],
        "xtick.labelsize": 15,
        "ytick.labelsize": 15
    })    
    
    #interpolate functions
    # fR_stt_1,m_stt_1,fR_stt_2,m_stt_2, fR_gr,m_gr1,m_gr2  =interpolated(radius,mass_ADM,A_r,radiusGR_stable,mass_ADMGR_stable)
    fR_stt_1,m_stt_1, fR_gr,m_gr1,m_gr2  =interpolated(radius,mass_ADM,A_r,radiusGR_stable,mass_ADMGR_stable)
    
    if (Madm_R==1):
        # maxR = [9.47651,11.08725,9.691275,9.90604,10.7651,10.65772,11.125,9.125,11.375,9.875,11.875,9.625,10.436,10.832]
        # maxMass = [1.592617,1.884564,1.631544,1.67047,1.826174,1.787248,1.1625,0.9375,1.2875,1.0125,0.9875,1.0875,0.5,1.46]
        # plt.figure()
        # plt.grid()
        # fig = plt.gcf()
        # fig.set_size_inches(13, 9)
        # plt.xlim([3,22])
        if (GR_on==1):
            plt.plot (radiusGR*r_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
            plt.plot (radiusGR_stable*r_scale,mass_ADMGR_stable/pp_mass_conv, 'k',linewidth=3)
        # plt.xlabel ("Radius (km)",fontsize=20)
        # plt.ylabel(r"$M_{ADM} / M_{sun}$",fontsize=20)
        # title = (r"$M_{ADM}$ vs R \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        # plt.title(title, fontsize=20)
        # plt.plot(radius*r_scale*A_r, mass_ADM/pp_mass_conv,'r.',linewidth=2)
        if beta_st==1:
            r_stable,r_unstable,mass_stable, mass_unstable = seperate(fR_stt_1,m_stt_1)
            plt.plot ((fR_stt_1*r_scale),m_stt_1/pp_mass_conv,color=clr,linewidth=3,linestyle="--",dashes=(1, 0.3))
            plt.plot ((fR_stt_1[:525]*r_scale),m_stt_1[:525]/pp_mass_conv,color=clr,linewidth=3,linestyle="-")
            plt.plot ((fR_stt_1[910:1000]*r_scale),m_stt_1[910:1000]/pp_mass_conv,color=clr,linewidth=3,linestyle="-")
        else:
            r_stable,r_unstable,mass_stable, mass_unstable = seperate(fR_stt_1,m_stt_1)
            plt.plot ((fR_stt_1*r_scale),m_stt_1/pp_mass_conv,color=clr,linewidth=3,linestyle="--",dashes=(1, 0.3))
            plt.plot ((r_stable*r_scale),mass_stable/pp_mass_conv,'-',color=clr,linewidth=3)
        # plt.plot ((fR_stt_2*r_scale),m_stt_2/pp_mass_conv, 'r--')
        # plt.scatter(maxR,maxMass)
        # plt.legend()
        # plt.savefig(savedir+"Madm_R_"+ datafile + ".png",dpi=300, bbox_inches='tight')
        # plt.close()



    if (rhoc_Madm==1):
        # plt.figure()
        # plt.grid()
        # fig = plt.gcf()
        # fig.set_size_inches(13, 9)
        # plt.xlim([6,15])
        if (GR_on==1):
            plt.plot (rho_cGR*rho_scale,mass_ADMGR/pp_mass_conv,'k--',linewidth=3)
            plt.plot (rho_cGR_stable*rho_scale,mass_ADMGR_stable/pp_mass_conv, 'k',linewidth=3)


        if beta_st==1:
            rho_stable,rho_unstable,mass_stable, mass_unstable = seperate(rho_c,mass_ADM)
            plt.plot(rho_stable[:200]*rho_scale,mass_stable[:200]/pp_mass_conv,'-',color=clr,linewidth=3,label= r"$\beta =$" +str(int(beta_st*2)))
            plt.plot(rho_c*rho_scale, mass_ADM/pp_mass_conv,'--',color=clr,linewidth=3,dashes=(1, 0.3))
            plt.plot(rho_stable[250:]*rho_scale,mass_stable[250:]/pp_mass_conv,'-',color=clr,linewidth=3)

        else:
            rho_stable,rho_unstable,mass_stable, mass_unstable = seperate(rho_c,mass_ADM)
            plt.plot(rho_stable*rho_scale,mass_stable/pp_mass_conv,'-',color=clr,linewidth=3,label= r"$\beta =$" +str(int(beta_st*2)))
            plt.plot(rho_c*rho_scale, mass_ADM/pp_mass_conv,'--',color=clr,linewidth=3,dashes=(1, 0.3))

        # plt.plot ((fR_stt_1*r_scale),m_stt_1/pp_mass_conv, 'r--', label = 'STT')
        # plt.plot ((fR_stt_2*r_scale),m_stt_2/pp_mass_conv, 'r--')
        # plt.scatter(maxR,maxMass)

def seperate(rho,mass):
    
    rho_stable = []
    rho_unstable = []
    mass_stable = []
    mass_unstable = []

    for i in range(len(mass)-1):
        if(mass[i]>mass[i+1]):
            mass_stable.append(mass[i])
            rho_stable.append(rho[i])
        else:
            rho_unstable.append(rho[i])
            mass_unstable.append(mass[i])
    return np.array(rho_stable),np.array(rho_unstable),np.array(mass_stable), np.array(mass_unstable)