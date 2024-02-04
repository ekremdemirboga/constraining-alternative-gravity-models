from cProfile import label
from pprint import pp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep,splev,CubicSpline,RectBivariateSpline,interp1d#,lagrange
from scipy.integrate import simpson,trapezoid

#conservation constants
pp_mass_conv =0.01185185
Mass_Sun = 1.989e30
G = 6.674e-11
c = 299792458
Goverc2 = 7.42591549e-28
r_scale = (Goverc2*Mass_Sun/1000)/pp_mass_conv
beta_conv =-2
rho0 = 1.66*10**14 #g/cm^3
rho_scale = ((Mass_Sun*1000)/(pp_mass_conv*r_scale**3*10**15))/rho0

def MRprobRead(directory): # reads the datafile in directory
    firstColumn = []
    secondColumn = []
    thirdColumn = []
    with open(directory) as f:
    ##reads the neutron stars mass-radius-probability data files in a given directory.
        for i in range(6):
            next(f)
        for line in f:
            this_line = line.split()
            if this_line[0] =="!":
                continue
            firstColumn.append(float(this_line[0]))
            secondColumn.append(float(this_line[1]))
            thirdColumn.append(float(this_line[2]))
    name = (directory)[8:-4]
    return firstColumn, secondColumn, thirdColumn, name 

def convertProb(Prob,n_radius,n_mass,radius,mass): ## makes probability vector into an array for 2d plot
    Prob2d = np.reshape(Prob,(n_radius,n_mass))
    spl= RectBivariateSpline(radius, mass, Prob2d)
    normalize= spl.integral(np.amin(radius), np.amax(radius), np.amin(mass), np.amax(mass))
    Prob2d = Prob2d/ normalize
    return Prob2d,spl,normalize
    
def MassRadius(datafile,type): #returns mass,radius,prob2d for each data file
    firstColumn,secondColumn,thirdColumn,name = MRprobRead("MRprob2"+"/"+datafile)
    print(name)
    if (type==1):
        Radius = np.asarray(np.sort(list(set(secondColumn))))
        n_radius = len(Radius)
        Mass = np.asarray(firstColumn[:n_radius])
        n_mass = len(Mass)
        Prob = np.asarray(thirdColumn)
    elif type==2:
        if name+".dat" in ["MRprob_M13.dat", "MRprob_M30.dat", "MRprob_NGC6304.dat", "MRprob_NGC6397.dat", "MRprob_OmCen.dat"]:
            n_mass= 99
            n_radius= 100
        if name+".dat" in ["MRprob_M28.dat"]:
            n_mass= 100
            n_radius= 100      
        if name+".dat" in ["MRprob_X5.dat" , "MRprob_X7.dat"]:
            n_mass= 51
            n_radius= 51
        Radius = np.asarray(np.sort(list(set(firstColumn))))
        Mass = np.asarray(secondColumn[:n_mass])
        Prob = np.asarray(thirdColumn)
        print(n_mass)
        print(n_radius)
    Prob2D,spl,normalize = convertProb(Prob,n_radius,n_mass,Radius,Mass)
    return Radius,Mass,Prob2D,name,spl,normalize

def contourplot(datafile,savedir,type,beta_st,mass_st=0): ## plots the datafile
    Radius, Mass, Prob2D, name,spl,normalize = MassRadius(datafile,type)
    import matplotlib.colors as mcolors
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
        if n == -1:
            n = cmap.N
        new_cmap = mcolors.LinearSegmentedColormap.from_list(
            'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap
    new_cmap = truncate_colormap(plt.get_cmap("Greys"), minval=0.0, maxval=0.6, n=-1)

    cs = plt.contourf(Radius,Mass, Prob2D.T,levels=[0.1,0.2,0.3,0.4,0.5],cmap =new_cmap)

    # bar = plt.colorbar(cs)
    # plt.savefig(savedir+"MR_diagram_for_"+name+"_beta_"+str(beta_st)+".png", dpi=300, bbox_inches='tight')
    # plt.close()
    return Radius,Mass,Prob2D,cs

    

### below is for relaxation data 


def Repeated(Mass_array,Mass_value):
    for i in range(len(Mass_array)):
        if (abs(Mass_value - Mass_array[i])<1e-10):
            return True
    return False


def nameArrange(directory):
    ## takes directory as string and gives the name of the datafile. directory = ..../datafile
    datafile = ""
    for i in directory[::-1]:
        if (i == '/'):
            datafile = datafile[::-1]
            break
        datafile += i
    
    return datafile

def interpolated(R_stt,M_stt,A_r,R_gr,M_gr,Model):
    #interpolate GR
    fR_gr = CubicSpline(M_gr,R_gr, bc_type='natural')
    fR_stt,m_stt= Interpolated_stt(R_stt*A_r, M_stt)

    if Model == 'linear/positive':
        m_stt,fR_stt = arrange_pos(m_stt,fR_stt)
    else:
        m_stt,fR_stt = arrange_neg(m_stt,fR_stt)
    # print(m_stt)
    print("------------------")
    if min(M_gr)<min(m_stt):
        m_gr1 = np.linspace(min(M_gr),min(m_stt), num=500,endpoint = False)
    else:
        m_gr1 = np.array([])
    
    if Model == 'linear/positive':
        m_gr2 = np.array([])
    else:
        m_gr2 = np.linspace(max(m_stt),max(M_gr), num=500,endpoint = False)

    

    # m_stt = np.flip(m_stt)
    # fR_stt = np.flip(fR_stt)

    return fR_stt,m_stt, fR_gr,m_gr1,m_gr2 


def interpolated2(R_stt,M_stt,A_r,R_gr,M_gr):
    if (M_stt[0]>M_stt[1]):
        M_stt = np.flip(M_stt)
        R_stt = np.flip(R_stt)
        A_r = np.flip(A_r)
    
    m_stt = np.linspace(min(M_stt),max(M_stt), num=1000)

    fR_gr = CubicSpline(M_gr,R_gr, bc_type='natural')

    M_stt = M_stt[0:-2]
    R_stt = R_stt[0:-2]
    A_r = A_r[0:-2]
    # interpolate 
    fR_stt = CubicSpline(M_stt,R_stt*A_r, bc_type='natural')
    m_stt = np.linspace(min(M_stt),max(M_stt), num=1000)

    print("------------------")
    m_gr1 = np.linspace(min(M_gr),min(M_stt), num=200,endpoint = False)
    m_gr2 = np.linspace(max(M_stt),max(M_gr), num=200)
    return fR_stt(m_stt),m_stt,fR_gr,m_gr1,m_gr2 


def Interpolated_stt(X,Y):
    # Define some points:
    points = np.array([X.tolist(),
                    Y.tolist()]).T  # a (nbre_points x nbre_dim) array
    # Linear length along the line:
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)/distance[-1]
    # Interpolation for different methods:
    alpha = np.linspace(0, 1, 2000)


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

def arrange_sgb(m_stt,fR_stt):
    new_m_stt = []
    new_fR_stt = []
    unstable_m_stt = []
    unstable_R_stt = []
    for i in range(len(m_stt)):
        if i == len(m_stt)-1:
            if (m_stt[i]>m_stt[i-1]):
                new_m_stt.append(m_stt[i])
                new_fR_stt.append(fR_stt[i])
            else:
                unstable_m_stt.append(m_stt[i])
                unstable_R_stt.append(fR_stt[i])
        else:
            if (m_stt[i+1]>m_stt[i]):
                new_m_stt.append(m_stt[i])
                new_fR_stt.append(fR_stt[i])
            else:
                unstable_m_stt.append(m_stt[i])
                unstable_R_stt.append(fR_stt[i])
    return np.asarray(new_m_stt), np.asarray(new_fR_stt),np.asarray(unstable_m_stt), np.asarray(unstable_R_stt)

def arrange_pos(m_stt,fR_stt):
    new_m_stt = []
    new_fR_stt = []
    for i in range(len(m_stt)-1):
        if (m_stt[i+1]>m_stt[i]):
            new_m_stt.append(m_stt[i])
            new_fR_stt.append(fR_stt[i])
        else:
            # new_m_stt.append(m_stt[i])
            # new_fR_stt.append(fR_stt[i])
            return np.asarray(new_m_stt), np.asarray(new_fR_stt)
    return np.asarray(new_m_stt), np.asarray(new_fR_stt)

def arrange_neg(m_stt,fR_stt):
    new_m_stt = []
    new_fR_stt = []
    for i in range(len(m_stt)-1):
        if (m_stt[i+1]>m_stt[i]):
            new_m_stt.append(m_stt[i])
            new_fR_stt.append(fR_stt[i])
        else:
            continue
    return np.asarray(new_m_stt), np.asarray(new_fR_stt)

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

    beta_st = float(datafile[11:15]) + float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29])+ float("0."+ datafile[30:33])
    eos = datafile[4]
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))
    
    with open(directory) as f:

        
        for line in f:
            if len(line.strip()) == 0 :
                continue
            this_line = line.split()
            if 'isScalarized' in this_line:
                continue
            
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
                phi_c_stable.append(float(this_line[8]))
                rho_c_stable.append(float(this_line[9]))
                ctr.append(float(this_line[10]))
            
            mass_Baryon.append(float(this_line[4]))
            mass_ADM.append(float(this_line[5]))
            A_r.append(float(this_line[6]))
            radius.append(float(this_line[7]))
            phi_c.append(float(this_line[8]))
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
        phi_c = np.asarray(phi_c)
        rho_c = np.asarray(rho_c)
        mass_Baryon_stable = np.asarray(mass_Baryon_stable)
        mass_ADM_stable = np.asarray(mass_ADM_stable)
        A_r_stable = np.asarray(A_r_stable)
        radius_stable = np.asarray(radius_stable)
        phi_c_stable = np.asarray(phi_c_stable)
        rho_c_stable = np.asarray(rho_c_stable)
        ctr = np.asarray(ctr)
        
        #fix the first star
        if (mass_ADM[1]>mass_ADM[0]):
            massIncreased[0] = 1
        else:
            massIncreased[0] = 0

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, \
            A_r_stable, radius_stable,phi_c_stable, rho_c_stable, ctr, beta_st, mass_st, eos

def dataReader2(directory):

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
    counter = 0
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

            # if(abs( (float(this_line[5])) -mas)<1e-15 or abs( float(this_line[7])-rad)<1e-15):
            #     continue

            if (Repeated(mass_ADM,float(this_line[5]))):
                continue

            # if (int(this_line[1]) ==1 and counter <= 5):
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
            #         mass_Baryon.append(float(this_line[4]))
            #         mass_ADM.append(float(this_line[5]))
            #         A_r.append(float(this_line[6]))
            #         radius.append(float(this_line[7]))
            #         phi_c.append(float(this_line[8]))
            #         rho_c.append(float(this_line[9]))
            #         ctr.append(float(this_line[10]))

            #         mass_Baryon = np.asarray(mass_Baryon)
            #         mass_ADM = np.asarray(mass_ADM)
            #         A_r = np.asarray(A_r)
            #         radius = np.asarray(radius)
            #         phi_c = np.asarray(phi_c)
            #         rho_c = np.asarray(rho_c)
            #         ctr = np.asarray(ctr)
            #         return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos
            #     else:
            #         if counter>= 5:
            #             stable_branch = False
            #         counter += 1
            #         continue
            mass_Baryon.append(float(this_line[4]))
            mass_ADM.append(float(this_line[5]))
            A_r.append(float(this_line[6]))
            radius.append(float(this_line[7]))
            phi_c.append(float(this_line[8]))
            rho_c.append(float(this_line[9]))
            ctr.append(float(this_line[10]))

            # rad = float(this_line[7])
            # mas = float(this_line[5])
        # turning to arrays
        isScalarized = np.asarray(isScalarized)
        massIncreased = np.asarray(massIncreased)
        isMonotoneInside = np.asarray(isMonotoneInside)
        isMonotoneOutside = np.asarray(isMonotoneOutside)
        
        mass_Baryon = np.asarray(mass_Baryon)
        mass_ADM = np.asarray(mass_ADM)
        A_r = np.asarray(A_r)
        radius = np.asarray(radius)
        phi_c = np.asarray(phi_c)
        rho_c = np.asarray(rho_c)
        ctr = np.asarray(ctr)
        
        

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos

def dataReaderMendes(directory):

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

    flag = 0
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))
    with open(directory) as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            this_line = line.split()
            if 'isScalarized' in this_line:
                continue

            # if(abs( (float(this_line[5])) -mas)<1e-15 or abs( float(this_line[7])-rad)<1e-15):
            #     continue
            mass_Baryon.append(float(this_line[4]))
            mass_ADM.append(float(this_line[5]))
            A_r.append(float(this_line[6]))
            radius.append(float(this_line[7]))
            phi_c.append(float(this_line[8]))
            rho_c.append(float(this_line[9]))
            ctr.append(float(this_line[10]))
            
            
            if (int(this_line[1]) ==0 and flag ==0):
                print("inside")
                flag =1
                continue

            if (int(this_line[1]) ==1 ):
                
                mass_Baryon = np.asarray(mass_Baryon)
                mass_ADM = np.asarray(mass_ADM)
                A_r = np.asarray(A_r)
                radius = np.asarray(radius)
                phi_c = np.asarray(phi_c)
                rho_c = np.asarray(rho_c)
                ctr = np.asarray(ctr)
                print(mass_ADM)
                print("end")
                
                return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos, 
            if (Repeated(mass_ADM,float(this_line[5]))):
                    continue
            # mass_Baryon.append(float(this_line[4]))
            # mass_ADM.append(float(this_line[5]))
            # A_r.append(float(this_line[6]))
            # radius.append(float(this_line[7]))
            # phi_c.append(float(this_line[8]))
            # rho_c.append(float(this_line[9]))
            # ctr.append(float(this_line[10]))


        
        # mass_Baryon = np.asarray(mass_Baryon)
        # mass_ADM = np.asarray(mass_ADM)
        # A_r = np.asarray(A_r)
        # radius = np.asarray(radius)
        # phi_c = np.asarray(phi_c)
        # rho_c = np.asarray(rho_c)
        # ctr = np.asarray(ctr)
        
        

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos

def dataReaderMendes_neg(directory):

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
    MADM_prev =0
    stable =0
    beta_st = float(datafile[11:15]) + float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29])+ float("0."+ datafile[30:33])
    eos = datafile[4]
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))
    
    with open(directory) as f:
        for line in f:
            if len(line.strip()) == 0 :
                continue
            this_line = line.split()
            if 'isScalarized' in this_line:
                continue
            
            if (Repeated(mass_ADM,float(this_line[5]))):
                continue
            # if(float(this_line[5])<MADM_prev):
            #     stable =1
            # else:
            #     stable =0

            # if (stable):
   
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
            # MADM_prev = float(this_line[5])
        # turning to arrays
        isScalarized = np.asarray(isScalarized)
        massIncreased = np.asarray(massIncreased)
        isMonotoneInside = np.asarray(isMonotoneInside)
        isMonotoneOutside = np.asarray(isMonotoneOutside)
        mass_Baryon = np.asarray(mass_Baryon)
        mass_ADM = np.asarray(mass_ADM)
        A_r = np.asarray(A_r)
        radius = np.asarray(radius)
        phi_c = np.asarray(phi_c)
        rho_c = np.asarray(rho_c)
        ctr = np.asarray(ctr)
        
        # #fix the first star
        # if (mass_ADM[1]>mass_ADM[0]):
        #     massIncreased[0] = 1
        # else:
        #     massIncreased[0] = 0

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos


def Prob_vs_Mass(m,fR,Prob2d,norm,probability_data,GR): # returns
    minM = 0
    if GR == 1:
        ProbvsMass = []
        for i in m:
            Mass = i/pp_mass_conv
            Radius = fR(i)*r_scale
            Probability = Prob2d.ev(Radius,Mass)/norm
            # if i/pp_mass_conv < 0.5 or i/pp_mass_conv>2.5:
            #     Probability = 0

            # if probability_data == "MRprob_X5.dat" or probability_data == "MRprob_X7.dat":
            #     if i/pp_mass_conv<0.5:
            #         Probability =0
            if probability_data == "MRprob_M28.dat" or probability_data == "MRprob_OmCen.dat" or probability_data == "MRprob_NGC6304.dat" or probability_data == "MRprob_NGC6397.dat" or probability_data == "MRprob_M30.dat" or probability_data =="MRprob_M13.dat":
                if Mass<0.5125:
                    Probability = 0
                minM = 0.5125   
            elif probability_data == "MRprob_X7.dat" or probability_data == "MRprob_X5.dat":
                if Radius >14:
                    Probability=0
                if Mass<0.5:
                    Probability =0
                minM = 0.5
                if Mass>2.1:
                    Probability =0
            else:
                if Mass<0.6:
                    Probability =0
                minM = 0.6           
            ProbvsMass.append(Probability)

        ProbvsMass = np.asarray(ProbvsMass)
        return ProbvsMass,minM
    if GR == 0:
        ProbvsMass = []
        for i in range(len(m)):
            Mass = m[i]/pp_mass_conv
            Radius = fR[i]*r_scale
            Probability = Prob2d.ev(Radius,Mass)/norm
            # if m[i]/pp_mass_conv < 0.5 or m[i]/pp_mass_conv > 2.5:
            #     Probability = 0
            # if probability_data == "MRprob_X5.dat" or probability_data == "MRprob_X7.dat":
            #     if m[i]/pp_mass_conv>2.1:
            #         Probability =0
            
            # if Mass>2.0:
            #     Probability = 0
            if probability_data == "MRprob_M28.dat" or probability_data == "MRprob_OmCen.dat" or probability_data == "MRprob_NGC6304.dat" or probability_data == "MRprob_NGC6397.dat" or probability_data == "MRprob_M30.dat" or probability_data =="MRprob_M13.dat":
                if Mass<0.5125:
                    Probability = 0
                minM = 0.5125
            elif probability_data == "MRprob_X7.dat" or probability_data == "MRprob_X5.dat":
                if Radius >14:
                    Probability=0
                if Mass<0.5:
                    Probability =0
                minM = 0.5
                if Mass>2.1:
                    Probability =0
            else:
                if Mass<0.6:
                    Probability =0
                minM = 0.6  
            ProbvsMass.append(Probability)

        ProbvsMass = np.asarray(ProbvsMass)
        return ProbvsMass,minM
        
# def integrate(fR_stt,m_stt,fR_gr,m_gr1,m_gr2,Prob2d):


prefix = "C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/"
GR_dir_eos3 = prefix + "GR/eos5/eos_5_beta_0000p000_mass_0000p000"
def plot(directory, directoryGR , Model,GR_on=1, Madm_R = 1):

    ### plots desired figures.

    if (GR_on==1):
        mass_BaryonGR, mass_ADMGR, _, radiusGR, _, rho_cGR, mass_BaryonGR_stable, mass_ADMGR_stable, _, \
        radiusGR_stable, _,rho_cGR_stable ,_,_,_,_\
        = dataReader(directoryGR)
    

        
    # mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, mass_Baryon_stable, mass_ADM_stable, A_r_stable, \
    # radius_stable,phi_c_stable, rho_c_stable, _, beta_st, mass_st, eos \
    # = dataReader(directory)

    # mass_Baryon_stable, mass_ADM_stable, A_r_stable, radius_stable,phi_c, rho_c_stable, ctr, beta_st, mass_st, eos = dataReaderMendes(directory)

    mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos = dataReaderMendes_neg(directory)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.serif": ["Palatino"],
        "xtick.labelsize": 15,
        "ytick.labelsize": 15
    })
    # if beta_st>50:
    #interpolate functions/ returns arrays
    fR_stt,m_stt,fR_gr,m_gr1,m_gr2 =interpolated(radius,mass_ADM,A_r,radiusGR_stable,mass_ADMGR_stable,Model)

    
    # else:
    #     fR_stt,m_stt,fR_gr,m_gr1,m_gr2 =interpolated2(radius_stable,mass_ADM_stable,A_r_stable,radiusGR_stable,mass_ADMGR_stable)

    if (Madm_R==1):
        # maxR = [9.47651,11.08725,9.691275,9.90604,10.7651,10.65772,11.125,9.125,11.375,9.875,11.875,9.625,10.436,10.832]
        # maxMass = [1.592617,1.884564,1.631544,1.67047,1.826174,1.787248,1.1625,0.9375,1.2875,1.0125,0.9875,1.0875,0.5,1.46]
        plt.figure()
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13, 11)
        plt.ylim([5,20])
        plt.xlim([0.1,3.1])



        if (GR_on==1):
            # plt.plot (mass_ADMGR/pp_mass_conv,radiusGR*r_scale,'k--', label = 'GR unstable')
            # plt.plot (mass_ADMGR_stable/pp_mass_conv, radiusGR_stable*r_scale, 'k', label = 'GR stable')
            if len(m_gr1)!=0:
                plt.plot(m_gr1/pp_mass_conv, fR_gr(m_gr1)*r_scale,'k',label = 'GR stable')

            if len(m_gr2)!=0:
                plt.plot(m_gr2/pp_mass_conv, fR_gr(m_gr2)*r_scale,'k')

        plt.ylabel ("Radius (km)",fontsize=20)
        plt.xlabel(r"$M_{ADM} / M_{sun}$",fontsize=20)
        title = (r"$M_{ADM}$ vs R \\ for $\beta=$"+str(beta_st*beta_conv)+r" $ m_S =$" + str(mass_st))
        plt.title(title, fontsize=20)
        # plt.plot (mass_ADM/pp_mass_conv,radius*A_r*r_scale, 'co', label = 'STT unstable')
        
        # plt.plot (mass_ADM_stable/pp_mass_conv,radius_stable*A_r_stable*r_scale, 'b', label = 'STT stable')
        plt.plot(m_stt/pp_mass_conv, fR_stt*r_scale, 'r',label='STT stable')

        # plt.plot(maxR,maxMass,'k.')

        # plt.vlines(x = (m_stt.min()/pp_mass_conv), ymin = fR_stt[-1]*r_scale, ymax = fR_gr((m_gr1[-1]))*r_scale, color = 'b',linestyle = '--')
        # plt.vlines(x = (m_stt.max()/pp_mass_conv),ymin = fR_stt[0]*r_scale, ymax = fR_gr((m_gr2[0]))*r_scale, color = 'b',linestyle = '--')
        plt.legend()
    return fR_stt,m_stt,fR_gr,m_gr1,m_gr2,beta_st, mass_st ,phi_c[0]


def main(relaxation_data,GRdata, probability_data,savedir,type,Model):

    fR_stt,m_stt,fR_gr,m_gr1,m_gr2,beta_st, mass_st,phi_c = plot(relaxation_data,GRdata,Model)
    Prob_spl,norm = contourplot(probability_data,savedir,type,beta_st,mass_st)

    plt.figure()
    if len(m_gr1)!=0:
        ProbvsMass_gr1,_ = Prob_vs_Mass(m_gr1,fR_gr,Prob_spl,norm, probability_data, GR=1)
        plt.plot(m_gr1/pp_mass_conv,ProbvsMass_gr1,'k',label ='GR')

    ProbvsMass,minM = Prob_vs_Mass(m_stt,fR_stt,Prob_spl,norm, probability_data, GR = 0)
    plt.plot(m_stt/pp_mass_conv,ProbvsMass,'r.',label ='STT')

    if len(m_gr2)!=0:
        ProbvsMass_gr2,_ = Prob_vs_Mass(m_gr2,fR_gr,Prob_spl,norm,probability_data,GR = 1)
        plt.plot(m_gr2/pp_mass_conv,ProbvsMass_gr2,'k',label ='GR')

    plt.title("probability for " + probability_data +" in theory beta "+str(beta_st)+" mass " + str(mass_st))
    plt.legend()
    plt.savefig(savedir+"prob_"+probability_data+"_beta_"+str(beta_st)+".png", dpi=300, bbox_inches='tight')
    plt.close()
    if Model == 'linear/positive':
        maxM = max(m_stt)
    else:
        if len(m_gr1)!=0 and len(m_gr2)!=0:
            max1 = max(m_gr1)
            max3 = max(m_gr2)
            max2 = max(m_stt)
            maxM = max(max1,max2,max3)
        elif len(m_gr1)==0 and len(m_gr2)!=0:
            max3 = max(m_gr2)
            max2 = max(m_stt)
            maxM = max(max2,max3)
        elif len(m_gr1)!=0 and len(m_gr2)==0:
            max1 = max(m_gr1)
            max2 = max(m_stt)
            maxM = max(max2,max1)
        else:
            maxM = max(m_stt)
    minM = 0.1
    prior = 1/(maxM/pp_mass_conv-minM)

    if len(m_gr1)!=0:
        Integral1 = simpson(ProbvsMass_gr1, m_gr1/pp_mass_conv)
    else:
        Integral1 =0
    Integral2 = simpson(ProbvsMass, m_stt/pp_mass_conv)
    if len(m_gr2)!=0:
        Integral3 = simpson(ProbvsMass_gr2, m_gr2/pp_mass_conv) ##or simpson should check
    else:
        Integral3 =0
    Integral = (Integral1+Integral2+Integral3)*prior
    
    if maxM/pp_mass_conv< 2.0:
        Integral = 0
    return Integral,beta_st, mass_st,phi_c