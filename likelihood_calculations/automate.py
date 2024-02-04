from likelihood_funcs import main
import matplotlib.pyplot as plt
import glob
import os
from scipy.integrate import simpson
from plot_functions2 import plot2
## Neutron Star Data
datafiles1= ["MRprob_1608.dat","MRprob_1724.dat","MRprob_1731.dat", "MRprob_1745.dat", "MRprob_1748.dat", "MRprob_1820.dat"]
datafiles2 = ["MRprob_M13.dat", "MRprob_M28.dat", "MRprob_M30.dat" , "MRprob_NGC6304.dat", "MRprob_NGC6397.dat","MRprob_OmCen.dat" ,"MRprob_X5.dat","MRprob_X7.dat"]

# relaxation data && GR data
prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/"

def auto(Model,eos):
    if Model == 'quadratic':
        if eos ==3:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/quadratic/eos3/relaxation_data/eos_3*"))
            save_dir = prefix + r"massless/quadratic/eos3/plots/"
            GR_dir = prefix + r"GR/eos3/eos_3_beta_0000p000_mass_0000p000"
        elif eos ==4:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/quadratic/eos4/relaxation_data/eos_4*"))
            save_dir = prefix + r"massless/quadratic/eos4/plots/"
            GR_dir = prefix + r"GR/eos4/eos_4_beta_0000p000_mass_0000p000"
        elif eos==5:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/quadratic/eos5/relaxation_data/eos_5*"))
            save_dir = prefix + r"massless/quadratic/eos5/plots/"
            GR_dir = prefix + r"GR/eos5/eos_5_beta_0000p000_mass_0000p000"
    elif Model == 'linear/positive':
        name='positive'
        if eos ==1:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/positive/eos1/relaxation_data/eos_1*"))
            save_dir = prefix + r"massless/linear/positive/eos1/plots/"
            GR_dir = prefix + r"GR/eos1/eos_1_beta_0000p000_mass_0000p000"
        if eos ==2:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/positive/eos2/relaxation_data/eos_2*"))
            save_dir = prefix + r"massless/linear/positive/eos2/plots/"
            GR_dir = prefix + r"GR/eos2/eos_2_beta_0000p000_mass_0000p000"
        if eos ==3:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/positive/eos3/relaxation_data/eos_3*"))
            save_dir = prefix + r"massless/linear/positive/eos3/plots/"
            GR_dir = prefix + r"GR/eos3/eos_3_beta_0000p000_mass_0000p000"
        elif eos ==4:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/positive/eos4/relaxation_data/eos_4*"))
            save_dir = prefix + r"massless/linear/positive/eos4/plots/"
            GR_dir = prefix + r"GR/eos4/eos_4_beta_0000p000_mass_0000p000"
        elif eos==5:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/positive/eos5/relaxation_data/eos_5*"))
            save_dir = prefix + r"massless/linear/positive/eos5/plots/"
            GR_dir = prefix + r"GR/eos5/eos_5_beta_0000p000_mass_0000p000"
    elif Model == 'linear/negative':
        name='negative'
        if eos ==1:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/negative/eos1/relaxation_data/eos_1*"))
            save_dir = prefix + r"massless/linear/negative/eos1/plots/"
            GR_dir = prefix + r"GR/eos1/eos_1_beta_0000p000_mass_0000p000"
        if eos ==2:
            print("here")
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/negative/eos2/relaxation_data/eos_2*"))
            save_dir = prefix + r"massless/linear/negative/eos2/plots/"
            GR_dir = prefix + r"GR/eos2/eos_2_beta_0000p000_mass_0000p000"
        if eos ==3:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/negative/eos3/relaxation_data/eos_3*"))
            save_dir = prefix + r"massless/linear/negative/eos3/plots/"
            GR_dir = prefix + r"GR/eos3/eos_3_beta_0000p000_mass_0000p000"
        elif eos ==4:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/negative/eos4/relaxation_data/eos_4*"))
            save_dir = prefix + r"massless/linear/negative/eos4/plots/"
            GR_dir = prefix + r"GR/eos4/eos_4_beta_0000p000_mass_0000p000"
        elif eos==5:
            relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/linear/negative/eos5/relaxation_data/eos_5*"))
            save_dir = prefix + r"massless/linear/negative/eos5/plots/"
            GR_dir = prefix + r"GR/eos5/eos_5_beta_0000p000_mass_0000p000"
    elif Model == 'GR':
        name='GR'
        if eos ==1:
            save_dir = prefix + r"massless/GR/eos1/plots/"
            GR_dir = prefix + r"GR/eos1/eos_1_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]
        if eos ==2:
            save_dir = prefix + r"massless/GR/eos2/plots/"
            GR_dir = prefix + r"GR/eos2/eos_2_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]
        if eos ==3:
            save_dir = prefix + r"massless/GR/eos1/plots/"
            GR_dir = prefix + r"GR/eos1/eos_1_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]
        if eos ==2:
            save_dir = prefix + r"massless/GR/eos2/plots/"
            GR_dir = prefix + r"GR/eos2/eos_2_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]
        if eos ==3:
            save_dir = prefix + r"massless/GR/eos3/plots/"
            GR_dir = prefix + r"GR/eos3/eos_3_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]
        elif eos ==4:
            save_dir = prefix + r"massless/GR/eos4/plots/"
            GR_dir = prefix + r"GR/eos4/eos_4_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]
        elif eos==5:
            save_dir = prefix + r"massless/GR/eos5/plots/"
            GR_dir = prefix + r"GR/eos5/eos_5_beta_0000p000_mass_0000p000"
            relaxation_data_eos3 =[GR_dir]

    #rearrange data
    if  Model != 'GR':
        for i in range(len(relaxation_data_eos3)):
            relaxation_data_eos3[i] = relaxation_data_eos3[i] + "/" +relaxation_data_eos3[i][-33:]
            relaxation_data_eos3[i] = relaxation_data_eos3[i].replace(os.sep,'/')
        

    totProb = 1
    p =[]
    beta =[]
    phi = []

    PLOT_FLAG = 0
    LIKELIHOOD_FLAG = 1
    for relax_data in relaxation_data_eos3:
        print("FOR " +relax_data+" .....")
        if (PLOT_FLAG):
            ###Plot Everything
            
            plot2(relax_data,prefix+"massless/"+Model+"/eos"+str(eos)+"/mr/",GR_dir)
        if(LIKELIHOOD_FLAG):
            ##Then Likelihood
            for prob_data in datafiles1:
                type =1
                print(prob_data)
                Integral,beta_st,mass_st,_ = main(relax_data,GR_dir,prob_data,save_dir,type,Model)
                print(Integral)
                totProb *= Integral
            for prob_data in datafiles2:
                type =2
                print(prob_data)
                Integral,beta_st,mass_st,_ = main(relax_data,GR_dir,prob_data,save_dir,type,Model)
                print(Integral)
                totProb *= Integral
            print("probability of data in theory with beta "+str(beta_st)+" is found to be "+str(totProb))
            beta.append(beta_st)
            p.append(totProb)
            totProb=1 #reset

    if (LIKELIHOOD_FLAG):
        ## normalize

        with open(prefix+'massless/'+Model +'/eos'+str(eos)+'/likelihood_'+name+'_eos'+str(eos)+'.txt', 'w') as f:
            f.write("beta   |   likelihood  |\n")
            for i in range(len(beta)):
                line = str(beta[i]) + " " + str(p[i])
                f.write(line)
                f.write('\n')
        norm =  simpson(p, beta)
        print("norm is "+str(norm))
        # p = p/norm
        plt.figure()
        plt.plot(beta,p/norm)
        plt.xlabel("beta")
        plt.ylabel("p")
        plt.grid()
        plt.savefig(prefix + "massless/"+Model+"/eos"+str(eos)+"/p_beta.png")

        
EOSs= [1,2,3,4,5]
Models = ['linear/negative', 'linear/positive', 'GR']

for eos in EOSs:
    for model in Models:
        auto(model,eos)