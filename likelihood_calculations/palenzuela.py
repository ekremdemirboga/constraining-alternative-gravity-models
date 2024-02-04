from plot_functions import nameArrange,GR_dir_eos3
import numpy as np
from main import prefix
import glob
import matplotlib.pyplot as plt
import os
print("inside palenzuela")

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
    mass_ADM = []
    A_r = []
    radius = []
    phi_c = []
    rho_c = []
    ctr = []

    beta_st = float(datafile[11:15]) + float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29])+ float("0."+ datafile[30:33])
    eos = datafile[4]
    print("beta_st: " + str(beta_st))
    print("mass_st: " + str(mass_st))
    print("eos: " + str(eos))
    
    with open(directory) as f:
        flag = 0
        for line in f:
            this_line = line.split()
            if 'isScalarized' in this_line:
                continue
            if ((this_line) ==["\n"] or this_line ==[]):
                print("EMPTY")
                if ((beta_st !=0)):
                    print("beta_st: " + str(beta_st))
                    print("mass_st: " + str(mass_st))
                    print("eos: " + str(eos))
                return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos, massIncreased
            isScalarized.append(int(this_line[0]))
            massIncreased.append(int(this_line[1]))
            isMonotoneInside.append(int(this_line[2]))
            isMonotoneOutside.append(int(this_line[3]))
            if ((int(this_line[1]) == 1)):
                flag =1

            if ((int(this_line[1]) == 0) and flag): # if stable 
                mass_Baryon.append(float(this_line[4]))
                mass_ADM.append(float(this_line[5]))
                A_r.append(float(this_line[6]))
                radius.append(float(this_line[7]))
                phi_c.append(abs(float(this_line[8])))
                rho_c.append(float(this_line[9]))
                ctr.append(float(this_line[10]))
                return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos, massIncreased
            
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
        ctr = np.asarray(ctr)

        return mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos, massIncreased


data = sorted(glob.glob(prefix+r"massless/quadratic/eos3/relaxation_data/eos_3*"))
data2 = sorted(glob.glob(prefix+r"massless/quadratic/eos3/relaxation_data/eos_3*"))

for i in range(len(data)):
    data[i] = data[i] + "/" +data[i][-33:]
    data[i] = data[i].replace(os.sep,'/')

for i in range(len(data2)):
    data2[i] = data2[i] + "/" +data2[i][-33:]
    data2[i] = data2[i].replace(os.sep,'/')

betas = []
first_star_compactness = []
betas2 = []
last_star_compactness = []

for directory in data:
    mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos,massIncreased = dataReader(directory)
    betas.append(beta_st)
    first_star_compactness.append(mass_ADM[0]/(radius[0]*A_r[0]))
    # plt.plot(mass_Baryon,np.divide(mass_ADM,radius))
    # plt.show()
for directory in data2:
    mass_Baryon, mass_ADM, A_r, radius,phi_c, rho_c, ctr, beta_st, mass_st, eos,massIncreased = dataReader(directory)
    # print(mass_ADM[-1],massIncreased[-1])
    betas2.append(beta_st)
    last_star_compactness.append(mass_ADM[-1]/(radius[-1]*A_r[-1]))

plt.semilogy(first_star_compactness,betas)
plt.semilogy(last_star_compactness,betas2,'o')
plt.xlabel("M/R")
plt.ylabel("Beta")
plt.title("First Scalarized Stars")
plt.show()