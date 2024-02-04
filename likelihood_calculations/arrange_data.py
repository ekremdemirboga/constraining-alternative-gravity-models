
import glob
import os

def Repeated(Mass_array,Mass_value):
    for i in range(len(Mass_array)):
        if (abs(Mass_value - Mass_array[i])<1e-10):
            print("repeat!!!")
            return True
    return False


def dataReader(directory):

    
    #Takes the input as "only" the name of the file and reads the data files 
    #with the names format: eos_x_beta_xxxxpxxx_mass_xxxxpxxx
    
    #Order: isScalarized | massIncreased | isMonotoneInside | isMonotoneOutside | 
    #mass_Baryon | mass_ADM | AST(radius) | radius  |  phi_c | rho_c | ctr |
    

    massIncreased = []
    mass_ADM = []
    old_mass = 9999
    flag = True
    with open(directory, "r") as f:
        lines = []
        for line in f:
            this_line = line.split()

            if 'isScalarized' in this_line:
                continue
            if ((this_line) ==["\n"] or this_line ==[]):
                print("EMPTY")
                continue
            if (Repeated(mass_ADM,float(this_line[5]))):
                continue
            
            current_mass = float(this_line[5])
            if (current_mass > old_mass):
                flag = False
            old_mass = current_mass
            mass_ADM.append(float(this_line[5]))
            ## categorizing according to instabilities
            if ((int(this_line[1]) == 1) and flag): # if stable 
                lines.append(line)
    with open(directory, "w") as f:
        for line in lines:
            f.write(line)
    return mass_ADM

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/"
relaxation_data_eos3 =sorted(glob.glob(prefix+r"massless/test_eos3/RELAXATION_DATA/eos_3*"))
#GR data
GR_dir_eos3 = prefix + r"GR/eos3/eos_3_beta_0000p000_mass_0000p000"
save_dir = prefix + r"massless/test_eos3/plots/"
for i in range(len(relaxation_data_eos3)):
    relaxation_data_eos3[i] = relaxation_data_eos3[i] + "/" +relaxation_data_eos3[i][-33:]
    relaxation_data_eos3[i] = relaxation_data_eos3[i].replace(os.sep,'/')

for relax_data in relaxation_data_eos3:
    dataReader(relax_data)