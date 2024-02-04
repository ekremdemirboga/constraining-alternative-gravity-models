import numpy as np
import glob
import os

prefix = r"C:/Users/ekrem/Desktop/PositiveBeta/likelihood/likelihood_calculations/"
relaxation_data_eos5 =sorted(glob.glob(prefix+r"massless/linear/positive/data/relaxation_data/eos_5*"))
extra =[]
for i in range(len(relaxation_data_eos5)):
    extra.append(relaxation_data_eos5[i] + "/_extra")
    relaxation_data_eos5[i] = relaxation_data_eos5[i] + "/" +relaxation_data_eos5[i][-33:]
    relaxation_data_eos5[i] = relaxation_data_eos5[i].replace(os.sep,'/')
    

for i in range(len(extra)):
    extra[i] = extra[i].replace(os.sep,'/')

def nameArrange(directory):
    ## takes directory as string and gives the name of the datafile. directory = ..../datafile
    datafile = ""
    for i in directory[::-1]:
        if (i == '/'):
            datafile = datafile[::-1]
            break
        datafile += i
    
    return datafile

## this will read the normal file.
def dataReader(directory,extra):
    datafile=nameArrange(directory)
    beta_st = float(datafile[11:15]) #+ float("0."+ datafile[16:19])
    mass_st = float(datafile[25:29]) #+ float("0."+ datafile[30:33])
    eos = datafile[4]
    if ((beta_st !=0)):
        print("beta_st: " + str(beta_st))
        print("mass_st: " + str(mass_st))
        print("eos: " + str(eos))


    with open(extra) as (e):
        extra_data =[]
        for line in reversed(e.readlines()):
            extra_data.append(line)
        
    with open(directory) as f:
        main_data = []
        for line in f:
            this_line = line.split()
            if 'isScalarized' in this_line:
                continue
            main_data.append(line)

    # data = extra_data + main_data  
    # print(data)
    with open(directory, 'w') as d:
        for line in extra_data:
            d.write(line)
        for line in main_data:
            d.write(line)

        
for i in range(len(relaxation_data_eos5)):
        
        dataReader(relaxation_data_eos5[i],extra[i])