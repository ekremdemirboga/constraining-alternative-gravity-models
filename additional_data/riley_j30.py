import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')



def read30(datafile):
    name = "J0730"
    Mass = []
    Rad = []
    prob = []
    weight = []
    with open(datafile) as f:
        ##reads the neutron stars mass-radius-probability data files in a given directory.
            for line in f:
                this_line = line.split()
                Mass.append(float(this_line[2]))
                Rad.append(float(this_line[3]))
                prob.append( (float(this_line[1])) )
                weight.append(float(this_line[0]))
    return np.array(Mass),np.array(Rad),np.array(prob),np.array(prob), name
    
Mass,Rad,prob,weight, name = read30("j30/riley/j30_riley.txt")
print(weight)
M_grid = np.linspace(min(Mass), max(Mass),40)
R_grid = np.linspace(min(Rad), max(Rad),40)
n = len(M_grid)
dx = abs(M_grid[1]-M_grid[0])
dy = abs(R_grid[1]-R_grid[0])
Z = np.zeros((len(M_grid),len(R_grid)))
for i in range(len(Mass)):
    x = int( (Mass[i] - min(M_grid)) /dx) ## index x in the grid
    y = int( (Rad[i] - min(R_grid) ) /dy) ## index y in the grid
    Z[x,y] = Z[x,y] + prob[i]*weight[i] ## add the weight to point Z(x,y)


plt.pcolormesh(R_grid,M_grid,Z) ## NOTE: X and Y is reversed for this function
plt.colorbar()
plt.title("riley j0030")
plt.savefig("mr_j30_riley.png",dpi=300, bbox_inches='tight')

### write the data 
with open('j0030_likelihood_riley.txt', 'w') as f:
    f.write('# =============================================\n')
    f.write('# Column 1: Equatorial circumferential radius in kilometers \n')
    f.write('# Column 2: Gravitational mass in solar masses \n')
    f.write('# Column 3: likelihood \n')
    f.write('# =============================================\n')
    f.write('# =============================================\n')
    for i in range(n):
        for j in range(n):
            f.write(str(M_grid[j])+ " " + str(R_grid[i])+ " " + str(Z[j][i])+"\n")