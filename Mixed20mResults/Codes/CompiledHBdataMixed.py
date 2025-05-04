import numpy as np



data = np.loadtxt('Plots/DataFiles/Combined.csv', dtype = 'int',delimiter = ',',skiprows=1)

HB_water = data[:,1]
HB_tfsi = data[:,2]
HB_total = data[:,3]
Cood_Li = data[:,4]
Cood_Zn = data[:,5]
Cood_total_cation = data[:,6]
HB_F = data[:,7]
HB_N = data[:,8]
HB_Ot = data[:,9]

# GENERAL H-BOND TABLE

row_combo = np.array([])

# HB=0,1,2,3,4 with other water
for i in range(0,5):   
    # Zn=0,1
    for j in range(0,2):
        sum=0
        # TFSI=0,1,2
        for k in range(0,3):
            #print(np.round(np.average((HB_water==i)*(Cood_Zn==j)*(HB_tfsi==k))*100,2), end = ' ')
            # sum+=np.round(np.average((HB_water==i)*(Cood_Li==j)*(HB_tfsi==k))*100,2)
            row_combo = np.append(row_combo,np.round(np.average((HB_water==i)*(Cood_total_cation==j)*(HB_tfsi==k))*100,2))
        #print(np.round(sum,2), end = ' ')
        # row_combo = np.append(row_combo,np.round(sum,2))
    #print('\n')
    

with open("Plots/DataFiles/HBabsValues.txt", "ab") as f:
    np.savetxt(f, row_combo.reshape(1, row_combo.shape[0]),fmt = '%.2f',delimiter=' ')
    




# TABLE FOR ISOLATED WATER

iso = (HB_water==0)


row_combo2 = np.array([])

# N=0,1
for i in range(0,2):
    # Otfsi=0,1,2
    for j in range(0,3):
        # Zn=0,1
        for k in range(0,2):
            # F=0,1,2
            for l in range(0,3):
                row_combo2 = np.append(row_combo2, np.round(np.average((HB_N[iso]==i)*(HB_Ot[iso]==j)*(Cood_total_cation[iso]==k)*(HB_F[iso]==l))*100,2))

                
with open("Plots/DataFiles/IsoHBabsValues.txt", "ab") as f2:
    np.savetxt(f2, row_combo2.reshape(1, row_combo2.shape[0]),fmt = '%.2f',delimiter=' ')




# TABLE FOR 1-HB WATER

one = (HB_water==1)


row_combo3 = np.array([])

# N=0,1
for i in range(0,2):
    # Otfsi=0,1,2
    for j in range(0,3):
        # Zn=0,1
        for k in range(0,2):
            # F=0,1,2
            for l in range(0,3):
                row_combo3 = np.append(row_combo3, np.round(np.average((HB_N[one]==i)*(HB_Ot[one]==j)*(Cood_total_cation[one]==k)*(HB_F[one]==l))*100,2))

                
with open("Plots/DataFiles/OneHBabsValues.txt", "ab") as f3:
    np.savetxt(f3, row_combo3.reshape(1, row_combo2.shape[0]),fmt = '%.2f',delimiter=' ')