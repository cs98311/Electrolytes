import numpy as np




data_avg = np.loadtxt('Averages/AvgHBabsValues.csv', dtype='float', delimiter = ',')

row_combo3 = data_avg.reshape(5,6)

lefter_rc3 = np.array([['0'],['1'],['2'],['3'],['4']])

rc_temp3 = np.append(lefter_rc3, row_combo3,axis =1)

header1_rc3 = np.array([[''],[''],['Cation (0)'],[''],[''],['Cation (1)'],['']])
header2_rc3 = np.array([[''],['TFSI (0)'],['TFSI (1)'],['TFSI (2)'],['TFSI (0)'],['TFSI (1)'],['TFSI (2)']])

rc3_semi = np.append(header2_rc3,rc_temp3.T,axis=1)

rc3_final = np.append(header1_rc3,rc3_semi,axis=1)

with open("Tables/TableAbsPercent.csv", "wb") as f:
    np.savetxt(f, rc3_final.T, delimiter=",",fmt = '%s')









### Isolated ###

data_avg4 = np.loadtxt('Averages/AvgIsoHBabsValues.csv', dtype='float', delimiter = ',')

row_combo4 = data_avg4.reshape(6,6)

lefter_rc4 = np.array([['O (0)'],['O (1)'],['O (2)'],['O (0)'],['O (1)'],['O (2)']])
leftest_rc4 = np.array([[''],['N (0)'],[''],[''],['N (1)'],['']])

rc_temp4 = np.append(lefter_rc4, row_combo4,axis =1)
rc_temp4 = np.append(leftest_rc4,rc_temp4,axis=1)

header1_rc4 = np.array([[''],[''],[''],['Cation (0)'],[''],[''],['Cation (1)'],['']])
header2_rc4 = np.array([[''],[''],['F (0)'],['F (1)'],['F (2)'],['F (0)'],['F (1)'],['F (2)']])

rc4_semi = np.append(header2_rc4,rc_temp4.T,axis=1)

rc4_final = np.append(header1_rc4,rc4_semi,axis=1)

with open("Tables/TableIsoAbsPercent.csv", "wb") as f2:
    np.savetxt(f2, rc4_final.T, delimiter=",",fmt = '%s')









### 1HB ###

data_avg5 = np.loadtxt('Averages/AvgOneHBabsValues.csv', dtype='float', delimiter = ',')

row_combo5 = data_avg5.reshape(6,6)

lefter_rc5 = np.array([['O (0)'],['O (1)'],['O (2)'],['O (0)'],['O (1)'],['O (2)']])
leftest_rc5 = np.array([[''],['N (0)'],[''],[''],['N (1)'],['']])

rc_temp5 = np.append(lefter_rc5, row_combo5,axis =1)
rc_temp5 = np.append(leftest_rc5,rc_temp5,axis=1)

header1_rc5 = np.array([[''],[''],[''],['Cation (0)'],[''],[''],['Cation (1)'],['']])
header2_rc5 = np.array([[''],[''],['F (0)'],['F (1)'],['F (2)'],['F (0)'],['F (1)'],['F (2)']])

rc5_semi = np.append(header2_rc5,rc_temp5.T,axis=1)

rc5_final = np.append(header1_rc5,rc5_semi,axis=1)

with open("Tables/TableOneAbsPercent.csv", "wb") as f3:
    np.savetxt(f3, rc5_final.T, delimiter=",",fmt = '%s')