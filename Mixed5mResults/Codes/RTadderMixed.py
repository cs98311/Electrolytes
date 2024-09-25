import numpy as np


z = np.loadtxt("RT/Zrt.txt").T[3].astype(int)
with open("RT/RtimesZ.txt","a") as f0:
	print(*z, file=f0)


l = np.loadtxt("RT/Lrt.txt").T[3].astype(int)
with open("RT/RtimesL.txt","a") as f1:
	print(*l, file=f1)

f = np.loadtxt("RT/Frt.txt").T[6].astype(int)
with open("RT/RtimesF.txt","a") as f2:
	print(*f, file=f2)

o = np.loadtxt("RT/Ort.txt").T[6].astype(int)
with open("RT/RtimesO.txt","a") as f3:
	print(*o, file=f3)

w = np.loadtxt("RT/Wrt.txt").T[5].astype(int)
with open("RT/RtimesW.txt","a") as f4:
	print(*w, file=f4)


dd1 = np.loadtxt("RT/Drt.txt").T[0].astype(int)
dd2 = np.loadtxt("RT/Drt.txt").T[1].astype(int)
with open("RT/RtimesD2.txt","a") as f5:
	print(*dd2[dd1==2], file=f5)
with open("RT/RtimesD1.txt","a") as f6:
	print(*dd2[dd1==1], file=f6)
with open("RT/RtimesD0.txt","a") as f7:
	print(*dd2[dd1==0], file=f7)


d0 = np.loadtxt("RT/RtimesD0.txt").astype(int)
d1 = np.loadtxt("RT/RtimesD1.txt").astype(int)
d2 = np.loadtxt("RT/RtimesD2.txt").astype(int)

with open("Averages/AvgDonorLifetimes.txt","a") as f8:
	print(f"{np.round(np.average(d0),2)},{np.round(np.average(d1),2)},{np.round(np.average(d2),2)}",file=f8)