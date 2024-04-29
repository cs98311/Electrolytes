#! usr/bin/env python3

import subprocess
from os.path import exists

# Emptying the folders
subprocess.run("rm Anion/*", shell=True)
subprocess.run("rm Averages/*", shell=True)
subprocess.run("rm Cation/*", shell=True)
subprocess.run("rm Combined/*", shell=True)
subprocess.run("rm Coordinates/*", shell=True)
subprocess.run("rm Exchange/*", shell=True)
subprocess.run("rm Orientation/*", shell=True)
subprocess.run("rm Plots/DataFiles/*", shell=True)
subprocess.run("rm RT/*", shell=True)
subprocess.run("rm Tables/*", shell=True)
subprocess.run("rm Water/*", shell=True)
subprocess.run("rm SysInfo.txt", shell=True)


### .TRR TO .XTC ###
"""
Checking whether .xtc file exists or not, and creating it in case not
since it is much faster to read iteratively than .trr file
"""
if not exists("MDSfiles/nvt3.xtc"):
    subprocess.run("gmx trjconv -f MDSfiles/nvt3.trr -o MDSfiles/nvt3.xtc", shell=True)


### WRITE HEADERS ROWS FOR CSV, TXT FILES ###
# Also creating files in advance that can be appended to by other codes.
subprocess.run("gcc Codes/Header.c -o Codes/Header -lm", shell=True)
subprocess.run("./Codes/Header")


### ITERATION LIMITS ###
"""
"start" is the timestep at which iteration begins
"end" is the timestep at which iteration ends +1
Add 1 to "end" since upper limit is not included in range(x,y)
Ex- start=20, end=101 if timesteps go from 20 to 100
'n' is the total number of timesteps being iterated
"""

start = 0
end = 1

n = end - start


# Total timesteps = 0 to 2000
# Exceptions: 20m Li and 10m Zn : both 10000 timesteps

for i in range(start, end):
    # Run trjconv command and generate coordinate text files (for each iteration)
    subprocess.run(
        "python3 Codes/Coordinates.py", shell=True, input="{}".format(i), text=True
    )

    # Do main calculations
    subprocess.run("gcc Codes/BondCheckerLi.c -o Codes/BondCheckerLi -lm", shell=True)
    subprocess.run(
        "./Codes/BondCheckerLi", input="{}".format(i), capture_output=True, text=True
    )

    # Print cumulative H-bond data (for each iteration) for table creation
    subprocess.run("python3 Codes/CompiledHBdataLi.py", shell=True)

    subprocess.run("python3 Codes/Cycles.py", shell=True)


subprocess.run("python3 Codes/RTadderLi.py", shell=True)


# Averaging
subprocess.run("gcc Codes/Averager.c -o Codes/Averager -lm", shell=True)

# Print percentage of water forming 0,1,2,3,4 H-bonds (averaged out over iterations)
m = 5
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nWater/HBpercent.txt\nAverages/AvgHBpercent.csv".format(n, m),
    capture_output=True,
    text=True,
)


# Print number of H-bond/Coordination with Water,Anion,Lithium,Zinc per water molecule (averaged out over iterations)
p = 4
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nWater/HBperWater.txt\nAverages/AvgHBperWater.csv".format(n, p),
    capture_output=True,
    text=True,
)


# Print out water cluster frequency distribution (averaged out over iterations)
o = 10
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nWater/FreqDistri.txt\nAverages/AvgFreqDistri1.csv".format(n, o),
    capture_output=True,
    text=True,
)


# Print out water cluster percentage distribution (averaged out over iterations)
q = 10
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nWater/ClusterPercent.txt\nAverages/AvgClusterPercent1.csv".format(
        n, q
    ),
    capture_output=True,
    text=True,
)


# Print cumulative H-bond data (averaged out over iterations) for table creation
r = 30
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nPlots/DataFiles/HBabsValues.txt\nAverages/AvgHBabsValues.csv".format(
        n, r
    ),
    capture_output=True,
    text=True,
)


# Print cumulative H-bond data of Isolated water(ie, 0 H-bonds with other water)(averaged out over iterations) for table creation
s = 36
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nPlots/DataFiles/IsoHBabsValues.txt\nAverages/AvgIsoHBabsValues.csv".format(
        n, s
    ),
    capture_output=True,
    text=True,
)


# Print cumulative H-bond data of water forming 1 H-bond with other water (averaged out over iterations) for table creation
s = 36
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nPlots/DataFiles/OneHBabsValues.txt\nAverages/AvgOneHBabsValues.csv".format(
        n, s
    ),
    capture_output=True,
    text=True,
)


# Print donor-only H-bond probablities for water(averaged out over iterations)
t = 3
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nWater/Donor.txt\nAverages/AvgDonorProbs.csv".format(n, t),
    capture_output=True,
    text=True,
)

# Print acceptor-only H-bond probablities for water(averaged out over iterations)
t = 3
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nWater/Acceptor.txt\nAverages/AvgAcceptorProbsOw.csv".format(n, t),
    capture_output=True,
    text=True,
)


# Print probablity of cation having 0,1,2,3,4,5,6 oxygen in coordination (averaged out over iterations)
u = 7
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nCation/Li_Ow.txt\nAverages/AvgLiCoodOw.csv".format(n, u),
    capture_output=True,
    text=True,
)
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nCation/Li_Ot.txt\nAverages/AvgLiCoodOt.csv".format(n, u),
    capture_output=True,
    text=True,
)


# Print acceptor-only H-bond probablities for Otfsi (averaged out over iterations)
t = 3
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nAnion/Ot_HB_accepted.txt\nAverages/AvgAcceptorProbsOt.csv".format(
        n, t
    ),
    capture_output=True,
    text=True,
)

# Print probabilities of Otfsi having 0,1,2 Cations in coordination (averaged out over iterations)
t = 3
subprocess.run(
    "./Codes/Averager",
    input="{}\n{}\nAnion/Ot_Li_Coodn.txt\nAverages/AvgOtCoodLi.csv".format(n, t),
    capture_output=True,
    text=True,
)


### TRANSPOSE FOR EASIER PLOTTING AND ANALYSIS ###
subprocess.run(
    ["python3", "Codes/Transposer.py"],
    input="Water/MyFile.csv\nWater/ClusterSizes.csv\n",
    text=True,
)
subprocess.run(
    ["python3", "Codes/Transposer.py"],
    input="Averages/AvgFreqDistri1.csv\nAverages/AvgFreqDistri.csv\n",
    text=True,
)
subprocess.run(
    ["python3", "Codes/Transposer.py"],
    input="Water/ClusterPercent1.csv\nWater/ClusterPercent.csv\n",
    text=True,
)
subprocess.run(
    ["python3", "Codes/Transposer.py"],
    input="Averages/AvgClusterPercent1.csv\nAverages/AvgClusterPercent.csv\n",
    text=True,
)


### CREATE H-BONDING INFO TABLES ###
subprocess.run("python3 Codes/TableCreateLi.py", shell=True)


# 3D Plotting
def Plot():
    # Print oxygen to oxygen H-bond vectors file ([Starting Point],[Direction])
    subprocess.run("gcc Codes/HBvectors.c -o Codes/HBvectors -lm", shell=True)
    subprocess.run("./Codes/HBvectors")

    # Print Center of Mass coordinates of TFSI (Xcom, Ycom, Zcom)
    subprocess.run("gcc Codes/TFSIcom.c -o Codes/TFSIcom -lm", shell=True)
    subprocess.run("./Codes/TFSIcom")

    # Save the 3D plots as images (.svg, .eps, .jpg)
    subprocess.run("python3 Codes/3DplotsLi.py", shell=True)


# Plot()


# Remove unnecessary files
subprocess.run("rm Water/MyFile.csv", shell=True)
subprocess.run("rm Water/FreqDistri.txt", shell=True)
subprocess.run("rm Averages/AvgFreqDistri1.csv", shell=True)
subprocess.run("rm Codes/Averager", shell=True)
subprocess.run("rm Codes/Header", shell=True)
if exists("Codes/HBvectors") == True:
    subprocess.run("rm Codes/HBvectors", shell=True)
if exists("Codes/TFSIcom") == True:
    subprocess.run("rm Codes/TFSIcom", shell=True)
subprocess.run("rm Codes/BondCheckerLi", shell=True)
subprocess.run("rm Water/ClusterPercent1.csv", shell=True)
subprocess.run("rm Averages/AvgClusterPercent1.csv", shell=True)
subprocess.run("rm Averages/AvgHBabsValues.csv", shell=True)
subprocess.run("rm Averages/AvgIsoHBabsValues.csv", shell=True)
subprocess.run("rm Averages/AvgOneHBabsValues.csv", shell=True)
subprocess.run("rm Cation/Zn_Ow.txt", shell=True)
subprocess.run("rm Cation/Zn_Ot.txt", shell=True)
subprocess.run("rm Averages/AvgZnCoodOw.csv", shell=True)
subprocess.run("rm Averages/AvgZnCoodOt.csv", shell=True)
subprocess.run("rm Anion/Ot_Zn_Coodn.txt", shell=True)
subprocess.run("rm Averages/AvgOtCoodZn.csv", shell=True)
subprocess.run("rm RT/RtimesZ.txt", shell=True)
subprocess.run("rm RT/Zrt.txt", shell=True)
