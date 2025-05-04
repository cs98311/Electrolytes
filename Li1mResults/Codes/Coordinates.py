import subprocess

# Open files for coordinates at each iteration

# File stores system info like box length, no. of water, TFSI, Li, Zn
file0 = open("SysInfo.txt", "w")

file1 = open("Coordinates/H1.txt", "w")
file2 = open("Coordinates/H2.txt", "w")
file3 = open("Coordinates/OW.txt", "w")

file4 = open("Coordinates/LI.txt", "w")
file5 = open("Coordinates/ZN.txt", "w")

# Combined files for O and F
file6 = open("Coordinates/F.txt", "w")
file7 = open("Coordinates/Ot.txt", "w")

file8 = open("Coordinates/N1.txt", "w")

fileF1 = open("Coordinates/F1.txt", "w")
fileF2 = open("Coordinates/F2.txt", "w")
fileF3 = open("Coordinates/F3.txt", "w")
fileF4 = open("Coordinates/F4.txt", "w")
fileF5 = open("Coordinates/F5.txt", "w")
fileF6 = open("Coordinates/F6.txt", "w")

fileO1 = open("Coordinates/O1.txt", "w")
fileO2 = open("Coordinates/O2.txt", "w")
fileO3 = open("Coordinates/O3.txt", "w")
fileO4 = open("Coordinates/O4.txt", "w")

fileS1 = open("Coordinates/S1.txt", "w")
fileS2 = open("Coordinates/S2.txt", "w")

fileC1 = open("Coordinates/C1.txt", "w")
fileC2 = open("Coordinates/C2.txt", "w")


# Takes input the iteration number from main file
# Print progress of iterations
num = input()
print(num)


# Generate .gro file for H2O, TFSI, Li, Zn at each iteration

subprocess.run(
    "gmx trjconv -f MDSfiles/nvt3.xtc -s MDSfiles/nvt3.gro -b {} -e {} -o all.gro".format(
        num, num
    ),
    shell=True,
    input="0\n",
    capture_output=True,
    text=True,
)


# Generate coordinate files for each element from the .gro files
# Also keep count of water, TFSI, cations for the system info file

countW = 0
countA = 0
countL = 0
countZ = 0

f = open("all.gro", "r")
for line in f:
    if line[12:15] == "HW1":
        countW += 1
        print(line[23:44], file=file1)
    if line[12:15] == "HW2":
        print(line[23:44], file=file2)
    if line[13:15] == "OW":
        print(line[23:44], file=file3)
    if line[4] == ".":
        print(line[3:11], file=file0)

    if (
        line[13:15] == "F1"
        or line[13:15] == "F2"
        or line[13:15] == "F3"
        or line[13:15] == "F4"
        or line[13:15] == "F5"
        or line[13:15] == "F6"
    ):
        print(line[23:44], file=file6)

    if line[13:15] == "F1":
        print(line[23:44], file=fileF1)
    if line[13:15] == "F2":
        print(line[23:44], file=fileF2)
    if line[13:15] == "F3":
        print(line[23:44], file=fileF3)
    if line[13:15] == "F4":
        print(line[23:44], file=fileF4)
    if line[13:15] == "F5":
        print(line[23:44], file=fileF5)
    if line[13:15] == "F6":
        print(line[23:44], file=fileF6)

    if line[13:15] == "S1":
        print(line[23:44], file=fileS1)

    if line[13:15] == "S2":
        print(line[23:44], file=fileS2)

    if line[13:15] == "C1":
        print(line[23:44], file=fileC1)

    if line[13:15] == "C2":
        print(line[23:44], file=fileC2)

    if (
        line[13:15] == "O1"
        or line[13:15] == "O2"
        or line[13:15] == "O3"
        or line[13:15] == "O4"
    ):
        print(line[23:44], file=file7)

    if line[13:15] == "O1":
        print(line[23:44], file=fileO1)
    if line[13:15] == "O2":
        print(line[23:44], file=fileO2)
    if line[13:15] == "O3":
        print(line[23:44], file=fileO3)
    if line[13:15] == "O4":
        print(line[23:44], file=fileO4)

    if line[13:15] == "N1":
        countA += 1
        print(line[23:44], file=file8)

    if line[13:15] == "LI":
        countL += 1
        print(line[23:44], file=file4)

    if line[13:15] == "ZN":
        countZ += 1
        print(line[23:44], file=file5)


print(countW, file=file0)
print(countA, file=file0)
print(countL, file=file0)
print(countZ, file=file0)


# Remove the .gro file at the end of iteration to avoid pileup

subprocess.run("rm all.gro", shell=True)


# Close the opened files

file0.close()
file1.close()
file2.close()
file3.close()
file4.close()
file5.close()
file6.close()
file7.close()
file8.close()
f.close()
