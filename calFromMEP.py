import makeInputFile as makeFile
import math
import subprocess

inputData = [[1.83973621968, 3.29141438637],
             [1.85099262972, 3.14718333462],
             [1.87642105079, 2.9226516802],
             [2.01783193557, 2.68092235024],
             [2.26569401743, 2.48734916481],
             [2.48734916576, 2.26569401648],
             [2.68092235118, 2.01783193466],
             [2.92265168028, 1.87642105074],
             [3.14718333458, 1.85099262972],
             [3.29141438637, 1.83973621968],]

f = open("calFromMEP","w")
for i in range(12):
    makeFile.makeInputFile(inputData[i][0],inputData[i][1])
    subprocess.call(["zsh", "start_program.sh"])
    fileFreeEnergy = open("freeEnergy.dat", "r")
    freeEnergy = float(fileFreeEnergy.read())
    fileFreeEnergy.close()

    f.write("{0:5.4f}".format(abs(inputData[i][0]-inputData[i][1])))
    f.write("   " + "{0:3.2f}".format(freeEnergy) + "\n")