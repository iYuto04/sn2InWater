import makeInputFile as makeFile
import subprocess
import matplotlib.pyplot as plt

freeEnergyZeropoint = 75.98 - 5.582
# f = open("interpNewData.dat","r")
f = open("logPhi20.dat", "r")
x = []
y = []
while True:
    readLine = f.readline()
    if readLine == "":
        break
    else:
        a,b,c = map(float,readLine.split())
        x.append((b-a) - (c-b))
        makeFile.makeInputFile(b-a, c-b)
        subprocess.call(["zsh", "start_program.sh"])
        fileFreeEnergy = open("freeEnergy.dat", "r")
        freeEnergy = float(fileFreeEnergy.read()) + freeEnergyZeropoint
        fileFreeEnergy.close()
        y.append(freeEnergy)

plt.plot(x,y,"o")
plt.plot(x,y,linewidth = 1,color = "m")
plt.savefig("plotEffectivePotential.png")
plt.show()

f.close()
