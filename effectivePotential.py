import potentialOfSN2 as sn2
import time
import  variables as var

numberOfGrid = var.numberOfGrit
fileFreeEnergy = open("freeEnergyContour.dat","r")
fileEffective = open("effectivePotential.dat","w")
numberOfPoints = (numberOfGrid + 1)*(numberOfGrid + 1)
print(numberOfPoints)
for i in range(numberOfPoints):
    x, y, freeEnergy = map(float,fileFreeEnergy.readline().split())
    theta = sn2.getPotential(x,y)
    theta += freeEnergy
    fileEffective.write("{0:4.3f}".format(x))
    fileEffective.write("   " + "{0:4.3f}".format(y))
    fileEffective.write("   " + str(theta) + "\n")
