import potentialOfSN2 as sn2
import time

t0 = time.clock()

numberOfGrid = 30 #main.pyのnumberOfGridの数
#numberOfGrid = 2
fileFreeEnergy = open("forContour.dat","r")
#fileFreeEnergy = open("test.dat","r")
fileEffective = open("effectivePotential.dat","w")
numberOfPoints = (numberOfGrid + 1)*(numberOfGrid + 1)
print(numberOfPoints)
for i in range(numberOfPoints):
    x, y, freeEnergy = map(float,fileFreeEnergy.readline().split())
    theta = sn2.getPotential(x,y)
    theta += (freeEnergy - 63.27)
    fileEffective.write("{0:4.3f}".format(x))
    fileEffective.write("   " + "{0:4.3f}".format(y))
    fileEffective.write("   " + str(theta) + "\n")

t1 = time.clock()
diff = t1 - t0
h = 0
m = 0
if diff >= 3600:
    h = int(diff / 3600)
    diff %=  3600
if diff > 60:
    m = int(diff / 60)
    diff %= 60
print("Time  "+ str(h) + "[h]"  + str(m) + "[m]" + str(diff) + "[s]")

