import  matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import  numpy as np
import variables as var

f = open("effectivePotential.dat","r")
filePhi = open("logPhi20.dat")

x1d = []
y1d = []
x2d = np.zeros((var.numberOfGrit + 1,var.numberOfGrit + 1))
y2d = np.zeros((var.numberOfGrit + 1,var.numberOfGrit + 1))
effective1d = []
effective2d = np.zeros((var.numberOfGrit + 1,var.numberOfGrit + 1))

aArray = []
bArray = []
cArray = []

a5Array = []
b5Array = []
c5Array = []


while True:
    readLine = f.readline()
    if readLine == "":
        break
    else:
        x, y, effective = map(float,readLine.split())
        x1d.append(x)
        y1d.append(y)
        effective1d.append(effective)

while True:
    readLine = filePhi.readline()
    if readLine == "":
        break
    else:
        a,b,c = map(float,readLine.split())
        aArray.append(a)
        bArray.append(b)
        cArray.append(c)



for i in range(len(aArray)):
    cArray[i] = cArray[i] - bArray[i]
    bArray[i] = bArray[i] - aArray[i]


count = 0
for i in range(var.numberOfGrit + 1):
    for j in range(var.numberOfGrit + 1):
        x2d[i][j] = x1d[count]
        y2d[i][j] = y1d[count]
        effective2d[i][j] = effective1d[count]
        count += 1

#0.14112089  2.50681886  4.87366601
#print(effective2d)
interval = np.arange(-10,50,1)
CS = plt.contour(x2d, y2d, effective2d, interval)
#plt.clabel(CS, inline = 1, fontsize = 10)
plt.plot(bArray,cArray,"o")
plt.plot(bArray,cArray,linewidth = 1,color = "m")
plt.plot(2.506-0.141,4.873-2.5068 - 2.51197913,"*",markersize = 8)
plt.xlim(1.5,5)
plt.ylim(1.5,5)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("stringInFreeEnergy.png")
plt.show()


