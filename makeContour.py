import  matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import  numpy as np
import variables as var

f = open("effectivePotential.dat","r")

x1d = []
y1d = []
x2d = np.zeros((var.numberOfGrit + 1,var.numberOfGrit + 1))
y2d = np.zeros((var.numberOfGrit + 1,var.numberOfGrit + 1))
effective1d = []
effective2d = np.zeros((var.numberOfGrit + 1,var.numberOfGrit + 1))


while True:
    readLine = f.readline()
    if readLine == "":
        break
    else:
        x, y, effective = map(float,readLine.split())
        x1d.append(x)
        y1d.append(y)
        effective1d.append(effective)
print(len(effective1d))
count = 0
for i in range(var.numberOfGrit + 1):
    for j in range(var.numberOfGrit + 1):
        x2d[i][j] = x1d[count]
        y2d[i][j] = y1d[count]
        effective2d[i][j] = effective1d[count]
        count += 1


print(effective2d)
interval = np.arange(-10,30,1)
CS = plt.contour(x2d, y2d, effective2d, interval)
plt.clabel(CS, inline = 1, fontsize = 10)
plt.xlim(1,6)
plt.ylim(1,6)
plt.savefig("effectivePotential.png")
plt.show()


