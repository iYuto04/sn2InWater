import  matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import  numpy as np
import variables as var
print(var.N_str)
input()
f = open("effectivePotential.dat","r")
filePhi = open("logLast25.dat")

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
f.close()

while True:
    readLine = filePhi.readline()
    if readLine == "":
        break
    else:
        a,b,c = map(float,readLine.split())
        aArray.append(a)
        bArray.append(b)
        cArray.append(c)
filePhi.close()
print(aArray)
input()
alpha = []
s_i = [0 for i in range(var.N_str + 1)]
newArrayOfA = []
newArrayOfB = []
newArrayOfC = []
s_i[0] = 0.0
for i in range(1,var.N_str+1):
    for j in range(3):
        s_i[i] += (aArray[i] - aArray[i-1])**2.0
        s_i[i] += (bArray[i] - bArray[i-1])**2.0
        s_i[i] += (cArray[i] - cArray[i-1])**2.0
    s_i[i] = s_i[i]**0.50 + s_i[i-1]

s_sum = s_i[var.N_str]

for i in range(var.N_str + 1):
    alpha.append(s_i[i]/s_sum)
print(alpha)
addAlpha = []
beta = []
start = 0.30
dt = 0.02
step = int((1 - 2*start) / dt)
# print(step)
# input()
for i in range(1,step + 1):
    addAlpha.append(start + dt*i)
for i in range(var.N_str + 1):
    beta.append(i/var.N_str)
# print(beta)
addAlpha.extend(beta)
addAlpha.sort()
# print(addAlpha)
# print(len(addAlpha))
# input()

print(addAlpha)
from scipy.interpolate import interp1d
for i in range(3):
    if i == 0:
        cubic_interp = interp1d(alpha, aArray, kind = 'cubic')
        for j in range(len(addAlpha)):
            newArrayOfA.append(cubic_interp(addAlpha[j]))
    elif i == 1:
        cubic_interp = interp1d(alpha, bArray, kind='cubic')
        for j in range(len(addAlpha)):
            newArrayOfB.append(cubic_interp(addAlpha[j]))
    elif i== 2:
        cubic_interp = interp1d(alpha, cArray, kind='cubic')
        for j in range(len(addAlpha)):
            newArrayOfC.append(cubic_interp(addAlpha[j]))
    else:
        print("error!")
        input()

f = open("interpNewData.dat","w")
for i in range(len(addAlpha)):
    f.write(str(newArrayOfA[i]) + "  ")
    f.write(str(newArrayOfB[i]) + "  ")
    f.write(str(newArrayOfC[i]) + "   \n")

for i in range(len(aArray)):
    cArray[i] = cArray[i] - bArray[i]
    bArray[i] = bArray[i] - aArray[i]
for i in range(len(addAlpha)):
    newArrayOfC[i] = newArrayOfC[i] - newArrayOfB[i]
    newArrayOfB[i] = newArrayOfB[i] - newArrayOfA[i]


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

plt.plot(newArrayOfB,newArrayOfC,linewidth = 1,color = "m")
plt.plot(newArrayOfB,newArrayOfC,"o")

#plt.plot(bArray,cArray,"o")
#plt.plot(bArray,cArray,linewidth = 1,color = "m")

plt.plot(2.506-0.141,4.873-2.5068 - 2.51197913,"*",markersize = 8)
plt.xlim(1.5,5)
plt.ylim(1.5,5)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("stringInFreeEnergy.png")
plt.show()





