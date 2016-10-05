#JCP 86,1354 "Molecular dynamis of a model S N 2 reaction in water"のモデルを使いました
#モデルは1次元上に並んでいて間の長さのみでポテンシャルが決まる
import numpy as np
import scipy as sp
import math

def getPotential(ab,bc):
    D = [[234.524674, 234.524674,  64.925971],
        [220.244820, 220.244820, 284.999867]]
    beta = [[0.929968, 0.929968, 0.432955],
        [4.822681, 4.822681, 1.016811]]
    r_0 =[[1.776382, 1.776382, 2.094857],
        [1.785014, 1.785014, 2.186060]]
    rVal = [ab,bc,ab + bc]
    E = np.zeros((2,3)) #(singlet及びtriplet)*(rab,rbc,rac)
    Q = np.zeros(3)
    J = np.zeros(3)
    for i in [0,1]:
        for j in [0,1,2]:
            E[i][j] = D[i][j]*(1.0 - (-1)**(i)* math.exp(-beta[i][j]*(rVal[j] - r_0[i][j])))**2.0 - D[i][j]
    for i in [0,1,2]:
        Q[i] = (E[0][i] + E[1][i])/2.0
        J[i] = (E[0][i] - E[1][i])/2.0
    THETA = Q[0] + Q[1] + Q[2] - (J[0]*J[0]+ J[1]*J[1] + J[2]*J[2] - J[0]*J[1] - J[1]*J[2] - J[2]*J[0])**0.50 + 234.524671
    return THETA

#ポテンシャル地形がうまくできているか確認するためのテスト
def test():
    startPosition = 1.0
    endPosition = 4.0
    dt = 100
    dtOfLine = (endPosition - startPosition) / dt
    dtOfRow = dtOfLine
    f = open("potentialOfAir.dat","w")
    for i in range(dt + 1):
        for j in range(dt +1):
            x = startPosition + dtOfLine * i
            y = startPosition + dtOfRow * j
            potential= getPotential(x,y)
            f.write("{0:2.1f}".format(x))
            f.write("   ")
            f.write("{0:2.1f}".format(y))
            f.write("   ")
            f.write(str(potential))
            f.write("\n")
    f.close()

#test()
#print("成功")
