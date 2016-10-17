import numpy as np
import scipy as sp
import copy
import matplotlib.pyplot as plt
import math
import subprocess
import csv
import sys
from mpl_toolkits.mplot3d import Axes3D
#import  makeContour as cont
import variables as var

import makeInputFile as makeFile

N =100 
N_str = var.N_str #ストリングの区間の数
d_old =10000
count = 0

freeEnergyZeropoint = 75.98 - 5.582


#initialize
ab = np.linspace(0.3,4,N)
bc = np.linspace(0.3,4,N)


D = [[234.524674, 234.524674,  64.925971],
        [220.244820, 220.244820, 284.999867]]


beta = [[0.929968, 0.929968, 0.432955],
        [4.822681, 4.822681, 1.016811]]
print(beta)

r_0 =[[1.776382, 1.776382, 2.094857],
        [1.785014, 1.785014, 2.186060]]



E1_ab = np.zeros(N)
E3_ab = np.zeros(N)
E1_bc = np.zeros(N)
E3_bc = np.zeros(N)
Q_ab = np.zeros(N)
Q_bc = np.zeros(N)
J_ab = np.zeros(N)
J_bc = np.zeros(N)
theta_leps = np.zeros((N,N))
r = np.zeros((N_str+1,3))
#seibun = np.zeros((N_str+1,3,3))
nabla_r = np.zeros((3,3,3))
nabla_E = np.zeros((3,3,3,2))
nabla_theta = np.zeros((3,3))
test_nabla_theta = np.zeros(3)
E = np.zeros((3,2))
E_plus = np.zeros((3,2))
E_minus = np.zeros((3,2))
test_E = np.zeros((3,2))
dE_dr = np.zeros((3,2))
fx = np.zeros(3)
fy = np.zeros(3)
fz = np.zeros(3)
Q = np.zeros(3)
Q_plus = np.zeros(3)
Q_minus = np.zeros(3)
J = np.zeros(3)
J_plus = np.zeros(3)
J_minus = np.zeros(3)
nabla_Q = np.zeros((3,3,3))
nabla_J = np.zeros((3,3,3))
dQ_dr = np.zeros(3)
dJ_dr = np.zeros(3)
dTHETA_dr = np.zeros(3)
test_dTHETA_dr = np.zeros(3)
force = np.zeros((3,3))
rk_k = np.zeros((3,3,4))
pote_val = np.zeros(N_str+1)
plot_phi = np.zeros(N_str+1)

str_num = 0
h = 0.001 #数値微分のテスト
dt_oiler = 0.001
h = 0.0001 #数値微分のテスト 

def get_norm(x1,y1,z1,x2,y2,z2):
    norm = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    norm = norm**0.50
    return norm

def get_potential(x1,y1,z1,x2,y2,z2,x3,y3,z3):
#def get_potential(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    r1 = get_norm(x1,y1,z1,x2,y2,z2)
    r2 = get_norm(x2,y2,z2,x3,y3,z3)
    r3 = get_norm(x1,y2,z3,x3,y3,z3)
    r_val = [r1,r2,r3]
    E = np.zeros((2,3))
    for i in range(2): #singlet,triplet
        for j in range(3): #r1,r2,r3
            E[i][j] = D[i][j]*(1.0 - (-1)**(i)* math.exp(-beta[i][j]*(r_val[j] - r_0[i][j])))**2.0 - D[i][j]

    for i in range(3):
        Q[i] = (E[0][i] + E[1][i])/2.0
        J[i] = (E[0][i] - E[1][i])/2.0

    THETA = Q[0] + Q[1] + Q[2] - (J[0]*J[0]+ J[1]*J[1] + J[2]*J[2] - J[0]*J[1] - J[1]*J[2] - J[2]*J[0])**0.50 + 234.524671
    return THETA

def get_nabla_THETA(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    get_potential(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    for i in range(3):
        for j in range(3):
            nabla_theta[i][j] =0.0
    cl1 = [x1,y1,z1]
    ch3 = [x2,y2,z2]
    cl2 = [x3,y3,z3]
    r1 = get_norm(x1,y1,z1,x2,y2,z2)
    r2 = get_norm(x2,y2,z2,x3,y3,z3)
    r3 = get_norm(x1,y2,z3,x3,y3,z3)
    get_potential(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    r_val = [r1,r2,r3]
    for i in range(3): #r1,r2,r3
        for j in range(2): #singlet,triplet
            dE_dr[i][j] = (-1)**j*2.0*beta[j][i]*D[j][i]*(1.0 -(-1)**j* math.exp(-beta[j][i]*(r_val[i] - r_0[j][i])))*(math.exp(-beta[j][i]*(r_val[i] - r_0[j][i])))
        dQ_dr[i] = (dE_dr[i][0] + dE_dr[i][1])/2.0
        dJ_dr[i] = (dE_dr[i][0] - dE_dr[i][1])/2.0
    
    for i in range(3): #x,y,zに関して dr/dxなどを計算
        nabla_r[0][0][i] = -(ch3[i] - cl1[i])/r_val[0]
        nabla_r[0][2][i] = -(cl2[i] - cl1[i])/r_val[2] #この2行はClに関する微分,そのためr2に関する微分はない

        nabla_r[1][0][i] =  (ch3[i] - cl1[i])/r_val[0]
        nabla_r[1][1][i] = -(cl2[i] - ch3[i])/r_val[1] #この2行はCH3に関する微分

        nabla_r[2][1][i] =  (cl2[i] - ch3[i])/r_val[1]
        nabla_r[2][2][i] =  (cl2[i] - cl1[i])/r_val[2]

    common = 0.50*(J[0]*J[0] + J[1]*J[1] + J[2]*J[2] -J[0]*J[1] - J[1]*J[2] - J[2]*J[0])**(-0.50)

    for j in range(3):#r1,r2,r3
        dTHETA_dr[j] = dQ_dr[j] - common*dJ_dr[j]*(3.0*J[j] - J[0] - J[1] - J[2]) 

    for j in range(3):#Cl,CH3,Cl'
            for k in range(3):#x,y,z
                for l in range(3):#r1,r2,r3
                    nabla_theta[j][k] += dTHETA_dr[l]*nabla_r[j][l][k]


def plot():
    plt_x = np.zeros(N_str+1)
    plt_y = np.zeros(N_str+1)
    for i in range(N_str+1):
        plt_x[i] = get_norm(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5])
        plt_y[i] = get_norm(phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
    plt.plot(plt_x,plt_y,"o")
    plt.plot(plt_x,plt_y ,linewidth = 3,color = "m")
    plt.legend()
    #plt.savefig("test.png")
    #plt.show()
    plt_x = []
    plt_y = []

    for i in range(N_str + 1 ):
        plt_x.append(get_norm(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5]))
        plt_y.append(get_norm(phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8]))
        # plt.xlim([0,5])
        # plt.ylim([0,5])
        plt.plot(plt_x,plt_y,"o")
        plt.plot(plt_x,plt_y ,linewidth = 3,color = "m")
        plt.legend()

    #plt.savefig("test.png")

    plt.show()

def plot_saddle(count):
    plt_x = np.zeros(N_str+1)
    plt_y = np.zeros(N_str+1)
    for i in range(N_str+1):
        plt_x[i] = get_norm(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5])
        plt_y[i] = get_norm(phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
    saddle_r1 = get_norm(phi_s[0],phi_s[1],phi_s[2],phi_s[3],phi_s[4],phi_s[5])
    saddle_r2 = get_norm(phi_s[3],phi_s[4],phi_s[5],phi_s[6],phi_s[7],phi_s[8])
    saddle_plot = [saddle_r1,saddle_r2]
    print(saddle_plot)
    plt.plot(plt_x,plt_y,"o",label = "%3d count" %count)
    plt.plot(plt_x,plt_y,label = "%3d count" %count)
    plt.plot(saddle_r1,saddle_r2,"kD",label = "saddle point")
    plt.clabel(CS, inline = 1, fontsize = 10)
    plt.legend()
    plt.savefig("saddle.png")
    plt.show()


def oiler():
    for j in range(3):
        phi[i][3*j] = phi[i][3*j] - dt_oiler*test_nabla_theta[j]


def runge_kutta():
    get_nabla_THETA(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
    for j in range(3):
        for k in range(3):
            phi[i][3*j + k] = phi[i][3*j + k] - dt*nabla_theta[j][k]



def runge_kutta():
    #get_nabla_THETA(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][0] = dt*nabla_theta[j][k]
    get_nabla_THETA(phi[i][0]+0.5*rk_k[0][0][0], phi[i][1]+0.5*rk_k[0][1][0], phi[i][2]+0.5*rk_k[0][2][0],\
            phi[i][3]+0.5*rk_k[1][0][0], phi[i][4]+0.5*rk_k[1][1][0], phi[i][5]+0.5*rk_k[1][2][0],\
            phi[i][6]+0.5*rk_k[2][0][0], phi[i][7]+0.5*rk_k[2][1][0], phi[i][8]+0.5*rk_k[2][2][0])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][1] = dt*nabla_theta[j][k]  
    get_nabla_THETA(phi[i][0]+0.5*rk_k[0][0][1],phi[i][1]+0.5*rk_k[0][1][1],phi[i][2]+0.5*rk_k[0][2][1],\
            phi[i][3]+0.5*rk_k[1][0][1],phi[i][4]+0.5*rk_k[1][1][1],phi[i][5]+0.5*rk_k[1][2][1],\
            phi[i][6]+0.5*rk_k[2][0][1],phi[i][7]+0.5*rk_k[2][1][1],phi[i][8]+0.5*rk_k[2][2][1])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][2] = dt*nabla_theta[j][k]
    get_nabla_THETA(phi[i][0]+rk_k[0][0][2],phi[i][1]+rk_k[0][1][2],phi[i][2]+rk_k[0][2][2],phi[i][3]+rk_k[1][0][2],phi[i][4]+rk_k[1][1][2],phi[i][5]+rk_k[1][2][2],\
            phi[i][6]+rk_k[2][0][2],phi[i][7]+rk_k[2][1][2],phi[i][8]+rk_k[2][2][2])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][3] = dt*nabla_theta[j][k]

    for j in range(3):#cl,ch3,cl'
        for k in range(3):#x,y,z
            phi[i][3*j+k] = phi[i][3*j+k] - rk_k[j][k][0]/6.0 - rk_k[j][k][1]/3.0 - rk_k[j][k][2]/3.0 - rk_k[j][k][3]/6.0


def runge_kutta_saddle(a,b,c,d,e,f,g,h,i):
    tau_r = np.zeros(9)
    tau_r = [a,b,c,d,e,f,g,h,i]
    get_nabla_THETA(phi_s[0],phi_s[1],phi[i][2],phi_s[3],phi_s[4],phi_s[5],phi_s[6],phi_s[7],phi_s[8])
    print("tau_r")
    print(tau_r)
    for j in range(3):
        for k in range(3):
            rk_k[j][k][0] = dt*(nabla_theta[j][k] - 2.0*(nabla_theta[j][k]*tau_r[3*j+k])*tau_r[3*j+k])
    get_nabla_THETA(phi_s[0]+0.5*rk_k[0][0][0], phi_s[1]+0.5*rk_k[0][1][0], phi_s[2]+0.5*rk_k[0][2][0],\
            phi_s[3]+0.5*rk_k[1][0][0], phi_s[4]+0.5*rk_k[1][1][0], phi_s[5]+0.5*rk_k[1][2][0],\
            phi_s[6]+0.5*rk_k[2][0][0], phi_s[7]+0.5*rk_k[2][1][0], phi_s[8]+0.5*rk_k[2][2][0])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][1] = dt*(nabla_theta[j][k] - 2.0*(nabla_theta[j][k]*tau_r[3*j+k])*tau_r[3*j+k] )
    get_nabla_THETA(phi_s[0]+0.5*rk_k[0][0][1],phi_s[1]+0.5*rk_k[0][1][1],phi_s[2]+0.5*rk_k[0][2][1],\
            phi_s[3]+0.5*rk_k[1][0][1],phi_s[4]+0.5*rk_k[1][1][1],phi_s[5]+0.5*rk_k[1][2][1],\
            phi_s[6]+0.5*rk_k[2][0][1],phi_s[7]+0.5*rk_k[2][1][1],phi_s[8]+0.5*rk_k[2][2][1])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][2] = dt*(nabla_theta[j][k] - 2.0*(nabla_theta[j][k]*tau_r[3*j+k])*tau_r[3*j+k])
    get_nabla_THETA(phi_s[0]+rk_k[0][0][2],phi_s[1]+rk_k[0][1][2],phi_s[2]+rk_k[0][2][2],phi_s[3]+rk_k[1][0][2],phi_s[4]+rk_k[1][1][2],phi_s[5]+rk_k[1][2][2],\
            phi_s[6]+rk_k[2][0][2],phi_s[7]+rk_k[2][1][2],phi_s[8]+rk_k[2][2][2])
    for j in range(3):
        for k in range(3):
            rk_k[j][k][3] = dt*(nabla_theta[j][k] - 2.0*(nabla_theta[j][k]*tau_r[3*j+k])*tau_r[3*j+k])

    for j in range(3):#cl,ch3,cl'
        for k in range(3):#x,y,z
            phi_s[3*j+k] = phi_s[3*j+k] - rk_k[j][k][0]/6.0 - rk_k[j][k][1]/3.0 - rk_k[j][k][2]/3.0 - rk_k[j][k][3]/6.0
    print(phi_s)



def isou():
    for i in range(N_str+1):
        pote_val[i] = get_potential(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])

def saddle():
    naiseki = 0.0
    #count_s = 0
    #difference_s = 0.0
    #d_s_old = 10000
    V_phi = np.zeros(N_str+1)
    phi_s_old = np.zeros(9)
    tau = np.zeros(9)
    l_s = 0.0
    for i in range(N_str+1):
        V_phi[i] = get_potential(phi[i][0], phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
        phi_s_V = max(V_phi)
    
    for i in range(N_str+1):
        if V_phi[i] == phi_s_V:
            s = i+1 #ポテンシャルが最大となるインデックス
    l = s-1 #sのひとつ左の点
    r = s+1 #sのひとつ右の点
    print("s,r,l")
    print(s,r,l)
    for i in range(9):
        l_s += (phi[r][i] - phi[l][i])*(phi[r][i] - phi[l][i])
    l_s = l_s**0.50
    print("l_s is ", l_s)
    for i in range(9):
        tau[i] = (phi[r][i] - phi[l][i])/l_s
    print("in saddle function")
    print(tau)
    n_step = int(math.log(N_str**4))
    print("n_step is")
    print(n_step)
    #n_step = 12 
    for i in range(9):
        phi_s[i] = phi[s][i]
    for step in range(n_step):
    #print("phi_s")
    #print(phi_s)
    #get_nabla_THETA(phi_s[0], phi_s[1], phi_s[2], phi_s[3], phi_s[4], phi_s[5], phi_s[6], phi_s[7], phi_s[8])
    #print("saddle point nabla_theta")
    #print(nabla_theta)
    #for step in range(n_step):
        #runge_kutta_saddle(tau[0],tau[1],tau[2],tau[3],tau[4],tau[5],tau[6],tau[7],tau[8])
        get_nabla_THETA(phi_s[0], phi_s[1], phi_s[2], phi_s[3], phi_s[4], phi_s[5], phi_s[6], phi_s[7], phi_s[8])
    # #while True:
        for i in range(3): #cl,ch3,cl'
            for j in range(3): #x,y,z
                naiseki += nabla_theta[i][j]*tau[3*i+j]
        for i in range(3):
            for j in range(3):
                phi_s[3*i+j] = phi_s[3*i+j] +dt*( -nabla_theta[i][j] + 2.0*naiseki*tau[3*i+j])
                #print("変異は", dt*(-nabla_theta[i][j] + 2.0*(nabla_theta[i][j]*tau[3*i+j])*tau[3*i+j]))
        naiseki = 0.0
        #print( get_potential(phi_s[0], phi_s[1], phi_s[2], phi_s[3], phi_s[4], phi_s[5], phi_s[6], phi_s[7], phi_s[8]))
        # for i in range(9):
        #     difference_s += (phi_s[i] - phi_s_old[i])**2.0
        # difference_s = difference_s**0.50
        # d_s= difference_s/dt
        # print(d_s)

        # if  d_s > d_s_old:
        #     break
        # d_s_old = d_s
        # phi_s_old = copy.deepcopy(phi_s)
        # difference_s = 0.0
        # count_s += 1
            # print("i is" ,i)
            # print("j is",j)
            # print(nabla_theta[i][j])
    print("the potential at saddle point")
    theta_val = get_potential(phi_s[0], phi_s[1], phi_s[2], phi_s[3], phi_s[4], phi_s[5], phi_s[6], phi_s[7], phi_s[8]) 
    theta_hashi = get_potential(phi[0][0],phi[0][1],phi[0][2],phi[0][3],phi[0][4],phi[0][5],phi[0][6],phi[0][7],phi[0][8])
    print(theta_val)
    print(phi_s)
    print("barrier height is")
    print(theta_val - theta_hashi)
    # print("count_s is ")
    # print(count_s)

def plot_pote():
    plt_x = np.zeros(N_str+1)
    plt_y = np.zeros(N_str+1)
    plt_s = np.zeros(2)
    all_plot_x = np.zeros(N_str+2)
    all_plot_y = np.zeros(N_str+2)

    for i in range(N_str+1):
        plot_pote_r1 = get_norm(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5])
        plot_pote_r2 = get_norm(phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
        plt_x[i] = plot_pote_r1 - plot_pote_r2
        plt_y[i] = get_potential(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
    plot_s_r1 = get_norm(phi_s[0],phi_s[1],phi_s[2],phi_s[3],phi_s[4],phi_s[5])
    plot_s_r2 = get_norm(phi_s[3],phi_s[4],phi_s[5],phi_s[6],phi_s[7],phi_s[8])
    plt_s[0] = plot_s_r1 - plot_s_r2
    plt_s[1] = get_potential(phi_s[0],phi_s[1],phi_s[2],phi_s[3],phi_s[4],phi_s[5],phi_s[6],phi_s[7],phi_s[8])
    for i in range(N_str+1):
        all_plot_x[i] = plt_x[i]
        all_plot_y[i] = plt_y[i]
    all_plot_x[N_str+1] = plt_s[0]
    all_plot_y[N_str+1] = plt_s[1]
    print("-----------all_plot_x------------")
    print(all_plot_x)


    for i in range(N_str+1):
        for j in range(i+1,N_str+2):
            if all_plot_x[i] > all_plot_x[j]:
                change_x = all_plot_x[i]
                change_y = all_plot_y[i]
                all_plot_x[i] = all_plot_x[j]
                all_plot_y[i] = all_plot_y[j]
                all_plot_x[j] = change_x
                all_plot_y[j] = change_y
    print("-----------after sort all_plot_x------------")
    print(all_plot_x)
    plt.plot(plt_x,plt_y,"o",markersize = 13,color = "green")



    #print("plt_s_x is",plt_s_x)
    #print("plt_s_y is",plt_s_y)
    plt.plot(plt_x,plt_y,"o",markersize = 13,color = "green")
    #plt.plot(plt_x,plt_y,linewidth = 2)
    #plt.plot(plt_x,plt_y)
    plt.scatter(plt_s[0],plt_s[1],marker = "D",color = "black",s = 50)
    plt.plot(all_plot_x,all_plot_y,linewidth = 2,color = "m")
    plt.legend()
    plt.savefig("plot_pote.png")
    plt.show()

def initEndPoint():
    phi[0][3] = 1.776
    phi[0][6] = phi[0][3] + 100
    phi[N_str][3] = 100
    phi[N_str][6] = phi[N_str][3] + 1.776

def initNextPoint():
    phi[1][3] = 1.776
    phi[1][6] = phi[1][3] + 10
    phi[N_str - 1][3] = 10
    phi[N_str - 1][6] = phi[N_str - 1][3] + 1.776




# print("length")
# print(get_norm(0,0,0,1,1,0))
# print("end")
#-----------------------------------------------------------------------------------------
for i in range(N):
    E1_ab[i] = D[0][0]*(1.0 - math.exp(-beta[0][0]*(ab[i]-r_0[0][0])))**2.0 -D[0][0]
    E3_ab[i] = D[1][0]*(1.0 + math.exp(-beta[1][0]*(ab[i]-r_0[1][0])))**2.0 -D[1][0]
    Q_ab[i] = (E1_ab[i] + E3_ab[i])/2.0
    J_ab[i] = (E1_ab[i] - E3_ab[i])/2.0

    E1_bc[i] = D[0][1]*(1.0 - math.exp(-beta[0][1]*(bc[i]-r_0[0][1])))**2.0 -D[0][1]
    E3_bc[i] = D[1][1]*(1.0 + math.exp(-beta[1][1]*(bc[i]-r_0[1][1])))**2.0 -D[1][1]
    Q_bc[i] = (E1_bc[i] + E3_bc[i])/2.0
    J_bc[i] = (E1_bc[i] - E3_bc[i])/2.0

f = open("output.dat","w")

for i in range(N):
    for j in range(N):
        ca = ab[i] + bc[j]

        E1_ca = D[0][2]*(1.0 - math.exp(-beta[0][2]*(ca - r_0[0][2])))**2.0 -D[0][2]
        E3_ca = D[1][2]*(1.0 + math.exp(-beta[1][2]*(ca - r_0[1][2])))**2.0 -D[1][2]
        Q_ca = (E1_ca + E3_ca)/2.0
        J_ca = (E1_ca - E3_ca)/2.0

        theta_leps[i][j] = Q_ab[i] + Q_bc[j] + Q_ca \
                -(J_ab[i]**2.0 + J_bc[j]**2.0 + J_ca**2.0\
                -J_ab[i]*J_bc[j] - J_bc[j]*J_ca - J_ca*J_ab[i])**0.50 +234.524671


        # f.write(str(ab[i]))
        # f.write(str(" "))
        # f.write(str(bc[j]))
        # f.write(str(" "))
        # f.write(str(theta_leps[i][j]))
        # f.write("\n")
#-------------------------------------------------------------------------------------------
#print(theta_leps)



# interval = np.arange(-10,30,5)
# CS = plt.contour(ab,bc,theta_leps,interval)
# plt.clabel(CS, inline = 1, fontsize = 10)
# plt.xlim(1,4)
# plt.ylim(1,4)

# print(theta_leps)


#等高線を書くのに必要な部分---------------
interval = np.arange(-10,30,5)
CS = plt.contour(ab,bc,theta_leps,interval)
plt.clabel(CS, inline = 1, fontsize = 10)
plt.xlim(1,4)
plt.ylim(1,4)
#------------------------------------------
# plt.savefig("contour.png")
#plt.show()



#--------------------------------------------------------------------------------------------------------
#ここからストリング法
#--------------------------------------------------------------------------------------------------------
dt = 0.001
dt = 0.0001
#dt = max(N_str**(-1),0.2)*0.05
print("dt is")
print(dt)
s_sum = 0.0 #地点間の距離の計算に使用
tol = N_str**(-4)
print("tol is")
print(tol)
count = 0

#phi = [[0 for i in range (2)] for j in range(N+1)]
phi = np.zeros((N_str+1,9))
test_phi_plus = np.zeros((N_str+1,9))
#new_plus = ((N_str+1,9))
test_phi_minus = np.zeros((N_str+1,9))
#new_minus = np.zeros((N_str+1,9))
test_r_plus = np.zeros((N_str+1,3))
test_r_minus = np.zeros((N_str+1,3))
phi_old = np.zeros((N_str+1,9))
s_i = np.zeros(N_str+1)
alpha = np.zeros(N_str+1)
difference = np.zeros(N_str+1)
new_plus = ((N_str+1,9))
test_phi_minus = np.zeros((N_str+1,9))
new_minus = np.zeros((N_str+1,9))
test_r_plus = np.zeros((N_str+1,3))
test_r_minus = np.zeros((N_str+1,3))
phi_old = np.zeros((N_str+1,9))
s_i = np.zeros(N_str + 1)
#s_i = []
alpha = np.zeros(N_str+1)

difference = np.zeros(N_str+1)
#ifference = []
phi_plus  = np.zeros((N_str+1,9))
phi_minus = np.zeros((N_str+1,9))
phi_s = np.zeros(9)


phi[0][3] = 1.79
phi[0][6] = 4.0
phi[N_str][3] = 4.0
phi[N_str][6] = 1.79

for i in range(1,N_str):
    phi[i][3] = phi[0][3] + i*(phi[N_str][3] - phi[0][3])/N_str
    phi[i][6] = phi[i][3] + phi[0][6] + i*(phi[N_str][6] - phi[0][6])/N_str

phi[0][6] += phi[0][3]
phi[N_str][6] += phi[N_str][3]

# print("phi")
# print(phi)
# plot()


#for test in range(1):
while True:
#while False:

    for j in range(N_str+1):
        s_i[j] = 0.0
        difference[j] = 0.0
    s_sum = 0.0
    count += 1
#
#

    for i in range (N_str+1):
        #get_nabla_THETA(phi[i][0],phi[i][1],phi[i][2],phi[i][3],phi[i][4],phi[i][5],phi[i][6],phi[i][7],phi[i][8])
        #print(nabla_theta)

#--------------------数値微分を記述----------------------------------------------------------------------
        for j in range(3): #Cl, CH3,CL'
            new_plus = copy.deepcopy(phi)
            new_minus = copy.deepcopy(phi)

            for switch in [-1,1]:
                if switch == -1:
                    new_minus[i][3 * j] += switch * h
                    # print(new_minus)
                    theta_minus = get_potential(new_minus[i][0], new_minus[i][1], new_minus[i][2], \
                                                new_minus[i][3], new_minus[i][4], new_minus[i][5], \
                                                new_minus[i][6], new_minus[i][7],new_minus[i][8])

                    makeFile.makeInputFile(new_minus[i][3] - new_minus[i][0],new_minus[i][6] - new_minus[i][3])
                    subprocess.call(["zsh", "start_program.sh"])
                    fileFreeEnergy = open("freeEnergy.dat", "r")
                    freeEnergy = float(fileFreeEnergy.read()) + freeEnergyZeropoint
                    fileFreeEnergy.close()
                    theta_minus += freeEnergy



                elif switch == 1:
                    new_plus[i][3 * j] += switch * h

                    theta_plus = get_potential(new_plus[i][0], new_plus[i][1], new_plus[i][2], \
                                               new_plus[i][3], new_plus[i][4], new_plus[i][5], \
                                               new_plus[i][6], new_plus[i][7],new_plus[i][8])

                    makeFile.makeInputFile(new_plus[i][3] - new_plus[i][0], new_plus[i][6] - new_plus[i][3])
                    subprocess.call(["zsh", "start_program.sh"])
                    fileFreeEnergy = open("freeEnergy.dat", "r")
                    freeEnergy = float(fileFreeEnergy.read()) + freeEnergyZeropoint
                    fileFreeEnergy.close()
                    theta_plus += freeEnergy

            test_nabla_theta[j] = (theta_plus - theta_minus) / 2.0 / h




        # print("--------------------------------------------------")
        # print(test_nabla_theta)
        # input()

        oiler()
        #runge_kutta()
#     #print(phi)

    if count % 1 == 0:
        filename = "logPhi" + str(count) + ".dat"
        f = open(filename,"w")
        for i in range(N_str + 1):
            for j in range(3):
                f.write(str(phi[i][3*j]) + "  ")
            f.write(("\n"))



    s_i[0] = 0.0
    for i in range(1,N_str+1):
        for j in range(9):
            s_i[i] += (phi[i][j] - phi[i-1][j])**2.0
        s_i[i] = s_i[i]**0.50 + s_i[i-1]

    s_sum = s_i[N_str]

    for i in range(N_str + 1):
        alpha[i] = s_i[i]/s_sum



    from scipy.interpolate import interp1d
    for i in range(9):
        measured = phi[:,i]
        cubic_interp = interp1d(alpha, measured, kind = 'cubic')
        for j in range(N_str+1):
            phi[j][i] = cubic_interp(j/N_str)

    for i in range(N_str+1):
        for j in range(9):
            difference[i] += (phi[i][j] - phi_old[i][j])**2.0
        difference[i] = difference[i]**0.50

    d = max(difference)/dt




    #print("d is")
    #print(d)
    # print("after oiler method, phi")
    # print(phi)

#plot()
# # #-------------------------------------プログラム作成中はコメントアウト---------------------------------------
# #     # plt.xlim(-1.3,1.3)
# #     # plt.ylim(-0.5,1)
# #     # plt.plot(phi[:,0],phi[:,1],label = "%03step" %count)
# #     # plt.legend()
# #     # plt.savefig("graph_%03d.png" %count)
# #     # plt.clf()
# # #-- HEAD

    if d_old < d:
        print("d_old < d")
        break
    d_old = d
    if d < tol:
        print("d < tol")
        break

    print(d)

    phi_old = copy.deepcopy(phi)


    for i in range(N_str + 1):
        s_i[i] =0.0
        difference[i] = 0.0

# print("search zero-point")
# print(get_potential(0,0,0,1.776282,0,0,1000,0,0))

filename = "logLast" + str(count) + ".dat"
f = open(filename, "w")
for i in range(N_str + 1):
    for j in range(3):
        f.write(str(phi[i][3 * j]) + "  ")
    f.write(("\n"))


#plot()

isou()
print("potential")
print(pote_val)
saddle()
plot()

#isou()
print("potential")
print(pote_val)
#saddle()
#plot_saddle(count)
#plot_pote()
# print ("dt is ")
# print (dt)
# print("\n")
print("\n")
print(count)
# print("parameter for rism")
# for i in range(N_str+1):
#     print(phi[i][3]-phi[i][0], phi[i][6]-phi[i][3], phi[i][6]-phi[i][0])
print("parameter for rism")
for i in range(N_str+1):
    print(phi[i][3]-phi[i][0], phi[i][6]-phi[i][3], phi[i][6]-phi[i][0])
